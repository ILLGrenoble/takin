#!/bin/env python3
#
# hdf5 tas file converter
# @author Tobias Weber <tweber@ill.fr>
# @author Gabriel Rebolini <rebolinig@ill.fr>
# @date 18-oct-2022
# @license GPLv2
#

import h5py
import numpy as np
import tabulate as tab
import re
import os
import sys
import argparse

print_retro = True
print_statistics = False


class H5Loader:
	#
	# get data out of an hdf5 entry
	#
	@staticmethod
	def get_dat(entry, path):
		try:
			return entry[path][0]
		except KeyError:
			return None


	#
	# get string data out of an hdf5 entry
	#
	@staticmethod
	def get_str(entry, path):
		dat = H5Loader.get_dat(entry, path)
		if dat == None:
			return ""
		return dat.decode("utf-8")


	#
	# load a TAS nexus file
	#
	def __init__(self, filename):
		file = h5py.File(filename, "r")
		entry = file["entry0"]

		# get scan data
		self.data = entry["data_scan/scanned_variables/data"][:]
		self.data = np.transpose(self.data)

		try:
			self.columns = entry["data_scan/scanned_variables/variables_names/label"][:]
		except KeyError:
			axes = entry["data_scan/scanned_variables/variables_names/axis"][:]
			names = entry["data_scan/scanned_variables/variables_names/name"][:]
			properties = entry["data_scan/scanned_variables/variables_names/property"][:]
			self.columns = [names[i] if axes[i]!=0 else properties[i] for i in range(axes.size)]
		self.columns = np.array([str.decode("utf-8") for str in self.columns])
		try:
			# select scanned columns
			scanned_cols = entry["data_scan/scanned_variables/variables_names/scanned"][:]
			self.selected_columns = [self.columns[idx] for idx in range(scanned_cols.size) if scanned_cols[idx] != 0]
		except KeyError:
			# select all columns
			self.selected_columns = self.columns

		# add data row index
		if not "PNT" in self.columns:
			num_rows = len(self.data)
			row_indices = np.linspace(1, num_rows, num_rows)
			self.data = np.append(self.data, row_indices.reshape(num_rows, 1), axis=1)
			self.columns = np.append(self.columns, "PNT")
			self.selected_columns.insert(0, "PNT")

		# add detector, monitor, tim, TT and TRT columns
		re_det = re.compile("([A-Za-z0-9]*)(Detector|Monitor)([A-Za-z0-9]*)|^Time$|^TT$|^TRT$")
		for col_name in self.columns:
			if col_name in self.selected_columns:
				continue
			if re_det.match(col_name) == None:
				continue
			self.selected_columns.append(col_name)

		# get instrument variables
		self.varias = {}
		self.zeros = {}
		self.targets = {}

		# find the instrument group
		instr_name = "instrument"
		if "instrument_name" in entry:
			real_instr_name = self.get_str(entry, "instrument_name")

			# check if there's an instrument group with that name
			if real_instr_name in entry:
				instr_name = real_instr_name

		# no instrument with the given name available?
		if not instr_name in entry:
			# get first group that is marked with "NXinstrument"
			for cur_entry in entry:
				nx_cls = entry[cur_entry].attrs.get("NX_class")
				if nx_cls != None and nx_cls.decode("utf-8") == "NXinstrument":
					instr_name = cur_entry
					break

		instr = entry[instr_name]
		self.instrname = self.get_str(instr, "name")
		self.commandline = self.get_str(instr, "command_line/actual_command")
		self.palcmd = self.get_str(instr, "pal/pal_contents")
		self.instrmode = self.get_str(entry, "instrument_mode")
		self.mono_d = self.get_dat(instr, "Monochromator/d_spacing")
		self.mono_k = self.get_dat(instr, "Monochromator/ki")
		self.mono_sense = self.get_dat(instr, "Monochromator/sens")
		self.mono_mosaic = self.get_dat(instr, "Monochromator/mosaic")
		self.reactor = self.get_dat(instr, "source/power")

		if self.get_dat(instr, "Monochromator/automatic_curvature"):
			self.mono_autocurve = "auto"
		else:
			self.mono_autocurve = "manu"
		self.ana_d = self.get_dat(instr, "Analyser/d_spacing")
		self.ana_k = self.get_dat(instr, "Analyser/kf")
		self.ana_sense = self.get_dat(instr, "Analyser/sens")
		self.ana_mosaic = self.get_dat(instr, "Analyser/mosaic")
		if self.get_dat(instr, "Analyser/automatic_curvature"):
			self.ana_autocurve = "auto"
		else:
			self.ana_autocurve = "manu"
		for key in instr.keys():
			varia_path = key + "/value"
			offs_path = key + "/offset_value"
			target_path = key + "/target_value"

			if varia_path in instr:
				self.varias[key] = self.get_dat(instr, varia_path)
			if offs_path in instr:
				self.zeros[key] = self.get_dat(instr, offs_path)
			if target_path in instr:
				self.targets[key] = self.get_dat(instr, target_path)

		self.colli_h = [
			self.get_dat(instr, "Distance/alf1"),
			self.get_dat(instr, "Distance/alf2"),
			self.get_dat(instr, "Distance/alf3"),
			self.get_dat(instr, "Distance/alf4"),
		]

		self.colli_v = [
			self.get_dat(instr, "Distance/bet1"),
			self.get_dat(instr, "Distance/bet2"),
			self.get_dat(instr, "Distance/bet3"),
			self.get_dat(instr, "Distance/bet4"),
		]

		# get user infos
		user = entry["user"]
		self.username = self.get_str(user, "name")
		self.localname = self.get_str(user, "namelocalcontact")
		self.expnumber = self.get_str(user, "proposal")

		# get experiment infos
		self.exptitle = self.get_str(entry, "title")
		self.starttime = self.get_str(entry, "start_time")
		self.endtime = self.get_str(entry, "end_time")
		self.numor = self.get_dat(entry, "run_number")

		# get and calculate steps infos
		qh = entry["data_scan/scanned_variables/data"][0]
		qh_st = len(qh) -1
		self.qh_step = (qh[qh_st] - qh[0]) / qh_st
		qk = entry["data_scan/scanned_variables/data"][1]
		qk_st = len(qk) -1
		self.qk_step = (qk[qk_st] - qk[0]) / qk_st

		# get sample infos
		sample = entry["sample"]
		self.gonio = self.get_dat(sample, "automatic_gonio")
		self.posqe = (
			self.get_dat(sample, "qh"),
			self.get_dat(sample, "qk"),
			self.get_dat(sample, "ql"),
			self.get_dat(sample, "en") )
		self.lattice = (
			self.get_dat(sample, "unit_cell_a"),
			self.get_dat(sample, "unit_cell_b"),
			self.get_dat(sample, "unit_cell_c") )
		self.angles = (
			self.get_dat(sample, "unit_cell_alpha"),
			self.get_dat(sample, "unit_cell_beta"),
			self.get_dat(sample, "unit_cell_gamma") )
		self.plane0 = ( self.get_dat(sample, "ax"), self.get_dat(sample, "ay"), self.get_dat(sample, "az") )
		self.plane1 = ( self.get_dat(sample, "bx"), self.get_dat(sample, "by"), self.get_dat(sample, "bz") )
		self.sample_sense = self.get_dat(sample, "sens")
		self.sample_mosaic = self.get_dat(sample, "mosaic")

		self.kfix_which = self.get_dat(sample, "fx")
		if self.kfix_which == 2:
			self.kfix = self.ana_k
		else:
			self.kfix = self.mono_k

		self.temp_t = self.get_dat(sample, "temperature")
		self.temp_r = self.get_dat(sample, "regulation_temperature")
		self.mag_field = self.get_dat(sample, "additional_environment/MagneticField/field")
		if self.temp_t == None:
			self.temp_t = -1.
		if self.temp_r == None:
			self.temp_r = -1.
		if self.mag_field == None:
			self.mag_field = -1.



	#
	# write a table of the scanned variables
	#
	def write_table(self, outfile = sys.stdout, table_format = "rounded_grid"):
		indices = np.array([np.where(self.columns == selected_column)[0][0] \
			for selected_column in self.selected_columns])
		outfile.write(tab.tabulate(self.data[:,indices], self.columns[indices],
			numalign = "right", tablefmt = table_format))


	#
	# write the old-style TAS text file format
	#
	def write_retro(self, f = sys.stdout):
		# write header variable
		def write_var(var, name, f):
			ctr = 0
			for key in var:
				if ctr % 4 == 0 and ctr != 0:
					f.write("\n")
				if ctr % 4 == 0:
					f.write("%s: " % name)
				val = float(var[key])
				f.write("%-8s = %6.2f, " % (key, val))
				ctr += 1

		# write header
		f.write("INSTR: %s\n" % self.instrname)
		f.write("EXPNO: %s\n" % self.expnumber)
		f.write("USER_: %s\n" % self.username)
		f.write("LOCAL: %s\n" % self.localname)
		f.write("FILE_: %d\n" % self.numor)
		f.write("DATE_: %s\n" % self.starttime)
		f.write("TITLE: %s\n" % self.exptitle)
		f.write("TYPE_: %s\n" % self.instrmode)
		f.write("COMND: %s\n" % self.commandline)
		f.write("POSQE: QH = %.4f, QK = %.4f, QL = %.4f, EN = %.4f, UN=meV\n" % self.posqe)
		f.write("CURVE: MONO = %s, ANA = %s\n" % (self.mono_autocurve, self.ana_autocurve))
		f.write("STEPS: QH = %.4f, QK = %.4f\n" % (self.qh_step, self.qk_step))
		f.write("PARAM: GONIO = %s\n" % self.gonio)
		f.write("PARAM: DM = %.5f, DA = %.5f, KFIX = %.5f\n" % (self.mono_d, self.ana_d, self.kfix))
		f.write("PARAM: SM = %d, SS = %d, SA = %d, FX = %d\n" % (self.mono_sense, self.sample_sense, self.ana_sense, self.kfix_which))
		if self.colli_h[0] != None:
			f.write("PARAM: ALF1 = %.2f, ALF2 = %.2f, ALF3 = %.2f, ALF4 = %.2f\n" % (self.colli_h[0], self.colli_h[1], self.colli_h[2], self.colli_h[3]))
		if self.colli_v[0] != None:
			f.write("PARAM: BET1 = %.2f, BET2 = %.2f, BET3 = %.2f, BET4 = %.2f\n" % (self.colli_v[0], self.colli_v[1], self.colli_v[2], self.colli_v[3]))
		f.write("PARAM: ETAM = %.2f, ETAS = %.5f, ETAA = %.2f\n" % (self.mono_mosaic, self.sample_mosaic, self.ana_mosaic))
		f.write("PARAM: AS = %.5f, BS = %.5f, CS = %.5f\n" % self.lattice)
		f.write("PARAM: AA = %.5f, BB = %.5f, CC = %.5f\n" % self.angles)
		f.write("PARAM: AX = %.3f, AY = %.3f, AZ = %.3f\n" % self.plane0)
		f.write("PARAM: BX = %.3f, BY = %.3f, BZ = %.3f\n" % self.plane1)
		f.write("PARAM: TT = %.4f, RT = %.4f, MAG = %.6f\n" % (self.temp_t, self.temp_r, self.mag_field))
		f.write("PARAM: REACTOR = %s\n" % self.reactor)
		write_var(self.varias, "VARIA", f)
		write_var(self.zeros, "ZEROS", f)
		write_var(self.targets, "TARGET", f)
		f.write("\n")
		for polcmd in self.palcmd.split("|"):
			polcmd = polcmd.strip()
			if polcmd != "":
				f.write("POLAN: %s" % polcmd)

		# write data
		f.write("FORMT:\n")  # TODO
		f.write("DATA_:\n")
		self.write_table(f, table_format = "plain")
		f.write("\n")


	#
	# get some statistics about the measurement
	#
	def get_statistics(self):
		indices = np.array([np.where(self.columns == "Time")])

		count_time = 0.
		for time in self.data[:,indices]:
			count_time += time[0][0][0]

		import datetime as dt

		start = dt.datetime.strptime(self.starttime, "%d-%b-%y %H:%M:%S")
		end = dt.datetime.strptime(self.endtime, "%d-%b-%y %H:%M:%S")
		scan_duration = end - start

		scan_time = scan_duration.total_seconds()
		move_time = scan_time - count_time

		return [ scan_duration, count_time, move_time ]


	#
	# write some statistics about the measurement
	#
	def write_statistics(self, title, scan_duration, count_time, move_time, f = sys.stdout):
		scan_time = scan_duration.total_seconds()

		f.write("\nTotal time needed for %s:\n" % title)
		if self != None:
			f.write("\tScan start time:          %s\n" % self.starttime)
			f.write("\tScan stop time:           %s\n" % self.endtime)
		f.write("\tScan time:                %d s = %.2f min = %.2f h = %s\n" % \
			(scan_time, scan_time / 60., scan_time / 3600., str(scan_duration)))
		f.write("\tActual counting time:     %.2f s = %.3f min = %.4f h = %.2f %%\n" % \
			(count_time, count_time / 60., count_time / 3600., count_time / scan_time * 100.))
		f.write("\tInstrument movement time: %.2f s = %.3f min = %.4f h = %.2f %%\n" % \
			(move_time, move_time / 60., move_time / 3600., move_time / scan_time * 100.))
		f.write("\n")


#
# loads TAS files from the command line and converts them
#
def main(argv):
	total_scan_duration = None
	total_count_time = 0.
	total_move_time = 0.
	files = []

	parser = argparse.ArgumentParser(description=".nxs files treatment")
	parser.add_argument("input", nargs="+", help="input paths of .nxs files and folders")
	parser.add_argument("-p", "--print", action="store_true", help="print converted file to console")
	parser.add_argument("-s", "--statistics", action="store_true", help="show statistics")
	args = parser.parse_args()
	input_arg = args.input
	print_statistics = args.statistics

	for argname in input_arg:
		nxs_file = re.compile(".*\\.nxs$")
		if os.path.isdir(argname):
			for file_name in os.listdir(argname):
				if nxs_file.match(file_name):
					files.append(os.path.join(argname, file_name))
		else:
			if nxs_file.match(argname):
				files.append(argname)

	def convert_file(h5, outfile = sys.stdout):
		if print_retro:
			#h5.selected_columns = [ "QH", "QK", "QL", "EN" ]
			h5.write_retro(outfile)

		if print_statistics:
			[scan_duration, count_time, move_time] = h5.get_statistics()
			h5.write_statistics("scan %s" % h5.numor, scan_duration, count_time, move_time)

			nonlocal total_scan_duration, total_count_time, total_move_time
			if total_scan_duration == None:
				total_scan_duration = scan_duration
			else:
				total_scan_duration += scan_duration
			total_count_time += count_time
			total_move_time += move_time

	for filename in files:
		try:
			h5 = H5Loader(filename)

			if args.print:
				convert_file(h5)
			else:
				output_name = filename[0:-3] + "dat"
				with open(output_name, "w") as outfile:
					print(filename + " -> " + output_name)
					convert_file(h5, outfile)
				
		except FileNotFoundError as err:
			print(err, file = sys.stderr)

	if print_statistics and len(files) > 1:
		H5Loader.write_statistics(None, "all scans" , total_scan_duration, total_count_time, total_move_time)



if __name__ == "__main__":
	import sys
	main(sys.argv)
