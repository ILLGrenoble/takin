#!/bin/env python3
#
# hdf5 tas file converter
# @author Tobias Weber <tweber@ill.fr>
# @date 18-oct-2022
# @license GPLv2
#

import h5py
import numpy as np
import tabulate as tab
import re


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

		# add detector and monitor columns
		re_det = re.compile("([A-Za-z0-9]*)(Detector|Monitor)([A-Za-z0-9]*)")
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
		if instr["Monochromator/automatic_curvature"]:
			self.mono_autocurve = "auto"
		else:
			self.mono_autocurve = "manu"
		self.ana_d = self.get_dat(instr, "Analyser/d_spacing")
		self.ana_k = self.get_dat(instr, "Analyser/kf")
		self.ana_sense = self.get_dat(instr, "Analyser/sens")
		self.ana_mosaic = self.get_dat(instr, "Analyser/mosaic")
		if instr["Analyser/automatic_curvature"]:
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

		# get sample infos
		sample = entry["sample"]
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


	#
	# prints a table of the scanned variables
	#
	def print_table(self, table_format = "rounded_grid"):
		indices = np.array([np.where(self.columns == selected_column)[0][0] \
			for selected_column in self.selected_columns])
		print(tab.tabulate(self.data[:,indices], self.columns[indices],
			numalign = "right", tablefmt = table_format))


	#
	# prints the old-style TAS text file format
	#
	def print_retro(self):
		# print header variable
		def print_var(var, name):
			ctr = 0
			for key in var:
				if ctr % 4 == 0:
					print("\n%s: " % name, end="")
				val = float(var[key])
				print("%-8s = %6.2f, " % (key, val), end="")
				ctr += 1

		# print header
		print("INSTR: %s" % self.instrname)
		print("EXPNO: %s" % self.expnumber)
		print("USER_: %s" % self.username)
		print("LOCAL: %s" % self.localname)
		print("FILE_: %d" % self.numor)
		print("DATE_: %s" % self.starttime)
		print("TITLE: %s" % self.exptitle)
		print("TYPE_: %s" % self.instrmode)
		print("COMND: %s" % self.commandline)
		print("POSQE: QH = %.4f, QK = %.4f, QL = %.4f, EN = %.4f, UN=meV" % self.posqe)
		print("CURVE: MONO = %s, ANA = %s" % (self.mono_autocurve, self.ana_autocurve))
		print("PARAM: DM = %.5f, DA = %.5f, KFIX = %.5f" % (self.mono_d, self.ana_d, self.kfix))
		print("PARAM: SM = %d, SS = %d, SA = %d, FX = %d" % (self.mono_sense, self.sample_sense, self.ana_sense, self.kfix_which))
		if self.colli_h[0] != None:
			print("PARAM: ALF1 = %.2f, ALF2 = %.2f, ALF3 = %.2f, ALF4 = %.2f" % (self.colli_h[0], self.colli_h[1], self.colli_h[2], self.colli_h[3]))
		if self.colli_v[0] != None:
			print("PARAM: BET1 = %.2f, BET2 = %.2f, BET3 = %.2f, BET4 = %.2f" % (self.colli_v[0], self.colli_v[1], self.colli_v[2], self.colli_v[3]))
		print("PARAM: ETAM = %.2f, ETAS = %.5f, ETAA = %.2f" % (self.mono_mosaic, self.sample_mosaic, self.ana_mosaic))
		print("PARAM: AS = %.5f, BS = %.5f, CS = %.5f" % self.lattice)
		print("PARAM: AA = %.5f, BB = %.5f, CC = %.5f" % self.angles)
		print("PARAM: AX = %.3f, AY = %.3f, AZ = %.3f" % self.plane0)
		print("PARAM: BX = %.3f, BY = %.3f, BZ = %.3f" % self.plane1)
		print_var(self.varias, "VARIA")
		print_var(self.zeros, "ZEROS")
		print_var(self.targets, "TARGET")
		print()
		for polcmd in self.palcmd.split("|"):
			polcmd = polcmd.strip()
			if polcmd != "":
				print("POLAN: %s" % polcmd)

		# print data
		print("FORMT:")  # TODO
		print("DATA_:")
		self.print_table(table_format = "plain")


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
	# prints some statistics about the measurement
	#
	def print_statistics(self, title, scan_duration, count_time, move_time):
		scan_time = scan_duration.total_seconds()

		print("\nTotal time needed for %s:" % title)
		if self != None:
			print("\tScan start time:          %s" % self.starttime)
			print("\tScan stop time:           %s" % self.endtime)
		print("\tScan time:                %d s = %.2f min = %.2f h = %s" % \
			(scan_time, scan_time / 60., scan_time / 3600., str(scan_duration)))
		print("\tActual counting time:     %.2f s = %.3f min = %.4f h = %.2f %%" % \
			(count_time, count_time / 60., count_time / 3600., count_time / scan_time * 100.))
		print("\tInstrument movement time: %.2f s = %.3f min = %.4f h = %.2f %%" % \
			(move_time, move_time / 60., move_time / 3600., move_time / scan_time * 100.))


#
# loads TAS files from the command line and converts them
#
def main(argv):

	total_scan_duration = None
	total_count_time = 0.
	total_move_time = 0.
	files = argv[1:]

	for filename in files:
		try:
			h5 = H5Loader(filename)

			if print_retro:
				#h5.selected_columns = [ "QH", "QK", "QL", "EN" ]
				h5.print_retro()

			if print_statistics:
				[scan_duration, count_time, move_time] = h5.get_statistics()
				h5.print_statistics("scan %s" % h5.numor, scan_duration, count_time, move_time)

				if total_scan_duration == None:
					total_scan_duration = scan_duration
				else:
					total_scan_duration += scan_duration
				total_count_time += count_time
				total_move_time += move_time
		except FileNotFoundError as err:
			print(err, file = sys.stderr)

	if print_statistics and len(files) > 1:
		H5Loader.print_statistics(None, "all scans" , total_scan_duration, total_count_time, total_move_time)


if __name__ == "__main__":
	import sys
	main(sys.argv)
