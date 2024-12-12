#
# magdyn py interface demo -- calculating a dispersion surface
# @author Tobias Weber <tweber@ill.fr>
# @date 12-dec-2024
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

import numpy

try:
	# try importing packaged version
	from takin import magdyn
except ModuleNotFoundError:
	# try importing unpackaged version
	import magdyn


# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
print_dispersion       = False  # write dispersion to console
plot_dispersion        = True   # show dispersion plot
plot_file              = ""     # file to save plot to
only_positive_energies = True   # ignore magnon annihilation?
num_Q_points           = 32     # number of Qs per direction to calculate on a dispersion direction

S_scale                = 64.    # weight scaling and clamp factors
S_clamp_min            = 1.     #
S_clamp_max            = 500.   #
S_filter_min           = -1.    # don't filter

modelfile              = "model.magdyn"
max_threads            = 0      # number of worker threads, 0: automatic determination

# dispersion plotting range
hkl_start              = numpy.array([ 0., 0., 0.0 ])
hkl_end1               = numpy.array([ 1., 1., 0.0 ])
hkl_end2               = numpy.array([ 0., 0., 0.5 ])
Q_idx1                 = 0
Q_idx2                 = 2

# plot options
alpha                  = 0xc0   # transparency
xlabel                 = "hh0 (rlu)"
ylabel                 = None
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# load the magnetic model
# -----------------------------------------------------------------------------
# create the magdyn object
mag = magdyn.MagDyn()


# load the model file
print("Loading {}...".format(modelfile))
if mag.Load(modelfile):
	print("Loaded {}.".format(modelfile))
else:
	print("Failed loading {}.".format(modelfile))
	exit(-1)


# minimum energy
print("\nEnergy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate the 2d dispersion surface
# -----------------------------------------------------------------------------
def calc_disp():
	print("\nCalculating dispersion surface...")
	if print_dispersion:
		print("{:>15} {:>15} {:>15} {:>15} {:>15} {:>8}".format(
			"h", "k", "l", "E", "S(Q,E)", "index"))

	data_h =   []
	data_k =   []
	data_l =   []
	data_E =   []
	data_S =   []
	data_idx = []

	#
	# add a data point
	#
	def append_data(h, k, l, E, S, idx):
		weight = S * S_scale

		if weight < S_filter_min:
			return
		if weight < S_clamp_min:
			weight = S_clamp_min
		elif weight > S_clamp_max:
			weight = S_clamp_max

		data_h.append(h)
		data_k.append(k)
		data_l.append(l)
		data_E.append(E)
		data_S.append(weight)
		data_idx.append(idx)


	# calculate the dispersion
	data_disp = mag.CalcDispersion(
		hkl_start[0], hkl_start[1], hkl_start[2],
		hkl_end1[0], hkl_end1[1], hkl_end1[2],
		hkl_end2[0], hkl_end2[1], hkl_end2[2],
		num_Q_points, max_threads)

	for S in data_disp:
		branch_idx = 0

		for data_EandS in S.E_and_S:
			if only_positive_energies and data_EandS.E < 0.:
				continue

			append_data(
				magdyn.get_h(S), magdyn.get_k(S), magdyn.get_l(S),
				data_EandS.E, data_EandS.weight, branch_idx)

			if print_dispersion:
				print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g} {:8d}".format(
					magdyn.get_h(S), magdyn.get_k(S), magdyn.get_l(S),
					data_EandS.E, data_EandS.weight, branch_idx))

			branch_idx += 1

	return [ data_h, data_k, data_l, data_E, data_S, data_idx ]
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the 2d dispersion surfaces
# -----------------------------------------------------------------------------
def plot_disp_2d(data, Q_idx1 = -1, Q_idx2 = -1):
	import matplotlib.pyplot as pyplot
	pyplot.rcParams.update({
		"font.sans-serif" : "DejaVu Sans",
		"font.family" : "sans-serif",
		"font.size" : 12,
	})

	data = data.T

	# find x and y axis to plot
	if Q_idx1 < 0 or Q_idx2 < 0:
		# Q start and end points
		b1 = ( data[:, 0][0], data[:, 1][0], data[:, 2][0] )
		b2 = ( data[:, 0][-1], data[:, 1][-1], data[:, 2][-1] )

		# find scan axis
		Q_diff = [
			numpy.abs(b1[0] - b2[0]),
			numpy.abs(b1[1] - b2[1]),
			numpy.abs(b1[2] - b2[2]) ]

		# first scan axis
		if Q_idx1 < 0:
			if Q_diff[1] > 0. and Q_idx2 != 0:
				Q_idx1 = 0
			if Q_diff[1] > Q_diff[Q_idx1] and Q_idx2 != 1:
				Q_idx1 = 1
			elif Q_diff[2] > Q_diff[Q_idx1] and Q_idx2 != 2:
				Q_idx1 = 2

		# second scan axis
		if Q_idx2 < 0:
			if Q_diff[1] > 0. and Q_idx1 != 0:
				Q_idx2 = 0
			if Q_diff[1] > Q_diff[Q_idx2] and Q_idx1 != 1:
				Q_idx2 = 1
			elif Q_diff[2] > Q_diff[Q_idx2] and Q_idx1 != 2:
				Q_idx2 = 2

	labels = [ "h (rlu)", "k (rlu)", "l (rlu)" ]

	(plt, axis) = pyplot.subplots(nrows = 1, ncols = 1,
		width_ratios = None, sharey = True,
		subplot_kw = { "projection" : "3d" })

	E_branch_idx = 0
	E_branch_max = int(numpy.max(data[:,5]))

	# iterate energy branches
	for E_branch_idx in range(0, E_branch_max + 1):
		# filter data for given branch
		data_Q = [
			[ row[Q_idx1] for row in data if row[5] == E_branch_idx ],
			[ row[Q_idx2] for row in data if row[5] == E_branch_idx ]
		]
		data_E = [ row[3] for row in data if row[5] == E_branch_idx ]
		data_S = [ row[4] for row in data if row[5] == E_branch_idx ]

		if(len(data_E) < 1):
			continue

		# choose branch colour
		r = int(0xff - 0xff * (E_branch_idx / E_branch_max))
		b = int(0x00 + 0xff * (E_branch_idx / E_branch_max))

		axis.plot_trisurf(data_Q[0], data_Q[1], data_E,
			color = "#%02x00%02x%02x" % (r, b, alpha),
			antialiased = True)

	if xlabel != None:
		axis.set_xlabel(xlabel)
	else:
		axis.set_xlabel(labels[Q_idx1])
	if ylabel != None:
		axis.set_ylabel(ylabel)
	else:
		axis.set_ylabel(labels[Q_idx2])
	axis.set_zlabel("E (meV)")

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0)

	if plot_file != "":
		pyplot.savefig(plot_file)
	pyplot.show()
# -----------------------------------------------------------------------------


if __name__ == "__main__":
	data = numpy.array(calc_disp())

	if plot_dispersion:
		plot_disp_2d(data, Q_idx1, Q_idx2)
