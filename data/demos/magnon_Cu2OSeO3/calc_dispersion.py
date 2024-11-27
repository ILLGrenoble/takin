#
# magdyn py interface demo -- plotting multiple dispersion branches
# @author Tobias Weber <tweber@ill.fr>
# @date 18-sep-2024
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

import time
import numpy
import numpy.linalg
from takin import magdyn


# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
print_dispersion       = False  # write dispersion to console
only_positive_energies = True   # ignore magnon annihilation?
max_threads            = 0      # number of worker threads, 0: automatic determination
num_Q_points           = 256    # number of Qs on a dispersion branch
show_dividers          = False  # show vertical bars between dispersion branches
use_custom_labels      = True   # show custom dispersion labels (branch_labels)
use_colours            = True   # use dispersion colours (branch_colours)
S_scale                = 64.    # weight scaling and clamp factors
S_clamp_min            = 1.     #
S_clamp_max            = 500.   #
S_filter_min           = -1.    # don't filter
modelfile              = "model.magdyn"
plotfile               = ""     # file to save plot to

# cubic high-symmetry points
pt_G  = numpy.array([ 0.0, 0.0, 0.0 ])
pt_X1 = numpy.array([ 0.5, 0.0, 0.0 ])
pt_X2 = numpy.array([ 0.0, 0.5, 0.0 ])
pt_X3 = numpy.array([ 0.0, 0.0, 0.5 ])
pt_M1 = numpy.array([ 0.5, 0.5, 0.0 ])
pt_M2 = numpy.array([ 0.5, 0.0, 0.5 ])
pt_M3 = numpy.array([ 0.0, 0.5, 0.5 ])
pt_R  = numpy.array([ 0.5, 0.5, 0.5 ])

# dispersion branches to plot
dispersion = [ pt_G, pt_X1, pt_M1, pt_R ]

branch_labels = [ "Γ", "X", "M", "R" ]
branch_colours = [ "#ffffff", "#eeeeee", "#ffffff" ]

width_ratios = []                   # lengths from one dispersion point to the next
num_branches = len(dispersion) - 1  # number of dispersion braches
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# load the magnetic model
# -----------------------------------------------------------------------------
mag = None


#
# calculate the energies and weights for a Q point
#
def calc_Es(h, k, l):
	Es_dicts = mag.CalcEnergies(h, k, l, False)

	Es = []
	for Es_dict in Es_dicts:
		Es.append(( h, k, l, Es_dict.E, Es_dict.weight ))

	return Es


def setup_struct():
	# create the magdyn object:
	global mag
	mag = magdyn.MagDyn()

	# load the model file
	print("Loading {}...".format(modelfile))
	if not mag.Load(modelfile):
		print("Failed loading {}.".format(modelfile))
		exit(-1)

	# minimum energy
	print("\nEnergy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
	print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate the dispersion branches
# -----------------------------------------------------------------------------
def calc_disp():
	# calculate the dispersion branches
	print("\nCalculating %d dispersion branches..." % num_branches)
	if print_dispersion:
		print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

	dispersion_plot_indices = []  # Q components to use for plotting branches
	data = []                     # data for all branches

	for branch_idx in range(num_branches):
		hkl_start = dispersion[branch_idx]
		hkl_end = dispersion[branch_idx + 1]

		width_ratios.append(numpy.linalg.norm(hkl_end - hkl_start))

		print("[%d/%d] Calculating %s  ->  %s branch..." %
		   (branch_idx + 1, num_branches, hkl_start, hkl_end))

		# find scan axis
		Q_diff = [
			numpy.abs(hkl_start[0] - hkl_end[0]),
			numpy.abs(hkl_start[1] - hkl_end[1]),
			numpy.abs(hkl_start[2] - hkl_end[2]) ]

		axis_idx = 0
		if Q_diff[1] > Q_diff[axis_idx]:
			axis_idx = 1
		elif Q_diff[2] > Q_diff[axis_idx]:
			axis_idx = 2
		dispersion_plot_indices.append(axis_idx)


		data_h = []
		data_k = []
		data_l = []
		data_E = []
		data_S = []

		#
		# add a data point
		#
		def append_data(h, k, l, E, S):
			#print(f"S(Q = ({h} {k} {l}) rlu, E = {E} meV) = {S}")
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


		# calculate the dispersion
		data_disp = mag.CalcDispersion(
			hkl_start[0], hkl_start[1], hkl_start[2],
			hkl_end[0], hkl_end[1], hkl_end[2],
			num_Q_points, max_threads)
		for S in data_disp:
			for data_EandS in S.E_and_S:
				if only_positive_energies and data_EandS.E < 0.:
					continue

				append_data(magdyn.get_h(S), magdyn.get_k(S), magdyn.get_l(S),
					data_EandS.E, data_EandS.weight)

				if print_dispersion:
					print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
						magdyn.get_h(S), magdyn.get_k(S), magdyn.get_l(S),
						data_EandS.E, data_EandS.weight))

		data.append([ branch_idx, data_h, data_k, data_l, data_E, data_S ])

	return (data, dispersion_plot_indices)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the results
# -----------------------------------------------------------------------------
def plot_disp(data, dispersion_plot_indices):
	import matplotlib.pyplot as plot
	plot.rcParams.update({
		"font.sans-serif" : "DejaVu Sans",
		"font.family" : "sans-serif",
		"font.size" : 12,
	})

	print("\nPlotting dispersion branches...")

	(plt, axes) = plot.subplots(nrows = 1, ncols = num_branches,
		width_ratios = width_ratios, sharey = True)

	# in case there's only one sub-plot
	if type(axes) != numpy.ndarray:
		axes = [ axes ]

	for ( branch_idx, data_h, data_k, data_l, data_E, data_S ) in data:
		# dispersion branch start and end points
		b1 = dispersion[branch_idx]
		b2 = dispersion[branch_idx + 1]

		plot_idx = dispersion_plot_indices[branch_idx]
		if plot_idx == 0:
			data_x = data_h
		elif plot_idx == 1:
			data_x = data_k
		elif plot_idx == 2:
			data_x = data_l

		# ticks and labels
		axes[branch_idx].set_xlim(data_x[0], data_x[-1])

		if use_colours:
			axes[branch_idx].set_facecolor(branch_colours[branch_idx])

		if use_custom_labels:
			tick_labels = [
				branch_labels[branch_idx],
				branch_labels[branch_idx + 1] ]
		else:
			tick_labels = [
				"(%.4g %.4g %.4g)" % (b1[0], b1[1], b1[2]),
				"(%.4g %.4g %.4g)" % (b2[0], b2[1], b2[2]) ]

		if branch_idx == 0:
			axes[branch_idx].set_ylabel("E (meV)")
		else:
			axes[branch_idx].get_yaxis().set_visible(False)
			if not show_dividers:
				axes[branch_idx].spines["left"].set_visible(False)

			tick_labels[0] = ""

		if not show_dividers and branch_idx != num_branches - 1:
			axes[branch_idx].spines["right"].set_visible(False)

		axes[branch_idx].set_xticks([data_x[0], data_x[-1]], labels = tick_labels)

		if branch_idx == num_branches / 2 - 1:
			axes[branch_idx].set_xlabel("Q (rlu)")

		# plot the dispersion branch
		axes[branch_idx].scatter(data_x, data_E, marker = '.', s = data_S)

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0)

	if plotfile != "":
		plot.savefig(plotfile)
	plot.show()
# -----------------------------------------------------------------------------


if __name__ == "__main__":
	setup_struct()

	tick = time.time()
	(data, dispersion_plot_indices) = calc_disp()
	time_needed = time.time() - tick
	print("Calculation took %.4g s." % time_needed)

	plot_disp(data, dispersion_plot_indices)
