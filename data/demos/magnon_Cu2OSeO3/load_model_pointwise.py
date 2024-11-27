#
# magdyn py interface demo -- loading the model from a file
#                             and query the dispersion point by point
# @author Tobias Weber <tweber@ill.fr>
# @date 12-oct-2023
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
from takin import magdyn


# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
save_dispersion        = False  # write dispersion to file
print_dispersion       = False  # write dispersion to console
plot_dispersion        = True   # show dispersion plot
only_positive_energies = True   # ignore magnon annihilation?
use_threadpool         = True   # parallelise calculation
max_threads            = 4      # number of worker threads
num_Q_points           = 256    # number of Qs to calculate on a dispersion direction
S_scale                = 64.    # weight scaling and clamp factors
S_clamp_min            = 1.     #
S_clamp_max            = 500.   #
S_filter_min           = -1.    # don't filter
modelfile              = "model.magdyn"
dispfile               = "disp.dat"

# dispersion plotting range
hkl_start              = numpy.array([ 0., 0., 0.5 ])
hkl_end                = numpy.array([ 1., 1., 0.5 ])
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


if save_dispersion:
	# directly calculate a dispersion and write it to a file
	print("\nSaving dispersion to {}...".format(dispfile))
	mag.SaveDispersion(dispfile,
		hkl_start[0], hkl_start[1], hkl_start[2],
		hkl_end[0], hkl_end[1], hkl_end[2],
		num_Q_points)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# manually calculate the dispersion
# -----------------------------------------------------------------------------
print("\nManually calculating dispersion...")
if print_dispersion:
	print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

data_h = []
data_k = []
data_l = []
data_E = []
data_S = []


#
# add a data point
#
def append_data(h, k, l, E, S):
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


if use_threadpool:
	import concurrent.futures as fut
	executor = fut.ThreadPoolExecutor # fut.ProcessPoolExecutor

	print(f"Using {max_threads} threads.")


	#
	# calculate the energies and weights for a Q point
	#
	def calc_Es(h, k, l):
		Es_dicts = mag.CalcEnergies(h, k, l, False)

		Es = []
		for Es_dict in Es_dicts:
			Es.append(( h, k, l, Es_dict.E, Es_dict.weight ))

		return Es


	with executor(max_workers = max_threads) as exe:
		Es_futures = []

		# submit tasks
		for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
			Es_futures.append(exe.submit(calc_Es, hkl[0], hkl[1], hkl[2]))

		# get results from tasks
		for Es_future in Es_futures:
			for (h, k, l, E, weight) in Es_future.result():
				if only_positive_energies and E < 0.:
					continue

				append_data(h, k, l, E, weight)

				if print_dispersion:
					print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
						h, k, l, E, weight))

else:  # no thread pool
	for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
		for S in mag.CalcEnergies(hkl[0], hkl[1], hkl[2], False):
			if only_positive_energies and S.E < 0.:
				continue

			append_data(hkl[0], hkl[1], hkl[2], S.E, S.weight)

			if print_dispersion:
				print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
					hkl[0], hkl[1], hkl[2], S.E, S.weight))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the results
# -----------------------------------------------------------------------------
if plot_dispersion:
	import matplotlib.pyplot as plot
	print("Plotting dispersion...")

	fig = plot.figure()

	plt = fig.add_subplot(1, 1, 1)
	plt.set_xlabel("h (rlu)")
	plt.set_ylabel("E (meV)")
	plt.scatter(data_h, data_E, marker='.', s=data_S)

	plot.tight_layout()
	plot.show()
# -----------------------------------------------------------------------------
