#
# magdyn py interface demo -- creating and plotting a model
# @author Tobias Weber <tweber@ill.fr>
# @date march-2024
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
save_config_file       = False  # save magdyn file
save_dispersion        = False  # write dispersion to file
print_dispersion       = False  # write dispersion to console
plot_dispersion        = True   # show dispersion plot
only_positive_energies = True   # ignore magnon annihilation?
max_threads            = 4      # number of worker threads, 0: automatic determination
num_Q_points           = 256    # number of Qs to calculate on a dispersion direction
S_scale                = 64.    # weight scaling and clamp factors
S_clamp_min            = 1.     #
S_clamp_max            = 500.   #
S_filter_min           = -1.    # don't filter
dispfile               = "disp.dat"

# dispersion plotting range
hkl_start              = numpy.array([ 0., 0., 0.5 ])
hkl_end                = numpy.array([ 1., 1., 0.5 ])
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Create the magnetic model
# -----------------------------------------------------------------------------
# The given magnetic model and its parameters are from this paper:
#     https://doi.org/10.1103/PhysRevB.101.144411
#     (which is also available here: https://arxiv.org/abs/2002.06283).
#

# create the magdyn object
mag = magdyn.MagDyn()

#
# add variables
#
magdyn.add_variable(mag, "J1", -0.58)
magdyn.add_variable(mag, "J2", -0.93)

magdyn.add_variable(mag, "D1x", +0.15)
magdyn.add_variable(mag, "D1y", -0.24)
magdyn.add_variable(mag, "D1z", -0.05)

magdyn.add_variable(mag, "D2x", +0.16)
magdyn.add_variable(mag, "D2y", +0.1)
magdyn.add_variable(mag, "D2z", +0.36)


#
# set the sample environment
#
magdyn.set_temperature(mag, 50)
magdyn.set_field(mag,  0, 0, 1,  1,  False)


#
# add magnetic sites
#                      name     x    y    z  sx sy sz   S
magdyn.add_site(mag, "Cu 1",  0.5,   0, 0.5,  0, 0, 1,  1)
magdyn.add_site(mag, "Cu 2",    0, 0.5, 0.5,  0, 0, 1,  1)
magdyn.add_site(mag, "Cu 3",  0.5, 0.5,   0,  0, 0, 1,  1)
magdyn.add_site(mag, "Cu 4",    0,   0,   0,  0, 0, 1,  1)


#
# add couplings between sites
#                                name   site 1  site 2    dx   dy   dz      J    dmix   dmiy   dmiz
magdyn.add_coupling(mag, "Coupling A",  "Cu 2", "Cu 1",  "0", "0", "0",  "J1",  "D1x", "D1y", "D1z")
magdyn.add_coupling(mag, "Coupling B",  "Cu 4", "Cu 1",  "0", "0", "0",  "J2",  "D2x", "D2y", "D2z")


# calculate the rest of the couplings by symmetry
print("Calculating symmetry-equivalent couplings...")
magdyn.symmetrise_couplings(mag, "P 21 3")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate sites and couplings
# -----------------------------------------------------------------------------
print("Calculating sites and couplings...")
magdyn.calc(mag)


# minimum energy
print("\nEnergy minimum at Q = (000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))


if save_config_file:
	# save the magdyn configuration file for use in the gui program or the convolution fitter
	print("\nSaving configuration to {}...".format("model.magdyn"))
	if not mag.Save("model.magdyn"):
		print("Error saving configuration.")


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
