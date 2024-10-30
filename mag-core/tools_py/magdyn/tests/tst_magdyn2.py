#
# magdyn py interface test
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
import magdyn


# options
save_config_file       = False  # save magdyn file
save_dispersion        = False  # write dispersion to file
print_dispersion       = False  # write dispersion to console
plot_dispersion        = True   # show dispersion plot
only_positive_energies = True   # ignore magnon annihilation?

# weight scaling and clamp factors
S_scale = 64.
S_min   = 1.
S_max   = 500.

# dispersion plotting range
hkl_start = numpy.array([ 0., 0., 0. ])
hkl_end   = numpy.array([ 1., 0., 0. ])

# number of Qs to calculate on a dispersion direction
num_Q_points = 256


# create the magdyn object
mag = magdyn.MagDyn()


# add variables (either register them with magdyn or use them directly)
#magdyn.add_variable(mag, "J", -1.)
J = -1


# add magnetic sites
# the numeric values have to be passed as either all strings or all doubles
# (the strings can contain the variables registered above)
#magdyn.add_site(mag, "site_1",  "0", "0", "0",  "0", "0", "1",  "1")
magdyn.add_site(mag, "site_1",  0, 0, 0,  0, 0, 1,  1)


# add couplings between the sites
# the numeric values have to be passed as either all strings or all doubles
# (the strings can contain the variables registered above)
#magdyn.add_coupling(mag, "coupling_1",  "site_1", "site_1",  "1", "0", "0",  "J")
magdyn.add_coupling(mag, "coupling_1",  "site_1", "site_1",  1, 0, 0,  J)


# explicitly re-calculate all sites and couplings if needed
# (e.g., when changing the external field after adding the sites)
#magdyn.calc(mag)


# minimum energy
print("Energy minimum at Q=(000): {:.4f} meV".format(mag.CalcMinimumEnergy()))
print("Ground state energy: {:.4f} meV".format(mag.CalcGroundStateEnergy()))


if save_dispersion:
	# directly calculate a dispersion and write it to a file
	print("\nSaving dispersion to {}...".format("disp.dat"))
	mag.SaveDispersion("disp.dat",
		hkl_start[0], hkl_start[1], hkl_start[2],
		hkl_end[0], hkl_end[1], hkl_end[2],
		num_Q_points)


# manually calculate the same dispersion
print("\nManually calculating dispersion...")
if print_dispersion:
	print("{:>15} {:>15} {:>15} {:>15} {:>15}".format("h", "k", "l", "E", "S(Q,E)"))

data_h = []
data_k = []
data_l = []
data_E = []
data_S = []

for hkl in numpy.linspace(hkl_start, hkl_end, num_Q_points):
	for S in mag.CalcEnergies(hkl[0], hkl[1], hkl[2], False):
		if only_positive_energies and S.E < 0.:
			continue

		weight = S.weight * S_scale
		if weight < S_min:
			weight = S_min
		elif weight > S_max:
			weight = S_max

		data_h.append(hkl[0])
		data_k.append(hkl[1])
		data_l.append(hkl[2])
		data_E.append(S.E)
		data_S.append(weight)

		if print_dispersion:
			print("{:15.4f} {:15.4f} {:15.4f} {:15.4f} {:15.4g}".format(
				hkl[0], hkl[1], hkl[2], S.E, S.weight))


# plot the results
if plot_dispersion:
	print("Plotting dispersion...")

	import matplotlib.pyplot as plot

	fig = plot.figure()

	plt = fig.add_subplot(1, 1, 1)
	plt.set_xlabel("h (rlu)")
	plt.set_ylabel("E (meV)")
	plt.scatter(data_h, data_E, marker='.', s=data_S)

	plot.tight_layout()
	plot.show()
