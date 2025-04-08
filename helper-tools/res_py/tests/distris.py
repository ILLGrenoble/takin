#
# compares random distributions
#
# @author Tobias Weber <tweber@ill.fr>
# @date 8-apr-2025
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

import numpy as np


#
# normal distribution
# see: https://en.wikipedia.org/wiki/Normal_distribution
#
def gauss(x, x0, sig):
	norm = (np.sqrt(2.*np.pi) * sig)
	return np.exp(-0.5*((x - x0)/sig)**2.) / norm


# see: https://en.wikipedia.org/wiki/Continuous_uniform_distribution
# see: https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
len = 1.
sig_unif = len / np.sqrt(12.)
sig_circ = len / np.sqrt(16.)
print("Length of uniform distribution: %.4f, sigma = %.4f." % (len, sig_unif))
print("  -> Replacement gaussian HWHM: %.4f." % (sig_unif*np.sqrt(2.*np.log(2.))))
print("Length of circular distribution: %.4f, sigma = %.4f." % (len, sig_circ))
print("  -> Replacement gaussian HWHM: %.4f." % (sig_circ*np.sqrt(2.*np.log(2.))))


# plotting
import matplotlib.pyplot as plt
import matplotlib.patches as pat

xs = np.linspace(-sig_unif*5., sig_unif*5., 128)
plt.plot(xs, gauss(xs, 0., sig_unif), label = "gaussian replacement for uniform distribution")
#plt.plot(xs, gauss(xs, 0., sig_circ), label = "gaussian replacement for circular distribution")
plt.gca().add_patch(pat.Rectangle((-0.5*len, 0.), len, 1./len, fill = False))
#plt.gca().add_patch(pat.Ellipse((0., 0.), len, 2.5/len, fill = False))

plt.ylim([0., gauss(0., 0., sig_circ)])
plt.legend()
plt.show()
