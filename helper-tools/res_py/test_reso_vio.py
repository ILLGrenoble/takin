#
# tests the resolution calculation
#
# @author Mecoli Victor <mecoli.ill.fr>
# @date feb-2025
# @license GPLv2
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

# requires numpy version >= 1.10
import numpy as np
import helpers
import tas
import reso
import vio
import params_in5


np.set_printoptions(floatmode = "fixed",  precision = 4)

reso_method = "vio"    # "vio", "eck", "pop", or "cn"
verbose = True
params = params_in5.params

# calculate resolution
res = vio.calc(params)
if not res["ok"]:
    print("RESOLUTION CALCULATION FAILED!")
    exit(-1)

if(verbose and reso_method != "vio"):
    print("R0 = %g, Vol = %g" % (res["r0"], res["res_vol"]))
    print("Resolution matrix:\n%s" % res["reso"])
    print("Resolution vector: %s" % res["reso_v"])
    print("Resolution scalar: %g" % res["reso_s"])


# describe and plot ellipses
ellipses = reso.calc_ellipses(res["reso"], verbose = verbose)
reso.plot_ellipses(ellipses, verbose = verbose)
