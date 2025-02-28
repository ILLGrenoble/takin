#
# TODO: implementation of the violini algo
# (this is just a stub file for now)
#
# @author
# @date 2025
# @license see 'LICENSE' file
#
# @desc for algorithm: [vio14] N. Violini et al., NIM A 736 (2014) pp. 31-39, doi: 10.1016/j.nima.2013.10.042
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
import numpy.linalg as la
import reso
import helpers


#
# resolution algorithm
#
def calc(param, pointlike = False):
    # instrument position
    ki = param["ki"]
    kf = param["kf"]
    E = param["E"]
    Q = param["Q"]

    lam = helpers.k2lam(ki)

    # angles
    twotheta = helpers.get_scattering_angle(ki, kf, Q) * param["sample_sense"]
    thetas = twotheta * 0.5
    Q_ki = helpers.get_angle_Q_ki(ki, kf, Q) * param["sample_sense"]
    Q_kf = helpers.get_angle_Q_kf(ki, kf, Q) * param["sample_sense"]

    if param["verbose"]:
        print("2theta = %g deg, Q_ki = %g deg, Q_kf = %g deg.\n" %
            (twotheta*helpers.rad2deg, Q_ki*helpers.rad2deg, Q_kf*helpers.rad2deg))


    # TODO: resolution matrix
    R = np.eye(4)
    R0 = 1.


    # dict with results
    res = {}

    res["Q_avg"] = np.array([ Q, 0., 0., E ])
    res["ki"] = ki
    res["kf"] = kf
    res["Q_ki"] = Q_ki
    res["Q_kf"] = Q_kf
    res["twotheta"] = twotheta

    res["reso"] = R
    res["reso_v"] = np.array([0., 0., 0., 0.])
    res["reso_s"] = 0.
    res["r0"] = R0
    res["res_vol"] = reso.ellipsoid_volume(R)

    if np.isnan(res["r0"]) or np.isinf(res["r0"]) or np.isnan(res["reso"].any()) or np.isinf(res["reso"].any()):
        res["ok"] = False

    res["ok"] = True
    return res
