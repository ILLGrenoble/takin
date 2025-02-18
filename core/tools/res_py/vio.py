#
# TODO: implementation of the violini algo
# (this is just a stub file for now)
#
# @author 
# @date 2025
# @license GPLv2
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
import cov


#
# resolution algorithm
#

def calc(paramUser, paramInstr):
    #paramInstr are read from a file (not yet implemented)

    # Storage of informations given by users
    vi = paramUser["vi"]
    vf = paramUser["vf"]
    rot_speedP = paramUser["rot_speedP"]
    rot_speedM = paramUser["rot_speedM"]
    shape = paramUser["shape"]

    # Information on the instrument
    dict_geo = {"dist_PM":paramInstr["dist_PM"], "dist_MS":paramInstr["dist_MS"], "dist_SD":paramInstr["dist_SD"], "angles":paramInstr["angles"]}
    param_chopP = np.zeros(5)
    param_chopP[:4] = paramInstr["chopperP"]
    param_chopP[4] = rot_speedP
    param_chopM = np.zeros(5)
    param_chopM[:4] = paramInstr["chopperM"]
    param_chopM[4] = rot_speedM
    dict_choppers = {"chopperP":param_chopP, "chopperM":param_chopM}

    # Energy, Q vector and Covariance matrix
    E = cov.energy(vi, vf)
    Q = np.array([0, 0, 0])
    covQhw = np.eye(4)
    if(shape == "Sph"):
        covQhw = cov.covSph(dict_geo, dict_choppers, vi, vf)
        Q = cov.vecQSph(vi, vf, paramInstr["angles"][0], paramInstr["angles"][2], paramInstr["angles"][4], paramInstr["angles"][6])
    elif(shape == "Hcyl"):
        covQhw = cov.covHcyl(dict_geo, dict_choppers, vi, vf)
        cos_phi_f = np.divide(paramInstr["dist_SD"][2], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
        sin_phi_f = np.divide(paramInstr["dist_SD"][0], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
        Q = cov.vecQHcyl(vi, vf, paramInstr["angles"][0], paramInstr["angles"][2], paramInstr["angles"][4], cos_phi_f, sin_phi_f)
    elif(shape == "Vcyl"):
        covQhw = cov.covVcyl(dict_geo, dict_choppers, vi, vf)
        cos_phi_f = np.divide(paramInstr["dist_SD"][0], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
        sin_phi_f = np.divide(paramInstr["dist_SD"][2], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
        Q = cov.vecQVcyl(vi, vf, paramInstr["angles"][0], paramInstr["angles"][2], paramInstr["angles"][4], cos_phi_f, sin_phi_f)
    covQhwInv = la.inv(covQhw)
    if(paramUser["verbose"]):
        print("E =", E, "; Q =", Q)
        print("covQhw =", covQhw)
        print("covQhwInv =", covQhwInv)


pU = {"vi":2000, "vf":1000, "rot_speedP":9000, "rot_speedM":8000, "shape":"Vcyl", "verbose": True}
pI = {"dist_PM":[20, 0.001, 3, 0.002], "dist_MS":[3, 0.005], "dist_SD":[4, 0.01, 0, 0.008], "angles":[0, 0.2, 0, 0.2, 0, 0.15], "chopperP":[0.8, 10, 7000, 1700], "chopperM":[0.8, 8, 7000, 1700]}
calc(pU, pI)


# def calc(param, pointlike = False):
#     # instrument position
#     ki = param["ki"]
#     kf = param["kf"]
#     E = param["E"]
#     Q = param["Q"]

#     lam = helpers.k2lam(ki)

#     # angles
#     twotheta = helpers.get_scattering_angle(ki, kf, Q) * param["sample_sense"]
#     thetas = twotheta * 0.5
#     Q_ki = helpers.get_angle_Q_ki(ki, kf, Q) * param["sample_sense"]
#     Q_kf = helpers.get_angle_Q_kf(ki, kf, Q) * param["sample_sense"]

#     if param["verbose"]:
#         print("2theta = %g deg, Q_ki = %g deg, Q_kf = %g deg.\n" %
#             (twotheta*helpers.rad2deg, Q_ki*helpers.rad2deg, Q_kf*helpers.rad2deg))


#     # TODO: resolution matrix
#     R = np.eye(4)
#     R0 = 1.


#     # dict with results
#     res = {}

#     res["Q_avg"] = np.array([ Q, 0., 0., E ])
#     res["ki"] = ki
#     res["kf"] = kf
#     res["Q_ki"] = Q_ki
#     res["Q_kf"] = Q_kf
#     res["twotheta"] = twotheta

#     res["reso"] = R
#     res["reso_v"] = np.array([0., 0., 0., 0.])
#     res["reso_s"] = 0.
#     res["r0"] = R0
#     res["res_vol"] = reso.ellipsoid_volume(R)

#     if np.isnan(res["r0"]) or np.isinf(res["r0"]) or np.isnan(res["reso"].any()) or np.isinf(res["reso"].any()):
#         res["ok"] = False

#     res["ok"] = True
#     return res
