#
# TODO: implementation of the violini algo
# 
#
# @author Mecoli Victor <mecoli@ill.fr>
# @date feb-2025
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
import cr_json as js


#
# resolution algorithm
#


# settings
ki = 1  # 2.5 in A
kf = 1
Q = 1

verbose = True


params = {
    "ki":ki,
    "kf":kf,
    "Q":Q,
    "E":helpers.get_E(ki, kf),

    "v_i":cov.k2v(ki),
    "v_f":cov.k2v(kf),

    "verbose" : verbose,

    "det_shape": 'SPHERE', #'VCYL', #'HCYL', #'SPHERE'

    "z":0
}

p_geo = {'dist_PM':[25000, 0], 'dist_MS':[1300, 6], 'dist_SD':[3000, 6], 'angles':[0, 0.2, 0, 0.2, 10, 0.2, 90, 0.2], 'delta_time_detector':0.006}
p_chop = {'chopperP':[1, 1, 1, 4.386], 'chopperM':[1, 1, 1, 27.778]}

instr = {
    "L_PM":25000,
    "L_MS":1300,
    "rad":3000,
    "delta_PM":0,
    "delta_MS":6,
    "delta_SD":6,
    "theta_i":0,
    "phi_i":0,
    "delta_theta_i":0.2,
    "delta_phi_i":0.2,
    "rad_chop_P":1,
    "RPM_min_chop_P":1,
    "RPM_max_chop_P":1000,
    "RPM_chop_P":4.386,
    "rad_chop_M":1,
    "RPM_min_chop_M":1,
    "RPM_max_chop_M":1000,
    "RPM_chop_M":27.778,
}

def calc(paramUser, paramInstr):
    #paramInstr are read from a file (not yet implemented)

    # Storage of informations given by users
    # Normes of ki, kf and Q
    ki = paramUser["ki"]
    kf = paramUser["kf"]
    Q = paramUser["Q"]
    # Angular velocity of choppers
    rot_speedP = paramUser["rot_speedP"]
    rot_speedM = paramUser["rot_speedM"]
    # Shape of the detector
    shape = paramUser["shape"]
    # Velocity
    vi = cov.k2v(ki)
    vf = cov.k2v(kf)
    # Angles
    Q_ki = helpers.get_angle_Q_ki(ki, kf, Q)
    print()

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
#    Q = np.array([0, 0, 0])
    covQhw = np.eye(4)
    if(shape == "Sph"):
        covQhw = cov.covSph(dict_geo, dict_choppers, vi, vf)
#        Q = cov.vecQSph(vi, vf, paramInstr["angles"][0], paramInstr["angles"][2], paramInstr["angles"][4], paramInstr["angles"][6])
    elif(shape == "Hcyl"):
        covQhw = cov.covHcyl(dict_geo, dict_choppers, vi, vf)
#        cos_phi_f = np.divide(paramInstr["dist_SD"][2], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
#        sin_phi_f = np.divide(paramInstr["dist_SD"][0], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
#        Q = cov.vecQHcyl(vi, vf, paramInstr["angles"][0], paramInstr["angles"][2], paramInstr["angles"][4], cos_phi_f, sin_phi_f)
    elif(shape == "Vcyl"):
        covQhw = cov.covVcyl(dict_geo, dict_choppers, vi, vf)
#        cos_phi_f = np.divide(paramInstr["dist_SD"][0], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
#        sin_phi_f = np.divide(paramInstr["dist_SD"][2], np.sqrt(np.square(paramInstr["dist_SD"][0]) + np.square(paramInstr["dist_SD"][2])))
#        Q = cov.vecQVcyl(vi, vf, paramInstr["angles"][0], paramInstr["angles"][2], paramInstr["angles"][4], cos_phi_f, sin_phi_f)
    covQhwInv = la.inv(covQhw)
    if(paramUser["verbose"]):
        print("E =", E) #, "; Q =", Q)
        print("covQhw =", covQhw)
        print("covQhwInv =", covQhwInv)
    # Going from ki, kf, Qz to Qpara, Qperp, Qz :
    rot = helpers.rotation_matrix_4d_zE(-Q_ki)
    covQhwInv = np.dot(rot.T, np.dot(covQhwInv, rot))
    if(paramUser["verbose"]):
        print("In the base (Qpara, Qperp, Qz) :")
        print("rot =", rot,'\ncovQhwInv =', covQhwInv)
    res={}
    res["ki"] = ki
    res["kf"] = kf
    res["Q"] = Q
    res["E"] = E
    res["reso"] = covQhwInv
    return res

#param_geo is a dictionary: {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[ (if HCYL: x, sigma_x), radius, sigma_r, (if VCYL: z, sigma_z)], 
#                            angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f, (if SPHERE: phi_f, sigma_phi_f)], delta_time_detector:value (0 by default)},
#param_choppers is a dictionary: {chopperP:[window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[window_angle, min_rot_speed, max_rot_speed, rot_speed]}, dist in mm, angles in degree, rot_speed in RPM,
#v_i, v_f: velocity of the incident and scattered neutron m/s,
#shape = SPHERE, VCYL, HCYL: shape of the detector (sphere, vertical cylinder or horizontal cylinder)"""
p_geo = {'dist_PM':[25000, 0], 'dist_MS':[1300, 6], 'dist_SD':[3000, 6], 'angles':[0, 0.2, 0, 0.2, 10, 0.2, 90, 0.2], 'delta_time_detector':0.006}
p_chop = {'chopperP':[1, 1, 1, 4.386], 'chopperM':[1, 1, 1, 27.778]}
covQhw = cov.cov(p_geo, p_chop, 620, 620, 'SPHERE', True)
aff = la.inv(covQhw)
print('aff =', aff)
ellipses = reso.calc_ellipses(aff, True)
reso.plot_ellipses(ellipses, True)
print(620*cov.mohSI/cov.m2A)

#distances in mm, angle in degree
#IN5 = {"pos_P1":9000, "delta_p1":40, "pos_P2":9000, "delta_P2":40, "pos_M1":1000, "delta_M1":10, "pos_M2":1000, "delta_M2":10, "pos_S":0, "delat_S":0,
#       "det_R":4000, "delta_R":12.7, "det_h":3000, "delat_h":10, "det_theta":146.73, "delta_theta":0.19} #h_min = -1500, h_max = 1500, theta_min, theta_max
#js.create(IN5, '/home/mecoli/Git/IN5')
#param = js.read('/home/mecoli/Git/IN5.json')

# pU = {"ki":16000000000, "kf":20000000000, "Q":7000000000, "rot_speedP":9000, "rot_speedM":8000, "shape":"Vcyl", "verbose": True}
# pI = {"dist_PM":[20, 0.001, 3, 0.002], "dist_MS":[3, 0.005], "dist_SD":[4, 0.01, 0, 0.008], "angles":[0, 0.2, 0, 0.2, 0, 0.15], "chopperP":[0.8, 10, 7000, 1700], "chopperM":[0.8, 8, 7000, 1700]}
# aff = calc(pU, pI)

# describe and plot ellipses
# ellipses = reso.calc_ellipses(aff["reso"], True)
# reso.plot_ellipses(ellipses, True)

#print(param)

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
