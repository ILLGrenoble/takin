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

import libs.tas as tas
import libs.helpers as helpers

import algos.vio_cov as vio_cov

import numpy as np
import numpy.linalg as la


#
# resolution algorithm
#

# Qup = Qz



def calc(param):
    # Storage of informations
    verbose = param["verbose"]
    # Normes of ki, kf and Q
    ki = param["ki"]
    kf = param["kf"]
    Q = param["Q"]
    ###########################################################################################
    if(verbose):
        print("ki =", ki, "; kf =", kf, "; Q =", Q)

    # Shape of the detector
    det_shape = param["det_shape"]
    ###########################################################################################
    if(verbose):
        print("det_shape =", det_shape)

    # Velocity in m/s
    vi = vio_cov.k2v(ki)
    vf = vio_cov.k2v(kf)
    ###########################################################################################
    if(verbose):
        print("vi =", vi, "; vf =", vf)

    # Angles
    theta_i = param["angles"][0]
    phi_i = param["angles"][2]
    theta_f = param["angles"][4]
    phi_f = 0
    if(det_shape == 'SPHERE'):
        phi_f = param["angles"][6]
    if(det_shape == 'HCYL'):
        phi_f = np.rad2deg( np.arctan( np.divide(param["dist_SD"][0], param["dist_SD"][2]) ) )
    if(det_shape == 'VCYL'):
        phi_f = np.rad2deg( np.arctan( np.divide(param["dist_SD"][2], param["dist_SD"][0]) ) )
    ###########################################################################################
    if(verbose):
        print("theta_i =", theta_i, "; phi_i =", phi_i,"; theta_f =", theta_f, "; phi_f =", phi_f)

    ki_xy = ki*np.cos(np.deg2rad(phi_i))
    ki_z = ki*np.sin(np.deg2rad(phi_i))
    kf_xy = kf*np.cos(np.deg2rad(phi_f))
    kf_z = kf*np.sin(np.deg2rad(phi_f))
    Q_x = ki_xy*np.cos(np.deg2rad(theta_i)) - kf_xy*np.cos(np.deg2rad(theta_f))
    Q_y = ki_xy*np.sin(np.deg2rad(theta_i)) - kf_xy*np.sin(np.deg2rad(theta_f))
    Q_z = ki_z - kf_z
    Q_xy = np.sqrt( np.square(Q_x) + np.square(Q_y) )
    ###########################################################################################
    if(verbose):
        print("ki_xy =", ki_xy, "; ki_z =", ki_z,"; kf_xy =", kf_xy, "; kf_z =", kf_z, "; Q_x =", Q_x, "; Q_y =", Q_y, "; Q_xy =", Q_xy, "; Q_z =", Q_z)

    # Information on the instrument
    dict_geo = {"dist_PM":param["dist_PM"], "dist_MS":param["dist_MS"], "dist_SD":param["dist_SD"], "angles":param["angles"], "delta_time_detector":param["delta_time_det"]}
    dict_choppers = {"chopperP":param["chopperP"], "chopperM":param["chopperM"]}
    ###########################################################################################
    if(verbose):
        print("\ndict_geo =", dict_geo, "\ndict_choppers =", dict_choppers, "\n")

    # Energy transfer, Q vector and Covariance matrix
    E = tas.get_E(ki, kf)
    vec_Q = np.array([Q_x, Q_y, Q_z])
    covQhw = vio_cov.cov(dict_geo, dict_choppers, vi, vf, det_shape, 'SI', verbose)
    covQhwInv = la.inv(covQhw)
    ###########################################################################################
    if(param["verbose"]):
        print("E =", E, "; vec_Q =", vec_Q)
        print("covQhw =", covQhw)
        print("covQhwInv =", covQhwInv)

    # Going from ki, kf, Qz to Qpara, Qperp, Qz :
    Q_ki = tas.get_psi(ki_xy, kf_xy, Q_xy, param["sample_sense"])
    rot = helpers.rotation_matrix_nd(-Q_ki, 4)
    covQhwInv = np.dot(rot.T, np.dot(covQhwInv, rot))
    ###########################################################################################
    if(verbose):
        print("In the base (Qpara, Qperp, Qz) :")
        print("rot =", rot,'\ncovQhwInv =', covQhwInv)

    res={}
    res["ki"] = ki
    res["kf"] = kf
    res["Q"] = Q
    res["vec_Q"] = vec_Q
    res["E"] = E
    res["reso"] = covQhwInv
    res["ok"] = True
    return res
