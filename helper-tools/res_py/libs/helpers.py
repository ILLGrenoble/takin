#
# helper functions
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license see 'LICENSE' file
#
# @desc for reso algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984)
# @desc for vertical scattering modification: [eck20] G. Eckold, personal communication, 2020.
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

import numpy as np
import numpy.linalg as la



#--------------------------------------------------------------------------
# constants
#--------------------------------------------------------------------------
sig2fwhm = 2.*np.sqrt(2.*np.log(2.))
cm2A = 1e8
min2rad = 1./ 60. / 180.*np.pi
rad2deg = 180. / np.pi
deg2rad = np.pi / 180.
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
# general helper functions
#--------------------------------------------------------------------------
#
# x rotation matrix
#
def rotation_matrix_3d_x(angle):
    s = np.sin(angle)
    c = np.cos(angle)

    return np.array([
        [ 1, 0,  0 ],
        [ 0, c, -s ],
        [ 0, s,  c ]])


def rotation_matrix_2d(angle):
    s = np.sin(angle)
    c = np.cos(angle)

    return np.array([
        [ c, -s ],
        [ s,  c ]])


def rotation_matrix_nd(angle, dims = 3):
    R = np.eye(dims)
    R[0:2, 0:2] = rotation_matrix_2d(angle)
    return R


def mirror_matrix(iSize, iComp):
    mat = np.identity(iSize)
    mat[iComp, iComp] = -1.

    return mat;


#
# thin lens equation: 1/f = 1/lenB + 1/lenA
#
def focal_len(lenBefore, lenAfter):
    f_inv = 1./lenBefore + 1./lenAfter
    return 1. / f_inv


#
# optimal mono/ana curvature,
# see e.g.
#    - (Shirane 2002) p. 66
#    - or nicos/nicos-core.git/tree/nicos/devices/tas/mono.py in nicos
#    - or Monochromator_curved.comp in McStas
#
def foc_curv(lenBefore, lenAfter, tt, vert):
    f = focal_len(lenBefore, lenAfter)
    s = np.abs(np.sin(0.5*tt))

    if vert:
        curv = 2. * f*s
    else:
        curv = 2. * f/s

    return curv


#
# optimal mono/ana curvature, using wavenumber k and crystal d
# see e.g.
#    - (Shirane 2002) p. 66
#    - or nicos/nicos-core.git/tree/nicos/devices/tas/mono.py in nicos
#    - or Monochromator_curved.comp in McStas
#
def foc_curv_2(f, k, d, vert):
    #f = focal_len(lenBefore, lenAfter)
    s = np.abs(np.pi / (d * k))

    if vert:
        curv = 2. * f*s
    else:
        curv = 2. * f/s

    return curv


#
# adjugate matrix
# see e.g.: https://en.wikipedia.org/wiki/Adjugate_matrix
#
def adjugate(mat):
    rows = len(mat)
    cols = len(mat[0])

    adj = np.zeros((rows, cols))

    for i in range(0, rows):
        for j in range(0, cols):
            submat = np.delete(np.delete(mat, i, axis = 0), j, axis = 1)
            adj[i, j] = (-1.)**(i + j) * la.det(submat)

    return np.transpose(adj)


#
# orthonormalise a given system
# see e.g.
#    - https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
#    - (Arens 2015), p. 744
#
def orthonormalise(mat):
    order = len(mat)
    mat_new = mat

    for idx in range(order):
        for idx2 in range(idx):
            # remove projections onto row vector with index 2
            len_sq = np.dot(mat_new[idx2, :], mat_new[idx2, :])
            proj = np.dot(mat[idx, :], mat_new[idx2, :]) / len_sq
            mat_new[idx, :] -= np.dot(proj, mat_new[idx2, :])

        # normalise
        mat_new[idx, :] /= la.norm(mat_new[idx, :])

    return mat_new

#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
# scattering triangle calculations
#--------------------------------------------------------------------------
def calc_triangle(param, eps = 1e-4):
    ki = param["ki"]
    kf = param["kf"]
    Q = param["Q"]

    # is Q given as (hkl) vector in rlu?
    if isinstance(Q, type([])) or isinstance(Q, type(np.array([]))):
        Q_rlu = param["Q"]
        orient_rlu = param["sample_plane_1"]
        B = tas.get_B(param["sample_lattice"], param["sample_angles"])
        orient_up_rlu = tas.cross(orient_rlu, param["sample_plane_2"], B)

        [a3, a4, dist_Q_plane] = tas.get_sample_angle(ki, kf, Q_rlu, orient_rlu, orient_up_rlu, B)
        param["theta"] = a3
        Q = param["Q"] = tas.get_Q(ki, kf, a4)

        if np.abs(dist_Q_plane) > eps:
            raise Warning("Q is not in the scattering plane.")

        param["sample_orient"] = orthonormalise(
            np.array([
                np.dot(B, orient_rlu),
                np.dot(B, param["sample_plane_2"]),
                np.dot(B, orient_up_rlu)
            ]))

        if param["verbose"]:
            print("Q = %g / A, theta = %g deg." % (Q, a3*rad2deg))
            print("Scattering plane orthonormal system [1/A]:\n%s" % param["sample_orient"])

    # is Q given as a scalar in 1/A?
    else:
        param["theta"] = 0.
        param["sample_orient"] = np.eye(3)

    # angles
    param["twotheta"] = tas.get_scattering_angle(ki, kf, Q) * param["sample_sense"]
    param["thetam"] = tas.get_mono_angle(ki, param["mono_xtal_d"], True) * param["mono_sense"]
    param["thetaa"] = tas.get_mono_angle(kf, param["ana_xtal_d"], True) * param["ana_sense"]
    param["Q_ki"] = tas.get_psi(ki, kf, Q, param["sample_sense"])
    param["Q_kf"] = tas.get_eta(ki, kf, Q, param["sample_sense"])

    if param["verbose"]:
        print("2theta = %g deg, thetam = %g deg, thetaa = %g deg, Q_ki = %g deg, Q_kf = %g deg.\n" % (
            param["twotheta"] * rad2deg,
            param["thetam"] * rad2deg, param["thetaa"] * rad2deg,
            param["Q_ki"] * rad2deg, param["Q_kf"] * rad2deg))
#--------------------------------------------------------------------------
