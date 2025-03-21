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
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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
# helpers
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

#--------------------------------------------------------------------------
