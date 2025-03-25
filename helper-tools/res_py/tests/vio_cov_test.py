#
# Test for every function of cov.py
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

# To execute the test : pytest helper-tools/res_py/tests/vio_cov_test.py

import sys
import os
sys.path.append(os.path.dirname(__file__) + "/..")

import algos.vio_cov as vio_cov

import numpy as np
import pytest


def test_constant():
    assert vio_cov.m_nSI == 1.67492749804e-27
    assert vio_cov.hSI == 6.62607015e-34
    assert vio_cov.hbarSI == pytest.approx(1.054571817e-34, 0.000001e-34)
    assert vio_cov.mohSI == pytest.approx(15882536.11573, 0.000001)
    assert vio_cov.eV2J == 1.602176634e-19
    assert vio_cov.meV2J == pytest.approx(1.602176634e-22, 0.0000001e-22)
    assert vio_cov.m2A == 1e10
    assert vio_cov.mm2m == 1e-3
    assert vio_cov.me2kg == 9.1093826e-31
    assert vio_cov.mu2kg == 1.66053906892e-27
    assert vio_cov.a02m == 5.2917721092e-11
    assert vio_cov.a02A == 5.2917721092e-1
    assert vio_cov.hartree2J == 4.35974417e-18
    assert vio_cov.hartree2meV == pytest.approx(27211.382799, 0.00001)
    assert vio_cov.timeAU2s == pytest.approx(2.418884e-17, 0.00001e-17)
    assert vio_cov.vAU2vSI == pytest.approx(2187691.263, 0.001)
    assert vio_cov.hbarAU == 1
    assert vio_cov.meAU == 1
    assert vio_cov.m_nAU == pytest.approx(1838.68388, 0.00001)
    assert vio_cov.mohAU == pytest.approx(1838.68388, 0.00001)

def test_k2v():
    assert vio_cov.k2v(10) == pytest.approx(100000000000/vio_cov.mohSI, 0.00001)

def test_calcDistLin():
    assert vio_cov.calcDistLin(np.array([1,2,3,4,5,6,7,8])) == 16

def test_reducedCoef():
    assert vio_cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'SPHERE') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 2*np.sqrt(6), 2*np.sqrt(2), 4*np.sqrt(2)), 0.00001)
    assert vio_cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'HCYL') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 4*np.sqrt(2), 2*np.sqrt(2), 2*np.sqrt(6)), 0.00001)
    assert vio_cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'VCYL') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 2*np.sqrt(6), 2*np.sqrt(2), 4*np.sqrt(2)), 0.00001)

def test_jacobTerms():
    dQxSPHERE, dQySPHERE, dQzSPHERE, dESPHERE = vio_cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'SPHERE', 'AU')
    assert dQxSPHERE == pytest.approx(np.array([20*vio_cov.mohAU, -32*vio_cov.mohAU, -16*vio_cov.mohAU, -20*vio_cov.mohAU, 32*vio_cov.mohAU, -0.75*vio_cov.mohAU, -0.25*vio_cov.mohAU, 0.5*np.sqrt(2)*vio_cov.mohAU, 0.5*np.sqrt(6)*vio_cov.mohAU]), 0.00001)
    assert dQySPHERE == pytest.approx(np.array([25*vio_cov.mohAU, -40*vio_cov.mohAU, -20*vio_cov.mohAU, -25*vio_cov.mohAU, 40*vio_cov.mohAU, 0.25*np.sqrt(3)*vio_cov.mohAU, -0.25*np.sqrt(3)*vio_cov.mohAU, -0.5*np.sqrt(6)*vio_cov.mohAU, 0.5*np.sqrt(2)*vio_cov.mohAU]), 0.00001)
    assert dQzSPHERE == pytest.approx(np.array([30*vio_cov.mohAU, -48*vio_cov.mohAU, -24*vio_cov.mohAU, -30*vio_cov.mohAU, 48*vio_cov.mohAU, 0, 0.5*np.sqrt(3)*vio_cov.mohAU, 0, -np.sqrt(2)*vio_cov.mohAU]), 0.00001)
    assert dESPHERE == pytest.approx(np.array([500*vio_cov.m_nAU, -800*vio_cov.m_nAU, -400*vio_cov.m_nAU, -500*vio_cov.m_nAU, 800*vio_cov.m_nAU, 0, 0, 0, 0]), 0.00001)

    dQxHCYL, dQyHCYL, dQzHCYL, dEHCYL = vio_cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'HCYL', 'AU')
    assert dQxHCYL == pytest.approx(np.array([20*vio_cov.mohAU, -32*vio_cov.mohAU, -0.4*vio_cov.mohAU, -0, -20*vio_cov.mohAU, 32*vio_cov.mohAU, -0.75*vio_cov.mohAU, -0.25*vio_cov.mohAU, 0]), 0.00001)
    assert dQyHCYL == pytest.approx(np.array([25*vio_cov.mohAU, -40*vio_cov.mohAU, 0, -50*vio_cov.mohAU, -25*vio_cov.mohAU, 40*vio_cov.mohAU, 0.25*np.sqrt(3)*vio_cov.mohAU, -0.25*np.sqrt(3)*vio_cov.mohAU, -0.5*np.sqrt(6)*vio_cov.mohAU]), 0.00001)
    assert dQzHCYL == pytest.approx(np.array([30*vio_cov.mohAU, -48*vio_cov.mohAU, 0, -60*vio_cov.mohAU, -30*vio_cov.mohAU, 48*vio_cov.mohAU, 0, 0.5*np.sqrt(3)*vio_cov.mohAU, 0.5*np.sqrt(2)*vio_cov.mohAU]), 0.00001)
    assert dEHCYL == pytest.approx(np.array([500*vio_cov.m_nAU, -800*vio_cov.m_nAU, -200*np.sqrt(2)*vio_cov.m_nAU, -200*np.sqrt(2)*vio_cov.m_nAU, -500*vio_cov.m_nAU, 800*vio_cov.m_nAU, 0, 0, 0]), 0.00001)

    dQxVCYL, dQyVCYL, dQzVCYL, dEVCYL = vio_cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'VCYL', 'AU')
    assert dQxVCYL == pytest.approx(np.array([20*vio_cov.mohAU, -32*vio_cov.mohAU, -40*vio_cov.mohAU, -0, -20*vio_cov.mohAU, 32*vio_cov.mohAU, -0.75*vio_cov.mohAU, -0.25*vio_cov.mohAU, 0.5*np.sqrt(2)*vio_cov.mohAU]), 0.00001)
    assert dQyVCYL == pytest.approx(np.array([25*vio_cov.mohAU, -40*vio_cov.mohAU, -50*vio_cov.mohAU, 0, -25*vio_cov.mohAU, 40*vio_cov.mohAU, 0.25*np.sqrt(3)*vio_cov.mohAU, -0.25*np.sqrt(3)*vio_cov.mohAU, -0.5*np.sqrt(6)*vio_cov.mohAU]), 0.00001)
    assert dQzVCYL == pytest.approx(np.array([30*vio_cov.mohAU, -48*vio_cov.mohAU, 0, -0.4*vio_cov.mohAU, -30*vio_cov.mohAU, 48*vio_cov.mohAU, 0, 0.5*np.sqrt(3)*vio_cov.mohAU, 0]), 0.00001)
    assert dEVCYL == pytest.approx(np.array([500*vio_cov.m_nAU, -800*vio_cov.m_nAU, -200*np.sqrt(2)*vio_cov.m_nAU, -200*np.sqrt(2)*vio_cov.m_nAU, -500*vio_cov.m_nAU, 800*vio_cov.m_nAU, 0, 0, 0]), 0.00001)

def test_getParamChopper():
    win1, rot1 = vio_cov.getParamChopper(np.array([1,2,3,4]))
    win2, rot2 = vio_cov.getParamChopper(np.array([5,6,7,-1]))
    assert (win1, rot1) == (1, 4)
    assert (win2, rot2) == (5, 6)

def test_listDeltaGeo():
    l_test = vio_cov.listDeltaGeo(np.array([10, 1, 20, 2, 30, 3, 40, 4]))
    assert (l_test[0], l_test[1], l_test[2], l_test[3]) == (1, 2, 3, 4)

def test_deltaTimeChopper():
    assert vio_cov.deltaTimeChopper(60, 0.5, 'AU') == 20/vio_cov.timeAU2s
    assert vio_cov.deltaTimeChopper(18, 1, 'AU') == 3/vio_cov.timeAU2s

def test_listDeltaTime():
    l_test = vio_cov.listDeltaTime(24, 1, 18, 1, 4/vio_cov.timeAU2s, 'AU')
    assert (l_test[0], l_test[1]) == (5/vio_cov.timeAU2s, 5/vio_cov.timeAU2s)

def test_getDeltas():
    pG = {'dist_PM_u':[10, 1], 'dist_MS_u':[20, 2], 'dist_SD_u':[30, 3], 'angles_rad':[1, 0.1, 2, 0.2, 3, 0.3, 4, 0.4], 'delta_time_detector_u':4/vio_cov.timeAU2s}
    pC = {'chopperP':[24, 1, 1, 1], 'chopperM':[18, 1, 1, 1]}
    s = [1, 1, 1, 2, 4, 9]
    l_test = vio_cov.getDeltas(pG, pC, s, 'AU')
    assert (l_test[0], l_test[1], l_test[2], l_test[3], l_test[4], l_test[5], l_test[6], l_test[7], l_test[8]) == (1, 2, 3, 5/vio_cov.timeAU2s, 5/vio_cov.timeAU2s, 0.1, 0.2, 0.3, 0.4)

def test_jacobianMatrix():
    assert np.array_equal( vio_cov.jacobianMatrix(np.array([0, 1, 2, 3]), np.array([4, 5, 6, 7]), np.array([8, 9, 10, 11]), np.array([12, 13, 14, 15])), np.array(([0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15])) )

def test_covxiMatrix():
    assert np.array_equal( vio_cov.covxiMatrix([1, 2, 3, 4, 5]),  np.array(([1, 0, 0, 0, 0], [0, 4, 0, 0, 0], [0, 0, 9, 0, 0], [0, 0, 0, 16, 0], [0, 0, 0, 0, 25])) )

def test_fromSItoAU():
    assert vio_cov.fromSItoAU(1, 'time') == 1/vio_cov.timeAU2s
    assert vio_cov.fromSItoAU(1, 'distance') == 1/vio_cov.a02m
    assert vio_cov.fromSItoAU(1, 'velocity') == 1/vio_cov.vAU2vSI
