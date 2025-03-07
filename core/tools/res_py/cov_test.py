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

# To execute the test : pytest core/tools/res_py/cov_test.py

import cov
import numpy as np
import pytest

def test_constant():
    assert cov.m_nSI == 1.67492749804e-27
    assert cov.hSI == 6.62607015e-34
    assert cov.hbarSI == pytest.approx(1.054571817e-34, 0.000001e-34)
    assert cov.mohSI == pytest.approx(15882536.11573, 0.000001)
    assert cov.eV2J == 1.602176634e-19
    assert cov.m2A == 1e10
    assert cov.mm2m == 1e-3
    assert cov.me2kg == 9.1093826e-31
    assert cov.mu2kg == 1.66053906892e-27
    assert cov.a02m == 5.2917721092e-11
    assert cov.a02A == 5.2917721092e-1
    assert cov.hartree2J == 4.35974417e-18
    assert cov.hartree2meV == pytest.approx(27211.382799, 0.00001)
    assert cov.timeAU2s == pytest.approx(2.418884e-17, 0.00001e-17)
    assert cov.vAU2vSI == pytest.approx(2187691.263, 0.001)
    assert cov.hbarAU == 1
    assert cov.meAU == 1
    assert cov.m_nAU == pytest.approx(1838.68388, 0.00001)
    assert cov.mohAU == pytest.approx(1838.68388, 0.00001)

def test_k2v():
    assert cov.k2v(10) == pytest.approx(100000000000/cov.mohSI, 0.00001)

def test_calcDistLin():
    assert cov.calcDistLin(np.array([1,2,3,4,5,6,7,8])) == 16

def test_reducedCoef():
    assert cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'SPHERE') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 2*np.sqrt(6), 2*np.sqrt(2), 4*np.sqrt(2)), 0.00001)
    assert cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'HCYL') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 4*np.sqrt(2), 2*np.sqrt(2), 2*np.sqrt(6)), 0.00001)
    assert cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'VCYL') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 2*np.sqrt(6), 2*np.sqrt(2), 4*np.sqrt(2)), 0.00001)

def test_jacobTerms():
    dQxSPHERE, dQySPHERE, dQzSPHERE, dESPHERE = cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'SPHERE', 'AU')
    assert dQxSPHERE == pytest.approx(np.array([20*cov.mohAU, -32*cov.mohAU, -16*cov.mohAU, -20*cov.mohAU, 32*cov.mohAU, -0.75*cov.mohAU, -0.25*cov.mohAU, 0.5*np.sqrt(2)*cov.mohAU, 0.5*np.sqrt(6)*cov.mohAU]), 0.00001)
    assert dQySPHERE == pytest.approx(np.array([25*cov.mohAU, -40*cov.mohAU, -20*cov.mohAU, -25*cov.mohAU, 40*cov.mohAU, 0.25*np.sqrt(3)*cov.mohAU, -0.25*np.sqrt(3)*cov.mohAU, -0.5*np.sqrt(6)*cov.mohAU, 0.5*np.sqrt(2)*cov.mohAU]), 0.00001)
    assert dQzSPHERE == pytest.approx(np.array([30*cov.mohAU, -48*cov.mohAU, -24*cov.mohAU, -30*cov.mohAU, 48*cov.mohAU, 0, 0.5*np.sqrt(3)*cov.mohAU, 0, -np.sqrt(2)*cov.mohAU]), 0.00001)
    assert dESPHERE == pytest.approx(np.array([500*cov.m_nAU, -800*cov.m_nAU, -400*cov.m_nAU, -500*cov.m_nAU, 800*cov.m_nAU, 0, 0, 0, 0]), 0.00001)

    dQxHCYL, dQyHCYL, dQzHCYL, dEHCYL = cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'HCYL', 'AU')
    assert dQxHCYL == pytest.approx(np.array([20*cov.mohAU, -32*cov.mohAU, -0.4*cov.mohAU, -0, -20*cov.mohAU, 32*cov.mohAU, -0.75*cov.mohAU, -0.25*cov.mohAU, 0]), 0.00001)
    assert dQyHCYL == pytest.approx(np.array([25*cov.mohAU, -40*cov.mohAU, 0, -50*cov.mohAU, -25*cov.mohAU, 40*cov.mohAU, 0.25*np.sqrt(3)*cov.mohAU, -0.25*np.sqrt(3)*cov.mohAU, -0.5*np.sqrt(6)*cov.mohAU]), 0.00001)
    assert dQzHCYL == pytest.approx(np.array([30*cov.mohAU, -48*cov.mohAU, 0, -60*cov.mohAU, -30*cov.mohAU, 48*cov.mohAU, 0, 0.5*np.sqrt(3)*cov.mohAU, 0.5*np.sqrt(2)*cov.mohAU]), 0.00001)
    assert dEHCYL == pytest.approx(np.array([500*cov.m_nAU, -800*cov.m_nAU, -200*np.sqrt(2)*cov.m_nAU, -200*np.sqrt(2)*cov.m_nAU, -500*cov.m_nAU, 800*cov.m_nAU, 0, 0, 0]), 0.00001)
    
    dQxVCYL, dQyVCYL, dQzVCYL, dEVCYL = cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'VCYL', 'AU')
    assert dQxVCYL == pytest.approx(np.array([20*cov.mohAU, -32*cov.mohAU, -40*cov.mohAU, -0, -20*cov.mohAU, 32*cov.mohAU, -0.75*cov.mohAU, -0.25*cov.mohAU, 0.5*np.sqrt(2)*cov.mohAU]), 0.00001)
    assert dQyVCYL == pytest.approx(np.array([25*cov.mohAU, -40*cov.mohAU, -50*cov.mohAU, 0, -25*cov.mohAU, 40*cov.mohAU, 0.25*np.sqrt(3)*cov.mohAU, -0.25*np.sqrt(3)*cov.mohAU, -0.5*np.sqrt(6)*cov.mohAU]), 0.00001)
    assert dQzVCYL == pytest.approx(np.array([30*cov.mohAU, -48*cov.mohAU, 0, -0.4*cov.mohAU, -30*cov.mohAU, 48*cov.mohAU, 0, 0.5*np.sqrt(3)*cov.mohAU, 0]), 0.00001)
    assert dEVCYL == pytest.approx(np.array([500*cov.m_nAU, -800*cov.m_nAU, -200*np.sqrt(2)*cov.m_nAU, -200*np.sqrt(2)*cov.m_nAU, -500*cov.m_nAU, 800*cov.m_nAU, 0, 0, 0]), 0.00001)

def test_getParamChopper():
    win1, rot1 = cov.getParamChopper(np.array([1,2,3,4]))
    win2, rot2 = cov.getParamChopper(np.array([5,6,7,-1]))
    assert (win1, rot1) == (1, 4)
    assert (win2, rot2) == (5, 6)

def test_listDeltaGeo():
    l_test = cov.listDeltaGeo(np.array([10, 1, 20, 2, 30, 3, 40, 4]))
    assert (l_test[0], l_test[1], l_test[2], l_test[3]) == (1, 2, 3, 4)

def test_deltaTimeChopper():
    assert cov.deltaTimeChopper(60, 0.5, 'AU') == 20/cov.timeAU2s
    assert cov.deltaTimeChopper(18, 1, 'AU') == 3/cov.timeAU2s

def test_listDeltaTime():
    l_test = cov.listDeltaTime(24, 1, 18, 1, 4/cov.timeAU2s, 'AU')
    assert (l_test[0], l_test[1]) == (5/cov.timeAU2s, 5/cov.timeAU2s)

def test_getDeltas():
    pG = {'dist_PM_u':[10, 1], 'dist_MS_u':[20, 2], 'dist_SD_u':[30, 3], 'angles_rad':[1, 0.1, 2, 0.2, 3, 0.3, 4, 0.4], 'delta_time_detector_u':4/cov.timeAU2s}
    pC = {'chopperP':[24, 1, 1, 1], 'chopperM':[18, 1, 1, 1]}
    s = [1, 1, 1, 2, 4, 9]
    l_test = cov.getDeltas(pG, pC, s, 'AU')
    assert (l_test[0], l_test[1], l_test[2], l_test[3], l_test[4], l_test[5], l_test[6], l_test[7], l_test[8]) == (1, 2, 3, 5/cov.timeAU2s, 5/cov.timeAU2s, 0.1, 0.2, 0.3, 0.4)

def test_jacobianMatrix():
    assert np.array_equal( cov.jacobianMatrix(np.array([0, 1, 2, 3]), np.array([4, 5, 6, 7]), np.array([8, 9, 10, 11]), np.array([12, 13, 14, 15])), np.array(([0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15])) )

def test_covxiMatrix():
    assert np.array_equal( cov.covxiMatrix([1, 2, 3, 4, 5]),  np.array(([1, 0, 0, 0, 0], [0, 4, 0, 0, 0], [0, 0, 9, 0, 0], [0, 0, 0, 16, 0], [0, 0, 0, 0, 25])) )

def test_fromSItoAU():
    assert cov.fromSItoAU(1, 'time') == 1/cov.timeAU2s
    assert cov.fromSItoAU(1, 'distance') == 1/cov.a02m
    assert cov.fromSItoAU(1, 'velocity') == 1/cov.vAU2vSI