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

# To execute the test : pytest PATH/cov_test.py

import cov
import numpy as np
import pytest


def test_k2v():
    assert int(cov.k2v(10)) == 6296

def test_calcDistLin():
    assert cov.calcDistLin(np.array([1,2,3,4,5,6,7,8])) == 16

def test_reducedCoef():
    assert cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'SPHERE') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 2*np.sqrt(6), 2*np.sqrt(2), 4*np.sqrt(2)), 0.00001)
    assert cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'HCYL') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 4*np.sqrt(2), 2*np.sqrt(2), 2*np.sqrt(6)), 0.00001)
    assert cov.reducedCoef(3, 4, 2, np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), 'VCYL') == pytest.approx((13.5, 1.125*np.sqrt(3), 3.375, 2.25, 32, 2*np.sqrt(6), 2*np.sqrt(2), 4*np.sqrt(2)), 0.00001)

def test_jacobTerms():
    dQxSPHERE, dQySPHERE, dQzSPHERE, dESPHERE = cov.jacobTerms(1, 2, np.array([20, 10, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'SPHERE')
    assert dQxSPHERE == pytest.approx(np.array([20*cov.mohAU, -32*cov.mohAU, -16*cov.mohAU, -20*cov.mohAU, 32*cov.mohAU, -0.75*cov.mohAU, -0.25*cov.mohAU, 0.5*np.sqrt(2)*cov.mohAU, 0.5*np.sqrt(6)*cov.mohAU]), 0.00001)
    # dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f, dQxdphi_f]) ))
    assert dQySPHERE == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f, dQydphi_f]) ))
    assert dQzSPHERE == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f, dQzdphi_f]) ))
    assert dESPHERE == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdrad, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f, dEdphi_f]) ))

    dQxHCYL, dQyHCYL, dQzHCYL, dEHCYL = cov.jacobTerms(1, 2, np.array([10, 20, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'HCYL')
    assert dQxHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdx, dQxdrad, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))
    assert dQyHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydx, dQydrad, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))
    assert dQzHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdx, dQzdrad, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))
    assert dEHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdx, dEdrad, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f]) ))
    
    dQxVCYL, dQyVCYL, dQzVCYL, dEVCYL = cov.jacobTerms(1, 2, np.array([10, 20, 5, 2, 2, 2]), np.array([np.pi/3, np.pi/6, np.pi/6, np.pi/4]), np.array([100, 4, 5, 6, 200, 8, 10, 12]), np.array([1, 1, 1, 2, 4, 9]), 'VCYL')
    assert dQxHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdz, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))
    assert dQyHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydz, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))
    assert dQzHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdz, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))
    assert dEHCYL == pytest.approx(np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]), 0.00001)
    # dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdrad, dEdz, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f]) ))