#
# Covariance matrix cov(Q,\hbar\omega) calculations for different shape of the detector
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

# Constant
m_n = 1.674927351e-27 #mass of the neutron in kg
h = 6.62607015e-34 #Plank constante in J.s
hbar = np.divide(h, 2*np.pi)
moh = np.divide(m_n, hbar)
covQhw = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

# Shape of the detector: vertical cylinder. The incident beam arrive align with the radius of the cylinder
def covVcyl(param_geo, param_choppers, v_i, v_f, verbose=False):
    """param_geo is a dictionary: {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[radius, sigma_r, heigh, sigma_h], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f]},
    param_choppers is a dictionary: {chopperP:[diameter, window_angle, min_rot_speed, max_rot_speed], chopperM:[diameter, window_angle, min_rot_speed, max_rot_speed]}, dist in mm, angle in degree, rot_speed in RPM,
    v_i, v_f: velocity of the incident and scattered neutron"""

    # Calcul of distances
    L_PM = 0
    L_MS = 0
    nb_PM = len(param_geo['dist_PM'])
    nb_MS = len(param_geo['dist_MS'])
    for i in range(0,nb_PM,2):
        L_PM += param_geo['dist_PM'][i]
    for i in range(0,nb_MS,2):
        L_MS += param_geo['dist_MS'][i]
    rad = param_geo['dist_SD'][0]
    hei = param_geo['dist_SD'][2]
    L_SD = np.sqrt(np.square(rad) + np.square(hei))
    if(verbose):
        print('param_geo =', param_geo)
        print('param_choppers =', param_choppers)
        print('v_i =', v_i, 'v_f =', v_f)
        print('L_PM =', L_PM, 'L_MS =', L_MS, 'L_SD =', L_SD)

    # Definition of variables
    theta_i = param_geo['angles'][0]
    phi_i = param_geo['angles'][2]
    theta_f = param_geo['angles'][4]
    cos_phi_f = np.divide(rad, L_SD)
    sin_phi_f = np.divide(hei, L_SD)
    Ae = np.divide(np.power(v_i,3), L_PM)
    Ax = np.divide(np.square(v_i), L_PM)*np.cos(theta_i)*np.cos(phi_i)
    Ay = np.divide(np.square(v_i), L_PM)*np.sin(theta_i)*np.cos(phi_i)
    Az = np.divide(np.square(v_i), L_PM)*np.sin(phi_i)
    Be = np.divide(np.power(v_f,3), L_PM)
    Bx = np.divide(np.square(v_i), L_PM)*np.cos(theta_f)*cos_phi_f
    By = np.divide(np.square(v_i), L_PM)*np.sin(theta_f)*cos_phi_f
    Bz = np.divide(np.square(v_i), L_PM)*sin_phi_f
    if(verbose):
        print('Ae =', Ae, 'Ax =', Ax, 'Ay =', Ay, 'Az =', Az)
        print('Be =', Be, 'Bx =', Bx, 'By =', By, 'Bz =', Bz)
    
    # Creation of the Jacobian and  the parameter's uncertainty matrices
    size_PM = np.divide(nb_PM, 2)
    size_MS = np.divide(nb_MS, 2)
    size_SD = np.divide(len(param_geo['dist_SD']), 2)
    size_tps = 2
    size_angles = np.divide(len(param_geo['angles']), 2)
    nb_col = int(size_PM + size_MS + size_SD + size_tps + size_angles)
    jacobian = np.zeros((4, nb_col))
    covxi = np.eye(nb_col)
    # Calcul of Jacobian's terms
    # x terms
    dQxdL_PMj = np.divide(moh, v_i)*(Ax + Bx*np.divide(L_MS, L_SD))
    dQxdL_MSj = -np.divide(moh, v_i)*Bx*np.divide(L_PM,L_SD)
    dQxdrad = -np.divide(moh, v_f)*Bx*np.divide(L_PM,rad)
    dQxdt_PM = -moh*(Ax + Bx*np.divide(L_MS, L_SD))
    dQxdt_MD = moh*Bx*np.divide(L_PM, L_SD)
    dQxdtheta_i = -moh*v_i*np.sin(theta_i)*np.cos(phi_i)
    dQxdphi_i = -moh*v_i*np.cos(theta_i)*np.sin(phi_i)
    dQxdtheta_f = moh*v_f*np.sin(theta_f)*cos_phi_f
    if(verbose):
        print('dQx/dL_PMj =', dQxdL_PMj, 'dQx/dL_MSj =', dQxdL_MSj, 'dQx/drad =', dQxdrad)
        print('dQx/dt_PM =', dQxdt_PM, 'dQx/dt_MD =', dQxdt_MD)
        print('dQx/dtheta_i =', dQxdtheta_i, 'dQx/dphi_i =', dQxdphi_i, 'dQx/dtheta_f =', dQxdtheta_f)
    # y terms
    dQydL_PMj = np.divide(moh, v_i)*(Ay + By*np.divide(L_MS, L_SD))
    dQydL_MSj = -np.divide(moh, v_i)*By*np.divide(L_PM, L_SD)
    dQydrad = -np.divide(moh, v_f)*By*np.divide(L_PM, rad)
    dQydt_PM = -moh*(Ay + By*np.divide(L_MS, L_SD))
    dQydt_MD = moh*By*np.divide(L_PM, L_SD)
    dQydtheta_i = moh*v_i*np.cos(theta_i)*np.cos(phi_i)
    dQydphi_i = -moh*v_i*np.sin(theta_i)*np.sin(phi_i)
    dQydtheta_f = -moh*v_f*np.cos(theta_f)*cos_phi_f
    if(verbose):
        print('dQy/dL_PMj =', dQydL_PMj, 'dQy/dL_MSj =', dQydL_MSj, 'dQy/drad =', dQydrad)
        print('dQy/dt_PM =', dQydt_PM, 'dQy/dt_MD =', dQydt_MD)
        print('dQy/dtheta_i =', dQydtheta_i, 'dQy/dphi_i =', dQydphi_i, 'dQy/dtheta_f =', dQydtheta_f)
    # z terms
    dQzdL_PMj = np.divide(moh, v_i)*(Az + Bz*np.divide(L_MS, L_SD))
    dQzdL_MSj = -np.divide(moh, v_i)*Bz*np.divide(L_PM, L_SD)
    dQzdz = -moh*np.divide(v_f, L_SD)
    dQzdt_PM = -moh*(Az + Bz*np.divide(L_MS, L_SD))
    dQzdt_MD = moh*Bz*np.divide(L_PM, L_SD)
    dQzdphi_i = moh*v_i*np.cos(phi_i)
    if(verbose):
        print('dQz/dL_PMj =', dQzdL_PMj, 'dQz/dL_MSj =', dQzdL_MSj, 'dQz/dz =', dQzdz)
        print('dQz/dt_PM =', dQzdt_PM, 'dQz/dt_MD =', dQzdt_MD)
        print('dQz/dphi_i =', dQzdphi_i)
    # energy terms
    dEdL_PMj = np.divide(m_n, v_i)*(Ae + Be*np.divide(L_MS,L_SD))
    dEdL_MSj = -np.divide(m_n, v_i)*Be*np.divide(L_PM,L_SD)
    dEdL_SDj = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)
    dEdt_PM = -m_n*(Ae + Be*np.divide(L_MS,L_SD))
    dEdt_MD = m_n*Be*np.divide(L_PM,L_SD)
    if(verbose):
        print('dE/dL_PMj =', dEdL_PMj, 'dE/dL_MSj =', dEdL_MSj, 'dE/dL_SDj =', dEdL_SDj)
        print('dE/dt_PM =', dEdt_PM, 'dE/dt_MD =', dEdt_MD)



#test perso

a = {"P":[1,2], "U":[3,4]}
print(len(a['P']))
print(np.square(5))
print(np.power(3,4))
print(np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]))
b = {'u':12, 'x':'f'}
print('b =', b['u'])
print(np.zeros((4,7)))
print(np.eye(7))