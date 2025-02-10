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
print(moh)

#Initialisation of the covariance Matrix
covQhw = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

# Calcul of the length of a segment cuts in small segments
def calcDistLin(dist):
    """dist in an array: [L1, delta1, L2, delta2, ...], L: distance, delta: uncerntainty"""
    L = 0
    nb = len(dist)
    for i in range(0,nb,2):
        L += dist[i]
    return L

# Storage of the given uncerntainty in a list
def listDeltaGeo(list_param):
    """list_param in an array: [P1, delta1, P2, delta2, ...], P: parameter, delta: uncerntainty"""
    dlt = []
    nb = len(list_param)
    for i in range(1,nb,2):
        dlt.append(list_param[i])
    return np.array(dlt)

# Calcul of the chopper's time uncertainty
def deltaTimeChopper(window_angle, rot_speed):
    """angle in degree, rot_speed in RPM"""
    return( np.divide(window_angle, 6*rot_speed) ) # result in seconds

# Return the time uncerntainty for t_PM and t_MD
def listDeltaTime(window_angleP, rot_speedP, window_angleM, rot_speedM):
    """angles in degree, rot_speed in RPM"""
    dltP = deltaTimeChopper(window_angleP, rot_speedP)
    dltM = deltaTimeChopper(window_angleM, rot_speedM)
    return np.array([np.sqrt(np.square(dltP) + np.square(dltM)), dltM])

# Creation and filling of the covxi matrix (uncerntainty on independant parameters xi)
def covxiMatrix(delta):
    """delta is a list of uncerntainties (for independant variables)"""
    nb_dlt = len(delta)
    covxi = np.eye(nb_dlt)
    for i in range(nb_dlt):
        covxi[i][i] = np.square(delta[i])
    return covxi

# Creation and filling of the Jacobian matrix
def jacobianMatrix(dQx, dQy, dQz, dE):
    nb_col = len(dQx)
    jacob = np.zeros((4, nb_col))
    jacob[0] = np.copy(dQx)
    jacob[1] = np.copy(dQy)
    jacob[2] = np.copy(dQz)
    jacob[3] = np.copy(dE)
    return jacob

# Shape of the detector: vertical cylinder. The incident beam arrive align with the radius of the cylinder
def covVcyl(param_geo, param_choppers, v_i, v_f, verbose=False):
    """param_geo is a dictionary: {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[radius, sigma_r, heigh, sigma_h], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f]},
    param_choppers is a dictionary: {chopperP:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed]}, dist in m, angle in degree, rot_speed in RPM,
    v_i, v_f: velocity of the incident and scattered neutron (m/s)"""

    # Calcul of distances
    rad = param_geo['dist_SD'][0]
    hei = param_geo['dist_SD'][2]
    L_PM = calcDistLin(param_geo['dist_PM'])
    L_MS = calcDistLin(param_geo['dist_MS'])
    L_SD = np.sqrt(np.square(rad) + np.square(hei))
    ###########################################################################################
    if(verbose):
        print('param_geo =', param_geo)
        print('param_choppers =', param_choppers)
        print('v_i =', v_i, '; v_f =', v_f)
        print('L_PM =', L_PM, '; L_MS =', L_MS, '; L_SD =', L_SD)

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
    ###########################################################################################
    if(verbose):
        print('Ae =', Ae, '; Ax =', Ax, '; Ay =', Ay, '; Az =', Az)
        print('Be =', Be, '; Bx =', Bx, '; By =', By, '; Bz =', Bz)
    
    # Number of variables for each set (distances, times, angles)
    size_PM = int(np.divide(len(param_geo['dist_PM']), 2))
    size_MS = int(np.divide(len(param_geo['dist_MS']), 2))
    size_SD = int(np.divide(len(param_geo['dist_SD']), 2))
    size_tps = 2
    size_angles = int(np.divide(len(param_geo['angles']), 2))
    nb_param = size_PM + size_MS + size_SD + size_tps + size_angles

    # Calcul of Jacobian's terms
    # x terms
    dQxdL_PMj = np.divide(moh, v_i)*(Ax + Bx*np.divide(L_MS, L_SD))
    dQxdL_MSj = -np.divide(moh, v_i)*Bx*np.divide(L_PM,L_SD)
    dQxdrad = -np.divide(moh, v_f)*Bx*np.divide(L_PM,rad)
    dQxdz = 0
    dQxdt_PM = -moh*(Ax + Bx*np.divide(L_MS, L_SD))
    dQxdt_MD = moh*Bx*np.divide(L_PM, L_SD)
    dQxdtheta_i = -moh*v_i*np.sin(theta_i)*np.cos(phi_i)
    dQxdphi_i = -moh*v_i*np.cos(theta_i)*np.sin(phi_i)
    dQxdtheta_f = moh*v_f*np.sin(theta_f)*cos_phi_f
    dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdz, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))
    ###########################################################################################
    if(verbose):
        print('x terms:')
        print('dQx/dL_PMj =', dQxdL_PMj, '; dQx/dL_MSj =', dQxdL_MSj, '; dQx/drad =', dQxdrad)
        print('dQx/dt_PM =', dQxdt_PM, '; dQx/dt_MD =', dQxdt_MD)
        print('dQx/dtheta_i =', dQxdtheta_i, '; dQx/dphi_i =', dQxdphi_i, '; dQx/dtheta_f =', dQxdtheta_f)
    # y terms
    dQydL_PMj = np.divide(moh, v_i)*(Ay + By*np.divide(L_MS, L_SD))
    dQydL_MSj = -np.divide(moh, v_i)*By*np.divide(L_PM, L_SD)
    dQydrad = -np.divide(moh, v_f)*By*np.divide(L_PM, rad)
    dQydz = 0
    dQydt_PM = -moh*(Ay + By*np.divide(L_MS, L_SD))
    dQydt_MD = moh*By*np.divide(L_PM, L_SD)
    dQydtheta_i = moh*v_i*np.cos(theta_i)*np.cos(phi_i)
    dQydphi_i = -moh*v_i*np.sin(theta_i)*np.sin(phi_i)
    dQydtheta_f = -moh*v_f*np.cos(theta_f)*cos_phi_f
    dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydz, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))
    ###########################################################################################
    if(verbose):
        print('y terms:')
        print('dQy/dL_PMj =', dQydL_PMj, '; dQy/dL_MSj =', dQydL_MSj, '; dQy/drad =', dQydrad)
        print('dQy/dt_PM =', dQydt_PM, '; dQy/dt_MD =', dQydt_MD)
        print('dQy/dtheta_i =', dQydtheta_i, '; dQy/dphi_i =', dQydphi_i, '; dQy/dtheta_f =', dQydtheta_f)
    # z terms
    dQzdL_PMj = np.divide(moh, v_i)*(Az + Bz*np.divide(L_MS, L_SD))
    dQzdL_MSj = -np.divide(moh, v_i)*Bz*np.divide(L_PM, L_SD)
    dQzdrad = 0
    dQzdz = -moh*np.divide(v_f, L_SD)
    dQzdt_PM = -moh*(Az + Bz*np.divide(L_MS, L_SD))
    dQzdt_MD = moh*Bz*np.divide(L_PM, L_SD)
    dQzdtheta_i = 0
    dQzdphi_i = moh*v_i*np.cos(phi_i)
    dQzdtheta_f = 0
    dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdz, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))
    ###########################################################################################
    if(verbose):
        print('z terms:')
        print('dQz/dL_PMj =', dQzdL_PMj, '; dQz/dL_MSj =', dQzdL_MSj, '; dQz/dz =', dQzdz)
        print('dQz/dt_PM =', dQzdt_PM, '; dQz/dt_MD =', dQzdt_MD)
        print('dQz/dphi_i =', dQzdphi_i)
    # energy terms
    dEdL_PMj = np.divide(m_n, v_i)*(Ae + Be*np.divide(L_MS,L_SD))
    dEdL_MSj = -np.divide(m_n, v_i)*Be*np.divide(L_PM,L_SD)
    dEdrad = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)*cos_phi_f
    dEdz = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)*sin_phi_f
    dEdt_PM = -m_n*(Ae + Be*np.divide(L_MS,L_SD))
    dEdt_MD = m_n*Be*np.divide(L_PM,L_SD)
    dEdtheta_i = 0
    dEdphi_i = 0
    dEdtheta_f = 0
    dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdrad, dEdz, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f]) ))
    ###########################################################################################
    if(verbose):
        print('E  terms:')
        print('dE/dL_PMj =', dEdL_PMj, '; dE/dL_MSj =', dEdL_MSj, '; dE/drad =', dEdrad, '; dE/dz =', dEdz)
        print('dE/dt_PM =', dEdt_PM, '; dE/dt_MD =', dEdt_MD)
    
    # List of uncertainty
    set1 = size_PM
    set2 = size_PM + size_MS
    set3 = size_PM + size_MS + size_SD
    set4 = size_PM + size_MS + size_SD + size_tps
    win_angleP = param_choppers['chopperP'][1]
    rot_speedP = param_choppers['chopperP'][4]
    if(rot_speedP == -1):
        rot_speedP = param_choppers['chopperP'][2]
    win_angleM = param_choppers['chopperM'][1]
    rot_speedM = param_choppers['chopperM'][4]
    if(rot_speedM == -1):
        rot_speedM = param_choppers['chopperM'][2]
    ###########################################################################################
    if(verbose):
        print('win_angleP =', win_angleP, '; rot_speedP = ', rot_speedP)
        print('win_angleM =', win_angleM, '; rot_speedM = ', rot_speedM)

    deltas = np.zeros(nb_param)
    deltas[:set1] = listDeltaGeo(param_geo['dist_PM'])
    deltas[set1:set2] = listDeltaGeo(param_geo['dist_MS'])
    deltas[set2:set3] = listDeltaGeo(param_geo['dist_SD'])
    deltas[set3:set4] = listDeltaTime(win_angleP, rot_speedP, win_angleM, rot_speedM)
    deltas[set4:] = listDeltaGeo(param_geo['angles'])
    ###########################################################################################
    if(verbose):
        print('deltas =', deltas)

    # Jacobian and parameter's uncertainty matrices
    jacobian = jacobianMatrix(dQx, dQy, dQz, dE)
    covxi = covxiMatrix(deltas)
    jacobian[0] = np.copy(dQx)
    jacobian[1] = np.copy(dQy)
    jacobian[2] = np.copy(dQz)
    jacobian[3] = np.copy(dE)
    ###########################################################################################
    if(verbose):
        print('Jacobian =', jacobian)
        print(covxi)

    jacobT = jacobian.T
    print(jacobian)
    print(jacobT)
    print(np.dot(jacobT, jacobian))
    global covQhw
    covQhw = np.dot(jacobian, np.dot(covxi, jacobT))
    print('ihikjdfkkahdfjlhadfklhjdafhjk')
    print(np.dot(jacobT, np.dot(covQhw, jacobian)))

print(covQhw)


#test perso

#a = {"P":[1,2], "U":[3,4]}
#print(len(a['P']))
#print(np.square(5))
#print(np.power(3,4))
#print(np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]))
#b = {'u':12, 'x':'f'}
#print('b =', b['u'])
#print(np.zeros((4,7)))
#print(np.eye(7))
#c = np.array([1,2])
#d = np.array([3,4])
#print(np.concatenate((c,d)))
#print(np.array([1]*3))
#e = np.array([5,6])
#print(np.concatenate((c,d,e)))


geo = {'dist_PM':[5, 0.001, 3, 0.001, 2, 0.001], 'dist_MS':[1.2, 0.002, 0.8, 0.002], 'dist_SD':[4, 0.003, 0, 0.003], 'angles':[0, 0.4, 0, 0.4, 0, 0.4]}
chop ={'chopperP':[500.0, 9.0, 7000.0, 17000.0, -1], 'chopperM':[800.0, 3.25, 7000.0, 17000.0, -1]}
covVcyl(geo,chop,2000,1000,True)

#print(covQhw)

print()