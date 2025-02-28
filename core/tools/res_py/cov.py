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


### See definition in my notes (Important to modify this line when a publication is done) /!\/!\/!\/!\/!\

import numpy as np
import numpy.linalg as la

# Conversion
eV2J = 1.602176634e-19 # eV to J
m2mm = 1e3 # m to mm
m2A = 1e10 # m to A
mm2A = 1e7 # mm to A

# Constant
m_n = 1.67492749804e-27 # mass of the neutron in kg
h_SI = 6.62607015e-34 # Plank constant in J.s
h_eV = np.divide(h_SI, eV2J) # Plank constant in eV.s
hbar_SI = np.divide(h_SI, 2*np.pi) # J.s
hbar_eV = np.divide(hbar_SI, eV2J) # eV.s
moh = np.divide(m_n, hbar_SI) # s/m^2
mohA = np.divide(moh, np.square(m2A)) # s/A^2
print(m_n, h_SI, hbar_SI, h_eV, hbar_eV, moh, mohA)

#Initialisation of the covariance Matrix
covQhw = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

# Give v in m/s using k in 1/A
def k2v(k):
    """k in 1/A"""
    km = k*m2A                  # 1/m
    return np.divide(k, moh)    # m/s

# Calcul of the length of a segment cuts in small segments
def calcDistLin(dist):
    """dist in an array: [L1, delta1, L2, delta2, ...], L: distance in mm, delta: uncerntainty"""
    L = 0
    nb = len(dist)
    for i in range(0,nb,2):
        L += dist[i]
    return L

def reducedCoef(v_i, v_f, L_PM, l_angles, shape, verbose=False):
    """v_i, v_f in m/s, L_PM in mm, l_angles is a list of angles in radian, shape in ("SPHERE", "VCYL", "HCYL")"""
    theta_i = l_angles[0]                                                   # radian
    phi_i = l_angles[1]                                                     # radian
    theta_f = l_angles[2]                                                   # radian
    phi_f = l_angles[3]                                                     # radian
    v_iA = v_i*m2A                                                          # A/s
    v_fA = v_f*m2A                                                          # A/s
    L_PMA = L_PM*mm2A                                                       # A

    Ae = np.divide(np.power(v_iA,3), L_PMA)                                 # A^2/s^3
    Ax = np.divide(np.square(v_iA), L_PMA)*np.cos(theta_i)*np.cos(phi_i)    # A/s^2
    Ay = np.divide(np.square(v_iA), L_PMA)*np.sin(theta_i)*np.cos(phi_i)    # A/s^2
    Az = np.divide(np.square(v_iA), L_PMA)*np.sin(phi_i)                    # A/s^2
    Be = np.divide(np.power(v_fA,3), L_PMA)                                 # A^2/s^3
    Bq = np.divide(np.square(v_fA), L_PMA)                                  # A/s^2
    if shape in ("SPHERE", "VCYL"):
        Bx = Bq*np.cos(theta_f)*np.cos(phi_f)                               # A/s^2
        By = Bq*np.sin(theta_f)*np.cos(phi_f)                               # A/s^2
        Bz = Bq*np.sin(phi_f)                                               # A/s^2
    elif(shape == "HCYL"):
        Bx = Bq*np.sin(phi_f)                                               # A/s^2
        By = Bq*np.sin(theta_f)*np.cos(phi_f)                               # A/s^2
        Bz = Bq*np.cos(theta_f)*np.cos(phi_f)                               # A/s^2
    if(verbose):
        print('Ae =', Ae, '; Ax =', Ax, '; Ay =', Ay, '; Az =', Az)
        print('Be =', Be, '; Bx =', Bx, '; By =', By, '; Bz =', Bz, '\n')
    return(Ae, Ax, Ay, Az, Be, Bx, By, Bz)

# Calcul of Jacobian's terms for every shape
def jacobTerms(v_i, v_f, l_dist, l_angles, l_AB, l_sizes, shape, verbose=False):
    '''v_i, v_f are velocities in m/s, l_dist is a list of 6 distances, l_AB is a list of 8 coef got from reducedCoef(), l_sizes is a list of 6 numbers, shape in (SPHERE, HCYL, VCYL)'''
    v_imm = v_i*m2mm                                                                     # mm/s
    v_fmm = v_f*m2mm                                                                     # mm/s
    v_iA = v_i*m2A                                                                       # A/s
    v_fA = v_f*m2A                                                                       # A/s
    L_PM = l_dist[0]                                                                     # mm
    L_MS = l_dist[1]                                                                     # mm
    L_SD = l_dist[2]                                                                     # mm
    rad = l_dist[4]                                                                      # mm
    theta_i = l_angles[0]                                                                # radian
    phi_i = l_angles[1]                                                                  # radian
    theta_f = l_angles[2]                                                                # radian
    phi_f = l_angles[3]                                                                  # radian
    Ae = l_AB[0]                                                                         # A^2/s^3
    Ax = l_AB[1]                                                                         # A/s^2
    Ay = l_AB[2]                                                                         # A/s^2
    Az = l_AB[3]                                                                         # A/s^2
    Be = l_AB[4]                                                                         # A^2/s^3
    Bx = l_AB[5]                                                                         # A/s^2
    By = l_AB[6]                                                                         # A/s^2
    Bz = l_AB[7]                                                                         # A/s^2
    Aem = np.divide(Ae,np.square(m2A))                                                   # m^2/s^3
    Bem = np.divide(Be,np.square(m2A))                                                   # m^2/s^3
    size_PM = l_sizes[0]
    size_MS = l_sizes[1]
    dQx = np.array([])
    dQy = np.array([])
    dQz = np.array([])
    dE = np.array([])

    # x terms in common
    dQxdL_PMj = np.divide(mohA, v_imm)*(Ax + Bx*np.divide(L_MS, L_SD))                   # 1/(A*mm)
    dQxdL_MSj = -np.divide(mohA, v_imm)*Bx*np.divide(L_PM,L_SD)                          # 1/(A*mm)
    dQxdt_PM = -mohA*(Ax + Bx*np.divide(L_MS, L_SD))                                     # 1/(A*s)
    dQxdt_MD = mohA*Bx*np.divide(L_PM, L_SD)                                             # 1/(A*s)
    dQxdtheta_i = -mohA*v_iA*np.sin(theta_i)*np.cos(phi_i)                               # 1/A
    dQxdphi_i = -mohA*v_iA*np.cos(theta_i)*np.sin(phi_i)                                 # 1/A
    # x terms that belong to a precise shape
    dQxdx = 0
    dQxdrad = 0
    dQxdz = 0
    dQxdtheta_f = 0
    dQxdphi_f = 0
    if(shape == 'SPHERE'):
        dQxdrad = -np.divide(mohA, v_fmm)*Bx*np.divide(L_PM,L_SD)                        # 1/(A*mm)
        dQxdtheta_f = mohA*v_fA*np.sin(theta_f)*np.cos(phi_f)                            # 1/A
        dQxdphi_f = mohA*v_fA*np.cos(theta_f)*np.sin(phi_f)                              # 1/A
        dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f, dQxdphi_f]) ))
    elif(shape == 'HCYL'):
        dQxdx = -moh*np.divide(v_fmm, L_SD)*np.divide(1,m2A*m2mm)                        # 1/(A*mm)
        dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdx, dQxdrad, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))
    elif(shape == 'VCYL'):
        dQxdrad = -np.divide(mohA, v_fmm)*Bx*np.divide(L_PM,rad)                         # 1/(A*mm)
        dQxdtheta_f = mohA*v_fA*np.sin(theta_f)*np.cos(phi_f)                            # 1/A
        dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdz, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))
    
    # y terms in common
    dQydL_PMj = np.divide(mohA, v_imm)*(Ay + By*np.divide(L_MS, L_SD))                   # 1/(A*mm)
    dQydL_MSj = -np.divide(mohA, v_imm)*By*np.divide(L_PM, L_SD)                         # 1/(A*mm)
    dQydt_PM = -mohA*(Ay + By*np.divide(L_MS, L_SD))                                     # 1/(A*s)
    dQydt_MD = mohA*By*np.divide(L_PM, L_SD)                                             # 1/(A*s)
    dQydtheta_i = mohA*v_iA*np.cos(theta_i)*np.cos(phi_i)                                # 1/A
    dQydphi_i = -mohA*v_iA*np.sin(theta_i)*np.sin(phi_i)                                 # 1/A
    dQydtheta_f = -mohA*v_fA*np.cos(theta_f)*np.cos(phi_f)                               # 1/A
    # y terms that belong to a precise shape
    dQydx = 0
    dQydrad = 0
    dQydz = 0
    dQydphi_f = 0
    if(shape == 'SPHERE'):
        dQydrad = -np.divide(mohA, v_fmm)*By*np.divide(L_PM, L_SD)                       # 1/(A*mm)
        dQydphi_f = mohA*v_fA*np.sin(theta_f)*np.sin(phi_f)                              # 1/A
        dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f, dQydphi_f]) ))
    elif(shape == 'HCYL'):
        dQydrad = -np.divide(mohA, v_fmm)*By*np.divide(L_PM, rad)                        # 1/(A*mm)
        dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydx, dQydrad, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))
    elif(shape == 'VCYL'):
        dQydrad = -np.divide(mohA, v_fmm)*By*np.divide(L_PM, rad)                        # 1/(A*mm)
        dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydz, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))

    # z terms in common
    dQzdL_PMj = np.divide(mohA, v_imm)*(Az + Bz*np.divide(L_MS, L_SD))                   # 1/(A*mm)
    dQzdL_MSj = -np.divide(mohA, v_imm)*Bz*np.divide(L_PM, L_SD)                         # 1/(A*mm)
    dQzdt_PM = -mohA*(Az + Bz*np.divide(L_MS, L_SD))                                     # 1/(A*s)
    dQzdt_MD = mohA*Bz*np.divide(L_PM, L_SD)                                             # 1/(A*s)
    dQzdtheta_i = 0
    dQzdphi_i = mohA*v_iA*np.cos(phi_i)                                                  # 1/A
    # z terms that belong to a precise shape
    dQzdx = 0
    dQzdrad = 0
    dQzdz = 0
    dQzdtheta_f = 0
    dQzdphi_f = 0
    if(shape == 'SPHERE'):
        dQzdrad = -np.divide(mohA, v_fmm)*Bz*np.divide(L_PM, L_SD)                       # 1/(A*mm)
        dQzdphi_f = -mohA*v_fA*np.cos(phi_f)                                             # 1/A
        dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f, dQzdphi_f]) ))
    elif(shape == 'HCYL'):
        dQzdrad = -np.divide(mohA, v_fmm)*Bz*np.divide(L_PM, rad)                        # 1/(A*mm)
        dQzdtheta_f = mohA*v_fA*np.sin(theta_f)*np.cos(phi_f)                            # 1/A
        dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdx, dQzdrad, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))
    elif(shape == 'VCYL'):
        dQzdz = -moh*np.divide(v_fmm, L_SD)*np.divide(1,m2A*m2mm)                        # 1/(A*mm)
        dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdz, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))

    # energy terms in common
    dEdL_PMj = np.divide(m_n, v_imm)*(Aem + Bem*np.divide(L_MS,L_SD))*np.divide(1,eV2J)  # eV/mm
    dEdL_MSj = -np.divide(m_n, v_imm)*Bem*np.divide(L_PM,L_SD)*np.divide(1,eV2J)         # eV/mm
    dEdt_PM = -m_n*(Aem + Bem*np.divide(L_MS,L_SD))*np.divide(1,eV2J)                    # eV/s
    dEdt_MD = m_n*Bem*np.divide(L_PM,L_SD)*np.divide(1,eV2J)                             # eV/s
    dEdtheta_i = 0
    dEdphi_i = 0
    dEdtheta_f = 0
    # energy terms that belong to a precise shape
    dEdx = 0
    dEdrad = 0
    dEdz = 0
    dEdphi_f = 0
    if(shape == 'SPHERE'):
        dEdrad = -np.divide(m_n, v_fmm)*Bem*np.divide(L_PM,L_SD)                         # eV/mm
        dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdrad, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f, dEdphi_f]) ))
    elif(shape == 'HCYL'):
        dEdx = -np.divide(m_n, v_fmm)*Bem*np.divide(L_PM,L_SD)*np.sin(phi_f)             # eV/mm
        dEdrad = -np.divide(m_n, v_fmm)*Bem*np.divide(L_PM,L_SD)*np.cos(phi_f)           # eV/mm
        dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdx, dEdrad, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f]) ))
    elif(shape == 'VCYL'):
        dEdrad = -np.divide(m_n, v_fmm)*Bem*np.divide(L_PM,L_SD)*np.cos(phi_f)           # eV/mm
        dEdz = -np.divide(m_n, v_fmm)*Bem*np.divide(L_PM,L_SD)*np.sin(phi_f)             # eV/mm
        dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdrad, dEdz, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f]) ))
    ###########################################################################################
    if(verbose):
        print('x terms:')
        print('dQx/dL_PMj =', dQxdL_PMj, '; dQx/dL_MSj =', dQxdL_MSj, '; dQx/dx =', dQxdx, '; dQx/drad =', dQxdrad)
        print('dQx/dt_PM =', dQxdt_PM, '; dQx/dt_MD =', dQxdt_MD)
        print('dQx/dtheta_i =', dQxdtheta_i, '; dQx/dphi_i =', dQxdphi_i, '; dQx/dtheta_f =', dQxdtheta_f, '; dQx/dphi_f =', dQxdphi_f, '\n')

        print('y terms:')
        print('dQy/dL_PMj =', dQydL_PMj, '; dQy/dL_MSj =', dQydL_MSj, '; dQy/drad =', dQydrad)
        print('dQy/dt_PM =', dQydt_PM, '; dQy/dt_MD =', dQydt_MD)
        print('dQy/dtheta_i =', dQydtheta_i, '; dQy/dphi_i =', dQydphi_i, '; dQy/dtheta_f =', dQydtheta_f, '; dQy/dphi_f =', dQydphi_f, '\n')

        print('z terms:')
        print('dQz/dL_PMj =', dQzdL_PMj, '; dQz/dL_MSj =', dQzdL_MSj, '; dQz/drad =', dQzdrad, '; dQz/dz =', dQzdz)
        print('dQz/dt_PM =', dQzdt_PM, '; dQz/dt_MD =', dQzdt_MD)
        print('dQz/dphi_i =', dQzdphi_i, 'dQz/dphi_f =', dQzdphi_f, '; dQzdtheta_f =', dQzdtheta_f, '\n')

        print('E  terms:')
        print('dE/dL_PMj =', dEdL_PMj, '; dE/dL_MSj =', dEdL_MSj, '; dE/dx =', dEdx, '; dE/drad =', dEdrad, '; dE/dz =', dEdz)
        print('dE/dt_PM =', dEdt_PM, '; dE/dt_MD =', dEdt_MD, '\n')
    return(dQx, dQy, dQz, dE)


# Getting parameters of choppers
def getParamChopper(l_param_chopper):
    '''l_param_chopper is a np.array [win_angle, min_rot_speed, max_rot_speed, rot_speed], win_angle in degree, *rot_speed in RPM'''
    win_angle = l_param_chopper[0]          # degree
    rot_speed = l_param_chopper[3]          # RPM
    if(rot_speed == -1):
        rot_speed = l_param_chopper[1]      # degree
    return win_angle, rot_speed

# Storage of the given uncerntainty in a list
def listDeltaGeo(list_param):
    """list_param is an array: [P1, delta1, P2, delta2, ...], P: parameter, delta: uncerntainty"""
    nb = len(list_param)
    dlt = np.zeros(int(nb/2))
    for i in range(1,nb,2):
        ind = int((i-1)/2)
        dlt[ind]=list_param[i]
    return dlt

# Calcul of the chopper's time uncertainty
def deltaTimeChopper(window_angle, rot_speed):
    """angle in degree, rot_speed in RPM"""
    return( np.divide(window_angle, 6*rot_speed) ) # result in seconds

# Return the time uncerntainty for t_PM and t_MD
def listDeltaTime(window_angleP, rot_speedP, window_angleM, rot_speedM, dltD = 0):
    """angles in degree, rot_speed in RPM, dltD in second (uncertainty for the detector)"""
    dltP = deltaTimeChopper(window_angleP, rot_speedP)      # s
    dltM = deltaTimeChopper(window_angleM, rot_speedM)      # s
    return np.array([ np.sqrt(np.square(dltP) + np.square(dltM)), np.sqrt(np.square(dltM) + np.square(dltD)) ])

# Getting list of uncerntainty
def getDeltas(param_geo, param_choppers, l_sizes, verbose = False):
    '''param_geo is a dictionary: {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[ (if HCYL: x, sigma_x), radius, sigma_r, (if VCYL: z, sigma_z)], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f, (if SPHERE: phi_f, sigma_phi_f)], delta_time_detector:value (0 by default)},
    param_choppers is a dictionary: {chopperP:[window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[window_angle, min_rot_speed, max_rot_speed, rot_speed]}, dist in m, angle in degree, rot_speed in RPM
    l_sizes is a list of 6 numbers'''
    set1 = l_sizes[0]
    set2 = l_sizes[0] + l_sizes[1]
    set3 = l_sizes[0] + l_sizes[1] + l_sizes[2]
    set4 = l_sizes[0] + l_sizes[1] + l_sizes[2] + l_sizes[3]
    win_angleP, rot_speedP = getParamChopper(param_choppers['chopperP'])
    win_angleM, rot_speedM = getParamChopper(param_choppers['chopperM'])

    deltas = np.zeros(l_sizes[5])
    deltas[:set1] = listDeltaGeo(param_geo['dist_PM'])
    deltas[set1:set2] = listDeltaGeo(param_geo['dist_MS'])
    deltas[set2:set3] = listDeltaGeo(param_geo['dist_SD'])
    deltas[set3:set4] = listDeltaTime(win_angleP, rot_speedP, win_angleM, rot_speedM, param_geo['delta_time_detector'])
    deltas[set4:] = listDeltaGeo(np.deg2rad(param_geo['angles']))
    ###########################################################################################
    if(verbose):
        print('win_angleP =', win_angleP, '; rot_speedP = ', rot_speedP)
        print('win_angleM =', win_angleM, '; rot_speedM = ', rot_speedM, '\n')
    return deltas

# Creation and filling of the Jacobian matrix
def jacobianMatrix(dQx, dQy, dQz, dE):
    """dQi and dE are lists of derivations in respect to every parameters"""
    nb_col = len(dQx)
    jacob = np.zeros((4, nb_col))
    jacob[0] = np.copy(dQx)
    jacob[1] = np.copy(dQy)
    jacob[2] = np.copy(dQz)
    jacob[3] = np.copy(dE)
    return jacob

# Creation and filling of the covxi matrix (uncerntainty on independant parameters xi)
def covxiMatrix(deltas):
    """delta is a list of uncerntainties (for independant variables)"""
    nb_dlt = len(deltas)
    covxi = np.eye(nb_dlt)
    for i in range(nb_dlt):
        covxi[i][i] = np.square(deltas[i])
    return covxi

def cov(param_geo, param_choppers, v_i, v_f, shape, verbose=False):
    """param_geo is a dictionary: {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[ (if HCYL: x, sigma_x), radius, sigma_r, (if VCYL: z, sigma_z)], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f, (if SPHERE: phi_f, sigma_phi_f)], delta_time_detector:value (0 by default)},
    param_choppers is a dictionary: {chopperP:[window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[window_angle, min_rot_speed, max_rot_speed, rot_speed]}, dist in mm, angles in degree, rot_speed in RPM,
    v_i, v_f: velocity of the incident and scattered neutron m/s,
    shape = SPHERE, VCYL, HCYL: shape of the detector (sphere, vertical cylinder or horizontal cylinder)"""

    if shape not in ('SPHERE', 'VCYL', 'HCYL'):
        print("this shape is not taken in account")
        return None
    
    # Calcul of distances and angles for each shape of detectors
    L_PM = calcDistLin(param_geo['dist_PM'])                    # mm
    L_MS = calcDistLin(param_geo['dist_MS'])                    # mm
    x = 0
    rad = 0
    z = 0
    theta_i = 0
    phi_i = 0
    theta_f = 0
    phi_f = 0
    if(shape == "SPHERE"):
        # Distances
        L_SD = param_geo['dist_SD'][0]                          # mm
        # Angles
        theta_i = np.deg2rad(param_geo['angles'][0])            # radian
        phi_i = np.deg2rad(param_geo['angles'][2])              # radian
        theta_f = np.deg2rad(param_geo['angles'][4])            # radian
        phi_f = np.deg2rad(param_geo['angles'][6])              # radian
    if(shape == "HCYL"):
        # Distances
        x = param_geo['dist_SD'][0]                             # mm
        rad = param_geo['dist_SD'][2]                           # mm
        L_SD = np.sqrt(np.square(x) + np.square(rad))           # mm
        # Angles
        theta_i = np.deg2rad(param_geo['angles'][0])            # radian
        phi_i = np.deg2rad(param_geo['angles'][2])              # radian
        theta_f = np.deg2rad(param_geo['angles'][4])            # radian
        phi_f = np.acos(np.divide(rad, L_SD))                   # radian
    if(shape == "VCYL"):
        # Distances
        rad = param_geo['dist_SD'][0]                           # mm
        z = param_geo['dist_SD'][2]                             # mm
        L_SD = np.sqrt(np.square(rad) + np.square(z))           # mm
        # Angles
        theta_i = np.deg2rad(param_geo['angles'][0])            # radian
        phi_i = np.deg2rad(param_geo['angles'][2])              # radian
        theta_f = np.deg2rad(param_geo['angles'][4])            # radian
        phi_f = np.acos(np.divide(rad, L_SD))                   # radian
    l_dist = np.array([L_PM, L_MS, L_SD, x, rad, z])
    l_angles = np.array([theta_i, phi_i, theta_f, phi_f])
    ###########################################################################################
    if(verbose):
        print('param_geo =', param_geo)
        print('param_choppers =', param_choppers)
        print('v_i =', v_i, '; v_f =', v_f)
        print('L_PM =', L_PM, '; L_MS =', L_MS, '; L_SD =', L_SD)
        print('x =', x, '; rad =', rad, '; z=', z, '; theta_i =', theta_i, '; phi_i =', phi_i, '; theta_f =', theta_f, '; phi_f =', phi_f, '\n')

    # Definition of variables depending on the shape of the detector
    Ae, Ax, Ay, Az, Be, Bx, By, Bz = reducedCoef(v_i, v_f, L_PM, l_angles, shape, verbose)
    l_AB = np.array([Ae, Ax, Ay, Az, Be, Bx, By, Bz])

    # Number of variables for each set of instrument's parameters (distances, times, angles)
    size_PM = int(np.divide(len(param_geo['dist_PM']), 2))
    size_MS = int(np.divide(len(param_geo['dist_MS']), 2))
    size_SD = int(np.divide(len(param_geo['dist_SD']), 2))
    size_tps = 2
    size_angles = int(np.divide(len(param_geo['angles']), 2))
    nb_param = size_PM + size_MS + size_SD + size_tps + size_angles
    l_sizes = np.array([size_PM, size_MS, size_SD, size_tps, size_angles, nb_param])
    ###########################################################################################
    if(verbose):
        print('size_PM =', size_PM, '; size_MS =', size_MS, '; size_SD =', size_SD, '; size_tps =', size_tps, '; size_angles =', size_angles, '\n')

    # Calcul of Jacobian's terms
    dQx, dQy, dQz, dE = jacobTerms(v_i, v_f, l_dist, l_angles, l_AB, l_sizes, shape, verbose)
    ###########################################################################################
    if(verbose):
        print('dQx =', dQx, '\ndQy =', dQy, '\ndQz =', dQz, '\ndE =', dE, '\n')
    
    # List of uncertainty
    deltas = getDeltas(param_geo, param_choppers, l_sizes, verbose)
    ###########################################################################################
    if(verbose):
        print('deltas =', deltas, '\n')

    # Jacobian and parameter's uncertainty matrices
    jacobian = jacobianMatrix(dQx, dQy, dQz, dE)
    covxi = covxiMatrix(deltas)
    ###########################################################################################
    if(verbose):
        print('Jacobian =', jacobian, '\n')
        print('covxi = ', covxi, '\n')

    jacobT = jacobian.T
    global covQhw
    covQhw = np.dot(jacobian, np.dot(covxi, jacobT))
    ###########################################################################################
    if(verbose):
        print('covQhw =', covQhw)

    return covQhw

#print('Vertical cylinder\n')
# Test for a vertical cylindrical detector
# covVcyl(param_geo, param_choppers, v_i, v_f, verbose=False):
# param_geo : {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[radius, sigma_r, z_heigh, sigma_z], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f]}
# param_choppers : {chopperP:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed]}

#geoV = {'dist_PM':[5, 0.001, 3, 0.001, 2, 0.001], 'dist_MS':[1.2, 0.002, 0.8, 0.002], 'dist_SD':[4, 0.003, 0, 0.003], 'angles':[0, 0.4, 0, 0.4, 0, 0.4]}
#chopV ={'chopperP':[500.0, 9.0, 7000.0, 17000.0, -1], 'chopperM':[800.0, 3.25, 7000.0, 17000.0, -1]}
#covVcyl(geoV,chopV,2000,1000,True)

#print('\nHorizontal cylinder\n')
# Test for a horizontal cylindrical detector
# covHcyl(param_geo, param_choppers, v_i, v_f, verbose=False):
# param_geo : {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[x_width, sigma_x, radius, sigma_r], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f]}
# param_choppers : {chopperP:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed]}

#geoH = {'dist_PM':[5, 0.001, 3, 0.001, 2, 0.001], 'dist_MS':[1.2, 0.002, 0.8, 0.002], 'dist_SD':[0, 0.003, 4, 0.003], 'angles':[0, 0.4, 0, 0.4, 0, 0.4]}
#chopH ={'chopperP':[500.0, 9.0, 7000.0, 17000.0, -1], 'chopperM':[800.0, 3.25, 7000.0, 17000.0, -1]}
#covHcyl(geoH,chopH,2000,1000,True)

#print('\nSphere\n')
# Test for a spherical detector
# covSph(param_geo, param_choppers, v_i, v_f, verbose=False):
# param_geo : {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[radius, sigma_r], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f, phi_f, sigma_phi_f]}
# param_choppers : {chopperP:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[diameter, window_angle, min_rot_speed, max_rot_speed, rot_speed]}

#geoS = {'dist_PM':[5, 0.001, 3, 0.001, 2, 0.001], 'dist_MS':[1.2, 0.002, 0.8, 0.002], 'dist_SD':[4, 0.003], 'angles':[0, 0.4, 0, 0.4, 0, 0.4, 0, 0.4]}
#chopS ={'chopperP':[500.0, 9.0, 7000.0, 17000.0, -1], 'chopperM':[800.0, 3.25, 7000.0, 17000.0, -1]}
#covSph(geoS,chopS,2000,1000,True)

# print(reducedCoef(1,1,3,[0,0,0,0],"VCYL"))
# print(np.zeros(5))
# a = np.array([])
# print(a)
# a = np.concatenate(( np.array([0,1]), np.array([2,3]) ))
# print(a)