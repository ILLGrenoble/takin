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

#Constant in international system SI
m_nSI = 1.67492749804e-27 # mass of the neutron in kg
hSI = 6.62607015e-34 # Plank constant in J.s
hbarSI = np.divide(hSI, 2*np.pi) # J.s
mohSI = np.divide(m_nSI, hbarSI) # s/m^2

# Conversion
eV2J = 1.602176634e-19 # eV to J
meV2J = eV2J*1e-3
m2A = 1e10 # m to A
mm2m = 1e-3 # mm to m
# Conversion from atomic units
me2kg = 9.1093826e-31 # electron rest mass (me) to kg
mu2kg = 1.66053906892e-27 # unified atomic mass unit to kg (not using it in the following functions)
a02m = 5.2917721092e-11 # atomic length in m
a02A = a02m*m2A
hartree2J = 4.35974417e-18 # hatree to J
hartree2meV = np.divide(hartree2J,eV2J)*1e3 # hartree to meV
timeAU2s = np.divide(hbarSI,hartree2J) # time in AU to s
vAU2vSI = np.divide(a02m,timeAU2s) # velocity in AU to m/s

#Constant in atomic units AU
hbarAU = 1
meAU = 1
m_nAU = np.divide(m_nSI,me2kg)
mohAU = np.divide(m_nAU, hbarAU)
#print(me2kg, mu2kg, a02m, hartree2J, hartree2meV, timeAU2s, vAU2vSI, hbarAU, meAU, m_nAU, mohAU)

#Initialisation of the covariance Matrix
covQhw = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

# Give v in m/s using k in 1/A
def k2v(k):
    """k in 1/A"""
    km = k*m2A                  # 1/m
    return np.divide(km, mohSI)    # m/s

# Calcul of the length of a segment cuts in small segments
def calcDistLin(dist):
    """dist in an array: [L1, delta1, L2, delta2, ...], L: distance in mm, delta: uncerntainty"""
    L = 0
    nb = len(dist)
    for i in range(0,nb,2):
        L += dist[i]
    return L

def reducedCoef(v_i, v_f, L_PM, l_angles, shape, verbose=False):
    """v_i, v_f in AU, L_PM in AU, l_angles is a list of angles in radian, shape in ("SPHERE", "VCYL", "HCYL")"""
    theta_i = l_angles[0]       # radian
    phi_i = l_angles[1]         # radian
    theta_f = l_angles[2]       # radian
    phi_f = l_angles[3]         # radian

    Ae = np.divide(np.power(v_i,3), L_PM)
    Ax = np.divide(np.square(v_i), L_PM)*np.cos(theta_i)*np.cos(phi_i)
    Ay = np.divide(np.square(v_i), L_PM)*np.sin(theta_i)*np.cos(phi_i)
    Az = np.divide(np.square(v_i), L_PM)*np.sin(phi_i)
    Be = np.divide(np.power(v_f,3), L_PM)
    Bq = np.divide(np.square(v_f), L_PM)
    if shape in ("SPHERE", "VCYL"):
        Bx = Bq*np.cos(theta_f)*np.cos(phi_f)
        By = Bq*np.sin(theta_f)*np.cos(phi_f)
        Bz = Bq*np.sin(phi_f)
    elif(shape == "HCYL"):
        Bx = Bq*np.sin(phi_f)
        By = Bq*np.sin(theta_f)*np.cos(phi_f)
        Bz = Bq*np.cos(theta_f)*np.cos(phi_f)
    if(verbose):
        print('Ae =', Ae, '; Ax =', Ax, '; Ay =', Ay, '; Az =', Az)
        print('Be =', Be, '; Bx =', Bx, '; By =', By, '; Bz =', Bz, '\n')
    return(Ae, Ax, Ay, Az, Be, Bx, By, Bz) # Ae, Be : [L]^2/[T]^3 ; Ax, Ay, Az, Bx, By, Bz : [L]/[T]^2

# Calcul of Jacobian's terms for every shape
def jacobTerms(v_i, v_f, l_dist, l_angles, l_AB, l_sizes, shape, unit = 'SI', verbose=False):
    '''v_i, v_f are velocities, l_dist is a list of 6 distances, l_angles is a list of 4 angles in radian, l_AB is a list of 8 coef got from reducedCoef(), l_sizes is a list of 6 numbers, shape in (SPHERE, HCYL, VCYL), unit in (SI, AU)'''
    if(unit == 'SI'):
        m_n = m_nSI
        moh = mohSI
    elif(unit == 'AU'):
        m_n = m_nAU
        moh = mohAU

    L_PM = l_dist[0]
    L_MS = l_dist[1]
    L_SD = l_dist[2]
    rad = l_dist[4]
    theta_i = l_angles[0]
    phi_i = l_angles[1]
    theta_f = l_angles[2]
    phi_f = l_angles[3]
    Ae = l_AB[0]
    Ax = l_AB[1]
    Ay = l_AB[2]
    Az = l_AB[3]
    Be = l_AB[4]
    Bx = l_AB[5]
    By = l_AB[6]
    Bz = l_AB[7]
    size_PM = l_sizes[0]
    size_MS = l_sizes[1]

    dQx = np.array([])
    dQy = np.array([])
    dQz = np.array([])
    dE = np.array([])

    # x terms in common
    dQxdL_PMj = np.divide(moh, v_i)*(Ax + Bx*np.divide(L_MS, L_SD))
    dQxdL_MSj = -np.divide(moh, v_i)*Bx*np.divide(L_PM,L_SD)
    dQxdt_PM = -moh*(Ax + Bx*np.divide(L_MS, L_SD))
    dQxdt_MD = moh*Bx*np.divide(L_PM, L_SD)
    dQxdtheta_i = -moh*v_i*np.sin(theta_i)*np.cos(phi_i)
    dQxdphi_i = -moh*v_i*np.cos(theta_i)*np.sin(phi_i)
    # x terms that belong to a precise shape
    dQxdx = 0
    dQxdrad = 0
    dQxdz = 0
    dQxdtheta_f = 0
    dQxdphi_f = 0
    if(shape == 'SPHERE'):
        dQxdrad = -np.divide(moh, v_f)*Bx*np.divide(L_PM,L_SD)
        dQxdtheta_f = moh*v_f*np.sin(theta_f)*np.cos(phi_f)
        dQxdphi_f = moh*v_f*np.cos(theta_f)*np.sin(phi_f)
        dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f, dQxdphi_f]) ))
    elif(shape == 'HCYL'):
        dQxdx = -moh*np.divide(v_f, L_SD)
        dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdx, dQxdrad, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))
    elif(shape == 'VCYL'):
        dQxdrad = -np.divide(moh, v_f)*Bx*np.divide(L_PM,rad)
        dQxdtheta_f = moh*v_f*np.sin(theta_f)*np.cos(phi_f)
        dQx = np.concatenate(( np.array([dQxdL_PMj]*size_PM), np.array([dQxdL_MSj]*size_MS), np.array([dQxdrad, dQxdz, dQxdt_PM, dQxdt_MD, dQxdtheta_i, dQxdphi_i, dQxdtheta_f]) ))

    # y terms in common
    dQydL_PMj = np.divide(moh, v_i)*(Ay + By*np.divide(L_MS, L_SD))
    dQydL_MSj = -np.divide(moh, v_i)*By*np.divide(L_PM, L_SD)
    dQydt_PM = -moh*(Ay + By*np.divide(L_MS, L_SD))
    dQydt_MD = moh*By*np.divide(L_PM, L_SD)
    dQydtheta_i = moh*v_i*np.cos(theta_i)*np.cos(phi_i)
    dQydphi_i = -moh*v_i*np.sin(theta_i)*np.sin(phi_i)
    dQydtheta_f = -moh*v_f*np.cos(theta_f)*np.cos(phi_f)
    # y terms that belong to a precise shape
    dQydx = 0
    dQydrad = 0
    dQydz = 0
    dQydphi_f = 0
    if(shape == 'SPHERE'):
        dQydrad = -np.divide(moh, v_f)*By*np.divide(L_PM, L_SD)
        dQydphi_f = moh*v_f*np.sin(theta_f)*np.sin(phi_f)
        dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f, dQydphi_f]) ))
    elif(shape == 'HCYL'):
        dQydrad = -np.divide(moh, v_f)*By*np.divide(L_PM, rad)
        dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydx, dQydrad, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))
    elif(shape == 'VCYL'):
        dQydrad = -np.divide(moh, v_f)*By*np.divide(L_PM, rad)
        dQy = np.concatenate(( np.array([dQydL_PMj]*size_PM), np.array([dQydL_MSj]*size_MS), np.array([dQydrad, dQydz, dQydt_PM, dQydt_MD, dQydtheta_i, dQydphi_i, dQydtheta_f]) ))

    # z terms in common
    dQzdL_PMj = np.divide(moh, v_i)*(Az + Bz*np.divide(L_MS, L_SD))
    dQzdL_MSj = -np.divide(moh, v_i)*Bz*np.divide(L_PM, L_SD)
    dQzdt_PM = -moh*(Az + Bz*np.divide(L_MS, L_SD))
    dQzdt_MD = moh*Bz*np.divide(L_PM, L_SD)
    dQzdtheta_i = 0
    dQzdphi_i = moh*v_i*np.cos(phi_i)
    # z terms that belong to a precise shape
    dQzdx = 0
    dQzdrad = 0
    dQzdz = 0
    dQzdtheta_f = 0
    dQzdphi_f = 0
    if(shape == 'SPHERE'):
        dQzdrad = -np.divide(moh, v_f)*Bz*np.divide(L_PM, L_SD)
        dQzdphi_f = -moh*v_f*np.cos(phi_f)
        dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f, dQzdphi_f]) ))
    elif(shape == 'HCYL'):
        dQzdrad = -np.divide(moh, v_f)*Bz*np.divide(L_PM, rad)
        dQzdtheta_f = moh*v_f*np.sin(theta_f)*np.cos(phi_f)
        dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdx, dQzdrad, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))
    elif(shape == 'VCYL'):
        dQzdz = -moh*np.divide(v_f, L_SD)
        dQz = np.concatenate(( np.array([dQzdL_PMj]*size_PM), np.array([dQzdL_MSj]*size_MS), np.array([dQzdrad, dQzdz, dQzdt_PM, dQzdt_MD, dQzdtheta_i, dQzdphi_i, dQzdtheta_f]) ))

    # energy terms in common
    dEdL_PMj = np.divide(m_n, v_i)*(Ae + Be*np.divide(L_MS,L_SD))
    dEdL_MSj = -np.divide(m_n, v_i)*Be*np.divide(L_PM,L_SD)
    dEdt_PM = -m_n*(Ae + Be*np.divide(L_MS,L_SD))
    dEdt_MD = m_n*Be*np.divide(L_PM,L_SD)
    dEdtheta_i = 0
    dEdphi_i = 0
    dEdtheta_f = 0
    # energy terms that belong to a precise shape
    dEdx = 0
    dEdrad = 0
    dEdz = 0
    dEdphi_f = 0
    if(shape == 'SPHERE'):
        dEdrad = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)
        dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdrad, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f, dEdphi_f]) ))
    elif(shape == 'HCYL'):
        dEdx = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)*np.sin(phi_f)
        dEdrad = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)*np.cos(phi_f)
        dE = np.concatenate(( np.array([dEdL_PMj]*size_PM), np.array([dEdL_MSj]*size_MS), np.array([dEdx, dEdrad, dEdt_PM, dEdt_MD, dEdtheta_i, dEdphi_i, dEdtheta_f]) ))
    elif(shape == 'VCYL'):
        dEdrad = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)*np.cos(phi_f)
        dEdz = -np.divide(m_n, v_f)*Be*np.divide(L_PM,L_SD)*np.sin(phi_f)
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

# Calcul of the chopper's time uncertainty for a pair of counter rotating choppers
def deltaTimeChopper(window_angle, rot_speed, unit='SI', verbose=False):
    """angle in degree, rot_speed in RPM, unit in (SI, AU)"""
    deltaTpsChop = np.divide(window_angle, 2*6*rot_speed) # s
    ###########################################################################################
    if(verbose):
        print('deltaTpsChop in s =', deltaTpsChop)

    if( unit == 'AU'):
        deltaTpsChop = fromSItoAU(deltaTpsChop, 'time') # AU
    return deltaTpsChop

# Return the time uncerntainty for t_PM and t_MD
def listDeltaTime(window_angleP, rot_speedP, window_angleM, rot_speedM, dltD = 0, unit='SI', verbose = False):
    """angles in degree, rot_speed in RPM, dltD in second (uncertainty for the detector), unit in (SI, AU)"""
    dltP = deltaTimeChopper(window_angleP, rot_speedP, unit, verbose)      # AU
    dltM = deltaTimeChopper(window_angleM, rot_speedM, unit, verbose)      # AU
    ###########################################################################################
    if(verbose):
        print('\ndltP =', dltP, 'dltM =', dltM, '\n')

    return np.array([ np.sqrt(np.square(dltP) + np.square(dltM)), np.sqrt(np.square(dltM) + np.square(dltD)) ])

# Getting list of uncerntainty
def getDeltas(param_geo_u, param_choppers, l_sizes, unit='SI', verbose = False):
    '''param_geo_u is a dictionary: {dist_PM_u:[PM1, sigma1, PM2, sigma2, ...], dist_MS_u:[MS1, sigma1, MS2, sigma2, ...], dist_SD_u:[ (if HCYL: x, sigma_x), radius, sigma_r, (if VCYL: z, sigma_z)],
            angles_rad:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f, (if SPHERE: phi_f, sigma_phi_f)], delta_time_detector_u:value (0 by default)}, distances and time in AU, angles in rad
    param_choppers is a dictionary: {chopperP:[window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[window_angle, min_rot_speed, max_rot_speed, rot_speed]}, angles in degree, rot_speed in RPM
    l_sizes is a list of 6 numbers, unit in (SI, AU)'''
    set1 = l_sizes[0]
    set2 = l_sizes[0] + l_sizes[1]
    set3 = l_sizes[0] + l_sizes[1] + l_sizes[2]
    set4 = l_sizes[0] + l_sizes[1] + l_sizes[2] + l_sizes[3]
    win_angleP, rot_speedP = getParamChopper(param_choppers['chopperP']) # degree and RPM
    win_angleM, rot_speedM = getParamChopper(param_choppers['chopperM']) # degree and RPM
    ###########################################################################################
    if(verbose):
        print('win_angleP =', win_angleP, '; rot_speedP = ', rot_speedP)
        print('win_angleM =', win_angleM, '; rot_speedM = ', rot_speedM, '\n')

    deltas = np.zeros(l_sizes[5])
    deltas[:set1] = listDeltaGeo(param_geo_u['dist_PM_u'])
    deltas[set1:set2] = listDeltaGeo(param_geo_u['dist_MS_u'])
    deltas[set2:set3] = listDeltaGeo(param_geo_u['dist_SD_u'])
    deltas[set3:set4] = listDeltaTime(win_angleP, rot_speedP, win_angleM, rot_speedM, param_geo_u['delta_time_detector_u'], unit, verbose)
    deltas[set4:] = listDeltaGeo(param_geo_u['angles_rad'])
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

# To convert from SI to AU
def fromSItoAU(object, type_of_unit):
    '''object : the object to cenvert, type_of_unit : "time", "distance", "velocity"'''
    if(type_of_unit == 'time'):
        return np.divide(object, timeAU2s)
    elif(type_of_unit == 'distance'):
        return np.divide(object, a02m)
    elif(type_of_unit == 'velocity'):
        return np.divide(object, vAU2vSI)

# To get covQhw in A and eV
def fromCovQhwAUToCovQhwAeV(covQhwAU):
    covAmeV = np.zeros((4,4))
    covAmeV[0][0] = covQhwAU[0][0]*np.divide(1,np.square(a02A))
    covAmeV[0][1] = covQhwAU[0][1]*np.divide(1,np.square(a02A))
    covAmeV[0][2] = covQhwAU[0][2]*np.divide(1,np.square(a02A))
    covAmeV[0][3] = covQhwAU[0][3]*np.divide(hartree2meV,a02A)
    covAmeV[1][0] = covQhwAU[1][0]*np.divide(1,np.square(a02A))
    covAmeV[1][1] = covQhwAU[1][1]*np.divide(1,np.square(a02A))
    covAmeV[1][2] = covQhwAU[1][2]*np.divide(1,np.square(a02A))
    covAmeV[1][3] = covQhwAU[1][3]*np.divide(hartree2meV,a02A)
    covAmeV[2][0] = covQhwAU[2][0]*np.divide(1,np.square(a02A))
    covAmeV[2][1] = covQhwAU[2][1]*np.divide(1,np.square(a02A))
    covAmeV[2][2] = covQhwAU[2][2]*np.divide(1,np.square(a02A))
    covAmeV[2][3] = covQhwAU[2][3]*np.divide(hartree2meV,a02A)
    covAmeV[3][0] = covQhwAU[3][0]*np.divide(hartree2meV,a02A)
    covAmeV[3][1] = covQhwAU[3][1]*np.divide(hartree2meV,a02A)
    covAmeV[3][2] = covQhwAU[3][2]*np.divide(hartree2meV,a02A)
    covAmeV[3][3] = covQhwAU[3][3]*np.square(hartree2meV)
    return covAmeV

def dQiFromSItoAeV(dQx, dQy, dQz, dE):
    return (np.divide(dQx,m2A), np.divide(dQy,m2A), np.divide(dQz,m2A), np.divide(dE, meV2J))

def cov(param_geo, param_choppers, v_i, v_f, shape, unit='SI', verbose=False):
    """param_geo is a dictionary: {dist_PM:[PM1, sigma1, PM2, sigma2, ...], dist_MS:[MS1, sigma1, MS2, sigma2, ...], dist_SD:[ (if HCYL: x, sigma_x), radius, sigma_r, (if VCYL: z, sigma_z)], angles:[theta_i, sigma_theta_i, phi_i, sigma_phi_i, theta_f, sigma_theta_f, (if SPHERE: phi_f, sigma_phi_f)], delta_time_detector:value (0 by default)},
    param_choppers is a dictionary: {chopperP:[window_angle, min_rot_speed, max_rot_speed, rot_speed], chopperM:[window_angle, min_rot_speed, max_rot_speed, rot_speed]}, dist in mm, angles in degree, rot_speed in RPM, delta_time in s
    v_i, v_f: velocity of the incident and scattered neutron m/s,
    shape = SPHERE, VCYL, HCYL: shape of the detector (sphere, vertical cylinder or horizontal cylinder)"""

    if shape not in ('SPHERE', 'VCYL', 'HCYL'):
        print("this shape is not taken in account")
        return None
    ###########################################################################################
    if(verbose):
        print('param_geo =', param_geo)
        print('param_choppers =', param_choppers)
        print('v_i =', v_i, '; v_f =', v_f)
        print('det_shape =', shape, '; unit =', unit, '\n')

    # Storage of values given by the user
    vi = v_i                                                    # m/s
    vf = v_f                                                    # m/s
    dist_PM = np.multiply(param_geo['dist_PM'], mm2m)           # m
    dist_MS = np.multiply(param_geo['dist_MS'], mm2m)           # m
    dist_SD = np.multiply(param_geo['dist_SD'], mm2m)           # m
    angles = param_geo['angles']                                # degree
    delta_time_detector = param_geo['delta_time_detector']      # s
    ###########################################################################################
    if(verbose):
        print('dist_PM =', dist_PM)
        print('dist_MS =', dist_MS)
        print('dist_SD =', dist_SD)
        print('angles =', angles)
        print('delta_time_detector =', delta_time_detector, '\n')

    #Convertion to AU and radian
    if(unit == 'AU'):
        dist_PM = fromSItoAU(dist_PM, 'distance')                # AU
        dist_MS = fromSItoAU(dist_MS, 'distance')                # AU
        dist_SD = fromSItoAU(dist_SD, 'distance')                # AU
        delta_time_detector = fromSItoAU(delta_time_detector, 'time') #AU
        vi = fromSItoAU(vi, 'velocity')                         #AU
        vf = fromSItoAU(vf, 'velocity')                         #AU
    angles_rad = np.deg2rad(angles)                             #radian

    #dictionnary in AU
    param_geo_u = {'dist_PM_u':dist_PM, 'dist_MS_u':dist_MS, 'dist_SD_u':dist_SD, 'angles_rad':angles_rad, 'delta_time_detector_u': delta_time_detector}
    ###########################################################################################
    if(verbose):
        print('dist_PM_u =', dist_PM)
        print('dist_MS_u =', dist_MS)
        print('dist_SD_u =', dist_SD)
        print('angles_rad =', angles_rad)
        print('delta_time_detector_u =', delta_time_detector)
        print('vi =', vi)
        print('vf =', vf)
        print('param_geo_u =', param_geo_u, '\n')

    # Calcul of distances and angles for each shape of detectors
    L_PM = calcDistLin(dist_PM)
    L_MS = calcDistLin(dist_MS)
    x = 0
    rad = 0
    z = 0
    theta_i = 0
    phi_i = 0
    theta_f = 0
    phi_f = 0
    if(shape == "SPHERE"):
        # Distances
        L_SD = dist_SD[0]
        # Angles
        theta_i = angles_rad[0]
        phi_i = angles_rad[2]
        theta_f = angles_rad[4]
        phi_f = angles_rad[6]
    if(shape == "HCYL"):
        # Distances
        x = dist_SD[0]
        rad = dist_SD[2]
        L_SD = np.sqrt(np.square(x) + np.square(rad))
        # Angles
        theta_i = angles_rad[0]
        phi_i = angles_rad[2]
        theta_f = angles_rad[4]
        phi_f = np.arccos(np.divide(rad, L_SD))
    if(shape == "VCYL"):
        # Distances
        rad = dist_SD[0]
        z = dist_SD[2]
        L_SD = np.sqrt(np.square(rad) + np.square(z))
        # Angles
        theta_i = angles_rad[0]
        phi_i = angles_rad[2]
        theta_f = angles_rad[4]
        phi_f = np.arccos(np.divide(rad, L_SD))
    l_dist = np.array([L_PM, L_MS, L_SD, x, rad, z])
    l_angles = np.array([theta_i, phi_i, theta_f, phi_f])
    ###########################################################################################
    if(verbose):
        print('L_PM =', L_PM, '; L_MS =', L_MS, '; L_SD =', L_SD)
        print('x =', x, '; rad =', rad, '; z=', z, '; theta_i =', theta_i, '; phi_i =', phi_i, '; theta_f =', theta_f, '; phi_f =', phi_f, '\n')

    # Definition of variables depending on the shape of the detector
    Ae, Ax, Ay, Az, Be, Bx, By, Bz = reducedCoef(vi, vf, L_PM, l_angles, shape, verbose)
    l_AB = np.array([Ae, Ax, Ay, Az, Be, Bx, By, Bz])

    # Number of variables for each set of instrument's parameters (distances, times, angles)
    size_PM = int(np.divide(len(dist_PM), 2))
    size_MS = int(np.divide(len(dist_MS), 2))
    size_SD = int(np.divide(len(dist_SD), 2))
    size_tps = 2
    size_angles = int(np.divide(len(angles), 2))
    nb_param = size_PM + size_MS + size_SD + size_tps + size_angles
    l_sizes = np.array([size_PM, size_MS, size_SD, size_tps, size_angles, nb_param])
    ###########################################################################################
    if(verbose):
        print('size_PM =', size_PM, '; size_MS =', size_MS, '; size_SD =', size_SD, '; size_tps =', size_tps, '; size_angles =', size_angles, '\nl_sizes =', l_sizes, '\n')

    # Calcul of Jacobian's terms
    dQx, dQy, dQz, dE = jacobTerms(vi, vf, l_dist, l_angles, l_AB, l_sizes, shape, unit, verbose)
    if(unit == 'SI'):
        dQx, dQy, dQz, dE = dQiFromSItoAeV(dQx, dQy, dQz, dE)
    ###########################################################################################
    if(verbose):
        print('dQx =', dQx, '\ndQy =', dQy, '\ndQz =', dQz, '\ndE =', dE, '\n')

    # List of uncertainty
    deltas = getDeltas(param_geo_u, param_choppers, l_sizes, unit, verbose)
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
    covQhw_u = np.dot(jacobian, np.dot(covxi, jacobT))
    ###########################################################################################
    if(verbose):
        print('covQhwAU =', covQhw_u, '\n')

    if(unit == 'AU'):
        covQhw = fromCovQhwAUToCovQhwAeV(covQhw_u)
    elif(unit == 'SI'):
        covQhw = np.copy(covQhw_u)
    ###########################################################################################
    if(verbose):
        print('covQhw =', covQhw, '\n')

    return covQhw
