#
# Sample Python S(q,w) module for simple acoustic phonons
#
# @author Tobias Weber <tweber@ill.fr>
# @license GPLv2
# @date nov-2019
#

import math as m

import numpy as np
import numpy.linalg as la
from numpy import array	# in global namespace so that Takin can access it

import scipy as sp
import scipy.constants as const



# -----------------------------------------------------------------------------
# dispersion
# -----------------------------------------------------------------------------

# kB in meV/K
kB = const.k / const.e * 1e3


# dispersion relations
def disp_phonon(q, amp, freq, E0):
	return np.abs(amp*np.sin(freq*q)) + E0


#
# power law
# see: https://en.wikipedia.org/wiki/Landau_theory
#
def power_law(T, Tc, crit_exp, amp):
	if T > Tc:
		return amp * np.abs(T - Tc)**crit_exp
	else:
		return 0.


# Bose factor
def bose(E, T):
	n = 1./(m.exp(abs(E)/(kB*T)) - 1.)
	if E >= 0.:
		n += 1.
	return n

# Bose factor which is cut off below Ecut
def bose_cutoff(E, T, Ecut=0.02):
	Ecut = abs(Ecut)

	if abs(E) < Ecut:
		b = bose(np.sign(E)*Ecut, T)
	else:
		b = bose(E, T)

	return b



# Gaussian peak
def gauss(x, x0, sig, amp):
	norm = (np.sqrt(2.*m.pi) * sig)
	return amp * np.exp(-0.5*((x-x0)/sig)**2.) / norm


#
# peak shape of a damped harmonic oscillator
# see: B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108
#
def DHO(E, T, E0, hwhm, amp):
	#E0 = np.sqrt(E0**2 - hwhm**2)
	return np.abs(bose(E, T)*amp/(E0*np.pi) *
		(hwhm/((E-E0)**2. + hwhm**2.) - hwhm/((E+E0)**2. + hwhm**2.)))

# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Takin interface
# -----------------------------------------------------------------------------

# global variables which can be accessed / changed by Takin
#g_G = np.array([4., 4., 0.])	# Bragg peak

g_h = 4.              #
g_k = 4.              # Bragg peak (hkl)
g_l = 0.              #

g_amp = 20.           # amplitude of sinusoidal dispersion
g_freq = np.pi/2.     # frequency of sinusoidal dispersion
g_offs = 0.           # energy gap
g_HWHM = 0.02         # linewidth
g_S0 = 1.             # intensity

g_inc_sig = 0.02      # incoherent width
g_inc_amp = 1.        # incoherent intensity

g_T = 300.            # temperature
g_bose_cut = 0.01     # lower cutoff energy for the Bose factor

#g_use_powerlaw = 0    # alternatively use a power law for the energy gap
#g_Tc = 50.            # critical temperature
#g_critexp = 0.5       # critical exponent
#g_amp_powerlaw = 1.0  # scaling factor of the power law


#
# the init function is called after Takin has changed a global variable (optional)
#
def TakinInit():
	print("Init: G=(%.2f %.2f %.2f), T=%.2f" % (g_h, g_k, g_l, g_T))


#
# dispersion E(Q) and weight factor (optional)
#
def TakinDisp(h, k, l):
	E_peak = 0.		# energy
	w_peak = 1.		# weight

	try:
		Q = np.array([h, k, l])
		G = np.array([g_h, g_k, g_l])
		q = la.norm(Q - G)

		E0 = g_offs
		#if g_use_powerlaw:
		#	E0 = power_law(g_T, g_Tc, g_critexp, g_amp_powerlaw)
		E_peak = disp_phonon(q, g_amp, g_freq, E0)
	except ZeroDivisionError:
		return [0., 0.]

	return [[E_peak, -E_peak], [w_peak, w_peak]]


#
# S(Q,E) function, called for every Monte-Carlo point
#
def TakinSqw(h, k, l, E):
	try:
#		print("h={0}, k={1}, l={2}, E={3}".format(h,k,l,E))

		[Ep_peak, Em_peak], [wp_peak, wm_peak] = TakinDisp(h,k,l)

#		S_p = DHO(E, g_T, Ep_peak, g_HWHM, g_S0*wp_peak)
#		S_m = DHO(E, g_T, Em_peak, g_HWHM, g_S0*wm_peak)
#
#		S = (S_p + S_m)*bose_cutoff(E, g_T, g_bose_cut) + incoh

		S = np.abs(DHO(E, g_T, Ep_peak, g_HWHM, g_S0*wp_peak))
		incoh = gauss(E, 0., g_inc_sig, g_inc_amp)

		S = S + incoh
#		print("S={0}".format(S))
		return S
	except ZeroDivisionError:
		return 0.

# -----------------------------------------------------------------------------



#
# this will be executed when the module loads
#
import os
print("Script working directory: " + os.getcwd())



#
# this python file can also be called directly for testing or plotting
#
if __name__ == "__main__":
	# testing the TakinSqw() function
	print(TakinSqw(4.1, 3.9, 0., 0.4))

	# plotting the transverse dispersion branch
	qs = np.linspace(-0.75, 0.75, 64)
	Es_plus = []
	Es_minus = []
	for q in qs:
		[[E_p, E_m], [w_p, w_m]] = TakinDisp(g_h+q, g_k-q, g_l)
		Es_plus.append(E_p)
		Es_minus.append(E_m)

	try:
		import matplotlib.pyplot as plt

		plt.xlabel("q (rlu)")
		plt.ylabel("E (meV)")
		plt.plot(qs, Es_plus)
		plt.plot(qs, Es_minus)
		plt.show()
	except ModuleNotFoundError:
		print("Could not plot dispersion.")
