#
# Testing the lswt algorithm (https://arxiv.org/abs/1402.6069) directly for comparison,
# 1d ferromagnetic chain
# @authors : Tobias Weber <tweber@ill.fr> & Victor Mecoli <mecoli@ill.fr>
# @date 31-oct-2024
# @license GPLv3, see 'LICENSE' file
#

import numpy as np
import numpy.linalg as la
from numpy import array	 # in global namespace so that Takin can access it

verbose_print = False    # print intermediate results


def print_infos(str):
	if verbose_print:
		print(str)



# sites and couplings are defined in TakinInit
site = {}
coupling = {}



#----------------------------------------------------------------------------------------------------
# Energy
#----------------------------------------------------------------------------------------------------

# Get the energies of the dispersion at the momentum transfer Qvec - correct when you have only one site
def get_energies(Qvec):
	print_infos("\n\nQ = %s" % Qvec)

	#--------------------------
	# Initialisation:
	#--------------------------

	# Fourier transform interaction matrices
	J_fourier = np.zeros((1, 1, 3, 3), dtype = complex)
	J0_fourier = np.zeros((1, 1, 3, 3), dtype = complex)

	# Parameters of the ferromagnetic chain
	S_i = site["S"]
	u_i = site["u"]
	v_i = site["v"]

	S_j = site["S"]
	u_j = site["u"]
	v_j = site["v"]

	S = 0.5 * np.sqrt(S_i * S_j)

	dist = np.array(coupling["dist"])
	J = coupling["J"]
	site1 = coupling["sites"][0]
	site2 = coupling["sites"][1]

	# Hamiltonian
	H = np.zeros((2*1, 2*1), dtype = complex)

	#--------------------------
	# Calculation:
	#--------------------------

	J_ft = J * np.exp(-1j * 2.*np.pi * np.dot(dist, Qvec))
	J_fourier[site1, site2] += J_ft
	J_fourier[site2, site1] += J_ft.transpose().conj()
	J0_fourier[site1, site2] += J
	J0_fourier[site2, site1] += J.transpose().conj()

	print_infos("\nJ_fourier =\n%s\n\nJ0_fourier =\n%s" % (J_fourier, J0_fourier))


	H[0, 0] += S * np.dot(u_i, np.dot(J_fourier[0, 0], u_j.conj()))
	H[0, 0] -= S_j * np.dot(v_i, np.dot(J0_fourier[0, 0], v_j))

	H[1, 1] += S * np.dot(u_i.conj(), np.dot(J_fourier[0, 0], u_j))
	H[1, 1] -= S_j * np.dot(v_i, np.dot(J0_fourier[0, 0], v_j))

	H[0, 1] += S * np.dot(u_i, np.dot(J_fourier[0, 0], u_j))
	H[1, 0] += (S * np.dot(u_j, np.dot(J_fourier[0, 0], u_i))).conj()

	print_infos("\nH =\n%s" % H)


	# Trafo
	C = la.cholesky(H)
	signs = np.array([[ 1, 0 ], [ 0, -1 ]])
	H_trafo = np.dot(C.transpose().conj(), np.dot(signs, C))

	print_infos("\nC =\n%s\n\nH_trafo =\n%s" % (C, H_trafo))


	# The eigenvalues of H give the energies
	Es = np.real(la.eigvals(H_trafo))

	print_infos("\nEs = %s" % Es)

	return Es




#----------------------------------------------------------------------------------------------------
# dispersion
#----------------------------------------------------------------------------------------------------

# kB in meV/K, value from scipy
kB = 0.08617333
# kB = const.k / const.e * 1e3


# Bose factor
def bose(E, T):
	n = 1./(np.exp(abs(E)/(kB*T)) - 1.)
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
	norm = (np.sqrt(2.*np.pi) * sig)
	return amp * np.exp(-0.5*((x-x0)/sig)**2.) / norm




#----------------------------------------------------------------------------------------------------
# Takin interface
#----------------------------------------------------------------------------------------------------

# global variables which can be accessed / changed by Takin
g_sig = 0.5           # magnon linewidth
g_amp = 2.            # magnon amplitude

g_inc_sig = 0.02      # incoherent width
g_inc_amp = 0.        # incoherent intensity

g_T = 300.            # temperature
g_bose_cut = 0.01     # lower cutoff energy for the Bose factor



#
# the init function is called after Takin has changed a global variable (optional)
#
def TakinInit():
	#
	# Magnetic site - ferromagnetic case
	# "S": spin magnitude
	# "u" and "v": spin rotations towards ferromagnetic order along [001]
	#
	global site
	site = { "S" : 1., "u" : np.array([ 1.+0.j, 0.+1.j, 0.+0.j ]),  "v" : np.array([ 0., 0., 1. ]) }

	print_infos("\nu = %s\nv = %s" % (site["u"], site["v"]))


	#
	# Magnetic coupling - ferromagnetic case
	# "sites": indices of the sites to couple
	# "J": (symmetric) exchange interaction
	# "dist": distance in reduced lattice unit (rlu) to the next unit cell for the coupling
	# "J": real interaction matrices
	#
	Jc = -1.
	global coupling
	coupling = { "sites" : [ 0, 0 ], "dist" : [ 1, 0, 0 ], "J" : np.diag([ Jc, Jc, Jc ]) }

	print("TakinInit: T = %g, J =\n%s" % (g_T, coupling["J"]))



#
# dispersion E(Q) and weight factor (optional)
#
def TakinDisp(h, k, l):
	E_peak = 0.		# energy
	w_peak = 1.		# weight

	try:
		Q = np.array([h, k, l])
		E_peaks = get_energies(Q)
	except (ZeroDivisionError, la.LinAlgError):
		return [[0.], [0.]]

	return [E_peaks, [w_peak] * len(E_peaks)]



#
# S(Q,E) function, called for every Monte-Carlo point
#
def TakinSqw(h, k, l, E):
	try:
		S = 0.

		# magnon peaks
		[E_peaks, w_peaks] = TakinDisp(h,k,l)
		for (E_peak, w_peak) in zip(E_peaks, w_peaks):
			S += gauss(E, E_peak, g_sig, g_amp) * w_peak * bose_cutoff(E_peak, g_T, g_bose_cut)

		# incoherent peak
		S += gauss(E, 0., g_inc_sig, g_inc_amp)

		return S
	except ZeroDivisionError:
		return 0.

#----------------------------------------------------------------------------------------------------




#
# this will be executed when the module loads
#
import os
print("Script working directory: " + os.getcwd())


#
# this python file can also be called directly for testing or plotting
#
if __name__ == "__main__":
	TakinInit()

	# plotting the dispersion
	qs = []
	Es = []
	for q in np.linspace(-1, 1, 128):
		for E in TakinDisp(q, 0, 0)[0]:
			if E < 0.:
				continue
			qs.append(q)
			Es.append(E)

	try:
		import matplotlib.pyplot as plt

		plt.xlabel("q (rlu)")
		plt.ylabel("E (meV)")
		plt.plot(qs, Es)
		plt.show()
	except ModuleNotFoundError:
		print("Could not plot dispersion.")
