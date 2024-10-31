#
# testing the lswt algorithm (https://arxiv.org/abs/1402.6069) directly for comparison,
# 1d ferromagnetic chain
# @authors : Tobias Weber <tweber@ill.fr> & Victor Mecoli <mecoli@ill.fr>
# @date 31-oct-2024
# @license GPLv3, see 'LICENSE' file
#

import math as m

import numpy as np
import numpy.linalg as la
from numpy import array	# in global namespace so that Takin can access it

import matplotlib.pyplot as plt

import scipy as sp
import scipy.constants as const



#is_ferromagnetic = True  # choose ferromagnetic or antiferromagnetic 1d spin chain
only_pos_E = True        # hide magnon annihilation?
verbose_print = False    # print intermediate results


def print_infos(str):
	if verbose_print:
		print(str)


#
# magnetic sites
# "S": spin magnitude
# "Sdir": spin direction
#
sites = [
	{ "S" : 1., "Sdir" : [ 0, 0, 1 ] },
]

#
# magnetic couplings
# "sites": indices of the sites to couple
# "J": (symmetric) exchange interaction
# "DMI": (antisymmetric) Dzyaloshinskii-Moryia interaction
# "dist": distance in rlu to the next unit cell for the coupling
#
couplings = [
	{ "sites" : [ 0, 0 ], "J" : -1., "dist" : [ 1, 0, 0 ] },
]


# calculate spin rotations towards ferromagnetic order along [001]
for site in sites:
	#Sdir = np.array(site["Sdir"]) / la.norm(site["Sdir"])
	c = 1.
	rot = np.diag([ c, c, c ])
	site["u"] = rot[0, :] + 1j * rot[1, :]
	site["v"] = rot[2, :]

#	print_infos(np.dot(rot, Sdir))
	print_infos("\nrot = \n%s\nu = %s\nv = %s" % (rot, site["u"], site["v"]))



# calculate real interaction matrices
for coupling in couplings:
	J = coupling["J"]
	coupling["J_real"] = np.diag([ J, J, J ])
	print_infos("\nJ_real =\n%s" % coupling["J_real"])


# get the energies of the dispersion at the momentum transfer Qvec
def get_energies(Qvec):
	print_infos("\n\nQ = %s" % Qvec)

	# fourier transform interaction matrices
	num_sites = len(sites)
	J_fourier = np.zeros((num_sites, num_sites, 3, 3), dtype = complex)
	J0_fourier = np.zeros((num_sites, num_sites, 3, 3), dtype = complex)

	for coupling in couplings:
		dist = np.array(coupling["dist"])
		J_real = coupling["J_real"]
		site1 = coupling["sites"][0]
		site2 = coupling["sites"][1]

		J_ft = J_real * np.exp(-1j * 2.*np.pi * np.dot(dist, Qvec))
		J_fourier[site1, site2] += J_ft
		J_fourier[site2, site1] += J_ft.transpose().conj()
		J0_fourier[site1, site2] += J_real
		J0_fourier[site2, site1] += J_real.transpose().conj()

	print_infos("\nJ_fourier =\n%s\n\nJ0_fourier =\n%s" % (J_fourier, J0_fourier))


	# hamiltonian
	H = np.zeros((2*num_sites, 2*num_sites), dtype = complex)

	for i in range(num_sites):
		S_i = sites[i]["S"]
		u_i = sites[i]["u"]
		v_i = sites[i]["v"]

		for j in range(num_sites):
			S_j = sites[j]["S"]
			u_j = sites[j]["u"]
			v_j = sites[j]["v"]
			S = 0.5 * np.sqrt(S_i * S_j)

			H[i, j] += S * np.dot(u_i, np.dot(J_fourier[i, j], u_j.conj()))
			H[i, i] -= S_j * np.dot(v_i, np.dot(J0_fourier[i, j], v_j))
			H[num_sites + i, num_sites + j] += \
				S * np.dot(u_i.conj(), np.dot(J_fourier[i, j], u_j))
			H[num_sites + i, num_sites + i] -= \
				S_j * np.dot(v_i, np.dot(J0_fourier[i, j], v_j))
			H[i, num_sites + j] += \
				S * np.dot(u_i, np.dot(J_fourier[i, j], u_j))
			H[num_sites + i, j] += \
				(S * np.dot(u_j, np.dot(J_fourier[j, i], u_i))).conj()

	print_infos("\nH =\n%s" % H)


	# trafo
	C = la.cholesky(H)
	signs = np.diag(np.concatenate((np.repeat(1, num_sites), np.repeat(-1, num_sites))))
	H_trafo = np.dot(C.transpose().conj(), np.dot(signs, C))
	print_infos("\nC =\n%s\n\nH_trafo =\n%s" % (C, H_trafo))


	# the eigenvalues of H give the energies
	Es = np.real(la.eigvals(H_trafo))
	print_infos("\nEs = %s" % Es)

	return Es

# -----------------------------------------------------------------------------
# dispersion
# -----------------------------------------------------------------------------

# kB in meV/K
kB = const.k / const.e * 1e3


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

# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Takin interface
# -----------------------------------------------------------------------------

# global variables which can be accessed / changed by Takin
#g_G = np.array([4., 4., 0.])	# Bragg peak

g_h = 4.              #
g_k = 4.              # Bragg peak (hkl)
g_l = 0.              #

g_inc_sig = 0.02      # incoherent width
g_inc_amp = 1.        # incoherent intensity

g_T = 300.            # temperature
g_bose_cut = 0.01     # lower cutoff energy for the Bose factor


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
		E_peak = get_energies(Q)
	except ZeroDivisionError:
		return [0., 0.]

	return [[E_peak, -E_peak], [w_peak, w_peak]]


#
# S(Q,E) function, called for every Monte-Carlo point
#
def TakinSqw(h, k, l, E):
	try:
		[Ep_peak, Em_peak], [wp_peak, wm_peak] = TakinDisp(h,k,l)

		incoh = gauss(E, 0., g_inc_sig, g_inc_amp)

		S = incoh
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
