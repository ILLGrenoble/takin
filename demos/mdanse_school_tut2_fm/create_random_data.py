#
# create some random ferromagnon data to test the fitter
# @author Tobias Weber <tweber@ill.fr>
# @date 30-oct-2024
# @license see 'LICENSE' file
#

import numpy as np
import numpy.random as rnd
import scipy.constants as const


q     = 0.2        # reduced momentum transfer
G     = [1, 1, 0]  # Bragg peak
T     = 10.        # temperature
Emin  = -5.        # energy range
Emax  = 5.         #
pts   = 32         # number of points to generate
sig   = 0.5        # magnon peak with
amp   = 1.         # magnon peak amplitude
S0    = 100.       # S scaling factor
noise = 0.05       # noise on the magnon peaks
bck   = 10.        # background noise
ferro = True       # ferromagnetic or antiferromagnetic dispersion?
plot  = True       # plot the data
plot2 = False      # plot dispersion


# ferromagnetic dispersion
def disp_fm(Q):
	J = -1.  # exchange constant
	S =  1.  # total spin
	d =  1.  # distance between spins

	return 2*J*S * (np.cos(Q*d*2.*np.pi) - 1.)


# antiferromagnetic dispersion
def disp_afm(Q):
	J = 1.  # exchange constant
	S = 1.  # total spin
	d = 1.  # distance between spins

	return 2*J*S * np.abs(np.sin(Q*d*2.*np.pi))


def disp(Q):
	if ferro:
		return disp_fm(Q)
	return disp_afm(Q)


# Gaussian peak
def gauss(E, E0, sig, amp):
	norm = (np.sqrt(2.*np.pi) * sig)

	return amp * np.exp(-0.5*((E - E0)/sig)**2.) / norm


# Bose factor
def bose(E, T):
	# kB in meV/K
	kB = const.k / const.e * 1e3

	n = 1. / (np.exp(abs(E)/(kB*T)) - 1.)
	if E >= 0.:
		n += 1.
	return n


# Bose factor which is cut off below Ecut
def bose_cutoff(E, T, Ecut = 0.02):
	Ecut = abs(Ecut)

	if abs(E) < Ecut:
		b = bose(np.sign(E)*Ecut, T)
	else:
		b = bose(E, T)

	return b


Es = np.linspace(Emin, Emax, pts)
E0 = disp(q)

# add magnons
S = gauss(Es, E0, sig, amp)*bose_cutoff(E0, T) + rnd.default_rng().random(pts)*noise \
	+ gauss(Es, -E0, sig, amp)*bose_cutoff(-E0, T) + rnd.default_rng().random(pts)*noise
S *= S0

# add background
S += rnd.default_rng().random(pts) * bck
S = S.round()

# header for the data file to be read by takin
hdr = """
this is randomly generated data to test the convolution fitter

sample_a     = 5
sample_b     = 5
sample_c     = 5

sample_alpha = 90
sample_beta  = 90
sample_gamma = 90

mono_d       = 3.355
ana_d        = 3.355

sense_m      = +1
sense_s      = -1
sense_a      = +1

orient1_x    = 1
orient1_y    = 0
orient1_z    = 0

orient2_x    = 0
orient2_y    = 1
orient2_z    = 0

is_ki_fixed  = 0
k_fix        = 2.662

col_h        = 1
col_k        = 2
col_l        = 3
col_E        = 4
cols_scanned = 4
col_ctr      = 5
col_ctr_err  = 6
col_mon      = 7
col_mon_err  = 8
"""

hs = [q + G[0]] * pts
ks = [G[1]] * pts
ls = [G[2]] * pts
mon = [1e4] * pts
np.savetxt("random.dat", np.array([ hs, ks, ls, Es, S, np.sqrt(S), mon, np.sqrt(mon) ]).T, fmt = "%.4g", header = hdr)


if plot:
	import matplotlib.pyplot as plt

	plt.plot()
	plt.xlabel("E (meV)")
	plt.ylabel("Counts")

	plt.errorbar(Es, S, yerr = np.sqrt(S), fmt = "o")

	plt.tight_layout()
	plt.show()

	if plot2:
		plt.plot()
		plt.xlabel("Q (rlu)")
		plt.ylabel("E (meV)")
		Qs = np.linspace(-1, 1, 1024)
		Es_fm = disp_fm(Qs)
		Es_afm = disp_afm(Qs)
		plt.plot(Qs, Es_fm, lw = 2)
		plt.plot(Qs, Es_afm, lw = 2)
		plt.show()
