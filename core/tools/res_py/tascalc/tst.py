#
# calculates TAS angles from rlu (part of in20tools)
# @author Tobias Weber <tweber@ill.fr>
# @date 1-aug-18
# @license see 'LICENSE' file
#

import tas
import numpy as np
import numpy.linalg as la



# ------------------------------------------------------------------------------
# example calculations
# ------------------------------------------------------------------------------
if __name__ == "__main__":
	# --------------------------------------------------------------------------
	# lattice input
	# --------------------------------------------------------------------------
	lattice = np.array([5, 5, 5])
	angles = np.array([90, 90, 90])
	orient_rlu = np.array([1, 0, 0])
	orient2_rlu = np.array([0, 1, 0])
	# --------------------------------------------------------------------------

	# --------------------------------------------------------------------------
	# measurement position and instrument configuration input
	# --------------------------------------------------------------------------
	Q_rlu = np.array([1, -2, 0])
	E = 0.5
	kf = 2.662
	dmono = 3.355
	dana = 3.355
	# --------------------------------------------------------------------------

	# --------------------------------------------------------------------------
	# lattice and TAS angle calculation
	# --------------------------------------------------------------------------
	B = tas.get_B(lattice, angles/180.*np.pi)
	orient_up_rlu = tas.cross(orient_rlu, orient2_rlu, B)	# up vector in rlu

	ki = tas.get_ki(kf, E)
	[a1, a2] = tas.get_a1a2(ki, dmono)
	[a5, a6] = tas.get_a1a2(kf, dana)
	[a3, a4, dist_Q_plane] = tas.get_a3a4(ki, kf, Q_rlu, orient_rlu, orient_up_rlu, B)
	# --------------------------------------------------------------------------

	# --------------------------------------------------------------------------
	# output
	# --------------------------------------------------------------------------
	np.set_printoptions(suppress=True, precision=4)

	print("B [rlu -> 1/A] = \n" + str(B))
	print("scattering plane normal = " + str(orient_up_rlu/la.norm(orient_up_rlu)) + " rlu")
	print("\na1 = %.4f deg, a2 = %.4f deg, a3 = %.4f deg, a4 = %.4f deg, a5 = %.4f deg, a6 = %.4f deg" \
		% (a1/np.pi*180., a2/np.pi*180., a3/np.pi*180., a4/np.pi*180., a5/np.pi*180., a6/np.pi*180.))
	# --------------------------------------------------------------------------


	# --------------------------------------------------------------------------
	# CHECK: reproducing input values
	# --------------------------------------------------------------------------
	ki = tas.get_monok(a1, dmono)
	kf = tas.get_monok(a5, dana)
	E = tas.get_E(ki, kf)
	Qlen = tas.get_Q(ki, kf, a4)
	Qvec = tas.get_hkl(ki, kf, a3, Qlen, orient_rlu, orient_up_rlu, B)
	# --------------------------------------------------------------------------

	# --------------------------------------------------------------------------
	# output
	# --------------------------------------------------------------------------
	print("ki = %.4f 1/A, kf = %.4f 1/A, E = %.4f meV, |Q| = %.4f 1/A, "\
		"Q = %s rlu" % (ki, kf, E, Qlen, Qvec))
	# --------------------------------------------------------------------------


	# --------------------------------------------------------------------------
	# angle between two reciprocal vectors
	metric = tas.get_metric(B)
	vec1  = tas.np.array([1., 0., 0.])
	vec2  = tas.np.array([0., 1., 0.])
	ang = tas.angle(vec1, vec2, metric)
	print("\nAngle between %s rlu and %s rlu: %.4f deg" % (str(vec1), str(vec2), ang/np.pi*180.))
	# --------------------------------------------------------------------------

# ------------------------------------------------------------------------------
