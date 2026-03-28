#
# Example Sunny interface module
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
# @date 2-sep-2024
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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


# -----------------------------------------------------------------------------
# helper functions and constants
# -----------------------------------------------------------------------------
# from scipy.constants
kB = 0.08617
ge = -2.00232


#
# Gaussian peak
#
function gauss(x, x0, sig, amp)
	norm = (sqrt(2.0*pi) * sig)
	return amp * exp(-0.5*((x-x0)/sig)^2.0) / norm
end


#
# Bose factor
#
function bose(E, T)
	n = 1.0/(exp(abs(E)/(kB*T)) - 1.0)
	if E >= 0.0
		n += 1.0
	end
	return n
end


#
# Bose factor which is cut off below Ecut
#
function bose_cutoff(E, T, Ecut=0.02)
	Ecut = abs(Ecut)

	b = 0.0
	if abs(E) < Ecut
		b = bose(sign(E)*Ecut, T)
	else
		b = bose(E, T)
	end

	return b
end
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# global variables which can be accessed / changed by Takin
# -----------------------------------------------------------------------------
g_sig      = 0.02       # magnon linewidth
g_S0       = 10.0       # S(Q, E) scaling factor

g_inc_sig  = 0.02       # incoherent width
g_inc_amp  = 10.0       # incoherent intensity

g_T        = 100.0      # temperature
g_bose_cut = 0.02       # Bose cutoff

g_S        = 1.0        # spin magnitudes
g_J        = -1.0       # coupling constant
g_negative = true       # repeat dispersion for -E
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# interface with takin and sunny
# -----------------------------------------------------------------------------
using Sunny

num_sites = 1
calc = nothing


#
# the init function is called after Takin has changed a global variable (optional)
#
function TakinInit()
	println("TakinInit: Setting up magnetic model...")

	# set up magnetic sites and xtal lattice
	magsites = Crystal(
		lattice_vectors(5, 5, 5, 90, 90, 90),
		[
			# site list
			[ 0, 0, 0 ],
		], 1)

	global num_sites = length(magsites.positions)


	# set up spin magnitudes and magnetic system
	magsys = System(magsites, ( 1, 1, 1 ),
		[
			SpinInfo(1, S = g_S, g = -[ ge 0 0; 0 ge 0; 0 0 ge ]),
		], :dipole)


	# spin directions
	set_dipole!(magsys, [ 0, 0, 1 ], ( 1, 1, 1, 1 ))


	# set up magnetic couplings
	set_exchange!(magsys,
		[
			g_J     0     0;
			  0   g_J     0;
			  0     0   g_J
		], Bond(1, 1, [ 1, 0, 0 ]))


	# set up spin-wave calculator
	global calc = SpinWaveTheory(magsys; apply_g = true)

	println("TakinInit: Finished setting up magnetic model.")
end


#
# dispersion E(Q) and weight factor (optional)
#
function TakinDisp(h::Float64, k::Float64, l::Float64)
	# momentum transfer
	Q = vec([h, k, l])

	momenta = collect(range(Q, Q, 1))
	energies, correlations = intensities_bands(calc, momenta,
		intensity_formula(calc, :trace; kernel = delta_function_kernel))

	return [ energies[1, :], correlations[1, :] / Float64(num_sites) ]
end


#
# called for every Monte-Carlo point
#
function TakinSqw(h::Float64, k::Float64, l::Float64, E::Float64)::Float64
	#println("Calling TakinSqw(", h, ", ", k, ", ", l, ", ", E, ") -> ", S)

	# magnon branches
	S_mag = 0.
	Es, ws = TakinDisp(h, k, l)
	for e_idx in 1 : length(Es)
		E_branch = Es[e_idx]
		w_branch = ws[e_idx]
		S_mag += gauss(E, E_branch, g_sig, g_S0*w_branch)

		if g_negative
			S_mag += gauss(E, -E_branch, g_sig, g_S0*w_branch)			
		end
	end

	# incoherent scattering
	incoh = gauss(E, 0.0, g_inc_sig, g_inc_amp)

	b = bose_cutoff(E, g_T, g_bose_cut)
	S = S_mag*b + incoh
	return Float64(S)
end


#
# background function, called for every nominal (Q, E) point (optional)
#
function TakinBackground(h::Float64, k::Float64, l::Float64, E::Float64)::Float64
	#println("Calling TakinBackground(", h, ", ", k, ", ", l, ", ", E, ") -> ", S)
	return Float64(0.)
end
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# test
# -----------------------------------------------------------------------------
#TakinInit()
#println(TakinDisp(1.1, 0.0, 0.0))
#println(TakinSqw(1.1, 0.0, 0.0, 0.4))
# -----------------------------------------------------------------------------
