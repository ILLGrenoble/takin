/**
 * tlibs2 -- magnetic dynamics -- hamiltonian
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_MAGDYN_HAM_H__
#define __TLIBS2_MAGDYN_HAM_H__

#include <tuple>
#include <string>
#include <iostream>
#include <iomanip>

#include "../maths.h"
#include "../units.h"

#include "magdyn.h"


// --------------------------------------------------------------------
// calculation functions for interaction matrix and hamiltonian
// --------------------------------------------------------------------

/**
 * calculate the real-space interaction matrix J of
 * equations (10) - (13) from (Toth 2015)
 */
MAGDYN_TEMPL
t_mat MAGDYN_INST::CalcRealJ(const MAGDYN_TYPE::ExchangeTerm& term) const
{
	// symmetric part of the exchange interaction matrix, see (Toth 2015) p. 2
	t_mat J = tl2::diag<t_mat>(
		tl2::create<t_vec>({ term.J_calc, term.J_calc, term.J_calc }));

	// dmi as anti-symmetric part of interaction matrix
	// using a cross product matrix, see (Toth 2015) p. 2
	if(term.dmi_calc.size() == 3)
		J += tl2::skewsymmetric<t_mat, t_vec>(-term.dmi_calc);

	// general J matrix
	if(term.Jgen_calc.size1() == 3 && term.Jgen_calc.size2() == 3)
		J += term.Jgen_calc;

	// incommensurate case: rotation wrt magnetic unit cell
	// equations (21), (6), (2) as well as section 10 from (Toth 2015)
	if(IsIncommensurate())
	{
		const t_real rot_UC_angle =
			s_twopi * tl2::inner<t_vec_real>(m_ordering, term.dist_calc);

		if(!tl2::equals_0<t_real>(rot_UC_angle, m_eps))
		{
			t_mat rot_UC = tl2::convert<t_mat>(
				tl2::rotation<t_mat_real, t_vec_real>(
					m_rotaxis, rot_UC_angle));
			J = J * rot_UC;

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
			std::cout << "Coupling " << term.name << ": rot_UC =\n";
			tl2::niceprint(std::cout, rot_UC, 1e-4, 4);
#endif
		}
	}

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
	std::cout << "Coupling " << term.name << ": J =\n";
	tl2::niceprint(std::cout, J, 1e-4, 4);
#endif
	return J;
}



/**
 * calculate the reciprocal interaction matrices J(Q) and J(-Q) of
 * equations (12) and (14) from (Toth 2015)
 */
MAGDYN_TEMPL
std::tuple<MAGDYN_TYPE::t_Jmap, MAGDYN_TYPE::t_Jmap>
MAGDYN_INST::CalcReciprocalJs(const t_vec_real& Qvec) const
{
	t_Jmap J_Q{}, J_Q0{};

	// iterate couplings to pre-calculate corresponding J matrices
	for(const ExchangeTerm& term : GetExchangeTerms())
	{
		// insert or add an exchange matrix at the given site indices
		auto insert_or_add = [](t_Jmap& J, const t_indices& indices, const t_mat& J33)
		{
			if(auto iter = J.find(indices); iter != J.end())
				iter->second += J33;
			else
				J.emplace(std::move(std::make_pair(indices, J33)));
		};

		if(!CheckMagneticSite(term.site1_calc) || !CheckMagneticSite(term.site2_calc))
			continue;

		const t_indices indices = std::make_pair(term.site1_calc, term.site2_calc);
		const t_indices indices_t = std::make_pair(term.site2_calc, term.site1_calc);

		const t_mat J = CalcRealJ(term);
		if(J.size1() == 0 || J.size2() == 0)
			continue;

		// get J in reciprocal space by fourier trafo
		// equations (14), (12), (11), and (52) from (Toth 2015)
		const t_mat J_fourier = J * std::exp(m_phase_sign * s_imag * s_twopi *
			tl2::inner<t_vec_real>(term.dist_calc, Qvec));

		// symmetry: equation (15) from (Toth 2015)
		insert_or_add(J_Q, indices, J_fourier);
		insert_or_add(J_Q, indices_t, tl2::herm(J_fourier));

		insert_or_add(J_Q0, indices, J);
		insert_or_add(J_Q0, indices_t, tl2::herm(J));
	}  // end of iteration over couplings

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
	for(const auto& pair : J_Q)
	{
		std::cout << "Coupling J_Q[" << pair.first.first << ", " << pair.first.second << "] =\n";
		tl2::niceprint(std::cout, pair.second, 1e-4, 4);
	}
	for(const auto& pair : J_Q0)
	{
		std::cout << "Coupling J_Q0[" << pair.first.first << ", " << pair.first.second << "] =\n";
		tl2::niceprint(std::cout, pair.second, 1e-4, 4);
	}
#endif

	return std::make_tuple(J_Q, J_Q0);
}



/**
 * get the hamiltonian at the given momentum
 * @note implements the formalism given by (Toth 2015)
 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
 */
MAGDYN_TEMPL
t_mat MAGDYN_INST::CalcHamiltonian(const t_vec_real& Qvec) const
{
	const t_size N = GetMagneticSitesCount();
	if(N == 0)
		return t_mat{};

	// build the interaction matrices J(Q) and J(-Q) of
	// equations (12) and (14) from (Toth 2015)
	const auto [J_Q, J_Q0] = CalcReciprocalJs(Qvec);

	// create the hamiltonian of equation (25) and (26) from (Toth 2015)
	t_mat H00 = tl2::zero<t_mat>(N, N);
	t_mat HNN = tl2::zero<t_mat>(N, N);  // = H00*(-Q)
	t_mat H0N = tl2::zero<t_mat>(N, N);

	bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
		&& m_field.dir.size() == 3;

	// iterate magnetic sites
	for(t_size i = 0; i < N; ++i)
	{
		const MagneticSite& s_i = GetMagneticSite(i);

		// get the pre-calculated u and v vectors for the commensurate case
		const t_vec& u_i  = s_i.trafo_plane_calc;
		const t_vec& uc_i = s_i.trafo_plane_conj_calc;  // = u*_i
		const t_vec& v_i  = s_i.trafo_z_calc;

		for(t_size j = 0; j < N; ++j)
		{
			const MagneticSite& s_j = GetMagneticSite(j);

			// get the pre-calculated u and v vectors for the commensurate case
			const t_vec& u_j  = s_j.trafo_plane_calc;
			const t_vec& uc_j = s_j.trafo_plane_conj_calc;  // = u*_j
			const t_vec& v_j  = s_j.trafo_z_calc;

			// get the pre-calculated exchange matrices for the (i, j) coupling
			const t_indices indices_ij = std::make_pair(i, j);
			const t_mat* J_Q33 = nullptr;
			const t_mat* J_Q033 = nullptr;
			if(auto iter = J_Q.find(indices_ij); iter != J_Q.end())
				J_Q33 = &iter->second;
			if(auto iter = J_Q0.find(indices_ij); iter != J_Q0.end())
				J_Q033 = &iter->second;

			if(J_Q33)
			{
				// equation (26) from (Toth 2015)
				const t_real S_mag = 0.5 * std::sqrt(s_i.spin_mag_calc * s_j.spin_mag_calc);

				H00(i, j) += S_mag * tl2::inner<t_vec>(uc_i, (*J_Q33) * uc_j);
				HNN(i, j) += S_mag * tl2::inner<t_vec>(u_i,  (*J_Q33) * u_j);
				H0N(i, j) += S_mag * tl2::inner<t_vec>(uc_i, (*J_Q33) * u_j);
			}

			if(J_Q033)
			{
				// equation (26) from (Toth 2015)
				t_cplx c = s_j.spin_mag_calc * tl2::inner_noconj<t_vec>(v_i, (*J_Q033) * v_j);

				H00(i, i) -= c;
				HNN(i, i) -= c;
			}
		}  // end of iteration over j sites

		// include external field, equation (28) from (Toth 2015)
		if(use_field)
		{
			const t_vec field = tl2::convert<t_vec>(-m_field.dir) * m_field.mag;
			const t_vec gv    = s_i.g_e * v_i;
			const t_cplx Bgv  = tl2::inner_noconj<t_vec>(field, gv);

			// bohr magneton in [meV/T]
			constexpr const t_real muB = tl2::mu_B<t_real>
				/ tl2::meV<t_real> * tl2::tesla<t_real>;

			H00(i, i) -= muB * Bgv;
			HNN(i, i) -= std::conj(muB * Bgv);
		}
	}  // end of iteration over i sites

	// equation (25) from (Toth 2015)
	const t_mat HN0 = tl2::herm(H0N);
	t_mat H = tl2::create<t_mat>(2*N, 2*N);
	tl2::set_submat(H, H00, 0, 0); tl2::set_submat(H, H0N, 0, N);
	tl2::set_submat(H, HN0, N, 0); tl2::set_submat(H, HNN, N, N);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
	using namespace tl2_ops;
	std::cout << "H[ Q = " << Qvec << " ] =\n";
	tl2::niceprint(std::cout, H, 1e-4, 4);
#endif
	return H;
}



/**
 * get the energies from a hamiltonian
 * @note implements the formalism given by (Toth 2015)
 */
MAGDYN_TEMPL
MAGDYN_TYPE::EnergiesAndWeights MAGDYN_INST::CalcEnergiesFromHamiltonian(
	t_mat _H, const t_vec_real& Qvec, bool only_energies) const
{
	using namespace tl2_ops;
	const t_size N = GetMagneticSitesCount();
	if(N == 0 || _H.size1() == 0 || _H.size2() == 0)
		return EnergiesAndWeights{};

	// equation (30) from (Toth 2015)
	t_mat g_sign = tl2::unit<t_mat>(2*N);
	for(t_size i = N; i < g_sign.size1(); ++i)
		g_sign(i, i) = -1.;

	// equation (31) from (Toth 2015)
	t_mat chol_mat;
	t_size chol_try = 0;
	bool chol_failed = false;
	for(; chol_try < m_tries_chol; ++chol_try)
	{
		// upper cholesky decomposition: _H = _C^H _C
		const auto [chol_ok, _C] = tl2_la::chol<t_mat>(_H);

		if(chol_ok)
		{
			chol_mat = std::move(_C);
			break;
		}

		if(chol_try >= m_tries_chol - 1)
		{
			CERR_OPT << "Magdyn error: Cholesky decomposition failed"
				<< " at Q = " << Qvec << "." << std::endl;

			chol_mat = std::move(_C);
			chol_failed = true;
			break;
		}

		// try forcing the hamilton to be positive definite
		for(t_size i = 0; i < _H.size1(); ++i)
			_H(i, i) += m_delta_chol;
	}

	if(chol_failed || chol_mat.size1() == 0 || chol_mat.size2() == 0)
	{
		CERR_OPT << "Magdyn error: Invalid Cholesky decomposition"
			<< " at Q = " << Qvec << "." << std::endl;

		return EnergiesAndWeights{};
	}

	if(m_perform_checks && chol_try > 0)
	{
		CERR_OPT << "Magdyn warning: Needed " << chol_try
			<< " correction(s) for Cholesky decomposition"
			<< " at Q = " << Qvec << "." << std::endl;
	}

	// see p. 5 in (Toth 2015)
	t_mat H_mat = chol_mat * g_sign * tl2::herm<t_mat>(chol_mat);

	const bool is_herm = tl2::is_symm_or_herm<t_mat, t_real>(H_mat, m_eps);
	if(m_perform_checks && !is_herm)
	{
		CERR_OPT << "Magdyn warning: Hamiltonian is not hermitian"
			<< " at Q = " << Qvec << "." << std::endl;
	}

	// eigenvalues of the hamiltonian correspond to the energies
	// eigenvectors correspond to the spectral weights
	const auto [evecs_ok, evals, evecs] =
		tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(
			H_mat, only_energies, is_herm, true);
	if(!evecs_ok)
	{
		CERR_OPT << "Magdyn error: Eigensystem calculation failed"
			<< " at Q = " << Qvec << "." << std::endl;

		return EnergiesAndWeights{};
	}

	EnergiesAndWeights Es_and_Ws{};
	Es_and_Ws.reserve(evals.size());

	// register energies
	for(const auto& eval : evals)
	{
		const EnergyAndWeight EandS { .E = eval.real(), };
		Es_and_Ws.emplace_back(std::move(EandS));
	}

	// weight factors
	if(!only_energies)
	{
		bool corr_ok = CalcCorrelationsFromHamiltonian(Es_and_Ws,
			H_mat, chol_mat, g_sign, Qvec, evecs);
		if(!corr_ok)
		{
			CERR_OPT << "Magdyn warning: Invalid correlations"
				<< " at Q = " << Qvec << "." << std::endl;
		}
	}

	return Es_and_Ws;
}


#endif
