/**
 * tlibs2 -- magnetic dynamics -- spin-spin correlation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (McClarty 2022) https://doi.org/10.1146/annurev-conmatphys-031620-104715
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

#ifndef __TLIBS2_MAGDYN_CORREL_H__
#define __TLIBS2_MAGDYN_CORREL_H__

#include <vector>
#include <iostream>
#include <iomanip>

#include "../algos.h"
#include "../maths.h"
#include "../phys.h"

#include "magdyn.h"


// --------------------------------------------------------------------
// calculation functions
// --------------------------------------------------------------------

/**
 * get the dynamical structure factor from a hamiltonian
 * @note implements the formalism given by (Toth 2015)
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CalcCorrelationsFromHamiltonian(MAGDYN_TYPE::SofQE& S) const
{
	using namespace tl2_ops;

	const t_size N = GetMagneticSitesCount();
	if(N == 0)
		return false;

	SortByEnergies(S);

	// create a matrix of eigenvectors
	std::vector<t_vec> evecs;
	evecs.reserve(S.E_and_S.size());
	for(t_size idx = 0; idx < S.E_and_S.size(); ++idx)
		evecs.push_back(S.E_and_S[idx].state);
	S.evec_mat = tl2::create<t_mat>(evecs);

	if(m_perform_checks)
	{
		// check commutator relations, see before equ. 7 in (McClarty 2022)
		//t_mat check_comm = S.evec_mat * S.comm * tl2::herm(S.evec_mat);
		t_mat check_comm = S.evec_mat * tl2::herm(S.evec_mat);
		if(!tl2::equals(/*S.comm*/ tl2::unit<t_mat>(2*N), check_comm, m_eps))
		{
			CERR_OPT << "Magdyn error: Wrong commutator at Q = "
				<< S.Q_rlu << ": " << check_comm
				<< "." << std::endl;
		}
	}

	// equation (32) from (Toth 2015)
	const t_mat energy_mat = tl2::herm(S.evec_mat) * S.H_comm * S.evec_mat;  // energies
	t_mat E_sqrt = S.comm * energy_mat;          // abs. energies
	for(t_size i = 0; i < E_sqrt.size1(); ++i)
		E_sqrt(i, i) = std::sqrt(E_sqrt(i, i));  // sqrt. of abs. energies

	// re-create energies, to be consistent with the weights
	if(energy_mat.size1() != S.E_and_S.size())
	{
		CERR_OPT << "Magdyn warning: Expected " << S.E_and_S.size() << " energies at Q = "
			<< S.Q_rlu << ", but got " << energy_mat.size1() << " energies"
			<< "." << std::endl;

		S.E_and_S.resize(energy_mat.size1());
	}

	for(t_size i = 0; i < energy_mat.size1(); ++i)
	{
		if(m_perform_checks && !tl2::equals_0(energy_mat(i, i).imag(), m_eps))
		{
			CERR_OPT << "Magdyn warning: Remaining imaginary energy component at Q = "
				<< S.Q_rlu << " and E = " << energy_mat(i, i)
				<< "." << std::endl;
		}

		if(m_perform_checks && !tl2::equals(energy_mat(i, i).real(), S.E_and_S[i].E, m_eps))
		{
			CERR_OPT << "Magdyn warning: Mismatching energy at Q = "
				<< S.Q_rlu << " and E = " << energy_mat(i, i).real()
				<< ", expected E = " << S.E_and_S[i].E
				<< "." << std::endl;
		}

		S.E_and_S[i].E = energy_mat(i, i).real();
		S.E_and_S[i].S = tl2::zero<t_mat>(3, 3);
		S.E_and_S[i].S_perp = tl2::zero<t_mat>(3, 3);
	}

	const auto [chol_inv, inv_ok] = tl2::inv(S.H_chol);
	if(!inv_ok)
	{
		CERR_OPT << "Magdyn error: Cholesky inversion failed"
			<< " at Q = " << S.Q_rlu << "." << std::endl;
		return false;
	}

	// equation (34) from (Toth 2015)
	const t_mat trafo = chol_inv * S.evec_mat * E_sqrt;
	const t_mat trafo_herm = tl2::herm(trafo);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
	t_mat D_mat = trafo_herm * S.H_comm * trafo;
	std::cout << "D =\n";
	tl2::niceprint(std::cout, D_mat, 1e-4, 4);
	std::cout << "E_sqrt =\n";
	tl2::niceprint(std::cout, E_sqrt, 1e-4, 4);
	std::cout << "L_energy =\n";
	tl2::niceprint(std::cout, energy_mat, 1e-4, 4);
#endif

	// building the spin correlation functions of equation (47) from (Toth 2015)
	for(std::uint8_t x_idx = 0; x_idx < 3; ++x_idx)
	for(std::uint8_t y_idx = 0; y_idx < 3; ++y_idx)
	{
		// equations (44) and (47) from (Toth 2015)
		t_mat M = tl2::create<t_mat>(2*N, 2*N);

		for(t_size i = 0; i < N; ++i)
		for(t_size j = 0; j < N; ++j)
		{
			// get the sites
			const MagneticSite& s_i = GetMagneticSite(i);
			const MagneticSite& s_j = GetMagneticSite(j);

			// get the pre-calculated u vectors
			const t_vec& u_i = s_i.ge_trafo_plane_calc;
			const t_vec& u_j = s_j.ge_trafo_plane_calc;
			const t_vec& uc_i = s_i.ge_trafo_plane_conj_calc;
			const t_vec& uc_j = s_j.ge_trafo_plane_conj_calc;

			// pre-factors of equation (44) from (Toth 2015)
			const t_real S_mag = std::sqrt(s_i.spin_mag_calc * s_j.spin_mag_calc);
			const t_cplx phase = std::exp(-m_phase_sign * s_imag * s_twopi *
				tl2::inner<t_vec_real>(s_j.pos_calc - s_i.pos_calc, S.Q_rlu));

			// matrix elements of equation (44) from (Toth 2015)
			M(    i,     j) = phase * S_mag * u_i[x_idx]  * uc_j[y_idx];
			M(    i, N + j) = phase * S_mag * u_i[x_idx]  * u_j[y_idx];
			M(N + i,     j) = phase * S_mag * uc_i[x_idx] * uc_j[y_idx];
			M(N + i, N + j) = phase * S_mag * uc_i[x_idx] * u_j[y_idx];
		} // end of site iteration

		const t_mat M_trafo = trafo_herm * M * trafo;

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
		std::cout << "M_trafo for x=" << int(x_idx) << ", y=" << int(y_idx) << ":\n";
		tl2::niceprint(std::cout, M_trafo, 1e-4, 4);
#endif

		for(t_size i = 0; i < S.E_and_S.size(); ++i)
		{
			S.E_and_S[i].S(x_idx, y_idx) +=
				M_trafo(i, i) / t_real(M.size1());
		}
	} // end of coordinate iteration

	return true;
}



/**
 * applies projectors, form and weight factors to get neutron intensities
 * @note implements the formalism given by (Toth 2015)
 */
MAGDYN_TEMPL
void MAGDYN_INST::CalcIntensities(MAGDYN_TYPE::SofQE& S) const
{
	using namespace tl2_ops;
	tl2::ExprParser<t_cplx> magffact = m_magffact;

	for(EnergyAndWeight& E_and_S : S.E_and_S)
	{
		// apply bose factor
		if(m_temperature >= 0.)
			E_and_S.S *= tl2::bose_cutoff(E_and_S.E, m_temperature, m_bose_cutoff);

		// apply form factor
		if(m_magffact_formula != "")
		{
			// get |Q| in units of A^(-1)
			t_vec_real Q_invA = m_xtalB * S.Q_rlu;
			t_real Q_abs = tl2::norm<t_vec_real>(Q_invA);

			// evaluate form factor expression
			magffact.register_var("Q", Q_abs);
			t_real ffact = magffact.eval_noexcept().real();
			E_and_S.S *= ffact;
		}

		// apply orthogonal projector for magnetic neutron scattering,
		// see (Shirane 2002), p. 37, equation (2.64)
		t_mat proj_neutron = tl2::ortho_projector<t_mat, t_vec>(S.Q_rlu, false);
		E_and_S.S_perp = proj_neutron * E_and_S.S * proj_neutron;

		CalcPolarisation(S.Q_rlu, E_and_S);

		// weights
		E_and_S.S_sum       = tl2::trace<t_mat>(E_and_S.S);
		E_and_S.S_perp_sum  = tl2::trace<t_mat>(E_and_S.S_perp);
		E_and_S.weight_full = std::abs(E_and_S.S_sum.real());
		E_and_S.weight      = std::abs(E_and_S.S_perp_sum.real());
	}
}
// --------------------------------------------------------------------

#endif
