/**
 * tlibs2 -- magnetic dynamics
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

#ifndef __TLIBS2_MAGDYN_DISP_H__
#define __TLIBS2_MAGDYN_DISP_H__

#include <numeric>
#include <vector>
#include <thread>

#include <boost/asio.hpp>

#include "../maths.h"
#include "magdyn.h"


// --------------------------------------------------------------------
// calculation functions
// --------------------------------------------------------------------

/**
 * unite degenerate energies and their corresponding eigenstates
 * TODO: handle case where the energy is the same but the state is different
 */
MAGDYN_TEMPL
MAGDYN_TYPE::SofQE MAGDYN_INST::UniteEnergies(const MAGDYN_TYPE::SofQE& S) const
{
	SofQE new_S = S;
	new_S.E_and_S.clear();
	new_S.E_and_S.reserve(S.E_and_S.size());

	for(const EnergyAndWeight& curState : S.E_and_S)
	{
		// find states with the same energy
		auto iter = std::find_if(new_S.E_and_S.begin(), new_S.E_and_S.end(),
			[&curState, this](const EnergyAndWeight& E_and_S) -> bool
		{
			t_real E = E_and_S.E;
			return tl2::equals<t_real>(E, curState.E, m_eps);
		});

		if(iter == new_S.E_and_S.end())
		{
			// energy not yet seen
			new_S.E_and_S.push_back(curState);
		}
		else
		{
			// energy already seen: add correlation matrices and weights
			iter->S           += curState.S;
			iter->S_perp      += curState.S_perp;
			iter->S_sum       += curState.S_sum;
			iter->S_perp_sum  += curState.S_perp_sum;
			iter->weight      += curState.weight;
			iter->weight_full += curState.weight_full;
			iter->degeneracy  += curState.degeneracy;
		}
	}

	return new_S;
}



/**
 * get the energies and the spin-correlation at the given momentum
 * (also calculates incommensurate contributions and applies weight factors)
 * @note implements the formalism given by (Toth 2015)
 */
MAGDYN_TEMPL
MAGDYN_TYPE::SofQE
MAGDYN_INST::CalcEnergies(const t_vec_real& Q_rlu, bool only_energies) const
{
	auto calc_EandS = [only_energies, this](const t_vec_real& Q) -> SofQE
	{
		const t_mat H = CalcHamiltonian(Q);
		return CalcEnergiesFromHamiltonian(H, Q, only_energies);
	};

	SofQE S;
	S.Q_rlu = Q_rlu;
	if(m_calc_H)
		S = calc_EandS(Q_rlu);

	if(IsIncommensurate())
	{
		// equations (39) and (40) from (Toth 2015)
		const t_mat proj_norm = tl2::convert<t_mat>(
			tl2::projector<t_mat_real, t_vec_real>(m_rotaxis, true));

		t_mat rot_incomm = tl2::unit<t_mat>(3);
		rot_incomm -= s_imag * m_phase_sign * tl2::skewsymmetric<t_mat, t_vec>(m_rotaxis);
		rot_incomm -= proj_norm;
		rot_incomm *= 0.5;

		// calculate additional hamiltonians for Q+-O
		SofQE S_p{}, S_m{};
		if(m_calc_Hp)
			S_p = calc_EandS(Q_rlu + m_ordering);
		if(m_calc_Hm)
			S_m = calc_EandS(Q_rlu - m_ordering);

		// move over additional hamiltonians for Q+-O
		S.H_p = std::move(S_p.H);
		S.H_chol_p = std::move(S_p.H_chol);
		S.H_comm_p = std::move(S_p.H_comm);
		S.evec_mat_p = std::move(S_p.evec_mat);
		S.H_m = std::move(S_m.H);
		S.H_chol_m = std::move(S_m.H_chol);
		S.H_comm_m = std::move(S_m.H_comm);
		S.evec_mat_m = std::move(S_m.evec_mat);

		if(!only_energies)
		{
			const t_mat rot_incomm_conj = tl2::conj(rot_incomm);

			// formula 40 from (Toth 2015)
			for(EnergyAndWeight& EandW : S.E_and_S)
				EandW.S = EandW.S * proj_norm;
			for(EnergyAndWeight& EandW : S_p.E_and_S)
				EandW.S = EandW.S * rot_incomm;
			for(EnergyAndWeight& EandW : S_m.E_and_S)
				EandW.S = EandW.S * rot_incomm_conj;
		}

		// merge energies and weights
		S.E_and_S.reserve(S.E_and_S.size() + S_p.E_and_S.size() + S_m.E_and_S.size());
		for(EnergyAndWeight& EandW : S_p.E_and_S)
			S.E_and_S.emplace_back(std::move(EandW));
		for(EnergyAndWeight& EandW : S_m.E_and_S)
			S.E_and_S.emplace_back(std::move(EandW));
	}

	if(!only_energies)
		CalcIntensities(S);

	if(m_unite_degenerate_energies)
		S = UniteEnergies(S);

	if(!only_energies)
		CheckImagWeights(S);

	return S;
}



MAGDYN_TEMPL
MAGDYN_TYPE::EnergiesAndWeights
MAGDYN_INST::CalcEnergies(t_real h, t_real k, t_real l, bool only_energies) const
{
	// momentum transfer
	const t_vec_real Qvec = tl2::create<t_vec_real>({ h, k, l });
	return CalcEnergies(Qvec, only_energies).E_and_S;
}



/**
 * generates the dispersion plot along the given Q path
 */
MAGDYN_TEMPL
MAGDYN_TYPE::SofQEs
MAGDYN_INST::CalcDispersion(t_real h_start, t_real k_start, t_real l_start,
	t_real h_end, t_real k_end, t_real l_end, t_size num_Qs,
	t_size num_threads, std::function<bool(int, int)> *progress_fkt) const
{
	// determine number of threads
	if(num_threads == 0)
		num_threads = std::max<t_size>(1, std::thread::hardware_concurrency() / 2);

	// thread pool and tasks
	using t_pool = boost::asio::thread_pool;
	using t_task = std::packaged_task<SofQE()>;
	using t_taskptr = std::shared_ptr<t_task>;

	t_pool pool{num_threads};
	std::vector<t_taskptr> tasks;
	tasks.reserve(num_Qs);

	// calculate dispersion
	for(t_size i = 0; i < num_Qs; ++i)
	{
		if(progress_fkt && !(*progress_fkt)(0, num_Qs))
			break;

		auto task = [this, i, num_Qs,
			h_start, k_start, l_start,
			h_end, k_end, l_end]() -> SofQE
		{
			// get Q
			const t_real h = std::lerp(h_start, h_end, t_real(i) / t_real(num_Qs - 1));
			const t_real k = std::lerp(k_start, k_end, t_real(i) / t_real(num_Qs - 1));
			const t_real l = std::lerp(l_start, l_end, t_real(i) / t_real(num_Qs - 1));
			const t_vec_real Q = tl2::create<t_vec_real>({ h, k, l });

			// get E and S(Q, E) for this Q
			return CalcEnergies(Q, false);
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		boost::asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	// collect results
	SofQEs results;
	results.reserve(tasks.size());

	t_size Qs_finished = 0;
	for(auto& task : tasks)
	{
		if(progress_fkt && !(*progress_fkt)(Qs_finished + 1, num_Qs))
			break;

		const SofQE& result = task->get_future().get();
		results.push_back(result);
		++Qs_finished;
	}

	return results;
}



/**
 * generates the dispersion plot along the given 2d Q surface
 */
MAGDYN_TEMPL
MAGDYN_TYPE::SofQEs
MAGDYN_INST::CalcDispersion(t_real h_start, t_real k_start, t_real l_start,
	t_real h_end1, t_real k_end1, t_real l_end1,
	t_real h_end2, t_real k_end2, t_real l_end2,
	t_size num_Qs_sqrt, t_size num_threads,
	std::function<bool(int, int)> *progress_fkt) const
{
	// determine number of threads
	if(num_threads == 0)
		num_threads = std::max<t_size>(1, std::thread::hardware_concurrency() / 2);

	// thread pool and tasks
	using t_pool = boost::asio::thread_pool;
	using t_task = std::packaged_task<SofQE()>;
	using t_taskptr = std::shared_ptr<t_task>;

	t_pool pool{num_threads};
	std::vector<t_taskptr> tasks;
	tasks.reserve(num_Qs_sqrt * num_Qs_sqrt);

	// calculate dispersion
	for(t_size i = 0; i < num_Qs_sqrt; ++i)
	for(t_size j = 0; j < num_Qs_sqrt; ++j)
	{
		if(progress_fkt && !(*progress_fkt)(0, num_Qs_sqrt * num_Qs_sqrt))
			break;

		auto task = [this, i, j, num_Qs_sqrt,
			h_start, k_start, l_start,
			h_end1, k_end1, l_end1,
			h_end2, k_end2, l_end2]() -> SofQE
		{
			// get Q
			const t_real h =
				std::lerp(h_start, h_end1, t_real(i) / t_real(num_Qs_sqrt - 1)) +
				std::lerp(h_start, h_end2, t_real(j) / t_real(num_Qs_sqrt - 1)) - h_start;
			const t_real k =
				std::lerp(k_start, k_end1, t_real(i) / t_real(num_Qs_sqrt - 1)) +
				std::lerp(k_start, k_end2, t_real(j) / t_real(num_Qs_sqrt - 1)) - k_start;
			const t_real l =
				std::lerp(l_start, l_end1, t_real(i) / t_real(num_Qs_sqrt - 1)) +
				std::lerp(l_start, l_end2, t_real(j) / t_real(num_Qs_sqrt - 1)) - l_start;
			const t_vec_real Q = tl2::create<t_vec_real>({ h, k, l });

			// get E and S(Q, E) for this Q
			return CalcEnergies(Q, false);
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		boost::asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	// collect results
	SofQEs results;
	results.reserve(tasks.size());

	t_size Qs_finished = 0;
	for(auto& task : tasks)
	{
		if(progress_fkt && !(*progress_fkt)(Qs_finished + 1, num_Qs_sqrt * num_Qs_sqrt))
			break;

		const SofQE& result = task->get_future().get();
		results.push_back(result);
		++Qs_finished;
	}

	return results;
}
// --------------------------------------------------------------------

#endif
