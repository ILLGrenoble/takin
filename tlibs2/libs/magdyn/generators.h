/**
 * tlibs2 -- magnetic dynamics -- generators
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
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

#ifndef __TLIBS2_MAGDYN_GEN_H__
#define __TLIBS2_MAGDYN_GEN_H__

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "../maths.h"
#include "../expr.h"

#include "magdyn.h"


// --------------------------------------------------------------------
// symmetrisation and generation functions
// --------------------------------------------------------------------
/**
 * generate symmetric positions based on the given symops
 */
MAGDYN_TEMPL
void MAGDYN_INST::SymmetriseMagneticSites(const std::vector<t_mat_real>& symops)
{
	CalcExternalField();
	CalcMagneticSites();

	MagneticSites newsites;
	newsites.reserve(GetMagneticSitesCount() * symops.size());

	for(const MagneticSite& site : GetMagneticSites())
	{
		// get symmetry-equivalent positions
		const auto positions = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			site.pos_calc, symops, m_eps);

		for(t_size idx = 0; idx < positions.size(); ++idx)
		{
			MagneticSite newsite = site;
			newsite.pos_calc = positions[idx];
			newsite.pos[0] = tl2::var_to_str(newsite.pos_calc[0], m_prec);
			newsite.pos[1] = tl2::var_to_str(newsite.pos_calc[1], m_prec);
			newsite.pos[2] = tl2::var_to_str(newsite.pos_calc[2], m_prec);
			newsite.name += "_" + tl2::var_to_str(idx + 1, m_prec);

			newsites.emplace_back(std::move(newsite));
		}
	}

	m_sites = std::move(newsites);
	RemoveDuplicateMagneticSites();
	CalcSymmetryIndices(symops);
	CalcMagneticSites();
}



/**
 * generate symmetric exchange terms based on the given symops
 */
MAGDYN_TEMPL
void MAGDYN_INST::SymmetriseExchangeTerms(const std::vector<t_mat_real>& symops)
{
	CalcExternalField();
	CalcMagneticSites();
	CalcExchangeTerms();

	ExchangeTerms newterms;
	newterms.reserve(GetExchangeTermsCount() * symops.size());

	tl2::ExprParser<t_cplx> parser = GetExprParser();

	// create unit cell site vectors
	std::vector<t_vec_real> sites_uc = GetMagneticSitePositions(true);

	for(const ExchangeTerm& term : GetExchangeTerms())
	{
		// check if the site indices are valid
		if(!CheckMagneticSite(term.site1_calc) || !CheckMagneticSite(term.site2_calc))
			continue;

		// super cell distance vector
		t_vec_real dist_sc = to_4vec<t_vec_real>(term.dist_calc, 0.);

		// generate new (possibly supercell) sites with symop
		auto sites1_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			sites_uc[term.site1_calc], symops, m_eps,
			false /*keep in uc*/, true /*ignore occupied*/,
			true /*return homogeneous*/);
		auto sites2_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			sites_uc[term.site2_calc] + dist_sc, symops, m_eps,
			false /*keep in uc*/, true /*ignore occupied*/,
			true /*return homogeneous*/);

		// generate new dmi vectors
		t_vec_real dmi = tl2::zero<t_vec_real>(4);

		for(std::uint8_t dmi_idx = 0; dmi_idx < 3; ++dmi_idx)
		{
			if(term.dmi[dmi_idx] == "")
				continue;

			if(parser.parse_noexcept(term.dmi[dmi_idx]))
			{
				dmi[dmi_idx] = parser.eval_noexcept().real();
			}
			else
			{
				CERR_OPT << "Magdyn error: Parsing DMI component "
					<< int(dmi_idx) << " of term \"" << term.name
					<< "\"." << std::endl;
			}
		}

		const auto newdmis = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
			dmi, symops, m_eps,
			false /*keep in uc*/, true /*ignore occupied*/,
			false /*return homogeneous*/, true /*pseudovector*/);

		// generate new general J matrices
		t_real Jgen_arr[3][3]{};

		for(std::uint8_t J_idx1 = 0; J_idx1 < 3; ++J_idx1)
		{
			for(std::uint8_t J_idx2 = 0; J_idx2 < 3; ++J_idx2)
			{
				if(term.Jgen[J_idx1][J_idx2] == "")
					continue;

				if(parser.parse_noexcept(term.Jgen[J_idx1][J_idx2]))
				{
					Jgen_arr[J_idx1][J_idx2] = parser.eval_noexcept().real();
				}
				else
				{
					CERR_OPT << "Magdyn error: Parsing general J component ("
						<< int(J_idx1) << ", " << int(J_idx2)
						<< ") of term \"" << term.name << "\"."
						<< std::endl;
				}
			}
		}

		t_mat_real Jgen = tl2::create<t_mat_real>({
			Jgen_arr[0][0], Jgen_arr[0][1], Jgen_arr[0][2], 0,
			Jgen_arr[1][0], Jgen_arr[1][1], Jgen_arr[1][2], 0,
			Jgen_arr[2][0], Jgen_arr[2][1], Jgen_arr[2][2], 0,
			0,              0,              0,              0 });

		const auto newJgens = tl2::apply_ops_hom<t_mat_real, t_real>(Jgen, symops);

		// iterate and insert generated couplings
		for(t_size op_idx = 0; op_idx < std::min(sites1_sc.size(), sites2_sc.size()); ++op_idx)
		{
			// get position of the site in the supercell
			const auto [sc1_ok, site1_sc_idx, sc1] = tl2::get_supercell(
				sites1_sc[op_idx], sites_uc, 3, m_eps);
			const auto [sc2_ok, site2_sc_idx, sc2] = tl2::get_supercell(
				sites2_sc[op_idx], sites_uc, 3, m_eps);

			if(!sc1_ok || !sc2_ok)
			{
				CERR_OPT << "Magdyn error: Could not find supercell"
					<< " for position generated from symop "
					<< op_idx << "." << std::endl;
			}

			ExchangeTerm newterm = term;
			newterm.site1_calc = site1_sc_idx;
			newterm.site2_calc = site2_sc_idx;
			newterm.site1 = GetMagneticSite(newterm.site1_calc).name;
			newterm.site2 = GetMagneticSite(newterm.site2_calc).name;
			newterm.dist_calc = to_3vec<t_vec_real>(sc2 - sc1);
			newterm.dist[0] = tl2::var_to_str(newterm.dist_calc[0], m_prec);
			newterm.dist[1] = tl2::var_to_str(newterm.dist_calc[1], m_prec);
			newterm.dist[2] = tl2::var_to_str(newterm.dist_calc[2], m_prec);

			for(std::uint8_t idx1 = 0; idx1 < 3; ++idx1)
			{
				newterm.dmi[idx1] = tl2::var_to_str(newdmis[op_idx][idx1], m_prec);
				for(std::uint8_t idx2 = 0; idx2 < 3; ++idx2)
					newterm.Jgen[idx1][idx2] = tl2::var_to_str(
						newJgens[op_idx](idx1, idx2), m_prec);
			}
			newterm.name += "_" + tl2::var_to_str(op_idx + 1, m_prec);

			newterms.emplace_back(std::move(newterm));
		}
	}

	m_exchange_terms = std::move(newterms);
	RemoveDuplicateExchangeTerms();
	CalcSymmetryIndices(symops);
	CalcExchangeTerms();
}



/**
 * generate possible couplings up to a certain distance
 */
MAGDYN_TEMPL
void MAGDYN_INST::GeneratePossibleExchangeTerms(
	t_real dist_max, t_size _sc_max,
	std::optional<t_size> couplings_max)
{
	if(GetMagneticSitesCount() == 0)
		return;

	// helper struct for finding possible couplings
	struct PossibleCoupling
	{
		// corresponding unit cell position
		t_vec_real pos1_uc{};
		t_vec_real pos2_uc{};

		// magnetic site position in supercell
		t_vec_real sc_vec{};
		t_vec_real pos2_sc{};

		// coordinates in orthogonal lab units
		t_vec_real pos1_uc_lab{};
		t_vec_real pos2_sc_lab{};

		// corresponding unit cell index
		t_size idx1_uc{};
		t_size idx2_uc{};

		// distance between the two magnetic sites
		t_real dist{};
	};

	CalcExternalField();
	CalcMagneticSites();
	CalcExchangeTerms();

	ExchangeTerms newterms;
	std::vector<PossibleCoupling> couplings;

	// generate a list of supercell vectors
	t_real sc_max = t_real(_sc_max);
	for(t_real sc_h = -sc_max; sc_h <= sc_max; sc_h += 1.)
	for(t_real sc_k = -sc_max; sc_k <= sc_max; sc_k += 1.)
	for(t_real sc_l = -sc_max; sc_l <= sc_max; sc_l += 1.)
	{
		// super cell vector
		const t_vec_real sc_vec = tl2::create<t_vec_real>({ sc_h, sc_k, sc_l });

		for(t_size idx1 = 0; idx1 < GetMagneticSitesCount() - 1; ++idx1)
		for(t_size idx2 = idx1 + 1; idx2 < GetMagneticSitesCount(); ++idx2)
		{
			PossibleCoupling coupling;

			coupling.idx1_uc = idx1;
			coupling.idx2_uc = idx2;

			coupling.pos1_uc = GetMagneticSite(idx1).pos_calc;
			coupling.pos2_uc = GetMagneticSite(idx2).pos_calc;

			coupling.sc_vec = sc_vec;
			coupling.pos2_sc = coupling.pos2_uc + sc_vec;

			// transform to lab units for correct distance
			coupling.pos1_uc_lab = m_xtalA * coupling.pos1_uc;
			coupling.pos2_sc_lab = m_xtalA * coupling.pos2_sc;

			coupling.dist = tl2::norm<t_vec_real>(
				coupling.pos2_sc_lab - coupling.pos1_uc_lab);
			if(coupling.dist <= dist_max && coupling.dist > m_eps)
				couplings.emplace_back(std::move(coupling));
		}
	}

	// sort couplings by distance
	std::stable_sort(couplings.begin(), couplings.end(),
		[](const PossibleCoupling& coupling1, const PossibleCoupling& coupling2) -> bool
	{
		return coupling1.dist < coupling2.dist;
	});

	// add couplings to list
	t_size coupling_idx = 0;
	for(const PossibleCoupling& coupling : couplings)
	{
		if(couplings_max && coupling_idx >= *couplings_max)
			break;

		ExchangeTerm newterm{};
		newterm.site1_calc = coupling.idx1_uc;
		newterm.site2_calc = coupling.idx2_uc;
		newterm.site1 = GetMagneticSite(newterm.site1_calc).name;
		newterm.site2 = GetMagneticSite(newterm.site2_calc).name;
		newterm.dist_calc = coupling.sc_vec;
		newterm.dist[0] = tl2::var_to_str(newterm.dist_calc[0], m_prec);
		newterm.dist[1] = tl2::var_to_str(newterm.dist_calc[1], m_prec);
		newterm.dist[2] = tl2::var_to_str(newterm.dist_calc[2], m_prec);
		newterm.length_calc = coupling.dist;
		newterm.J = "0";
		newterm.name += "coupling_" + tl2::var_to_str(coupling_idx + 1, m_prec);
		newterms.emplace_back(std::move(newterm));

		++coupling_idx;
	}

	m_exchange_terms = std::move(newterms);
	RemoveDuplicateExchangeTerms();
	//CalcSymmetryIndices(symops);
	CalcExchangeTerms();
}



/**
 * extend the magnetic structure
 */
MAGDYN_TEMPL
void MAGDYN_INST::ExtendStructure(t_size x_size, t_size y_size, t_size z_size,
	bool remove_duplicates, bool flip_spin)
{
	CalcExternalField();
	CalcMagneticSites();

	const t_size num_sites = GetMagneticSitesCount();
	const t_size num_terms = GetExchangeTermsCount();
	m_sites.reserve(num_sites * x_size * y_size * z_size);
	m_exchange_terms.reserve(num_terms * x_size * y_size * z_size);

	// iterate over extended structure
	for(t_size x_idx = 0; x_idx < x_size; ++x_idx)
	for(t_size y_idx = 0; y_idx < y_size; ++y_idx)
	for(t_size z_idx = 0; z_idx < z_size; ++z_idx)
	{
		// ignore sites in original cell
		if(x_idx == 0 && y_idx == 0 && z_idx == 0)
			continue;

		std::string ext_id = "_" + tl2::var_to_str(x_idx + 1, m_prec) +
			"_" + tl2::var_to_str(y_idx + 1, m_prec) +
			"_" + tl2::var_to_str(z_idx + 1, m_prec);

		// duplicate sites
		for(t_size site_idx = 0; site_idx < num_sites; ++site_idx)
		{
			MagneticSite newsite = GetMagneticSite(site_idx);

			newsite.name += ext_id;
			newsite.pos_calc += tl2::create<t_vec_real>({
				t_real(x_idx), t_real(y_idx), t_real(z_idx) });
			newsite.pos[0] = tl2::var_to_str(newsite.pos_calc[0], m_prec);
			newsite.pos[1] = tl2::var_to_str(newsite.pos_calc[1], m_prec);
			newsite.pos[2] = tl2::var_to_str(newsite.pos_calc[2], m_prec);

			if(flip_spin && (x_idx + y_idx + z_idx) % 2 != 0)
			{
				// flip spin
				newsite.spin_dir_calc = -newsite.spin_dir_calc;
				newsite.spin_dir[0] = tl2::var_to_str(newsite.spin_dir_calc[0], m_prec);
				newsite.spin_dir[1] = tl2::var_to_str(newsite.spin_dir_calc[1], m_prec);
				newsite.spin_dir[2] = tl2::var_to_str(newsite.spin_dir_calc[2], m_prec);
			}

			AddMagneticSite(std::move(newsite));
		}

		// duplicate couplings
		for(t_size term_idx = 0; term_idx < num_terms; ++term_idx)
		{
			ExchangeTerm newterm = GetExchangeTerm(term_idx);

			newterm.site1 += ext_id;
			newterm.site2 += ext_id;

			AddExchangeTerm(std::move(newterm));
		}
	}

	if(remove_duplicates)
	{
		RemoveDuplicateMagneticSites();
		RemoveDuplicateExchangeTerms();
	}

	FixExchangeTerms(x_size, y_size, z_size);
	CalcMagneticSites();
	CalcExchangeTerms();
}



/**
 * modify exchange term whose sites point to sc positions that are also available in the uc
 */
MAGDYN_TEMPL
void MAGDYN_INST::FixExchangeTerms(t_size x_size, t_size y_size, t_size z_size)
{
	for(ExchangeTerm& term : GetExchangeTerms())
	{
		// coupling within uc?
		if(tl2::equals_0<t_vec_real>(term.dist_calc, m_eps))
			continue;

		// find site 2
		const MagneticSite* site2_uc = FindMagneticSite(term.site2);
		if(!site2_uc)
			continue;

		// get site 2 sc vector
		t_vec_real site2_sc = site2_uc->pos_calc + term.dist_calc;

		// fix couplings that are now internal:
		// see if site 2's sc vector is now also available in the uc
		bool fixed_coupling = false;
		for(const MagneticSite& site : GetMagneticSites())
		{
			if(tl2::equals<t_vec_real>(site.pos_calc, site2_sc, m_eps))
			{
				// found the identical site
				term.site2 = site.name;
				term.dist[0] = term.dist[1] = term.dist[2] = "0";
				term.dist_calc = tl2::zero<t_vec_real>(3);
				fixed_coupling = true;
				break;
			}
		}
		if(fixed_coupling)
			continue;

		// fix external couplings
		t_vec_real site2_newsc = site2_sc;

		site2_newsc[0] = std::fmod(site2_newsc[0], t_real(x_size));
		site2_newsc[1] = std::fmod(site2_newsc[1], t_real(y_size));
		site2_newsc[2] = std::fmod(site2_newsc[2], t_real(z_size));

		if(site2_newsc[0] < t_real(0))
			site2_newsc[0] += t_real(x_size);
		if(site2_newsc[1] < t_real(0))
			site2_newsc[1] += t_real(y_size);
		if(site2_newsc[2] < t_real(0))
			site2_newsc[2] += t_real(z_size);

		for(const MagneticSite& site : GetMagneticSites())
		{
			if(tl2::equals<t_vec_real>(site.pos_calc, site2_newsc, m_eps))
			{
				// found the identical site
				term.site2 = site.name;
				term.dist_calc = site2_sc - site2_newsc;
				term.dist[0] = tl2::var_to_str(term.dist_calc[0], m_prec);
				term.dist[1] = tl2::var_to_str(term.dist_calc[1], m_prec);
				term.dist[2] = tl2::var_to_str(term.dist_calc[2], m_prec);
				break;
			}
		}
	}
}



/**
 * remove literal duplicate sites (not symmetry-equivalent ones)
 */
MAGDYN_TEMPL
void MAGDYN_INST::RemoveDuplicateMagneticSites()
{
	for(auto iter1 = m_sites.begin(); iter1 != m_sites.end(); std::advance(iter1, 1))
	for(auto iter2 = std::next(iter1, 1); iter2 != m_sites.end();)
	{
		if(tl2::equals<t_vec_real>(iter1->pos_calc, iter2->pos_calc, m_eps))
			iter2 = m_sites.erase(iter2);
		else
			std::advance(iter2, 1);
	}
}



/**
 * remove literal duplicate couplings (not symmetry-equivalent ones)
 */
MAGDYN_TEMPL
void MAGDYN_INST::RemoveDuplicateExchangeTerms()
{
	for(auto iter1 = m_exchange_terms.begin(); iter1 != m_exchange_terms.end(); std::advance(iter1, 1))
	for(auto iter2 = std::next(iter1, 1); iter2 != m_exchange_terms.end();)
	{
		// identical coupling
		bool same_uc = (iter1->site1 == iter2->site1 && iter1->site2 == iter2->site2);
		bool same_sc = tl2::equals<t_vec_real>(iter1->dist_calc, iter2->dist_calc, m_eps);

		// flipped coupling
		bool inv_uc = (iter1->site1 == iter2->site2 && iter1->site2 == iter2->site1);
		bool inv_sc = tl2::equals<t_vec_real>(iter1->dist_calc, -iter2->dist_calc, m_eps);

		if((same_uc && same_sc) || (inv_uc && inv_sc))
			iter2 = m_exchange_terms.erase(iter2);
		else
			std::advance(iter2, 1);
	}
}



/**
 * are two sites equivalent with respect to the given symmetry operators?
 */
MAGDYN_TEMPL
bool MAGDYN_INST::IsSymmetryEquivalent(
	const MAGDYN_TYPE::MagneticSite& site1, const MAGDYN_TYPE::MagneticSite& site2,
	const std::vector<t_mat_real>& symops) const
{
	// get symmetry-equivalent positions
	const auto positions = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
		site1.pos_calc, symops, m_eps);

	for(const auto& pos : positions)
	{
		// symmetry-equivalent site found?
		if(tl2::equals<t_vec_real>(site2.pos_calc, pos, m_eps))
			return true;
	}

	return false;
}



/**
 * are two couplings equivalent with respect to the given symmetry operators?
 */
MAGDYN_TEMPL
bool MAGDYN_INST::IsSymmetryEquivalent(
	const MAGDYN_TYPE::ExchangeTerm& term1, const MAGDYN_TYPE::ExchangeTerm& term2,
	const std::vector<t_mat_real>& symops) const
{
	// check if site indices are within bounds
	const t_size N = GetMagneticSitesCount();
	if(term1.site1_calc >= N || term1.site2_calc >= N ||
		term2.site1_calc >= N || term2.site2_calc >= N)
		return false;

	// create unit cell site vectors
	std::vector<t_vec_real> sites_uc = GetMagneticSitePositions(true);

	// super cell distance vector
	t_vec_real dist_sc = to_4vec<t_vec_real>(term1.dist_calc, 0.);

	// generate new (possibly supercell) sites with symop
	auto sites1_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
		sites_uc[term1.site1_calc], symops, m_eps,
		false /*keep in uc*/, true /*ignore occupied*/,
		true /*return homogeneous*/);
	auto sites2_sc = tl2::apply_ops_hom<t_vec_real, t_mat_real, t_real>(
		sites_uc[term1.site2_calc] + dist_sc, symops, m_eps,
		false /*keep in uc*/, true /*ignore occupied*/,
		true /*return homogeneous*/);

	for(t_size idx = 0; idx < std::min(sites1_sc.size(), sites2_sc.size()); ++idx)
	{
		// get position of the site in the supercell
		const auto [sc1_ok, site1_sc_idx, sc1] = tl2::get_supercell(
			sites1_sc[idx], sites_uc, 3, m_eps);
		const auto [sc2_ok, site2_sc_idx, sc2] = tl2::get_supercell(
			sites2_sc[idx], sites_uc, 3, m_eps);
		if(!sc1_ok || !sc2_ok)
			continue;

		// symmetry-equivalent coupling found?
		if(tl2::equals<t_vec_real>(to_3vec<t_vec_real>(sc2 - sc1), term2.dist_calc, m_eps)
			&& site1_sc_idx == term2.site1_calc && site2_sc_idx == term2.site2_calc)
			return true;
	}

	return false;
}



/**
 * assign symmetry group indices to sites and couplings
 */
MAGDYN_TEMPL
void MAGDYN_INST::CalcSymmetryIndices(const std::vector<t_mat_real>& symops)
{
	// iterate sites
	t_size site_sym_idx_ctr = 0;
	std::vector<const MagneticSite*> seen_sites;

	for(MagneticSite& site : GetMagneticSites())
	{
		const MagneticSite* equivalent_site = nullptr;
		for(const MagneticSite* seen_site : seen_sites)
		{
			// found symmetry-equivalent site
			if(IsSymmetryEquivalent(site, *seen_site, symops))
			{
				equivalent_site = seen_site;
				break;
			}
		}

		if(equivalent_site)
		{
			// assign symmetry index from equivalent site
			site.sym_idx = equivalent_site->sym_idx;
		}
		else
		{
			// assign new symmetry index
			site.sym_idx = ++site_sym_idx_ctr;
			seen_sites.push_back(&site);
		}
	}


	// iterate couplings
	std::vector<const ExchangeTerm*> seen_terms;
	t_size term_sym_idx_ctr = 0;

	for(ExchangeTerm& term : GetExchangeTerms())
	{
		const ExchangeTerm* equivalent_term = nullptr;
		for(const ExchangeTerm* seen_term : seen_terms)
		{
			// found symmetry-equivalent coupling
			if(IsSymmetryEquivalent(term, *seen_term, symops))
			{
				equivalent_term = seen_term;
				break;
			}
		}

		if(equivalent_term)
		{
			// assign symmetry index from equivalent coupling
			term.sym_idx = equivalent_term->sym_idx;
		}
		else
		{
			// assign new symmetry index
			term.sym_idx = ++term_sym_idx_ctr;
			seen_terms.push_back(&term);
		}
	}
}
// --------------------------------------------------------------------

#endif
