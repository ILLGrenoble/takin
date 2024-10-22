/**
 * tlibs2 -- magnetic dynamics -- getters / setters / cleanup functions
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

#ifndef __TLIBS2_MAGDYN_GETTERS_H__
#define __TLIBS2_MAGDYN_GETTERS_H__

#include <vector>
#include <tuple>
#include <string>
#include <iostream>
#include <iomanip>

#include "../maths.h"
#include "../expr.h"



// --------------------------------------------------------------------
// cleanup functions
// --------------------------------------------------------------------
/**
 * clear all
 */
MAGDYN_TEMPL void MAGDYN_INST::Clear()
{
	ClearVariables();
	ClearMagneticSites();
	ClearExchangeTerms();
	ClearExternalField();

	// clear temperature, -1: don't use
	m_temperature = -1.;

	// clear form factor
	m_magffact_formula = "";

	// clear ordering wave vector
	m_ordering = tl2::zero<t_vec_real>(3);
	m_is_incommensurate = false;

	// reset rotation axis
	m_rotaxis = tl2::create<t_vec_real>({ 1., 0., 0. });

	// clear crystal
	m_xtallattice[0] = m_xtallattice[1] = m_xtallattice[2] = 0.;
	m_xtalangles[0] = m_xtalangles[1] = m_xtalangles[2] = t_real(0.5) * tl2::pi<t_real>;
	m_xtalA = m_xtalB = tl2::unit<t_mat_real>(3);
	m_xtalUB = m_xtalUBinv = tl2::unit<t_mat_real>(3);

	// clear scattering plane
	m_scatteringplane[0] = tl2::create<t_vec_real>({ 1., 0., 0. });
	m_scatteringplane[1] = tl2::create<t_vec_real>({ 0., 1., 0. });
	m_scatteringplane[2] = tl2::create<t_vec_real>({ 0., 0., 1. });
}



/**
 * clear all parser variables
 */
MAGDYN_TEMPL void MAGDYN_INST::ClearVariables()
{
	m_variables.clear();
}



/**
 * clear all magnetic sites
 */
MAGDYN_TEMPL void MAGDYN_INST::ClearMagneticSites()
{
	m_sites.clear();
}



/**
 * clear all couplings
 */
MAGDYN_TEMPL void MAGDYN_INST::ClearExchangeTerms()
{
	m_exchange_terms.clear();
}



/**
 * clear the external field settings
 */
MAGDYN_TEMPL void MAGDYN_INST::ClearExternalField()
{
	m_field.dir.clear();
	m_field.mag = 0.;
	m_field.align_spins = false;
}
// --------------------------------------------------------------------



// --------------------------------------------------------------------
// getter
// --------------------------------------------------------------------
MAGDYN_TEMPL const MAGDYN_TYPE::Variables& MAGDYN_INST::GetVariables() const
{
	return m_variables;
}


MAGDYN_TEMPL const MAGDYN_TYPE::MagneticSites& MAGDYN_INST::GetMagneticSites() const
{
	return m_sites;
}


MAGDYN_TEMPL MAGDYN_TYPE::MagneticSites& MAGDYN_INST::GetMagneticSites()
{
	return m_sites;
}


MAGDYN_TEMPL t_size MAGDYN_INST::GetMagneticSitesCount() const
{
	return m_sites.size();
}


MAGDYN_TEMPL const MAGDYN_TYPE::ExchangeTerms& MAGDYN_INST::GetExchangeTerms() const
{
	return m_exchange_terms;
}


MAGDYN_TEMPL MAGDYN_TYPE::ExchangeTerms& MAGDYN_INST::GetExchangeTerms()
{
	return m_exchange_terms;
}


MAGDYN_TEMPL t_size MAGDYN_INST::GetExchangeTermsCount() const
{
	return m_exchange_terms.size();
}


MAGDYN_TEMPL const MAGDYN_TYPE::ExternalField& MAGDYN_INST::GetExternalField() const
{
	return m_field;
}


MAGDYN_TEMPL const t_vec_real& MAGDYN_INST::GetRotationAxis() const
{
	return m_rotaxis;
}


MAGDYN_TEMPL const t_vec_real& MAGDYN_INST::GetOrderingWavevector() const
{
	return m_ordering;
}


MAGDYN_TEMPL t_real MAGDYN_INST::GetTemperature() const
{
	return m_temperature;
}


MAGDYN_TEMPL t_real MAGDYN_INST::GetBoseCutoffEnergy() const
{
	return m_bose_cutoff;
}


MAGDYN_TEMPL const std::string& MAGDYN_INST::GetMagneticFormFactor() const
{
	return m_magffact_formula;
}


MAGDYN_TEMPL const t_mat_real& MAGDYN_INST::GetCrystalATrafo() const
{
	return m_xtalA;
}


MAGDYN_TEMPL const t_mat_real& MAGDYN_INST::GetCrystalBTrafo() const
{
	return m_xtalB;
}


MAGDYN_TEMPL const t_mat_real& MAGDYN_INST::GetCrystalUBTrafo() const
{
	return m_xtalUB;
}



MAGDYN_TEMPL const MAGDYN_TYPE::MagneticSite&
MAGDYN_INST::GetMagneticSite(t_size idx) const
{
	if(!CheckMagneticSite(idx))
	{
		static MagneticSite null_site{};
		return null_site;
	}

	return m_sites[idx];
}



MAGDYN_TEMPL const MAGDYN_TYPE::ExchangeTerm&
MAGDYN_INST::GetExchangeTerm(t_size idx) const
{
	if(!CheckExchangeTerm(idx))
	{
		static ExchangeTerm null_term{};
		return null_term;
	}

	return m_exchange_terms[idx];
}



MAGDYN_TEMPL bool MAGDYN_INST::IsIncommensurate() const
{
	return m_is_incommensurate || m_force_incommensurate;
}



MAGDYN_TEMPL bool MAGDYN_INST::GetSilent() const
{
	return m_silent;
}



/**
 * get number of magnetic sites with the given name (to check if the name is unique)
 */
MAGDYN_TEMPL std::vector<const MAGDYN_TYPE::MagneticSite*>
MAGDYN_INST::FindMagneticSites(const std::string& name) const
{
	std::vector<const MagneticSite*> sites;

	for(const MagneticSite& site : GetMagneticSites())
	{
		if(site.name == name)
			sites.push_back(&site);
	}

	return sites;
}



/**
 * get magnetic site with the given name
 */
MAGDYN_TEMPL const MAGDYN_TYPE::MagneticSite*
MAGDYN_INST::FindMagneticSite(const std::string& name) const
{
	for(const MagneticSite& site : GetMagneticSites())
	{
		if(site.name == name)
			return &site;
	}

	return nullptr;
}



/**
 * get the index of a magnetic site from its name
 */
MAGDYN_TEMPL t_size MAGDYN_INST::GetMagneticSiteIndex(const std::string& name) const
{
	// try to find the site index by name
	for(t_size idx = 0; idx < GetMagneticSitesCount(); ++idx)
	{
		if(GetMagneticSite(idx).name == name)
			return idx;
	}

	// alternatively try to parse the expression for the index
	tl2::ExprParser<t_size> parser;
	parser.SetInvalid0(false);
	parser.SetAutoregisterVariables(false);
	if(parser.parse_noexcept(name))
	{
		if(t_size idx = parser.eval_noexcept(); idx < GetMagneticSitesCount())
			return idx;
	}
	else
	{
		CERR_OPT << "Magdyn error: "
			<< "Invalid site name \"" << name << "\"."
			<< std::endl;
	}

	// nothing found: return invalid index
	return GetMagneticSitesCount();
}



/**
 * get the index of an exchange term from its name
 */
MAGDYN_TEMPL t_size MAGDYN_INST::GetExchangeTermIndex(const std::string& name) const
{
	// try to find the term index by name
	for(t_size idx = 0; idx < GetExchangeTermsCount(); ++idx)
	{
		if(GetExchangeTerm(idx).name == name)
			return idx;
	}

	// alternatively try to parse the expression for the index
	tl2::ExprParser<t_size> parser;
	parser.SetInvalid0(false);
	parser.SetAutoregisterVariables(false);
	if(parser.parse_noexcept(name))
	{
		if(t_size idx = parser.eval_noexcept(); idx < GetExchangeTermsCount())
			return idx;
	}
	else
	{
		CERR_OPT << "Magdyn error: Invalid coupling name \"" << name << "\"."
			<< std::endl;
	}

	// nothing found: return invalid index
	return GetExchangeTermsCount();
}



MAGDYN_TEMPL std::vector<t_vec_real>
MAGDYN_INST::GetMagneticSitePositions(bool homogeneous) const
{
	std::vector<t_vec_real> sites;
	sites.reserve(GetMagneticSitesCount());

	for(const MagneticSite& site : GetMagneticSites())
	{
		if(homogeneous)
			sites.push_back(to_4vec<t_vec_real>(site.pos_calc, 1.));
		else
			sites.push_back(to_3vec<t_vec_real>(site.pos_calc));
	}

	return sites;
}



MAGDYN_TEMPL t_vec_real MAGDYN_INST::GetCrystalLattice() const
{
	return tl2::create<t_vec_real>(
	{
		m_xtallattice[0], m_xtallattice[1], m_xtallattice[2],
		m_xtalangles[0], m_xtalangles[1], m_xtalangles[2],
	});
}



MAGDYN_TEMPL const t_vec_real* MAGDYN_INST::GetScatteringPlane() const
{
	return m_scatteringplane;
}



/**
 * get the needed supercell ranges from the exchange terms
 */
MAGDYN_TEMPL std::tuple<t_vec_real, t_vec_real>
MAGDYN_INST::GetSupercellMinMax() const
{
	t_vec_real min = tl2::zero<t_vec_real>(3);
	t_vec_real max = tl2::zero<t_vec_real>(3);

	for(const ExchangeTerm& term : GetExchangeTerms())
	{
		for(std::uint8_t i = 0; i < 3; ++i)
		{
			min[i] = std::min(min[i], term.dist_calc[i]);
			max[i] = std::max(max[i], term.dist_calc[i]);
		}
	}

	return std::make_tuple(min, max);
}
// --------------------------------------------------------------------



// --------------------------------------------------------------------
// setter
// --------------------------------------------------------------------
MAGDYN_TEMPL void MAGDYN_INST::SetEpsilon(t_real eps)
{
	m_eps = eps;
}


MAGDYN_TEMPL void MAGDYN_INST::SetPrecision(int prec)
{
	m_prec = prec;
}


MAGDYN_TEMPL void MAGDYN_INST::SetTemperature(t_real T)
{
	m_temperature = T;
}


MAGDYN_TEMPL void MAGDYN_INST::SetBoseCutoffEnergy(t_real E)
{
	m_bose_cutoff = E;
}


MAGDYN_TEMPL void MAGDYN_INST::SetUniteDegenerateEnergies(bool b)
{
	m_unite_degenerate_energies = b;
}


MAGDYN_TEMPL void MAGDYN_INST::SetForceIncommensurate(bool b)
{
	m_force_incommensurate = b;
}


MAGDYN_TEMPL void MAGDYN_INST::SetPerformChecks(bool b)
{
	m_perform_checks = b;
}


MAGDYN_TEMPL void MAGDYN_INST::SetSilent(bool b)
{
	m_silent = b;
}


MAGDYN_TEMPL void MAGDYN_INST::SetPhaseSign(t_real sign)
{
	m_phase_sign = sign;
}


MAGDYN_TEMPL void MAGDYN_INST::SetCholeskyMaxTries(t_size max_tries)
{
	m_tries_chol = max_tries;
}


MAGDYN_TEMPL void MAGDYN_INST::SetCholeskyInc(t_real delta)
{
	m_delta_chol = delta;
}



MAGDYN_TEMPL void MAGDYN_INST::SetMagneticFormFactor(const std::string& ffact)
{
	m_magffact_formula = ffact;
	if(m_magffact_formula == "")
		return;

	// parse the given formula
	m_magffact = GetExprParser();
	m_magffact.SetInvalid0(false);
	m_magffact.register_var("Q", 0.);

	if(!m_magffact.parse_noexcept(ffact))
	{
		m_magffact_formula = "";

		CERR_OPT << "Magdyn error: Magnetic form facor formula: \""
			<< ffact << "\" could not be parsed."
			<< std::endl;
	}
}



MAGDYN_TEMPL void MAGDYN_INST::SetExternalField(const MAGDYN_TYPE::ExternalField& field)
{
	m_field = field;

	// normalise direction vector
	const t_real len = tl2::norm<t_vec_real>(m_field.dir);
	if(!tl2::equals_0<t_real>(len, m_eps))
		m_field.dir /= len;
}



MAGDYN_TEMPL void MAGDYN_INST::RotateExternalField(const t_vec_real& axis, t_real angle)
{
	const t_mat_real rot = tl2::rotation<t_mat_real, t_vec_real>(
		axis, angle, false);
	m_field.dir = rot * m_field.dir;
}



MAGDYN_TEMPL void MAGDYN_INST::RotateExternalField(t_real x, t_real y, t_real z, t_real angle)
{
	RotateExternalField(tl2::create<t_vec_real>({ x, y, z }), angle);
}



/**
 * set the ordering wave vector (e.g., the helix pitch) for incommensurate structures
 */
MAGDYN_TEMPL void MAGDYN_INST::SetOrderingWavevector(const t_vec_real& ordering)
{
	m_ordering = ordering;
	m_is_incommensurate = !tl2::equals_0<t_vec_real>(m_ordering, m_eps);
}



/**
 * set the rotation axis for the ordering wave vector
 */
MAGDYN_TEMPL void MAGDYN_INST::SetRotationAxis(const t_vec_real& axis)
{
	m_rotaxis = axis;

	// normalise
	const t_real len = tl2::norm<t_vec_real>(m_rotaxis);
	if(!tl2::equals_0<t_real>(len, m_eps))
		m_rotaxis /= len;
}



MAGDYN_TEMPL void MAGDYN_INST::SetCalcHamiltonian(bool H, bool Hp, bool Hm)
{
	m_calc_H  = H;
	m_calc_Hp = Hp;
	m_calc_Hm = Hm;
}



MAGDYN_TEMPL void MAGDYN_INST::AddVariable(MAGDYN_TYPE::Variable&& var)
{
	m_variables.emplace_back(std::forward<Variable&&>(var));
}



MAGDYN_TEMPL void MAGDYN_INST::SetVariable(MAGDYN_TYPE::Variable&& var)
{
	// is a variable with the same name already registered?
	auto iter = std::find_if(m_variables.begin(), m_variables.end(),
		[&var](const Variable& thevar)
	{
		return thevar.name == var.name;
	});

	if(iter == m_variables.end())
	{
		// add a new variable
		AddVariable(std::forward<Variable&&>(var));
	}
	else
	{
		// replace the value of an existing variable
		iter->value = var.value;
	}
}



MAGDYN_TEMPL void MAGDYN_INST::AddMagneticSite(MAGDYN_TYPE::MagneticSite&& site)
{
	m_sites.emplace_back(std::forward<MagneticSite&&>(site));
}



MAGDYN_TEMPL void MAGDYN_INST::AddExchangeTerm(MAGDYN_TYPE::ExchangeTerm&& term)
{
	m_exchange_terms.emplace_back(std::forward<ExchangeTerm&&>(term));
}



/**
 * calculate the B matrix from the crystal lattice
 */
MAGDYN_TEMPL
void MAGDYN_INST::SetCrystalLattice(t_real a, t_real b, t_real c,
	t_real alpha, t_real beta, t_real gamma)
{
	try
	{
		m_xtallattice[0] = a;
		m_xtallattice[1] = b;
		m_xtallattice[2] = c;

		m_xtalangles[0] = alpha;
		m_xtalangles[1] = beta;
		m_xtalangles[2] = gamma;

		// crystal fractional coordinate trafo matrices
		m_xtalA = tl2::A_matrix<t_mat_real>(a, b, c, alpha, beta, gamma);
		m_xtalB = tl2::B_matrix<t_mat_real>(a, b, c, alpha, beta, gamma);
	}
	catch(const std::exception& ex)
	{
		m_xtalA = m_xtalB = tl2::unit<t_mat_real>(3);

		CERR_OPT << "Magdyn error: Could not calculate crystal matrices."
			<< std::endl;
	}
}



/**
 * calculate the UB matrix from the scattering plane and the crystal lattice
 * note: SetCrystalLattice() has to be called before this function
 */
MAGDYN_TEMPL
void MAGDYN_INST::SetScatteringPlane(t_real ah, t_real ak, t_real al,
	t_real bh, t_real bk, t_real bl)
{
	try
	{
		m_scatteringplane[0] = tl2::create<t_vec_real>({ ah, ak, al });
		m_scatteringplane[1] = tl2::create<t_vec_real>({ bh, bk, bl });
		m_scatteringplane[2] = tl2::cross(m_xtalB, m_scatteringplane[0], m_scatteringplane[1]);

		for(std::uint8_t i = 0; i < 3; ++i)
			m_scatteringplane[i] /= tl2::norm<t_vec_real>(m_scatteringplane[i]);

		m_xtalUB = tl2::UB_matrix(m_xtalB,
			m_scatteringplane[0], m_scatteringplane[1], m_scatteringplane[2]);
		bool inv_ok = false;
		std::tie(m_xtalUBinv, inv_ok) = tl2::inv(m_xtalUB);

		if(!inv_ok)
		{
			CERR_OPT << "Magdyn error: UB matrix is not invertible."
				<< std::endl;
		}
	}
	catch(const std::exception& ex)
	{
		m_xtalUB = m_xtalUBinv = tl2::unit<t_mat_real>(3);

		CERR_OPT << "Magdyn error: Could not calculate scattering plane matrices."
			<< std::endl;
	}
}
// --------------------------------------------------------------------



// --------------------------------------------------------------------
/**
 * get an expression parser object with registered variables
 */
MAGDYN_TEMPL
tl2::ExprParser<t_cplx> MAGDYN_INST::GetExprParser() const
{
	tl2::ExprParser<t_cplx> parser;

	// register all variables
	parser.SetAutoregisterVariables(false);
	for(const Variable& var : m_variables)
		parser.register_var(var.name, var.value);

	return parser;
}
// --------------------------------------------------------------------



// --------------------------------------------------------------------
// sanity checks
// --------------------------------------------------------------------
/**
 * check if the site index is valid
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CheckMagneticSite(t_size idx, bool print_error) const
{
	if(!m_perform_checks)
		return true;

	if(idx >= m_sites.size())
	{
		if(print_error)
		{
			CERR_OPT << "Magdyn error: Site index " << idx
				<< " is out of bounds."
				<< std::endl;
		}

		return false;
	}

	return true;
}



/**
 * check if the term index is valid
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CheckExchangeTerm(t_size idx, bool print_error) const
{
	if(!m_perform_checks)
		return true;

	if(idx >= m_exchange_terms.size())
	{
		if(print_error)
		{
			CERR_OPT << "Magdyn error: Coupling index " << idx
				<< " is out of bounds."
				<< std::endl;
		}

		return false;
	}

	return true;
}



/**
 * check if imaginary weights remain
 */
MAGDYN_TEMPL
bool MAGDYN_INST::CheckImagWeights(const t_vec_real& Q_rlu,
	const MAGDYN_TYPE::EnergiesAndWeights& Es_and_S) const
{
	if(!m_perform_checks)
		return true;

	using namespace tl2_ops;
	bool ok = true;

	for(const EnergyAndWeight& EandS : Es_and_S)
	{
		// imaginary parts should be gone after UniteEnergies()
		if(std::abs(EandS.S_perp_sum.imag()) > m_eps ||
			std::abs(EandS.S_sum.imag()) > m_eps)
		{
			ok = false;

			CERR_OPT << "Magdyn warning: Remaining imaginary S(Q, E) component at Q = "
				<< Q_rlu << " and E = " << EandS.E
				<< ": imag(S) = " << EandS.S_sum.imag()
				<< ", imag(S_perp) = " << EandS.S_perp_sum.imag()
				<< "." << std::endl;
		}
	}

	return ok;
}
// --------------------------------------------------------------------

#endif
