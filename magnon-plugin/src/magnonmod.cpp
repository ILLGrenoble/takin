/**
 * S(Q, E) module for magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv2, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

//#define MAGNONMOD_USE_CPLX
//#define MAGNONMOD_ALLOW_QSIGNS

#include "magnonmod.h"

#include "core/libs/version.h"
#include "tlibs/string/string.h"
#include "tlibs/math/math.h"
#include "tlibs/phys/neutrons.h"

using t_real = MagnonMod::t_real;


// ----------------------------------------------------------------------------
// constructors

MagnonMod::MagnonMod()
{
	SqwBase::m_bOk = false;
}


MagnonMod::MagnonMod(const std::string& cfg_file) : MagnonMod()
{
	if(cfg_file == "")
	{
		tl::log_info("No config file given for magnon module.");
		SqwBase::m_bOk = false;
		return;
	}

	tl::log_info("Magnon module config file: \"", cfg_file, "\".");

	// load config file
	SqwBase::m_bOk = m_dyn.Load(cfg_file);
}


MagnonMod::~MagnonMod()
{
	m_dyn.Clear();
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// dispersion, spectral weight and structure factor

std::tuple<std::vector<t_real>, std::vector<t_real>>
	MagnonMod::disp(t_real h, t_real k, t_real l) const
{
#ifdef MAGNONMOD_ALLOW_QSIGNS
	if(m_Qsigns.size() == 3)
	{
		h *= m_Qsigns[0];
		k *= m_Qsigns[1];
		l *= m_Qsigns[2];
	}
#endif

	// calculate dispersion relation
	auto modes = m_dyn.CalcEnergies(h, k, l, false);

	std::vector<t_real> energies;
	std::vector<t_real> weights;
	energies.reserve(modes.size());
	weights.reserve(modes.size());

	for(const auto& mode : modes)
	{
		energies.push_back(mode.E);

		if(m_channel >= 0 && m_channel < 3)
			weights.push_back(std::abs(mode.S_perp(m_channel, m_channel).real()));
		else
			weights.push_back(mode.weight);
	}

	return std::make_tuple(energies, weights);
}


t_real MagnonMod::operator()(t_real h, t_real k, t_real l, t_real E) const
{
	// bose factor
	t_real bose = 1.;
	if(!m_use_model_bose)
	{
		// calculate bose factor here (not in model)
		bose = tl::bose_cutoff(E, m_T, m_dyn.GetBoseCutoffEnergy());
	}

	std::vector<t_real> Es, Ws;
	std::tie(Es, Ws) = disp(h, k, l);

	// incoherent peak
	t_real incoh = 0.;
	if(!tl::float_equal(m_incoh_amp, t_real(0)))
		incoh = tl::gauss_model(E, t_real(0),
			m_incoh_sigma, m_incoh_amp, t_real(0));

	// magnon peaks
	t_real S = 0.;
	for(std::size_t iE = 0; iE < Es.size(); ++iE)
	{
		if(!tl::float_equal(Ws[iE], t_real(0)))
			S += tl::gauss_model(E, Es[iE], m_sigma, Ws[iE], t_real(0));
	}

	return m_S0*S*bose + incoh;
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// get & set variables

std::vector<MagnonMod::t_var> MagnonMod::GetVars() const
{
	std::vector<t_var> vars;

	// get external field
	const t_magdyn::ExternalField& field = m_dyn.GetExternalField();
	std::vector<t_real> B;
	if(field.dir.size() == 3)
	{
		B = std::vector<t_real>{{
			field.dir[0], field.dir[1], field.dir[2] }};
	}
	else
	{
		B = std::vector<t_real>{{ 0., 0., 1. }};
	}

	vars.push_back(SqwBase::t_var{
		"sigma", "real", tl::var_to_str(m_sigma)});
	vars.push_back(SqwBase::t_var{
		"inc_amp", "real", tl::var_to_str(m_incoh_amp)});
	vars.push_back(SqwBase::t_var{
		"inc_sigma", "real", tl::var_to_str(m_incoh_sigma)});
	vars.push_back(SqwBase::t_var{
		"S0", "real", tl::var_to_str(m_S0)});
	vars.push_back(SqwBase::t_var{
		"T", "real", tl::var_to_str(m_T /*m_dyn.GetTemperature()*/)});
	vars.push_back(SqwBase::t_var{
		"cutoff", "real", tl::var_to_str(m_dyn.GetBoseCutoffEnergy())});
	vars.push_back(SqwBase::t_var{
		"use_model_bose", "int", tl::var_to_str((int)m_use_model_bose)});
	vars.push_back(SqwBase::t_var{
		"channel", "int", tl::var_to_str(m_channel)});
	vars.push_back(SqwBase::t_var{
		"B_dir", "vector", vec_to_str(B)});
	vars.push_back(SqwBase::t_var{
		"B_mag", "real", tl::var_to_str(field.mag)});
	vars.push_back(SqwBase::t_var{
		"B_align_spins", "int", tl::var_to_str((int)field.align_spins)});
	vars.push_back(SqwBase::t_var{
		"silent", "real", tl::var_to_str((int)m_dyn.GetSilent())});
#ifdef MAGNONMOD_ALLOW_QSIGNS
	vars.push_back(SqwBase::t_var{
		"Q_signs", "vector", vec_to_str(m_Qsigns)});
#endif

	// get variables from the magdyn model
	for(const auto& modelvar : m_dyn.GetVariables())
	{
#ifdef MAGNONMOD_USE_CPLX
		vars.push_back(SqwBase::t_var{
			modelvar.name, "complex",
			tl::var_to_str(modelvar.value)});
#else
		vars.push_back(SqwBase::t_var{
			modelvar.name, "real",
			tl::var_to_str(modelvar.value.real())});
#endif
	}

	return vars;
}


void MagnonMod::SetVars(const std::vector<MagnonMod::t_var>& vars)
{
	if(!vars.size())
		return;

	bool set_model_temp = false;
	bool calc_sites = false;
	bool calc_terms = false;

	for(const SqwBase::t_var& var : vars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "sigma")
			m_sigma = tl::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp")
			m_incoh_amp = tl::str_to_var<decltype(m_incoh_amp)>(strVal);
		else if(strVar == "inc_sigma")
			m_incoh_sigma = tl::str_to_var<decltype(m_incoh_sigma)>(strVal);
		else if(strVar == "S0")
			m_S0 = tl::str_to_var<decltype(m_S0)>(strVal);
		else if(strVar == "T")
		{
			m_T = tl::str_to_var<t_real>(strVal);
			set_model_temp = true;
		}
		else if(strVar == "cutoff")
			m_dyn.SetBoseCutoffEnergy(tl::str_to_var<t_real>(strVal));
		else if(strVar == "use_model_bose")
		{
			m_use_model_bose = (tl::str_to_var<int>(strVal) != 0);
			set_model_temp = true;
		}
		else if(strVar == "channel")
			m_channel = tl::str_to_var<int>(strVal);
		else if(strVar == "B_dir")
		{
			std::vector<t_real> dir = str_to_vec<std::vector<t_real>>(strVal);
			if(dir.size() == 3)
			{
				t_magdyn::ExternalField field = m_dyn.GetExternalField();
				field.dir = tl2::create<t_vec_real>({dir[0], dir[1], dir[2]});

				m_dyn.SetExternalField(field);
				calc_sites = true;
			}
			else
			{
				tl::log_err("Invalid field direction.");
			}
		}
		else if(strVar == "B_mag")
		{
			t_magdyn::ExternalField field = m_dyn.GetExternalField();
			field.mag = tl::str_to_var<decltype(m_S0)>(strVal);

			m_dyn.SetExternalField(field);
			calc_sites = true;
		}
		else if(strVar == "B_align_spins")
		{
			t_magdyn::ExternalField field = m_dyn.GetExternalField();
			field.align_spins = (tl::str_to_var<int>(strVal) != 0);

			m_dyn.SetExternalField(field);
			calc_sites = true;
		}
		else if(strVar == "silent")
			m_dyn.SetSilent(tl::str_to_var<int>(strVal) != 0);
#ifdef MAGNONMOD_ALLOW_QSIGNS
		else if(strVar == "Q_signs")
		{
			std::vector<t_real> signs = str_to_vec<std::vector<t_real>>(strVal);
			if(signs.size() == 3)
			{
				m_Qsigns[0] = signs[0];
				m_Qsigns[1] = signs[1];
				m_Qsigns[2] = signs[2];
			}
		}
#endif
		else
		{
			// set model variables
			tl::log_info("Model variable: ", strVar, " = ",  strVal, ".");

			t_magdyn::Variable modelvar;
			modelvar.name = strVar;
#ifdef MAGNONMOD_USE_CPLX
			modelvar.value = tl::str_to_var<t_cplx>(strVal);
#else
			modelvar.value = tl::str_to_var<t_real>(strVal);
#endif
			m_dyn.SetVariable(std::move(modelvar));
			calc_terms = true;
		}
	}

	if(set_model_temp)
	{
		if(!m_use_model_bose)
			m_dyn.SetTemperature(-1.);
		else
			m_dyn.SetTemperature(m_T);
	}

	if(calc_sites)
	{
		m_dyn.CalcExternalField();
		m_dyn.CalcMagneticSites();
	}
	if(calc_terms)
		m_dyn.CalcExchangeTerms();
}


bool MagnonMod::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	return SqwBase::SetVarIfAvail(strKey, strNewVal);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// copy

SqwBase* MagnonMod::shallow_copy() const
{
	MagnonMod *mod = new MagnonMod();

	mod->m_sigma = this->m_sigma;
	mod->m_incoh_amp = this->m_incoh_amp;
	mod->m_incoh_sigma = this->m_incoh_sigma;
	mod->m_S0 = this->m_S0;
	mod->m_dyn = this->m_dyn;
	mod->m_channel = this->m_channel;
	mod->m_use_model_bose = this->m_use_model_bose;
	mod->m_T = this->m_T;
#ifdef MAGNONMOD_ALLOW_QSIGNS
	mod->m_Qsigns = this->m_Qsigns;
#endif

	return mod;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// SO interface

static const char* g_help = R"RAWSTR(Magnetic Dynamics Module.

This module serves as an interface between the magnon calculator ("Tools" -> "Magnetic Dynamics...") and the resolution-convolution simulator and fitter.

Please refer to the Takin help for more information and tutorials.)RAWSTR";

#ifndef __MINGW32__

#include <boost/dll/alias.hpp>

std::tuple<std::string, std::string, std::string, std::string> sqw_info()
{
	return std::make_tuple(TAKIN_VER, "magnonmod", "Magnetic Dynamics [Magpie]", g_help);
}

std::shared_ptr<SqwBase> sqw_construct(const std::string& cfg_file)
{
	return std::make_shared<MagnonMod>(cfg_file);
}

// exports from so file
BOOST_DLL_ALIAS(sqw_info, takin_sqw_info);
BOOST_DLL_ALIAS(sqw_construct, takin_sqw);

#else  // mingw exports

extern "C" std::tuple<std::string, std::string, std::string, std::string> takin_sqw_info()
{
	return std::make_tuple(TAKIN_VER, "magnonmod", "Magnetic Dynamics [Magpie]", g_help);
}

extern "C" std::shared_ptr<SqwBase> takin_sqw(const std::string& cfg_file)
{
	return std::make_shared<MagnonMod>(cfg_file);
}

#endif
// ----------------------------------------------------------------------------
