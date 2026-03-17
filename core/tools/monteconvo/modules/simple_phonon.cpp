/**
 * monte carlo convolution tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015 -- 2018
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "simple_phonon.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/phys/neutrons.h"
#include <fstream>
#include <list>

using t_real = t_real_reso;



t_real SqwPhonon::phonon_disp(t_real dq, t_real da, t_real df)
{
	return std::abs(da*std::sin(dq*df));
}


void SqwPhonon::create()
{
#ifdef USE_RTREE
	m_rt = std::make_shared<tl::Rt<t_real,3,RT_ELEMS>>();
#else
	m_kd = std::make_shared<tl::Kd<t_real>>();
#endif

	const bool bSaveOnlyIndices = 1;
	destroy();

	if(m_vecBragg.size()==0 || m_vecLA.size()==0 || m_vecTA1.size()==0 || m_vecTA2.size()==0)
	{
		// examples
		m_vecBragg = tl::make_vec({1., 0., 0.});
		m_vecLA = tl::make_vec({1., 0., 0.});
		m_vecTA1 = tl::make_vec({0., 1., 0.});
		m_vecTA2 = tl::make_vec({0., 0., 1.});
	}

	m_vecLA /= ublas::norm_2(m_vecLA);
	m_vecTA1 /= ublas::norm_2(m_vecTA1);
	m_vecTA2 /= ublas::norm_2(m_vecTA2);

	tl::log_info("LA: ", m_vecLA);
	tl::log_info("TA1: ", m_vecTA1);
	tl::log_info("TA2: ", m_vecTA2);

	std::list<std::vector<t_real>> lst;
	// generate qs for one brillouin zone
	for(t_real dq =- 1.; dq < 1.; dq += 1./t_real(m_iNumqs))
	{
		ublas::vector<t_real> vecQLA = dq*m_vecLA;
		ublas::vector<t_real> vecQTA1 = dq*m_vecTA1;
		ublas::vector<t_real> vecQTA2 = dq*m_vecTA2;

		t_real dELA = phonon_disp(dq, m_dLA_amp, m_dLA_freq);
		t_real dETA1 = phonon_disp(dq, m_dTA1_amp, m_dTA1_freq);
		t_real dETA2 = phonon_disp(dq, m_dTA2_amp, m_dTA2_freq);

		t_real dLA_E_HWHM = m_dLA_E_HWHM;
		t_real dLA_q_HWHM = m_dLA_q_HWHM;
		t_real dTA1_E_HWHM = m_dTA1_E_HWHM;
		t_real dTA1_q_HWHM = m_dTA1_q_HWHM;
		t_real dTA2_E_HWHM = m_dTA2_E_HWHM;
		t_real dTA2_q_HWHM = m_dTA2_q_HWHM;

		t_real dLA_S0 = m_dLA_S0;
		t_real dTA1_S0 = m_dTA1_S0;
		t_real dTA2_S0 = m_dTA2_S0;

		if(bSaveOnlyIndices)
		{
			dTA1_E_HWHM = dTA1_q_HWHM = dTA1_S0 = -1.;
			dTA2_E_HWHM = dTA2_q_HWHM = dTA2_S0 = -2.;
			dLA_E_HWHM = dLA_q_HWHM = dLA_S0 = -3.;
		}

		// only generate exact phonon branches, no arcs
		if(m_iNumArc == 0 || m_iNumArc == 1)
		{
			lst.push_back(std::vector<t_real>({
				vecQLA[0] + m_vecBragg[0],
				vecQLA[1] + m_vecBragg[1],
				vecQLA[2] + m_vecBragg[2],
				dELA, dLA_S0, dLA_E_HWHM, dLA_q_HWHM}));
			lst.push_back(std::vector<t_real>({
				vecQTA1[0] + m_vecBragg[0],
				vecQTA1[1] + m_vecBragg[1],
				vecQTA1[2] + m_vecBragg[2],
				dETA1, dTA1_S0, dTA1_E_HWHM, dTA1_q_HWHM}));
			lst.push_back(std::vector<t_real>({
				vecQTA2[0] + m_vecBragg[0],
				vecQTA2[1] + m_vecBragg[1],
				vecQTA2[2] + m_vecBragg[2],
				dETA2, dTA2_S0, dTA2_E_HWHM, dTA2_q_HWHM}));
		}
		else
		{
			const t_real dArcMax = std::abs(tl::d2r(m_dArcMax));
			for(t_real dph = -dArcMax; dph <= dArcMax; dph += 1./t_real(m_iNumArc))
			for(t_real dth = -dArcMax; dth <= dArcMax; dth += 1./t_real(m_iNumArc))
			{
				// ta2
				ublas::vector<t_real> vecArcTA2 = tl::sph_shell(vecQTA2, dph, dth) + m_vecBragg;;
				lst.push_back(std::vector<t_real>({
					vecArcTA2[0], vecArcTA2[1], vecArcTA2[2],
					dETA2, dTA2_S0, dTA2_E_HWHM, dTA2_q_HWHM}));

				// ta1
				ublas::vector<t_real> vecArcTA1 = tl::sph_shell(vecQTA1, dph, dth) + m_vecBragg;;
				lst.push_back(std::vector<t_real>({
					vecArcTA1[0], vecArcTA1[1], vecArcTA1[2],
					dETA1, dTA1_S0, dTA1_E_HWHM, dTA1_q_HWHM}));

				// la
				ublas::vector<t_real> vecArcLA = tl::sph_shell(vecQLA, dph, dth) + m_vecBragg;;
				lst.push_back(std::vector<t_real>({
					vecArcLA[0], vecArcLA[1], vecArcLA[2],
					dELA, dLA_S0, dLA_E_HWHM, dLA_q_HWHM}));
			}
		}
	}

	tl::log_info("Generated ", lst.size(), " S(Q, E) points.");
#ifdef USE_RTREE
	m_rt->Load(lst);
	tl::log_info("Generated R* tree.");
#else
	m_kd->Load(lst, 3);
	tl::log_info("Generated k-d tree.");
#endif

	m_bOk = true;
}


void SqwPhonon::destroy()
{
#ifdef USE_RTREE
	m_rt->Unload();
#else
	m_kd->Unload();
#endif
}


SqwPhonon::SqwPhonon(const ublas::vector<t_real>& vecBragg,
	const ublas::vector<t_real>& vecTA1,
	const ublas::vector<t_real>& vecTA2,
	t_real dLA_amp, t_real dLA_freq, t_real dLA_E_HWHM, t_real dLA_q_HWHM, t_real dLA_S0,
	t_real dTA1_amp, t_real dTA1_freq, t_real dTA1_E_HWHM, t_real dTA1_q_HWHM, t_real dTA1_S0,
	t_real dTA2_amp, t_real dTA2_freq, t_real dTA2_E_HWHM, t_real dTA2_q_HWHM, t_real dTA2_S0,
	t_real dT)
		: m_vecBragg(vecBragg), m_vecLA(vecBragg),
		m_vecTA1(vecTA1), m_vecTA2(vecTA2),
		m_dLA_amp(dLA_amp), m_dLA_freq(dLA_freq),
		m_dLA_E_HWHM(dLA_E_HWHM), m_dLA_q_HWHM(dLA_q_HWHM),
		m_dLA_S0(dLA_S0),

		m_dTA1_amp(dTA1_amp), m_dTA1_freq(dTA1_freq),
		m_dTA1_E_HWHM(dTA1_E_HWHM), m_dTA1_q_HWHM(dTA1_q_HWHM),
		m_dTA1_S0(dTA1_S0),

		m_dTA2_amp(dTA2_amp), m_dTA2_freq(dTA2_freq),
		m_dTA2_E_HWHM(dTA2_E_HWHM), m_dTA2_q_HWHM(dTA2_q_HWHM),
		m_dTA2_S0(dTA2_S0),

		m_dT(dT)
{
	create();
}


SqwPhonon::SqwPhonon(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr)
		tl::log_warn("Cannot open phonon config file \"", pcFile, "\".");

	// if a file is given, load the parameters
	if(!!ifstr)
	{
		std::string strLine;
		while(std::getline(ifstr, strLine))
		{
			std::vector<std::string> vecToks;
			tl::get_tokens<std::string>(strLine, std::string("=,"), vecToks);
			std::for_each(vecToks.begin(), vecToks.end(), [](std::string& str) {tl::trim(str); });

			if(vecToks.size() == 0) continue;

			if(vecToks[0] == "num_qs") m_iNumqs = tl::str_to_var<unsigned int>(vecToks[1]);
			if(vecToks[0] == "num_arc") m_iNumArc = tl::str_to_var<unsigned int>(vecToks[1]);
			if(vecToks[0] == "arc_max") m_dArcMax = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "G") m_vecLA = m_vecBragg = tl::make_vec({tl::str_to_var_parse<t_real>(vecToks[1]), tl::str_to_var_parse<t_real>(vecToks[2]), tl::str_to_var_parse<t_real>(vecToks[3])});
			else if(vecToks[0] == "TA1") m_vecTA1 = tl::make_vec({tl::str_to_var_parse<t_real>(vecToks[1]), tl::str_to_var_parse<t_real>(vecToks[2]), tl::str_to_var_parse<t_real>(vecToks[3])});
			else if(vecToks[0] == "TA2") m_vecTA2 = tl::make_vec({tl::str_to_var_parse<t_real>(vecToks[1]), tl::str_to_var_parse<t_real>(vecToks[2]), tl::str_to_var_parse<t_real>(vecToks[3])});

			else if(vecToks[0] == "LA_amp") m_dLA_amp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "LA_freq") m_dLA_freq = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "LA_E_HWHM") m_dLA_E_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "LA_q_HWHM") m_dLA_q_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "LA_S0") m_dLA_S0 = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "TA1_amp") m_dTA1_amp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA1_freq") m_dTA1_freq = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA1_E_HWHM") m_dTA1_E_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA1_q_HWHM") m_dTA1_q_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA1_S0") m_dTA1_S0 = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "TA2_amp") m_dTA2_amp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA2_freq") m_dTA2_freq = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA2_E_HWHM") m_dTA2_E_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA2_q_HWHM") m_dTA2_q_HWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "TA2_S0") m_dTA2_S0 = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "inc_amp") m_dIncAmp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "inc_sig") m_dIncSig = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "use_q_distr") m_bUseQDistr = (tl::str_to_var<unsigned int>(vecToks[1]) != 0);

			else if(vecToks[0] == "T") m_dT = tl::str_to_var_parse<t_real>(vecToks[1]);
		}
	}

	create();
}


/**
 * TODO: dispersion E(Q)
 */
/*std::tuple<std::vector<t_real>, std::vector<t_real>>
SqwPhonon::disp(t_real dh, t_real dk, t_real dl) const
{
}*/


t_real SqwPhonon::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	std::vector<t_real> vechklE = { dh, dk, dl, dE };
#ifdef USE_RTREE
	if(!m_rt->IsPointInGrid(vechklE))
		return 0.;
	std::vector<t_real> vec = m_rt->GetNearestNode(vechklE);
#else
	if(!m_kd->IsPointInGrid(vechklE))
		return 0.;
	std::vector<t_real> vec = m_kd->GetNearestNode(vechklE);
#endif

	t_real dE0 = vec[3];
	t_real dS = vec[4];
	t_real dT = m_dT;
	t_real dE_HWHM = vec[5];
	t_real dQ_HWHM = vec[6];

	// index, not value
	if(dE_HWHM < 0.)
	{
		if(tl::float_equal<t_real>(dE_HWHM, -1., 0.1))      // TA1
			dE_HWHM = m_dTA1_E_HWHM;
		else if(tl::float_equal<t_real>(dE_HWHM, -2., 0.1)) // TA2
			dE_HWHM = m_dTA2_E_HWHM;
		else if(tl::float_equal<t_real>(dE_HWHM, -3., 0.1)) // LA
			dE_HWHM = m_dLA_E_HWHM;
	}
	if(dQ_HWHM < 0.)
	{
		if(tl::float_equal<t_real>(dQ_HWHM, -1., 0.1))      // TA1
			dQ_HWHM = m_dTA1_q_HWHM;
		else if(tl::float_equal<t_real>(dQ_HWHM, -2., 0.1)) // TA2
			dQ_HWHM = m_dTA2_q_HWHM;
		else if(tl::float_equal<t_real>(dQ_HWHM, -3., 0.1)) // LA
			dQ_HWHM = m_dLA_q_HWHM;
	}
	if(dS < 0.)
	{
		if(tl::float_equal<t_real>(dS, -1., 0.1))           // TA1
			dS = m_dTA1_S0;
		else if(tl::float_equal<t_real>(dS, -2., 0.1))      // TA2
			dS = m_dTA2_S0;
		else if(tl::float_equal<t_real>(dS, -3., 0.1))      // LA
			dS = m_dLA_S0;
	}

	// distribution for the distances between the wanted and the available qs
	t_real dq_distr = 1.;
	if(m_bUseQDistr)
	{
		t_real dqDist = std::sqrt(
			std::pow(vec[0] - vechklE[0], 2.) +
			std::pow(vec[1] - vechklE[1], 2.) +
			std::pow(vec[2] - vechklE[2], 2.));
		dq_distr = tl::gauss_model/*_amp*/<t_real>(dqDist, 0., dQ_HWHM*tl::get_HWHM2SIGMA<t_real>(), 1., 0.);
	}

	t_real dInc = 0.;
	if(!tl::float_equal<t_real>(m_dIncAmp, 0.))
		dInc = tl::gauss_model<t_real>(dE, 0., m_dIncSig, m_dIncAmp, 0.);

	return dInc + dS * std::abs(tl::DHO_model<t_real>(dE, dT, dE0, dE_HWHM, 1., 0.)) * dq_distr;
}


std::vector<SqwBase::t_var> SqwPhonon::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"num_qs", "uint", tl::var_to_str(m_iNumqs)});
	vecVars.push_back(SqwBase::t_var{"num_arc", "uint", tl::var_to_str(m_iNumArc)});
	vecVars.push_back(SqwBase::t_var{"arc_max", "real", tl::var_to_str(m_dArcMax)});

	vecVars.push_back(SqwBase::t_var{"G", "vector", vec_to_str(m_vecBragg)});
	vecVars.push_back(SqwBase::t_var{"TA1", "vector", vec_to_str(m_vecTA1)});
	vecVars.push_back(SqwBase::t_var{"TA2", "vector", vec_to_str(m_vecTA2)});

	vecVars.push_back(SqwBase::t_var{"LA_amp", "real", tl::var_to_str(m_dLA_amp)});
	vecVars.push_back(SqwBase::t_var{"LA_freq", "real", tl::var_to_str(m_dLA_freq)});
	vecVars.push_back(SqwBase::t_var{"LA_E_HWHM", "real", tl::var_to_str(m_dLA_E_HWHM)});
	vecVars.push_back(SqwBase::t_var{"LA_q_HWHM", "real", tl::var_to_str(m_dLA_q_HWHM)});
	vecVars.push_back(SqwBase::t_var{"LA_S0", "real", tl::var_to_str(m_dLA_S0)});

	vecVars.push_back(SqwBase::t_var{"TA1_amp", "real", tl::var_to_str(m_dTA1_amp)});
	vecVars.push_back(SqwBase::t_var{"TA1_freq", "real", tl::var_to_str(m_dTA1_freq)});
	vecVars.push_back(SqwBase::t_var{"TA1_E_HWHM", "real", tl::var_to_str(m_dTA1_E_HWHM)});
	vecVars.push_back(SqwBase::t_var{"TA1_q_HWHM", "real", tl::var_to_str(m_dTA1_q_HWHM)});
	vecVars.push_back(SqwBase::t_var{"TA1_S0", "real", tl::var_to_str(m_dTA1_S0)});

	vecVars.push_back(SqwBase::t_var{"TA2_amp", "real", tl::var_to_str(m_dTA2_amp)});
	vecVars.push_back(SqwBase::t_var{"TA2_freq", "real", tl::var_to_str(m_dTA2_freq)});
	vecVars.push_back(SqwBase::t_var{"TA2_E_HWHM", "real", tl::var_to_str(m_dTA2_E_HWHM)});
	vecVars.push_back(SqwBase::t_var{"TA2_q_HWHM", "real", tl::var_to_str(m_dTA2_q_HWHM)});
	vecVars.push_back(SqwBase::t_var{"TA2_S0", "real", tl::var_to_str(m_dTA2_S0)});

	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sig", "real", tl::var_to_str(m_dIncSig)});
	vecVars.push_back(SqwBase::t_var{"use_q_distr", "uint", tl::var_to_str(m_bUseQDistr)});

	vecVars.push_back(SqwBase::t_var{"T", "real", tl::var_to_str(m_dT)});

	return vecVars;
}


void SqwPhonon::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(vecVars.size() == 0)
		return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "num_qs") m_iNumqs = tl::str_to_var<decltype(m_iNumqs)>(strVal);
		if(strVar == "num_arc") m_iNumArc = tl::str_to_var<decltype(m_iNumArc)>(strVal);
		if(strVar == "arc_max") m_dArcMax = tl::str_to_var<decltype(m_dArcMax)>(strVal);

		else if(strVar == "G") m_vecLA = m_vecBragg = str_to_vec<decltype(m_vecBragg)>(strVal);
		else if(strVar == "TA1") m_vecTA1 = str_to_vec<decltype(m_vecTA1)>(strVal);
		else if(strVar == "TA2") m_vecTA2 = str_to_vec<decltype(m_vecTA2)>(strVal);

		else if(strVar == "LA_amp") m_dLA_amp = tl::str_to_var<decltype(m_dLA_amp)>(strVal);
		else if(strVar == "LA_freq") m_dLA_freq = tl::str_to_var<decltype(m_dLA_freq)>(strVal);
		else if(strVar == "LA_E_HWHM") m_dLA_E_HWHM = tl::str_to_var<decltype(m_dLA_E_HWHM)>(strVal);
		else if(strVar == "LA_q_HWHM") m_dLA_q_HWHM = tl::str_to_var<decltype(m_dLA_q_HWHM)>(strVal);
		else if(strVar == "LA_S0") m_dLA_S0 = tl::str_to_var<decltype(m_dLA_S0)>(strVal);

		else if(strVar == "TA1_amp") m_dTA1_amp = tl::str_to_var<decltype(m_dTA1_amp)>(strVal);
		else if(strVar == "TA1_freq") m_dTA1_freq = tl::str_to_var<decltype(m_dTA1_freq)>(strVal);
		else if(strVar == "TA1_E_HWHM") m_dTA1_E_HWHM = tl::str_to_var<decltype(m_dTA1_E_HWHM)>(strVal);
		else if(strVar == "TA1_q_HWHM") m_dTA1_q_HWHM = tl::str_to_var<decltype(m_dTA1_q_HWHM)>(strVal);
		else if(strVar == "TA1_S0") m_dTA1_S0 = tl::str_to_var<decltype(m_dTA1_S0)>(strVal);

		else if(strVar == "TA2_amp") m_dTA2_amp = tl::str_to_var<decltype(m_dTA2_amp)>(strVal);
		else if(strVar == "TA2_freq") m_dTA2_freq = tl::str_to_var<decltype(m_dTA2_freq)>(strVal);
		else if(strVar == "TA2_E_HWHM") m_dTA2_E_HWHM = tl::str_to_var<decltype(m_dTA2_E_HWHM)>(strVal);
		else if(strVar == "TA2_q_HWHM") m_dTA2_q_HWHM = tl::str_to_var<decltype(m_dTA2_q_HWHM)>(strVal);
		else if(strVar == "TA2_S0") m_dTA2_S0 = tl::str_to_var<decltype(m_dTA2_S0)>(strVal);

		else if(strVar == "inc_amp") m_dIncAmp = tl::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sig") m_dIncSig = tl::str_to_var<decltype(m_dIncSig)>(strVal);
		else if(strVar == "use_q_distr") m_bUseQDistr = (tl::str_to_var<unsigned int>(strVal) != 0);

		else if(strVar == "T") m_dT = tl::str_to_var<decltype(m_dT)>(strVal);
	}

	bool bRecreateTree = false;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		if(strVar != "T" && strVar.find("HWHM") == std::string::npos &&
			strVar.find("inc") == std::string::npos &&
			strVar.find("S0") == std::string::npos)
			bRecreateTree = true;
	}

	if(bRecreateTree)
		create();
}


SqwBase* SqwPhonon::shallow_copy() const
{
	SqwPhonon *pCpy = new SqwPhonon();
	*static_cast<SqwBase*>(pCpy) = *static_cast<const SqwBase*>(this);

#ifdef USE_RTREE
	pCpy->m_rt = m_rt;
#else
	pCpy->m_kd = m_kd;
#endif
	pCpy->m_iNumqs = m_iNumqs;
	pCpy->m_iNumArc = m_iNumArc;
	pCpy->m_dArcMax = m_dArcMax;

	pCpy->m_vecBragg = m_vecBragg;
	pCpy->m_vecLA = m_vecLA;
	pCpy->m_vecTA1 = m_vecTA1;
	pCpy->m_vecTA2 = m_vecTA2;

	pCpy->m_dLA_amp = m_dLA_amp;
	pCpy->m_dLA_freq = m_dLA_freq;
	pCpy->m_dLA_E_HWHM = m_dLA_E_HWHM;
	pCpy->m_dLA_q_HWHM = m_dLA_q_HWHM;
	pCpy->m_dLA_S0 = m_dLA_S0;

	pCpy->m_dTA1_amp = m_dTA1_amp;
	pCpy->m_dTA1_freq = m_dTA1_freq;
	pCpy->m_dTA1_E_HWHM = m_dTA1_E_HWHM;
	pCpy->m_dTA1_q_HWHM = m_dTA1_q_HWHM;
	pCpy->m_dTA1_S0 = m_dTA1_S0;

	pCpy->m_dTA2_amp = m_dTA2_amp;
	pCpy->m_dTA2_freq = m_dTA2_freq;
	pCpy->m_dTA2_E_HWHM = m_dTA2_E_HWHM;
	pCpy->m_dTA2_q_HWHM = m_dTA2_q_HWHM;
	pCpy->m_dTA2_S0 = m_dTA2_S0;

	pCpy->m_dIncAmp = m_dIncAmp;
	pCpy->m_dIncSig = m_dIncSig;
	pCpy->m_bUseQDistr = m_bUseQDistr;

	pCpy->m_dT = m_dT;
	return pCpy;
}



//------------------------------------------------------------------------------



// for model compatibility testing
#define BRANCH_PREFIX ""
//#define BRANCH_PREFIX "TA2_"

t_real SqwPhononSingleBranch::phonon_disp(t_real dq, t_real da, t_real df)
{
	return std::abs(da*std::sin(dq*df));
}


SqwPhononSingleBranch::SqwPhononSingleBranch(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr)
		tl::log_warn("Cannot open phonon config file \"", pcFile, "\".");

	// if a file is given, load the parameters
	if(!!ifstr)
	{
		std::string strLine;
		while(std::getline(ifstr, strLine))
		{
			std::vector<std::string> vecToks;
			tl::get_tokens<std::string>(strLine, std::string("=,"), vecToks);
			std::for_each(vecToks.begin(), vecToks.end(), [](std::string& str) {tl::trim(str); });

			if(vecToks.size() == 0)
				continue;

			if(vecToks[0] == "G")
			{
				m_vecBragg = tl::make_vec({
					tl::str_to_var_parse<t_real>(vecToks[1]),
					tl::str_to_var_parse<t_real>(vecToks[2]),
					tl::str_to_var_parse<t_real>(vecToks[3])});
			}

			else if(vecToks[0] == BRANCH_PREFIX"amp") m_damp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == BRANCH_PREFIX"freq") m_dfreq = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == BRANCH_PREFIX"E_HWHM") m_dHWHM = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == BRANCH_PREFIX"S0") m_dS0 = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "inc_amp") m_dIncAmp = tl::str_to_var_parse<t_real>(vecToks[1]);
			else if(vecToks[0] == "inc_sig") m_dIncSig = tl::str_to_var_parse<t_real>(vecToks[1]);

			else if(vecToks[0] == "T") m_dT = tl::str_to_var_parse<t_real>(vecToks[1]);
		}
	}

	// just use (100) if nothing else is defined
	if(m_vecBragg.size() < 3)
		m_vecBragg = tl::make_vec({1., 0., 0.});

	m_bOk = true;
}


/**
 * dispersion E(Q)
 */
std::tuple<std::vector<t_real>, std::vector<t_real>>
SqwPhononSingleBranch::disp(t_real dh, t_real dk, t_real dl) const
{
	dh -= m_vecBragg[0];
	dk -= m_vecBragg[1];
	dl -= m_vecBragg[2];

	t_real dq = std::sqrt(dh*dh + dk*dk + dl*dl);
	t_real dE0 = phonon_disp(dq, m_damp, m_dfreq);
	t_real dWeight = t_real(1);

	return std::make_tuple(std::vector<t_real>({dE0, -dE0}),
		std::vector<t_real>({dWeight, dWeight}));
}


/**
 * dynamical structure factor S(Q,E)
 */
t_real SqwPhononSingleBranch::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	std::vector<t_real> vecE0, vecW;
	std::tie(vecE0, vecW) = disp(dh, dk, dl);

	t_real dInc = 0.;
	if(!tl::float_equal<t_real>(m_dIncAmp, 0.))
		dInc = tl::gauss_model<t_real>(dE, 0., m_dIncSig, m_dIncAmp, 0.);

	return std::abs(tl::DHO_model<t_real>(dE, m_dT, vecE0[0], m_dHWHM, m_dS0*vecW[0], 0.)) + dInc;
}


std::vector<SqwBase::t_var> SqwPhononSingleBranch::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"G", "vector", vec_to_str(m_vecBragg)});

	vecVars.push_back(SqwBase::t_var{BRANCH_PREFIX"amp", "real", tl::var_to_str(m_damp)});
	vecVars.push_back(SqwBase::t_var{BRANCH_PREFIX"freq", "real", tl::var_to_str(m_dfreq)});
	vecVars.push_back(SqwBase::t_var{BRANCH_PREFIX"E_HWHM", "real", tl::var_to_str(m_dHWHM)});
	vecVars.push_back(SqwBase::t_var{BRANCH_PREFIX"S0", "real", tl::var_to_str(m_dS0)});

	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sig", "real", tl::var_to_str(m_dIncSig)});

	vecVars.push_back(SqwBase::t_var{"T", "real", tl::var_to_str(m_dT)});

	return vecVars;
}


void SqwPhononSingleBranch::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(vecVars.size() == 0)
		return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "G") m_vecBragg = str_to_vec<decltype(m_vecBragg)>(strVal);

		else if(strVar == BRANCH_PREFIX"amp") m_damp = tl::str_to_var<decltype(m_damp)>(strVal);
		else if(strVar == BRANCH_PREFIX"freq") m_dfreq = tl::str_to_var<decltype(m_dfreq)>(strVal);
		else if(strVar == BRANCH_PREFIX"E_HWHM") m_dHWHM = tl::str_to_var<decltype(m_dHWHM)>(strVal);
		else if(strVar == BRANCH_PREFIX"S0") m_dS0 = tl::str_to_var<decltype(m_dS0)>(strVal);

		else if(strVar == "inc_amp") m_dIncAmp = tl::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sig") m_dIncSig = tl::str_to_var<decltype(m_dIncSig)>(strVal);

		else if(strVar == "T") m_dT = tl::str_to_var<decltype(m_dT)>(strVal);
	}
}


SqwBase* SqwPhononSingleBranch::shallow_copy() const
{
	SqwPhononSingleBranch *pCpy = new SqwPhononSingleBranch();
	*static_cast<SqwBase*>(pCpy) = *static_cast<const SqwBase*>(this);

	pCpy->m_vecBragg = m_vecBragg;

	pCpy->m_damp = m_damp;
	pCpy->m_dfreq = m_dfreq;
	pCpy->m_dHWHM = m_dHWHM;
	pCpy->m_dS0 = m_dS0;

	pCpy->m_dIncAmp = m_dIncAmp;
	pCpy->m_dIncSig = m_dIncSig;

	pCpy->m_dT = m_dT;
	return pCpy;
}
