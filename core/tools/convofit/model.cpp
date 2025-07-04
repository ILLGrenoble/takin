/**
 * Convolution fitting model
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include <fstream>

#include "model.h"
#include "tlibs/math/math.h"
#include "tlibs/math/stat.h"
#include "tlibs/log/log.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/array.h"
#include "libs/globals.h"
#include "../res/defs.h"
#include "../res/ellipse.h"
#include "convofit.h"


using t_real = t_real_mod;


SqwFuncModel::SqwFuncModel(std::shared_ptr<SqwBase> pSqw, const TASReso& reso)
	: m_pSqw(pSqw), m_vecResos({reso})
{}


SqwFuncModel::SqwFuncModel(std::shared_ptr<SqwBase> pSqw, const std::vector<TASReso>& vecResos)
	: m_pSqw(pSqw), m_vecResos(vecResos)
{}


TASReso* SqwFuncModel::GetTASReso()
{
	TASReso *pReso = nullptr;
	// multi-fits
	if(m_pScans && m_vecResos.size() > 1)
		pReso = &m_vecResos[m_iCurParamSet];
	else
		pReso = &m_vecResos[0];
	return pReso;
}

const TASReso* SqwFuncModel::GetTASReso() const
{
	return const_cast<SqwFuncModel*>(this)->GetTASReso();
}


bool SqwFuncModel::SetTASPos(t_real dPrincipalX, TASReso& reso) const
{
	const t_real xrange = t_real(m_dPrincipalAxisMax - m_dPrincipalAxisMin);
	const t_real xscale = (t_real(dPrincipalX) - t_real(m_dPrincipalAxisMin)) / xrange;

	const ublas::vector<t_real> vecScanPos = m_vecScanOrigin + xscale*m_vecScanDir;

	if(!reso.SetHKLE(vecScanPos[0], vecScanPos[1], vecScanPos[2], vecScanPos[3]))
	{
		std::ostringstream ostrErr;
		ostrErr << "Invalid crystal position: ("
			<< vecScanPos[0] << " " << vecScanPos[1] << " " << vecScanPos[2]
			<< ") rlu, " << vecScanPos[3] << " meV.";
		tl::log_err(ostrErr.str());
		return false;
	}

	return true;
}


tl::t_real_min SqwFuncModel::operator()(tl::t_real_min x_principal) const
{
	TASReso/*&*/ reso = *GetTASReso();
	if(!SetTASPos(t_real_mod(x_principal), reso))
		return 0.;

	const t_real xrange = t_real(m_dPrincipalAxisMax - m_dPrincipalAxisMin);
	const t_real xscale = (t_real(x_principal) - t_real(m_dPrincipalAxisMin)) / xrange;
	const ublas::vector<t_real> vecScanPos = m_vecScanOrigin + t_real(xscale)*m_vecScanDir;

	std::vector<ublas::vector<t_real_reso>> vecNeutrons;
	Ellipsoid4d<t_real_reso> elli;
	if(m_bUseThreadedMC)
	{
		elli = reso.GenerateMC(m_iNumNeutrons, vecNeutrons);
	}
	else
	{
		if(m_iSeed)
		{
			// reset the seed here in case chi^2 runs multi-threaded and neutron recycling is enabled
			tl::init_rand_seed(*m_iSeed, false);
		}
		elli = reso.GenerateMC_deferred(m_iNumNeutrons, vecNeutrons);
	}

	t_real dS = 0.;
	t_real dhklE_mean[4] = { 0., 0., 0., 0. };

	for(const ublas::vector<t_real_reso>& vecHKLE : vecNeutrons)
	{
		dS += t_real((*m_pSqw)(vecHKLE[0], vecHKLE[1], vecHKLE[2], vecHKLE[3]));

		for(int i = 0; i < 4; ++i)
			dhklE_mean[i] += t_real(vecHKLE[i]);
	}

	dS /= t_real(m_iNumNeutrons);
	dS += m_pSqw->GetBackground(vecScanPos[0], vecScanPos[1], vecScanPos[2], vecScanPos[3]);

	for(int i = 0; i < 4; ++i)
		dhklE_mean[i] /= t_real(m_iNumNeutrons);

	dS *= reso.GetResoResults().dR0 * reso.GetR0Scale();
	//if(reso.GetResoParams().flags & CALC_RESVOL)
	//	dS /= reso.GetResoResults().dResVol * tl::get_pi<t_real>() * t_real(3.);

	t_real dYVal = m_dScale*(dS + m_dSlope*x_principal) + m_dOffs;
	if(dYVal < 0.)
		dYVal = 0.;

	if(m_psigFuncResult)
	{
		(*m_psigFuncResult)(vecScanPos[0], vecScanPos[1], vecScanPos[2], vecScanPos[3],
			dYVal, m_iCurParamSet);
	}
	return tl::t_real_min(dYVal);
}


SqwFuncModel* SqwFuncModel::copy() const
{
	// cannot rebuild kd tree in phonon model with only a shallow copy
	SqwFuncModel* pMod = new SqwFuncModel(
		std::shared_ptr<SqwBase>(m_pSqw->shallow_copy()), m_vecResos);

	pMod->m_vecScanOrigin = this->m_vecScanOrigin;
	pMod->m_vecScanDir = this->m_vecScanDir;
	pMod->m_dPrincipalAxisMin = this->m_dPrincipalAxisMin;
	pMod->m_dPrincipalAxisMax = this->m_dPrincipalAxisMax;

	pMod->m_iNumNeutrons = this->m_iNumNeutrons;
	pMod->m_bUseThreadedMC = this->m_bUseThreadedMC;
	pMod->m_iSeed = this->m_iSeed;

	pMod->m_dScale = this->m_dScale;
	pMod->m_dSlope = this->m_dSlope;
	pMod->m_dOffs = this->m_dOffs;
	pMod->m_dScaleErr = this->m_dScaleErr;
	pMod->m_dSlopeErr = this->m_dSlopeErr;
	pMod->m_dOffsErr = this->m_dOffsErr;
	pMod->m_nonSQEParamNames = this->m_nonSQEParamNames;

	pMod->m_vecModelParamNames = this->m_vecModelParamNames;
	pMod->m_vecModelParams = this->m_vecModelParams;
	pMod->m_vecModelErrs = this->m_vecModelErrs;

	pMod->m_strTempParamName = this->m_strTempParamName;
	pMod->m_strFieldParamName = this->m_strFieldParamName;

	pMod->m_iCurParamSet = this->m_iCurParamSet;
	pMod->m_pScans = this->m_pScans;
	pMod->m_psigFuncResult = this->m_psigFuncResult;
	pMod->m_psigParamsChanged = this->m_psigParamsChanged;
	pMod->m_vecSqwParams = this->m_vecSqwParams;

	return pMod;
}


void SqwFuncModel::SetNonSQEParams(const std::vector<std::string>& nonSQEParams)
{
	m_nonSQEParamNames = nonSQEParams;
}


void SqwFuncModel::SetOtherParamNames(std::string strTemp, std::string strField)
{
	m_strTempParamName = strTemp;
	m_strFieldParamName = strField;
}


void SqwFuncModel::SetOtherParams(t_real dTemperature, t_real dField)
{
	std::vector<SqwBase::t_var> vecVars;
	if(m_strTempParamName != "")
		vecVars.push_back(std::make_tuple(m_strTempParamName, "double", tl::var_to_str(dTemperature, g_iPrec)));
	if(m_strFieldParamName != "")
		vecVars.push_back(std::make_tuple(m_strFieldParamName, "double", tl::var_to_str(dField, g_iPrec)));
	m_pSqw->SetVars(vecVars);
}


void SqwFuncModel::SetModelParams()
{
	const std::size_t iNumParams = m_vecModelParams.size();
	std::vector<SqwBase::t_var> vecVars;
	vecVars.reserve(iNumParams);

	for(std::size_t iParam = 0; iParam < iNumParams; ++iParam)
	{
		std::string strVal = tl::var_to_str(m_vecModelParams[iParam], g_iPrec);
		SqwBase::t_var var = std::make_tuple(m_vecModelParamNames[iParam], "double", strVal);
		vecVars.push_back(var);
	}

	m_pSqw->SetVars(vecVars);
}


std::size_t SqwFuncModel::GetNonSQEParamIdx(const std::string& param) const
{
	for(std::size_t i = 0; i < m_nonSQEParamNames.size(); ++i)
		if(m_nonSQEParamNames[i] == param)
			return i;

	// not found
	return m_nonSQEParamNames.size();
}


bool SqwFuncModel::SetParams(const std::vector<tl::t_real_min>& vecParams)
{
	// --------------------------------------------------------------------
	// prints changed model parameters
	std::vector<t_real> vecOldParams;

	// add parameters that are not part of the S(Q, E) model
	for(const std::string& param_name : m_nonSQEParamNames)
	{
		if(param_name == "scale")
			vecOldParams.push_back(m_dScale);
		else if(param_name == "slope")
			vecOldParams.push_back(m_dSlope);
		else if(param_name == "offs")
			vecOldParams.push_back(m_dOffs);
	}

	// add S(Q, E) model parameters
	vecOldParams.insert(vecOldParams.end(), m_vecModelParams.begin(), m_vecModelParams.end());

	std::vector<std::string> vecParamNames = GetParamNames();

	if(vecOldParams.size() == vecParams.size() && vecParamNames.size() == vecParams.size())
	{
		std::ostringstream ostrDebug;
		std::transform(vecParams.begin(), vecParams.end(), vecParamNames.begin(),
			std::ostream_iterator<std::string>(ostrDebug, ", "),
			[&vecOldParams, &vecParamNames](tl::t_real_min _dVal, const std::string& strParam) -> std::string
			{
				t_real dVal = t_real(_dVal);
				std::vector<std::string>::const_iterator iterParam =
					std::find(vecParamNames.begin(), vecParamNames.end(), strParam);
				if(iterParam == vecParamNames.end())
					return "";
				t_real dOldParam = vecOldParams[iterParam-vecParamNames.begin()];

				bool bChanged = !tl::float_equal(dVal, dOldParam);
				std::string strRet = strParam +
					std::string(" = ") +
					tl::var_to_str(dVal, g_iPrec);
				if(bChanged)
					strRet += " (old: " + tl::var_to_str(dOldParam, g_iPrec) + ")";
				return strRet;
			});

		if(m_psigParamsChanged)
			(*m_psigParamsChanged)(ostrDebug.str());
	}
	else
	{
		tl::log_err("Parameter size mismatch.");
	}
	// --------------------------------------------------------------------

	if(vecParams.size() < m_nonSQEParamNames.size())
	{
		tl::log_err("Invalid size of model parameters. Has to have at least three parameters: scale, slope, offs.");
		return false;
	}

	// parameters that are not part of the S(Q, E) model
	std::size_t idxScale = GetNonSQEParamIdx("scale");
	std::size_t idxSlope = GetNonSQEParamIdx("slope");
	std::size_t idxOffs = GetNonSQEParamIdx("offs");

	if(idxScale < m_nonSQEParamNames.size())
		m_dScale = t_real(vecParams[idxScale]);
	if(idxSlope < m_nonSQEParamNames.size())
		m_dSlope = t_real(vecParams[idxSlope]);
	if(idxOffs < m_nonSQEParamNames.size())
		m_dOffs = t_real(vecParams[idxOffs]);

	for(std::size_t iParam = m_nonSQEParamNames.size(); iParam < vecParams.size(); ++iParam)
		m_vecModelParams[iParam - m_nonSQEParamNames.size()] = t_real(vecParams[iParam]);

	SetModelParams();
	return true;
}


bool SqwFuncModel::SetErrs(const std::vector<tl::t_real_min>& vecErrs)
{
	if(vecErrs.size() < m_nonSQEParamNames.size())
	{
		tl::log_err("Invalid size of model parameter errors. Has to have at least three parameters: scale, slope, offs.");
		return false;
	}

	// parameters that are not part of the S(Q, E) model
	std::size_t idxScale = GetNonSQEParamIdx("scale");
	std::size_t idxSlope = GetNonSQEParamIdx("slope");
	std::size_t idxOffs = GetNonSQEParamIdx("offs");

	if(idxScale < m_nonSQEParamNames.size())
		m_dScaleErr = t_real(vecErrs[idxScale]);
	if(idxSlope < m_nonSQEParamNames.size())
		m_dSlopeErr = t_real(vecErrs[idxSlope]);
	if(idxOffs < m_nonSQEParamNames.size())
		m_dOffsErr = t_real(vecErrs[idxOffs]);

	for(std::size_t iParam = m_nonSQEParamNames.size(); iParam<vecErrs.size(); ++iParam)
		m_vecModelErrs[iParam - m_nonSQEParamNames.size()] = t_real(vecErrs[iParam]);

	//SetModelParams();
	return true;
}


std::vector<std::string> SqwFuncModel::GetParamNames() const
{
	// parameters that are not part of the S(Q, E) model
	std::vector<std::string> vecNames = m_nonSQEParamNames;

	for(const std::string& str : m_vecModelParamNames)
		vecNames.push_back(str);

	return vecNames;
}


std::vector<tl::t_real_min> SqwFuncModel::GetParamValues() const
{
	std::vector<tl::t_real_min> vecVals;

	// add parameters that are not part of the S(Q, E) model
	for(const std::string& param_name : m_nonSQEParamNames)
	{
		if(param_name == "scale")
			vecVals.push_back(m_dScale);
		else if(param_name == "slope")
			vecVals.push_back(m_dSlope);
		else if(param_name == "offs")
			vecVals.push_back(m_dOffs);
	}

	// add S(Q, E) model parameters
	for(t_real d : m_vecModelParams)
		vecVals.push_back(tl::t_real_min(d));

	return vecVals;
}


std::vector<tl::t_real_min> SqwFuncModel::GetParamErrors() const
{
	std::vector<tl::t_real_min> vecErrs;

	// add parameter errors that are not part of the S(Q, E) model
	for(const std::string& param_name : m_nonSQEParamNames)
	{
		if(param_name == "scale")
			vecErrs.push_back(m_dScaleErr);
		else if(param_name == "slope")
			vecErrs.push_back(m_dSlopeErr);
		else if(param_name == "offs")
			vecErrs.push_back(m_dOffsErr);
	}

	// add S(Q, E) model parameter errors
	for(t_real d : m_vecModelErrs)
		vecErrs.push_back(tl::t_real_min(d));

	return vecErrs;
}


void SqwFuncModel::SetMinuitParams(const minuit::MnUserParameters& state)
{
	std::vector<t_real> vecNewVals;
	std::vector<t_real> vecNewErrs;

	const std::vector<std::string> vecNames = GetParamNames();
	for(std::size_t iParam = 0; iParam < vecNames.size(); ++iParam)
	{
		const std::string& strName = vecNames[iParam];

		const t_real dVal = t_real(state.Value(strName));
		const t_real dErr = t_real(state.Error(strName));

		vecNewVals.push_back(dVal);
		vecNewErrs.push_back(dErr);
	}

	SetParams(tl::container_cast<tl::t_real_min, t_real, std::vector>()(vecNewVals));
	SetErrs(tl::container_cast<tl::t_real_min, t_real, std::vector>()(vecNewErrs));
}


minuit::MnUserParameters SqwFuncModel::GetMinuitParams() const
{
	minuit::MnUserParameters params;

	// parameters that are not part of the S(Q, E) model
	for(const std::string& paramname : m_nonSQEParamNames)
	{
		if(paramname == "scale")
			params.Add("scale", m_dScale, m_dScaleErr);
		if(paramname == "slope")
			params.Add("slope", m_dSlope, m_dSlopeErr);
		if(paramname == "offs")
			params.Add("offs", m_dOffs, m_dOffsErr);
	}

	for(std::size_t iParam = 0; iParam < m_vecModelParamNames.size(); ++iParam)
	{
		const std::string& strParam = m_vecModelParamNames[iParam];
		tl::t_real_min dHint = tl::t_real_min(m_vecModelParams[iParam]);
		tl::t_real_min dErr = tl::t_real_min(m_vecModelErrs[iParam]);

		params.Add(strParam, dHint, dErr);
	}

	return params;
}


bool SqwFuncModel::Save(const char *pcFile, std::size_t iNum,
	std::size_t iSkipBegin, std::size_t iSkipEnd, const char* comment) const
{
	if(iSkipBegin + iSkipEnd >= iNum)
	{
		tl::log_err("No points to plot");
		return 0;
	}

	try
	{
		std::ofstream ofstr;
		//ofstr.exceptions(ofstr.exceptions() | std::ios_base::failbit | std::ios_base::badbit);
		ofstr.open(pcFile, std::ios_base::out | std::ios_base::trunc);
		if(!ofstr)
		{
			tl::log_err("Cannot open model file \"", pcFile, "\" for writing.");
			return false;
		}

		ofstr.precision(g_iPrec);
		const std::vector<std::string> vecNames = GetParamNames();

		tl::container_cast<t_real, tl::t_real_min, std::vector> cst;
		const std::vector<t_real> vecVals = cst(GetParamValues());
		const std::vector<t_real> vecErrs = cst(GetParamErrors());

		ofstr << "#\n";
		for(std::size_t iParam = 0; iParam < vecNames.size(); ++iParam)
		{
			ofstr << "# " << vecNames[iParam] << " = "
				<< vecVals[iParam] << " +- "
				<< vecErrs[iParam] << "\n";
		}

		ofstr << "#\n";
		if(comment)
			ofstr << "# " << comment << "\n#\n";

		ofstr << "# Data columns: (1) scan axis, (2) intensity";
		ofstr << ", (3) Bragg Qx (rlu), (4) Bragg Qy (rlu), (5) Bragg Qz (rlu), (6) Bragg E (meV)\n";
		ofstr << "#\n";

		std::vector<t_real> vecFWHMs[4];

		tl::log_info("Number of plot points: ", iNum, ", skip: ", iSkipBegin, ", ", iSkipEnd, ".");

		for(std::size_t i = iSkipBegin; i < iNum - iSkipEnd; ++i)
		{
			t_real dX = tl::lerp(t_real(m_dPrincipalAxisMin), t_real(m_dPrincipalAxisMax), t_real(i)/t_real(iNum-1));
			t_real dY = (*this)(dX);

			ofstr << std::left << std::setw(g_iPrec*2.5) << dX << " "
				<< std::left << std::setw(g_iPrec*2.5) << dY << " ";

			// save bragg widths for error calculation
			// TODO: also save bragg width in rlu
			TASReso reso = *GetTASReso();
			SetTASPos(dX, reso);
			const auto& crysopts = reso.GetMCOpts();
			const ResoResults& resores = reso.GetResoResults();

			ublas::matrix<t_real> resoHKL;
			std::tie(resoHKL, std::ignore, std::ignore) =
			conv_lab_to_rlu<ublas::matrix<t_real>, ublas::vector<t_real>>
				(crysopts.dAngleQVec0, crysopts.matUB, crysopts.matUBinv,
				 resores.reso, resores.reso_v, resores.Q_avg);

			const std::vector<t_real> vecFwhms = calc_bragg_fwhms(resoHKL);
			for(int iFwhm = 0; iFwhm < 4; ++iFwhm)
				vecFWHMs[iFwhm].push_back(vecFwhms[iFwhm]);

			for(t_real dFwhm : vecFwhms)
				ofstr << std::left << std::setw(g_iPrec*2.5) << dFwhm << " ";
			ofstr << "\n";

			ofstr.flush();
		}

		for(int iFwhm = 0; iFwhm < 4; ++iFwhm)
		{
			t_real dSig = tl::mean_value(vecFWHMs[iFwhm]);
			t_real dDSig = tl::std_dev(vecFWHMs[iFwhm]);
			ofstr << "# " << "bragg_sig_" << iFwhm << " = "
				<< dSig*tl::get_FWHM2SIGMA<t_real>() << " +- "
				<< dDSig*tl::get_FWHM2SIGMA<t_real>() << "\n";
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err("Saving model failed: ", ex.what());
		return false;
	}

	return true;
}



// -----------------------------------------------------------------------------
// optional, for multi-fits
void SqwFuncModel::SetParamSet(std::size_t iSet)
{
	if(!m_pScans)
		return;

	if(iSet >= m_pScans->size())
	{
		tl::log_err("Requested invalid scan group ", iSet, ".");
		return;
	}

	m_iCurParamSet = iSet;

	// parameters that are not part of the S(Q, E) model
	static const std::unordered_set<std::string> ignored_parms{{ "scale", "slope", "offs" }};

	// set S(Q, E) parameters from scan file
	const Scan& sc = m_pScans->operator[](m_iCurParamSet);
	if(m_vecResos.size() > 1)
		set_tasreso_params_from_scan(m_vecResos[m_iCurParamSet], sc);
	else
		set_tasreso_params_from_scan(m_vecResos[0], sc);
	set_model_params_from_scan(*this, sc);


	// set S(Q, E) parameters from optional overrides
	std::string strParams;
	if(m_iCurParamSet < m_vecSqwParams.size())
	{
		bool params_ok = false;
		std::unordered_map<std::string, std::string> all_params;

		std::tie(params_ok, strParams) =
			m_pSqw->SetVars(m_vecSqwParams[m_iCurParamSet],
				false, &ignored_parms, &all_params);


		// change parameters that are not part of the S(Q, E) model
		auto iter_scale = all_params.find("scale");
		auto iter_slope = all_params.find("slope");
		auto iter_offs = all_params.find("offs");

		if(iter_scale != all_params.end())
		{
			m_dScale = tl::str_to_var<t_real>(iter_scale->second);

			if(strParams.length())
				strParams += ", ";
			strParams += "scale = " + tl::var_to_str(m_dScale, g_iPrec);
		}

		if(iter_slope != all_params.end())
		{
			m_dSlope = tl::str_to_var<t_real>(iter_slope->second);

			if(strParams.length())
				strParams += ", ";
			strParams += "slope = " + tl::var_to_str(m_dSlope, g_iPrec);
		}

		if(iter_offs != all_params.end())
		{
			m_dOffs = tl::str_to_var<t_real>(iter_offs->second);

			if(strParams.length())
				strParams += ", ";
			strParams += "offs = " + tl::var_to_str(m_dOffs, g_iPrec);
		}


		if(!params_ok)
		{
			tl::log_err("Could not override S(Q, E) model parameter(s) for scan group ",
				m_iCurParamSet, ".");
		}
	}

	if(m_psigParamsChanged)
	{
		std::ostringstream ostrDescr;
		ostrDescr << "Scan group " << m_iCurParamSet;
		if(strParams.length())
			ostrDescr << ": " << strParams;
		ostrDescr << ".";

		(*m_psigParamsChanged)(ostrDescr.str());
	}
}


std::size_t SqwFuncModel::GetParamSetCount() const
{
	// multi-fits
	if(m_pScans)
		return m_pScans->size();

	// single-fits
	return 1;
}


std::size_t SqwFuncModel::GetExpLen() const
{
	if(m_pScans)
		return m_pScans->operator[](m_iCurParamSet).vecX.size();
	return 0;
}


const t_real_mod* SqwFuncModel::GetExpX() const
{
	if(m_pScans)
		return m_pScans->operator[](m_iCurParamSet).vecX.data();
	return nullptr;
}


const t_real_mod* SqwFuncModel::GetExpY() const
{
	if(m_pScans)
		return m_pScans->operator[](m_iCurParamSet).vecCts.data();
	return nullptr;
}


const t_real_mod* SqwFuncModel::GetExpDY() const
{
	if(m_pScans)
		return m_pScans->operator[](m_iCurParamSet).vecCtsErr.data();
	return nullptr;
}
// -----------------------------------------------------------------------------
