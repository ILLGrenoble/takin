/**
 * Convolution fitting model
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __CONVOFIT_MOD_H__
#define __CONVOFIT_MOD_H__

#include <memory>
#include <vector>
#include <string>

#include <boost/optional.hpp>

#include "tlibs/fit/minuit.h"
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnPrint.h>

#include <boost/signals2.hpp>

#include "../monteconvo/sqwbase.h"
#include "../monteconvo/TASReso.h"
#include "../res/defs.h"
#include "scan.h"

namespace minuit = ROOT::Minuit2;
namespace sig = boost::signals2;

//using t_real_mod = tl::t_real_min;
using t_real_mod = t_real_reso;


class SqwFuncModel : public tl::MinuitMultiFuncModel<t_real_mod>
{
protected:
	std::shared_ptr<SqwBase> m_pSqw;
	std::vector<TASReso> m_vecResos;
	std::vector<std::string> m_vecSqwParams;
	unsigned int m_iNumNeutrons = 1000;
	bool m_bUseThreadedMC = true;
	boost::optional<unsigned int> m_iSeed;

	ublas::vector<t_real_mod> m_vecScanOrigin;	// hklE
	ublas::vector<t_real_mod> m_vecScanDir;		// hklE
	t_real_mod m_dPrincipalAxisMin, m_dPrincipalAxisMax;

	// parameters that are not part of the S(Q, E) model
	t_real_mod m_dScale = 1., m_dSlope = 0., m_dOffs = 0.;
	t_real_mod m_dScaleErr = 0.1, m_dSlopeErr = 0., m_dOffsErr = 0.;
	std::vector<std::string> m_nonSQEParamNames = { "scale", "slope", "offs" };

	std::vector<std::string> m_vecModelParamNames;
	std::vector<t_real_mod> m_vecModelParams;
	std::vector<t_real_mod> m_vecModelErrs;

	std::string m_strTempParamName = "T";
	std::string m_strFieldParamName = "";

	// -------------------------------------------------------------------------
	// optional, for multi-fits
	std::size_t m_iCurParamSet = 0;
	const std::vector<Scan>* m_pScans = nullptr;
	// -------------------------------------------------------------------------


protected:
	// -------------------------------------------------------------------------
	// signals
	public: using t_sigFuncResult = sig::signal<void(t_real_mod h, t_real_mod k, t_real_mod l,
		t_real_mod E, t_real_mod S, std::size_t scan_group)>;
	protected: std::shared_ptr<t_sigFuncResult> m_psigFuncResult;
	public: void AddFuncResultSlot(const t_sigFuncResult::slot_type& slot)
	{
		if(!m_psigFuncResult) m_psigFuncResult = std::make_shared<t_sigFuncResult>();
		m_psigFuncResult->connect(slot);
	}

	public: using t_sigParamsChanged = sig::signal<void(const std::string&)>;
	protected: std::shared_ptr<t_sigParamsChanged> m_psigParamsChanged;
	public: void AddParamsChangedSlot(const t_sigParamsChanged::slot_type& slot)
	{
		if(!m_psigParamsChanged) m_psigParamsChanged = std::make_shared<t_sigParamsChanged>();
		m_psigParamsChanged->connect(slot);
	}
	// -------------------------------------------------------------------------


protected:
	void SetModelParams();

	bool SetTASPos(t_real_mod dX, TASReso& reso) const;
	TASReso* GetTASReso();
	const TASReso* GetTASReso() const;

	std::size_t GetNonSQEParamIdx(const std::string& param) const;


public:
	SqwFuncModel(std::shared_ptr<SqwBase> pSqw, const TASReso& reso);
	SqwFuncModel(std::shared_ptr<SqwBase> pSqw, const std::vector<TASReso>& vecResos);
	SqwFuncModel() = delete;
	virtual ~SqwFuncModel() = default;

	virtual bool SetParams(const std::vector<tl::t_real_min>& vecParams) override;
	virtual bool SetErrs(const std::vector<tl::t_real_min>& vecErrs);
	virtual tl::t_real_min operator()(tl::t_real_min x) const override;

	virtual SqwFuncModel* copy() const override;

	virtual const char* GetModelName() const override { return "SqwFuncModel"; }
	virtual std::vector<std::string> GetParamNames() const override;
	virtual std::vector<tl::t_real_min> GetParamValues() const override;
	virtual std::vector<tl::t_real_min> GetParamErrors() const override;

	// -------------------------------------------------------------------------
	// optional, for multi-fits
	virtual void SetParamSet(std::size_t iSet) override;
	virtual std::size_t GetParamSetCount() const override;
	virtual std::size_t GetExpLen() const override;
	virtual const t_real_mod* GetExpX() const override;
	virtual const t_real_mod* GetExpY() const override;
	virtual const t_real_mod* GetExpDY() const override;

	void SetScans(const std::vector<Scan>* pScans) { m_pScans = pScans; }
	// -------------------------------------------------------------------------


	void SetNonSQEParams(const std::vector<std::string>& nonSQEParams);
	void SetScale(t_real_mod scale) { m_dScale = scale; }
	void SetSlope(t_real_mod slope) { m_dSlope = slope; }
	void SetOffs(t_real_mod offs) { m_dOffs = offs; }

	void SetOtherParamNames(std::string strTemp, std::string strField);
	void SetOtherParams(t_real_mod dTemperature, t_real_mod dField);

	void SetReso(const TASReso& reso) { m_vecResos = { reso }; }
	void SetResos(const std::vector<TASReso>& vecResos) { m_vecResos = vecResos; }
	void SetSqwParamOverrides(const std::vector<std::string>& params) { m_vecSqwParams = params; }
	void SetNumNeutrons(unsigned int iNum) { m_iNumNeutrons = iNum; }
	void SetThreadedMC(bool b) { m_bUseThreadedMC = b; }
	void SetSeed(unsigned int iSeed) { m_iSeed = iSeed; }

	void SetScanOrigin(t_real_mod h, t_real_mod k, t_real_mod l, t_real_mod E)
	{ m_vecScanOrigin = tl::make_vec({h,k,l,E}); }
	void SetScanDir(t_real_mod h, t_real_mod k, t_real_mod l, t_real_mod E)
	{ m_vecScanDir = tl::make_vec({h,k,l,E}); }
	void SetPrincipalScanAxisMinMax(t_real_mod dMin, t_real_mod dMax)
	{ m_dPrincipalAxisMin = dMin; m_dPrincipalAxisMax = dMax; }

	void AddModelFitParams(const std::string& strName, t_real_mod dInitValue=0., t_real_mod dErr=0.)
	{
		m_vecModelParamNames.push_back(strName);
		m_vecModelParams.push_back(dInitValue);
		m_vecModelErrs.push_back(dErr);
	}

	minuit::MnUserParameters GetMinuitParams() const;
	void SetMinuitParams(const minuit::MnUserParameters& params);
	void SetMinuitParams(const minuit::MnUserParameterState& state)
	{ SetMinuitParams(state.Parameters()); }

	bool Save(const char *pcFile, std::size_t iPts = 256,
		std::size_t iSkipBegin = 0, std::size_t iSkipEnd = 0) const;

	SqwBase* GetSqwBase() { return m_pSqw.get(); }
};


#endif
