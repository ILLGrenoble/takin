/**
 * scan viewer -- curve fitting
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2023
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

#include "scanviewer.h"

#include <QMessageBox>

#include <string>
#include <algorithm>

#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/stat.h"
#include "tlibs/log/log.h"
#include "tlibs/helper/misc.h"

#include "tlibs/fit/minuit.h"
#include "tlibs/fit/interpolation.h"
#include "tlibs/fit/swarm.h"

using tl::t_real_min;
using t_real = t_real_glob;


void ScanViewerDlg::ShowFitParams()
{
	focus_dlg(m_pFitParamDlg);
}


/**
 * fit a function to data points
 */
template<std::size_t iFuncArgs, class t_func>
bool ScanViewerDlg::Fit(t_func&& func,
	const std::vector<std::string>& vecParamNames,
	std::vector<t_real>& vecVals,
	std::vector<t_real>& vecErrs,
	const std::vector<bool>& vecFixed)
{
	bool bUseSwarm = m_settings.value("use_swarm", false).toBool();

	m_vecFitX.clear();
	m_vecFitY.clear();

	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return false;

	bool bOk = 0;
	try
	{
		std::vector<t_real_min>
			_vecVals = tl::container_cast<t_real_min, t_real, std::vector>()(vecVals),
			_vecErrs = tl::container_cast<t_real_min, t_real, std::vector>()(vecErrs);

		if(bUseSwarm)
		{
			bOk = tl::swarmfit<t_real, iFuncArgs>
				(func, m_vecX, m_vecY, m_vecYErr, vecParamNames, vecVals, vecErrs);
		}

		bOk = tl::fit<iFuncArgs>(func,
			tl::container_cast<t_real_min, t_real, std::vector>()(m_vecX),
			tl::container_cast<t_real_min, t_real, std::vector>()(m_vecY),
			tl::container_cast<t_real_min, t_real, std::vector>()(m_vecYErr),
			vecParamNames, _vecVals, _vecErrs, &vecFixed);
		vecVals = tl::container_cast<t_real, t_real_min, std::vector>()(_vecVals);
		vecErrs = tl::container_cast<t_real, t_real_min, std::vector>()(_vecErrs);
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	if(!bOk)
	{
		QMessageBox::critical(this, "Error", "Could not fit function. Please set or improve the initial parameters.");
		return false;
	}

	auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());

	m_vecFitX.reserve(GFX_NUM_POINTS);
	m_vecFitY.reserve(GFX_NUM_POINTS);

	std::vector<t_real> vecValsWithX = { 0. };
	for(t_real dVal : vecVals) vecValsWithX.push_back(dVal);

	for(std::size_t i=0; i<GFX_NUM_POINTS; ++i)
	{
		t_real dX = t_real(i)*(*minmaxX.second - *minmaxX.first)/t_real(GFX_NUM_POINTS-1) + *minmaxX.first;
		vecValsWithX[0] = dX;
		t_real dY = tl::call<iFuncArgs, decltype(func), t_real, std::vector>(func, vecValsWithX);

		m_vecFitX.push_back(dX);
		m_vecFitY.push_back(dY);
	}

	PlotScan();
	m_pFitParamDlg->UnsetAllBold();
	return true;
}


void ScanViewerDlg::FitLine()
{
	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return;

	auto func = [](t_real x, t_real m, t_real offs) -> t_real { return m*x + offs; };
	constexpr std::size_t iFuncArgs = 3;

	t_real_glob dSlope = m_pFitParamDlg->GetSlope(),	dSlopeErr = m_pFitParamDlg->GetSlopeErr();
	t_real_glob dOffs = m_pFitParamDlg->GetOffs(),		dOffsErr = m_pFitParamDlg->GetOffsErr();

	bool bSlopeFixed = m_pFitParamDlg->GetSlopeFixed();
	bool bOffsFixed = m_pFitParamDlg->GetOffsFixed();

	// automatic parameter determination
	if(!m_pFitParamDlg->WantParams())
	{
		auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
		auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());

		dSlope = (*minmaxY.second - *minmaxY.first) / (*minmaxX.second - *minmaxX.first);
		dOffs = *minmaxY.first;

		dSlopeErr = dSlope * 0.1;
		dOffsErr = dOffs * 0.1;

		bSlopeFixed = bOffsFixed = 0;
	}

	std::vector<std::string> vecParamNames = { "slope", "offs" };
	std::vector<t_real> vecVals = { dSlope, dOffs };
	std::vector<t_real> vecErrs = { dSlopeErr, dOffsErr };
	std::vector<bool> vecFixed = { bSlopeFixed, bOffsFixed };

	if(!Fit<iFuncArgs>(func, vecParamNames, vecVals, vecErrs, vecFixed))
		return;

	for(t_real &d : vecErrs)
		d = std::abs(d);

	m_pFitParamDlg->SetSlope(vecVals[0]);	m_pFitParamDlg->SetSlopeErr(vecErrs[0]);
	m_pFitParamDlg->SetOffs(vecVals[1]);	m_pFitParamDlg->SetOffsErr(vecErrs[1]);
}


/**
 * parabola: y = amp*(x-x0)^2 + offs
 */
void ScanViewerDlg::FitParabola()
{
	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return;

	const bool bUseSlope = checkSloped->isChecked();

	auto func = tl::parabola_model<t_real>;
	auto funcSloped = tl::parabola_model_slope<t_real>;
	constexpr std::size_t iFuncArgs = 4;
	constexpr std::size_t iFuncArgsSloped = iFuncArgs+1;

	t_real_glob dAmp = m_pFitParamDlg->GetAmp(),	dAmpErr = m_pFitParamDlg->GetAmpErr();
	t_real_glob dX0 = m_pFitParamDlg->GetX0(),	dX0Err = m_pFitParamDlg->GetX0Err();
	t_real_glob dOffs = m_pFitParamDlg->GetOffs(),	dOffsErr = m_pFitParamDlg->GetOffsErr();
	t_real_glob dSlope = m_pFitParamDlg->GetSlope(),dSlopeErr = m_pFitParamDlg->GetSlopeErr();

	bool bAmpFixed = m_pFitParamDlg->GetAmpFixed();
	bool bX0Fixed = m_pFitParamDlg->GetX0Fixed();
	bool bOffsFixed = m_pFitParamDlg->GetOffsFixed();
	bool bSlopeFixed = m_pFitParamDlg->GetSlopeFixed();

	// automatic parameter determination
	if(!m_pFitParamDlg->WantParams())
	{
		auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());

		dX0 = m_vecX[minmaxY.second - m_vecY.begin()];
		dAmp = std::abs(*minmaxY.second-*minmaxY.first);
		dOffs = *minmaxY.first;

		dX0Err = dX0 * 0.1;
		dAmpErr = dAmp * 0.1;
		dOffsErr = dOffs * 0.1;

		bAmpFixed = bX0Fixed = bOffsFixed = 0;
	}

	std::vector<std::string> vecParamNames = { "x0", "amp", "offs" };
	std::vector<t_real> vecVals = { dX0, dAmp, dOffs };
	std::vector<t_real> vecErrs = { dX0Err, dAmpErr, dOffsErr };
	std::vector<bool> vecFixed = { bX0Fixed, bAmpFixed, bOffsFixed };

	if(bUseSlope)
	{
		vecParamNames.push_back("slope");
		vecVals.push_back(dSlope);
		vecErrs.push_back(dSlopeErr);
		vecFixed.push_back(bSlopeFixed);
	}

	bool bOk = false;
	if(bUseSlope)
		bOk = Fit<iFuncArgsSloped>(funcSloped, vecParamNames, vecVals, vecErrs, vecFixed);
	else
		bOk = Fit<iFuncArgs>(func, vecParamNames, vecVals, vecErrs, vecFixed);
	if(!bOk)
		return;

	for(t_real &d : vecErrs) d = std::abs(d);
	vecVals[1] = std::abs(vecVals[1]);

	m_pFitParamDlg->SetX0(vecVals[0]);	m_pFitParamDlg->SetX0Err(vecErrs[0]);
	m_pFitParamDlg->SetAmp(vecVals[1]);	m_pFitParamDlg->SetAmpErr(vecErrs[1]);
	m_pFitParamDlg->SetOffs(vecVals[2]);	m_pFitParamDlg->SetOffsErr(vecErrs[2]);

	if(bUseSlope)
	{
		m_pFitParamDlg->SetSlope(vecVals[3]);
		m_pFitParamDlg->SetSlopeErr(vecErrs[3]);
	}
}


/**
 * sin(x + pi) = -sin(x)
 * sin(-x + phi) = -sin(x - phi)
 */
template<class t_real = double>
void sanitise_sine_params(t_real& amp, t_real& freq, t_real& phase, t_real& offs)
{
	if(freq < t_real(0))
	{
		freq = -freq;
		phase = -phase;
		amp = -amp;
	}

	if(amp < t_real(0))
	{
		amp = -amp;
		phase += tl::get_pi<t_real>();
	}

	if(phase < t_real(0))
	{
		int iNum = std::abs(int(phase / (t_real(2)*tl::get_pi<t_real>()))) + 1;
		phase += t_real(2*iNum) * tl::get_pi<t_real>();
	}

	phase = std::fmod(phase, t_real(2)*tl::get_pi<t_real>());
}


void ScanViewerDlg::FitSine()
{
	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return;

	const bool bUseSlope = checkSloped->isChecked();

	auto func = [](t_real x, t_real amp, t_real freq, t_real phase, t_real offs) -> t_real
		{ return amp*std::sin(freq*x + phase) + offs; };
	constexpr std::size_t iFuncArgs = 5;

	auto funcSloped = [](t_real x, t_real amp, t_real freq, t_real phase, t_real offs, t_real slope) -> t_real
	{ return amp*std::sin(freq*x + phase) + slope*x + offs; };
	constexpr std::size_t iFuncArgsSloped = 6;


	t_real_glob dAmp = m_pFitParamDlg->GetAmp(),	dAmpErr = m_pFitParamDlg->GetAmpErr();
	t_real_glob dFreq = m_pFitParamDlg->GetFreq(),	dFreqErr = m_pFitParamDlg->GetFreqErr();
	t_real_glob dPhase = m_pFitParamDlg->GetPhase(),dPhaseErr = m_pFitParamDlg->GetPhaseErr();
	t_real_glob dOffs = m_pFitParamDlg->GetOffs(),	dOffsErr = m_pFitParamDlg->GetOffsErr();
	t_real_glob dSlope = m_pFitParamDlg->GetSlope(),	dSlopeErr = m_pFitParamDlg->GetSlopeErr();

	bool bAmpFixed = m_pFitParamDlg->GetAmpFixed();
	bool bFreqFixed = m_pFitParamDlg->GetFreqFixed();
	bool bPhaseFixed = m_pFitParamDlg->GetPhaseFixed();
	bool bOffsFixed = m_pFitParamDlg->GetOffsFixed();
	bool bSlopeFixed = m_pFitParamDlg->GetSlopeFixed();

	// automatic parameter determination
	if(!m_pFitParamDlg->WantParams())
	{
		auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
		auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());

		dFreq = t_real(2.*M_PI) / (*minmaxX.second - *minmaxX.first);
		//dOffs = *minmaxY.first + (*minmaxY.second - *minmaxY.first)*0.5;
		dOffs = tl::mean_value(m_vecY);
		dAmp = (std::abs(*minmaxY.second - dOffs) + std::abs(dOffs - *minmaxY.first)) * 0.5;
		dPhase = 0.;

		dFreqErr = dFreq * 0.1;
		dOffsErr = tl::std_dev(m_vecY);;
		dAmpErr = dAmp * 0.1;
		dPhaseErr = M_PI;

		bAmpFixed = bFreqFixed = bPhaseFixed = bOffsFixed = 0;
	}

	std::vector<std::string> vecParamNames = { "amp", "freq", "phase",  "offs" };
	std::vector<t_real> vecVals = { dAmp, dFreq, dPhase, dOffs };
	std::vector<t_real> vecErrs = { dAmpErr, dFreqErr, dPhaseErr, dOffsErr };
	std::vector<bool> vecFixed = { bAmpFixed, bFreqFixed, bPhaseFixed, bOffsFixed };

	if(bUseSlope)
	{
		vecParamNames.push_back("slope");
		vecVals.push_back(dSlope);
		vecErrs.push_back(dSlopeErr);
		vecFixed.push_back(bSlopeFixed);
	}

	bool bOk = false;
	if(bUseSlope)
		bOk = Fit<iFuncArgsSloped>(funcSloped, vecParamNames, vecVals, vecErrs, vecFixed);
	else
		bOk = Fit<iFuncArgs>(func, vecParamNames, vecVals, vecErrs, vecFixed);
	if(!bOk)
		return;

	for(t_real &d : vecErrs)
		d = std::abs(d);

	sanitise_sine_params<t_real>(vecVals[0], vecVals[1], vecVals[2], vecVals[3]);


	m_pFitParamDlg->SetAmp(vecVals[0]);		m_pFitParamDlg->SetAmpErr(vecErrs[0]);
	m_pFitParamDlg->SetFreq(vecVals[1]);	m_pFitParamDlg->SetFreqErr(vecErrs[1]);
	m_pFitParamDlg->SetPhase(vecVals[2]);	m_pFitParamDlg->SetPhaseErr(vecErrs[2]);
	m_pFitParamDlg->SetOffs(vecVals[3]);	m_pFitParamDlg->SetOffsErr(vecErrs[3]);

	if(bUseSlope)
	{
		m_pFitParamDlg->SetSlope(vecVals[4]);
		m_pFitParamDlg->SetSlopeErr(vecErrs[4]);
	}
}


/**
 * fits a gaussian peak
 */
void ScanViewerDlg::FitGauss()
{
	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return;


	// prefit
	unsigned int iOrder = m_settings.value("spline_order", 6).toInt();
	std::vector<t_real> vecMaximaX, vecMaximaSize, vecMaximaWidth;
	tl::find_peaks<t_real>(m_vecX.size(), m_vecX.data(), m_vecY.data(), iOrder,
		vecMaximaX, vecMaximaSize, vecMaximaWidth, g_dEpsGfx);


	const bool bUseSlope = checkSloped->isChecked();

	auto func = tl::gauss_model_amp<t_real>;
	auto funcSloped = tl::gauss_model_amp_slope<t_real>;
	constexpr std::size_t iFuncArgs = 5;
	constexpr std::size_t iFuncArgsSloped = iFuncArgs+1;

	t_real_glob dAmp = m_pFitParamDlg->GetAmp(),	dAmpErr = m_pFitParamDlg->GetAmpErr();
	t_real_glob dSig = m_pFitParamDlg->GetSig(),	dSigErr = m_pFitParamDlg->GetSigErr();
	t_real_glob dX0 = m_pFitParamDlg->GetX0(),		dX0Err = m_pFitParamDlg->GetX0Err();
	t_real_glob dOffs = m_pFitParamDlg->GetOffs(),	dOffsErr = m_pFitParamDlg->GetOffsErr();
	t_real_glob dSlope = m_pFitParamDlg->GetSlope(),	dSlopeErr = m_pFitParamDlg->GetSlopeErr();

	bool bAmpFixed = m_pFitParamDlg->GetAmpFixed();
	bool bSigFixed = m_pFitParamDlg->GetSigFixed();
	bool bX0Fixed = m_pFitParamDlg->GetX0Fixed();
	bool bOffsFixed = m_pFitParamDlg->GetOffsFixed();
	bool bSlopeFixed = m_pFitParamDlg->GetSlopeFixed();

	// automatic parameter determination
	if(!m_pFitParamDlg->WantParams())
	{
		auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
		auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());

		// use prefitter values
		if(vecMaximaX.size())
		{
			dX0 = vecMaximaX[0];
			dAmp = vecMaximaSize[0];
			dSig = tl::get_FWHM2SIGMA<t_real>()*vecMaximaWidth[0] * 0.5;
			dOffs = *minmaxY.first;
		}
		// try to guess values
		else
		{
			dX0 = m_vecX[minmaxY.second - m_vecY.begin()];
			dSig = std::abs((*minmaxX.second-*minmaxX.first) * 0.5);
			dAmp = std::abs(*minmaxY.second-*minmaxY.first);
			dOffs = *minmaxY.first;
		}

		dX0Err = dX0 * 0.1;
		dSigErr = dSig * 0.1;
		dAmpErr = dAmp * 0.1;
		dOffsErr = dOffs * 0.1;

		bAmpFixed = bSigFixed = bX0Fixed = bOffsFixed = 0;
	}

	std::vector<std::string> vecParamNames = { "x0", "sig", "amp", "offs" };
	std::vector<t_real> vecVals = { dX0, dSig, dAmp, dOffs };
	std::vector<t_real> vecErrs = { dX0Err, dSigErr, dAmpErr, dOffsErr };
	std::vector<bool> vecFixed = { bX0Fixed, bSigFixed, bAmpFixed, bOffsFixed };

	if(bUseSlope)
	{
		vecParamNames.push_back("slope");
		vecVals.push_back(dSlope);
		vecErrs.push_back(dSlopeErr);
		vecFixed.push_back(bSlopeFixed);
	}

	bool bOk = false;
	if(bUseSlope)
		bOk = Fit<iFuncArgsSloped>(funcSloped, vecParamNames, vecVals, vecErrs, vecFixed);
	else
		bOk = Fit<iFuncArgs>(func, vecParamNames, vecVals, vecErrs, vecFixed);
	if(!bOk)
		return;

	for(t_real &d : vecErrs) d = std::abs(d);
	vecVals[1] = std::abs(vecVals[1]);

	m_pFitParamDlg->SetX0(vecVals[0]);		m_pFitParamDlg->SetX0Err(vecErrs[0]);
	m_pFitParamDlg->SetSig(vecVals[1]);		m_pFitParamDlg->SetSigErr(vecErrs[1]);
	m_pFitParamDlg->SetAmp(vecVals[2]);		m_pFitParamDlg->SetAmpErr(vecErrs[2]);
	m_pFitParamDlg->SetOffs(vecVals[3]);	m_pFitParamDlg->SetOffsErr(vecErrs[3]);

	if(bUseSlope)
	{
		m_pFitParamDlg->SetSlope(vecVals[4]);
		m_pFitParamDlg->SetSlopeErr(vecErrs[4]);
	}
}


/**
 * fits a lorentzian peak
 */
void ScanViewerDlg::FitLorentz()
{
	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return;


	// prefit
	unsigned int iOrder = m_settings.value("spline_order", 6).toInt();
	std::vector<t_real> vecMaximaX, vecMaximaSize, vecMaximaWidth;
	tl::find_peaks<t_real>(m_vecX.size(), m_vecX.data(), m_vecY.data(), iOrder,
		vecMaximaX, vecMaximaSize, vecMaximaWidth, g_dEpsGfx);


	const bool bUseSlope = checkSloped->isChecked();

	auto func = tl::lorentz_model_amp<t_real>;
	auto funcSloped = tl::lorentz_model_amp_slope<t_real>;
	constexpr std::size_t iFuncArgs = 5;
	constexpr std::size_t iFuncArgsSloped = iFuncArgs+1;

	t_real_glob dAmp = m_pFitParamDlg->GetAmp(),	dAmpErr = m_pFitParamDlg->GetAmpErr();
	t_real_glob dHWHM = m_pFitParamDlg->GetHWHM(),	dHWHMErr = m_pFitParamDlg->GetHWHMErr();
	t_real_glob dX0 = m_pFitParamDlg->GetX0(),		dX0Err = m_pFitParamDlg->GetX0Err();
	t_real_glob dOffs = m_pFitParamDlg->GetOffs(),	dOffsErr = m_pFitParamDlg->GetOffsErr();
	t_real_glob dSlope = m_pFitParamDlg->GetSlope(),	dSlopeErr = m_pFitParamDlg->GetSlopeErr();

	bool bAmpFixed = m_pFitParamDlg->GetAmpFixed();
	bool bHWHMFixed = m_pFitParamDlg->GetHWHMFixed();
	bool bX0Fixed = m_pFitParamDlg->GetX0Fixed();
	bool bOffsFixed = m_pFitParamDlg->GetOffsFixed();
	bool bSlopeFixed = m_pFitParamDlg->GetSlopeFixed();

	// automatic parameter determination
	if(!m_pFitParamDlg->WantParams())
	{
		auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
		auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());

		// use prefitter values
		if(vecMaximaX.size())
		{
			dX0 = vecMaximaX[0];
			dAmp = vecMaximaSize[0];
			dHWHM = 0.5*vecMaximaWidth[0] * 0.5;
			dOffs = *minmaxY.first;
		}
		// try to guess values
		else
		{
			dX0 = m_vecX[minmaxY.second - m_vecY.begin()];
			dHWHM = std::abs((*minmaxX.second-*minmaxX.first) * 0.5);
			dAmp = std::abs(*minmaxY.second-*minmaxY.first);
			dOffs = *minmaxY.first;
		}

		dX0Err = dX0 * 0.1;
		dHWHMErr = dHWHM * 0.1;
		dAmpErr = dAmp * 0.1;
		dOffsErr = dOffs * 0.1;

		bAmpFixed = bHWHMFixed = bX0Fixed = bOffsFixed = 0;
	}

	std::vector<std::string> vecParamNames = { "x0", "hwhm", "amp", "offs" };
	std::vector<t_real> vecVals = { dX0, dHWHM, dAmp, dOffs };
	std::vector<t_real> vecErrs = { dX0Err, dHWHMErr, dAmpErr, dOffsErr  };
	std::vector<bool> vecFixed = { bX0Fixed, bHWHMFixed, bAmpFixed, bOffsFixed };

	if(bUseSlope)
	{
		vecParamNames.push_back("slope");
		vecVals.push_back(dSlope);
		vecErrs.push_back(dSlopeErr);
		vecFixed.push_back(bSlopeFixed);
	}

	bool bOk = false;
	if(bUseSlope)
		bOk = Fit<iFuncArgsSloped>(funcSloped, vecParamNames, vecVals, vecErrs, vecFixed);
	else
		bOk = Fit<iFuncArgs>(func, vecParamNames, vecVals, vecErrs, vecFixed);
	if(!bOk)
		return;

	for(t_real &d : vecErrs) d = std::abs(d);
	vecVals[1] = std::abs(vecVals[1]);

	m_pFitParamDlg->SetX0(vecVals[0]);		m_pFitParamDlg->SetX0Err(vecErrs[0]);
	m_pFitParamDlg->SetHWHM(vecVals[1]);	m_pFitParamDlg->SetHWHMErr(vecErrs[1]);
	m_pFitParamDlg->SetAmp(vecVals[2]);		m_pFitParamDlg->SetAmpErr(vecErrs[2]);
	m_pFitParamDlg->SetOffs(vecVals[3]);	m_pFitParamDlg->SetOffsErr(vecErrs[3]);

	if(bUseSlope)
	{
		m_pFitParamDlg->SetSlope(vecVals[4]);
		m_pFitParamDlg->SetSlopeErr(vecErrs[4]);
	}
}


#ifndef HAS_COMPLEX_ERF

void ScanViewerDlg::FitVoigt() {}

#else


/**
 * fits a voigtian peak
 */
void ScanViewerDlg::FitVoigt()
{
	if(std::min(m_vecX.size(), m_vecY.size()) == 0)
		return;


	// prefit
	unsigned int iOrder = m_settings.value("spline_order", 6).toInt();
	std::vector<t_real> vecMaximaX, vecMaximaSize, vecMaximaWidth;
	tl::find_peaks<t_real>(m_vecX.size(), m_vecX.data(), m_vecY.data(), iOrder,
		vecMaximaX, vecMaximaSize, vecMaximaWidth, g_dEpsGfx);


	const bool bUseSlope = checkSloped->isChecked();

	auto func = tl::voigt_model_amp<t_real>;
	auto funcSloped = tl::voigt_model_amp_slope<t_real>;
	constexpr std::size_t iFuncArgs = 6;
	constexpr std::size_t iFuncArgsSloped = iFuncArgs+1;

	t_real_glob dAmp = m_pFitParamDlg->GetAmp(),	dAmpErr = m_pFitParamDlg->GetAmpErr();
	t_real_glob dSig = m_pFitParamDlg->GetSig(),	dSigErr = m_pFitParamDlg->GetSigErr();
	t_real_glob dHWHM = m_pFitParamDlg->GetHWHM(),	dHWHMErr = m_pFitParamDlg->GetHWHMErr();
	t_real_glob dX0 = m_pFitParamDlg->GetX0(),		dX0Err = m_pFitParamDlg->GetX0Err();
	t_real_glob dOffs = m_pFitParamDlg->GetOffs(),	dOffsErr = m_pFitParamDlg->GetOffsErr();
	t_real_glob dSlope = m_pFitParamDlg->GetSlope(),	dSlopeErr = m_pFitParamDlg->GetSlopeErr();

	bool bAmpFixed = m_pFitParamDlg->GetAmpFixed();
	bool bSigFixed = m_pFitParamDlg->GetSigFixed();
	bool bHWHMFixed = m_pFitParamDlg->GetHWHMFixed();
	bool bX0Fixed = m_pFitParamDlg->GetX0Fixed();
	bool bOffsFixed = m_pFitParamDlg->GetOffsFixed();
	bool bSlopeFixed = m_pFitParamDlg->GetSlopeFixed();

	// automatic parameter determination
	if(!m_pFitParamDlg->WantParams())
	{
		auto minmaxX = std::minmax_element(m_vecX.begin(), m_vecX.end());
		auto minmaxY = std::minmax_element(m_vecY.begin(), m_vecY.end());

		// use prefitter values
		if(vecMaximaX.size())
		{
			dX0 = vecMaximaX[0];
			dAmp = vecMaximaSize[0];
			dSig = tl::get_FWHM2SIGMA<t_real>()*vecMaximaWidth[0] * 0.5*0.5;
			dHWHM = 0.5*vecMaximaWidth[0] * 0.5*0.5;
			dOffs = *minmaxY.first;
		}
		// try to guess values
		else
		{
			dX0 = m_vecX[minmaxY.second - m_vecY.begin()];
			dHWHM = std::abs((*minmaxX.second-*minmaxX.first) * 0.25);
			dSig = std::abs((*minmaxX.second-*minmaxX.first) * 0.25);
			dAmp = std::abs(*minmaxY.second-*minmaxY.first);
			dOffs = *minmaxY.first;
		}

		dX0Err = dX0 * 0.1;
		dHWHMErr = dHWHM * 0.1;
		dSigErr = dSig * 0.1;
		dAmpErr = dAmp * 0.1;
		dOffsErr = dOffs * 0.1;

		bAmpFixed = bHWHMFixed = bSigFixed = bX0Fixed = bOffsFixed = 0;
	}

	std::vector<std::string> vecParamNames = { "x0", "sig", "hwhm", "amp", "offs" };
	std::vector<t_real> vecVals = { dX0, dSig, dHWHM, dAmp, dOffs };
	std::vector<t_real> vecErrs = { dX0Err, dSigErr, dHWHMErr, dAmpErr, dOffsErr };
	std::vector<bool> vecFixed = { bX0Fixed, bSigFixed, bHWHMFixed, bAmpFixed, bOffsFixed };

	if(bUseSlope)
	{
		vecParamNames.push_back("slope");
		vecVals.push_back(dSlope);
		vecErrs.push_back(dSlopeErr);
		vecFixed.push_back(bSlopeFixed);
	}

	bool bOk = false;
	if(bUseSlope)
		bOk = Fit<iFuncArgsSloped>(funcSloped, vecParamNames, vecVals, vecErrs, vecFixed);
	else
		bOk = Fit<iFuncArgs>(func, vecParamNames, vecVals, vecErrs, vecFixed);
	if(!bOk)
		return;

	for(t_real &d : vecErrs) d = std::abs(d);
	vecVals[1] = std::abs(vecVals[1]);
	vecVals[2] = std::abs(vecVals[2]);

	m_pFitParamDlg->SetX0(vecVals[0]);		m_pFitParamDlg->SetX0Err(vecErrs[0]);
	m_pFitParamDlg->SetSig(vecVals[1]);		m_pFitParamDlg->SetSigErr(vecErrs[1]);
	m_pFitParamDlg->SetHWHM(vecVals[2]);	m_pFitParamDlg->SetHWHMErr(vecErrs[2]);
	m_pFitParamDlg->SetAmp(vecVals[3]);		m_pFitParamDlg->SetAmpErr(vecErrs[3]);
	m_pFitParamDlg->SetOffs(vecVals[4]);	m_pFitParamDlg->SetOffsErr(vecErrs[4]);

	if(bUseSlope)
	{
		m_pFitParamDlg->SetSlope(vecVals[5]);
		m_pFitParamDlg->SetSlopeErr(vecErrs[5]);
	}
}
#endif  // HAS_COMPLEX_ERF
