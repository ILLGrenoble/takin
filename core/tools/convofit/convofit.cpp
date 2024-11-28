/**
 * convolution fitting
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

#include "tlibs/file/prop.h"
#include "tlibs/file/loaddat.h"
#include "tlibs/log/log.h"
#include "tlibs/log/debug.h"
#include "tlibs/math/rand.h"

#include <iostream>
#include <fstream>
#include <locale>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "convofit.h"
#include "convofit_import.h"
#include "scan.h"
#include "model.h"
#include "../monteconvo/monteconvo_common.h"
#include "../monteconvo/sqwfactory.h"
#include "../res/defs.h"

using t_real = t_real_reso;


// ----------------------------------------------------------------------------
// default plotter

#define DEFAULT_TERM "x11 noraise"
//#define DEFAULT_TERM "qt noraise"


/**
 * plotter init callback function
 */
static void* init_convofit_plot(const std::string& strTerm)
{
	tl::GnuPlot<t_real> *pPlt = new tl::GnuPlot<t_real>();
	if(!pPlt->Init())
		return nullptr;

	pPlt->SetTerminal(0, strTerm.c_str());
	return pPlt;
}


/**
 * plotter deinit callback function
 */
static void deinit_convofit_plot(void *&_pPlt)
{
	if(_pPlt)
	{
		tl::GnuPlot<t_real> *pPlt = reinterpret_cast<tl::GnuPlot<t_real>*>(_pPlt);
		delete pPlt;
		_pPlt = nullptr;
	}
}


/**
 * plotter callback function
 */
static void convofit_plot(void* _pPlt, const char* pcX, const char* pcY, const char *pcTitle,
	const tl::PlotObj<t_real>& pltMeas, const tl::PlotObj<t_real>& pltMod, bool bIsFinal)
{
	if(!_pPlt)
		return;
	tl::GnuPlot<t_real> *pPlt = reinterpret_cast<tl::GnuPlot<t_real>*>(_pPlt);
	pPlt->StartPlot();

	if(pcTitle)
		pPlt->SetTitle(pcTitle);
	if(pcX)
		pPlt->SetXLabel(pcX);
	if(pcY)
		pPlt->SetYLabel(pcY);

	pPlt->AddLine(pltMod);
	pPlt->AddLine(pltMeas);

	pPlt->FinishPlot();
}

// ----------------------------------------------------------------------------


Convofit::Convofit(bool bUseDefaultPlotter)
{
	if(bUseDefaultPlotter)
	{
		addsig_initplotter(&::init_convofit_plot);
		addsig_deinitplotter(&::deinit_convofit_plot);
		addsig_plot(&::convofit_plot);
	}
}


Convofit::~Convofit()
{
	m_sigDeinitPlotter(m_pPlt);
}


// ----------------------------------------------------------------------------
// global command line overrides
bool g_bVerbose = false;
bool g_bSkipFit = false;
bool g_bUseValuesFromModel = false;
unsigned int g_iNumNeutrons = 0;
std::string g_strSetParams;
std::string g_strOutFileSuffix;
unsigned int g_iPlotPoints = 0;
unsigned int g_iPlotSkipBegin = 0;
unsigned int g_iPlotSkipEnd = 0;
// ----------------------------------------------------------------------------


bool Convofit::run_job(const std::string& _strJob)
{
	// --------------------------------------------------------------------
	// set working directory for job
	bool bChangedCWD = 0;
	fs::path pathProg = fs::system_complete(_strJob).remove_filename();
	if(pathProg.filename().string() == ".")		// remove "./"
		pathProg.remove_filename();
	fs::path pathCWD = fs::system_complete(fs::current_path());
	if(pathProg != pathCWD)
	{
		fs::current_path(pathProg);
		bChangedCWD = 1;
		tl::log_debug("Working directory: ", pathProg.string(), ".");
	}
	// --------------------------------------------------------------------


	std::string strJob = _strJob;

	if(bChangedCWD)
		strJob = fs::path(strJob).filename().string();


	// if a monteconvo file is given, convert it to a convofit job file
	tl::Prop<std::string> propMC;
	if(propMC.Load(strJob.c_str(), tl::PropType::XML) &&
		propMC.Exists("taz/monteconvo"))
	{
		tl::log_info("Importing monteconvo file \"", strJob, "\".");
		strJob = convert_monteconvo(propMC);
		if(strJob == "")
			return false;
		tl::log_info("Converted convofit file is \"", strJob, "\".");
	}


	const unsigned iSeed = tl::get_rand_seed();
	tl::init_rand_seed(iSeed);

	// Parameters
	tl::Prop<std::string> prop;
	if(!prop.Load(strJob.c_str(), tl::PropType::INFO))
	{
		tl::log_err("Cannot load job file \"", strJob, "\".");
		return false;
	}

	std::string strScFile = prop.Query<std::string>("input/scan_file");
	if(strScFile == "")	// "scan_file_0" is synonymous to "scan_file"
		strScFile = prop.Query<std::string>("input/scan_file_0");

	boost::optional<unsigned> iMainScanAxis = prop.QueryOpt<unsigned>("input/scan_axis");  // 1-based scan axis
	if(!iMainScanAxis)	// "scan_axis_0" is synonymous to "scan_axis"
		iMainScanAxis = prop.Query<unsigned>("input/scan_axis_0", 0);
	if(!iMainScanAxis)
		iMainScanAxis = 0;	// automatic selection

	std::string strTempCol = prop.Query<std::string>("input/temp_col");
	std::string strFieldCol = prop.Query<std::string>("input/field_col");
	bool bTempOverride = prop.Exists("input/temp_override");
	bool bFieldOverride = prop.Exists("input/field_override");
	t_real dTempOverride = prop.QueryAndParse<t_real>("input/temp_override");
	t_real dFieldOverride = prop.QueryAndParse<t_real>("input/field_override");
	std::string strCntCol = prop.Query<std::string>("input/counts_col");
	std::string strMonCol = prop.Query<std::string>("input/monitor_col");
	std::string strCntErrCol = prop.Query<std::string>("input/counts_err_col");
	std::string strMonErrCol = prop.Query<std::string>("input/monitor_err_col");
	std::string strResFile = prop.Query<std::string>("input/instrument_file");
	if(strResFile == "")	// "instrument_file_0" is synonymous to "instrument_file"
		strResFile = prop.Query<std::string>("input/instrument_file_0");

	std::string strSqwMod = prop.Query<std::string>("input/sqw_model");
	std::string strSqwFile = prop.Query<std::string>("input/sqw_file");
	std::string strTempVar = prop.Query<std::string>("input/sqw_temp_var", "T");
	std::string strFieldVar = prop.Query<std::string>("input/sqw_field_var", "");
	std::string strSetParams = prop.Query<std::string>("input/sqw_set_params", "");

	bool bNormToMon = prop.Query<bool>("input/norm_to_monitor", true);
	bool bFlipCoords = prop.Query<bool>("input/flip_lhs_rhs", false);
	bool bUseFirstAndLastScanPt = prop.Query<bool>("input/use_first_last_pt", false);
	bool bAllowScanMerging = prop.Query<bool>("input/allow_scan_merging", false);

	if(g_strSetParams != "")
	{
		if(strSetParams != "")
			strSetParams += "; ";
		strSetParams += g_strSetParams;
	}

	Filter filter;
	if(prop.Exists("input/filter_lower"))
		filter.dLower = prop.QueryAndParse<t_real>("input/filter_lower", 0);
	if(prop.Exists("input/filter_upper"))
		filter.dUpper = prop.QueryAndParse<t_real>("input/filter_upper", 0);


	// --------------------------------------------------------------------
	// files in inner vector will be merged
	// files in outer vector will be used for multi-function fitting
	std::vector<std::vector<std::string>> vecvecScFiles;
	std::vector<unsigned> vecScanAxes;

	// primary scan file(s)
	{
		std::vector<std::string> vecScFiles;
		tl::get_tokens<std::string, std::string>(strScFile, ";", vecScFiles);
		std::for_each(vecScFiles.begin(), vecScFiles.end(), [](std::string& str){ tl::trim(str); });
		vecvecScFiles.emplace_back(std::move(vecScFiles));

		vecScanAxes.push_back(*iMainScanAxis);
	}

	// get secondary scan files for multi-function fitting
	for(std::size_t iSecFile=1; 1; ++iSecFile)
	{
		// scan file names
		std::string strSecFile = "input/scan_file_" + tl::var_to_str(iSecFile);
		std::string strSecScFile = prop.Query<std::string>(strSecFile, "");
		if(strSecScFile == "")
			break;

		std::vector<std::string> vecSecScFiles;
		tl::get_tokens<std::string, std::string>(strSecScFile, ";", vecSecScFiles);
		std::for_each(vecSecScFiles.begin(), vecSecScFiles.end(), [](std::string& str){ tl::trim(str); });
		vecvecScFiles.emplace_back(std::move(vecSecScFiles));

		// scan axes
		unsigned iScanAxis = prop.Query<unsigned>("input/scan_axis_" + tl::var_to_str(iSecFile), *iMainScanAxis);
		vecScanAxes.push_back(iScanAxis);
	}

	// fall-back default
	if(!vecScanAxes.size())
		vecScanAxes.push_back(*iMainScanAxis);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// primary resolution file
	std::vector<std::string> vecResFiles({strResFile});

	// get secondary resolution files for multi-function fitting
	for(std::size_t iSecFile=1; 1; ++iSecFile)
	{
		std::string strSecFile = "input/instrument_file_" + tl::var_to_str(iSecFile);
		std::string strSecResFile = prop.Query<std::string>(strSecFile, "");
		if(strSecResFile == "")
			break;
		tl::trim(strSecResFile);
		vecResFiles.emplace_back(std::move(strSecResFile));
	}

	if(vecResFiles.size()!=1 && vecResFiles.size()!=vecvecScFiles.size())
	{
		tl::log_err("Number of resolution files has to be either one or match the number of scan file groups.");
		tl::log_err("Number of scan file groups: ", vecvecScFiles.size(), ", number of resolution files: ", vecResFiles.size(), ".");
		return false;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// optional s(Q, E) parameter overrides
	std::vector<std::string> vecSetParams;

	// get secondary resolution files for multi-function fitting
	for(std::size_t iGroup = 0; iGroup < vecvecScFiles.size(); ++iGroup)
	{
		std::string strSecParams = "input/sqw_set_params_" + tl::var_to_str(iGroup);
		std::string strSetParamsOverrides = prop.Query<std::string>(strSecParams, "");
		tl::trim(strSetParamsOverrides);
		vecSetParams.emplace_back(std::move(strSetParamsOverrides));
	}

	if(vecSetParams.size() != vecvecScFiles.size())
	{
		tl::log_err("Number of S(Q, E) parameter overrides has to match the number of scan file groups.");
		tl::log_err("Number of scan file groups: ", vecvecScFiles.size(), ", number of parameter files: ", vecSetParams.size(), ".");
		return false;
	}
	// --------------------------------------------------------------------

	t_real eps_plane_dist = prop.Query<t_real>("tolerances/plane_dist", EPS_PLANE);

	unsigned iNumNeutrons = prop.Query<unsigned>("montecarlo/neutrons", 1000);
	unsigned iNumSample = prop.Query<unsigned>("montecarlo/sample_positions", 1);
	int iRecycleMC = prop.Query<int>("montecarlo/recycle_neutrons", 1);

	// global override
	if(g_iNumNeutrons > 0)
		iNumNeutrons = g_iNumNeutrons;

	std::string strResAlgo = prop.Query<std::string>("resolution/algorithm", "pop");

	// -1: unchanged (use curvature value from reso file), 0: flat, 1: optimal
	int iResFocMonoV = prop.Query<int>("resolution/focus_mono_v", -1);
	int iResFocMonoH = prop.Query<int>("resolution/focus_mono_h", -1);
	int iResFocAnaV = prop.Query<int>("resolution/focus_ana_v", -1);
	int iResFocAnaH = prop.Query<int>("resolution/focus_ana_h", -1);

	std::string strMinimiser = prop.Query<std::string>("fitter/minimiser");
	int iStrat = prop.Query<int>("fitter/strategy", 0);
	t_real dSigma = prop.Query<t_real>("fitter/sigma", 1.);
	unsigned iNumThreads = prop.Query<unsigned>("fitter/num_threads", g_iMaxThreads);

	bool bDoFit = prop.Query<bool>("fitter/do_fit", true);
	if(g_bSkipFit) bDoFit = 0;

	unsigned int iMaxFuncCalls = prop.Query<unsigned>("fitter/max_funccalls", 0);
	t_real dTolerance = prop.Query<t_real>("fitter/tolerance", 0.5);

	std::string strScOutFile = prop.Query<std::string>("output/scan_file");
	std::string strModOutFile = prop.Query<std::string>("output/model_file");
	std::string strLogOutFile = prop.Query<std::string>("output/log_file");
	bool bPlot = prop.Query<bool>("output/plot", false);
	bool bPlotIntermediate = prop.Query<bool>("output/plot_intermediate", false);

	unsigned int iPlotPoints = prop.Query<unsigned>("output/plot_points", 128);
	unsigned int iPlotPointsSkipBegin = prop.Query<unsigned>("output/plot_points_skip_begin", 0);
	unsigned int iPlotPointsSkipEnd = prop.Query<unsigned>("output/plot_points_skip_end", 0);
	if(g_iPlotPoints) iPlotPoints = g_iPlotPoints;
	if(g_iPlotSkipBegin) iPlotPointsSkipBegin = g_iPlotSkipBegin;
	if(g_iPlotSkipEnd) iPlotPointsSkipEnd = g_iPlotSkipEnd;

	if(bPlot || bPlotIntermediate)
	{
		std::string strTerm = prop.Query<std::string>("output/plot_term", DEFAULT_TERM);
		boost::optional<void*> optPlt = m_sigInitPlotter(strTerm);
		if(optPlt)
			m_pPlt = *optPlt;

		if(!m_pPlt)
			tl::log_err("Could not initialise plotter. Is gnuplot (correctly) installed?");
	}


	if(g_strOutFileSuffix != "")
		strLogOutFile += g_strOutFileSuffix;

	// thread-local debug log
	std::unique_ptr<std::ostream> ofstrLog;
	if(strLogOutFile != "")
	{
		ofstrLog.reset(new std::ofstream(strLogOutFile));

		for(tl::Log* plog : { &tl::log_info, &tl::log_warn, &tl::log_err, &tl::log_crit, &tl::log_debug })
			plog->AddOstr(ofstrLog.get(), 0, 1);
	}

	if(strScOutFile=="" || strModOutFile=="")
	{
		tl::log_err("Not output files selected.");
		return false;
	}


	std::string strFitParams = prop.Query<std::string>("fit_parameters/params");
	std::string strFitValues = prop.Query<std::string>("fit_parameters/values");
	std::string strFitErrors = prop.Query<std::string>("fit_parameters/errors");
	std::string strFitFixed = prop.Query<std::string>("fit_parameters/fixed");

	// parameter using either strModInFile if it is defined or strModOutFile if
	// reuse_values_from_model_file is set to 1
	bool bUseValuesFromModel = prop.Query<bool>("fit_parameters/reuse_values_from_model_file", false);
	std::string strModInFile = prop.Query<std::string>("input/model_file");
	if(g_bUseValuesFromModel || strModInFile != "")
		bUseValuesFromModel = 1;


	std::vector<std::string> vecFitParams;
	tl::get_tokens<std::string, std::string>(strFitParams, " \t\n,;", vecFitParams);
	std::vector<t_real> vecFitValues;
	tl::parse_tokens<t_real, std::string>(strFitValues, " \t\n,;", vecFitValues);
	std::vector<t_real> vecFitErrors;
	tl::parse_tokens<t_real, std::string>(strFitErrors, " \t\n,;", vecFitErrors);
	std::vector<bool> vecFitFixed;
	tl::get_tokens<bool, std::string>(strFitFixed, " \t\n,;", vecFitFixed);

	if(vecFitParams.size() != vecFitValues.size() ||
		vecFitParams.size() != vecFitErrors.size() ||
		vecFitParams.size() != vecFitFixed.size())
	{
		tl::log_err("Fit parameter size mismatch.");
		return false;
	}


	if(bUseValuesFromModel)
	{
		const std::string *pModOverrideFile = &strModOutFile;
		if(strModInFile != "")
		{
			pModOverrideFile = &strModInFile;
			tl::log_info("Overriding parameters with model input file \"", strModInFile, "\".");
		}
		else
		{
			tl::log_info("Overriding parameters with model output file \"", strModOutFile, "\".");
		}

		tl::DatFile<t_real, char> datMod;
		if(datMod.Load(*pModOverrideFile))
		{
			const auto& mapHdr = datMod.GetHeader();
			//for(const auto& pair : mapHdr)
			//	std::cout << pair.first << ", " << pair.second << std::endl;

			for(std::size_t iParam=0; iParam<vecFitParams.size(); ++iParam)
			{
				auto iterParam = mapHdr.find(vecFitParams[iParam]);
				if(iterParam != mapHdr.end())
				{
					std::pair<std::string, std::string> pairVal =
						tl::split_first<std::string>(iterParam->second, "+-", 1, 1);

					const std::string& strNewVal = pairVal.first;
					const std::string& strNewErr = pairVal.second;

					vecFitValues[iParam] = tl::str_to_var<t_real>(strNewVal);
					vecFitErrors[iParam] = tl::str_to_var<t_real>(strNewErr);

					tl::log_info("Overriding parameter \"", iterParam->first,
						"\" with model value: ", strNewVal, " +- ", strNewErr, ".");
				}
				else
				{
					tl::log_warn("Requested override parameter \"",
						vecFitParams[iParam], "\" is not available in model file.");
				}
			}
		}
		else
		{
			tl::log_err("Parameter override using model file requested, but model file \"",
				*pModOverrideFile, "\" is invalid.");
		}
	}

	if(g_strOutFileSuffix != "")
	{
		strScOutFile += g_strOutFileSuffix;
		strModOutFile += g_strOutFileSuffix;
	}



	// --------------------------------------------------------------------
	// Scan files
	std::vector<Scan> vecSc;
	for(std::size_t iSc = 0; iSc < vecvecScFiles.size(); ++iSc)
	{
		Scan sc;
		if(strTempCol != "")
			sc.strTempCol = strTempCol;
		if(strFieldCol != "")
			sc.strFieldCol = strFieldCol;
		sc.strCntCol = strCntCol;
		sc.strMonCol = strMonCol;
		sc.strCntErrCol = strCntErrCol;
		sc.strMonErrCol = strMonErrCol;

		if(vecvecScFiles.size() > 1)
			tl::log_info("Loading scan group ", iSc, ".");
		if(!load_file(vecvecScFiles[iSc], sc, bNormToMon, filter,
			bFlipCoords, bAllowScanMerging, bUseFirstAndLastScanPt,
			vecScanAxes[iSc], g_bVerbose))
		{
			tl::log_err("Cannot load scan files of group ", iSc, ".");
			continue;
		}

		// read back the determined scan axis
		vecScanAxes[iSc] = sc.m_iScIdx + 1;

		vecSc.emplace_back(std::move(sc));
	}
	if(!vecSc.size())
	{
		tl::log_err("No scans loaded.");
		return false;
	}

	tl::log_info("Number of scan groups: ", vecSc.size(), ".");

	// scan plot objects
	std::vector<tl::PlotObj<t_real>> pltMeas;
	if(bPlot || bPlotIntermediate)
	{
		pltMeas.reserve(vecSc.size());

		// data points for plotter
		for(std::size_t scan_group=0; scan_group<vecSc.size(); ++scan_group)
		{
			tl::PlotObj<t_real> thescan;

			thescan.vecX = vecSc[scan_group].vecX;
			thescan.vecY = vecSc[scan_group].vecCts;
			thescan.vecErrY = vecSc[scan_group].vecCtsErr;
			thescan.linestyle = tl::STYLE_POINTS;
			thescan.oiColor = 0xff0000;

			pltMeas.emplace_back(std::move(thescan));
		}
	}
	// --------------------------------------------------------------------




	// --------------------------------------------------------------------
	// resolution files
	std::vector<TASReso> vecResos;
	for(std::size_t iGroup = 0; iGroup < vecResFiles.size(); ++iGroup)
	{
		const std::string& strCurResFile = vecResFiles[iGroup];

		TASReso reso;
		reso.SetPlaneDistTolerance(eps_plane_dist);

		tl::log_info("Loading instrument resolution file \"", strCurResFile, "\" for scan group ", iGroup, ".");
		if(!reso.LoadRes(strCurResFile.c_str()))
			return false;

		if(strResAlgo == "pop")
			reso.SetAlgo(ResoAlgo::POP);
		else if(strResAlgo == "pop_cn")
			reso.SetAlgo(ResoAlgo::POP_CN);
		else if(strResAlgo == "cn")
			reso.SetAlgo(ResoAlgo::CN);
		else if(strResAlgo == "eck")
			reso.SetAlgo(ResoAlgo::ECK);
		else if(strResAlgo == "vio" || strResAlgo == "viol")
			reso.SetAlgo(ResoAlgo::VIO);
		else
		{
			tl::log_err("Invalid resolution algorithm selected: \"", strResAlgo, "\".");
			return false;
		}

		{
			// TODO: distinguish between flat_v and flat_h
			unsigned ifocMode = unsigned(ResoFocus::FOC_UNCHANGED);

			if(iResFocMonoH == 0 && iResFocMonoV == 0) ifocMode |= unsigned(ResoFocus::FOC_MONO_FLAT);		// flat
			if(iResFocMonoH == 1) ifocMode |= unsigned(ResoFocus::FOC_MONO_H);		// horizontal
			if(iResFocMonoV == 1) ifocMode |= unsigned(ResoFocus::FOC_MONO_V);		// vertical

			if(iResFocAnaH == 0 && iResFocAnaV == 0) ifocMode |= unsigned(ResoFocus::FOC_ANA_FLAT);			// flat
			if(iResFocAnaH == 1) ifocMode |= unsigned(ResoFocus::FOC_ANA_H);	// horizontal
			if(iResFocAnaV == 1) ifocMode |= unsigned(ResoFocus::FOC_ANA_V);	// vertical

			reso.SetOptimalFocus(ResoFocus(ifocMode));
		}

		reso.SetRandomSamplePos(iNumSample);
		vecResos.emplace_back(std::move(reso));
	}

	// base parameter set for single-fits
	set_tasreso_params_from_scan(vecResos[0], vecSc[0]);
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// Model file
	tl::log_info("Loading S(Q, E) file \"", strSqwFile, "\".");
	std::shared_ptr<SqwBase> pSqw = construct_sqw(strSqwMod, strSqwFile);

	if(!pSqw)
	{
		tl::log_err("Invalid S(Q, E) model selected: \"", strSqwMod, "\".");
		return false;
	}

	if(!pSqw->IsOk())
	{
		tl::log_err("S(Q, E) model cannot be initialised.");
		return false;
	}
	SqwFuncModel mod(pSqw, vecResos);
	mod.SetSqwParamOverrides(vecSetParams);


	// temporary data buffer for plotting
	std::vector<t_real> vecModPlotX, vecModPlotY;
	std::mutex mtxPlot;

	// callback for outputting results
	mod.AddFuncResultSlot(
	[this, &pltMeas, &vecModPlotX, &vecModPlotY, bPlotIntermediate, &vecScanAxes, &mtxPlot]
		(t_real h, t_real k, t_real l, t_real E, t_real S, std::size_t scan_group)
	{
		if(g_bVerbose)
			tl::log_info("Q = (", h, ", ", k, ", ", l, ") rlu, E = ", E, " meV -> S = ", S);

		if(bPlotIntermediate)
		{
			t_real x = 0.;
			std::string scan_axis = "";
			unsigned iScanAxis = vecScanAxes[scan_group];

			switch(iScanAxis)
			{
				case 1: x = h; scan_axis = "h (rlu)"; break;
				case 2: x = k; scan_axis = "k (rlu)"; break;
				case 3: x = l; scan_axis = "l (rlu)"; break;
				case 4: x = E; scan_axis = "E (meV)"; break;
				default: x = E; scan_axis = "E (meV)"; break;
			}

			std::lock_guard<std::mutex> lock(mtxPlot);
			vecModPlotX.push_back(x);
			vecModPlotY.push_back(S);
			tl::sort_2(vecModPlotX.begin(), vecModPlotX.end(), vecModPlotY.begin());

			// convolution curve
			tl::PlotObj<t_real> pltMod;
			pltMod.vecX = vecModPlotX;
			pltMod.vecY = vecModPlotY;
			pltMod.linestyle = tl::STYLE_LINES_SOLID;
			pltMod.oiColor = 0x0000ff;
			pltMod.odSize = 1.5;

			if(scan_group < pltMeas.size())
			{
				std::string title = "Takin/Convofit, scan group #" + tl::var_to_str(scan_group);
				m_sigPlot(this->m_pPlt, scan_axis.c_str(), "Intensity", title.c_str(),
					pltMeas[scan_group], pltMod, false);
			}
		}
	});

	// callback for changed parameters
	mod.AddParamsChangedSlot(
	[&vecModPlotX, &vecModPlotY, bPlotIntermediate, iSeed, iRecycleMC](const std::string& strDescr)
	{
		tl::log_info("Changed model parameters: ", strDescr);

		if(bPlotIntermediate)
		{
			vecModPlotX.clear();
			vecModPlotY.clear();
		}

		// do we use the same MC neutrons again?
		if(iRecycleMC)
		{
			tl::init_rand_seed(iSeed);
			tl::log_debug("Resetting random seed to ", iSeed, ".");
		}
	});


	// only needed for multi-fits
	if(vecSc.size() > 1)
		mod.SetScans(&vecSc);

	mod.SetNumNeutrons(iNumNeutrons);
	// execution has to be in a determined order to recycle the same neutrons
	mod.SetUseThreads(iRecycleMC == 0);

	// if threads are used in the fitter or in the chi^2 function, we need to
	// even more aggressively recycle neutrons before the mc generation step
	if(iNumThreads && iRecycleMC == 2)
		mod.SetSeed(iSeed);

	if(bTempOverride)
	{
		for(Scan& sc : vecSc)
		{
			sc.dTemp = dTempOverride;
			sc.dTempErr = 0.;
		}
	}
	if(bFieldOverride)
	{
		for(Scan& sc : vecSc)
		{
			sc.dField = dFieldOverride;
			sc.dFieldErr = 0.;
		}
	}
	mod.SetOtherParamNames(strTempVar, strFieldVar);

	// base parameter set for single-fits
	set_model_params_from_scan(mod, vecSc[0]);

	// set the given individual global model parameters
	mod.GetSqwBase()->SetVars(strSetParams);

	tl::log_info("Number of neutrons: ", iNumNeutrons, ".");
	tl::log_info("Number of threads: ", iNumThreads, ".");
	tl::log_info("Model temperature variable: \"", strTempVar, "\", value: ", vecSc[0].dTemp);
	tl::log_info("Model field variable: \"", strFieldVar, "\", value: ", vecSc[0].dField);
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// fitting

	// non-S(Q, E) model parameters to include in fitting
	std::vector<std::string> vecNonSQEParms;

	for(std::size_t iParam = 0; iParam < vecFitParams.size(); ++iParam)
	{
		const std::string& strParam = vecFitParams[iParam];
		const t_real dValue = vecFitValues[iParam];

		// only include "scale", "slope", or "offs" if they are not fixed
		if(strParam == "scale")
		{
			mod.SetScale(dValue);
			if(!vecFitFixed[iParam])
				vecNonSQEParms.push_back("scale");
		}
		else if(strParam == "slope")
		{
			mod.SetSlope(dValue);
			if(!vecFitFixed[iParam])
				vecNonSQEParms.push_back("slope");
		}
		else if(strParam == "offs")
		{
			mod.SetOffs(dValue);
			if(!vecFitFixed[iParam])
				vecNonSQEParms.push_back("offs");
		}

		// ignore non-S(Q, E) model parameters
		if(strParam == "scale" || strParam == "slope" || strParam == "offs")
			continue;

		mod.AddModelFitParams(strParam, dValue, vecFitErrors[iParam]);
	}

	mod.SetNonSQEParams(vecNonSQEParms);


	tl::Chi2Function_mult<t_real_sc, std::vector> chi2fkt;
	// the vecSc[0] data sets are the default data set (will not be used if scan groups are defined)
	chi2fkt.AddFunc(&mod, vecSc[0].vecX.size(), vecSc[0].vecX.data(), vecSc[0].vecCts.data(), vecSc[0].vecCtsErr.data());
	chi2fkt.SetDebug(true);
	chi2fkt.SetSigma(dSigma);
	chi2fkt.SetNumThreads(iNumThreads);


	minuit::MnUserParameters params = mod.GetMinuitParams();
	for(std::size_t iParam = 0; iParam < vecFitParams.size(); ++iParam)
	{
		const std::string& strParam = vecFitParams[iParam];

		// ignore non-included parameters
		if(strParam == "scale" && std::find(vecNonSQEParms.begin(), vecNonSQEParms.end(), strParam) == vecNonSQEParms.end())
			continue;
		else if(strParam == "slope" && std::find(vecNonSQEParms.begin(), vecNonSQEParms.end(), strParam) == vecNonSQEParms.end())
			continue;
		else if(strParam == "offs" && std::find(vecNonSQEParms.begin(), vecNonSQEParms.end(), strParam) == vecNonSQEParms.end())
			continue;

		const t_real dVal = vecFitValues[iParam];
		const t_real dErr = vecFitErrors[iParam];

		params.SetValue(strParam, dVal);
		params.SetError(strParam, dErr);
		if(vecFitFixed[iParam])
			params.Fix(strParam);
	}
	// set initials
	mod.SetMinuitParams(params);


	minuit::MnStrategy strat(iStrat);

	std::unique_ptr<minuit::MnApplication> pmini;
	if(strMinimiser == "simplex")
		pmini.reset(new minuit::MnSimplex(chi2fkt, params, strat));
	else if(strMinimiser == "migrad")
		pmini.reset(new minuit::MnMigrad(chi2fkt, params, strat));
	else
	{
		tl::log_err("Invalid minimiser selected: \"", strMinimiser, "\".");
		return false;
	}

	bool bValidFit = 0;
	if(bDoFit)
	{
		tl::log_info("Performing fit.");
		minuit::FunctionMinimum mini = (*pmini)(iMaxFuncCalls, dTolerance);
		const minuit::MnUserParameterState& state = mini.UserState();
		bValidFit = mini.IsValid() && mini.HasValidParameters() && state.IsValid();
		mod.SetMinuitParams(state);

		std::ostringstream ostrMini;
		ostrMini << "Final fit results: " << mini << "\n";
		tl::log_info(ostrMini.str(), "Fit valid: ", bValidFit);
	}
	else
	{
		tl::log_info("Skipping fit, keeping initial values.");
	}


	tl::log_info("Saving results.");

	for(std::size_t iSc = 0; iSc < vecSc.size(); ++iSc)
	{
		const Scan& sc = vecSc[iSc];

		std::string strCurModOutFile = strModOutFile;
		std::string strCurScOutFile = strScOutFile;

		// append scan group number if this is a multi-fit
		if(vecSc.size() > 1)
		{
			strCurModOutFile += tl::var_to_str(iSc);
			strCurScOutFile += tl::var_to_str(iSc);
		}

		mod.SetParamSet(iSc);
		mod.Save(strCurModOutFile.c_str(), iPlotPoints, iPlotPointsSkipBegin, iPlotPointsSkipEnd);
		save_file(strCurScOutFile.c_str(), sc);
	}
	// --------------------------------------------------------------------



	// --------------------------------------------------------------------
	// final plotting of results
	if(bPlot && vecSc.size() <= 1)
	{
		tl::DatFile<t_real, char> datMod;
		if(datMod.Load(strModOutFile))
		{
			// convolution curve
			tl::PlotObj<t_real> pltMod;
			pltMod.vecX = datMod.GetColumn(0);
			pltMod.vecY = datMod.GetColumn(1);
			pltMod.linestyle = tl::STYLE_LINES_SOLID;
			pltMod.oiColor = 0x0000ff;
			pltMod.odSize = 1.5;

			if(pltMeas.size())
				m_sigPlot(m_pPlt, "", "Intensity", "Takin/Convofit result", pltMeas[0], pltMod, true);
		}
		else
		{
			tl::log_err("Cannot open model file \"", strModOutFile, "\" for plotting.");
		}
	}
	// --------------------------------------------------------------------


	// remove thread-local loggers
	if(!!ofstrLog)
	{
		for(tl::Log* plog : { &tl::log_info, &tl::log_warn, &tl::log_err, &tl::log_crit, &tl::log_debug })
			plog->RemoveOstr(ofstrLog.get());
	}


	if(!bDoFit)
		return true;
	return bValidFit;
}
