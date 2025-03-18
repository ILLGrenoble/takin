/**
 * Converts monteconvo format to convofit
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date mar-2017
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

#include "convofit_import.h"
#include "tlibs/file/tmp.h"
#include "tlibs/log/log.h"
#include "../res/defs.h"
#include "../monteconvo/monteconvo_common.h"

#include <map>
#include <boost/filesystem.hpp>

using t_real = t_real_reso;
namespace fs = boost::filesystem;


/**
 * if possible re-express strFile relative to pathBase
 */
static void trymake_relpath(std::string& strFile, const fs::path& pathBase)
{
	fs::path pathFile(strFile);
	pathFile.remove_filename();

	if(pathFile == pathBase)
	{
		fs::path pathOnlyFile = fs::path(strFile).filename();
		strFile = (fs::path(".") / pathOnlyFile).string();
	}
}


/**
 * converts from monteconvo input file format
 */
std::string convert_monteconvo(
	const tl::Prop<std::string>& propMC, std::string strFile, bool bRelPaths)
{
	// write temporary job file if none given
	if(strFile == "")
		strFile = tl::create_tmp_file<char>("convofit");


	// paths
	fs::path pathFile = fs::path(strFile);
	fs::path pathFileBase = pathFile.filename().replace_extension("");
	fs::path pathBase = pathFile;
	pathBase.remove_filename();

	if(bRelPaths)
		pathFile = ".";
	else
		pathFile = pathBase;


	// inputs
	std::map<std::string, std::string> mapJob;
	mapJob["input/scan_file"] = propMC.Query<std::string>("taz/monteconvo/scanfile");
	mapJob["input/instrument_file"] = propMC.Query<std::string>("taz/monteconvo/instr");
	mapJob["input/sqw_model"] = propMC.Query<std::string>("taz/monteconvo/sqw");
	mapJob["input/sqw_file"] = propMC.Query<std::string>("taz/monteconvo/sqw_conf");
	mapJob["input/counts_col"] = propMC.Query<std::string>("taz/convofit/counter");
	mapJob["input/monitor_col"] = propMC.Query<std::string>("taz/convofit/monitor");
	mapJob["input/flip_lhs_rhs"] = propMC.Query<std::string>("taz/convofit/flip_coords", "0");
	mapJob["input/allow_scan_merging"] = propMC.Query<std::string>("taz/monteconvo/allow_scan_merging", "0");
	mapJob["input/scan_axis"] = propMC.Query<std::string>("taz/convofit/scanaxis", "0");

	if(bRelPaths)
	{
		trymake_relpath(mapJob["input/scan_file"], pathBase);
		trymake_relpath(mapJob["input/instrument_file"], pathBase);
		trymake_relpath(mapJob["input/sqw_file"], pathBase);
	}


	// input overrides
	std::string strTOver = propMC.Query<std::string>("taz/convofit/temp_override");
	std::string strBOver = propMC.Query<std::string>("taz/convofit/field_override");
	std::string strSqwOver = propMC.Query<std::string>("taz/convofit/sqw_params");
	for(std::string *pstr : {&strTOver, &strBOver, &strSqwOver})
		tl::trim(*pstr);

	if(strTOver != "") mapJob["input/temp_override"] = strTOver;
	if(strBOver != "") mapJob["input/field_override"] = strBOver;


	fs::path pathScan = pathFile; pathScan /= "sc"; pathScan += pathFileBase; pathScan += ".dat";
	fs::path pathModel = pathFile; pathModel /= "mod"; pathModel += pathFileBase; pathModel += ".dat";
	fs::path pathLog = pathFile; pathLog /= "log"; pathLog += pathFileBase; pathLog += ".dat";

	mapJob["output/scan_file"] = pathScan.string();
	mapJob["output/model_file"] = pathModel.string();
	mapJob["output/log_file"] = pathLog.string();
	mapJob["output/plot"] = "1";
	mapJob["output/plot_intermediate"] = "1";


	// neutrons
	mapJob["montecarlo/neutrons"] = propMC.Query<std::string>("taz/monteconvo/neutron_count");
	mapJob["montecarlo/sample_positions"] =
		propMC.Query<std::string>("taz/monteconvo/sample_step_count", "1");
	mapJob["montecarlo/recycle_neutrons"] =
		propMC.Query<std::string>("taz/convofit/recycle_neutrons", "1");

	// fitting
	std::string strMin = "simplex";
	switch(propMC.Query<int>("taz/convofit/minimiser"))
	{
		case 0: strMin = "simplex"; break;
		case 1: strMin = "migrad"; break;
	}
	mapJob["fitter/minimiser"] = strMin;
	mapJob["fitter/strategy"] = propMC.Query<std::string>("taz/convofit/strategy", "1");
	mapJob["fitter/max_funccalls"] = propMC.Query<std::string>("taz/convofit/max_calls", "250");
	mapJob["fitter/tolerance"] = propMC.Query<std::string>("taz/convofit/tolerance", "25.");
	mapJob["fitter/sigma"] = "1.";


	// resolution
	std::string strAlgo = propMC.Query<std::string>("taz/monteconvo/algo");

	// if the algo id string is unrecognised, try the algo index
	if(strAlgo != "cn" && strAlgo != "pop_cn" && strAlgo != "pop"
		&& strAlgo != "eck" && strAlgo != "eck_ext" && strAlgo != "vio")
	{
		switch(propMC.Query<int>("taz/monteconvo/algo_idx"))
		{
			case 0: strAlgo = "cn"; break;
			case 1: strAlgo = "pop_cn"; break;
			case 2: strAlgo = "pop"; break;
			case 3: strAlgo = "eck"; break;
			case 4: strAlgo = "eck_ext"; break;
			case 5: strAlgo = "vio"; break;
		}
	}
	mapJob["resolution/algorithm"] = strAlgo;

	std::string strFocMonoH, strFocMonoV;
	std::string strFocAnaH, strFocAnaV;
	switch(propMC.Query<int>("taz/monteconvo/mono_foc"))
	{
		case 0: strFocMonoH = "-1"; strFocMonoV = "-1"; break;  // unchanged
		case 1: strFocMonoH = "0"; strFocMonoV = "0"; break;    // flat
		case 2: strFocMonoH = "1"; strFocMonoV = "0"; break;    // horizontal
		case 3: strFocMonoH = "0"; strFocMonoV = "1"; break;    // vertical
		case 4: strFocMonoH = "1"; strFocMonoV = "1"; break;    // both
	}
	switch(propMC.Query<int>("taz/monteconvo/ana_foc"))
	{
		case 0: strFocAnaH = "-1"; strFocAnaV = "-1"; break;    // unchanged
		case 1: strFocAnaH = "0"; strFocAnaV = "0"; break;      // flat
		case 2: strFocAnaH = "1"; strFocAnaV = "0"; break;      // horizontal
		case 3: strFocAnaH = "0"; strFocAnaV = "1"; break;      // vertical
		case 4: strFocAnaH = "1"; strFocAnaV = "1"; break;      // both
	}

	mapJob["resolution/focus_mono_v"] = strFocMonoV;
	mapJob["resolution/focus_mono_h"] = strFocMonoH;
	mapJob["resolution/focus_ana_v"] = strFocAnaV;
	mapJob["resolution/focus_ana_h"] = strFocAnaH;


	// epsilons
	t_real eps_plane = propMC.Query<t_real>("taz/monteconvo/eps_plane_dist", t_real(EPS_PLANE));
	mapJob["tolerances/plane_dist"] = tl::var_to_str(eps_plane);


	// parameters
	t_real dScale = propMC.Query<t_real>("taz/monteconvo/S_scale", t_real(1));
	t_real dSlope = propMC.Query<t_real>("taz/monteconvo/S_slope", t_real(0));
	t_real dOffs = propMC.Query<t_real>("taz/monteconvo/S_offs", t_real(0));

	mapJob["fit_parameters/params"] = "scale slope offs ";
	mapJob["fit_parameters/fixed"] = "1 1 1 ";
	mapJob["fit_parameters/values"] =
		tl::var_to_str(dScale) + " " + tl::var_to_str(dSlope) + " " + tl::var_to_str(dOffs) + " ";
	mapJob["fit_parameters/errors"] =
		tl::var_to_str(dScale*0.1) + " " + tl::var_to_str(dSlope*0.1) + " " + tl::var_to_str(dOffs*0.1) + " ";

	std::vector<std::string> vecParams = propMC.GetChildNodes("taz/monteconvo/sqw_params/");
	for(const std::string& strParam : vecParams)
	{
		t_real dVal = propMC.Query<t_real>(
			"taz/monteconvo/sqw_params/" + strParam, t_real(0));
		t_real dErr = propMC.Query<t_real>(
			"taz/monteconvo/sqw_errors/" + strParam, t_real(0));
		bool bFit = propMC.Query<bool>(
			"taz/monteconvo/sqw_fitvar/" + strParam, 0);

		if(bFit)	// open fit parameter
		{
			mapJob["fit_parameters/params"] += strParam + " ";
			mapJob["fit_parameters/values"] += tl::var_to_str(dVal) + " ";
			mapJob["fit_parameters/errors"] += tl::var_to_str(dErr) + " ";
			mapJob["fit_parameters/fixed"] += (bFit ? "0 " : "1 ");
		}
		else	// fixed sqw parameter
		{
			// get the full value string which might also contain code
			std::string strVal = propMC.Query<std::string>(
				"taz/monteconvo/sqw_params/" + strParam, "0");
			mapJob["input/sqw_set_params"] += strParam + " = " + strVal + "; ";
		}
	}

	// additional user-given sqw params
	if(strSqwOver != "")
		mapJob["input/sqw_set_params"] += strSqwOver + "; ";


	tl::Prop<std::string> job;
	job.Add(mapJob);
	if(!job.Save(strFile, tl::PropType::INFO))
	{
		tl::log_err("Could not convert to convofit file format.");
		return "";
	}

	return strFile;
}
