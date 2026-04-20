/**
 * subtract scan files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
 * @date 2014, apr-2026
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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
// g++ -I../.. -I. -DNO_IOSTR -o subscans ../../tlibs/file/loadinstr.cpp ../../tlibs/log/log.cpp subscans.cpp -std=c++17 -lboost_system -lboost_filesystem

#include <iostream>
#include <fstream>
#include <iomanip>
#include "tlibs/file/loadinstr.h"
#include "tlibs/log/log.h"
#include "tlibs/phys/neutrons.h"
#include "libs/globals.h"

using t_real = t_real_glob;


int main(int argc, char **argv)
{
	if(argc < 3)
	{
		tl::log_err("Usage:\n\t", argv[0], "<file_+> <file_-> <file_out>");
		return -1;
	}


	// original file
	tl::log_info("Loading file ", argv[1], "...");
	tl::FileInstrBase<t_real> *dat_add = tl::FileInstrBase<t_real>::LoadInstr(argv[1]);

	if(!dat_add)
	{
		tl::log_err("Cannot load data file ", argv[1], ".");
		return -1;
	}


	std::string strCnt = dat_add->GetCountVar();
	std::string strMon = dat_add->GetMonVar();
	std::string strTime = dat_add->GetTimerVar();
	tl::log_info("Count var: ", strCnt, ", monitor var: ", strMon, ", timer var: ", strTime);

	t_real Q_eps = 0.01;
	t_real E_eps = 0.05;
	std::array<t_real, 4> pos_add = dat_add->GetPosHKLE();

	tl::FileInstrBase<>::t_vecVals& vecCntAdd = dat_add->GetCol(strCnt);
	std::vector<t_real> vecCntAdd_err;;
	tl::FileInstrBase<>::t_vecVals& vecMonAdd = dat_add->GetCol(strMon);
	tl::FileInstrBase<>::t_vecVals& vecTimeAdd = dat_add->GetCol(strTime);


	// file to subtract
	const char* pcFile = argv[2];
	tl::log_info("Loading file ", pcFile, "...");
	tl::FileInstrBase<t_real> *dat_sub = tl::FileInstrBase<t_real>::LoadInstr(pcFile);
	if(!dat_sub)
	{
		tl::log_err("Cannot load data file ", pcFile, ".");
		return -1;
	}
	std::array<t_real, 4> pos_sub = dat_sub->GetPosHKLE();
	const tl::FileInstrBase<>::t_vecVals& vecCntSub = dat_sub->GetCol(strCnt);
	const tl::FileInstrBase<>::t_vecVals& vecMonSub = dat_sub->GetCol(strMon);
	const tl::FileInstrBase<>::t_vecVals& vecTimeSub = dat_sub->GetCol(strTime);

	if(vecCntSub.size() != vecCntAdd.size() || vecMonSub.size() != vecMonAdd.size() || vecTimeSub.size() != vecTimeAdd.size())
	{
		tl::log_err("Size mismatch in file ", pcFile, ".");
		return -1;
	}

	// TODO
	t_real T_add = 300.;
	t_real T_sub = 300.;

	for(unsigned int i = 0; i < vecCntAdd.size(); ++i)
	{
		if(!tl::float_equal<t_real>(std::get<0>(pos_add), std::get<0>(pos_sub), Q_eps))
			tl::log_warn("Mismatching h in index ", i, ".");
		if(!tl::float_equal<t_real>(std::get<1>(pos_add), std::get<1>(pos_sub), Q_eps))
			tl::log_warn("Mismatching k in index ", i, ".");
		if(!tl::float_equal<t_real>(std::get<2>(pos_add), std::get<2>(pos_sub), Q_eps))
			tl::log_warn("Mismatching l in index ", i, ".");
		if(!tl::float_equal<t_real>(std::get<3>(pos_add), std::get<3>(pos_sub), E_eps))
			tl::log_warn("Mismatching E in index ", i, ".");

		t_real bose_add = tl::bose<t_real>(std::get<3>(pos_add), T_add);
		t_real bose_sub = tl::bose<t_real>(std::get<3>(pos_sub), T_sub);

		t_real err_add = std::sqrt(vecCntAdd[i]);
		if(tl::float_equal<t_real>(err_add, 0.))
			err_add = 1.;

		t_real err_sub = std::sqrt(vecCntSub[i]);
		if(tl::float_equal<t_real>(err_sub, 0.))
			err_sub = 1.;

		vecCntAdd[i] -= vecCntSub[i] * vecMonAdd[i]/vecMonSub[i] * bose_add/bose_sub;

		err_sub *= vecMonAdd[i]/vecMonSub[i];
		vecCntAdd_err.push_back(std::sqrt(err_add*err_add + err_sub*err_sub));
	}
	delete dat_sub;


	tl::log_info("Saving file ", argv[argc - 1]);
	std::ofstream ofstr(argv[argc-1]);
	if(!ofstr.is_open())
	{
		tl::log_err("Cannot save data file ", argv[argc - 1], ".");
		return -1;
	}

	int w = 15;
	ofstr.precision(8);
	const tl::FileInstrBase<>::t_vecColNames& vecColNames = dat_add->GetColNames();
	for(const std::string& strColName : vecColNames)
		ofstr << std::setw(w) << std::right << strColName << " ";
	ofstr << std::setw(w) << std::right << "cnt_err" << " ";
	ofstr << "\n";

	for(unsigned int i = 0; i < vecCntAdd.size(); ++i)
	{
		for(const std::string& strColName : vecColNames)
		{
			tl::FileInstrBase<>::t_vecVals& vecCol = dat_add->GetCol(strColName);
			ofstr << std::setw(w) << std::right << vecCol[i] << " ";
		}
		ofstr << std::setw(w) << std::right << vecCntAdd_err[i] << " ";
		ofstr << "\n";
	}

	delete dat_add;
	return 0;
}
