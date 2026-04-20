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

// g++ -I../.. -I. -I/opt/homebrew/include -L/opt/homebrew/lib -DNO_IOSTR -o subscans ../../tlibs/file/loadinstr.cpp ../../tlibs/log/log.cpp subscans.cpp -std=c++17 -lboost_system -lboost_filesystem

#include <iostream>
#include <fstream>
#include <iomanip>

#include "tlibs/file/loadinstr.h"
#include "tlibs/log/log.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/math/stat.h"

#include "libs/globals.h"


using t_real = t_real_glob;


int main(int argc, char **argv)
{
	if(argc < 3)
	{
		tl::log_err("Usage:\n\t", argv[0], "<file_add> <file_sub> <file_out>");
		return -1;
	}


	// original file
	const char* pcFileAdd = argv[1];
	tl::log_info("Loading file ", pcFileAdd, "...");
	tl::FileInstrBase<t_real> *dat_add = tl::FileInstrBase<t_real>::LoadInstr(pcFileAdd);
	if(!dat_add)
	{
		tl::log_err("Cannot load data file ", pcFileAdd, ".");
		return -1;
	}

	std::string strCnt = dat_add->GetCountVar();
	std::string strMon = dat_add->GetMonVar();
	std::string strTemp = "TT";
	tl::log_info("Count var: ", strCnt, ", monitor var: ", strMon, ".");

	t_real Q_eps = 0.01;
	t_real E_eps = 0.05;
	std::array<t_real, 4> pos_add = dat_add->GetPosHKLE();

	tl::FileInstrBase<>::t_vecVals& vecCntAdd = dat_add->GetCol(strCnt);
	tl::FileInstrBase<>::t_vecVals& vecMonAdd = dat_add->GetCol(strMon);
	tl::FileInstrBase<>::t_vecVals& vecTempAdd = dat_add->GetCol(strTemp);
	t_real T_add = tl::mean_value(vecTempAdd);
	std::vector<t_real> vecCntAdd_err;


	// file to subtract
	const char* pcFileSub = argv[2];
	tl::log_info("Loading file ", pcFileSub, "...");
	tl::FileInstrBase<t_real> *dat_sub = tl::FileInstrBase<t_real>::LoadInstr(pcFileSub);
	if(!dat_sub)
	{
		tl::log_err("Cannot load data file ", pcFileSub, ".");
		return -1;
	}

	std::array<t_real, 4> pos_sub = dat_sub->GetPosHKLE();

	const tl::FileInstrBase<>::t_vecVals& vecCntSub = dat_sub->GetCol(strCnt);
	const tl::FileInstrBase<>::t_vecVals& vecMonSub = dat_sub->GetCol(strMon);
	const tl::FileInstrBase<>::t_vecVals& vecTempSub = dat_sub->GetCol(strTemp);
	t_real T_sub = tl::mean_value(vecTempSub);
	tl::log_info("T_add: ", T_add, "K, T_sub: ", T_sub, " K.");

	if(vecCntSub.size() != vecCntAdd.size() || vecMonSub.size() != vecMonAdd.size())
	{
		tl::log_err("Size mismatch in file ", pcFileSub, ".");
		return -1;
	}

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


	// output file
	const char* pcFileOut = argv[3];
	tl::log_info("Saving file ", pcFileOut);
	std::ofstream ofstr(pcFileOut);
	if(!ofstr.is_open())
	{
		tl::log_err("Cannot save data file ", pcFileOut, ".");
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
