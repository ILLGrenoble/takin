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


static t_real get_err(t_real cnt)
{
	t_real err = std::sqrt(cnt);
	if(tl::float_equal<t_real>(err, 0.))
		err = 1.;
	return err;
}


int main(int argc, char **argv)
{
	if(argc < 3)
	{
		tl::log_err("Usage:\n\t", argv[0], " <file_add> <file_sub> <file_out>");
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
	std::vector<t_real> vecCntAdd_err, vecMonAdd_err;


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

		t_real err_add = get_err(vecCntAdd[i]);
		t_real err_sub = get_err(vecCntSub[i]);
		t_real err_mon_add = get_err(vecMonAdd[i]);
		t_real err_mon_sub = get_err(vecMonSub[i]);

		vecCntAdd[i] -= vecCntSub[i] * vecMonAdd[i]/vecMonSub[i] * bose_add/bose_sub;

		err_sub *= vecMonAdd[i]/vecMonSub[i];
		vecCntAdd_err.push_back(std::sqrt(err_add*err_add + err_sub*err_sub));
		vecMonAdd_err.push_back(std::sqrt(err_mon_add*err_mon_add + err_mon_sub*err_mon_sub));
	}
	delete dat_sub;


	// output file
	const char* pcFileOut = argv[3];
	tl::log_info("Saving file ", pcFileOut);
	std::ofstream ofstr(pcFileOut);
	if(!ofstr.is_open())
	{
		tl::log_err("Cannot save output data file ", pcFileOut, ".");
		return -1;
	}

	int w = 15;
	ofstr.precision(8);


	// header
	std::array<t_real, 3> lattice = dat_add->GetSampleLattice();
	ofstr << "#\n";
	ofstr << "# sample_a     = " << std::get<0>(lattice) << "\n";
	ofstr << "# sample_b     = " << std::get<1>(lattice) << "\n";
	ofstr << "# sample_c     = " << std::get<2>(lattice) << "\n";

	std::array<t_real, 3> angles = dat_add->GetSampleAngles();
	ofstr << "# sample_alpha = " << tl::r2d(std::get<0>(angles)) << "\n";
	ofstr << "# sample_beta  = " << tl::r2d(std::get<1>(angles)) << "\n";
	ofstr << "# sample_gamma = " << tl::r2d(std::get<2>(angles)) << "\n";

	std::array<t_real, 3> vec0 = dat_add->GetScatterPlane0();
	ofstr << "# orient1_x    = " << std::get<0>(vec0) << "\n";
	ofstr << "# orient1_y    = " << std::get<1>(vec0) << "\n";
	ofstr << "# orient1_z    = " << std::get<2>(vec0) << "\n";

	std::array<t_real, 3> vec1 = dat_add->GetScatterPlane0();
	ofstr << "# orient2_x    = " << std::get<0>(vec1) << "\n";
	ofstr << "# orient2_y    = " << std::get<1>(vec1) << "\n";
	ofstr << "# orient2_z    = " << std::get<2>(vec1) << "\n";

	std::array<t_real, 2> ds = dat_add->GetMonoAnaD();
	ofstr << "# mono_d       = " << std::get<0>(ds) << "\n";
	ofstr << "# ana_d        = " << std::get<1>(ds) << "\n";

	std::array<bool, 3> senses = dat_add->GetScatterSenses();
	ofstr << "# sense_m      = " << (std::get<0>(senses) ? "+1": "-1") << "\n";
	ofstr << "# sense_s      = " << (std::get<1>(senses) ? "+1": "-1") << "\n";
	ofstr << "# sense_a      = " << (std::get<2>(senses) ? "+1": "-1") << "\n";

	ofstr << "# is_ki_fixed  = " << (dat_add->IsKiFixed() ? "1" : "0") << "\n";
	ofstr << "# k_fix        = " << dat_add->GetKFix() << "\n";

	ofstr << "# col_h        = 1" << "\n";
	ofstr << "# col_k        = 2" << "\n";
	ofstr << "# col_l        = 3" << "\n";
	ofstr << "# col_E        = 4" << "\n";
	ofstr << "# col_ctr      = 5" << "\n";
	ofstr << "# col_ctr_err  = 6" << "\n";
	ofstr << "# col_mon      = 7" << "\n";
	ofstr << "# col_mon_err  = 8" << "\n";
	ofstr << "# col_timer    = -1" << "\n";
	ofstr << "# cols_scanned = 4" << "\n";
	ofstr << "#\n";


	// write data columns
	for(unsigned int i = 0; i < vecCntAdd.size(); ++i)
	{
		/*std::array<t_real, 5> hklkikf = dat_add->GetScanHKLKiKf(i);
		t_real E = (tl::k2E(std::get<3>(hklkikf) / tl::get_one_angstrom<t_real>())
			- tl::k2E(std::get<4>(hklkikf) / tl::get_one_angstrom<t_real>())) / tl::get_one_meV<t_real>();
		std::cout << E << " " << std::get<3>(hklE) << std::endl;*/
		std::array<t_real, 4> hklE = dat_add->GetScanHKLE(i);

		ofstr << std::setw(w) << std::right << std::get<0>(hklE) << " ";
		ofstr << std::setw(w) << std::right << std::get<1>(hklE) << " ";
		ofstr << std::setw(w) << std::right << std::get<2>(hklE) << " ";
		ofstr << std::setw(w) << std::right << std::get<3>(hklE) << " ";
		ofstr << std::setw(w) << std::right << vecCntAdd[i] << " ";
		ofstr << std::setw(w) << std::right << vecCntAdd_err[i] << " ";
		ofstr << std::setw(w) << std::right << vecMonAdd[i] << " ";
		ofstr << std::setw(w) << std::right << vecMonAdd_err[i] << "\n";
	}

	/*const tl::FileInstrBase<>::t_vecColNames& vecColNames = dat_add->GetColNames();
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
	}*/

	delete dat_add;
	return 0;
}
