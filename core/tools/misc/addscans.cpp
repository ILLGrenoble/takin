/**
 * add scan files
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
// g++ -I../.. -I. -DNO_IOSTR -o addscans ../../tlibs/file/loadinstr.cpp ../../tlibs/log/log.cpp addscans.cpp -std=c++17 -lboost_system -lboost_filesystem
// e.g. ./addscans /home/tweber/Auswertungen/MnSi-Mira-15/data3/11009_00016851.dat /home/tweber/Auswertungen/MnSi-Mira-15/data3/11009_00016867.dat merged.dat


#include <iostream>
#include <fstream>
#include <iomanip>
#include "tlibs/file/loadinstr.h"
#include "tlibs/log/log.h"
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
	bool write_raw_format = true;
	t_real Q_eps = 0.01;
	t_real E_eps = 0.05;
	int w = 15;    // output column width
	int prec = 8;  // output precision


	if(argc < 3)
	{
		tl::log_err("Usage:\n\t", argv[0], " <file1> <file2> ... <fileN> <file_out>");
		return -1;
	}

	tl::log_info("Loading file ", argv[1], "...");
	tl::FileInstrBase<t_real> *dat0 = tl::FileInstrBase<t_real>::LoadInstr(argv[1]);
	if(!dat0)
	{
		tl::log_err("Cannot load data file ", argv[1], ".");
		return -1;
	}


	std::string strCnt = dat0->GetCountVar();
	std::string strMon = dat0->GetMonVar();
	std::string strTime = dat0->GetTimerVar();
	// additional columns to add
	std::vector<std::string> strAdditional = { /*"M2"*/ };
	//std::string strMon = "M2";
	tl::log_info("Count var: ", strCnt, ", monitor var: ", strMon, ", timer var: ", strTime);

	tl::FileInstrBase<>::t_vecVals& vecCnt0 = dat0->GetCol(strCnt);
	tl::FileInstrBase<>::t_vecVals& vecMon0 = dat0->GetCol(strMon);
	tl::FileInstrBase<>::t_vecVals& vecTime0 = dat0->GetCol(strTime);
	std::vector<tl::FileInstrBase<>::t_vecVals*> vecAdditional0;
	for(const std::string& strAdd : strAdditional)
		vecAdditional0.push_back(&dat0->GetCol(strAdd));


	for(int iArg = 2; iArg < argc - 1; ++iArg)
	{
		const char* pcFile = argv[iArg];

		tl::log_info("Loading file ", pcFile, "...");
		tl::FileInstrBase<t_real> *dat = tl::FileInstrBase<t_real>::LoadInstr(pcFile);
		if(!dat)
		{
			tl::log_err("Cannot load data file ", pcFile, ".");
			return -1;
		}

		const tl::FileInstrBase<>::t_vecVals& vecCnt = dat->GetCol(strCnt);
		const tl::FileInstrBase<>::t_vecVals& vecMon = dat->GetCol(strMon);
		const tl::FileInstrBase<>::t_vecVals& vecTime = dat->GetCol(strTime);
		std::vector<const tl::FileInstrBase<>::t_vecVals*> vecAdditional;
		for(const std::string& strAdd : strAdditional)
			vecAdditional.push_back(&dat->GetCol(strAdd));

		if(vecCnt.size() != vecCnt0.size() || vecMon.size() != vecMon0.size() || vecTime.size() != vecTime0.size())
		{
			tl::log_err("Size mismatch in file ", pcFile, ".");
			return -1;
		}

		for(unsigned int i = 0; i < vecCnt0.size(); ++i)
		{
			std::array<t_real, 4> pos0 = dat0->GetScanHKLE(i);
			std::array<t_real, 4> pos = dat->GetScanHKLE(i);

			if(!tl::float_equal<t_real>(std::get<0>(pos0), std::get<0>(pos), Q_eps))
				tl::log_warn("Mismatching h in index ", i, ": ", std::get<0>(pos0), " != ", std::get<0>(pos), ".");
			if(!tl::float_equal<t_real>(std::get<1>(pos0), std::get<1>(pos), Q_eps))
				tl::log_warn("Mismatching k in index ", i, ": ", std::get<1>(pos0), " != ", std::get<1>(pos), ".");
			if(!tl::float_equal<t_real>(std::get<2>(pos0), std::get<2>(pos), Q_eps))
				tl::log_warn("Mismatching l in index ", i, ": ", std::get<2>(pos0), " != ", std::get<2>(pos), ".");
			if(!tl::float_equal<t_real>(std::get<3>(pos0), std::get<3>(pos), E_eps))
				tl::log_warn("Mismatching E in index ", i, ": ", std::get<3>(pos0), " != ", std::get<3>(pos), ".");

			vecCnt0[i] += vecCnt[i];
			vecMon0[i] += vecMon[i];
			vecTime0[i] += vecTime[i];

			for(unsigned int iAdd = 0; iAdd < vecAdditional.size(); ++iAdd)
				(*vecAdditional0[iAdd])[i] += (*vecAdditional[iAdd])[i];
		}

		delete dat;
	}


	tl::log_info("Saving file ", argv[argc - 1]);
	std::ofstream ofstr(argv[argc - 1]);
	if(!ofstr.is_open())
	{
		tl::log_err("Cannot save data file ", argv[argc - 1], ".");
		return -1;
	}

	ofstr.precision(prec);


	if(write_raw_format)
	{
		// header
		std::array<t_real, 3> lattice = dat0->GetSampleLattice();
		ofstr << "#\n";
		ofstr << "# sample_a     = " << std::get<0>(lattice) << "\n";
		ofstr << "# sample_b     = " << std::get<1>(lattice) << "\n";
		ofstr << "# sample_c     = " << std::get<2>(lattice) << "\n";

		std::array<t_real, 3> angles = dat0->GetSampleAngles();
		ofstr << "# sample_alpha = " << tl::r2d(std::get<0>(angles)) << "\n";
		ofstr << "# sample_beta  = " << tl::r2d(std::get<1>(angles)) << "\n";
		ofstr << "# sample_gamma = " << tl::r2d(std::get<2>(angles)) << "\n";

		std::array<t_real, 3> vec0 = dat0->GetScatterPlane0();
		ofstr << "# orient1_x    = " << std::get<0>(vec0) << "\n";
		ofstr << "# orient1_y    = " << std::get<1>(vec0) << "\n";
		ofstr << "# orient1_z    = " << std::get<2>(vec0) << "\n";

		std::array<t_real, 3> vec1 = dat0->GetScatterPlane1();
		ofstr << "# orient2_x    = " << std::get<0>(vec1) << "\n";
		ofstr << "# orient2_y    = " << std::get<1>(vec1) << "\n";
		ofstr << "# orient2_z    = " << std::get<2>(vec1) << "\n";

		std::array<t_real, 2> ds = dat0->GetMonoAnaD();
		ofstr << "# mono_d       = " << std::get<0>(ds) << "\n";
		ofstr << "# ana_d        = " << std::get<1>(ds) << "\n";

		std::array<bool, 3> senses = dat0->GetScatterSenses();
		ofstr << "# sense_m      = " << (std::get<0>(senses) ? "+1": "-1") << "\n";
		ofstr << "# sense_s      = " << (std::get<1>(senses) ? "+1": "-1") << "\n";
		ofstr << "# sense_a      = " << (std::get<2>(senses) ? "+1": "-1") << "\n";

		ofstr << "# is_ki_fixed  = " << (dat0->IsKiFixed() ? "1" : "0") << "\n";
		ofstr << "# k_fix        = " << dat0->GetKFix() << "\n";

		ofstr << "# col_h        = 1" << "\n";
		ofstr << "# col_k        = 2" << "\n";
		ofstr << "# col_l        = 3" << "\n";
		ofstr << "# col_E        = 4" << "\n";
		ofstr << "# col_ctr      = 5" << "\n";
		ofstr << "# col_ctr_err  = 6" << "\n";
		ofstr << "# col_mon      = 7" << "\n";
		ofstr << "# col_mon_err  = 8" << "\n";
		ofstr << "# col_timer    = 9" << "\n";
		ofstr << "# col_timer    = -1" << "\n";
		ofstr << "# cols_scanned = 4" << "\n";
		ofstr << "#\n";

		// column header
		ofstr << "#" << std::setw(w - 1) << std::right << "h" << " ";
		ofstr << std::setw(w) << std::right << "k" << " ";
		ofstr << std::setw(w) << std::right << "l" << " ";
		ofstr << std::setw(w) << std::right << "E" << " ";
		ofstr << std::setw(w) << std::right << "ctr" << " ";
		ofstr << std::setw(w) << std::right << "ctr_err" << " ";
		ofstr << std::setw(w) << std::right << "mon" << " ";
		ofstr << std::setw(w) << std::right << "mon_err" << " ";
		ofstr << std::setw(w) << std::right << "time" << "\n";

		// write data columns
		for(unsigned int i = 0; i < vecCnt0.size(); ++i)
		{
			std::array<t_real, 4> hklE = dat0->GetScanHKLE(i);

			ofstr << std::setw(w) << std::right << std::get<0>(hklE) << " ";
			ofstr << std::setw(w) << std::right << std::get<1>(hklE) << " ";
			ofstr << std::setw(w) << std::right << std::get<2>(hklE) << " ";
			ofstr << std::setw(w) << std::right << std::get<3>(hklE) << " ";
			ofstr << std::setw(w) << std::right << vecCnt0[i] << " ";
			ofstr << std::setw(w) << std::right << get_err(vecCnt0[i]) << " ";
			ofstr << std::setw(w) << std::right << vecMon0[i] << " ";
			ofstr << std::setw(w) << std::right << get_err(vecMon0[i]) << " ";
			ofstr << std::setw(w) << std::right << vecTime0[i] << "\n";
		}
	}
	else
	{
		const tl::FileInstrBase<>::t_vecColNames& vecColNames = dat0->GetColNames();
		for(const std::string& strColName : vecColNames)
			ofstr << std::setw(w) << std::right << strColName << " ";
		ofstr << "\n";

		for(unsigned int i = 0; i < vecCnt0.size(); ++i)
		{
			for(const std::string& strColName : vecColNames)
			{
				tl::FileInstrBase<>::t_vecVals& vecCol = dat0->GetCol(strColName);
				ofstr << std::setw(w) << std::right << vecCol[i] << " ";
			}
			ofstr << "\n";
		}
	}

	delete dat0;
	return 0;
}
