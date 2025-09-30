/**
 * print scan files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2
 * @date oct-2017
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

// clang++ -I../.. -I. -I/usr/local/include -DNO_IOSTR -DUSE_HDF5 -o printscans ../../tlibs/file/loadinstr.cpp ../../tlibs/log/log.cpp printscans.cpp -std=c++17 -L/usr/local/lib -lboost_filesystem -lhdf5_cpp -lhdf5

#include <iostream>
#include <fstream>
#include <iomanip>

#include "tlibs/file/loadinstr.h"
#include "tlibs/log/log.h"
#include "tlibs/phys/neutrons.h"

using t_real = double;
int g_iPrec = 6;

const auto angs = tl::get_one_angstrom<t_real>();
const auto meV = tl::get_one_meV<t_real>();

int main(int argc, char **argv)
{
	std::cout.precision(g_iPrec);

	if(argc < 2)
	{
		tl::log_err("Usage:\n\t", argv[0], "<file1> <file2> ... <fileN>");
		return -1;
	}


	std::cout 
		<< std::left << std::setw(g_iPrec+2) << "Number" << " "
		<< std::left << std::setw(g_iPrec+2) << "h" << " "
		<< std::left << std::setw(g_iPrec+2) << "k" << " "
		<< std::left << std::setw(g_iPrec+2) << "l" << " "
		<< std::left << std::setw(g_iPrec+2) << "E" << " "
		<< std::left << std::setw(g_iPrec+2) << "ki" << " "
		<< std::left << std::setw(g_iPrec+2) << "kf" << " "
		<< std::left << std::setw(g_iPrec+2) << "Var" << " "
		<< std::left << std::setw(g_iPrec+2) << "T" << " "
		<< std::endl;

	for(int iArg=1; iArg<argc; ++iArg)
	{
		const char* pcFile = argv[iArg];
		//tl::log_info("Loading file ", pcFile);

		tl::FileInstrBase<t_real> *pInstr = tl::FileInstrBase<t_real>::LoadInstr(pcFile);
	        if(!pInstr)
		{
			tl::log_err("Cannot load data file ", pcFile, ".");
			continue;
		}

		std::string strT = "n/a";
		//const auto& mapParams = pInstr->GetAllParams();
		//const auto iter = mapParams.find("TT");
		//if(iter != mapParams.end())
		//	strT = iter->second;

		const tl::FileInstrBase<t_real>::t_vecVals& vecTT =
			pInstr->GetCol("TT");
		if(vecTT.size() != 0)
			strT = tl::var_to_str(vecTT[0], g_iPrec);


		const tl::FileInstrBase<t_real>::t_vecVals& vecCnt =
			pInstr->GetCol(pInstr->GetCountVar());
		if(vecCnt.size() == 0)
			continue;

		auto vecVars = pInstr->GetScannedVars();
		if(vecVars.size() == 0)
			continue;

		std::array<t_real, 5> arrPos = pInstr->GetScanHKLKiKf(0);
		if(std::get<3>(arrPos) == 0 && std::get<4>(arrPos) == 0)
			continue;

		auto E = tl::get_energy_transfer(std::get<3>(arrPos)/angs, std::get<4>(arrPos)/angs);

		std::cout 
			<< std::left << std::setw(g_iPrec+2) << pInstr->GetScanNumber() << " "
			<< std::left << std::setw(g_iPrec+2) << std::get<0>(arrPos) << " "
			<< std::left << std::setw(g_iPrec+2) << std::get<1>(arrPos) << " "
			<< std::left << std::setw(g_iPrec+2) << std::get<2>(arrPos) << " "
			<< std::left << std::setw(g_iPrec+2) << t_real(E/meV) << " "
			<< std::left << std::setw(g_iPrec+2) << std::get<3>(arrPos) << " "
			<< std::left << std::setw(g_iPrec+2) << std::get<4>(arrPos) << " "
			<< std::left << std::setw(g_iPrec+2) << vecVars[0] << " "
			<< std::left << std::setw(g_iPrec+2) << strT << " "
			<< std::endl;
	}

	return 0;
}
