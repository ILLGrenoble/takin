/**
 * saves instrument-specific data files
 * @author Tobias Weber <tweber@ill.fr>
 * @date may-2025
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS_SAVEINST_IMPL_H__
#define __TLIBS_SAVEINST_IMPL_H__

#include "loadinstr.h"
#include <iostream>


namespace tl{


/**
 * convert scan file to the psi text file format
 * TODO: also add the other formats
 */
template<class t_real>
bool FileInstrBase<t_real>::Save(std::ostream& ostr)
{
	// output header
	ostr << "USER_: " << GetUser() << "\n";
	ostr << "LOCAL: " << GetLocalContact() << "\n";
	ostr << "FILE_: " << GetScanNumber() << "\n";
	ostr << "DATE_: " << GetTimestamp() << "\n";
	ostr << "TITLE: " << GetTitle() << "\n";
	ostr << "POSQE: QH = " << std::get<0>(GetPosHKLE())
		<< ", QK = " << std::get<1>(GetPosHKLE())
		<< ", QL = " << std::get<2>(GetPosHKLE())
		<< ", EN = " << std::get<3>(GetPosHKLE())
		<< ", UN = meV" << "\n";
	ostr << "COMND: " << GetScanCommand() << "\n";

	ostr << "PARAM: DM = " << std::get<0>(GetMonoAnaD())
		<< ", DA = " << std::get<1>(GetMonoAnaD())
		<< ", KFIX = " << GetKFix() << "\n";
	ostr << "PARAM: SM = " << (std::get<0>(GetScatterSenses()) ? "1" : "-1")
		<< ", SS = " << (std::get<1>(GetScatterSenses()) ? "1" : "-1")
		<< ", SA = " << (std::get<2>(GetScatterSenses()) ? "1" : "-1")
		<< ", FX = " << (IsKiFixed() ? "1" : "2") << "\n";
	ostr << "PARAM: AS = " << std::get<0>(GetSampleLattice())
		<< ", BS = " << std::get<1>(GetSampleLattice())
		<< ", CS = " << std::get<2>(GetSampleLattice()) << "\n";
	ostr << "PARAM: AA = " << tl::r2d(std::get<0>(GetSampleAngles()))
		<< ", BB = " << tl::r2d(std::get<1>(GetSampleAngles()))
		<< ", CC = " << tl::r2d(std::get<2>(GetSampleAngles())) << "\n";
	ostr << "PARAM: AX = " << std::get<0>(GetScatterPlane0())
		<<", AY = " << std::get<1>(GetScatterPlane0())
		<<", AZ = " << std::get<2>(GetScatterPlane0()) << "\n";
	ostr << "PARAM: BX = " << std::get<0>(GetScatterPlane1())
		<<", BY = " << std::get<1>(GetScatterPlane1())
		<<", BZ = " << std::get<2>(GetScatterPlane1()) << "\n";

	// output data columns
	const t_vecColNames& colnames = GetColNames();

	ostr << "FORMT:\n";
	ostr << "DATA_:\n";

	const int colwidth = std::max<int>(ostr.precision() * 2.5, 15);

	// column names
	for(const std::string& colname : colnames)
	{
		if(colname == "Point_Index")
			ostr << std::setw(colwidth) << std::left << "PNT" << " ";
		else
			ostr << std::setw(colwidth) << std::left << colname << " ";
	}
	ostr << "\n";

	// data
	const std::size_t rows = GetScanCount();
	const std::size_t cols = colnames.size();
	const t_vecDat& data = GetData();

	for(std::size_t row = 0; row < rows; ++row)
	{
		for(std::size_t col = 0; col < cols; ++col)
			ostr << std::setw(colwidth) << std::left << data[col][row] << " ";

		ostr << "\n";
	}

	return true;
}


}

#endif
