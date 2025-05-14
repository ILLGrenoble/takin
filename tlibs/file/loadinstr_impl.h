/**
 * Loads instrument-specific data files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015 - 2023
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS_LOADINST_IMPL_H__
#define __TLIBS_LOADINST_IMPL_H__

#include "../log/log.h"
#include "../helper/py.h"
#include "../file/file.h"
#include "../math/math.h"
#include "../math/stat.h"
#include "../phys/neutrons.h"
#include "../phys/lattice.h"
#if !defined NO_IOSTR
	#include "../file/comp.h"
#endif

#include "loadinstr.h"
#include "saveinstr.h"
#include "loadinstr_psi.h"
#include "loadinstr_frm.h"
#include "loadinstr_raw.h"
#include "loadinstr_tax.h"
#include "loadinstr_macs.h"
#include "loadinstr_trisp.h"
#ifdef USE_HDF5
	#include "loadinstr_h5.h"
#endif

#ifndef USE_BOOST_REX
	#include <regex>
	namespace rex = ::std;
#else
	#include <boost/tr1/regex.hpp>
	namespace rex = ::boost;
#endif

#include <numeric>
#include <fstream>


namespace tl{

template<class t_real>
void FileInstrBase<t_real>::RenameDuplicateCols()
{
	using t_mapCols = std::unordered_map<std::string, std::size_t>;
	t_mapCols mapCols;

	t_vecColNames& vecCols = const_cast<t_vecColNames&>(this->GetColNames());
	for(std::string& strCol : vecCols)
	{
		t_mapCols::iterator iter = mapCols.find(strCol);
		if(iter == mapCols.end())
		{
			mapCols.emplace(std::make_pair(strCol, 0));
		}
		else
		{
			log_warn("Column \"", strCol, "\" is duplicate, renaming it.");

			++iter->second;
			strCol += "_" + var_to_str(iter->second);
		}
	}
}


// automatically choose correct instrument
template<class t_real>
FileInstrBase<t_real>* FileInstrBase<t_real>::LoadInstr(const char* pcFile)
{
	FileInstrBase<t_real>* pDat = nullptr;

	std::string filename = pcFile;
	std::string ext = str_to_lower(get_fileext(filename));

	// binary hdf5 files
	if(ext == "nxs" || ext == "hdf")
	{
#ifdef USE_HDF5
		pDat = new FileH5<t_real>();
#else
		return nullptr;
#endif
	}

	// text files
	else
	{
		std::ifstream ifstr(pcFile);
		if(!ifstr.is_open())
			return nullptr;

#if !defined NO_IOSTR
		std::shared_ptr<std::istream> ptrIstr = create_autodecomp_istream(ifstr);
		if(!ptrIstr) return nullptr;
		std::istream* pIstr = ptrIstr.get();
#else
		std::istream* pIstr = &ifstr;
#endif

		std::string strLine, strLine2, strLine3;
		std::getline(*pIstr, strLine);
		std::getline(*pIstr, strLine2);
		std::getline(*pIstr, strLine3);


		trim(strLine);
		trim(strLine2);
		trim(strLine3);
		strLine = str_to_lower(strLine);
		strLine2 = str_to_lower(strLine2);
		strLine3 = str_to_lower(strLine3);

		if(strLine == "")
			return nullptr;

		const std::string strNicos("nicos data file");
		const std::string strMacs("ice");
		const std::string strPsi("tas data");
		const std::string strPsiOld("instr:");
		const std::string strTax("scan =");

		if(strLine.find(strNicos) != std::string::npos)
		{ // frm file
			//log_debug(pcFile, " is an frm file.");
			pDat = new FileFrm<t_real>();
		}
		else if(strLine.find('#') != std::string::npos &&
			strLine.find(strMacs) != std::string::npos &&
			strLine2.find('#') != std::string::npos)
		{ // macs file
			//log_debug(pcFile, " is a macs file.");
			pDat = new FileMacs<t_real>();
		}
		else if(strLine2.find("scan start") != std::string::npos)
		{ // trisp file
			//log_debug(pcFile, " is a trisp file.");
			pDat = new FileTrisp<t_real>();
		}
		else if(strLine.find('#') == std::string::npos &&
			strLine2.find('#') == std::string::npos &&
			(strLine3.find(strPsi) != std::string::npos ||
			strLine.find(strPsiOld) != std::string::npos))
		{ // psi or ill file
			//log_debug(pcFile, " is an ill or psi file.");
			pDat = new FilePsi<t_real>();
		}
		else if(strLine.find('#') != std::string::npos &&
			strLine.find(strTax) != std::string::npos &&
			strLine2.find('#') != std::string::npos &&
			strLine3.find('#') != std::string::npos)
		{ // tax file
			//log_debug(pcFile, " is a tax file.");
			pDat = new FileTax<t_real>();
		}
		else
		{ // raw file
			log_warn("\"", pcFile, "\" is of unknown type, falling back to raw loader.");
			pDat = new FileRaw<t_real>();
		}
	}

	if(pDat && !pDat->Load(pcFile))
	{
		delete pDat;
		return nullptr;
	}

	return pDat;
}


template<class t_real>
std::array<t_real, 5> FileInstrBase<t_real>::GetScanHKLKiKf(const char* pcH, const char* pcK,
	const char* pcL, const char* pcE, std::size_t i) const
{
	// zero position to fallback if no position is given in scan rows
	const std::array<t_real, 4> arrZeroPos = GetPosHKLE();

	const t_vecVals& vecH = GetCol(pcH);
	const t_vecVals& vecK = GetCol(pcK);
	const t_vecVals& vecL = GetCol(pcL);
	const t_vecVals& vecE = GetCol(pcE);

	t_real h = i < vecH.size() ? vecH[i] : std::get<0>(arrZeroPos);
	t_real k = i < vecK.size() ? vecK[i] : std::get<1>(arrZeroPos);
	t_real l = i < vecL.size() ? vecL[i] : std::get<2>(arrZeroPos);
	t_real E = i < vecE.size() ? vecE[i] : std::get<3>(arrZeroPos);

	bool bKiFix = IsKiFixed();
	t_real kfix = GetKFix();
	t_real kother = get_other_k<units::si::system, t_real>
		(E*get_one_meV<t_real>(), kfix/get_one_angstrom<t_real>(), bKiFix) *
			get_one_angstrom<t_real>();

	return std::array<t_real,5>{{h,k,l, bKiFix?kfix:kother, bKiFix?kother:kfix}};
}


template<class t_real>
bool FileInstrBase<t_real>::MatchColumn(const std::string& strRegex,
	std::string& strColName, bool bSortByCounts, bool bFilterEmpty) const
{
	const FileInstrBase<t_real>::t_vecColNames& vecColNames = GetColNames();
	rex::regex rx(strRegex, rex::regex::ECMAScript | rex::regex_constants::icase);

	using t_pairCol = std::pair<std::string, t_real>;
	std::vector<t_pairCol> vecMatchedCols;

	for(const std::string& strCurColName : vecColNames)
	{
		rex::smatch m;
		if(rex::regex_match(strCurColName, m, rx))
		{
			const typename FileInstrBase<t_real>::t_vecVals& vecVals = GetCol(strCurColName);

			t_real dSum = std::accumulate(vecVals.begin(), vecVals.end(), 0.,
				[](t_real t1, t_real t2) -> t_real { return t1+t2; });
			if(!bFilterEmpty || !float_equal<t_real>(dSum, 0.))
				vecMatchedCols.push_back(t_pairCol{strCurColName, dSum});
		}
	}

	if(bSortByCounts)
	{
		std::sort(vecMatchedCols.begin(), vecMatchedCols.end(),
		[](const t_pairCol& pair1, const t_pairCol& pair2) -> bool
			{ return pair1.second > pair2.second; });
	}

	if(vecMatchedCols.size())
	{
		strColName = vecMatchedCols[0].first;
		return true;
	}
	return false;
}


template<class t_real>
bool FileInstrBase<t_real>::IsCompatible(const FileInstrBase<t_real>* pDat)
{
	if(this->GetColNames().size() != pDat->GetColNames().size())
		return false;

	for(const std::string& strCol : GetColNames())
	{
		const t_vecVals& col1 = this->GetCol(strCol);
		const t_vecVals& col2 = pDat->GetCol(strCol);

		if(col1.size() == 0 || col2.size() == 0)
			return false;
	}

	return true;
}


template<class t_real>
bool FileInstrBase<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	if(!allow_col_mismatch && this->GetColNames().size() != pDat->GetColNames().size())
	{
		log_err("Cannot merge: Mismatching number of columns.");
		return false;
	}

	for(const std::string& strCol : GetColNames())
	{
		t_vecVals& col1 = this->GetCol(strCol);
		const t_vecVals& col2 = pDat->GetCol(strCol);

		if(!allow_col_mismatch && (col1.size() == 0 || col2.size() == 0))
		{
			log_err("Cannot merge: Column \"", strCol, "\" is empty.");
			return false;
		}

		col1.insert(col1.end(), col2.begin(), col2.end());

		if(allow_col_mismatch)
		{
			// fill up to match number of rows
			for(std::size_t i = col2.size(); i < pDat->GetScanCount(); ++i)
				col1.push_back(0.);
		}
	}

	return true;
}


template<class t_real>
void FileInstrBase<t_real>::SmoothData(const std::string& strCol,
	t_real dEps, bool bIterate)
{
	std::size_t iIdxCol;
	this->GetCol(strCol, &iIdxCol);		// get column index
	if(iIdxCol == GetColNames().size())	// no such column?
	{
		log_err("No such data column: \"", strCol, "\".");
		return;
	}

	while(1)
	{
		t_vecDat& vecDatOld = this->GetData();
		const std::size_t iNumCols = vecDatOld.size();
		const std::size_t iNumRows = vecDatOld[0].size();
		t_vecDat vecDatNew(iNumCols);
		std::vector<bool> vecValidRows(iNumRows, 1);

		for(std::size_t iPt1=0; iPt1<iNumRows; ++iPt1)
		{
			if(!vecValidRows[iPt1]) continue;

			t_vecVals vecVals(iNumCols, 0);
			std::size_t iNumUnited = 0;
			for(std::size_t iPt2=iPt1; iPt2<iNumRows; ++iPt2)
			{
				if(!vecValidRows[iPt2]) continue;
				if(std::abs(vecDatOld[iIdxCol][iPt1]-vecDatOld[iIdxCol][iPt2]) <= dEps)
				{
					for(std::size_t iCol=0; iCol<iNumCols; ++iCol)
						vecVals[iCol] += vecDatOld[iCol][iPt2];
					++iNumUnited;
					vecValidRows[iPt2] = 0;
				}
			}

			for(std::size_t iCol=0; iCol<iNumCols; ++iCol)
			{
				vecVals[iCol] /= t_real(iNumUnited);
				vecDatNew[iCol].push_back(vecVals[iCol]);
			}
		}

		vecDatOld = vecDatNew;
		// if no iteration requested or no more change -> break
		if(!bIterate || iNumRows==vecDatNew[0].size())
			break;
	}
}


template<class t_real>
bool FileInstrBase<t_real>::HasCol(const std::string& strName) const
{
	const t_vecColNames& cols = GetColNames();

	for(std::size_t i=0; i<cols.size(); ++i)
	{
		if(str_to_lower(cols[i]) == str_to_lower(strName))
			return true;
	}

	return false;
}


template<class t_real>
void FileInstrBase<t_real>::ParsePolData()
{}


template<class t_real>
void FileInstrBase<t_real>::SetPolNames(const char* pVec1, const char* pVec2,
	const char* pCur1, const char* pCur2)
{}


template<class t_real>
void FileInstrBase<t_real>::SetLinPolNames(const char* pFlip1, const char* pFlip2,
	const char* pXYZ)
{}


template<class t_real>
std::size_t FileInstrBase<t_real>::NumPolChannels() const
{
	return 0;
}


template<class t_real>
const std::vector<std::array<t_real, 6>>& FileInstrBase<t_real>::GetPolStates() const
{
	static const std::vector<std::array<t_real, 6>> vecNull;
	return vecNull;
}


template<class t_real>
std::string FileInstrBase<t_real>::GetCountErr() const
{
	return "";
}


template<class t_real>
std::string FileInstrBase<t_real>::GetMonErr() const
{
	return "";
}

}

#endif
