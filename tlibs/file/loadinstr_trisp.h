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

#ifndef __TLIBS_LOADINST_TRISP_IMPL_H__
#define __TLIBS_LOADINST_TRISP_IMPL_H__

#include "loadinstr.h"

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

#include <numeric>
#include <fstream>


namespace tl{

template<class t_real>
void FileTrisp<t_real>::ReadHeader(std::istream& istr)
{
	bool bInVarSection = 0;
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);

		trim(strLine);
		if(strLine.length()==0)
			continue;

		if(str_contains<std::string>(strLine, "----", 0))	// new variable section beginning
		{
			bInVarSection = 1;

			if(str_contains<std::string>(strLine, "steps", 0))
				break;
			continue;
		}

		if(bInVarSection)
		{
			std::pair<std::string, std::string> pairLine =
					split_first<std::string>(strLine, " \t", 1);

			if(pairLine.first == "")
				continue;

			typename t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
		else
		{
			if(begins_with<std::string>(str_to_lower(strLine), "scan start:"))
				m_mapParams["scan_start_timestamp"] = trimmed(strLine.substr(11));
			else if(begins_with<std::string>(str_to_lower(strLine), "sc"))
				m_mapParams["scan_command"] = strLine;
		}
	}
}


template<class t_real>
void FileTrisp<t_real>::ReadData(std::istream& istr)
{
	bool bAtStepsBeginning = 0;
	bool bInFooter = 0;

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);

		if(!bAtStepsBeginning)
		{
			if(begins_with<std::string>(str_to_lower(strLine), "pnt"))
			{
				get_tokens<std::string, std::string>(strLine, " \t", m_vecQuantities);
				FileInstrBase<t_real>::RenameDuplicateCols();

				bAtStepsBeginning = 1;
				m_vecData.resize(m_vecQuantities.size());
			}
			continue;
		}

		if(strLine.length()==0 || strLine[0]=='#')
			continue;


		// character in scan data -> beginning of footer
		for(typename std::string::value_type c : split_first<std::string>(strLine, " \t", 1).first)
		{
			if(std::isalpha(c))
			{
				if(begins_with<std::string>(str_to_lower(strLine), "scan end:"))
					m_mapParams["scan_finish_timestamp"] = trimmed(strLine.substr(9));
				else if(begins_with<std::string>(str_to_lower(strLine), "scan"))
				{
					std::pair<std::string, std::string> pairLine =
						split_first<std::string>(strLine, " \t", 1);

					m_mapParams["scan_vars"] = trimmed(pairLine.second);
				}

				bInFooter = 1;
			}
		}


		if(!bInFooter)
		{
			std::vector<t_real> vecToks;
			get_tokens<t_real, std::string>(strLine, " \t", vecToks);

			if(vecToks.size() != m_vecQuantities.size())
			{
				log_warn("Loader: Line size mismatch.");

				// add zeros
				while(m_vecQuantities.size() > vecToks.size())
					vecToks.push_back(0.);
			}

			for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
				m_vecData[iTok].push_back(vecToks[iTok]);
		}
	}
}


template<class t_real>
bool FileTrisp<t_real>::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

#if !defined NO_IOSTR
	std::shared_ptr<std::istream> ptrIstr = create_autodecomp_istream(ifstr);
	if(!ptrIstr) return false;
	std::istream *pIstr = ptrIstr.get();
#else
	std::istream *pIstr = &ifstr;
#endif

	ReadHeader(*pIstr);
	ReadData(*pIstr);

	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileTrisp<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileTrisp*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileTrisp<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	for(std::size_t i=0; i<m_vecQuantities.size(); ++i)
	{
		if(m_vecQuantities[i] == strName)
		{
			if(pIdx)
				*pIdx = i;
			return m_vecData[i];
		}
	}

	if(pIdx)
		*pIdx = m_vecQuantities.size();

	log_err("Column \"", strName, "\" does not exist.");
	return vecNull;
}


template<class t_real>
std::array<t_real, 3> FileTrisp<t_real>::GetSampleLattice() const
{
	typename t_mapParams::const_iterator iterA = m_mapParams.find("AS");
	typename t_mapParams::const_iterator iterB = m_mapParams.find("BS");
	typename t_mapParams::const_iterator iterC = m_mapParams.find("CS");

	t_real a = (iterA!=m_mapParams.end() ? str_to_var<t_real>(iterA->second) : 0.);
	t_real b = (iterB!=m_mapParams.end() ? str_to_var<t_real>(iterB->second) : 0.);
	t_real c = (iterC!=m_mapParams.end() ? str_to_var<t_real>(iterC->second) : 0.);

	return std::array<t_real,3>{{a,b,c}};
}


template<class t_real>
std::array<t_real, 3> FileTrisp<t_real>::GetSampleAngles() const
{
	typename t_mapParams::const_iterator iterA = m_mapParams.find("AA");
	typename t_mapParams::const_iterator iterB = m_mapParams.find("BB");
	typename t_mapParams::const_iterator iterC = m_mapParams.find("CC");

	t_real alpha = (iterA!=m_mapParams.end() ? d2r(str_to_var<t_real>(iterA->second)) : get_pi<t_real>()/2.);
	t_real beta = (iterB!=m_mapParams.end() ? d2r(str_to_var<t_real>(iterB->second)) : get_pi<t_real>()/2.);
	t_real gamma = (iterC!=m_mapParams.end() ? d2r(str_to_var<t_real>(iterC->second)) : get_pi<t_real>()/2.);

	return std::array<t_real,3>{{alpha, beta, gamma}};
}


template<class t_real>
std::array<t_real, 2> FileTrisp<t_real>::GetMonoAnaD() const
{
	typename t_mapParams::const_iterator iterM = m_mapParams.find("DM");
	typename t_mapParams::const_iterator iterA = m_mapParams.find("DA");

	t_real m = (iterM!=m_mapParams.end() ? str_to_var<t_real>(iterM->second) : 3.355);
	t_real a = (iterA!=m_mapParams.end() ? str_to_var<t_real>(iterA->second) : 3.355);

	return std::array<t_real,2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FileTrisp<t_real>::GetScatterSenses() const
{
	typename t_mapParams::const_iterator iterM = m_mapParams.find("SM");
	typename t_mapParams::const_iterator iterS = m_mapParams.find("SS");
	typename t_mapParams::const_iterator iterA = m_mapParams.find("SA");

	bool m = (iterM!=m_mapParams.end() ? (str_to_var<int>(iterM->second)>0) : 0);
	bool s = (iterS!=m_mapParams.end() ? (str_to_var<int>(iterS->second)>0) : 1);
	bool a = (iterA!=m_mapParams.end() ? (str_to_var<int>(iterA->second)>0) : 0);

	return std::array<bool,3>{{m, s, a}};
}


template<class t_real>
std::array<t_real, 3> FileTrisp<t_real>::GetScatterPlane0() const
{
	typename t_mapParams::const_iterator iterX = m_mapParams.find("AX");
	typename t_mapParams::const_iterator iterY = m_mapParams.find("AY");
	typename t_mapParams::const_iterator iterZ = m_mapParams.find("AZ");

	t_real x = (iterX!=m_mapParams.end() ? str_to_var<t_real>(iterX->second) : 1.);
	t_real y = (iterY!=m_mapParams.end() ? str_to_var<t_real>(iterY->second) : 0.);
	t_real z = (iterZ!=m_mapParams.end() ? str_to_var<t_real>(iterZ->second) : 0.);

	return std::array<t_real,3>{{x,y,z}};
}


template<class t_real>
std::array<t_real, 3> FileTrisp<t_real>::GetScatterPlane1() const
{
	typename t_mapParams::const_iterator iterX = m_mapParams.find("BX");
	typename t_mapParams::const_iterator iterY = m_mapParams.find("BY");
	typename t_mapParams::const_iterator iterZ = m_mapParams.find("BZ");

	t_real x = (iterX!=m_mapParams.end() ? str_to_var<t_real>(iterX->second) : 0.);
	t_real y = (iterY!=m_mapParams.end() ? str_to_var<t_real>(iterY->second) : 1.);
	t_real z = (iterZ!=m_mapParams.end() ? str_to_var<t_real>(iterZ->second) : 0.);

	return std::array<t_real,3>{{x,y,z}};
}


template<class t_real>
std::array<t_real, 4> FileTrisp<t_real>::GetPosHKLE() const
{
	t_real h = 0;
	t_real k = 0;
	t_real l = 0;
	t_real E = 0;

	// get the first position from the scan if available
	const t_vecVals& vecH = GetCol("QH");
	const t_vecVals& vecK = GetCol("QK");
	const t_vecVals& vecL = GetCol("QL");
	const t_vecVals& vecE = GetCol("E");

	if(vecH.size()) h = vecH[0];
	if(vecK.size()) k = vecK[0];
	if(vecL.size()) l = vecL[0];
	if(vecE.size()) E = vecE[0];

	return std::array<t_real, 4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileTrisp<t_real>::GetKFix() const
{
	const std::string strKey = IsKiFixed() ? "KI" : "KF";

	typename t_mapParams::const_iterator iter = m_mapParams.find(strKey);
	if(iter==m_mapParams.end())
	{
		log_err("Cannot determine kfix.");
		return 0.;
	}

	return str_to_var<t_real>(iter->second);
}


template<class t_real>
bool FileTrisp<t_real>::IsKiFixed() const
{
	return false;		// assume ckf
}


template<class t_real>
std::size_t FileTrisp<t_real>::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}


template<class t_real>
std::array<t_real, 5> FileTrisp<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstrBase<t_real>::GetScanHKLKiKf("QH", "QK", "QL", "E", i);
}


template<class t_real>
bool FileTrisp<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	return FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch);
}

template<class t_real> std::string FileTrisp<t_real>::GetTitle() const { return ""; }
template<class t_real> std::string FileTrisp<t_real>::GetProposal() const { return ""; }
template<class t_real> std::string FileTrisp<t_real>::GetUser() const { return ""; }
template<class t_real> std::string FileTrisp<t_real>::GetInstrument() const { return "trisp"; }
template<class t_real> std::string FileTrisp<t_real>::GetLocalContact() const { return ""; }
template<class t_real> std::string FileTrisp<t_real>::GetScanNumber() const { return ""; }
template<class t_real> std::string FileTrisp<t_real>::GetSampleName() const { return ""; }
template<class t_real> std::string FileTrisp<t_real>::GetSpacegroup() const { return ""; }


template<class t_real>
std::vector<std::string> FileTrisp<t_real>::GetScannedVars() const
{
	std::vector<std::string> vecScan;

	typename t_mapParams::const_iterator iter = m_mapParams.find("scan_vars");
	if(iter != m_mapParams.end())
		get_tokens<std::string, std::string>(iter->second, " \t", vecScan);

	if(!vecScan.size())
	{
		log_warn("Could not determine scan variable.");
		if(m_vecQuantities.size() >= 1)
		{
			log_warn("Using first column: \"", m_vecQuantities[0], "\".");
			vecScan.push_back(m_vecQuantities[0]);
		}
	}

	return vecScan;
}


template<class t_real>
std::string FileTrisp<t_real>::GetCountVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(c[0-9])REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(c[0-9])REX", strRet, false, false))
		return strRet;
	return "";
}


template<class t_real>
std::string FileTrisp<t_real>::GetMonVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(mon[a-z0-9]*)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(mon[a-z0-9]*)REX", strRet, false, false))
		return strRet;
	return "";
}


template<class t_real>
std::string FileTrisp<t_real>::GetScanCommand() const
{
	std::string str;
	typename t_mapParams::const_iterator iter = m_mapParams.find("scan_command");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}


template<class t_real>
std::string FileTrisp<t_real>::GetTimestamp() const
{
	std::string str;
	typename t_mapParams::const_iterator iter = m_mapParams.find("scan_start_timestamp");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}

}


#endif
