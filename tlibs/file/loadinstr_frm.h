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

#ifndef __TLIBS_LOADINST_FRM_IMPL_H__
#define __TLIBS_LOADINST_FRM_IMPL_H__

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

#include <regex>
#include <numeric>
#include <fstream>


namespace tl{

template<class t_real>
void FileFrm<t_real>::ReadHeader(std::istream& istr)
{
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);
		if(strLine.length()==0 || strLine[0]!='#')
			continue;
		if(strLine.length()>=3 && strLine[0]=='#' && strLine[1]=='#' && strLine[2]=='#')
		{
			std::string strCreatedAt("created at");
			std::size_t iPosCreated = strLine.find(strCreatedAt);
			if(iPosCreated != std::string::npos)
			{
				iPosCreated += strCreatedAt.length();
				std::string strDate = strLine.substr(iPosCreated);
				trim(strDate);

				m_mapParams["file_timestamp"] = strDate;
			}

			continue;
		}

		strLine = strLine.substr(1);

		std::pair<std::string, std::string> pairLine =
			split_first<std::string>(strLine, ":", true);

		if(pairLine.first == "")
		{
			continue;
		}
		else
		{
			typename t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;

			// try to find instrument name
			if(m_instrname == "")
			{
				try
				{
					const std::string strRegex = R"REX(([a-z0-9]+)_responsible)REX";
					std::regex rx(strRegex, std::regex::ECMAScript|std::regex_constants::icase);
					std::smatch m;
					if(std::regex_search(pairLine.first, m, rx) && m.size()>=2)
						m_instrname = m[1];
				}
				catch(const std::exception& ex)
				{
					log_warn(ex.what());
				}
			}
		}
	}
}


template<class t_real>
void FileFrm<t_real>::ReadData(std::istream& istr)
{
	skip_after_line<char>(istr, "### scan data", true, false);

	// column headers
	skip_after_char<char>(istr, '#');
	std::string strLineQuantities;
	std::getline(istr, strLineQuantities);
	get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t", m_vecQuantities);
	for(std::string& _str : m_vecQuantities)
	{
		find_all_and_replace<std::string>(_str, "\n", "");
		find_all_and_replace<std::string>(_str, "\r", "");
	}

	skip_after_char<char>(istr, '#');
	std::string strLineUnits;
	std::getline(istr, strLineUnits);
	get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t", m_vecUnits);


	m_vecData.resize(m_vecQuantities.size());

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);
		if(strLine.length()==0 || strLine[0]=='#')
			continue;

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

	FileInstrBase<t_real>::RenameDuplicateCols();
}


template<class t_real>
bool FileFrm<t_real>::Load(const char* pcFile)
{
	for(int iStep : {0,1})
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

		if(iStep==0)
			ReadHeader(*pIstr);
		else if(iStep==1)
			ReadData(*pIstr);
	}

	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileFrm<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileFrm*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileFrm<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
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
std::array<t_real, 3> FileFrm<t_real>::GetSampleLattice() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Sample_lattice");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{0.,0.,0.}};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid lattice array size.");
		return std::array<t_real,3>{{0.,0.,0.}};
	}

	return std::array<t_real,3>{{vec[0],vec[1],vec[2]}};
}


template<class t_real>
std::array<t_real, 3> FileFrm<t_real>::GetSampleAngles() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Sample_angles");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{0.,0.,0.}};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid angle array size.");
		return std::array<t_real,3>{{0.,0.,0.}};
	}

	return std::array<t_real,3>{{d2r(vec[0]), d2r(vec[1]), d2r(vec[2])}};
}


template<class t_real>
std::array<t_real, 2> FileFrm<t_real>::GetMonoAnaD() const
{
	typename t_mapParams::const_iterator iterM = m_mapParams.find("mono_dvalue");
	typename t_mapParams::const_iterator iterA = m_mapParams.find("ana_dvalue");

	t_real m = (iterM!=m_mapParams.end() ? str_to_var<t_real>(iterM->second) : 3.355);
	t_real a = (iterA!=m_mapParams.end() ? str_to_var<t_real>(iterA->second) : 3.355);

	return std::array<t_real,2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FileFrm<t_real>::GetScatterSenses() const
{
	std::vector<int> vec;

	typename t_mapParams::const_iterator iter;
	for(iter=m_mapParams.begin(); iter!=m_mapParams.end(); ++iter)
	{
		if(iter->first.find("scatteringsense") != std::string::npos)
		{
			vec = get_py_array<std::string, std::vector<int>>(iter->second);
			break;
		}
	}

	if(vec.size() != 3)
	{
		vec.resize(3);
		vec[0] = 0; vec[1] = 1; vec[2] = 0;
	}

	return std::array<bool,3>{{vec[0]>0, vec[1]>0, vec[2]>0}};
}


template<class t_real>
std::array<t_real, 3> FileFrm<t_real>::GetScatterPlane0() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Sample_orient1");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{0.,0.,0.}};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid sample peak 1 array size.");
		return std::array<t_real,3>{{0.,0.,0.}};
	}
	return std::array<t_real,3>{{vec[0],vec[1],vec[2]}};
}


template<class t_real>
std::array<t_real, 3> FileFrm<t_real>::GetScatterPlane1() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("Sample_orient2");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{0.,0.,0.}};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid sample peak 2 array size.");
		return std::array<t_real,3>{{0.,0.,0.}};
	}
	return std::array<t_real,3>{{-vec[0],-vec[1],-vec[2]}};	// LH -> RH
}


template<class t_real>
std::array<t_real, 4> FileFrm<t_real>::GetPosHKLE() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find(
		m_instrname + "_value");
	if(iter == m_mapParams.end())
		return std::array<t_real,4>{{0, 0, 0, 0}};

	std::vector<t_real> vecPos = get_py_array<std::string, std::vector<t_real>>(
		iter->second);
	if(vecPos.size() < 4)
		return std::array<t_real,4>{{0, 0, 0, 0}};

	return std::array<t_real,4>{{vecPos[0], vecPos[1], vecPos[2], vecPos[3]}};
}


template<class t_real>
t_real FileFrm<t_real>::GetKFix() const
{
	std::string strKey = (IsKiFixed() ? "ki_value" : "kf_value");

	typename t_mapParams::const_iterator iter = m_mapParams.find(strKey);
	return (iter!=m_mapParams.end() ? str_to_var<t_real>(iter->second) : 0.);
}


template<class t_real>
bool FileFrm<t_real>::IsKiFixed() const
{
	std::string strScanMode = "ckf";

	typename t_mapParams::const_iterator iter;
	for(iter=m_mapParams.begin(); iter!=m_mapParams.end(); ++iter)
	{
		if(iter->first.find("scanmode") != std::string::npos)
		{
			strScanMode = str_to_lower(iter->second);
			trim(strScanMode);
			break;
		}
	}

	if(strScanMode == "cki")
		return true;
	return false;
}


template<class t_real>
std::size_t FileFrm<t_real>::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}


template<class t_real>
std::array<t_real, 5> FileFrm<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstrBase<t_real>::GetScanHKLKiKf("h", "k", "l", "E", i);
}


template<class t_real>
bool FileFrm<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	if(!FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch))
		return false;

	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		typename t_mapParams::iterator iter = m_mapParams.find("number");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;
}


template<class t_real>
std::string FileFrm<t_real>::GetTitle() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Exp_title");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileFrm<t_real>::GetProposal() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Exp_proposal");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileFrm<t_real>::GetUser() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Exp_users");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FileFrm<t_real>::GetInstrument() const
{
	return m_instrname;
}


template<class t_real>
std::string FileFrm<t_real>::GetLocalContact() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Exp_localcontact");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FileFrm<t_real>::GetScanNumber() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("number");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileFrm<t_real>::GetSampleName() const
{
	std::string strName;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Sample_samplename");
	if(iter != m_mapParams.end())
		strName = iter->second;
	return strName;
}


template<class t_real>
std::string FileFrm<t_real>::GetSpacegroup() const
{
	std::string strSG;
	typename t_mapParams::const_iterator iter = m_mapParams.find("Sample_spacegroup");
	if(iter != m_mapParams.end())
		strSG = iter->second;
	return strSG;
}


template<class t_real>
std::vector<std::string> FileFrm<t_real>::GetScannedVars() const
{
	std::vector<std::string> vecVars;

	// scan command
	typename t_mapParams::const_iterator iter = m_mapParams.find("info");
	if(iter != m_mapParams.end())
	{
		const std::string& strInfo = iter->second;

		// try qscan/qcscan
		const std::string strRegex = R"REX((qscan|qcscan)\((\[.*\])[, ]+(\[.*\]).*\))REX";
		std::regex rx(strRegex, std::regex::ECMAScript|std::regex_constants::icase);
		std::smatch m;
		if(std::regex_search(strInfo, m, rx) && m.size()>3)
		{
			const std::string& strSteps = m[3];
			std::vector<t_real> vecSteps = get_py_array<std::string, std::vector<t_real>>(strSteps);

			if(vecSteps.size()>0 && !float_equal<t_real>(vecSteps[0], 0.))
				vecVars.push_back("h");
			if(vecSteps.size()>1 && !float_equal<t_real>(vecSteps[1], 0.))
				vecVars.push_back("k");
			if(vecSteps.size()>2 && !float_equal<t_real>(vecSteps[2], 0.))
				vecVars.push_back("l");
			if(vecSteps.size()>3 && !float_equal<t_real>(vecSteps[3], 0.))
				vecVars.push_back("E");
		}


		if(vecVars.size() == 0)
		{
			// try scan/cscan
			const std::string strRegexDevScan = R"REX((scan|cscan)\(([a-z0-9_\.]+)[, ]+.*\))REX";
			std::regex rxDev(strRegexDevScan, std::regex::ECMAScript|std::regex_constants::icase);
			std::smatch mDev;
			if(std::regex_search(strInfo, mDev, rxDev) && mDev.size()>2)
			{
				std::string strDev = mDev[2];
				if(std::find(m_vecQuantities.begin(), m_vecQuantities.end(), strDev) != m_vecQuantities.end())
					vecVars.push_back(strDev);
			}
		}
	}

	if(!vecVars.size())
	{
		log_warn("Could not determine scan variable.");
		if(m_vecQuantities.size() >= 1)
		{
			log_warn("Using first column: \"", m_vecQuantities[0], "\".");
			vecVars.push_back(m_vecQuantities[0]);
		}
	}

	return vecVars;
}


template<class t_real>
std::string FileFrm<t_real>::GetCountVar() const
{
	std::string strRet;

	// find counters that have counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX((det[a-z]*[0-9])|(ctr[0-9])|(counter[0-9])|([a-z0-9\.]*roi))REX", strRet, true))
		return strRet;

	// also include counters that have no counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX((det[a-z]*[0-9])|(ctr[0-9])|(counter[0-9])|([a-z0-9\.]*roi))REX", strRet, true, false))
		return strRet;

	return "";
}


template<class t_real>
std::string FileFrm<t_real>::GetMonVar() const
{
	std::string strRet;

	// find counters that have counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX((mon[a-z]*[0-9]))REX", strRet, true))
		return strRet;

	// also include counters that have no counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX((mon[a-z]*[0-9]))REX", strRet, true, false))
		return strRet;

	return "";
}


template<class t_real>
std::string FileFrm<t_real>::GetScanCommand() const
{
	std::string strCmd;
	typename t_mapParams::const_iterator iter = m_mapParams.find("info");
	if(iter != m_mapParams.end())
		strCmd = iter->second;
	return strCmd;
}


template<class t_real>
std::string FileFrm<t_real>::GetTimestamp() const
{
	std::string strDate;
	typename t_mapParams::const_iterator iter = m_mapParams.find("file_timestamp");
	if(iter != m_mapParams.end())
		strDate = iter->second;
	return strDate;
}

}

#endif
