/**
 * Loads instrument-specific data files
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-2015 - 2026
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS_LOADINST_FRMOLD_IMPL_H__
#define __TLIBS_LOADINST_FRMOLD_IMPL_H__

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


/**
 * loader for the frm text file format
 */
template<class t_real>
bool FileFrmOld<t_real>::ReadHeader(std::istream& istr)
{
	std::string section;

	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);

		if(strLine.length() == 0)
			continue;
		if(strLine.length() > 0 && strLine[0] == '*')
			continue;

		std::pair<std::string, std::string> pairLine =
			split_first<std::string>(strLine, ":", true);

		if(pairLine.first == "")
			continue;

		pairLine.first = str_to_lower(pairLine.first);

		// start of a new section in the file
		if(pairLine.second == "")
		{
			section = pairLine.first;

			if(section == "scan data")
				return true;
			continue;
		}

		// add file section name to field
		if(section != "")
			pairLine.first = section + " / " + pairLine.first;

		typename t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

		if(iter == m_mapParams.end())
			m_mapParams.insert(pairLine);
		else
			iter->second += ", " + pairLine.second;
	}

	return false;
}


template<class t_real>
void FileFrmOld<t_real>::ReadData(std::istream& istr)
{
	//skip_after_line<char>(istr, "scan data:", true, false);

	// scan command
	std::getline(istr, m_command);

	// column headers
	std::string strLineQuantities;
	std::getline(istr, strLineQuantities);
	get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t;", m_vecQuantities);

	std::string strLineUnits;
	std::getline(istr, strLineUnits);
	get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t;", m_vecUnits);


	m_vecData.resize(m_vecQuantities.size());

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);
		if(strLine.length() == 0 || strLine[0] == '#')
			continue;
		if(strLine.length() > 0 & strLine[0] == '*')
			break;

		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(strLine, " \t;", vecToks);

		if(vecToks.size() != m_vecQuantities.size())
		{
			log_warn("Loader: Line size mismatch in line: \"", strLine, "\".");

			// add zeros
			while(m_vecQuantities.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok = 0; iTok < vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}

	FileInstrBase<t_real>::RenameDuplicateCols();
}


template<class t_real>
bool FileFrmOld<t_real>::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

#if !defined NO_IOSTR
	std::shared_ptr<std::istream> ptrIstr = create_autodecomp_istream(ifstr);
	if(!ptrIstr)
		return false;
	std::istream *pIstr = ptrIstr.get();
#else
	std::istream *pIstr = &ifstr;
#endif

	if(!ReadHeader(*pIstr))
	{
		log_warn("Loader: No scan data found.");
		return false;
	}
	ReadData(*pIstr);

	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileFrmOld<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileFrmOld*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileFrmOld<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	for(std::size_t i = 0; i < m_vecQuantities.size(); ++i)
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
std::array<t_real, 3> FileFrmOld<t_real>::GetSampleLattice() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("sample information / a,b,c (a)");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second, false);
	if(vec.size() != 3)
	{
		log_err("Invalid lattice array size.");
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}

	return std::array<t_real,3>{{vec[0],vec[1],vec[2]}};
}


template<class t_real>
std::array<t_real, 3> FileFrmOld<t_real>::GetSampleAngles() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("sample information / alpha,beta,gamma (deg)");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second, false);
	if(vec.size() != 3)
	{
		log_err("Invalid angle array size.");
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}

	return std::array<t_real,3>{{d2r(vec[0]), d2r(vec[1]), d2r(vec[2])}};
}


template<class t_real>
std::array<t_real, 2> FileFrmOld<t_real>::GetMonoAnaD() const
{
	typename t_mapParams::const_iterator iterM = m_mapParams.find("monochromator / d- value (A)");
	typename t_mapParams::const_iterator iterA = m_mapParams.find("analyser / d- value (A)");

	t_real m = (iterM != m_mapParams.end() ? str_to_var<t_real>(iterM->second) : 3.355);
	t_real a = (iterA != m_mapParams.end() ? str_to_var<t_real>(iterA->second) : 3.355);

	return std::array<t_real,2>{{ m, a }};
}


template<class t_real>
std::array<bool, 3> FileFrmOld<t_real>::GetScatterSenses() const
{
	std::vector<int> vec;

	typename t_mapParams::const_iterator iter;
	for(iter = m_mapParams.begin(); iter != m_mapParams.end(); ++iter)
	{
		if(iter->first.find("scattering sense") != std::string::npos)
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

	return std::array<bool,3>{{ vec[0] > 0, vec[1] > 0, vec[2] > 0 }};
}


template<class t_real>
std::array<t_real, 3> FileFrmOld<t_real>::GetScatterPlane0() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("sample information / 1st orientation reflection");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second, true, true);
	if(vec.size() != 3)
	{
		log_err("Invalid sample peak 1 array size.");
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}
	return std::array<t_real,3>{{ vec[0], vec[1], vec[2] }};
}


template<class t_real>
std::array<t_real, 3> FileFrmOld<t_real>::GetScatterPlane1() const
{
	typename t_mapParams::const_iterator iter = m_mapParams.find("sample information / 2nd orientation reflection");
	if(iter == m_mapParams.end())
		return std::array<t_real,3>{{ 0., 0., 0. }};

	std::vector<t_real> vec = get_py_array<std::string, std::vector<t_real>>(iter->second, true, true);
	if(vec.size() != 3)
	{
		log_err("Invalid sample peak 2 array size.");
		return std::array<t_real,3>{{ 0., 0., 0. }};
	}
	return std::array<t_real,3>{{ -vec[0], -vec[1], -vec[2] }};  // LH -> RH
}


template<class t_real>
std::array<t_real, 4> FileFrmOld<t_real>::GetPosHKLE() const
{
	const t_vecVals& hs = GetCol("h");
	const t_vecVals& ks = GetCol("k");
	const t_vecVals& ls = GetCol("l");
	const t_vecVals& Es = GetCol("E");

	// get first position
	t_real h = hs.size() ? hs[0] : 0.;
	t_real k = ks.size() ? ks[0] : 0.;
	t_real l = ls.size() ? ls[0] : 0.;
	t_real E = Es.size() ? Es[0] : 0.;

	return std::array<t_real,4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileFrmOld<t_real>::GetKFix() const
{
	if(IsKiFixed())
	{
		const t_vecVals& kis = GetCol("mono");
		return tl::mean_value(kis);
	}
	else
	{
		const t_vecVals& kfs = GetCol("ana");
		return tl::mean_value(kfs);
	}
}


template<class t_real>
bool FileFrmOld<t_real>::IsKiFixed() const
{
	const t_vecVals& kis = GetCol("mono");
	const t_vecVals& kfs = GetCol("ana");
	if(kis.size() == 0 || kfs.size() == 0)
		return false;  // assume kf = const.

	t_real dki = tl::std_dev(kis);
	t_real dkf = tl::std_dev(kfs);

	return dki < dkf;
}


template<class t_real>
std::size_t FileFrmOld<t_real>::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}


template<class t_real>
std::array<t_real, 5> FileFrmOld<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstrBase<t_real>::GetScanHKLKiKf("h", "k", "l", "E", i);
}


template<class t_real>
bool FileFrmOld<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	if(!FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch))
		return false;

	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan filename
		typename t_mapParams::iterator iter = m_mapParams.find("filename");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetTitle() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("title");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetProposal() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("proposal");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetUser() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("user");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetInstrument() const
{
	std::string strInst;
	typename t_mapParams::const_iterator iterInst = m_mapParams.find("installation");
	if(iterInst != m_mapParams.end())
		strInst = iterInst->second;

	std::string strInstr;
	typename t_mapParams::const_iterator iterInstr = m_mapParams.find("instrument");
	if(iterInstr != m_mapParams.end())
		strInstr = iterInstr->second;

	if(strInst != "")
		return strInst + " / " + strInstr;
	return strInstr;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetLocalContact() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("responsable");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetScanNumber() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("filename");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetSampleName() const
{
	std::string strName;
	typename t_mapParams::const_iterator iter = m_mapParams.find("name");
	if(iter != m_mapParams.end())
		strName = iter->second;
	return strName;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetSpacegroup() const
{
	std::string strSG;
	typename t_mapParams::const_iterator iter = m_mapParams.find("sample information / spacegroup");
	if(iter != m_mapParams.end())
		strSG = iter->second;
	return strSG;
}


template<class t_real>
std::vector<std::string> FileFrmOld<t_real>::GetScannedVars() const
{
	std::vector<std::string> vecVars;

	// try qscan/qcscan command
	const std::string strRegex = R"REX((qscan|qcscan)\((.*)\))REX";
	std::regex rx(strRegex, std::regex::ECMAScript|std::regex_constants::icase);
	std::smatch m;
	if(std::regex_search(m_command, m, rx) && m.size() > 2)
	{
		const std::string& strSteps = m[2];
		std::vector<t_real> vecSteps = get_py_array<std::string, std::vector<t_real>>(strSteps, false);

		if(vecSteps.size() > 4 && !float_equal<t_real>(vecSteps[4], 0.))
			vecVars.push_back("h");
		if(vecSteps.size() > 5 && !float_equal<t_real>(vecSteps[5], 0.))
			vecVars.push_back("k");
		if(vecSteps.size() > 6 && !float_equal<t_real>(vecSteps[6], 0.))
			vecVars.push_back("l");
		if(vecSteps.size() > 7 && !float_equal<t_real>(vecSteps[7], 0.))
			vecVars.push_back("E");
	}

	if(vecVars.size() == 0)
	{
		// try scan/cscan
		const std::string strRegexDevScan = R"REX((scan|cscan)\(([a-z0-9_\.]+)[, ]+.*\))REX";
		std::regex rxDev(strRegexDevScan, std::regex::ECMAScript|std::regex_constants::icase);
		std::smatch mDev;
		if(std::regex_search(m_command, mDev, rxDev) && mDev.size() > 2)
		{
			std::string strDev = mDev[2];
			if(std::find(m_vecQuantities.begin(), m_vecQuantities.end(), strDev) != m_vecQuantities.end())
				vecVars.push_back(strDev);
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
std::string FileFrmOld<t_real>::GetCountVar() const
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
std::string FileFrmOld<t_real>::GetTimerVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(time)REX", strRet, true))
		return strRet;
	return "";
}


template<class t_real>
std::string FileFrmOld<t_real>::GetMonVar() const
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
std::string FileFrmOld<t_real>::GetScanCommand() const
{
	return m_command;
}


template<class t_real>
std::string FileFrmOld<t_real>::GetTimestamp() const
{
	std::string strDate;
	typename t_mapParams::const_iterator iter = m_mapParams.find("created at");
	if(iter != m_mapParams.end())
		strDate = iter->second;
	return strDate;
}

}

#endif
