/**
 * loads instrument-specific data files
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

#ifndef __TLIBS_LOADINST_PSI_IMPL_H__
#define __TLIBS_LOADINST_PSI_IMPL_H__

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
std::string FilePsi<t_real>::ReadData(std::istream& istr)
{
	std::size_t iLine = 0;

	// header
	std::string strHdr;
	std::getline(istr, strHdr);
	trim(strHdr);
	++iLine;

	get_tokens<std::string, std::string, t_vecColNames>(strHdr, " \t", m_vecColNames);
	for(std::string& _str : m_vecColNames)
	{
		find_all_and_replace<std::string>(_str, "\n", "");
		find_all_and_replace<std::string>(_str, "\r", "");
	}

	m_vecData.resize(m_vecColNames.size());
	FileInstrBase<t_real>::RenameDuplicateCols();


	// data
	while(!istr.eof())
	{
		// seekg/tellg doesn't work with compressed iostreams, return line instead!
		//std::streampos pos = istr.tellg();
		std::string line;
		std::getline(istr, line);

		// begin of another section (and end of data section)?
		if(line.find(':') != std::string::npos)
		{
			//istr.seekg(pos);
			return line;
		}

		trim(line);
		++iLine;

		if(line.length() == 0 || line[0] == '#')
			continue;

		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(line, " \t", vecToks);

		if(vecToks.size() != m_vecColNames.size())
		{
			log_warn("Loader: Column size mismatch in data line ", iLine,
				": Expected ", m_vecColNames.size(), ", got ", vecToks.size(), ".");

			// add zeros
			while(m_vecColNames.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok = 0; iTok < vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}

	// no leftover line
	return "";
}


template<class t_real>
std::string FilePsi<t_real>::ReadMultiData(std::istream& istr)
{
	// skip over multi-analyser data for the moment

	while(!istr.eof())
	{
		// seekg/tellg doesn't work with compressed iostreams, return line instead!
		//std::streampos pos = istr.tellg();
		std::string line;
		std::getline(istr, line);

		// begin of another section (and end of data section)?
		if(line.find(':') != std::string::npos)
		{
			//istr.seekg(pos);
			return line;
		}
	}

	// no leftover line
	return "";
}


template<class t_real>
void FilePsi<t_real>::GetInternalParams(const std::string& strAll,
	FilePsi<t_real>::t_mapIParams& mapPara, bool fix_broken)
{
	std::vector<std::string> vecToks;
	get_tokens<std::string, std::string>(strAll, ",\n", vecToks);

	for(const std::string& strTok : vecToks)
	{
		if(fix_broken)
		{
			// ignore broken parameter strings
			if(std::count(strTok.begin(), strTok.end(), '=') > 1)
				continue;
		}

		std::pair<std::string, std::string> pair =
			split_first<std::string>(strTok, "=", 1);

		if(pair.first == "")
			continue;

		t_real dVal = str_to_var<t_real>(pair.second);
		mapPara.emplace(typename t_mapIParams::value_type(pair.first, dVal));
	}

	// sometimes the "steps" parameters are written without commas
	if(fix_broken)
	{
		const std::string strRegex = R"REX(([a-zA-Z0-9]+)[ \t]*\=[ \t]*([+\-]?[0-9\.]+))REX";
		std::regex rx(strRegex, std::regex::ECMAScript|std::regex_constants::icase);

		std::smatch m;
		std::string data = strAll;
		while(std::regex_search(data, m, rx))
		{
			if(m.size() >= 3)
			{
				const std::string& strkey = m[1];
				const std::string& strval = m[2];
				if(strkey == "")
					continue;

				t_real val = str_to_var<t_real>(strval);
				mapPara.emplace(typename t_mapIParams::value_type(strkey, val));
			}

			data = m.suffix();
		}
	}
}


/**
 * check if a string describes a number
 */
static inline bool is_num(const std::string& str)
{
	for(typename std::string::value_type c : str)
	{
		if(c!='0' && c!='1' && c!='2' && c!='3' && c!='4'
			&& c!='5' && c!='6' && c!='7' && c!='8' && c!='9'
			&& c!='+'&& c!='-'&& c!='.')
			return false;
	}

	return true;
};


/**
 * parse the incoming and outgoing polarisation states from a .pal script file
 */
template<typename t_real = double>
static std::vector<std::array<t_real, 6>> parse_pol_states(const std::string& polScript,
	const std::string& strPolVec1, const std::string& strPolVec2,
	const std::string& strPolCur1, const std::string& strPolCur2,
	const std::string& strXYZ,
	const std::string& strFlip1, const std::string& strFlip2)
{
	std::vector<std::array<t_real, 6>> vecPolStates;

	std::vector<std::string> vecLines;
	get_tokens<std::string, std::string>(polScript, ",", vecLines);

	// initial and final polarisation states
	t_real Pix = t_real(0), Piy = t_real(0), Piz = t_real(0);
	t_real Pfx = t_real(0), Pfy = t_real(0), Pfz = t_real(0);
	t_real Pi_sign = t_real(1);
	t_real Pf_sign = t_real(1);

	bool bIsSphericalPA = true;
	bool bSwitchOn = false;

	// iterate command lines
	for(std::string& strLine : vecLines)
	{
		trim(strLine);
		strLine = str_to_lower(strLine);

		std::vector<std::string> vecLine;
		get_tokens<std::string, std::string>(strLine, " \t", vecLine);

		if(vecLine.size() == 0)
			continue;

		// first token: dr command
		// polarisation vector or current driven
		if(vecLine[0] == "dr")
		{
			std::string strCurDev = "";
			std::size_t iCurComp = 0;

			// scan next tokens in drive command
			for(std::size_t iDr = 1; iDr < vecLine.size(); ++iDr)
			{
				const std::string& strWord = vecLine[iDr];

				if(is_num(strWord))	// value to drive to
				{
					t_real dNum = str_to_var<t_real>(strWord);

					// --------------------------------------------------
					// for spherical polarisation analysis
					// --------------------------------------------------
					if(strCurDev == strPolVec1)
					{
						// incoming polarisation vector changed
						switch(iCurComp)
						{
							case 0: Pix = dNum; break;
							case 1: Piy = dNum; break;
							case 2: Piz = dNum; break;
						}
					}
					else if(strCurDev == strPolVec2)
					{	// outgoing polarisation vector changed
						switch(iCurComp)
						{
							case 0: Pfx = dNum; break;
							case 1: Pfy = dNum; break;
							case 2: Pfz = dNum; break;
						}
					}
					else if(strCurDev == strPolCur1)
					{	// sign of polarisation vector 1 changed
						if(iCurComp == 0)
						{
							if(dNum >= t_real(0))
								Pi_sign = t_real(1);
							else
								Pi_sign = t_real(-1);
						}
					}
					else if(strCurDev == strPolCur2)
					{	// sign of polarisation vector 2 changed
						if(iCurComp == 0)
						{
							if(dNum >= t_real(0))
								Pf_sign = t_real(1);
							else
								Pf_sign = t_real(-1);
						}
					}
					// --------------------------------------------------

					// --------------------------------------------------
					// for linear polarisation analysis
					// --------------------------------------------------
					else if(strCurDev == strXYZ)
					{
						bIsSphericalPA = 0;

						// sample guide field vector changed
						switch(iCurComp)
						{
							case 0: Pfx = Pix = dNum; break;
							case 1: Pfy = Piy = dNum; break;
							case 2: Pfz = Piz = dNum; break;
						}
					}
					// --------------------------------------------------

					++iCurComp;
				}
				else	// (next) device to drive
				{
					iCurComp = 0;
					strCurDev = strWord;
				}
			}
		}
		// --------------------------------------------------
		// for linear polarisation analysis
		// --------------------------------------------------
		else if(vecLine.size()>1 && (vecLine[0] == "on" || vecLine[0] == "off" || vecLine[0] == "of"))
		{
			if(vecLine[0] == "on")
				bSwitchOn = 1;
			else if(vecLine[0] == "off" || vecLine[0] == "of")
				bSwitchOn = 0;

			std::string strFlip = vecLine[1];

			if(strFlip == strFlip1)
			{
				bIsSphericalPA = 0;
				Pi_sign = bSwitchOn ? t_real(-1) : t_real(1);
			}
			else if(strFlip == strFlip2)
			{
				bIsSphericalPA = 0;
				Pf_sign = bSwitchOn ? t_real(-1) : t_real(1);
			}
		}
		// --------------------------------------------------


		// count command issued -> save current spin states
		else if(vecLine[0] == "co")
		{
			// for linear PA, we only have principal directions
			// and diagonal elements of the P matrix
			// values in the other components are correction currents
			if(!bIsSphericalPA)
			{
				int maxComp = 0;
				if(std::abs(Piy) > std::abs(Pix) && std::abs(Piy) > std::abs(Piz))
					maxComp = 1;
				if(std::abs(Piz) > std::abs(Piy) && std::abs(Piz) > std::abs(Pix))
					maxComp = 2;

				Pix = Piy = Piz = Pfx = Pfy = Pfz = 0;
				if(maxComp == 0)
					Pix = Pfx = 1;
				else if(maxComp == 1)
					Piy = Pfy = 1;
				else if(maxComp == 2)
					Piz = Pfz = 1;
			}


			vecPolStates.push_back(std::array<t_real,6>({{
				Pi_sign*Pix, Pi_sign*Piy, Pi_sign*Piz,
				Pf_sign*Pfx, Pf_sign*Pfy, Pf_sign*Pfz }}));
		}
	}

	// cleanup
	for(std::size_t iPol = 0; iPol < vecPolStates.size(); ++iPol)
	{
		for(unsigned iComp = 0; iComp < 6; ++iComp)
			set_eps_0(vecPolStates[iPol][iComp]);
	}

	return vecPolStates;
}


template<class t_real>
void FilePsi<t_real>::ParsePolData()
{
	m_vecPolStates.clear();
	typename t_mapParams::const_iterator iter = m_mapParams.find("POLAN");
	if(iter == m_mapParams.end())
		return;

	m_vecPolStates = parse_pol_states<t_real>(iter->second,
		m_strPolVec1, m_strPolVec2,
		m_strPolCur1, m_strPolCur2,
		m_strXYZ,
		m_strFlip1, m_strFlip2);
}


template<class t_real>
void FilePsi<t_real>::SetAutoParsePolData(bool b)
{
	m_bAutoParsePol = b;
}


template<class t_real>
const std::vector<std::array<t_real, 6>>& FilePsi<t_real>::GetPolStates() const
{
	return m_vecPolStates;
}


template<class t_real>
void FilePsi<t_real>::SetPolNames(const char* pVec1, const char* pVec2,
	const char* pCur1, const char* pCur2)
{
	m_strPolVec1 = pVec1; m_strPolVec2 = pVec2;
	m_strPolCur1 = pCur1; m_strPolCur2 = pCur2;
}


template<class t_real>
void FilePsi<t_real>::SetLinPolNames(const char* pFlip1, const char* pFlip2,
	const char* pXYZ)
{
	m_strFlip1 = pFlip1; m_strFlip2 = pFlip2;
	m_strXYZ = pXYZ;
}


template<class t_real>
bool FilePsi<t_real>::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

#if !defined NO_IOSTR
	std::shared_ptr<std::istream> ptrIstr = create_autodecomp_istream(ifstr);
	if(!ptrIstr) return false;
	std::istream* pIstr = ptrIstr.get();
#else
	std::istream* pIstr = &ifstr;
#endif

	std::string nextLine;
	while(!pIstr->eof())
	{
		std::string strLine;

		// see if a line is left over from the previous loop
		if(nextLine == "")
		{
			// read a new line
			std::getline(*pIstr, strLine);
		}
		else
		{
			strLine = nextLine;
			nextLine = "";
		}

		if(strLine.substr(0,4) == "RRRR")
			skip_after_line<char>(*pIstr, "VVVV", true);

		std::pair<std::string, std::string> pairLine =
			split_first<std::string>(strLine, ":", 1);
		if(pairLine.first == "DATA_")
			nextLine = ReadData(*pIstr);
		else if(pairLine.first == "MULTI")
			nextLine = ReadMultiData(*pIstr);
		else if(pairLine.first == "")
			continue;
		else
		{
			typename t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
	}

	typename t_mapParams::const_iterator
		iterParams = m_mapParams.find("PARAM"),
		iterZeros = m_mapParams.find("ZEROS"),
		iterVars = m_mapParams.find("VARIA"),
		iterPos = m_mapParams.find("POSQE"),
		iterSteps = m_mapParams.find("STEPS");

	if(iterParams!=m_mapParams.end())
		GetInternalParams(iterParams->second, m_mapParameters);
	if(iterZeros!=m_mapParams.end())
		GetInternalParams(iterZeros->second, m_mapZeros);
	if(iterVars!=m_mapParams.end())
		GetInternalParams(iterVars->second, m_mapVariables);
	if(iterPos!=m_mapParams.end())
		GetInternalParams(iterPos->second, m_mapPosHkl);
	if(iterSteps!=m_mapParams.end())
		GetInternalParams(iterSteps->second, m_mapScanSteps, true);

	if(m_bAutoParsePol)
		ParsePolData();
	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FilePsi<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FilePsi*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FilePsi<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	for(std::size_t i = 0; i < m_vecColNames.size(); ++i)
	{
		if(str_to_lower(m_vecColNames[i]) == str_to_lower(strName))
		{
			if(pIdx)
				*pIdx = i;
			return m_vecData[i];
		}
	}

	// set invalid index if the column was not found
	if(pIdx)
		*pIdx = m_vecColNames.size();

	log_err("Column \"", strName, "\" does not exist.");
	return vecNull;
}


template<class t_real>
void FilePsi<t_real>::PrintParams(std::ostream& ostr) const
{
	for(const typename t_mapParams::value_type& val : m_mapParams)
	{
		ostr << "Param: " << val.first
			<< ", Val: " << val.second << "\n";
	}
}


template<class t_real>
std::array<t_real,3> FilePsi<t_real>::GetSampleLattice() const
{
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("AS");
	typename t_mapIParams::const_iterator iterB = m_mapParameters.find("BS");
	typename t_mapIParams::const_iterator iterC = m_mapParameters.find("CS");

	t_real a = (iterA!=m_mapParameters.end() ? iterA->second : 0.);
	t_real b = (iterB!=m_mapParameters.end() ? iterB->second : 0.);
	t_real c = (iterC!=m_mapParameters.end() ? iterC->second : 0.);

	return std::array<t_real,3>{{a,b,c}};
}


template<class t_real>
std::array<t_real,3> FilePsi<t_real>::GetSampleAngles() const
{
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("AA");
	typename t_mapIParams::const_iterator iterB = m_mapParameters.find("BB");
	typename t_mapIParams::const_iterator iterC = m_mapParameters.find("CC");

	t_real alpha = (iterA!=m_mapParameters.end() ? d2r(iterA->second) : get_pi<t_real>()/2.);
	t_real beta = (iterB!=m_mapParameters.end() ? d2r(iterB->second) : get_pi<t_real>()/2.);
	t_real gamma = (iterC!=m_mapParameters.end() ? d2r(iterC->second) : get_pi<t_real>()/2.);

	return std::array<t_real,3>{{alpha, beta, gamma}};
}


template<class t_real>
std::array<t_real,2> FilePsi<t_real>::GetMonoAnaD() const
{
	typename t_mapIParams::const_iterator iterM = m_mapParameters.find("DM");
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("DA");

	t_real m = (iterM!=m_mapParameters.end() ? iterM->second : 3.355);
	t_real a = (iterA!=m_mapParameters.end() ? iterA->second : 3.355);

	return std::array<t_real,2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FilePsi<t_real>::GetScatterSenses() const
{
	typename t_mapIParams::const_iterator iterM = m_mapParameters.find("SM");
	typename t_mapIParams::const_iterator iterS = m_mapParameters.find("SS");
	typename t_mapIParams::const_iterator iterA = m_mapParameters.find("SA");

	bool m = (iterM!=m_mapParameters.end() ? (iterM->second>0) : 0);
	bool s = (iterS!=m_mapParameters.end() ? (iterS->second>0) : 1);
	bool a = (iterA!=m_mapParameters.end() ? (iterA->second>0) : 0);

	return std::array<bool,3>{{m, s, a}};
}


template<class t_real>
std::array<t_real, 3> FilePsi<t_real>::GetScatterPlane0() const
{
	typename t_mapIParams::const_iterator iterX = m_mapParameters.find("AX");
	typename t_mapIParams::const_iterator iterY = m_mapParameters.find("AY");
	typename t_mapIParams::const_iterator iterZ = m_mapParameters.find("AZ");

	t_real x = (iterX!=m_mapParameters.end() ? iterX->second : 1.);
	t_real y = (iterY!=m_mapParameters.end() ? iterY->second : 0.);
	t_real z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<t_real,3>{{x,y,z}};
}


template<class t_real>
std::array<t_real, 3> FilePsi<t_real>::GetScatterPlane1() const
{
	typename t_mapIParams::const_iterator iterX = m_mapParameters.find("BX");
	typename t_mapIParams::const_iterator iterY = m_mapParameters.find("BY");
	typename t_mapIParams::const_iterator iterZ = m_mapParameters.find("BZ");

	t_real x = (iterX!=m_mapParameters.end() ? iterX->second : 0.);
	t_real y = (iterY!=m_mapParameters.end() ? iterY->second : 1.);
	t_real z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<t_real,3>{{x,y,z}};
}


template<class t_real>
t_real FilePsi<t_real>::GetKFix() const
{
	typename t_mapIParams::const_iterator iterK = m_mapParameters.find("KFIX");
	t_real k = (iterK!=m_mapParameters.end() ? iterK->second : 0.);

	return k;
}


template<class t_real>
std::array<t_real, 4> FilePsi<t_real>::GetPosHKLE() const
{
	typename t_mapIParams::const_iterator iterH = m_mapPosHkl.find("QH");
	typename t_mapIParams::const_iterator iterK = m_mapPosHkl.find("QK");
	typename t_mapIParams::const_iterator iterL = m_mapPosHkl.find("QL");
	typename t_mapIParams::const_iterator iterE = m_mapPosHkl.find("EN");

	t_real h = (iterH!=m_mapPosHkl.end() ? iterH->second : 0.);
	t_real k = (iterK!=m_mapPosHkl.end() ? iterK->second : 0.);
	t_real l = (iterL!=m_mapPosHkl.end() ? iterL->second : 0.);
	t_real E = (iterE!=m_mapPosHkl.end() ? iterE->second : 0.);

	return std::array<t_real,4>{{ h, k, l, E }};
}


template<class t_real>
std::array<t_real, 4> FilePsi<t_real>::GetDeltaHKLE() const
{
	typename t_mapIParams::const_iterator iterH = m_mapScanSteps.find("DQH");
	if(iterH==m_mapScanSteps.end()) iterH = m_mapScanSteps.find("QH");

	typename t_mapIParams::const_iterator iterK = m_mapScanSteps.find("DQK");
	if(iterK==m_mapScanSteps.end()) iterK = m_mapScanSteps.find("QK");

	typename t_mapIParams::const_iterator iterL = m_mapScanSteps.find("DQL");
	if(iterL==m_mapScanSteps.end()) iterL = m_mapScanSteps.find("QL");

	typename t_mapIParams::const_iterator iterE = m_mapScanSteps.find("DEN");
	if(iterE==m_mapScanSteps.end()) iterE = m_mapScanSteps.find("EN");

        t_real h = (iterH!=m_mapScanSteps.end() ? iterH->second : 0.);
        t_real k = (iterK!=m_mapScanSteps.end() ? iterK->second : 0.);
        t_real l = (iterL!=m_mapScanSteps.end() ? iterL->second : 0.);
        t_real E = (iterE!=m_mapScanSteps.end() ? iterE->second : 0.);

        return std::array<t_real,4>{{ h, k, l, E }};
}


template<class t_real>
bool FilePsi<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	if(!FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch))
		return false;

	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		typename t_mapParams::iterator iter = m_mapParams.find("FILE_");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;
}


template<class t_real>
bool FilePsi<t_real>::IsKiFixed() const
{
	typename t_mapIParams::const_iterator iter = m_mapParameters.find("FX");
	t_real val = (iter!=m_mapParameters.end() ? iter->second : 2.);

	return float_equal<t_real>(val, 1., 0.25);
}


template<class t_real>
std::size_t FilePsi<t_real>::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}


template<class t_real>
std::array<t_real, 5> FilePsi<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	// default column names
	const char *h = "QH";
	const char *k = "QK";
	const char *l = "QL";
	const char *E = "EN";

	// alternate column names
	if(!FileInstrBase<t_real>::HasCol("QH") && FileInstrBase<t_real>::HasCol("H"))
		h = "H";
	if(!FileInstrBase<t_real>::HasCol("QK") && FileInstrBase<t_real>::HasCol("K"))
		k = "K";
	if(!FileInstrBase<t_real>::HasCol("QL") && FileInstrBase<t_real>::HasCol("L"))
		l = "L";
	if(!FileInstrBase<t_real>::HasCol("EN") && FileInstrBase<t_real>::HasCol("E"))
		E = "E";

	return FileInstrBase<t_real>::GetScanHKLKiKf(h, k, l, E, i);
}


template<class t_real>
std::string FilePsi<t_real>::GetTitle() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("TITLE");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FilePsi<t_real>::GetProposal() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("EXPNO");
	if(iter != m_mapParams.end())
		strTitle = iter->second;

	if(strTitle == "")
	{
		typename t_mapIParams::const_iterator iter = m_mapParameters.find("ProposalID");
		if(iter != m_mapParameters.end())
			strTitle = tl::var_to_str(unsigned(iter->second));
	}

	return strTitle;
}


template<class t_real>
std::string FilePsi<t_real>::GetUser() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("USER_");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FilePsi<t_real>::GetInstrument() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("INSTR");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FilePsi<t_real>::GetLocalContact() const
{
	std::string strUser;
	typename t_mapParams::const_iterator iter = m_mapParams.find("LOCAL");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}


template<class t_real>
std::string FilePsi<t_real>::GetScanNumber() const
{
	std::string strTitle;
	typename t_mapParams::const_iterator iter = m_mapParams.find("FILE_");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}


template<class t_real>
std::string FilePsi<t_real>::GetSampleName() const
{
	return "";
}


template<class t_real>
std::string FilePsi<t_real>::GetSpacegroup() const
{
	return "";
}


template<class t_real>
std::vector<std::string> FilePsi<t_real>::GetScannedVarsFromCommand(const std::string& cmd)
{
	std::vector<std::string> vecVars;

	// get tokens
	std::vector<std::string> vecToks;
	get_tokens<std::string, std::string>(cmd, " \t", vecToks);
	for(std::string& strTok : vecToks)
		trim(strTok);
	std::transform(vecToks.begin(), vecToks.end(), vecToks.begin(), str_to_lower<std::string>);

	// try to find the deltas for a Q/E scan
	typename std::vector<std::string>::iterator iterTok
		= std::find(vecToks.begin(), vecToks.end(), "dqh");

	if(iterTok != vecToks.end())
	{
		t_real dh = str_to_var<t_real>(*(++iterTok));
		t_real dk = str_to_var<t_real>(*(++iterTok));
		t_real dl = str_to_var<t_real>(*(++iterTok));
		t_real dE = str_to_var<t_real>(*(++iterTok));

		if(!float_equal<t_real>(dh, 0.)) vecVars.push_back("QH");
		if(!float_equal<t_real>(dk, 0.)) vecVars.push_back("QK");
		if(!float_equal<t_real>(dl, 0.)) vecVars.push_back("QL");
		if(!float_equal<t_real>(dE, 0.)) vecVars.push_back("EN");
	}

	// still nothing found, try regex
	if(!vecVars.size())
	{
		const std::string strRegex = R"REX((sc|scan|bs)[ \t]+([a-z0-9]+)[ \t]+[0-9\.-]+[ \t]+[d|D]([a-z0-9]+).*)REX";
		std::regex rx(strRegex, std::regex::ECMAScript|std::regex_constants::icase);
		std::smatch m;
		if(std::regex_search(cmd, m, rx) && m.size() > 3)
		{
			const std::string& strSteps = m[3];
			vecVars.push_back(str_to_upper(strSteps));
		}
	}

	return vecVars;
}


template<class t_real>
std::vector<std::string> FilePsi<t_real>::GetScannedVars() const
{
	std::vector<std::string> vecVars;

	// steps parameter
	for(const typename t_mapIParams::value_type& pair : m_mapScanSteps)
	{
		if(!float_equal<t_real>(pair.second, 0.) && pair.first.length())
		{
			if(std::tolower(pair.first[0]) == 'd')
				vecVars.push_back(pair.first.substr(1));
			else
				vecVars.push_back(pair.first);
		}
	}

	// nothing found yet -> try scan command instead
	if(!vecVars.size())
	{
		typename t_mapParams::const_iterator iter = m_mapParams.find("COMND");
		if(iter != m_mapParams.end())
			vecVars = GetScannedVarsFromCommand(iter->second);
	}

	if(!vecVars.size())
	{
		log_warn("Could not determine scan variable.");
		if(m_vecColNames.size() >= 1)
		{
			log_warn("Using first column: \"", m_vecColNames[0], "\".");
			vecVars.push_back(m_vecColNames[0]);
		}
	}

	return vecVars;
}


template<class t_real>
std::string FilePsi<t_real>::GetCountVar() const
{
	// try to match names
	std::string strRet;

	// find counters that have counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX(cnts)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(det)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX((Single)?Detector)REX", strRet))
		return strRet;

	// also include counters that have no counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX(cnts)REX", strRet, false, false))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(det)REX", strRet, false, false))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX((Single)?Detector)REX", strRet, false, false))
		return strRet;

	return "";
}


template<class t_real>
std::string FilePsi<t_real>::GetMonVar() const
{
	std::string strRet;

	// find counters that have counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX(m[0-9]+)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(Monitor[0-9]*)REX", strRet))
		return strRet;

	// also include counters that have no counts
	if(FileInstrBase<t_real>::MatchColumn(R"REX(m[0-9]+)REX", strRet, false, false))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(Monitor[0-9]*)REX", strRet, false, false))
		return strRet;

	return "";
}


template<class t_real>
std::string FilePsi<t_real>::GetScanCommand() const
{
	std::string strCmd;
	typename t_mapParams::const_iterator iter = m_mapParams.find("COMND");
	if(iter != m_mapParams.end())
		strCmd = iter->second;
	return strCmd;
}


template<class t_real>
std::string FilePsi<t_real>::GetTimestamp() const
{
	std::string strDate;
	typename t_mapParams::const_iterator iter = m_mapParams.find("DATE_");
	if(iter != m_mapParams.end())
		strDate = iter->second;
	return strDate;
}


template<class t_real>
std::size_t FilePsi<t_real>::NumPolChannels() const
{
	return m_vecPolStates.size();
}

}


#endif
