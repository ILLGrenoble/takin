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

#ifndef __TLIBS_LOADINST_RAW_IMPL_H__
#define __TLIBS_LOADINST_RAW_IMPL_H__

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


/**
 * loader for a simple column-based text format
 */
template<class t_real>
bool FileRaw<t_real>::Load(const char* pcFile)
{
	bool bOk = m_dat.Load(pcFile);
	m_vecCols.clear();
	for(std::size_t iCol=0; iCol<m_dat.GetColumnCount(); ++iCol)
		m_vecCols.emplace_back(var_to_str(iCol+1));
	return bOk;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileRaw<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileRaw*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileRaw<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	std::size_t iCol = str_to_var<std::size_t>(strName)-1;
	if(iCol < m_dat.GetColumnCount())
	{
		if(pIdx) *pIdx = iCol;
		return m_dat.GetColumn(iCol);
	}

	static std::vector<t_real> vecNull;
	if(pIdx)
		*pIdx = m_dat.GetColumnCount();

	log_err("Column \"", strName, "\" does not exist.");
	return vecNull;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecDat&
FileRaw<t_real>::GetData() const
{
	return m_dat.GetData();
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecDat&
FileRaw<t_real>::GetData()
{
	return m_dat.GetData();
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecColNames&
FileRaw<t_real>::GetColNames() const
{
	return m_vecCols;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_mapParams&
FileRaw<t_real>::GetAllParams() const
{
	return m_dat.GetHeader();
}


template<class t_real>
std::array<t_real, 3> FileRaw<t_real>::GetSampleLattice() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	{
		typename t_map::const_iterator iter = params.find("sample_a");
		if(iter != params.end())
			a = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sample_b");
		if(iter != params.end())
			b = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sample_c");
		if(iter != params.end())
			c = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,3>{{a, b, c}};
}


template<class t_real>
std::array<t_real, 3> FileRaw<t_real>::GetSampleAngles() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	{
		typename t_map::const_iterator iter = params.find("sample_alpha");
		if(iter != params.end())
			a = d2r(str_to_var<t_real>(iter->second));
	}
	{
		typename t_map::const_iterator iter = params.find("sample_beta");
		if(iter != params.end())
			b = d2r(str_to_var<t_real>(iter->second));
	}
	{
		typename t_map::const_iterator iter = params.find("sample_gamma");
		if(iter != params.end())
			c = d2r(str_to_var<t_real>(iter->second));
	}

	return std::array<t_real,3>{{a, b, c}};
}


template<class t_real>
std::array<t_real, 2> FileRaw<t_real>::GetMonoAnaD() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real m{0}, a{0};

	{
		typename t_map::const_iterator iter = params.find("mono_d");
		if(iter != params.end())
			m = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("ana_d");
		if(iter != params.end())
			a = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FileRaw<t_real>::GetScatterSenses() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real m{0}, s{1}, a{0};

	{
		typename t_map::const_iterator iter = params.find("sense_m");
		if(iter != params.end())
			m = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sense_s");
		if(iter != params.end())
			s = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("sense_a");
		if(iter != params.end())
			a = str_to_var<t_real>(iter->second);
	}

	return std::array<bool,3>{{m>0., s>0., a>0.}};
}


template<class t_real>
std::array<t_real, 3> FileRaw<t_real>::GetScatterPlane0() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real x{0}, y{0}, z{0};

	{
		typename t_map::const_iterator iter = params.find("orient1_x");
		if(iter != params.end())
			x = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient1_y");
		if(iter != params.end())
			y = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient1_z");
		if(iter != params.end())
			z = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,3>{{x, y, z}};
}


template<class t_real>
std::array<t_real, 3> FileRaw<t_real>::GetScatterPlane1() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real x{0}, y{0}, z{0};

	{
		typename t_map::const_iterator iter = params.find("orient2_x");
		if(iter != params.end())
			x = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient2_y");
		if(iter != params.end())
			y = str_to_var<t_real>(iter->second);
	}
	{
		typename t_map::const_iterator iter = params.find("orient2_z");
		if(iter != params.end())
			z = str_to_var<t_real>(iter->second);
	}

	return std::array<t_real,3>{{x, y, z}};
}


template<class t_real>
std::string FileRaw<t_real>::GetColNameFromParam(
	const std::string& paramName, const std::string& defaultVal) const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find(paramName);
	if(iter != params.end())
		return iter->second;

	return defaultVal;
}


template<class t_real>
std::array<t_real, 4> FileRaw<t_real>::GetPosHKLE() const
{
	// get column names
	std::string strColH = GetColNameFromParam("col_h", "1");
	std::string strColK = GetColNameFromParam("col_k", "2");
	std::string strColL = GetColNameFromParam("col_l", "3");
	std::string strColE = GetColNameFromParam("col_E", "4");

	// get the first position from the scan if available
	const t_vecVals& vecH = GetCol(strColH);
	const t_vecVals& vecK = GetCol(strColK);
	const t_vecVals& vecL = GetCol(strColL);
	const t_vecVals& vecE = GetCol(strColE);

	// get values
	t_real h = 0;
	t_real k = 0;
	t_real l = 0;
	t_real E = 0;

	if(vecH.size()) h = vecH[0];
	if(vecK.size()) k = vecK[0];
	if(vecL.size()) l = vecL[0];
	if(vecE.size()) E = vecE[0];

	return std::array<t_real, 4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileRaw<t_real>::GetKFix() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real k{0};

	typename t_map::const_iterator iter = params.find("k_fix");
	if(iter != params.end())
		k = str_to_var<t_real>(iter->second);

	return k;
}


template<class t_real>
bool FileRaw<t_real>::IsKiFixed() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	bool b{0};

	typename t_map::const_iterator iter = params.find("is_ki_fixed");
	if(iter != params.end())
		b = (str_to_var<int>(iter->second) != 0);

	return b;
}


template<class t_real>
std::size_t FileRaw<t_real>::GetScanCount() const
{
	if(m_dat.GetColumnCount() != 0)
		return m_dat.GetRowCount();
	return 0;
}


template<class t_real>
std::array<t_real, 5> FileRaw<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	// get column names
	std::string strColH = GetColNameFromParam("col_h", "1");
	std::string strColK = GetColNameFromParam("col_k", "2");
	std::string strColL = GetColNameFromParam("col_l", "3");
	std::string strColE = GetColNameFromParam("col_E", "4");

	return FileInstrBase<t_real>::GetScanHKLKiKf(
		strColH.c_str(), strColK.c_str(), strColL.c_str(),
		strColE.c_str(), i);
}


template<class t_real> std::vector<std::string> FileRaw<t_real>::GetScannedVars() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColVars;

	{
		typename t_map::const_iterator iter = params.find("cols_scanned");
		if(iter != params.end())
			strColVars = iter->second;
	}

	std::vector<std::string> vecVars;
	get_tokens<std::string, std::string>(strColVars, ",;", vecVars);

	// if nothing is given, default to E
	if(!vecVars.size())
		vecVars.push_back("4");

	return vecVars;
}


template<class t_real> std::string FileRaw<t_real>::GetCountVar() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr = "5";

	{
		typename t_map::const_iterator iter = params.find("col_ctr");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real> std::string FileRaw<t_real>::GetMonVar() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr = "6";

	{
		typename t_map::const_iterator iter = params.find("col_mon");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real> std::string FileRaw<t_real>::GetCountErr() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr;

	{
		typename t_map::const_iterator iter = params.find("col_ctr_err");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real> std::string FileRaw<t_real>::GetMonErr() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr;

	{
		typename t_map::const_iterator iter = params.find("col_mon_err");
		if(iter != params.end())
			strColCtr = iter->second;
	}

	return strColCtr;
}


template<class t_real>
bool FileRaw<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	return FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch);
}


template<class t_real> std::string FileRaw<t_real>::GetTitle() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetProposal() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetUser() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetInstrument() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetLocalContact() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetScanNumber() const { return "0"; }
template<class t_real> std::string FileRaw<t_real>::GetSampleName() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetSpacegroup() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetScanCommand() const { return ""; }
template<class t_real> std::string FileRaw<t_real>::GetTimestamp() const { return ""; }


}


#endif
