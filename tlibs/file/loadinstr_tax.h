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

#ifndef __TLIBS_LOADINST_TAX_IMPL_H__
#define __TLIBS_LOADINST_TAX_IMPL_H__

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
bool FileTax<t_real>::Load(const char* pcFile)
{
	// load data columns
	m_dat.SetCommentChar('#');
	m_dat.SetSeparatorChars("=");
	bool bOk = m_dat.Load(pcFile);

	// get column header names
	m_vecCols.clear();

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

#if !defined NO_IOSTR
	std::shared_ptr<std::istream> ptrIstr = create_autodecomp_istream(ifstr);
	if(!ptrIstr)
		return false;
	std::istream& istr = *ptrIstr.get();
#else
	std::istream& istr = ifstr;
#endif

	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim<std::string>(strLine);
		if(strLine.length() == 0)
			continue;

		if(strLine == "# col_headers =")
		{
			// column headers in next line
			std::getline(istr, strLine);
			trim<std::string>(strLine);
			strLine = strLine.substr(1);

			get_tokens<std::string, std::string>(strLine, " \t", m_vecCols);
			break;
		}
	}

	if(m_vecCols.size() != m_dat.GetColumnCount())
	{
		log_warn("Mismatch between the number of data columns (",
			m_dat.GetColumnCount(), ") and column headers (",
			m_vecCols.size(), ").");
	}

	// fill the rest with dummy column names
	for(std::size_t iCol=m_vecCols.size(); iCol<m_dat.GetColumnCount(); ++iCol)
		m_vecCols.emplace_back(var_to_str(iCol+1));

	return bOk;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileTax<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileTax*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileTax<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	auto iter = std::find(m_vecCols.begin(), m_vecCols.end(), strName);
	if(iter == m_vecCols.end())
	{
		log_err("Column \"", strName, "\" does not exist.");
		return vecNull;
	}

	std::size_t iCol = iter - m_vecCols.begin();
	if(iCol < m_dat.GetColumnCount())
	{
		if(pIdx)
			*pIdx = iCol;
		return m_dat.GetColumn(iCol);
	}

	if(pIdx)
		*pIdx = m_dat.GetColumnCount();

	log_err("Column \"", strName, "\" does not exist.");
	return vecNull;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecDat&
FileTax<t_real>::GetData() const
{
	return m_dat.GetData();
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecDat&
FileTax<t_real>::GetData()
{
	return m_dat.GetData();
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecColNames&
FileTax<t_real>::GetColNames() const
{
	return m_vecCols;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_mapParams&
FileTax<t_real>::GetAllParams() const
{
	return m_dat.GetHeader();
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetSampleLattice() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	typename t_map::const_iterator iter = params.find("latticeconstants");
	if(iter != params.end())
	{
		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(iter->second, ",", vecToks);

		if(vecToks.size() > 0)
			a = vecToks[0];
		if(vecToks.size() > 1)
			b = vecToks[1];
		if(vecToks.size() > 2)
			c = vecToks[2];
	}

	return std::array<t_real, 3>{{a, b, c}};
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetSampleAngles() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	t_real a{0}, b{0}, c{0};

	typename t_map::const_iterator iter = params.find("latticeconstants");
	if(iter != params.end())
	{
		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(iter->second, ",", vecToks);

		if(vecToks.size() > 3)
			a = d2r(vecToks[3]);
		if(vecToks.size() > 4)
			b = d2r(vecToks[4]);
		if(vecToks.size() > 5)
			c = d2r(vecToks[5]);
	}

	return std::array<t_real, 3>{{a, b, c}};
}


template<class t_real>
std::array<t_real, 2> FileTax<t_real>::GetMonoAnaD() const
{
	t_real m{0}, a{0};

	// TODO

	return std::array<t_real, 2>{{m, a}};
}


template<class t_real>
std::array<bool, 3> FileTax<t_real>::GetScatterSenses() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	bool m{0}, s{1}, a{0};

	typename t_map::const_iterator iter = params.find("sense");
	if(iter != params.end())
	{
		if(iter->second.length() > 0)
			m = (iter->second[0] == '+');
		if(iter->second.length() > 1)
			s = (iter->second[1] == '+');
		if(iter->second.length() > 2)
			a = (iter->second[2] == '+');
	}

	return std::array<bool, 3>{{m, s, a}};
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetScatterPlaneVector(int i) const
{
	using t_mat = tl::ublas::matrix<t_real>;
	using t_vec = tl::ublas::vector<t_real>;

	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	const t_vecVals& colGl = GetCol("sgl");
	const t_vecVals& colGu = GetCol("sgu");
	t_real gl = d2r(mean_value(colGl));
	t_real gu = d2r(mean_value(colGu));

	t_real x{0}, y{0}, z{0};
	typename t_map::const_iterator iter = params.find("ubmatrix");
	if(iter != params.end())
	{
		std::vector<t_real> vecToks;
		get_tokens<t_real, std::string>(iter->second, ",", vecToks);
		t_mat UB = tl::make_mat<t_mat>({
			{ vecToks[0], vecToks[3], vecToks[6] },
			{ vecToks[1], vecToks[4], vecToks[7] },
			{ vecToks[2], vecToks[5], vecToks[8] } });

		std::array<t_real, 3> lattice = GetSampleLattice();
		std::array<t_real, 3> angles = GetSampleAngles();

		tl::Lattice<t_real> latt(
			lattice[0], lattice[1], lattice[2],
			angles[0], angles[1], angles[2]);
		latt.RotateEuler(-gl, -gu, 0.);

		t_mat B = tl::get_B(latt, true);
		t_mat B_inv;
		tl::inverse(B, B_inv);

		t_mat U = prod_mm(UB, B_inv);
		t_vec vec = tl::make_vec<t_vec>({ U(0, i), U(1, i), U(2, i) });
		vec = prod_mv(B_inv, vec);
		vec /= tl::veclen(vec);

		x = vec[0]; y = vec[1]; z = vec[2];
	}

	return std::array<t_real, 3>{{x, y, z}};
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetScatterPlane0() const
{
	return GetScatterPlaneVector(0);
}


template<class t_real>
std::array<t_real, 3> FileTax<t_real>::GetScatterPlane1() const
{
	return GetScatterPlaneVector(1);
}


template<class t_real>
std::array<t_real, 4> FileTax<t_real>::GetPosHKLE() const
{
	t_real h = 0;
	t_real k = 0;
	t_real l = 0;
	t_real E = 0;

	// get the first position from the scan if available
	const t_vecVals& vecH = GetCol("h");
	const t_vecVals& vecK = GetCol("k");
	const t_vecVals& vecL = GetCol("l");
	const t_vecVals& vecE = GetCol("e");

	if(vecH.size()) h = vecH[0];
	if(vecK.size()) k = vecK[0];
	if(vecL.size()) l = vecL[0];
	if(vecE.size()) E = vecE[0];

	return std::array<t_real, 4>{{ h, k, l, E }};
}


template<class t_real>
t_real FileTax<t_real>::GetKFix() const
{
	bool ki_fixed = IsKiFixed();
	const t_vecVals& colEfix = GetCol(ki_fixed ? "ei" : "ef");

	t_real E = mean_value(colEfix);
	bool imag;
	t_real k = E2k<units::si::system, t_real>(
		E * get_one_meV<t_real>(), imag)
			* get_one_angstrom<t_real>();
	return k;
}


template<class t_real>
bool FileTax<t_real>::IsKiFixed() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();
	bool b{0};

	typename t_map::const_iterator iter = params.find("mode");
	if(iter != params.end())
		b = (str_to_var<int>(iter->second) != 0);

	return b;
}


template<class t_real>
std::size_t FileTax<t_real>::GetScanCount() const
{
	if(m_dat.GetColumnCount() != 0)
		return m_dat.GetRowCount();
	return 0;
}


template<class t_real>
std::array<t_real, 5> FileTax<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstrBase<t_real>::GetScanHKLKiKf("h", "k", "l", "e", i);
}


template<class t_real> std::vector<std::string> FileTax<t_real>::GetScannedVars() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColVars;
	typename t_map::const_iterator iter = params.find("def_x");
	if(iter != params.end())
		strColVars = iter->second;

	std::vector<std::string> vecVars;
	get_tokens<std::string, std::string>(strColVars, ",;", vecVars);

	// if nothing is given, default to E
	if(!vecVars.size())
		vecVars.push_back("e");

	return vecVars;
}


template<class t_real> std::string FileTax<t_real>::GetCountVar() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string strColCtr = "detector";
	typename t_map::const_iterator iter = params.find("def_y");
	if(iter != params.end())
		strColCtr = iter->second;

	return strColCtr;
}


template<class t_real> std::string FileTax<t_real>::GetMonVar() const
{
	return "monitor";
}


template<class t_real>
bool FileTax<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	return FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch);
}


template<class t_real> std::string FileTax<t_real>::GetTitle() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("experiment");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetProposal() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("proposal");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetUser() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("users");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetInstrument() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("instrument");
	if(iter != params.end())
		return iter->second;
	return "tax";
}


template<class t_real> std::string FileTax<t_real>::GetLocalContact() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("local_contact");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetScanNumber() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("scan");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetSampleName() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("samplename");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetSpacegroup() const
{
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetScanCommand() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	typename t_map::const_iterator iter = params.find("command");
	if(iter != params.end())
		return iter->second;
	return "";
}


template<class t_real> std::string FileTax<t_real>::GetTimestamp() const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	const t_map& params = GetAllParams();

	std::string timestamp;

	typename t_map::const_iterator iterDate = params.find("date");
	if(iterDate != params.end())
		timestamp = iterDate->second;

	typename t_map::const_iterator iterTime = params.find("time");
	if(iterTime != params.end())
	{
		if(timestamp.length() > 0)
			timestamp += ", ";
		timestamp += iterTime->second;
	}

	return timestamp;
}

}


#endif
