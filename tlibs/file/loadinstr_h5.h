/**
 * Loads instrument-specific data files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015 - 2023
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

#ifndef __TLIBS_LOADINST_H5_IMPL_H__
#define __TLIBS_LOADINST_H5_IMPL_H__

#include "loadinstr.h"

#include "../log/log.h"
#include "../helper/py.h"
#include "../file/file.h"
#include "../math/math.h"
#include "../math/stat.h"
#include "../phys/neutrons.h"
#include "../phys/lattice.h"
#include "../file/h5.h"
#if !defined NO_IOSTR
	#include "../file/comp.h"
#endif

#include <numeric>
#include <fstream>


namespace tl{


template<class t_real>
bool FileH5<t_real>::Load(const char* pcFile)
{
	const t_real eps = 1e-6;
	const int prec = 6;

	try
	{
		H5::H5File h5file = H5::H5File(pcFile, H5F_ACC_RDONLY);

		m_data.clear();
		m_vecCols.clear();
		m_params.clear();
		m_scanned_vars.clear();

		std::vector<std::string> entries;
		if(!get_h5_entries(h5file, "/", entries) || entries.size() == 0)
		{
			log_err("No entries in hdf5 file.");
			return false;
		}

		if(entries.size() > 1)
			log_warn(entries.size(), " root entries in hdf5 file, expected a single one.");

		const std::string& entry = entries[0];

		// get data matrix
		if(!get_h5_matrix(h5file, entry + "/data_scan/scanned_variables/data", m_data))
		{
			log_err("Cannot load count data.");
			return false;
		}

		// get column names
		if(!get_h5_string_vector(h5file, entry + "/data_scan/scanned_variables/variables_names/label", m_vecCols))
		{
			log_err("Cannot load column names.");
			return false;
		}

		// get scanned variables
		std::vector<int> scanned;
		if(!get_h5_vector(h5file, entry + "/data_scan/scanned_variables/variables_names/scanned", scanned))
		{
			log_err("Cannot load scanned variables.");
			return false;
		}

		std::vector<t_real> scanned_stddevs;
		for(std::size_t idx = 0; idx<std::min(m_vecCols.size(), scanned.size()); ++idx)
		{
			const std::string& col_name = m_vecCols[idx];
			const t_vecVals& col_vec = GetCol(col_name);

			// add variable to parameter map
			t_real dStd = t_real(0);
			if(col_vec.size())
			{
				t_real dMean = mean_value(col_vec);
				dStd = std_dev(col_vec);

				std::string col_val = var_to_str(dMean, prec);
				if(!float_equal(dStd, t_real(0), eps))
					col_val += " +- " + var_to_str(dStd, prec);
				m_params.emplace(std::make_pair("var_" + col_name, col_val));
			}

			if(scanned[idx])
			{
				m_scanned_vars.push_back(col_name);
				scanned_stddevs.push_back(dStd);
			}
		}

		// add index column
		m_vecCols.insert(m_vecCols.begin(), "Point_Index");
		t_vecVals vals_idx;
		vals_idx.reserve(GetScanCount());
		for(std::size_t idx=0; idx<GetScanCount(); ++idx)
			vals_idx.push_back(idx);
		m_data.emplace(m_data.begin(), std::move(vals_idx));

		// if Q, E coordinates are among the scan variables, move them to the front
		auto iterQL = std::find(m_scanned_vars.begin(), m_scanned_vars.end(), "QL");
		if(iterQL != m_scanned_vars.end())
		{
			m_scanned_vars.erase(iterQL);
			m_scanned_vars.insert(m_scanned_vars.begin(), "QL");
		}
		auto iterQK = std::find(m_scanned_vars.begin(), m_scanned_vars.end(), "QK");
		if(iterQK != m_scanned_vars.end())
		{
			m_scanned_vars.erase(iterQK);
			m_scanned_vars.insert(m_scanned_vars.begin(), "QK");
		}
		auto iterQH = std::find(m_scanned_vars.begin(), m_scanned_vars.end(), "QH");
		if(iterQH != m_scanned_vars.end())
		{
			m_scanned_vars.erase(iterQH);
			m_scanned_vars.insert(m_scanned_vars.begin(), "QH");
		}
		auto iterEN = std::find(m_scanned_vars.begin(), m_scanned_vars.end(), "EN");
		if(iterEN != m_scanned_vars.end())
		{
			m_scanned_vars.erase(iterEN);
			m_scanned_vars.insert(m_scanned_vars.begin(), "EN");
		}

		// move the first scan variable with non-zero deviation to the front
		for(std::size_t i=0; i<m_scanned_vars.size(); ++i)
		{
			if(tl::float_equal<t_real>(scanned_stddevs[i], t_real(0), eps))
				continue;

			if(i > 0)
			{
				std::string cur_var = m_scanned_vars[i];
				m_scanned_vars.erase(m_scanned_vars.begin() + i);
				m_scanned_vars.insert(m_scanned_vars.begin(), cur_var);
				break;
			}
		}

		// get the name of the instrument if available
		bool instr_found = get_h5_string(h5file, entry + "/instrument_name", m_instrname);

		if(instr_found)
		{
			// check if the instrument group exists
			instr_found = h5file.nameExists(entry + "/" + m_instrname);
		}

		if(!instr_found)
		{
			// get the first group marked with "NXinstrument"
			std::vector<std::string> main_entries;
			get_h5_entries(h5file, entry, main_entries);
			for(const std::string& main_entry : main_entries)
			{
				std::string nx_class = get_h5_attr<std::string>(h5file, entry + "/" + main_entry, "NX_class", true);
				if(nx_class == "NXinstrument")
				{
					// found an instrument entry
					m_instrname = main_entry;
					instr_found = true;
					break;
				}
			}
		}

		if(!instr_found)
		{
			m_instrname = "instrument";
			log_err("No instrument group found, defaulting to \"", m_instrname, "\".");
		}

		// get experiment infos
		std::string timestamp_end;
		get_h5_string(h5file, entry + "/title", m_title);
		get_h5_string(h5file, entry + "/start_time", m_timestamp);
		get_h5_string(h5file, entry + "/end_time", timestamp_end);
		get_h5_scalar(h5file, entry + "/run_number", m_scannumber);
		get_h5_string(h5file, entry + "/" + m_instrname + "/command_line/actual_command", m_scancommand);

		// get polarisation infos
		get_h5_string(h5file, entry + "/" + m_instrname + "/pal/pal_contents", m_palcommand);
		if(!get_h5_scalar(h5file, entry + "/data_scan/pal_steps", m_numPolChannels))
			m_numPolChannels = 0;

		// get user infos
		get_h5_string(h5file, entry + "/user/name", m_username);
		get_h5_string(h5file, entry + "/user/namelocalcontact", m_localname);
		get_h5_string(h5file, entry + "/user/proposal", m_proposal);

		// get instrument infos
		t_real mono_sense = -1., ana_sense = -1., sample_sense = 1.;
		t_real ki = 0., kf = 0.;
		int fx = 2;

		get_h5_scalar(h5file, entry + "/" + m_instrname + "/Monochromator/d_spacing", m_dspacings[0]);
		if(!get_h5_scalar(h5file, entry + "/" + m_instrname + "/Monochromator/sense", mono_sense))
			get_h5_scalar(h5file, entry + "/" + m_instrname + "/Monochromator/sens", mono_sense);
		if(!get_h5_scalar(h5file, entry + "/" + m_instrname + "/Monochromator/ki", ki))
		{
			// if ki does not exist, try to convert from Ei
			t_real Ei = 0.;
			if(get_h5_scalar(h5file, entry + "/" + m_instrname + "/Monochromator/ei", Ei))
				ki = std::sqrt(get_E2KSQ<t_real>() * Ei);
		}
		get_h5_scalar(h5file, entry + "/" + m_instrname + "/Analyser/d_spacing", m_dspacings[1]);
		if(!get_h5_scalar(h5file, entry + "/" + m_instrname + "/Analyser/sense", ana_sense))
			get_h5_scalar(h5file, entry + "/" + m_instrname + "/Analyser/sens", ana_sense);
		if(!get_h5_scalar(h5file, entry + "/" + m_instrname + "/Analyser/kf", kf))
		{
			// if kf does not exist, try to convert from Ef
			t_real Ef = 0.;
			if(get_h5_scalar(h5file, entry + "/" + m_instrname + "/Analyser/ef", Ef))
				kf = std::sqrt(get_E2KSQ<t_real>() * Ef);
		}
		if(!get_h5_scalar(h5file, entry + "/sample/sense", sample_sense))
			get_h5_scalar(h5file, entry + "/sample/sens", sample_sense);
		get_h5_scalar(h5file, entry + "/sample/fx", fx);

		m_senses[0] = mono_sense > 0.;
		m_senses[1] = sample_sense > 0.;
		m_senses[2] = ana_sense > 0.;

		m_iskifixed = (fx == 1);
		m_kfix = (m_iskifixed ? ki : kf);

		// get sample infos
		get_h5_scalar(h5file, entry + "/sample/unit_cell_a", m_lattice[0]);
		get_h5_scalar(h5file, entry + "/sample/unit_cell_b", m_lattice[1]);
		get_h5_scalar(h5file, entry + "/sample/unit_cell_c", m_lattice[2]);

		get_h5_scalar(h5file, entry + "/sample/unit_cell_alpha", m_angles[0]);
		get_h5_scalar(h5file, entry + "/sample/unit_cell_beta", m_angles[1]);
		get_h5_scalar(h5file, entry + "/sample/unit_cell_gamma", m_angles[2]);

		get_h5_scalar(h5file, entry + "/sample/ax", m_plane[0][0]);
		get_h5_scalar(h5file, entry + "/sample/ay", m_plane[0][1]);
		get_h5_scalar(h5file, entry + "/sample/az", m_plane[0][2]);

		get_h5_scalar(h5file, entry + "/sample/bx", m_plane[1][0]);
		get_h5_scalar(h5file, entry + "/sample/by", m_plane[1][1]);
		get_h5_scalar(h5file, entry + "/sample/bz", m_plane[1][2]);

		get_h5_scalar(h5file, entry + "/sample/qh", m_initialpos[0]);
		get_h5_scalar(h5file, entry + "/sample/qk", m_initialpos[1]);
		get_h5_scalar(h5file, entry + "/sample/ql", m_initialpos[2]);
		get_h5_scalar(h5file, entry + "/sample/en", m_initialpos[3]);

		m_angles[0] = d2r(m_angles[0]);
		m_angles[1] = d2r(m_angles[1]);
		m_angles[2] = d2r(m_angles[2]);

		ublas::vector<t_real> plane_1 = make_vec<ublas::vector<t_real>>({ m_plane[0][0],  m_plane[0][1],  m_plane[0][2] });
		ublas::vector<t_real> plane_2 = make_vec<ublas::vector<t_real>>({ m_plane[1][0],  m_plane[1][1],  m_plane[1][2] });
		ublas::vector<t_real> plane_n = cross_3(plane_1, plane_2);

		// try to determine scanned variables from scan command
		std::vector<std::string> scanned_vars = FilePsi<t_real>::GetScannedVarsFromCommand(m_scancommand);
		for(auto iterScVar = scanned_vars.rbegin(); iterScVar != scanned_vars.rend(); ++iterScVar)
		{
			// bring the scanned variable to the front
			const std::string& scanned_var = *iterScVar;

			auto iterVar = std::find(m_scanned_vars.begin(), m_scanned_vars.end(), scanned_var);
			if(iterVar != m_scanned_vars.end())
			{
				m_scanned_vars.erase(iterVar);
				m_scanned_vars.insert(m_scanned_vars.begin(), scanned_var);
			}
		}

		// use first column in case no scan variables are given
		if(!m_scanned_vars.size())
		{
			log_warn("Could not determine scan variable.");
			if(m_vecCols.size() >= 1)
			{
				log_warn("Using first column: \"", m_vecCols[0], "\".");
				m_scanned_vars.push_back(m_vecCols[0]);
			}
		}

		// check consistency with respect to the number of scan steps
		std::size_t scan_steps = 0;
		std::size_t pal_steps = m_numPolChannels ? m_numPolChannels : 1;
		if(get_h5_scalar(h5file, entry + "/data_scan/actual_step", scan_steps)
			&& GetScanCount() != scan_steps*pal_steps)
		{
			log_warn("Determined ", GetScanCount(),
				" scan steps, but file reports ", scan_steps*pal_steps, ".");
		}

		h5file.close();

		// add parameters to metadata map
		m_params.emplace(std::make_pair("exp_title", m_title));
		m_params.emplace(std::make_pair("exp_user", m_username));
		m_params.emplace(std::make_pair("exp_localcontact", m_localname));

		m_params.emplace(std::make_pair("scan_time_start", m_timestamp));
		m_params.emplace(std::make_pair("scan_time_end", timestamp_end));
		m_params.emplace(std::make_pair("scan_number", var_to_str(m_scannumber)));
		m_params.emplace(std::make_pair("scan_command", m_scancommand));
		m_params.emplace(std::make_pair("scan_command_pol", m_palcommand));

		m_params.emplace(std::make_pair("sample_lattice",
			var_to_str(m_lattice[0], prec) + ", " +
			var_to_str(m_lattice[1], prec) + ", " +
			var_to_str(m_lattice[2], prec)));
		m_params.emplace(std::make_pair("sample_angles",
			var_to_str(r2d(m_angles[0]), prec) + ", " +
			var_to_str(r2d(m_angles[1]), prec) + ", " +
			var_to_str(r2d(m_angles[2]), prec)));
		m_params.emplace(std::make_pair("sample_plane_vec1",
			var_to_str(m_plane[0][0], prec) + ", " +
			var_to_str(m_plane[0][1], prec) + ", " +
			var_to_str(m_plane[0][2], prec)));
		m_params.emplace(std::make_pair("sample_plane_vec2",
			var_to_str(m_plane[1][0], prec) + ", " +
			var_to_str(m_plane[1][1], prec) + ", " +
			var_to_str(m_plane[1][2], prec)));
		m_params.emplace(std::make_pair("sample_plane_norm",
			var_to_str(plane_n[0], prec) + ", " +
			var_to_str(plane_n[1], prec) + ", " +
			var_to_str(plane_n[2], prec)));
		m_params.emplace(std::make_pair("sample_hklE",
			var_to_str(m_initialpos[0], prec) + ", " +
			var_to_str(m_initialpos[1], prec) + ", " +
			var_to_str(m_initialpos[2], prec) + ", " +
			var_to_str(m_initialpos[3], prec)));

		m_params.emplace(std::make_pair("instr_senses",
			var_to_str(m_senses[0]) + ", " +
			var_to_str(m_senses[1]) + ", " +
			var_to_str(m_senses[2])));
		m_params.emplace(std::make_pair("instr_ki", var_to_str(ki, prec)));
		m_params.emplace(std::make_pair("instr_kf", var_to_str(kf, prec)));
		m_params.emplace(std::make_pair("instr_ki_fixed", var_to_str(m_iskifixed)));
		m_params.emplace(std::make_pair("instr_kf_fixed", var_to_str(!m_iskifixed)));
		m_params.emplace(std::make_pair("instr_dspacings",
			var_to_str(m_dspacings[0], prec) + ", " +
			var_to_str(m_dspacings[1], prec)));

		if(m_bAutoParsePol)
			ParsePolData();
	}
	catch(const H5::Exception& ex)
	{
		log_err(ex.getDetailMsg());
		return false;
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
		return false;
	}

	return true;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecVals&
FileH5<t_real>::GetCol(const std::string& strName, std::size_t *pIdx) const
{
	return const_cast<FileH5*>(this)->GetCol(strName, pIdx);
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecVals&
FileH5<t_real>::GetCol(const std::string& strName, std::size_t *pIdx)
{
	static std::vector<t_real> vecNull;

	for(std::size_t i=0; i<m_vecCols.size(); ++i)
	{
		if(str_to_lower(m_vecCols[i]) == str_to_lower(strName))
		{
			if(pIdx) *pIdx = i;
			return m_data[i];
		}
	}

	if(pIdx)
		*pIdx = m_vecCols.size();

	log_err("Column \"", strName, "\" does not exist.");
	return vecNull;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecDat&
FileH5<t_real>::GetData() const
{
	return m_data;
}


template<class t_real>
typename FileInstrBase<t_real>::t_vecDat&
FileH5<t_real>::GetData()
{
	return m_data;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_vecColNames&
FileH5<t_real>::GetColNames() const
{
	return m_vecCols;
}


template<class t_real>
const typename FileInstrBase<t_real>::t_mapParams&
FileH5<t_real>::GetAllParams() const
{
	return m_params;
}


template<class t_real>
std::array<t_real,3> FileH5<t_real>::GetSampleLattice() const
{
	return m_lattice;
}


template<class t_real>
std::array<t_real,3> FileH5<t_real>::GetSampleAngles() const
{
	return m_angles;
}


template<class t_real>
std::array<t_real,2> FileH5<t_real>::GetMonoAnaD() const
{
	return m_dspacings;
}


template<class t_real>
std::array<bool, 3> FileH5<t_real>::GetScatterSenses() const
{
	return m_senses;
}


template<class t_real>
std::array<t_real, 3> FileH5<t_real>::GetScatterPlane0() const
{
	return m_plane[0];
}


template<class t_real>
std::array<t_real, 3> FileH5<t_real>::GetScatterPlane1() const
{
	return m_plane[1];
}


template<class t_real>
std::array<t_real, 4> FileH5<t_real>::GetPosHKLE() const
{
	return m_initialpos;
}


template<class t_real>
t_real FileH5<t_real>::GetKFix() const
{
	return m_kfix;
}


template<class t_real>
bool FileH5<t_real>::IsKiFixed() const
{
	return m_iskifixed;
}


template<class t_real>
std::size_t FileH5<t_real>::GetScanCount() const
{
	if(m_data.size() != 0)
		return m_data[0].size();
	return 0;
}


template<class t_real>
std::array<t_real, 5> FileH5<t_real>::GetScanHKLKiKf(std::size_t i) const
{
	using t_map = typename FileInstrBase<t_real>::t_mapParams;
	//const t_map& params = GetAllParams();

	return FileInstrBase<t_real>::GetScanHKLKiKf("QH", "QK", "QL", "EN", i);
}


template<class t_real> std::vector<std::string> FileH5<t_real>::GetScannedVars() const
{
	return m_scanned_vars;
}


template<class t_real>
bool FileH5<t_real>::MergeWith(const FileInstrBase<t_real>* pDat, bool allow_col_mismatch)
{
	return FileInstrBase<t_real>::MergeWith(pDat, allow_col_mismatch);
}


template<class t_real> std::string FileH5<t_real>::GetCountVar() const
{
	// try to match names
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX((Single)?Detector)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(cnts)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(det)REX", strRet))
		return strRet;

	return "";
}


template<class t_real>
void FileH5<t_real>::ParsePolData()
{
	std::string palcommand = m_palcommand;
	find_all_and_replace<std::string>(palcommand, "|", ",");

	m_vecPolStates = parse_pol_states<t_real>(palcommand,
		m_strPolVec1, m_strPolVec2,
		m_strPolCur1, m_strPolCur2,
		m_strXYZ,
		m_strFlip1, m_strFlip2);

	if(m_vecPolStates.size() && m_numPolChannels != m_vecPolStates.size())
	{
		log_warn("Determined ", m_vecPolStates.size(),
			" polarisation channels, but file reports ", m_numPolChannels, ".");
		m_numPolChannels = m_vecPolStates.size();
	}
}


template<class t_real>
std::size_t FileH5<t_real>::NumPolChannels() const
{
	return m_numPolChannels;
	//return m_vecPolStates.size();
}


template<class t_real>
const std::vector<std::array<t_real, 6>>& FileH5<t_real>::GetPolStates() const
{
	return m_vecPolStates;
}


template<class t_real>
void FileH5<t_real>::SetPolNames(const char* pVec1, const char* pVec2,
	const char* pCur1, const char* pCur2)
{
	m_strPolVec1 = pVec1; m_strPolVec2 = pVec2;
	m_strPolCur1 = pCur1; m_strPolCur2 = pCur2;
}


template<class t_real>
void FileH5<t_real>::SetLinPolNames(const char* pFlip1, const char* pFlip2,
	const char* pXYZ)
{
	m_strFlip1 = pFlip1; m_strFlip2 = pFlip2;
	m_strXYZ = pXYZ;
}


template<class t_real> void FileH5<t_real>::SetAutoParsePolData(bool b)
{
	m_bAutoParsePol = b;
}


template<class t_real> std::string FileH5<t_real>::GetMonVar() const
{
	std::string strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(Monitor[0-9]*)REX", strRet))
		return strRet;
	if(FileInstrBase<t_real>::MatchColumn(R"REX(m[0-9]+)REX", strRet))
		return strRet;

	return "";
}


template<class t_real> std::string FileH5<t_real>::GetTitle() const
{
	return m_title;
}


template<class t_real> std::string FileH5<t_real>::GetProposal() const
{
	return m_proposal;
}


template<class t_real> std::string FileH5<t_real>::GetUser() const
{
	return m_username;
}


template<class t_real> std::string FileH5<t_real>::GetInstrument() const
{
	return m_instrname;
}


template<class t_real> std::string FileH5<t_real>::GetLocalContact() const
{
	return m_localname;
}


template<class t_real> std::string FileH5<t_real>::GetScanNumber() const
{
	return var_to_str(m_scannumber);
}


template<class t_real> std::string FileH5<t_real>::GetScanCommand() const
{
	return m_scancommand;
}


template<class t_real> std::string FileH5<t_real>::GetTimestamp() const
{
	return m_timestamp;
}


template<class t_real> std::string FileH5<t_real>::GetSampleName() const { return ""; }
template<class t_real> std::string FileH5<t_real>::GetSpacegroup() const { return ""; }

}


#endif
