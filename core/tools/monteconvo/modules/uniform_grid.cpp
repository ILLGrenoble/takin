/**
 * S(Q, E) module for a pre-calculated grid (file format version 2)
 * @author Tobias Weber <tweber@ill.fr>
 * @date 06-jan-2020
 * @license GPLv2
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

#include "uniform_grid.h"

#include "libs/version.h"
#include "tlibs/string/string.h"
#include "tlibs/math/math.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/file/file.h"

#include <fstream>

using t_real = typename SqwUniformGrid::t_real;


// ----------------------------------------------------------------------------
// constructors

SqwUniformGrid::SqwUniformGrid()
{
	SqwBase::m_bOk = 0;
}


SqwUniformGrid::SqwUniformGrid(const std::string& strDatFile) : m_strDataFile(strDatFile)
{
	SqwBase::m_bOk = 0;
	tl::log_info("Loading grid version 2 data file: \"", strDatFile, "\".");


	if(!tl::file_exists(strDatFile.c_str()))
	{
		tl::log_err("Grid data file \"", strDatFile, "\" does not exist.");
		return;
	}

	std::ifstream ifstr(strDatFile);
	if(!ifstr)
	{
		tl::log_err("Grid data file \"", strDatFile, "\" cannot be opened.");
		return;
	}

	auto _blockOffs = tl::get_file_mem<std::size_t>(ifstr, 0, 1);
	auto _dims = tl::get_file_mem<t_real>(ifstr, 8, 9);

	if(!_blockOffs.first || !_blockOffs.second)
	{
		tl::log_err("Grid data file \"", strDatFile, "\" cannot be mapped.");
		return;
	}

	m_indexBlockOffset = *_blockOffs.second.get();

	m_hmin = _dims.second.get()[0];
	m_hmax = _dims.second.get()[1];
	m_hstep = _dims.second.get()[2];
	m_hsize = std::size_t(std::round((m_hmax - m_hmin) / m_hstep));

	m_kmin = _dims.second.get()[3];
	m_kmax = _dims.second.get()[4];
	m_kstep = _dims.second.get()[5];
	m_ksize = std::size_t(std::round((m_kmax - m_kmin) / m_kstep));

	m_lmin = _dims.second.get()[6];
	m_lmax = _dims.second.get()[7];
	m_lstep = _dims.second.get()[8];
	m_lsize = std::size_t(std::round((m_lmax - m_lmin) / m_lstep));

	std::size_t numEntries =
		std::size_t(((m_hmax-m_hmin) / m_hstep)) *
		std::size_t(((m_kmax-m_kmin) / m_kstep)) *
		std::size_t(((m_lmax-m_lmin) / m_lstep));

	tl::log_info("Data block dimensions: h=", m_hmin, "..", m_hmax, " (delta=", m_hstep, "), ",
		"k=", m_kmin, "..", m_kmax, " (delta=", m_kstep, "), ",
		"l=", m_lmin, "..", m_lmax, " (delta=", m_lstep, ").");
	tl::log_info("Number of entries: ", numEntries, ".");
	//tl::log_info("Offset of index block: ", m_indexBlockOffset, ".");

	SqwBase::m_bOk = 1;
}


SqwUniformGrid::~SqwUniformGrid()
{
}


// ----------------------------------------------------------------------------
// dispersion, spectral weight and structure factor

std::tuple<std::vector<t_real>, std::vector<t_real>>
	SqwUniformGrid::disp(t_real dh, t_real dk, t_real dl) const
{
	/**
	 * calculate file index based on coordinates
	 */
	auto hkl_to_idx = [this](t_real h, t_real k, t_real l) -> std::size_t
	{
		// clamp values to boundaries
		if(h < m_hmin) h = m_hmin;
		if(k < m_kmin) k = m_kmin;
		if(l < m_lmin) l = m_lmin;
		if(h >= m_hmax) h = m_hmax - m_hstep;
		if(k >= m_kmax) k = m_kmax - m_kstep;
		if(l >= m_lmax) l = m_lmax - m_lstep;

		// position indices
		std::size_t iH = std::size_t(std::round(((h - m_hmin) / m_hstep)));
		std::size_t iK = std::size_t(std::round(((k - m_kmin) / m_kstep)));
		std::size_t iL = std::size_t(std::round(((l - m_lmin) / m_lstep)));

                // clamp again
		if(iH >= m_hsize) iH = m_hsize-1;
		if(iK >= m_ksize) iK = m_ksize-1;
		if(iL >= m_lsize) iL = m_lsize-1;

		return iH*m_ksize*m_lsize + iK*m_lsize + iL;
	};



	if(!tl::file_exists(m_strDataFile.c_str()))
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" does not exist.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	std::ifstream ifstr(m_strDataFile);

	if(!ifstr)
	{
		tl::log_err("Data file \"", m_strDataFile, "\" cannot be opened.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}


	// ------------------------------------------------------------------------
	// the index block the offsets into the data block
	std::size_t idx_file_offs = hkl_to_idx(dh, dk, dl);
	auto _dat_file_offs = tl::get_file_mem<std::size_t>(ifstr,
		m_indexBlockOffset + idx_file_offs*sizeof(std::size_t), 1);

	if(!_dat_file_offs.first)
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" cannot be mapped.");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	std::size_t dat_file_offs = *_dat_file_offs.second.get();
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// the data block holds the energies and spectral weights of the dispersion branches

	auto _num_branches = tl::get_file_mem<unsigned int>(ifstr, dat_file_offs, 1);
	if(!_num_branches.first)
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" cannot be mapped (1).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	// number of dispersion branches and weights
	unsigned int iNumBranches = *_num_branches.second.get();


	// map actual (E, w) data
	auto _branches = tl::get_file_mem<t_real>(ifstr,
		dat_file_offs+sizeof(iNumBranches), iNumBranches*2);
	if(!_branches.first)
	{
		tl::log_err("Grid data file \"", m_strDataFile, "\" cannot be mapped (2).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	const t_real *pBranches = _branches.second.get();

	std::vector<t_real> vecE, vecw;
	for(unsigned int iBranch=0; iBranch<iNumBranches; ++iBranch)
	{
		if(!tl::float_equal(pBranches[iBranch*2 + 1], t_real(0)))
		{
			vecE.push_back(pBranches[iBranch*2 + 0]);	// energy
			vecw.push_back(pBranches[iBranch*2 + 1]);	// weight
		}
	}

	// ------------------------------------------------------------------------


	return std::make_tuple(vecE, vecw);
}


/**
 * S(q,E)
 */
t_real SqwUniformGrid::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	std::vector<t_real> vecE, vecW;
	std::tie(vecE, vecW) = disp(dh, dk, dl);

	t_real dInc = 0;
	if(!tl::float_equal(m_dIncAmp, t_real(0)))
		dInc = tl::gauss_model(dE, t_real(0), m_dIncSigma, m_dIncAmp, t_real(0));

	t_real dS = 0;
	for(std::size_t iE=0; iE<vecE.size(); ++iE)
		dS += tl::gauss_model(dE, vecE[iE], m_dSigma, vecW[iE], t_real(0));

	return m_dS0*dS * tl::bose_cutoff(dE, m_dT, m_dcut) + dInc;
}



// ----------------------------------------------------------------------------
// get & set variables

std::vector<SqwUniformGrid::t_var> SqwUniformGrid::GetVars() const
{
	std::vector<t_var> vecVars;

	vecVars.push_back(SqwBase::t_var{"T", "real", tl::var_to_str(m_dT)});
	vecVars.push_back(SqwBase::t_var{"bose_cutoff", "real", tl::var_to_str(m_dcut)});
	vecVars.push_back(SqwBase::t_var{"sigma", "real", tl::var_to_str(m_dSigma)});
	vecVars.push_back(SqwBase::t_var{"inc_amp", "real", tl::var_to_str(m_dIncAmp)});
	vecVars.push_back(SqwBase::t_var{"inc_sigma", "real", tl::var_to_str(m_dIncSigma)});
	vecVars.push_back(SqwBase::t_var{"S0", "real", tl::var_to_str(m_dS0)});

	return vecVars;
}


void SqwUniformGrid::SetVars(const std::vector<SqwUniformGrid::t_var>& vecVars)
{
	if(!vecVars.size()) return;

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		if(strVar == "T") m_dT = tl::str_to_var<decltype(m_dT)>(strVal);
		else if(strVar == "bose_cutoff") m_dcut = tl::str_to_var<decltype(m_dcut)>(strVal);
		else if(strVar == "sigma") m_dSigma = tl::str_to_var<t_real>(strVal);
		else if(strVar == "inc_amp") m_dIncAmp = tl::str_to_var<decltype(m_dIncAmp)>(strVal);
		else if(strVar == "inc_sigma") m_dIncSigma = tl::str_to_var<decltype(m_dIncSigma)>(strVal);
		else if(strVar == "S0") m_dS0 = tl::str_to_var<decltype(m_dS0)>(strVal);
	}
}


bool SqwUniformGrid::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	return SqwBase::SetVarIfAvail(strKey, strNewVal);
}



// ----------------------------------------------------------------------------
// copy

SqwBase* SqwUniformGrid::shallow_copy() const
{
	SqwUniformGrid *pMod = new SqwUniformGrid();

	pMod->m_dT = this->m_dT;
	pMod->m_dcut = this->m_dcut;
	pMod->m_dSigma = this->m_dSigma;
	pMod->m_dIncAmp = this->m_dIncAmp;
	pMod->m_dIncSigma = this->m_dIncSigma;
	pMod->m_dS0 = this->m_dS0;

	pMod->m_strDataFile = this->m_strDataFile;
	pMod->m_indexBlockOffset = this->m_indexBlockOffset;

	pMod->m_hmin = this->m_hmin;
	pMod->m_hmax = this->m_hmax;
	pMod->m_hstep = this->m_hstep;
	pMod->m_hsize = this->m_hsize;

	pMod->m_kmin = this->m_kmin;
	pMod->m_kmax = this->m_kmax;
	pMod->m_kstep = this->m_kstep;
	pMod->m_ksize = this->m_ksize;

	pMod->m_lmin = this->m_lmin;
	pMod->m_lmax = this->m_lmax;
	pMod->m_lstep = this->m_lstep;
	pMod->m_lsize = this->m_lsize;

	return pMod;
}
