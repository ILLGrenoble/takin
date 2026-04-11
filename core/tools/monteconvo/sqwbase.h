/**
 * interface for S(Q, E) models
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015, 2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MCONV_SQW_BASE_H__
#define __MCONV_SQW_BASE_H__

#include <string>
#include <tuple>
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>

#ifndef __MCONV_OWN_DEFS__
	#include "../res/defs.h"
#else
	using t_real_reso = double;
#endif

#include "tlibs/string/string.h"


/**
 * base class for S(Q, E) models
 */
class SqwBase
{
public:
	// basic fields: [ ident, type, value ]
	using t_var = std::tuple<std::string, std::string, std::string>;

	// extended fields: [ ident, error, is fit var?, range ]
	using t_var_fit = std::tuple<std::string, std::string, bool, std::string>;


protected:
	bool m_bOk = false;
	std::vector<t_var_fit> m_vecFit{};


public:
	/**
	 * E(Q) dispersion and spectral weight function (optional)
	 * return [energies, weights]
	 */
	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso /*dh*/, t_real_reso /*dk*/, t_real_reso /*dl*/) const
	{
		return std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>({}, {});
	}

	// S(Q,E) dynamical structure factor function which is queried for every mc point
	virtual t_real_reso operator()(
		t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const = 0;

	// background which is queried for every nominal (Q, E) point
	virtual t_real_reso GetBackground(
		t_real_reso /*dh*/, t_real_reso /*dk*/, t_real_reso /*dl*/, t_real_reso /*dE*/) const
	{
		return 0.;
	}

	virtual bool IsOk() const
	{
		return m_bOk;
	}

	// return model variables
	virtual std::vector<t_var> GetVars() const = 0;
	virtual const std::vector<t_var_fit>& GetFitVars() const
	{
		return m_vecFit;
	}

	// updates model variables
	virtual std::tuple<bool, std::string> SetVars(const std::string& sqw_params,
		bool print_messages = true,
		const std::unordered_set<std::string> *ingored_params = nullptr,
		std::unordered_map<std::string, std::string> *all_params = nullptr);
	virtual void SetVars(const std::vector<t_var>&) = 0;
	virtual void SetFitVars(const std::vector<t_var_fit>&);

	// this replaces the current fit variables with the given ones
	virtual void InitFitVars(const std::vector<t_var_fit>&);
	virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal);

	virtual bool SetErrIfAvail(const std::string& strKey, const std::string& strNewErr);
	virtual bool SetRangeIfAvail(const std::string& strKey, const std::string& strNewRange);

	SqwBase() = default;
	virtual ~SqwBase() = default;

	virtual const SqwBase& operator=(const SqwBase& sqw);
	SqwBase(const SqwBase& sqw)
	{
		this->operator=(sqw);
	}

	virtual SqwBase* shallow_copy() const = 0;
};


// ----------------------------------------------------------------------------


template<class t_vec = std::vector<double>>
std::string vec_to_str(const t_vec& vec)
{
	std::ostringstream ostr;
	for(const typename t_vec::value_type& t : vec)
		ostr << t << " ";

	return ostr.str();
}


template<class t_vec = std::vector<double>>
t_vec str_to_vec(const std::string& str)
{
	using T = typename t_vec::value_type;

	std::vector<T> vec0;
	tl::get_tokens<T, std::string, std::vector<T>>(str, " \t", vec0);

	t_vec vec(vec0.size());
	for(std::size_t idx = 0; idx < vec0.size(); ++idx)
		vec[idx] = vec0[idx];

	return vec;
}

// ----------------------------------------------------------------------------

#endif
