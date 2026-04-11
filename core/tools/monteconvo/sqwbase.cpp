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

#include "sqwbase.h"
#include "tlibs/log/log.h"


/**
 * copy member variables
 */
const SqwBase& SqwBase::operator=(const SqwBase& sqw)
{
	this->m_bOk = sqw.m_bOk;
	this->m_vecFit = sqw.m_vecFit;

	return *this;
}


/**
 * set or override S(Q, E) model parameters from a string
 */
std::tuple<bool, std::string> SqwBase::SetVars(
  const std::string& sqw_params, bool print_messages,
	const std::unordered_set<std::string> *ingored_params,
	std::unordered_map<std::string, std::string> *all_params)
{
	bool ok = true;

	// set the given individual global model parameters
	if(sqw_params == "")
		return std::make_pair(ok, "");

	std::ostringstream ostrStatus;

	// get all parameters
	std::vector<std::string> vecSetParams;
	tl::get_tokens<std::string, std::string>(sqw_params, ";", vecSetParams);

	for(const std::string& strModParam : vecSetParams)
	{
		// get key = value pairs for current parameter
		std::vector<std::string> vecModParam;
		tl::get_tokens<std::string, std::string>(strModParam, "=", vecModParam);
		if(vecModParam.size() < 2)
			continue;

		std::string key = vecModParam[0];
		std::string val = vecModParam[1];
		tl::trim(key);
		tl::trim(val);

		// save all parameters to a map
		if(all_params)
			all_params->insert(std::make_pair(key, val));

		// skip ignored parameters
		if(ingored_params && ingored_params->find(key) != ingored_params->end())
			continue;

		if(SetVarIfAvail(key, val))
		{
			// list changed parameters in a readable form
			if(ostrStatus.str().length())
				ostrStatus << ", ";
			ostrStatus << key << " = " << val;

			if(print_messages)
				tl::log_info("Setting S(Q, E) model parameter \"", key, "\" to \"", val, "\".");
		}
		else
		{
			ok = false;

			tl::log_err("No parameter named \"", key, "\" available in S(Q, E) model, ",
				"cannot set new value \"", val, "\".");
		}
	}

	return std::make_pair(ok, ostrStatus.str());
}


/**
 * if the variable "strKey" is known, update it with the value "strNewVal"
 */
bool SqwBase::SetVarIfAvail(const std::string& strKey, const std::string& strNewVal)
{
	std::vector<t_var> vecVars = GetVars();
	for(const t_var& var : vecVars)
	{
		if(strKey != std::get<0>(var))
			continue;

		t_var varNew = var;
		std::get<2>(varNew) = strNewVal;
		SetVars(std::vector<t_var>({ varNew }));
		return true;
	}

	return false;
}


/**
 * if the variable "strKey" is known, update its error with the value "strNewErr"
 */
bool SqwBase::SetErrIfAvail(const std::string& strKey, const std::string& strNewErr)
{
	for(t_var_fit& var : m_vecFit)
	{
		if(strKey != std::get<0>(var))
			continue;

		std::get<1>(var) = strNewErr;
		return true;
	}

	return false;
}


/**
 * if the variable "strKey" is known, update its range with the value "strNewRange"
 */
bool SqwBase::SetRangeIfAvail(const std::string& strKey, const std::string& strNewRange)
{
	for(t_var_fit& var : m_vecFit)
	{
		if(strKey != std::get<0>(var))
			continue;

		std::get<3>(var) = strNewRange;
		return true;
	}

	return false;
}


/**
 * replaces all fit variables
 */
void SqwBase::InitFitVars(const std::vector<t_var_fit>& vecFit)
{
	m_vecFit = vecFit;
}


/**
 * keeps the current fit variable vector and only replaces changed values
 */
void SqwBase::SetFitVars(const std::vector<t_var_fit>& vecFitVars)
{
	using t_iter = typename decltype(m_vecFit)::iterator;

	for(const t_var_fit& var : vecFitVars)
	{
		const std::string& name = std::get<0>(var);
		t_iter iter = std::find_if(m_vecFit.begin(), m_vecFit.end(), [&name](const t_var_fit& fitvar)
		{
			return std::get<0>(fitvar) == name;
		});

		if(iter != m_vecFit.end())
		{
			std::get<1>(*iter) = std::get<1>(var);
			std::get<2>(*iter) = std::get<2>(var);
			std::get<3>(*iter) = std::get<3>(var);
		}
		else
		{
			tl::log_err("Tried to update non-existing fit variable \"", name, "\".");
		}
	}
}
