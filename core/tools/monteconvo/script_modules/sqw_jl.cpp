/**
 * S(Q, E) julia interface
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "sqw_jl.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"
#include "tlibs/ext/jl.h"

using t_real = t_real_reso;


static const char* pcModIdent = "jl";
static const char* pcModName = "Julia Scripting";
static const char* pcModHelp = R"RAWSTR(Julia Scripting Module.

This module allows for the inclusion of user-provided S(Q, E) Julia scripts for use with the convolution simulator or fitter.

Please refer to the Takin help for more information and tutorials.)RAWSTR";


#define MAX_PARAM_VAL_SIZE 128


//extern "C" void jl_init__threading();


SqwJl::SqwJl(const std::string& strFile) : m_pmtx(std::make_shared<std::mutex>())
{
	if(!tl::file_exists(strFile.c_str()))
	{
		tl::log_err("Could not find Julia script file: \"", strFile, "\".");
		m_bOk = false;
		return;
	}

	std::string strDir = tl::get_dir(strFile);
	const bool bSetScriptCWD = true;

	// init interpreter
	static bool bInited = false;
	if(!bInited)
	{
		//jl_init__threading();
		jl_init();
		std::string strJl = jl_ver_string();
		tl::log_debug("Initialised Julia interpreter version ", strJl, ".");
		bInited = true;
	}

	// get include function
	jl_function_t *pInc = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_base_module, jl_symbol("include")));
	if(!pInc)
	{
		m_bOk = false;
		tl::log_err("Cannot get Julia include() function.");
		return;
	}

	// include module
	jl_value_t* pModIncRet = jl_call2(pInc, (jl_value_t*)jl_main_module,
		jl_cstr_to_string(strFile.c_str()));
	if(!pModIncRet)
	{
		tl::log_err("Cannot load Julia script: \"", strFile, "\".");
		PrintExceptions();
	}


	// working dir
	if(bSetScriptCWD)
	{
		jl_function_t *pCwd = reinterpret_cast<jl_function_t*>(
			jl_get_global(jl_base_module, jl_symbol("cd")));
		jl_value_t *pDir = jl_cstr_to_string(strDir.c_str());

		if(pCwd && pDir)
			jl_call1(pCwd, pDir);
		else
			tl::log_err("Cannot get Julia cd() function.");
	}


	// import takin functions
	m_pInit = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_main_module, jl_symbol("TakinInit")));
	m_pSqw = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_main_module, jl_symbol("TakinSqw")));
	m_pDisp = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_main_module, jl_symbol("TakinDisp")));
	m_pBackground = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_main_module, jl_symbol("TakinBackground")));

	PrintExceptions();


	if(m_pInit)
		tl::log_info("TakinInit function was found in \"", strFile, "\".");
	else
		tl::log_warn("No TakinInit function was found in \"", strFile, "\".");

	if(m_pSqw)
		tl::log_info("TakinSqw function was found in \"", strFile, "\".");
	else
		tl::log_err("No TakinSqw function was found in \"", strFile, "\".");

	if(m_pDisp)
		tl::log_info("TakinDisp function was found in \"", strFile, "\".");
	else
		tl::log_warn("No TakinDisp function was found in \"", strFile, "\".");

	if(m_pBackground)
		tl::log_info("TakinBackground function was found in \"", strFile, "\".");
	else
		tl::log_warn("No TakinBackground function was found in \"", strFile, "\".");

	// does the module have a TakinSqw function?
	if(!m_pSqw)
	{
		m_bOk = false;
		return;
	}
	else
	{
		m_bOk = true;
	}

	if(m_pInit)
		jl_call0((jl_function_t*)m_pInit);
}


SqwJl::~SqwJl()
{
	// TODO: cleanup...
}


/**
 * converts a jl_value_t to a string
 */
std::string SqwJl::GetJlString(void* _pVal) const
{
	jl_value_t* pVal = (jl_value_t*)_pVal;
	if(!pVal)
		return "";

	std::string str;
	jl_value_t *pStr = nullptr;

	jl_function_t *pToStr = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_base_module, jl_symbol("string")));
	if(!pToStr)
	{
		tl::log_err("Cannot get Julia string() function.");
		return "";
	}

	if(pToStr)
		pStr = jl_call1(pToStr, pVal);
	if(pStr)
		str = jl_string_ptr(pStr);

	return str;
}


/**
 * Checks for and prints possible exceptions
 */
void SqwJl::PrintExceptions() const
{
	jl_value_t* pEx = jl_exception_occurred();

	if(pEx)
	{
		std::string strEx = GetJlString(pEx);

		if(strEx != "")
			tl::log_err("Julia error: ", strEx, ".");
		else
			tl::log_err("Unknown Julia error occurred.");
	}

	jl_exception_clear();
}


/**
 * dispersion relation, E(Q)
 */
std::tuple<std::vector<t_real>, std::vector<t_real>>
	SqwJl::disp(t_real dh, t_real dk, t_real dl) const
{
	if(!m_bOk)
	{
		tl::log_err("Julia interpreter has not initialised, cannot query S(Q, E).");
		return std::make_tuple(std::vector<t_real>(), std::vector<t_real>());
	}

	std::lock_guard<std::mutex> lock(*m_pmtx);
	std::vector<t_real> vecE, vecW;

	if(m_pDisp)
	{
		jl_value_t *phkl[3] =
		{
			tl::jl_traits<t_real>::box(dh),
			tl::jl_traits<t_real>::box(dk),
			tl::jl_traits<t_real>::box(dl)
		};
		jl_array_t *pEW = reinterpret_cast<jl_array_t*>(
			jl_call((jl_function_t*)m_pDisp, phkl, 3));

		if(jl_array_len(pEW) != 2)
		{
			tl::log_err("TakinDisp has to return [arrEnergies, arrWeights].");
			return std::make_tuple(vecE, vecW);
		}

		jl_array_t *parrE = jl_array_data(pEW, jl_array_t*)[0];
		jl_array_t *parrW = jl_array_data(pEW, jl_array_t*)[1];

		std::size_t iSizeE = jl_array_len(parrE);
		std::size_t iSizeW = jl_array_len(parrW);

		if(iSizeE != iSizeW)
			tl::log_warn("Size mismatch between energies and weights array in Julia script.");

		const t_real *Es = jl_array_data(parrE, t_real);
		const t_real *Ws = jl_array_data(parrW, t_real);

		for(std::size_t iElem = 0; iElem < std::min(iSizeE, iSizeW); ++iElem)
		{
			vecE.push_back(Es[iElem]);
			vecW.push_back(Ws[iElem]);
		}
	}

	PrintExceptions();
	return std::make_tuple(vecE, vecW);
}


/**
 * dynamical structure factor, S(Q,E)
 */
t_real SqwJl::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	if(!m_bOk)
	{
		tl::log_err("Julia interpreter has not initialised, cannot query S(Q, E).");
		return t_real(0);
	}

	std::lock_guard<std::mutex> lock(*m_pmtx);

	jl_value_t *phklE[4] =
	{
		tl::jl_traits<t_real>::box(dh),
		tl::jl_traits<t_real>::box(dk),
		tl::jl_traits<t_real>::box(dl),
		tl::jl_traits<t_real>::box(dE)
	};
	jl_value_t *pSqw = jl_call((jl_function_t*)m_pSqw, phklE, 4);

	PrintExceptions();
	return t_real(tl::jl_traits<t_real>::unbox(pSqw));
}


/**
 * background function
 */
t_real SqwJl::GetBackground(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	if(!m_bOk)
	{
		tl::log_err("Julia interpreter has not initialised, cannot query background.");
		return t_real(0);
	}

	std::lock_guard<std::mutex> lock(*m_pmtx);
	t_real val = 0.;

	if(m_pBackground)
	{
		jl_value_t *phklE[4] =
		{
			tl::jl_traits<t_real>::box(dh),
			tl::jl_traits<t_real>::box(dk),
			tl::jl_traits<t_real>::box(dl),
			tl::jl_traits<t_real>::box(dE)
		};
		jl_value_t *pVal = jl_call((jl_function_t*)m_pBackground, phklE, 4);

		val = t_real(tl::jl_traits<t_real>::unbox(pVal));
	}

	PrintExceptions();
	return val;
}


std::vector<SqwBase::t_var> SqwJl::GetVars() const
{
	std::vector<SqwBase::t_var> vecVars;
	if(!m_bOk)
	{
		tl::log_err("Julia interpreter has not initialised, cannot get variables.");
		return vecVars;
	}

	jl_function_t *pNames = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_base_module, jl_symbol("names")));
	jl_function_t *pGetField = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_base_module, jl_symbol("getfield")));
	jl_function_t *pPrint = reinterpret_cast<jl_function_t*>(
		jl_get_global(jl_base_module, jl_symbol("string")));

	if(!pNames || !pGetField || !pPrint)
	{
		tl::log_err("Required Julia functions not available.");
		return vecVars;
	}

	jl_array_t* pArrNames = (jl_array_t*)jl_call1(pNames, (jl_value_t*)jl_main_module);
	if(!pArrNames)
		return vecVars;

	std::size_t iSyms = jl_array_len(pArrNames);
	for(std::size_t iSym=0; iSym<iSyms; ++iSym)
	{
		jl_sym_t* pSym = (jl_sym_t*)jl_array_ptr_ref(pArrNames, iSym);
		if(!pSym) continue;

		// name
		std::string strName = jl_symbol_name(pSym);
		if(strName.length() == 0) continue;

		// filter out non-prefixed variables
		if(m_strVarPrefix.size() && strName.substr(0,m_strVarPrefix.size()) != m_strVarPrefix)
			continue;

		// type
		jl_value_t* pFld = jl_call2(pGetField, (jl_value_t*)jl_main_module, (jl_value_t*)pSym);
		if(!pFld) continue;
		std::string strType = jl_typeof_str(pFld);
		if(strType.length() == 0) continue;
		if(strType[0] == '#' || strType == "Module") continue;	// filter funcs and mods

		// value
		jl_value_t* pFldPr = jl_call1(pPrint, pFld);
		if(!pFldPr) continue;
		std::string strValue = jl_string_ptr(pFldPr);

		SqwBase::t_var var;
		std::get<0>(var) = std::move(strName);
		std::get<1>(var) = std::move(strType);
		std::get<2>(var) = std::move(strValue);

		vecVars.push_back(var);
	}

	PrintExceptions();
	return vecVars;
}


void SqwJl::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(!m_bOk)
	{
		tl::log_err("Julia interpreter has not initialised, cannot set variables.");
		return;
	}

	std::ostringstream ostrEval;
	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strName = std::get<0>(var);
		const std::string& strType = std::get<1>(var);
		const std::string& strValue = std::get<2>(var);

		if(!strName.length()) continue;

		ostrEval << strName << " = ";

		if(strType.length())
		{
			// if a type is given, filter out some names
			if(strType[0] == '#' || strType == "Module")
				continue;

			// with cast
			ostrEval << strType << "(" << strValue << ");\n";
		}
		else
		{
			//without cast
			ostrEval << strValue << ";\n";
		}
	}
	jl_eval_string(ostrEval.str().c_str());

	// TODO: check for changed parameters and if reinit is needed
	if(m_pInit)
		jl_call0((jl_function_t*)m_pInit);

	PrintExceptions();
}


SqwBase* SqwJl::shallow_copy() const
{
	SqwJl* pSqw = new SqwJl();
	*static_cast<SqwBase*>(pSqw) = *static_cast<const SqwBase*>(this);

	pSqw->m_pInit = this->m_pInit;
	pSqw->m_pSqw = this->m_pSqw;
	pSqw->m_pDisp = this->m_pDisp;
	pSqw->m_pmtx = this->m_pmtx;

	return pSqw;
}



// ----------------------------------------------------------------------------
// SO interface
#ifdef BUILD_PLUGIN

#include <boost/dll/alias.hpp>
#include "../sqw_proc.h"
#include "../sqw_proc_impl.h"
#include "libs/version.h"


std::tuple<std::string, std::string, std::string, std::string> sqw_info()
{
	//tl::log_info("In ", __func__, ".");
	return std::make_tuple(TAKIN_VER, pcModIdent, pcModName, pcModHelp);
}

std::shared_ptr<SqwBase> sqw_construct(const std::string& strCfgFile)
{
	//tl::log_info("In ", __func__, ".");
	return std::make_shared<SqwProc<SqwJl>>(strCfgFile);
}


SqwBase* sqw_construct_raw(const std::string& strCfgFile)
{
	return new SqwProc<SqwJl>(strCfgFile);
}


void sqw_destruct_raw(SqwBase *sqw)
{
	if(sqw)
		delete sqw;
}


// exports from so file
BOOST_DLL_ALIAS(sqw_info, takin_sqw_info);

// construction interface
BOOST_DLL_ALIAS(sqw_construct, takin_sqw);

// alternate raw construction interface
//BOOST_DLL_ALIAS(sqw_construct_raw, takin_sqw_new);
//BOOST_DLL_ALIAS(sqw_destruct_raw, takin_sqw_del);

#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// external program
#ifdef BUILD_APPLI

#include "../sqw_proc.h"
#include "../sqw_proc_impl.h"
#include "libs/version.h"


int main(int argc, char** argv)
{
	if(argc <= 2)
	{
		std::cout << "#\n# This is a Takin plugin module.\n#\n";
		std::cout << "module_ident: " << pcModIdent << "\n";
		std::cout << "module_name: " << pcModName << "\n";
		std::cout << "module_type: sqw\n";
		std::cout << "required_takin_version: " << TAKIN_VER << "\n";
		std::cout << "module_help: begin\n" << pcModHelp << "\nmodule_help: end\n";
		std::cout.flush();
		return 0;
	}

	const char* pcCfgFile = argv[1];
	const char* pcSharedMem = argv[2];
	SqwProc<SqwJl> proc(pcCfgFile, SqwProcStartMode::START_CHILD, pcSharedMem);

	return 0;
}

#endif
// ----------------------------------------------------------------------------
