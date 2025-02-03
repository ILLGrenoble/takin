/**
 * factory and plugin interface for S(Q, E) models
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2016 -- 2020
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "sqwfactory.h"

#include "modules/simple_phonon.h"
#include "modules/simple_magnon.h"
#include "modules/kdtree.h"
#include "modules/table1d.h"
#include "modules/elast.h"
#include "modules/uniform_grid.h"

#include "sqw_proc.h"
#include "sqw_proc_impl.h"
#include "sqwrawdelegate.h"
#include "sqwnull.h"

#include "tlibs/log/log.h"
#include "tlibs/file/file.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/proc.h"
#include "libs/globals.h"
#include "libs/version.h"

#include <algorithm>
#include <functional>
#include <unordered_map>


// sqw info function: "takin_sqw_info"
// returns: [takin ver, ident, long name]
using t_pfkt_info = std::tuple<std::string, std::string, std::string, std::string>(*)();
using t_fkt_info = typename std::remove_pointer<t_pfkt_info>::type;

// sqw module creation function: "takin_sqw"
// old interface returning a shared pointer, which might be dangerous for so files,
// see: https://www.boost.org/doc/libs/1_72_0/doc/html/boost_dll/missuses.html
using t_pfkt = std::shared_ptr<SqwBase>(*)(const std::string&);
using t_fkt = typename std::remove_pointer<t_pfkt>::type;

// new raw pointer constructor interface
using t_pfkt_raw_new = SqwBase*(*)(const std::string&);
using t_fkt_raw_new = typename std::remove_pointer<t_pfkt_raw_new>::type;
using t_pfkt_raw_del = void(*)(SqwBase*);
using t_fkt_raw_del = typename std::remove_pointer<t_pfkt_raw_del>::type;


// key: identifier, value: [func, long name, help text]
using t_mapSqw = std::unordered_map<std::string, std::tuple<t_pfkt, std::string, std::string>>;
using t_mapSqwRaw = std::unordered_map<std::string, std::tuple<t_pfkt_raw_new, t_pfkt_raw_del, std::string, std::string>>;

// key: identifier, value: [long name, binary file name, help text]
using t_mapSqwExt = std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>>;



// shared_ptr constructors
static t_mapSqw g_mapSqw =
{
	{ "kd", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwKdTree>(strCfgFile.c_str()); },
			"4D Nearest-Point Raster",
			"This modules uses a 4-d nearest-point raster of the Form (h, k, l, E, S). "
			"\n\nIt is based on querying a k-d tree for the nearest point, "
			"so please make sure to include zeros around the dispersion branches."
		}
	},
	{ "table_1d", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwTable1d>(strCfgFile.c_str()); },
			"1D Nearest-Point Raster",
			"This module uses a 1-d nearest-point raster of the Form (q, E, S). "
			"\n\nIt is based on querying a k-d tree for the nearest point, "
			"so please make sure to include zeros around the dispersion branches."
		}
	},
	{ "uniform_grid", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwUniformGrid>(strCfgFile.c_str()); },
			"Uniform Grid",
			"A uniformly-spaced grid in Q saving the individual dispersion branches per Q. "
			"\n\nSee the Takin help for more information.",
		}
	},
	{ "phonon", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwPhonon>(strCfgFile.c_str()); },
			"Simple Phonon Model",
			"A simple model for longitudinal and transversal acoustic phonons."
		}
	},
	{ "phonon_single", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwPhononSingleBranch>(strCfgFile.c_str()); },
			"Simple Single-Branch Phonon Model",
			"This module creates a simple sine-shaped dispersion for acoustic phonons."
		}
	},
	{ "magnon", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwMagnon>(strCfgFile.c_str()); },
			"Simple Magnon Model",
			"This module creates simple ferromagnetic and antiferromagnetic dispersions. "
			"\n\nPlease consider using the \"Magnetic Dynamics\" module instead for more advanced cases."
		}
	},
	{ "elastic", t_mapSqw::mapped_type
		{
			[](const std::string& strCfgFile) -> std::shared_ptr<SqwBase>
			{ return std::make_shared<SqwElast>(strCfgFile.c_str()); },
			"Elastic Model",
			"This model creates a collection of Bragg peaks."
		}
	},
};


// raw pointer constructors
static t_mapSqwRaw g_mapSqwRaw;

// external process plugins
static t_mapSqwExt g_mapSqwExt;


/**
 * gets the short and long names of the installed plugin modules
 * as well as the help texts
 */
std::vector<std::tuple<std::string, std::string, std::string>> get_sqw_names()
{
	using t_tup = std::tuple<std::string, std::string, std::string>;
	std::vector<t_tup> vec;
	vec.reserve(g_mapSqw.size());

	for(const t_mapSqw::value_type& val : g_mapSqw)
	{
		t_tup tup;
		std::get<0>(tup) = val.first;
		std::get<1>(tup) = std::get<1>(val.second);
		std::get<2>(tup) = std::get<2>(val.second);

		vec.emplace_back(std::move(tup));
	}

	for(const t_mapSqwRaw::value_type& val : g_mapSqwRaw)
	{
		t_tup tup;
		std::get<0>(tup) = val.first;
		std::get<1>(tup) = std::get<2>(val.second);
		std::get<2>(tup) = std::get<3>(val.second);

		vec.emplace_back(std::move(tup));
	}

	for(const t_mapSqwExt::value_type& val : g_mapSqwExt)
	{
		t_tup tup;
		std::get<0>(tup) = val.first;
		std::get<1>(tup) = std::get<0>(val.second);
		std::get<2>(tup) = std::get<2>(val.second);

		vec.emplace_back(std::move(tup));
	}

	std::sort(vec.begin(), vec.end(), [](const t_tup& tup0, const t_tup& tup1) -> bool
	{
		const std::string& str0 = std::get<1>(tup0);
		const std::string& str1 = std::get<1>(tup1);

		return std::lexicographical_compare(str0.begin(), str0.end(), str1.begin(), str1.end());
	});

	return vec;
}


std::shared_ptr<SqwBase> construct_sqw(const std::string& strName,
	const std::string& strConfigFile)
{
	typename t_mapSqw::const_iterator iter = g_mapSqw.find(strName);
	typename t_mapSqwRaw::const_iterator iterRaw = g_mapSqwRaw.find(strName);
	typename t_mapSqwExt::const_iterator iterExt = g_mapSqwExt.find(strName);

	if(iter != g_mapSqw.end())
	{
		t_pfkt pFkt = std::get<0>(iter->second);
		if(!pFkt)
		{
			tl::log_err("Invalid constructor function for S(Q, E) model.");
			return nullptr;
		}

		tl::log_debug("Constructing \"", iter->first, "\" S(Q, E) module.");
		return (*pFkt)(strConfigFile);
	}
	else if(iterRaw != g_mapSqwRaw.end())
	{
		t_pfkt_raw_new pFktNew = std::get<0>(iterRaw->second);
		t_pfkt_raw_del pFktDel = std::get<1>(iterRaw->second);

		if(!pFktNew || !pFktDel)
		{
			tl::log_err("Invalid constructor function for S(Q, E) model.");
			return nullptr;
		}

		tl::log_debug("Constructing \"", iterRaw->first, "\" S(Q, E) module via raw interface.");
		return std::make_shared<SqwRawDelegate>(pFktNew(strConfigFile));
	}
	else if(iterExt != g_mapSqwExt.end())
	{
		tl::log_debug("Constructing \"", iterExt->first, "\" S(Q, E) module via external interface.");

		// number of spawned child processes
		std::size_t iNumProcesses = get_max_processes();

		return std::make_shared<SqwProc<SqwNull>>(strConfigFile.c_str(),
			SqwProcStartMode::START_PARENT_CREATE_CHILD, nullptr,
			std::get<1>(iterExt->second).c_str(), iNumProcesses);
	}

	tl::log_err("No S(Q, E) model of name \"", strName, "\" found.");
	return nullptr;
}




// --------------------------------------------------------------------------------
// external process plugins
// --------------------------------------------------------------------------------

// tracking modules for refcounting
static std::vector<std::string> g_vecExtMods;

void unload_sqw_ext_plugins()
{
	g_mapSqwExt.clear();
	g_vecExtMods.clear();

	tl::log_debug("Unloaded all external plugins.");
}


void load_sqw_ext_plugins()
{
	static bool bPluginsLoaded = 0;
	if(!bPluginsLoaded)
	{
		tl::log_info("Loading external plugins from directory: ", g_strApp, ".");

		std::vector<std::string> vecPlugins = tl::get_all_files(g_strApp.c_str());
		for(const std::string& strPlugin : vecPlugins)
		{
			std::string strPluginNoDir = tl::get_file_nodir<std::string>(strPlugin);

			// file names have to start with "takinmod_"
			if(!tl::begins_with<std::string>(strPluginNoDir, "takinmod_", false))
				continue;


			// get module infos
			tl::PipeProc<char> proc(("\"" + strPlugin + "\"").c_str(), false);
			if(!proc.IsReady())
			{
				tl::log_err("Cannot query external plugin infos for \"", strPlugin, "\".");
				continue;
			}

			std::string strModIdent, strModName, strTakVer, strModHelp;
			bool in_help_text = false;
			// get module descriptor strings
			while(!proc.GetIstr().eof())
			{
				std::string line;
				std::getline(proc.GetIstr(), line);

				if(in_help_text)
				{
					if(line == "module_help: end" || line == "module_help:end")
					{
						// end of help text
						in_help_text = false;
					}
					else
					{
						// add line of help text
						if(!strModHelp.empty())
							strModHelp += "\n";
						strModHelp += line;
					}

					continue;
				}

				std::vector<std::string> vecTokens;
				tl::get_tokens<std::string, std::string>(line, std::string(":"), vecTokens);
				if(vecTokens.size() != 2)
					continue;

				if(vecTokens[0] == "module_ident")
					strModIdent = tl::trimmed(vecTokens[1]);
				else if(vecTokens[0] == "module_name")
					strModName = tl::trimmed(vecTokens[1]);
				else if(vecTokens[0] == "required_takin_version")
					strTakVer = tl::trimmed(vecTokens[1]);
				else if(vecTokens[0] == "module_help" && tl::trimmed(vecTokens[1]) == "begin")
					in_help_text = true;
				else if(vecTokens[0] == "module_help" && tl::trimmed(vecTokens[1]) == "end")
					in_help_text = false;
			}

			if(strTakVer == "")
			{
				tl::log_err("Skipping S(Q, E) plugin \"", strPlugin,
					"\" as it is not responding.");
				continue;
			}
			if(strTakVer != TAKIN_VER)
			{
				tl::log_err("Skipping external S(Q, E) plugin \"", strPlugin,
					"\" as it was compiled for Takin version ", strTakVer,
					", but this is version ", TAKIN_VER, ".");
				continue;
			}

			g_vecExtMods.push_back(strModIdent);
			g_mapSqwExt.insert(std::make_pair(strModIdent, std::make_tuple(strModName, strPlugin, strModHelp)));
			tl::log_info("Loaded plugin: ", strPlugin, " -> ", strModIdent, " (\"", strModName, "\").");
		}

		tl::log_debug("Loaded all exernal plugins.");
		bPluginsLoaded = 1;
	}
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// plugins
// --------------------------------------------------------------------------------
#ifdef USE_PLUGINS

#include <boost/dll/shared_library.hpp>
#include <boost/dll/import.hpp>

namespace so = boost::dll;


// tracking modules for refcounting
static std::vector<std::shared_ptr<so::shared_library>> g_vecMods;


void unload_sqw_plugins()
{
	for(auto& pMod : g_vecMods)
	{
		if(!pMod)
			continue;

		std::function<t_fkt_info> fktInfo =
#ifndef __MINGW32__
			pMod->get<t_pfkt_info>("takin_sqw_info");
#else
			pMod->get<t_fkt_info>("takin_sqw_info");
#endif

		if(fktInfo)
		{
			auto tupInfo = fktInfo();
			const std::string& strModIdent = std::get<1>(tupInfo);
			tl::log_debug("Unloading plugin: ", strModIdent, ".");

			// remove module from map
			auto iterMod = g_mapSqw.find(strModIdent);
			if(iterMod != g_mapSqw.end())
				g_mapSqw.erase(iterMod);

			// also look in the map for modules using the raw interface
			auto iterModRaw = g_mapSqwRaw.find(strModIdent);
			if(iterModRaw != g_mapSqwRaw.end())
				g_mapSqwRaw.erase(iterModRaw);
		}

		pMod->unload();
		pMod.reset();
	}

	g_vecMods.clear();
	tl::log_debug("Unloaded all plugins.");

	// also unload the external process plugins
	unload_sqw_ext_plugins();
}


void load_sqw_plugins()
{
	static bool bPluginsLoaded = 0;
	if(!bPluginsLoaded)
	{
		// look in the directories "plugins" and "takin_plugins"
		std::vector<std::string> vecPlugins = find_resource_dirs("plugins", false);
		for(const std::string& plugin : find_resource_dirs("takin_plugins", false))
			vecPlugins.push_back(plugin);

		for(const std::string& strPlugins : vecPlugins)
		{
			tl::log_info("Loading plugins from directory: ", strPlugins, ".");

			std::vector<std::string> vecPlugins = tl::get_all_files(strPlugins.c_str());
			for(const std::string& strPlugin : vecPlugins)
			{
				try
				{
					// TODO: libjulia.so needs rtld_global, but cannot be used here as the takin_sqw_info functions are named the same in all so files...
					std::shared_ptr<so::shared_library> pmod =
						std::make_shared<so::shared_library>(strPlugin,
							so::load_mode::rtld_lazy | so::load_mode::rtld_local);
					if(!pmod || !*pmod)
						continue;

					// import info function
					if(!pmod->has("takin_sqw_info"))
					{
						tl::log_err(strPlugin, " has no takin_sqw_info function.");
						continue;
					}
					std::function<t_fkt_info> fktInfo =
#ifndef __MINGW32__
						pmod->get<t_pfkt_info>("takin_sqw_info");
#else
						pmod->get<t_fkt_info>("takin_sqw_info");
#endif
					if(!fktInfo)
					{
						tl::log_err(strPlugin, " has no valid takin_sqw_info function.");
						continue;
					}

					auto tupInfo = fktInfo();
					const std::string& strTakVer = std::get<0>(tupInfo);
					const std::string& strModIdent = std::get<1>(tupInfo);
					const std::string& strModLongName = std::get<2>(tupInfo);
					const std::string& strModHelp = std::get<3>(tupInfo);

					// module already registered?
					if(g_mapSqw.find(strModIdent) != g_mapSqw.end())
					{
						tl::log_warn("Module \"", strModLongName, "\" (id=", strModIdent, ") is already registered."
							" Plugin: ", strPlugin, ".");
						pmod->unload();
						continue;
					}
					if(strTakVer == "")
					{
						tl::log_err("Skipping S(Q, E) plugin \"", strPlugin,
							"\" as it is not responding.");
						pmod->unload();
						continue;
					}
					if(strTakVer != TAKIN_VER)
					{
						tl::log_err("Skipping S(Q, E) plugin \"", strPlugin,
							"\" as it was compiled for Takin version ", strTakVer,
							", but this is version ", TAKIN_VER, ".");
						pmod->unload();
						continue;
					}


					// import factory function
					if(pmod->has("takin_sqw_new") && pmod->has("takin_sqw_del"))
					{
#ifndef __MINGW32__
						t_pfkt_raw_new pFktNew = pmod->get<t_pfkt_raw_new>("takin_sqw_new");
						t_pfkt_raw_del pFktDel = pmod->get<t_pfkt_raw_del>("takin_sqw_del");
#else
						t_pfkt_raw_new pFktNew = pmod->get<t_fkt_raw_new>("takin_sqw_new");
						t_pfkt_raw_del pFktDel = pmod->get<t_fkt_raw_del>("takin_sqw_del");
#endif
						if(!pFktNew || !pFktDel)
						{
							pmod->unload();
							continue;
						}

						// use the raw new/delete interface if it exists
						g_mapSqwRaw.insert( t_mapSqwRaw::value_type
						{
							strModIdent,
							t_mapSqwRaw::mapped_type
								{ pFktNew, pFktDel, strModLongName, strModHelp }
						});
					}
					else if(pmod->has("takin_sqw"))
					{
						// if raw interface does not exist, try the old shared_ptr one
#ifndef __MINGW32__
						t_pfkt pFkt = pmod->get<t_pfkt>("takin_sqw");
#else
						t_pfkt pFkt = pmod->get<t_fkt>("takin_sqw");
#endif
						if(!pFkt)
						{
							pmod->unload();
							continue;
						}

						g_mapSqw.insert( t_mapSqw::value_type
						{
							strModIdent,
							t_mapSqw::mapped_type
								{ pFkt, strModLongName, strModHelp }
						});
					}
					else
					{
						tl::log_err("No valid constructor interface found in \"", strPlugin, "\".");
						continue;
					}


					g_vecMods.emplace_back(std::move(pmod));
					tl::log_info("Loaded plugin: ", strPlugin,
						" -> ", strModIdent, " (\"", strModLongName, "\").");
				}
				catch(const std::exception& ex)
				{
					tl::log_err("Could not load ", strPlugin, ". Reason: ", ex.what());
				}
			}
		}

		tl::log_debug("Loaded all plugins.");
		bPluginsLoaded = 1;
	}


	// also load the external process plugins
	load_sqw_ext_plugins();
}
// --------------------------------------------------------------------------------


#else


void unload_sqw_plugins()
{
	unload_sqw_ext_plugins();
}

void load_sqw_plugins()
{
	// only load the external process plugins, skip the SO files
	load_sqw_ext_plugins();
}

#endif



// ----------------------------------------------------------------------------
// saving and loading of parameters

bool save_sqw_params(const SqwBase* pSqw,
    std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	if(!pSqw) return 0;

	const std::vector<SqwBase::t_var> vecVars = pSqw->GetVars();
	const std::vector<SqwBase::t_var_fit>& vecFitVars = pSqw->GetFitVars();

	for(const SqwBase::t_var& var : vecVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strVal = std::get<2>(var);

		mapConf[strXmlRoot + "sqw_params/" + strVar] = strVal;
	}

	for(const SqwBase::t_var_fit& var : vecFitVars)
	{
		const std::string& strVar = std::get<0>(var);
		const std::string& strErr = std::get<1>(var);
		const std::string& strRange = std::get<3>(var);
		const bool bFit = std::get<2>(var);

		mapConf[strXmlRoot + "sqw_errors/" + strVar] = strErr;
		mapConf[strXmlRoot + "sqw_fitvar/" + strVar] = bFit ? "1" : "0";
		mapConf[strXmlRoot + "sqw_ranges/" + strVar] = strRange;
	}

	return 1;
}


bool load_sqw_params(SqwBase* pSqw,
	const tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	if(!pSqw) return 0;

	std::vector<std::string> vecChildren =
		xml.GetChildNodes(strXmlRoot + "sqw_params/");

	std::vector<SqwBase::t_var> vecVars;
	std::vector<SqwBase::t_var_fit> vecVarsFit;
	vecVars.reserve(vecChildren.size());
	vecVarsFit.reserve(vecChildren.size());

	for(const std::string& strChild : vecChildren)
	{
		boost::optional<std::string> opVal =
			xml.QueryOpt<std::string>(strXmlRoot + "sqw_params/" + strChild);
		if(opVal)
		{
			SqwBase::t_var var;
			std::get<0>(var) = strChild;
			std::get<2>(var) = *opVal;
			vecVars.emplace_back(std::move(var));
		}

		boost::optional<std::string> opErr =
			xml.QueryOpt<std::string>(strXmlRoot + "sqw_errors/" + strChild);
		boost::optional<std::string> opRange =
			xml.QueryOpt<std::string>(strXmlRoot + "sqw_ranges/" + strChild);
		boost::optional<bool> opFit =
			xml.QueryOpt<bool>(strXmlRoot + "sqw_fitvar/" + strChild);
		SqwBase::t_var_fit varFit;
		std::get<0>(varFit) = strChild;
		std::get<1>(varFit) = opErr ? *opErr : "0";
		std::get<3>(varFit) = opRange ? *opRange : "open : open";
		std::get<2>(varFit) = opFit ? *opFit : 0;
		vecVarsFit.emplace_back(std::move(varFit));
	}

	pSqw->SetVars(vecVars);
	pSqw->InitFitVars(vecVarsFit);

	return 1;
}
// ----------------------------------------------------------------------------
