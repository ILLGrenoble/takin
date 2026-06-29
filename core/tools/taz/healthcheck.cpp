/**
 * program integrity check
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date mar-2026
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

#include "libs/globals.h"
#include "libs/version.h"
#include "tools/monteconvo/sqwfactory_defs.h"

#include "tlibs/log/log.h"
#include "tlibs/file/prop.h"
#include "tlibs/file/file.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/proc.h"

#include <unordered_set>
#include <sstream>
#include <string>
#include <cstdlib>


static bool check_tools()
{
	// test external tools
	std::string strTools = find_resource("res/conf/tools.xml");
	tl::Prop<std::string> propTools;
	if(strTools == "" || !propTools.Load(strTools.c_str(), tl::PropType::XML))
	{
		tl::log_err("Cannot load tool configuration file \"", strTools, "\".");
		return false;
	}

	// iterate all external tools
	bool ok = true;
	for(std::size_t entry = 0; true; ++entry)
	{
		std::ostringstream _xmlpath;
		_xmlpath << "tools/entry_" << entry;
		std::string xmlpath = _xmlpath.str();
		if(!propTools.Exists(xmlpath))
			break;
		std::string tooltype = propTools.Query<std::string>(xmlpath + "/type", "");
		if(tooltype != "program")
			continue;

		std::string toolname = propTools.Query<std::string>(xmlpath + "/name", "");
		std::string toolprog = propTools.Query<std::string>(xmlpath + "/program", "");
		bool toolmand = propTools.Query<bool>(xmlpath + "/mandatory", false);
		if(toolname == "")
		{
			tl::log_err("Invalid tool name");
			ok = false;
			continue;
		}
		std::string toolbin = find_program_binary(toolprog);
		if(toolbin == "")
		{
			if(toolmand)
			{
				ok = false;
				tl::log_err("Tool binary \"", toolprog, "\" was not found.");
			}
			else
			{
				tl::log_warn("Tool binary \"", toolprog, "\" was not found.");
			}
			continue;
		}

		// run exernal tool process
		tl::log_debug("Running integrity check for process \"", toolbin, "\"...");
		if(std::system(("\"" + toolbin + "\" --healthcheck").c_str()) != 0)
		{
			tl::log_err("Integrity check for process \"", toolbin, "\" failed.");
			ok = false;
		}
	}

	return ok;
}


static bool check_ext_plugins()
{
	bool ok = true;
	bool pymod_seen = false;
	bool jlmod_seen = false;

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
			ok = false;
			continue;
		}

		std::string strModIdent, strTakVer;
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
				continue;
			}

			std::vector<std::string> vecTokens;
			tl::get_tokens<std::string, std::string>(line, std::string(":"), vecTokens);
			if(vecTokens.size() != 2)
				continue;

			if(vecTokens[0] == "module_ident")
				strModIdent = tl::trimmed(vecTokens[1]);
			else if(vecTokens[0] == "required_takin_version")
				strTakVer = tl::trimmed(vecTokens[1]);
			else if(vecTokens[0] == "module_help" && tl::trimmed(vecTokens[1]) == "begin")
				in_help_text = true;
			else if(vecTokens[0] == "module_help" && tl::trimmed(vecTokens[1]) == "end")
				in_help_text = false;
		}

		if(strTakVer == "")
		{
			tl::log_err("External S(Q, E) plugin \"", strPlugin, "\" is not responding.");
			ok = false;
			continue;
		}

		if(strTakVer != TAKIN_VER)
		{
			tl::log_err("External S(Q, E) plugin \"", strPlugin,
				"\" was compiled for Takin version ", strTakVer,
				", but this is version ", TAKIN_VER, ".");
			ok = false;
			continue;
		}

		if(strPluginNoDir == "takinmod_py" && strModIdent == "py")
			pymod_seen = true;
		if(strPluginNoDir == "takinmod_jl" && strModIdent == "jl0")
			jlmod_seen = true;
	}

	if(!pymod_seen)
	{
		tl::log_err("Py scripting module, takinmod_py, doesn't seem to be available or working.");
		ok = false;
	}

	if(!jlmod_seen)
	{
		tl::log_warn("Jl scripting module, takinmod_jl, doesn't seem to be available or working.");
		//ok = false;
	}

	return ok;
}


#ifdef USE_PLUGINS

#include <boost/dll/shared_library.hpp>
#include <boost/dll/import.hpp>

namespace so = boost::dll;


static bool check_plugins()
{
	bool ok = true;
	bool seen_magnonmod = false;

	// look in the directories "plugins" and "takin_plugins"
	std::vector<std::string> vecPlugins = find_resource_dirs("plugins", false);
	for(const std::string& plugin : find_resource_dirs("takin_plugins", false))
		vecPlugins.push_back(plugin);

	std::unordered_set<std::string> seen_modules;

	for(const std::string& strPlugins : vecPlugins)
	{
		std::vector<std::string> vecPlugins = tl::get_all_files(strPlugins.c_str());
		for(const std::string& strPlugin : vecPlugins)
		{
			try
			{
				std::shared_ptr<so::shared_library> pmod =
					std::make_shared<so::shared_library>(strPlugin,
						so::load_mode::rtld_lazy | so::load_mode::rtld_local);
				if(!pmod || !*pmod)
					continue;

				// import info function
				if(!pmod->has("takin_sqw_info"))
					continue;

				std::function<t_fkt_info> fktInfo =
#ifndef __MINGW32__
					pmod->get<t_pfkt_info>("takin_sqw_info");
#else
					pmod->get<t_fkt_info>("takin_sqw_info");
#endif
				if(!fktInfo)
				{
					pmod->unload();
					continue;
				}

				auto tupInfo = fktInfo();
				const std::string& strTakVer = std::get<0>(tupInfo);
				const std::string& strModIdent = std::get<1>(tupInfo);
				const std::string& strModLongName = std::get<2>(tupInfo);

				// module already registered?
				if(seen_modules.find(strModIdent) != seen_modules.end())
				{
					tl::log_warn("Module \"", strModLongName, "\" (id=", strModIdent, ") is already registered."
						" Plugin: ", strPlugin, ".");
					pmod->unload();
					continue;
				}
				if(strTakVer == "")
				{
					tl::log_err("S(Q, E) plugin \"", strPlugin, "\" is not responding.");
					pmod->unload();
					ok = false;
					continue;
				}
				if(strTakVer != TAKIN_VER)
				{
					tl::log_err("S(Q, E) plugin \"", strPlugin,
						"\" was compiled for Takin version ", strTakVer,
						", but this is version ", TAKIN_VER, ".");
					pmod->unload();
					ok = false;
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
					seen_modules.insert(strModIdent);
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

					seen_modules.insert(strModIdent);
					pmod->unload();
				}
				else
				{
					tl::log_err("No valid constructor interface found in S(Q, E) plugin \"", strPlugin, "\".");
					ok = false;
					continue;
				}

				if(strModIdent == "magnonmod")
					seen_magnonmod = true;
			}
			catch(const std::exception& ex)
			{
				tl::log_warn("Could not load S(Q, E) plugin \"", strPlugin, "\". Reason: ", ex.what());
				//ok = false;
			}
		}
	}

	if(!seen_magnonmod)
	{
		tl::log_err("Magnon module doesn't seem to be available or working.");
		ok = false;
	}

	return ok;
}

#else

static bool check_plugins()
{
	return true;
}

#endif


bool healthcheck()
{
	bool tools_ok = check_tools();
	bool ext_plugins_ok = check_ext_plugins();
	bool plugins_ok = check_plugins();

	return tools_ok && ext_plugins_ok && plugins_ok;
}
