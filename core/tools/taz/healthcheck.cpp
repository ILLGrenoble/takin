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
#include "tlibs/log/log.h"
#include "tlibs/file/prop.h"

#include <sstream>
#include <string>
#include <cstdlib>


bool healthcheck()
{
	bool ok = true;

	// test external tools
	std::string strTools = find_resource("res/conf/tools.xml");
	tl::Prop<std::string> propTools;
	if(strTools != "" && propTools.Load(strTools.c_str(), tl::PropType::XML))
	{
		// to avoid two separators in a row
		bool bJustAddedSeparator = false;
		// add all menu entries
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
	}
	else
	{
		tl::log_err("Cannot load tool configuration file \"", strTools, "\".");
		ok = false;
	}

	return ok;
}
