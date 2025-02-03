/**
 * globals
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 20-mar-2015
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

#include "globals.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/proc.h"

#include <thread>

// -----------------------------------------------------------------------------

unsigned int g_iPrec = 6;
unsigned int g_iPrecGfx = 4;

t_real_glob g_dEps = 1e-6;
t_real_glob g_dEpsGfx = 1e-4;

bool g_bHasElements = false;
bool g_bHasFormfacts = false;
bool g_bHasMagFormfacts = false;
bool g_bHasScatlens = false;
bool g_bHasSpaceGroups = false;
bool g_bShowFsq = true;
bool g_b3dBZ = true;
bool g_bUseGlobalPaths = true;

std::string g_strApp;
std::string g_strHome;

std::string g_strGplTool = "gnuplot";

std::size_t GFX_NUM_POINTS = 512;
std::size_t g_iMaxNN = 4;

t_real_glob g_dFontSize = 10.;


// -----------------------------------------------------------------------------

unsigned int g_iMaxThreads = std::thread::hardware_concurrency();
unsigned int g_iMaxProcesses = std::thread::hardware_concurrency();

unsigned int get_max_threads()
{
	unsigned int iMaxThreads = std::thread::hardware_concurrency();
	return std::min(iMaxThreads, g_iMaxThreads);
}

unsigned int get_max_processes()
{
	unsigned int iMaxProcesses = std::thread::hardware_concurrency();
	return std::min(iMaxProcesses, g_iMaxProcesses);
}

// -----------------------------------------------------------------------------


static std::vector<std::string> s_vecInstallPaths =
{
	".",				// local resource dir
#ifdef INSTALL_PREFIX
	INSTALL_PREFIX "/share/takin",	// resource dir from install path
#endif
	"/usr/local/share/takin",	// some default fallback paths
	"/usr/share/takin",
};


void add_resource_path(const std::string& strPath, bool bToBack)
{
	if(bToBack)
		s_vecInstallPaths.push_back(strPath);
	else	// insert at beginning, after "."
		s_vecInstallPaths.insert(s_vecInstallPaths.begin()+1, strPath);
}


std::string find_resource(const std::string& strFile, bool bLogErr)
{
	// relative paths
	for(const std::string& strPrefix : s_vecInstallPaths)
	{
		std::string _strFile = strPrefix + "/" + strFile;
		//tl::log_debug("Looking for file: ", _strFile);
		if(tl::file_exists(_strFile.c_str()))
			return _strFile;
		else if(tl::file_exists((_strFile + ".gz").c_str()))
			return _strFile + ".gz";
		else if(tl::file_exists((_strFile + ".bz2").c_str()))
			return _strFile + ".bz2";
	}

	// absolute path
	if(tl::file_exists(strFile.c_str()))
		return strFile;

	if(bLogErr)
		tl::log_err("Could not find resource file \"", strFile, "\".");
	return "";
}


std::vector<std::string> find_resource_dirs(const std::string& strDir, bool bLogErr)
{
	std::vector<std::string> vecDirs;

	for(const std::string& strPrefix : s_vecInstallPaths)
	{
		std::string _strDir = strPrefix + "/" + strDir;
		if(tl::dir_exists(_strDir.c_str()))
			vecDirs.push_back(_strDir);
	}

	if(bLogErr && vecDirs.size()==0)
		tl::log_err("Could not find resource directory \"", strDir, "\".");

	return vecDirs;
}

// -----------------------------------------------------------------------------


static std::vector<std::string> g_vecPaths;

void add_global_path(const std::string& strPath, bool bToBack)
{
	if(bToBack)
		g_vecPaths.push_back(strPath);
	else
		g_vecPaths.insert(g_vecPaths.begin()+1, strPath);
}


const std::vector<std::string>& get_global_paths()
{
	static const std::vector<std::string> vecEmpty;
	return g_bUseGlobalPaths ? g_vecPaths : vecEmpty;
}


std::string find_file_in_global_paths(const std::string& strFile, bool bAlsoTryFileOnly)
{
	// no file given
	if(strFile == "")
		return "";

	// if the file exists, use it
	if(tl::file_exists(strFile.c_str()))
		return strFile;


	const std::vector<std::string>& vecGlobPaths = get_global_paths();

	// add full path of "strFile" to global paths
	for(const std::string& strGlobPath : vecGlobPaths)
	{
		std::string strNewFile = strGlobPath + "/" + strFile;
		if(tl::file_exists(strNewFile.c_str()))
			return strNewFile;
	}

	if(bAlsoTryFileOnly)
	{
		// add only the file name in "strFile" to global paths
		const std::string strFileOnly = tl::get_file_nodir(strFile);
		for(const std::string& strGlobPath : vecGlobPaths)
		{
			std::string strNewFile = strGlobPath + "/" + strFileOnly;
			if(tl::file_exists(strNewFile.c_str()))
				return strNewFile;
		}
	}

	// nothing found
	return "";
}


void clear_global_paths()
{
	g_vecPaths.clear();
}


// -----------------------------------------------------------------------------



std::string find_program_binary(const std::string& strExe, bool log_messages)
{
	// if the given binary file exists, directly use it
	if(tl::file_exists(strExe.c_str()))
	{
		if(log_messages)
			tl::log_info("Found external tool: \"", strExe, "\".");
		return strExe;
	}

	// try prefixing it with possible application paths
	std::vector<std::string> vecPaths =
	{
		g_strApp + "/" + strExe,
		g_strApp + "/" + strExe + ".exe",
		"/usr/local/bin/" + strExe,
		"/usr/bin/" + strExe,
	};

	for(const std::string& strPath : vecPaths)
	{
		//tl::log_debug("Trying to find ", strPath, "...");
		if(tl::file_exists(strPath.c_str()))
		{
			if(log_messages)
				tl::log_info("Found external tool in program path: \"", strPath, "\".");
			return strPath;
		}
	}

	return "";
}


// -----------------------------------------------------------------------------



std::string get_gpltool_version()
{
	tl::PipeProc<char> proc(("\"" + g_strGplTool + "\" 2>/dev/null --version").c_str(), false);
	if(!proc.IsReady())
		return "";

	std::string strVer;
	std::getline(proc.GetIstr(), strVer);
	tl::trim(strVer);
	return strVer;
}
