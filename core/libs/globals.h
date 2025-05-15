/**
 * globals
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 20-mar-2015
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

#ifndef __TAKIN_GLOBALS_H__
#define __TAKIN_GLOBALS_H__

#include <vector>
#include <string>


//using t_real_glob = float;
using t_real_glob = double;


extern unsigned int g_iPrec;
extern unsigned int g_iPrecGfx;

extern t_real_glob g_dEps;
extern t_real_glob g_dEpsGfx;

extern unsigned int g_iMaxThreads;
extern unsigned int g_iMaxProcesses;

extern std::size_t GFX_NUM_POINTS;
extern std::size_t g_iMaxNN;

extern bool g_bHasElements;
extern bool g_bHasFormfacts;
extern bool g_bHasMagFormfacts;
extern bool g_bHasScatlens;
extern bool g_bHasSpaceGroups;
extern bool g_bShowFsq;
extern bool g_b3dBZ;

extern std::string g_strApp;	// application dir
extern std::string g_strHome;	// home dir

extern std::string g_strGplTool;

extern t_real_glob g_dFontSize;

extern void add_resource_path(const std::string& strPath, bool bToBack = true);
extern std::string find_resource(const std::string& strFile, bool bLogErr = true);
extern std::vector<std::string> find_resource_dirs(const std::string& strDir, bool bLogErr = true);

extern void add_global_path(const std::string& strPath, bool bToBack=1);
extern const std::vector<std::string>& get_global_paths();
extern void clear_global_paths();
extern std::string find_file_in_global_paths(const std::string& strFile, bool bAlsoTryFileOnly = true);

extern std::string find_program_binary(const std::string& strExe, bool log_messages = true);

extern unsigned int get_max_threads();
extern unsigned int get_max_processes();

extern std::string get_gpltool_version();

#endif
