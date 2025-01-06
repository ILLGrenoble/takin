/**
 * magnetic dynamics -- type definitions and setting variables
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

#ifndef DONT_USE_QT
	#include <QtCore/QSettings>
#endif

#include <thread>

#include "defs.h"



// ----------------------------------------------------------------------------
// global settings variables
// ----------------------------------------------------------------------------
// number of threads for calculation
unsigned int g_num_threads = std::max<unsigned int>(
	1, std::thread::hardware_concurrency()/2);

// maximum number of recent files
unsigned int g_maxnum_recents = 16;

// epsilons and precisions
int g_prec = 6;
int g_prec_gui = 3;
t_real g_eps = 1e-6;

// delta for numerical differentiation
t_real g_delta_diff = 1e-12;

// bose cutoff energy
t_real g_bose_cutoff = 0.025;

// settings for cholesky decomposition
unsigned int g_cholesky_maxtries = 50;
t_real g_cholesky_delta = 0.0025;

// optional features
int g_allow_ortho_spin = 0;
int g_allow_general_J = 1;
int g_evecs_ortho = 1;

// console messages
int g_silent = 1;
int g_checks = 0;

#ifndef DONT_USE_QT
// gui theme
QString g_theme = "Fusion";

// gui font
QString g_font = "";
#endif

// use native menu bar?
int g_use_native_menubar = 0;

// use native dialogs?
int g_use_native_dialogs = 0;

// plot colour
std::string g_colPlot = "#0000ff";
std::string g_colPlotDegen = "#ff0000";

// structure plotter settings
t_real g_structplot_site_rad = 0.05;
t_real g_structplot_term_rad = 0.01;
t_real g_structplot_dmi_rad = 0.015;
t_real g_structplot_dmi_len = 0.25;
t_real g_structplot_fov = 90.;
// ----------------------------------------------------------------------------



#ifndef DONT_USE_QT
/**
 * transfer the setting from the takin core program
 */
void get_settings_from_takin_core()
{
	QSettings sett_core("takin", "core");

	if(sett_core.contains("main/max_threads"))
	{
		g_num_threads = sett_core.value("main/max_threads").toUInt();
	}

	if(sett_core.contains("main/font_gen"))
	{
		g_font = sett_core.value("main/font_gen").toString();
	}

	if(sett_core.contains("main/prec"))
	{
		g_prec = sett_core.value("main/prec").toInt();
		g_eps = std::pow(t_real(10), -t_real(g_prec));
	}

	if(sett_core.contains("main/prec_gfx"))
	{
		g_prec_gui = sett_core.value("main/prec_gfx").toInt();
	}

	if(sett_core.contains("main/gui_style_value"))
	{
		g_theme = sett_core.value("main/gui_style_value").toString();
	}
}
#endif
