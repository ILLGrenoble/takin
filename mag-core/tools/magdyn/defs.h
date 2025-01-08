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

#ifndef __MAGDYN_DEFS__
#define __MAGDYN_DEFS__

#ifndef DONT_USE_QT
	#include <QtCore/QString>
#endif

#include <string>
#include <array>
#include <variant>
#include <optional>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/magdyn.h"
#include "tlibs2/libs/qt/gl.h"

#include "libs/defs.h"



// ----------------------------------------------------------------------------
// type definitions
// ----------------------------------------------------------------------------
using t_vec_real = tl2::vec<t_real, std::vector>;
using t_mat_real = tl2::mat<t_real, std::vector>;

using t_vec = tl2::vec<t_cplx, std::vector>;
using t_mat = tl2::mat<t_cplx, std::vector>;

using t_real_gl = tl2::t_real_gl;
using t_vec2_gl = tl2::t_vec2_gl;
using t_vec3_gl = tl2::t_vec3_gl;
using t_vec_gl = tl2::t_vec_gl;
using t_mat_gl = tl2::t_mat_gl;

// magnon calculation core
using t_magdyn = tl2_mag::MagDyn<t_mat, t_vec, t_mat_real, t_vec_real, t_cplx, t_real, t_size>;
using t_site = typename t_magdyn::MagneticSite;
using t_term = typename t_magdyn::ExchangeTerm;;
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// global settings variables
// ----------------------------------------------------------------------------
// number of threads for calculation
extern unsigned int g_num_threads;

// maximum number of recent files
extern unsigned int g_maxnum_recents;

// number precisions
extern int g_prec, g_prec_gui;

// epsilon
extern t_real g_eps;

// delta for numerical differentiation
extern t_real g_delta_diff;

// bose cutoff energy
extern t_real g_bose_cutoff;

// settings for cholesky decomposition
extern unsigned int g_cholesky_maxtries;
extern t_real g_cholesky_delta;

// optional features
extern int g_allow_ortho_spin, g_allow_general_J;
extern int g_evecs_ortho;

// console messages
extern int g_silent, g_checks;

// use native menubar and dialogs?
extern int g_use_native_menubar, g_use_native_dialogs;

// plot colour
extern std::string g_colPlot, g_colPlotDegen;

// structure plotter settings
extern t_real g_structplot_site_rad;
extern t_real g_structplot_term_rad;
extern t_real g_structplot_dmi_rad;
extern t_real g_structplot_dmi_len;
extern t_real g_structplot_fov;

// fraction of points to check if the stop button was pressed
extern unsigned int g_stop_check_fraction;


#ifndef DONT_USE_QT
	// gui theme and font
	extern QString g_theme, g_font;

	// transfer the setting from the takin core program
	void get_settings_from_takin_core();
#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// variables register
// ----------------------------------------------------------------------------
#include "dialogs/settings.h"

constexpr std::array<SettingsVariable, 22> g_settingsvariables
{{
	// threads
	{
		.description = "Number of threads for calculation.",
		.key = "num_threads",
		.value = &g_num_threads,
	},

	// epsilons and precisions
	{
		.description = "Calculation epsilon.",
		.key = "eps",
		.value = &g_eps,
	},
	{
		.description = "Number precision.",
		.key = "prec",
		.value = &g_prec,
	},
	{
		.description = "GUI number precision.",
		.key = "prec_gui",
		.value = &g_prec_gui,
	},
	// deltas for numerical calculations
	{
		.description = "Differentiation delta.",
		.key = "delta_diff",
		.value = &g_delta_diff,
	},
	{
		.description = "Bose cutoff energy.",
		.key = "bose_cutoff",
		.value = &g_bose_cutoff,
	},

	// settings for cholesky decomposition
	{
		.description = "Cholesky maximum tries.",
		.key = "cholesky_maxtries",
		.value = &g_cholesky_maxtries,
	},
	{
		.description = "Cholesky delta per trial.",
		.key = "cholesky_delta",
		.value = &g_cholesky_delta,
	},

	// file options
	{
		.description = "Maximum number of recent files.",
		.key = "maxnum_recents",
		.value = &g_maxnum_recents,
	},

	// colours
	{
		.description = "Plot colour.",
		.key = "plot_colour",
		.value = &g_colPlot,
	},
	{
		.description = "Plot colour for degeneracies.",
		.key = "plot_colour_degen",
		.value = &g_colPlotDegen,
	},

	// structure plotter settings
	{
		.description = "Site radius in 3d structure plotter.",
		.key = "structplot_site_radius",
		.value = &g_structplot_site_rad,
	},
	{
		.description = "Coupling radius in 3d structure plotter.",
		.key = "structplot_term_radius",
		.value = &g_structplot_term_rad,
	},
	{
		.description = "DMI vector radius in 3d structure plotter.",
		.key = "structplot_dmi_radius",
		.value = &g_structplot_dmi_rad,
	},
	{
		.description = "DMI vector length in 3d structure plotter.",
		.key = "structplot_dmi_length",
		.value = &g_structplot_dmi_len,
	},
	{
		.description = "Camera field-of-view in 3d structure plotter.",
		.key = "structplot_fov",
		.value = &g_structplot_fov,
	},

	// optional features
	{
		.description = "Allow setting of orthogonal spins.",
		.key = "ortho_spins",
		.value = &g_allow_ortho_spin,
		.editor = SettingsVariableEditor::YESNO,
	},
	{
		.description = "Allow setting of general exchange matrix J.",
		.key = "allow_gen_J",
		.value = &g_allow_general_J,
		.editor = SettingsVariableEditor::YESNO,
	},
	{
		.description = "Eigenstates of H are always orthogonal.",
		.key = "evecs_ortho",
		.value = &g_evecs_ortho,
		.editor = SettingsVariableEditor::YESNO,
	},
	{
		.description = "Number of times to check for stop button press.",
		.key = "stop_check_fraction",
		.value = &g_stop_check_fraction,
	},

	// console messages
	{
		.description = "Silence output of error messages on console.",
		.key = "output_silent",
		.value = &g_silent,
		.editor = SettingsVariableEditor::YESNO,
	},
	{
		.description = "Perform extra sanity checks.",
		.key = "sanity_checks",
		.value = &g_checks,
		.editor = SettingsVariableEditor::YESNO,
	},
}};
// ----------------------------------------------------------------------------


#endif
