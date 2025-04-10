/**
 * type definitions for reso
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date mar-2016
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

#ifndef __RESO_DEFS_H__
#define __RESO_DEFS_H__

#include "libs/globals.h"
#include "tlibs/helper/boost_hacks.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace ublas = boost::numeric::ublas;

using t_real_reso = ::t_real_glob;


/**
 * flags for resolution calculation
 */
enum ResoFlags : std::size_t
{
	CALC_KI3        = 1<<0,  // mono efficiency factor
	CALC_KF3        = 1<<1,  // ana efficiency factor

	CALC_KFKI       = 1<<2,  // kf/ki factor
	CALC_MONKI      = 1<<3,  // monitor 1/ki factor
	CALC_MON        = 1<<4,  // monitor in ki axis

	CALC_GENERAL_R0	= 1<<5,  // alternative R0 normalisation

	NORM_TO_SAMPLE  = 1<<6,  // normalise by sample volume
};


struct ResoResults
{
	bool bOk;
	std::string strErr;

	ublas::matrix<t_real_reso> reso;    // quadratic part of quadric
	ublas::vector<t_real_reso> reso_v;  // linear part of quadric
	t_real_reso reso_s;                 // constant part of quadric

	ublas::vector<t_real_reso> Q_avg;
	t_real_reso dR0;                    // resolution prefactor
	t_real_reso dResVol;                // resolution volume in 1/A^3 * meV

	t_real_reso dBraggFWHMs[4];
};


/**
 * all resolution calculation methods
 */
enum class ResoAlgo
{
	CN         = 1,
	POP_CN     = 2,

	POP        = 3,
	ECK        = 4,
	ECK_EXT    = 5,

	VIO        = 6,
	SIMPLE     = 100,

	MC         = 500,

	UNKNOWN    = -1,
};


#endif
