/**
 * simple resolution calculation including only ki and kf errors
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2016
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

#ifndef __SIMPLERESO_H__
#define __SIMPLERESO_H__

#include "../defs.h"
#include "tlibs/phys/neutrons.h"


namespace units = boost::units;
namespace codata = boost::units::si::constants::codata;


struct SimpleResoParams
{
	// values
	tl::t_wavenumber_si<t_real_reso> ki, kf, Q;
	tl::t_energy_si<t_real_reso> E;

	tl::t_angle_si<t_real_reso> twotheta,
		angle_ki_Q, angle_kf_Q;

	// sigmas
	tl::t_wavenumber_si<t_real_reso> sig_ki, sig_kf;
	tl::t_wavenumber_si<t_real_reso> sig_ki_perp, sig_kf_perp;
	tl::t_wavenumber_si<t_real_reso> sig_ki_z, sig_kf_z;
};


extern ResoResults calc_simplereso(const SimpleResoParams& params);

#endif
