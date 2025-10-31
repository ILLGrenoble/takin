/**
 * implementation of the Violini TOF reso algorithm
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2016
 * @license GPLv2
 *
 * @desc for algo, see: [vio14] N. Violini et al., NIM A 736 (2014) pp. 31-39, doi: 10.1016/j.nima.2013.10.042
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

#ifndef __TOFRESO_H__
#define __TOFRESO_H__

#include "../defs.h"
#include "tlibs/phys/neutrons.h"


namespace units = boost::units;
namespace codata = boost::units::si::constants::codata;


enum class TofDetShape { SPH, CYL, UNKNOWN };


/**
 * TOF parameters in sigma
 */
struct VioParams
{
	// scattering triangle
	tl::t_wavenumber_si<t_real_reso> ki, kf, Q;
	tl::t_energy_si<t_real_reso> E;

	tl::t_angle_si<t_real_reso> twotheta,
		angle_ki_Q, angle_kf_Q;

	tl::t_angle_si<t_real_reso> angle_outplane_i, angle_outplane_f;
	tl::t_angle_si<t_real_reso> twotheta_i;


	// instrument lengths
	tl::t_length_si<t_real_reso> len_pulse_mono,
		len_mono_sample, len_sample_det;


	// instrument sigmas
	tl::t_length_si<t_real_reso> sig_len_pulse_mono,
		sig_len_mono_sample, sig_len_sample_det;

	tl::t_time_si<t_real_reso> sig_pulse, sig_mono, sig_det;

	tl::t_angle_si<t_real_reso> sig_twotheta_f, sig_outplane_f,
		sig_twotheta_i, sig_outplane_i;


	TofDetShape det_shape = TofDetShape::SPH;
};


extern ResoResults calc_vio(const VioParams& params);


#endif
