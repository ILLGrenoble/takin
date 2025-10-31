/**
 * re-implementation of V. Mecoli's extension of Violini's TOF reso algorithm [mec25]
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2025
 * @license GPLv2
 *
 * @desc for extended method, see: [mec25] V. Mecoli, PhD thesis in preparation
 * @desc for original method, see: [vio14] N. Violini et al., NIM A 736 (2014) pp. 31-39, doi: 10.1016/j.nima.2013.10.042
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TOFRESO_EXT_H__
#define __TOFRESO_EXT_H__

#include "vio.h"


/**
 * TOF parameters in sigma
 */
struct VioExtParams : public VioParams
{
	// instrument lengths
	tl::t_length_si<t_real_reso> len_ch_pulse_guide,
		len_ch_mono_guide, len_guide_sample;
	tl::t_length_si<t_real_reso> len_sample_det2;  // TODO: unite with len_sample_det

	// pulse chopper
	tl::t_angle_si<t_real_reso> ch_pulse_angle_win, ch_pulse_angle_beam;
	tl::t_length_si<t_real_reso> ch_pulse_width;
	t_real_reso ch_pulse_rpm;
	bool ch_pulse_counterrot;

	// monochromatising chopper
	tl::t_angle_si<t_real_reso> ch_mono_angle_win, ch_mono_angle_beam;
	tl::t_length_si<t_real_reso> ch_mono_width;
	t_real_reso ch_mono_rpm;
	bool ch_mono_counterrot;

	// guide
	tl::t_length_si<t_real_reso> endguide_width, endguide_height;

	// detector
	tl::t_length_si<t_real_reso> det_tube_width, det_height, det_z;

	// sample
	tl::t_length_si<t_real_reso> sample_width, sample_height;

	// mc points for length calculations
	unsigned int mc_lengths;
};


extern ResoResults calc_vio_ext(const VioExtParams& params);


#endif
