/**
 * implementation of the eckold-sobolev algo
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
 * @license GPLv2
 *
 * @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
 * @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
 * @desc for vertical scattering modification: [eck20] G. Eckold, personal communication, 2020.
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

#ifndef __TAKIN_ECK_H__
#define __TAKIN_ECK_H__

#include "pop.h"

/**
 * TAS parameters in fwhm
 */
struct EckParams : public PopParams
{
	tl::t_length_si<t_real_reso> pos_x, pos_y, pos_z;

	// vertical scattering in k_f: angle_kf = 90 deg
	tl::t_angle_si<t_real_reso> angle_kf;
};


extern ResoResults calc_eck(const EckParams& eck);


/**
 * general R0 normalisation factor from [mit84], equ. A.57
 */
template<class t_real = double>
t_real mitch_R0(bool norm_to_ki_vol,
	t_real dmono_refl, t_real dana_effic,
	t_real dKiVol, t_real dKfVol, t_real dResVol,
	bool bNormToResVol = false)
{
	t_real dR0 = dana_effic * dKfVol;
	if(!norm_to_ki_vol)
		dR0 *= dmono_refl * dKiVol;

	// not needed for MC simulations, because the gaussian generated
	// with std::normal_distribution is already normalised
	// see: tools/test/tst_norm.cpp
	if(bNormToResVol)
		dR0 /= (dResVol * tl::get_pi<t_real>() * t_real{3});

	return dR0;
}


#endif
