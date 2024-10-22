/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#ifndef __TLIBS2_MAGDYN_HELPERS_H__
#define __TLIBS2_MAGDYN_HELPERS_H__

#include <string>

#include "../maths.h"



namespace tl2_mag {

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

/**
 * rotate spin vector for incommensurate structures, i.e. helices
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
void rotate_spin_incommensurate(t_vec& spin_vec,
	const t_vec& sc_vec, const t_vec& ordering, const t_vec& rotaxis,
	t_real eps = std::numeric_limits<t_real>::epsilon())
{
	const t_real sc_angle = t_real(2) * tl2::pi<t_real> *
		tl2::inner<t_vec>(ordering, sc_vec);

	if(!tl2::equals_0<t_real>(sc_angle, t_real(eps)))
	{
		t_mat sc_rot = tl2::rotation<t_mat, t_vec>(rotaxis, sc_angle);
		spin_vec = sc_rot * spin_vec;
	}
}




/**
 * polarisation matrix for the helical case
 * @see https://doi.org/10.1088/1361-6463/aa7573
 */
template<class t_mat, class t_cplx = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
t_mat get_polarisation_incommensurate(int channel = 0, bool in_chiral_basis = true)
{
	if(in_chiral_basis)
	{
		t_mat pol = tl2::zero<t_mat>(3);

		// just pick the selected component on the diagonal
		if(channel >=0 && channel < 3)
			pol(channel, channel) = t_cplx(1);

		return pol;
	}
	else
	{
		constexpr const t_cplx halfi = t_cplx(0, 0.5);
		constexpr const t_cplx half = t_cplx(0.5, 0);

		// TODO: check coordinate system
		switch(channel)
		{
			case 0: return tl2::create<t_mat>({
				{   half, +halfi,  0 },
				{ -halfi,   half,  0 },
				{      0,      0,  0 } });
			case 1: return tl2::create<t_mat>({
				{   half, -halfi,  0 },
				{ +halfi,   half,  0 },
				{      0,      0,  0 } });
			case 2: return tl2::create<t_mat>({
				{      0,      0,  0 },
				{      0,      0,  0 },
				{      0,      0,  1 } });
		}
	}

	return tl2::zero<t_mat>(3);
}



/**
 * create a 3-vector from a homogeneous 4-vector
 */
template<class t_vec>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec>
#endif
t_vec to_3vec(const t_vec& vec)
{
	return tl2::create<t_vec>({ vec[0], vec[1], vec[2] });
}



/**
 * create a (homogeneous) 4-vector from a 3-vector
 */
template<class t_vec, class t_val = typename t_vec::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec>
#endif
t_vec to_4vec(const t_vec& vec, t_val w = 0.)
{
	return tl2::create<t_vec>({ vec[0], vec[1], vec[2], w });
}
// ----------------------------------------------------------------------------


}
#endif
