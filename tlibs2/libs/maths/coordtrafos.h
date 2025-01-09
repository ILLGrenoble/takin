/**
 * tlibs2 maths library -- coordinate trafos
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
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

#ifndef __TLIBS2_MATHS_COORDTRAFOS_H__
#define __TLIBS2_MATHS_COORDTRAFOS_H__

#include <cmath>
#include <tuple>

#include "decls.h"
#include "constants.h"



namespace tl2 {
// -----------------------------------------------------------------------------
// coordinate trafos
// -----------------------------------------------------------------------------

/**
 * cartesian -> spherical
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T, T, T> cart_to_sph(T x, T y, T z)
{
	T rho = std::sqrt(x*x + y*y + z*z);
	T phi = std::atan2(y, x);      // range: [-pi, pi]
	T theta = std::acos(z / rho);  // range: [0, pi]

	return std::make_tuple(rho, phi, theta);
}


/**
 * spherical -> cartesian
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T, T, T> sph_to_cart(T rho, T phi, T theta)
{
	T x = rho * std::cos(phi)*std::sin(theta);
	T y = rho * std::sin(phi)*std::sin(theta);
	T z = rho * std::cos(theta);

	return std::make_tuple(x, y, z);
}


/**
 * uv parameter -> spherical (for uniform point distribution on sphere)
 * u in [0, 1], v in [0, 1]
 * @see https://mathworld.wolfram.com/SpherePointPicking.html
 */
template<class T = double>
std::tuple<T, T> uv_to_sph(T u, T v)
{
	T phi = T(2)*pi<T>*u - pi<T>;

	// "reflect back" out-of-range parameter
	if(v < 0.)
		v = -v;
	else if(v > 1.)
		v -= v - 1.;

	T c = T(2)*v - T(1);
	c = std::clamp<T>(c, -1., 1.);
	T theta = std::acos(c);

	return std::make_tuple(phi, theta);
}


/**
 * spherical -> uv parameter
 * u in [0, 1], v in [0, 1]
 * @see https://mathworld.wolfram.com/SpherePointPicking.html
 */
template<class T = double>
std::tuple<T, T> sph_to_uv(T phi, T theta)
{
	T u = (phi + pi<T>) / (T(2) * pi<T>);
	T v = (std::cos(theta) + T(1)) / T(2);

	return std::make_tuple(u, v);
}


/**
 * cylindrical -> spherical
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T, T, T> cyl_to_sph(T rho_cyl, T phi_cyl, T z_cyl)
{
	T rho = std::sqrt(rho_cyl*rho_cyl + z_cyl*z_cyl);
	T theta = std::acos(z_cyl / rho);

	return std::make_tuple(rho, phi_cyl, theta);
}


/**
 * spherical -> cylindrical
 * @see https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
template<class T = double>
std::tuple<T, T, T> sph_to_cyl(T rho_sph, T phi_sph, T theta_sph)
{
	T rho = rho_sph * std::sin(theta_sph);
	T z = rho_sph * std::cos(theta_sph);

	return std::make_tuple(rho, phi_sph, z);
}


/**
 * cylindrical -> cartesian
 * @see https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
 */
template<class T = double>
std::tuple<T, T, T> cyl_to_cart(T rho, T phi, T z)
{
	T x = rho * std::cos(phi);
	T y = rho * std::sin(phi);

	return std::make_tuple(x, y, z);
}


/**
 * cartesian -> cylindrical
 * @see https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
 */
template<class T = double>
std::tuple<T, T, T> cart_to_cyl(T x, T y, T z)
{
	T rho = std::sqrt(x*x + y*y);
	T phi = std::atan2(y, x);

	return std::make_tuple(rho, phi, z);
}


/**
 * gnomonic projection (similar to perspective projection with fov=90°)
 * @return [x,y]
 * @see see http://mathworld.wolfram.com/GnomonicProjection.html
 */
template<class T = double>
std::tuple<T, T> gnomonic_proj(T twophi_crys, T twotheta_crys)
{
	T x = -std::tan(twophi_crys);
	T y = std::tan(twotheta_crys) / std::cos(twophi_crys);

	return std::make_tuple(x, y);
}


/**
 * stereographic projection
 * @return [x,y]
 * @see http://mathworld.wolfram.com/StereographicProjection.html
 */
template<class T = double>
std::tuple<T, T> stereographic_proj(T twophi_crys, T twotheta_crys, T rad)
{
	const T sth = std::sin(twotheta_crys);
	const T cth = std::cos(twotheta_crys);
	const T sph = std::sin(twophi_crys);
	const T cph = std::cos(twophi_crys);

	T x = -T(2) * rad * sph * cth / (T(1) + cth*cph);
	T y = T(2) * rad * sth / (T(1) + cth*cph);

	return std::make_tuple(x, y);
}

}

#endif
