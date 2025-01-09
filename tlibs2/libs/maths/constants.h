/**
 * tlibs2 maths library -- constants
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

#ifndef __TLIBS2_MATHS_CONSTS_H__
#define __TLIBS2_MATHS_CONSTS_H__

#include <cmath>

#ifdef __TLIBS2_USE_NUMBERS__
	#include <numbers>
#endif

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

#ifdef __TLIBS2_USE_NUMBERS__
	template<typename T = double> constexpr T golden{std::numbers::phi_v<T>};
	template<typename T = double> constexpr T pi{std::numbers::pi_v<T>};
#else
	// see: https://en.wikipedia.org/wiki/Golden_ratio
	template<typename T = double> constexpr T golden{1.618033988749895};
	template<typename T = double> constexpr T pi{M_PI};
#endif

template<typename INT = int> bool is_even(INT i) { return (i%2 == 0); }
template<typename INT = int> bool is_odd(INT i) { return !is_even<INT>(i); }

template<class T = double> constexpr T r2d(T rad) { return rad/pi<T>*T(180); }     // rad -> deg
template<class T = double> constexpr T d2r(T deg) { return deg/T(180)*pi<T>; }     // deg -> rad
template<class T = double> constexpr T r2m(T rad) { return rad/pi<T>*T(180*60); }  // rad -> min
template<class T = double> constexpr T m2r(T min) { return min/T(180*60)*pi<T>; }  // min -> rad

/**
 * Gaussian around 0: f(x) = exp(-1/2 * (x/sig)^2)
 * at hwhm: f(x_hwhm) = 1/2
 *          exp(-1/2 * (x_hwhm/sig)^2) = 1/2
 *          -1/2 * (x_hwhm/sig)^2 = ln(1/2)
 *          (x_hwhm/sig)^2 = -2*ln(1/2)
 *          x_hwhm^2 = sig^2 * 2*ln(2)
 */
template<class T = double> static constexpr T SIGMA2FWHM = T(2)*std::sqrt(T(2)*std::log(T(2)));
template<class T = double> static constexpr T SIGMA2HWHM = std::sqrt(T(2)*std::log(T(2)));
template<class T = double> static constexpr T FWHM2SIGMA = T(1)/SIGMA2FWHM<T>;
template<class T = double> static constexpr T HWHM2SIGMA = T(1)/SIGMA2HWHM<T>;

}

#endif
