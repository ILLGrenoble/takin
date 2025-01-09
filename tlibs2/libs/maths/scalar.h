/**
 * tlibs2 maths library -- scalar algorithms
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

#ifndef __TLIBS2_MATHS_SCALAR_H__
#define __TLIBS2_MATHS_SCALAR_H__

#include <limits>
#include <vector>
#include <random>
#include <cmath>

#include "decls.h"
#include "constants.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// scalar algorithms
// ----------------------------------------------------------------------------

/**
 * are two scalars equal within an epsilon range?
 */
template<class T>
bool equals(T t1, T t2, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	return std::abs(t1 - t2) <= eps;
}


/**
 * is the given value an integer?
 */
template<class T = double>
bool is_integer(T val, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	T val_diff = val - std::round(val);
	return tl2::equals<T>(val_diff, T(0), eps);
}


/**
 * get next multiple of the given granularity
 */
template<typename t_num = unsigned int>
t_num next_multiple(t_num num, t_num granularity)
requires is_scalar<t_num>
{
	t_num div = num / granularity;
	bool rest_is_0 = 1;

	if constexpr(std::is_floating_point_v<t_num>)
	{
		div = std::floor(div);
		t_num rest = std::fmod(num, granularity);
		rest_is_0 = tl2::equals(rest, t_num{0});
	}
	else
	{
		t_num rest = num % granularity;
		rest_is_0 = (rest==0);
	}

	return rest_is_0 ? num : (div+1) * granularity;
}


/**
 * mod operation, keeping result positive
 */
template<class t_real>
t_real mod_pos(t_real val, t_real tomod = t_real{2}*pi<t_real>)
requires is_scalar<t_real>
{
	val = std::fmod(val, tomod);
	if(val < t_real(0))
		val += tomod;

	return val;
}


/**
 * mod operation, keeping result in range [ -tomod/2, tomod/2 ]
 */
template<class t_real>
t_real mod_posneg(t_real val, t_real tomod = t_real{2}*pi<t_real>)
requires is_scalar<t_real>
{
	t_real halfrange = tomod/t_real{2};
	return mod_pos<t_real>(val + halfrange, tomod) - halfrange;
}


/**
 * are two angles equal within an epsilon range?
 */
template<class T>
bool angle_equals(T t1, T t2,
	T eps = std::numeric_limits<T>::epsilon(),
	T tomod = T{2}*pi<T>)
requires is_scalar<T>
{
	t1 = mod_pos<T>(t1, tomod);
	t2 = mod_pos<T>(t2, tomod);

	return std::abs(t1 - t2) <= eps;
}


/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T> requires is_complex<T>
bool equals(const T& t1, const T& t2,
	typename T::value_type eps = std::numeric_limits<typename T::value_type>::epsilon())
{
	return (std::abs(t1.real() - t2.real()) <= eps) &&
		(std::abs(t1.imag() - t2.imag()) <= eps);
}


/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T> requires is_complex<T>
bool equals(const T& t1, const T& t2, const T& eps)
{
	return tl2::equals<T>(t1, t2, eps.real());
}


template<typename T> T sign(T t)
{
	if(t < 0.)
		return -T(1);
	return T(1);
}


template<typename T> T cot(T t)
{
	return std::tan(T(0.5)*pi<T> - t);
}


template<typename T> T coth(T t)
{
	return T(1) / std::tanh(t);
}


template<typename T=double>
T log(T tbase, T tval)
{
	return T(std::log(tval)/std::log(tbase));
}


template<typename T=double>
T nextpow(T tbase, T tval)
{
	return T(std::pow(tbase, std::ceil(log(tbase, tval))));
}


/**
 * unsigned angle between two vectors
 * <q1|q2> / (|q1| |q2|) = cos(alpha)
 */
template<class t_vec>
typename t_vec::value_type angle_unsigned(const t_vec& q1, const t_vec& q2)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	if(q1.size() != q2.size())
		return t_real(0);

	t_real dot = t_real(0);
	t_real len1 = t_real(0);
	t_real len2 = t_real(0);

	for(std::size_t i = 0; i < q1.size(); ++i)
	{
		dot += q1[i]*q2[i];

		len1 += q1[i]*q1[i];
		len2 += q2[i]*q2[i];
	}

	len1 = std::sqrt(len1);
	len2 = std::sqrt(len2);

	dot /= len1;
	dot /= len2;

	return std::acos(dot);
}


/**
 * x = 0..1, y = 0..1
 * @see https://en.wikipedia.org/wiki/Bilinear_interpolation
 */
template<typename T=double>
T bilinear_interp(T x0y0, T x1y0, T x0y1, T x1y1, T x, T y)
{
	T top = std::lerp(x0y1, x1y1, x);
	T bottom = std::lerp(x0y0, x1y0, x);

	return std::lerp(bottom, top, y);
}


template<typename T=double, typename REAL=double,
template<class...> class t_vec = std::vector>
t_vec<T> linspace(const T& tmin, const T& tmax, std::size_t iNum)
{
	t_vec<T> vec;
	vec.reserve(iNum);

	if(iNum == 1)
	{
		// if just one point is requested, use the lower limit
		vec.push_back(tmin);
		return vec;
	}

	for(std::size_t i = 0; i < iNum; ++i)
		vec.push_back(std::lerp(tmin, tmax, REAL(i)/REAL(iNum-1)));
	return vec;
}


template<typename T=double, typename REAL=double,
template<class...> class t_vec = std::vector>
t_vec<T> logspace(const T& tmin, const T& tmax, std::size_t iNum, T tBase = T(10))
{
	t_vec<T> vec = linspace<T, REAL>(tmin, tmax, iNum);

	for(T& t : vec)
		t = std::pow(tBase, t);
	return vec;
}


template<class T>
bool is_in_range(T val, T centre, T pm)
{
	pm = std::abs(pm);

	if(val < centre-pm)
		return false;
	if(val > centre+pm)
		return false;
	return true;
}


/**
 * point contained in linear range?
 */
template<class T = double>
bool is_in_linear_range(T dStart, T dStop, T dPoint)
{
	if(dStop < dStart)
		std::swap(dStart, dStop);

	return (dPoint >= dStart) && (dPoint <= dStop);
}


/**
 * angle contained in angular range?
 */
template<class T = double>
bool is_in_angular_range(T dStart, T dRange, T dAngle)
{
	if(dStart < T(0)) dStart += T(2)*pi<T>;
	if(dAngle < T(0)) dAngle += T(2)*pi<T>;

	dStart = std::fmod(dStart, T(2)*pi<T>);
	dAngle = std::fmod(dAngle, T(2)*pi<T>);

	T dStop = dStart + dRange;


	// if the end point is contained in the circular range
	if(dStop < T(2)*pi<T>)
	{
		return is_in_linear_range<T>(dStart, dStop, dAngle);
	}
	// else end point wraps around
	else
	{
		return is_in_linear_range<T>(dStart, T(2)*pi<T>, dAngle) ||
			is_in_linear_range<T>(T(0), dRange-(T(2)*pi<T>-dStart), dAngle);
	}
}


/**
 * get a random number in the given range
 */
template<class t_num>
t_num get_rand(t_num min = 1, t_num max = -1)
{
	static std::mt19937 rng{std::random_device{}()};

	if(max <= min)
	{
		min = std::numeric_limits<t_num>::lowest() / 10.;
		max = std::numeric_limits<t_num>::max() / 10.;
	}

	if constexpr(std::is_integral_v<t_num>)
		return std::uniform_int_distribution<t_num>(min, max)(rng);
	else
		return std::uniform_real_distribution<t_num>(min, max)(rng);
}

}

#endif
