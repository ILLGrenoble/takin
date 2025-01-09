/**
 * tlibs2 maths library -- special functions
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

#ifndef __TLIBS2_MATHS_PEAKS_H__
#define __TLIBS2_MATHS_PEAKS_H__

#include <cmath>
#include <complex>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "decls.h"
#include "constants.h"


#ifdef __TLIBS2_USE_FADDEEVA__
	#include <Faddeeva.hh>
	using t_real_fadd = double;
#else
	//#pragma message("tlibs2: Disabling Faddeeva library (not found).")
#endif



namespace tl2 {
// -----------------------------------------------------------------------------
// peak model functions
// -----------------------------------------------------------------------------
/**
 * gaussian
 * @see https://en.wikipedia.org/wiki/Gaussian_function
 */
template<class T = double>
T gauss_model(T x, T x0, T sigma, T amp, T offs)
{
	T norm = T(1)/(std::sqrt(T(2)*pi<T>) * sigma);
	return amp * norm * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}


template<class T = double>
T gauss_model_amp(T x, T x0, T sigma, T amp, T offs)
{
	return amp * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}


template<class T = double>
T gauss_model_amp_slope(T x, T x0, T sigma, T amp, T offs, T slope)
{
	return amp * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + (x-x0)*slope + offs;
}


/**
 * lorentzian
 * @see https://en.wikipedia.org/wiki/Cauchy_distribution
 */
template<class T = double>
T lorentz_model_amp(T x, T x0, T hwhm, T amp, T offs)
{
	return amp*hwhm*hwhm / ((x-x0)*(x-x0) + hwhm*hwhm) + offs;
}


template<class T = double>
T lorentz_model_amp_slope(T x, T x0, T hwhm, T amp, T offs, T slope)
{
	return amp*hwhm*hwhm / ((x-x0)*(x-x0) + hwhm*hwhm) + (x-x0)*slope + offs;
}


template<class T = double>
T parabola_model(T x, T x0, T amp, T offs)
{
	return amp*(x-x0)*(x-x0) + offs;
}


template<class T = double>
T parabola_model_slope(T x, T x0, T amp, T offs, T slope)
{
	return amp*(x-x0)*(x-x0) + (x-x0)*slope + offs;
}


// -----------------------------------------------------------------------------
template<class t_real_to, class t_real_from,
bool bIsEqu = std::is_same<t_real_from, t_real_to>::value>
struct complex_cast
{
	const std::complex<t_real_to>& operator()(const std::complex<t_real_from>& c) const
	{ return c; }
};

template<class t_real_to, class t_real_from>
struct complex_cast<t_real_to, t_real_from, 0>
{
	std::complex<t_real_to> operator()(const std::complex<t_real_from>& c) const
	{ return std::complex<t_real_to>(t_real_to(c.real()), t_real_to(c.imag())); }
};
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// Faddeeva function
// -----------------------------------------------------------------------------
#ifdef __TLIBS2_USE_FADDEEVA__

/**
 * Complex error function
 */
template<class T=double>
std::complex<T> erf(const std::complex<T>& z)
{
	complex_cast<t_real_fadd, T> cst;
	complex_cast<T, t_real_fadd> inv_cst;
	return inv_cst(::Faddeeva::erf(cst(z)));
}

/**
 * Complex complementary error function
 */
template<class T=double>
std::complex<T> erfc(const std::complex<T>& z)
{
	complex_cast<t_real_fadd, T> cst;
	complex_cast<T, t_real_fadd> inv_cst;
	return inv_cst(::Faddeeva::erfc(cst(z)));
}

/**
 * Faddeeva function
 * @see https://en.wikipedia.org/wiki/Faddeeva_function
 */
template<class T=double>
std::complex<T> faddeeva(const std::complex<T>& z)
{
	std::complex<T> i(0, 1.);
	return std::exp(-z*z) * erfc(-i*z);
}

/**
 * Voigt profile
 * @see https://en.wikipedia.org/wiki/Voigt_profile
 */
template<class T=double>
T voigt_model(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	T norm = T(1)/(std::sqrt(T(2)*pi<T>) * sigma);
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));

	return amp*norm * faddeeva<T>(z).real() + offs;
}

template<class T=double>
T voigt_model_amp(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));
	return amp * faddeeva<T>(z).real() + offs;
}

template<class T=double>
T voigt_model_amp_slope(T x, T x0, T sigma, T gamma, T amp, T offs, T slope)
{
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));
	return amp * faddeeva<T>(z).real() + (x-x0)*slope + offs;
}

#endif
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// physics-related functions
// -----------------------------------------------------------------------------

/**
 * wrapper for boost's Y function
 */
template<class T=double>
std::complex<T> Ylm(int l /*0..i*/, int m /*-l..l*/, T th /*0..pi*/, T ph /*0..2pi*/)
{
	return boost::math::spherical_harmonic<T,T>(l,m, th, ph);
}



/**
 * CG coefficients
 * @see (Arfken 2013), p. 790 for the formula
 *
 * e.g. two e- spins: s1 = s2 = 0.5, ms[1,2] = 0.5 (up) or -0.5 (down), S = 0 (sing.) or 1 (trip.)
 */
template<class T = double>
T CG_coeff(T S, T s1, T s2, T ms1, T ms2)
{
	T (*fak)(T) = [](T t) -> T { return boost::math::factorial<T>(t); };

	T tCG = fak(S + s1 - s2)*fak(S - s1 + s2)*fak(-S + s1 + s2);
	tCG *= (T(2)*S + T(1));
	tCG *= fak(S + ms1 + ms2) * fak(S - (ms1 + ms2));
	tCG *= fak(s1 + ms1) * fak(s1 - ms1);
	tCG *= fak(s2 + ms2) * fak(s2 - ms2);
	tCG /= fak(S + s1 + s2 + T(1));
	tCG = std::sqrt(tCG);

	auto k_fkt = [&](T k) -> T
	{
		T t = std::pow(T(-1), k);
		t /= fak(k);
		t /= fak(-S + s1 + s2 - k)*fak(S - s1 - ms2 + k)*fak(S - s2 + ms1 + k);
		t /= fak(s2 + ms2 - k)*fak(s1 - ms1 - k);
		return t;
	};

	auto k_minmax = [&]() -> std::pair<T,T>
	{
		T kmax = s1 - ms1;
		kmax = std::min(kmax, s2 + ms2);
		kmax = std::min(kmax, -S + s1 + s2);

		T kmin = -(S - s1 - ms2);
		kmin = std::max(kmin, -(S - s2 + ms1));
		kmin = std::max(kmin, T(0));

		return std::make_pair(kmin, kmax);
	};

	T kmin, kmax;
	std::tie(kmin, kmax) = k_minmax();
	T kfact = T(0);
	for(T k=kmin; k<=kmax; k+=T(1))
		kfact += k_fkt(k);
	tCG *= kfact;

	return tCG;
}
// -----------------------------------------------------------------------------

}

#endif
