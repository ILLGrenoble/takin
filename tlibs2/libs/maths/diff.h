/**
 * tlibs2 maths library -- numerical differentiation and integration
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

#ifndef __TLIBS2_MATHS_DIFF_H__
#define __TLIBS2_MATHS_DIFF_H__

#include <cmath>
#include <vector>
#include <limits>

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// numerical differentiation and integration
// ----------------------------------------------------------------------------
/**
 * numerical differentiation
 */
template<class t_cont_out = std::vector<double>, class t_cont_in = std::vector<double>>
t_cont_out diff(const t_cont_in& xs, const t_cont_in& ys)
{
	using t_real = typename t_cont_in::value_type;

	const std::size_t N = xs.size();
	t_cont_out new_ys{};
	new_ys.reserve(N);

	for(std::size_t i = 0; i < N-1; ++i)
	{
		t_real diff_val = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i]);
		new_ys.push_back(diff_val);
	}

	// repeat last value
	new_ys[N-1] = new_ys[N-2];
	return new_ys;
}


/**
 * trapezoid rule
 * @see https://en.wikipedia.org/wiki/Trapezoidal_rule
 */
template<class R=double, class A=double>
R numint_trap(const std::function<R(A)>& fkt,
	A x0, A x1)
{
	return R(0.5)*R(x1-x0) * (fkt(x0) + fkt(x1));
}


/**
 * trapezoid rule
 * @see https://en.wikipedia.org/wiki/Trapezoidal_rule
 */
template<class R=double, class A=double>
R numint_trapN(const std::function<R(A)>& fkt,
	A x0, A x1, std::size_t N)
{
	const A xstep = A(x1-x0)/A(N);

	R xsum = fkt(x0) + fkt(x1);
	for(std::size_t i = 1; i < N; ++i)
		xsum += R(2)*fkt(x0 + A(i)*xstep);

	xsum *= R(0.5)*R(xstep);
	return R(xsum);
}


/**
 * rectangle rule
 * @see https://en.wikipedia.org/wiki/Rectangle_method
 */
template<class R=double, class A=double>
R numint_rect(const std::function<R(A)>& fkt,
	A x0, A x1, std::size_t N)
{
	const A xstep = (x1-x0)/A(N);

	R xsum = R(0);
	for(std::size_t i = 0; i < N; ++i)
		xsum += fkt(x0 + A(i)*xstep);

	xsum *= R(xstep);
	return xsum;
}


/**
 * Simpson's rule
 * @see https://en.wikipedia.org/wiki/Simpson%27s_rule
 */
template<class R=double, class A=double>
R numint_simp(const std::function<R(A)>& fkt, A x0, A x1)
{
	return (fkt(x0) + 4.*fkt(0.5*(x0+x1)) + fkt(x1)) * (x1-x0)/6.;
}


/**
 * Simpson's rule
 * @see https://en.wikipedia.org/wiki/Simpson%27s_rule
 */
template<class R=double, class A=double>
R numint_simpN(const std::function<R(A)>& fkt, A x0, A x1, std::size_t N)
{
	const A xstep = (x1-x0)/A(N);
	R xsum = fkt(x0) + fkt(x1);

	for(std::size_t i = 1; i <= N/2; ++i)
	{
		xsum += R(2) * fkt(x0 + A(2*i)*xstep);
		xsum += R(4) * fkt(x0 + A(2*i-1)*xstep);
	}
	xsum -= R(2)*fkt(x0 + A(2*N/2)*xstep);

	xsum *= R(xstep)/R(3);
	return xsum;
}


/**
 * convolution integral of fkt0 and fkt1
 * @see https://en.wikipedia.org/wiki/Convolution
 */
template<class R=double, class A=double>
R convolute(const std::function<R(A)>& fkt0,
	const std::function<R(A)>& fkt1,
	A x, A x0, A x1, std::size_t N)
{
	std::function<R(A,A)> fkt = [&fkt0, &fkt1](A t, A tau) -> R
	{
		return fkt0(tau) * fkt1(t-tau);
	};

	// ... at fixed arg x
	std::function<R(A)> fktbnd = [&fkt, &x](A t) -> R { return fkt(x, t); };

	return numint_simpN(fktbnd, x0, x1, N);
}


/**
 * discrete convolution
 * @see https://en.wikipedia.org/wiki/Convolution#Discrete_convolution
 */
template<class cont_type = std::vector<double>>
cont_type convolute_discrete(const cont_type& f, const cont_type& g)
{
	const std::size_t M = f.size();
	const std::size_t N = g.size();

	cont_type conv;
	conv.reserve(M+N-1);

	for(std::size_t n = 0; n < M+N-1; ++n)
	{
		typename cont_type::value_type val = 0.;

		for(std::size_t m=0; m<M; ++m)
			if(n>=m && n-m<N)
				val += f[m] * g[n-m];

		conv.push_back(val);
	}

	return conv;
}


/**
 * Newton iteration
 * @see https://en.wikipedia.org/wiki/Newton%27s_method
 */
template<class T = double>
T newton(const std::function<T(T)>& fkt, const std::function<T(T)>& diff,
	T x, std::size_t imax = 128, T eps = std::numeric_limits<T>::epsilon())
{
	T xnew = x;
	for(std::size_t i = 0; i < imax; ++i)
	{
		xnew = x - fkt(x)/diff(x);

		if(std::abs(xnew - x) < eps)
			break;

		x = xnew;
	}

	return xnew;
}
// ----------------------------------------------------------------------------


}

#endif
