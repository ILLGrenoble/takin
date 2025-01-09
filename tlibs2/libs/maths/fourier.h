/**
 * tlibs2 maths library -- Fourier transform
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

#ifndef __TLIBS2_MATHS_FOURIER_H__
#define __TLIBS2_MATHS_FOURIER_H__

#include <cmath>
#include <complex>
#include <vector>

#include "../algos.h"

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// Fourier transform
// @see (Scarpino 2011), ch. 14 for infos.
// ----------------------------------------------------------------------------

/**
 * dft
 * @see http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cplx dft_coeff(int k, const t_cont<t_cplx>& invec, bool bInv = false)
{
	const std::size_t N = invec.size();

	t_cplx imag(0., 1.);
	t_cplx f(0., 0.);

	for(std::size_t j=0; j<N; ++j)
	{
		T dv = T(-2)*pi < T>*T(j)*T(k)/T(N);
		if(bInv) dv = -dv;
		f += invec[j] * (std::cos(dv) + imag*std::sin(dv));
	}

	return f;
}


/**
 * dft
 * @see http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> dft(const t_cont<t_cplx>& invec,
	bool bInv = false, bool bNorm = false)
{
	const std::size_t N = invec.size();
	t_cont<t_cplx> outvec;
	outvec.resize(N);

	for(std::size_t k=0; k<N; ++k)
	{
		outvec[k] = dft_coeff<T, t_cplx, t_cont>(k, invec, bInv);

		if(bNorm && bInv)
			outvec[k] /= N;
	}

	return outvec;
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>>
requires is_complex<t_cplx>
t_cplx fft_factor(T N, T k, bool bInv = false)
{
	T ph = bInv ? -1 : 1.;

	T c = std::cos(T(2)*pi < T>*k/N * ph);
	T s = std::sin(T(2)*pi < T>*k/N * ph);

	return t_cplx{c, -s};
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> fft_reorder(const t_cont<t_cplx>& vecIn)
{
	t_cont<std::size_t> vecIdx =
		bit_reverse_indices<std::size_t, t_cont>(vecIn.size());

	t_cont<t_cplx> vecInRev;
	vecInRev.reserve(vecIn.size());

	for(std::size_t i = 0; i < vecIn.size(); ++i)
		vecInRev.push_back(vecIn[vecIdx[i]]);

	return vecInRev;
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> fft_merge(const t_cont<t_cplx>& vecIn, bool bInv = false)
{
	const std::size_t N = vecIn.size();
	const std::size_t N2 = N/2;

	if(N==0 || N==1)
		return vecIn;

	auto split_vec = [](const t_cont<t_cplx>& vec)
		-> std::pair<t_cont<t_cplx>, t_cont<t_cplx>>
	{
		std::size_t N = vec.size();

		t_cont<t_cplx> vec1, vec2;
		vec1.reserve(N/2);
		vec2.reserve(N/2);

		for(std::size_t i = 0; i < N/2; ++i)
		{
			vec1.push_back(vec[i]);
			vec2.push_back(vec[N/2 + i]);
		}

		return std::make_pair(std::move(vec1), std::move(vec2));
	};

	auto pair = split_vec(vecIn);
	t_cont<t_cplx> vec1 = fft_merge<T, t_cplx, t_cont>(pair.first, bInv);
	t_cont<t_cplx> vec2 = fft_merge<T, t_cplx, t_cont>(pair.second, bInv);

	t_cont<t_cplx> vecOut;
	vecOut.resize(N);

	for(std::size_t i = 0; i < N2; ++i)
	{
		vecOut[i] = vec1[i] + vec2[i]*fft_factor<T, t_cplx>(N, i, bInv);
		vecOut[N2+i] = vec1[i] + vec2[i]*fft_factor<T, t_cplx>(N, N2+i, bInv);
	}

	return vecOut;
}


/**
 * fft
 * @see (Scarpino 2011), ch. 14.
 */
template<typename T = double, class t_cplx = std::complex<T>,
	template<class...> class t_cont = std::vector>
requires is_complex<t_cplx>
t_cont<t_cplx> fft(const t_cont<t_cplx>& vecIn,
	bool bInv = false, bool bNorm = false)
{
	const std::size_t n = vecIn.size();
	t_cont<t_cplx> vecOut;
	vecOut.resize(n);

	vecOut = fft_reorder<T, t_cplx, t_cont>(vecIn);
	vecOut = fft_merge<T, t_cplx, t_cont>(vecOut, bInv);

	if(bInv && bNorm)
		for(t_cplx& c : vecOut)
			c /= n;

	return vecOut;
}
// ----------------------------------------------------------------------------


}

#endif
