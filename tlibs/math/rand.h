/**
 * random numbers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
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

#ifndef __TLIBS_RAND_H__
#define __TLIBS_RAND_H__

#include <random>
#include <vector>
#include <array>
#include <initializer_list>
#include <future>
#include <type_traits>
#include <limits>

#include <boost/type_traits/function_traits.hpp>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>

#include "../helper/traits.h"


namespace tl {

extern std::mt19937& get_randeng();


// ----------------------------------------------------------------------------
// initialisers

extern unsigned int get_rand_seed();
extern void init_rand();
extern void init_rand_seed(unsigned int uiSeed, bool debug_msg = true);



// ----------------------------------------------------------------------------
// very simple distributions

extern unsigned int simple_rand(unsigned int iMax);



// ----------------------------------------------------------------------------
// simple distributions

/**
 * generates a random integer between iMin and iMax
 */
template<typename INT>
INT rand_int(INT iMin, INT iMax)
{
	if(iMin > iMax)
		std::swap(iMin, iMax);

	std::uniform_int_distribution<INT> dist(iMin, iMax);
	return dist(get_randeng());
}


/**
 * generates a random real between iMin and iMax
 */
template<typename REAL>
REAL rand_real(REAL dMin, REAL dMax)
{
	if(dMin > dMax)
		std::swap(dMin, dMax);

	std::uniform_real_distribution<REAL> dist(dMin, dMax);
	return dist(get_randeng());
}


/**
 * generates a random real between 0 and 1
 */
template<typename T>
T rand01()
{
	return rand_real<T>(T(0), T(1));
}


/**
 * returns 'true' with probability p
 */
template<typename T>
bool rand_prob(T p)
{
	T rnd = rand01<T>();
	return rnd <= p;
}


template<typename INT> struct _rand_int
{
	INT operator()(INT iMin, INT iMax) const
	{ return rand_int<INT>(iMin, iMax); }
};


template<typename REAL> struct _rand_real
{
	REAL operator()(REAL dMin, REAL dMax) const
	{ return rand_real<REAL>(dMin, dMax); }
};


/**
 * generates a random real or integer between tMin and tMax
 */
template<typename T>
static T rand_minmax(T tMin, T tMax)
{
	using t_rand = typename
		std::conditional<std::is_integral<T>::value,
		_rand_int<T>, _rand_real<T>>::type;

	t_rand rand;
	return rand(tMin, tMax);
}


/**
 * specialisation for 'bool'
 */
template<> bool rand_minmax(bool tMin, bool tMax);


// ----------------------------------------------------------------------------
// continuous distributions

/**
 * generates Gaussian-distributed random numbers
 */
template<typename REAL>
REAL rand_norm(REAL dMu, REAL dSigma)
{
	std::normal_distribution<REAL> dist(dMu, dSigma);
	return dist(get_randeng());
}


/**
 * generates Gaussian-distributed random numbers, including the normalising prefactor
 */
template<typename REAL>
REAL rand_norm_normalised(REAL dMu, REAL dSigma)
{
	std::normal_distribution<REAL> dist(dMu, dSigma);
	REAL prob = dist(get_randeng());
	prob /= std::sqrt(REAL(2) * REAL(M_PI)) * dSigma;
	return prob;
}


/**
 * generates lognormally-distributed random numbers
 */
template<typename REAL>
REAL rand_lognorm(REAL dMu, REAL dSigma)
{
	std::lognormal_distribution<REAL> dist(dMu, dSigma);
	return dist(get_randeng());
}


/**
 * generates Cauchy-distributed random numbers
 */
template<typename REAL>
REAL rand_cauchy(REAL dMu, REAL dSigma)
{
	std::cauchy_distribution<REAL> dist(dMu, dSigma);
	return dist(get_randeng());
}


/**
 * generates chi^2-distributed random numbers
 */
template<typename REAL=double>
REAL rand_chi2(REAL dDof)
{
	std::chi_squared_distribution<REAL> dist(dDof);
	return dist(get_randeng());
}


/**
 * generates t-distributed random numbers
 */
template<typename REAL=double>
REAL rand_student(REAL dDof)
{
	std::student_t_distribution<REAL> dist(dDof);
	return dist(get_randeng());
}


/**
 * generates f-distributed random numbers
 */
template<typename REAL=double>
REAL rand_fisher(REAL dDof1, REAL dDof2)
{
	std::fisher_f_distribution<REAL> dist(dDof1, dDof2);
	return dist(get_randeng());
}


/**
 * generates Gamma-distributed random numbers
 */
template<typename REAL=double>
REAL rand_gamma(REAL dMu, REAL dSigma)
{
	std::gamma_distribution<REAL> dist(dMu, dSigma);
	return dist(get_randeng());
}


/**
 * generates exp-distributed random numbers
 */
template<typename REAL=double>
REAL rand_exp(REAL dLam)
{
	std::exponential_distribution<REAL> dist(dLam);
	return dist(get_randeng());
}



// ----------------------------------------------------------------------------
// discrete distributions

/**
 * generates Poisson-distributed random numbers
 */
template<typename INT, typename REAL=double>
INT rand_poisson(REAL dLam)
{
	std::poisson_distribution<INT> dist(dLam);
	return dist(get_randeng());
}


/**
 * generates binomially-distributed random numbers
 */
template<typename INT, typename REAL=double>
INT rand_binomial(INT n, REAL p)
{
	std::binomial_distribution<INT> dist(n, p);
	return dist(get_randeng());
}



// ----------------------------------------------------------------------------
// multi-dimensional distributions

template<class _t_func,
	template<class...> class t_vec = std::vector>
struct _rand_nd
{
	static constexpr std::size_t iNumArgs = boost::function_traits<_t_func>::arity;

	using t_real = typename boost::function_traits<_t_func>::result_type;
	using t_func = std::function<_t_func>;

	t_func m_func;
	bool m_bUseThreads = 0;

	_rand_nd(const t_func& func) : m_func(func) {}

	t_vec<t_real> operator()(const t_vec<t_vec<t_real>>& vecParams) const
	{
		const std::size_t iDim = vecParams[0].size();
		t_vec<t_real> vecRet(iDim);

		std::vector<std::future<t_real>> vecFut;
		vecFut.reserve(iDim);

		using t_iter = typename t_vec<t_real>::const_iterator;
		t_vec<t_iter> vecIters;
		for(const t_vec<t_real>& vec : vecParams)
			vecIters.push_back(vec.begin());

		std::launch lpol = std::launch::deferred;
		if(m_bUseThreads) lpol |= std::launch::async;

		while(1)
		{
			t_vec<t_real> vecArgs;
			for(t_iter iter : vecIters)
				vecArgs.push_back(*iter);

			auto fkt = [this, vecArgs]() -> t_real
			{
				return tl::call<iNumArgs, t_func, t_real, t_vec>(m_func, vecArgs);
			};

			vecFut.emplace_back(std::async(lpol, fkt));

			for(t_iter& iter : vecIters)
				++iter;
			if(vecIters[0] == vecParams[0].end())
				break;
		}

		for(std::size_t i=0; i<iDim; ++i)
			vecRet[i] = vecFut[i].get();

		return vecRet;
	}
};


/**
 * generates n-dimensional uniformly distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_minmax_nd(const t_vec<t_real>& vecMin, const t_vec<t_real>& vecMax)
{
	_rand_nd<decltype(rand_minmax<t_real>), t_vec> rnd(rand_minmax<t_real>);
	return rnd(t_vec<t_vec<t_real>>({vecMin, vecMax}));
}


/**
 * generates n-dimensional Gaussian-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_norm_nd(const t_vec<t_real>& vecMu, const t_vec<t_real>& vecSigma)
{
	_rand_nd<decltype(rand_norm<t_real>), t_vec> rnd(rand_norm<t_real>);
	return rnd(t_vec<t_vec<t_real>>({vecMu, vecSigma}));
}


/**
 * generates n-dimensional lognorm-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_lognorm_nd(const t_vec<t_real>& vecMu, const t_vec<t_real>& vecSigma)
{
	_rand_nd<decltype(rand_lognorm<t_real>), t_vec> rnd(rand_lognorm<t_real>);
	return rnd({vecMu, vecSigma});
}


/**
 * generates n-dimensional Cauchy-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_cauchy_nd(const t_vec<t_real>& vecMu, const t_vec<t_real>& vecSigma)
{
	_rand_nd<decltype(rand_cauchy<t_real>), t_vec> rnd(rand_cauchy<t_real>);
	return rnd({vecMu, vecSigma});
}


/**
 * generates n-dimensional Gamma-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_gamma_nd(const t_vec<t_real>& vecMu, const t_vec<t_real>& vecSigma)
{
	_rand_nd<decltype(rand_gamma<t_real>), t_vec> rnd(rand_gamma<t_real>);
	return rnd({vecMu, vecSigma});
}


/**
 * generates n-dimensional f-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_fisher_nd(const t_vec<t_real>& vecDof1, const t_vec<t_real>& vecDof2)
{
	_rand_nd<decltype(rand_fisher<t_real>), t_vec> rnd(rand_fisher<t_real>);
	return rnd({vecDof1, vecDof2});
}


/**
 * generates n-dimensional t-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_student_nd(const t_vec<t_real>& vecDof)
{
	_rand_nd<decltype(rand_student<t_real>), t_vec> rnd(rand_student<t_real>);
	return rnd({vecDof});
}


/**
 * generates n-dimensional chi^2-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_chi2_nd(const t_vec<t_real>& vecDof)
{
	_rand_nd<decltype(rand_chi2<t_real>), t_vec> rnd(rand_chi2<t_real>);
	return rnd({vecDof});
}


/**
 * generates n-dimensional exp-distributed random numbers
 */
template<class t_real=double, template<class...> class t_vec=std::vector>
t_vec<t_real> rand_exp_nd(const t_vec<t_real>& vecDof)
{
	_rand_nd<decltype(rand_exp<t_real>), t_vec> rnd(rand_exp<t_real>);
	return rnd({vecDof});
}



// ----------------------------------------------------------------------------
// other distributions

/**
 * generates a random character sequence usable for temporary file names
 */
template<class t_str=std::string>
t_str rand_name(std::size_t iLen=8)
{
	static const t_str strChars = "abcdefghijklmnopqrstuvwxyz1234567890_";
	static const std::size_t iLenChars = strChars.length();

	t_str strRnd;
	strRnd.reserve(iLen);
	for(std::size_t iRnd=0; iRnd<iLen; ++iRnd)
		strRnd.push_back(strChars[tl::simple_rand(iLenChars)]);

	return strRnd;
}



// ----------------------------------------------------------------------------
// random array


/**
 * generates a random array
 */
template<class T, std::size_t DIM,
	template<class, std::size_t, class...> class t_arr_1d = boost::array,
	template<class, std::size_t, class...> class t_arr_nd = boost::multi_array>
t_arr_nd<T, DIM> rand_array(
	const t_arr_1d<typename t_arr_nd<T,DIM>::index, DIM>& arrDims)
{
	using t_arr = t_arr_nd<T, DIM>;

	T tMin = std::numeric_limits<T>::min();
	T tMax = std::numeric_limits<T>::max();

	t_arr arrRet(arrDims);
	std::size_t iLen = arrRet.num_elements();
	T* pRaw = arrRet.data();

	for(std::size_t i=0; i<iLen; ++i)
		pRaw[i] = rand_minmax<T>(tMin, tMax);

	return arrRet;
}


/**
 * generates random array indices in range [arrMin, arrMax[
 */
template<class T = std::size_t, std::size_t DIM,
	template<class, std::size_t, class...> class t_arr = std::array>
t_arr<T, DIM> rand_idx(
	const t_arr<T, DIM>& arrMin, const t_arr<T, DIM>& arrMax)
{
	t_arr<T, DIM> arrRet;
	for(std::size_t i=0; i<DIM; ++i)
		arrRet[i] = rand_int<T>(arrMin[i], arrMax[i]-1);

	return arrRet;
}

// ----------------------------------------------------------------------------

}
#endif
