/**
 * tlibs2 maths library -- statistical functions
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

#ifndef __TLIBS2_MATHS_STATISTICS_H__
#define __TLIBS2_MATHS_STATISTICS_H__

#include <cmath>
#include <vector>
#include <limits>

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// statistical functions
// ----------------------------------------------------------------------------

/**
 * mean value
 */
template<class t_elem, template<class...> class t_cont /*= std::vector*/>
t_elem mean(const t_cont<t_elem>& vec)
requires is_basic_vec<t_cont<t_elem>>
{
	if(vec.size()==0) return t_elem{};
	else if(vec.size()==1) return *vec.begin();

	using namespace tl2_ops;

	//t_elem meanvec = std::accumulate(std::next(vec.begin(), 1), vec.end(), *vec.begin());

	t_elem meanvec = *vec.begin();
	auto iter = std::next(vec.begin(), 1);
	for(; iter!=vec.end(); iter=std::next(iter, 1))
		meanvec += *iter;

	meanvec /= vec.size();
	return meanvec;
}


/**
 * mean value with given probability
 */
template<class t_vec_prob, class t_vec>
typename t_vec::value_type mean(const t_vec_prob& vecP, const t_vec& vec)
requires is_basic_vec<t_vec> && is_basic_vec<t_vec_prob>
{
	typedef typename t_vec::value_type T;
	typedef typename t_vec_prob::value_type Tprob;
	std::size_t iSize = std::min(vecP.size(), vec.size());

	if(iSize==0) return T(0);

	T tMean = vecP[0]*vec[0];
	Tprob tProbTotal = vecP[0];
	for(std::size_t i=1; i<iSize; ++i)
	{
		tMean += vecP[i]*vec[i];
		tProbTotal += vecP[i];
	}
	tMean /= tProbTotal;

	return tMean;
}


/**
 * standard deviation of mean value, with correction factor
 * @see https://en.wikipedia.org/wiki/Bessel%27s_correction
 */
template<class t_vec>
typename t_vec::value_type std_dev(const t_vec& vec, bool bCorr=1)
requires is_basic_vec<t_vec>
{
	typedef typename t_vec::value_type T;
	if(vec.size()<=1) return T(0);

	T tProb = T(vec.size());
	if(bCorr) tProb -= T(1);

	T tMean = mean(vec);
	T t = T(0);
	for(const T& tval : vec)
		t += (tval-tMean) * (tval-tMean);
	t /= tProb;

	return std::sqrt(t);
}


/**
 * standard deviation with given probability
 * @see https://en.wikipedia.org/wiki/Standard_deviation
 */
template<class t_vec_prob, class t_vec>
typename t_vec::value_type std_dev(const t_vec_prob& vecP, const t_vec& vec)
requires is_basic_vec<t_vec> && is_basic_vec<t_vec_prob>
{
	typedef typename t_vec::value_type T;
	std::size_t iSize = std::min(vecP.size(), vec.size());
	if(iSize<=1) return T(0);

	T tMean = mean<t_vec_prob, t_vec>(vecP, vec);
	T t = T(0);
	T tProbTotal = T(0);

	for(std::size_t iIdx = 0; iIdx<iSize; ++iIdx)
	{
		t += (vec[iIdx]-tMean)*(vec[iIdx]-tMean) * vecP[iIdx];
		tProbTotal += vecP[iIdx];
	}
	t /= tProbTotal;

	return std::sqrt(t);
}


/**
 * calculates minimum and maximum components of a collection of vectors
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, t_vec> minmax(const t_cont<t_vec>& verts)
requires is_vec<t_vec>
{
	using namespace tl2_ops;
	using t_real = typename t_vec::value_type;

	if(!verts.size())
		return std::make_tuple(t_vec{}, t_vec{});

	// set to limit values
	t_vec vecmin = zero<t_vec>(verts.begin()->size());
	t_vec vecmax = zero<t_vec>(verts.begin()->size());

	for(std::size_t i=0; i<vecmin.size(); ++i)
	{
		vecmin[i] = std::numeric_limits<t_real>::max();
		vecmax[i] = std::numeric_limits<t_real>::lowest();
	}

	// iterate components
	for(std::size_t i=0; i<vecmin.size(); ++i)
	{
		// iterate vectors
		for(const t_vec& vec : verts)
		{
			vecmin[i] = std::min(vecmin[i], vec[i]);
			vecmax[i] = std::max(vecmax[i], vec[i]);
		}
	}

	return std::make_tuple(vecmin, vecmax);
}


/**
 * calculates the covariance and the correlation matrices
 * covariance: C_ij = cov(X_i, X_j) = < (X_i - <X_i>) * (X_j - <X_j>) >
 * correlation: K_ij = C_ij / (sigma_i sigma_j)
 *
 * @see http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
 * @see (Arfken 2013) pp. 1142-1144
 * @see (Arens 2015), p. 795 and p. 1372
 */
template<class t_mat, class t_vec, class T=typename t_vec::value_type>
std::tuple<t_mat, t_mat>
covariance(const std::vector<t_vec>& vecVals, const std::vector<T>* pProb = nullptr)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using t_vecvec = typename std::remove_reference<decltype(vecVals)>::type;
	using t_innervec_org = decltype(vecVals[0]);
	using t_innervec = typename std::remove_const<
		typename std::remove_reference<t_innervec_org>::type>::type;

	if(vecVals.size() == 0)
		return std::make_tuple(t_mat(), t_mat());

	// mean vector <X_i>
	t_innervec vecMean;
	if(pProb)
		vecMean = mean<std::vector<T>, t_vecvec>(*pProb, vecVals);
	else
		vecMean = mean<t_vec>(vecVals);

	t_mat matCov = zero<t_mat>(vecVals[0].size(), vecVals[0].size());
	T tSum = T{0};
	const std::size_t N = vecVals.size();

	for(std::size_t i=0; i<N; ++i)
	{
		T tprob = T{1};

		// X_i - <X_i>
		t_innervec vec = vecVals[i] - vecMean;

		// matrix elements, AA^t
		t_mat matOuter = outer<t_mat, t_vec>(vec, vec);

		// probabilities for final averaging, <...>
		if(pProb)
		{
			tprob = (*pProb)[i];
			matOuter *= tprob;
		}

		matCov += matOuter;
		tSum += tprob;
	}

	// average, sometimes defined as C /= (N-1)
	matCov /= tSum /*-T(1)*/;


	// --------------------------------------------------------------------------------
	// correlation matrix
	t_innervec vecVar = diag_vec<t_vec, t_mat>(matCov);
	t_innervec vecStdDev(vecVar.size());

	std::transform(vecVar.begin(), vecVar.end(), vecStdDev.begin(),
		[](typename t_innervec::value_type d) -> typename t_innervec::value_type
		{ return std::sqrt(d); });

	t_mat matStdDev = outer<t_mat, t_vec>(vecStdDev, vecStdDev);
	t_mat matCorr = div_perelem<t_mat>(matCov, matStdDev);
	// --------------------------------------------------------------------------------

	return std::make_tuple(matCov, matCorr);
}


/**
 * calculates chi^2 distance of a function model to data points
 * chi^2 = sum( (y_i - f(x_i))^2 / sigma_i^2 )
 *
 * @see (Arfken 2013), p. 1170
 */
template<class T, class t_func, class t_iter_dat=T*>
T chi2(const t_func& func, std::size_t N,
	const t_iter_dat x, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T{0};

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - func(T(x[i]));
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}


/**
 * chi^2 for vector types
 *
 * @see (Merziger 2006), p. 185
 */
template<class t_vec, class t_func>
typename t_vec::value_type chi2(const t_func& func,
	const t_vec& x, const t_vec& y, const t_vec& dy)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	return chi2<T, t_func, T*>(func, x.size(), x.data(), y.data(),
		dy.size() ? dy.data() : nullptr);
}


/**
 * chi^2 which doesn't use an x value, but an index instead: y[idx] - func(idx)
 * @see (Arfken 2013), p. 1170
 */
template<class T, class t_func, class t_iter_dat=T*>
T chi2_idx(const t_func& func, std::size_t N, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - func(i);
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}


/**
 * direct chi^2 calculation with a model array instead of a model function
 * @see (Arfken 2013), p. 1170
 */
template<class T, class t_iter_dat=T*>
T chi2_direct(std::size_t N, const t_iter_dat func_y, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - T(func_y[i]);
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}



/**
 * multi-dimensional chi^2 function
 * @see (Arfken 2013), p. 1170
 */
template<class T, class T_dat, class t_func, template<class...> class t_vec=std::vector>
T chi2_nd(const t_func& func,
	const t_vec<t_vec<T_dat>>& vecvecX, const t_vec<T_dat>& vecY, const t_vec<T_dat>& vecDY)
{
	T tchi2 = T(0);

	for(std::size_t i=0; i<vecvecX.size(); ++i)
	{
		T td = T(vecY[i]) - func(vecvecX[i]);
		T tdy = vecDY[i];

		if(std::abs(tdy) < std::numeric_limits<T_dat>::min())
			tdy = std::numeric_limits<T_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}
// ----------------------------------------------------------------------------

}

#endif
