/**
 * data point interpolation
 *
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 25-04-2013
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

#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <cmath>
#include <boost/math/special_functions/binomial.hpp>
#include <vector>
#include <algorithm>
#include <limits>

#include "../math/geo.h"
#include "../math/math.h"
#include "../helper/misc.h"
#include "../log/log.h"
#include "funcmod.h"

namespace tl {

namespace ublas = boost::numeric::ublas;


/**
 * see: http://mathworld.wolfram.com/BernsteinPolynomial.html
 */
template<typename T> T bernstein(int i, int n, T t)
{
	T bino = boost::math::binomial_coefficient<T>(n, i);
	return bino * pow(t, i) * pow(1-t, n-i);
}


/**
 * see: http://mathworld.wolfram.com/BezierCurve.html
 */
template<typename T>
ublas::vector<T> bezier(const ublas::vector<T>* P, std::size_t N, T t)
{
	if(N==0)
		return ublas::vector<T>(0);
	const int n = N-1;

	ublas::vector<T> vec(P[0].size());
	for(std::size_t i=0; i<vec.size(); ++i) vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bernstein(i, n, t);

	return vec;
}


/**
 * see: http://mathworld.wolfram.com/B-Spline.html
 */
template<typename T>
T bspline_base(int i, int j, T t, const std::vector<T>& knots)
{
	if(j==0)
	{
		if((knots[i] <= t) && (t < knots[i+1]) && (knots[i]<knots[i+1]))
			return 1.;
		return 0.;
	}

	T val11 = (t - knots[i]) / (knots[i+j]-knots[i]);
	T val12 = bspline_base(i, j-1, t, knots);
	T val1 = val11 * val12;

	T val21 = (knots[i+j+1]-t) / (knots[i+j+1]-knots[i+1]);
	T val22 = bspline_base(i+1, j-1, t, knots);
	T val2 = val21 * val22;

	T val = val1 + val2;
	return val;
}


/**
 * see: http://mathworld.wolfram.com/B-Spline.html
 */
template<typename T>
ublas::vector<T> bspline(const ublas::vector<T>* P, std::size_t N, T t, const std::vector<T>& knots)
{
	if(N==0)
		return ublas::vector<T>(0);
	const int n = N-1;
	const int m = knots.size()-1;
	const int degree = m-n-1;

	ublas::vector<T> vec(P[0].size());
	for(std::size_t i=0; i<vec.size(); ++i)
		vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bspline_base(i, degree, t, knots);

	return vec;
}


// ----------------------------------------------------------------------------


template<typename T=double>
class Bezier : public FunctionModel_param<ublas::vector<T>>
{
protected:
	std::unique_ptr<ublas::vector<T>[]> m_pvecs;
	std::size_t m_iN;

public:
	Bezier(std::size_t N, const T *px, const T *py);
	virtual ~Bezier() = default;

	virtual ublas::vector<T> operator()(T t) const override;
	virtual const char* GetModelName() const override { return "bezier"; };
};


template<typename T=double>
class BSpline : public FunctionModel_param<ublas::vector<T>>
{
protected:
	std::unique_ptr<ublas::vector<T>[]> m_pvecs;
	std::size_t m_iN, m_iDegree;
	std::vector<T> m_vecKnots;

public:
	BSpline(std::size_t N, const T *px, const T *py, unsigned int iDegree=3);
	virtual ~BSpline() = default;

	virtual ublas::vector<T> operator()(T t) const override;
	virtual const char* GetModelName() const override { return "bspline"; };
};


template<typename T=double>
class LinInterp : public FunctionModel<T>
{
protected:
	std::unique_ptr<ublas::vector<T>[]> m_pvecs;
	std::size_t m_iN;

public:
	LinInterp(std::size_t N, const T *px, const T *py);
	virtual ~LinInterp() = default;

	virtual T operator()(T x) const override;
	virtual const char* GetModelName() const override { return "linear_interpolation"; };
};



template<class T>
Bezier<T>::Bezier(std::size_t N, const T *px, const T *py) : m_iN(N)
{
	m_pvecs.reset(new ublas::vector<T>[m_iN]);

	for(std::size_t i=0; i<m_iN; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	//auto MinMax = boost::minmax_element(px, px+N);
	//m_dMin = *MinMax.first;
	//m_dMax = *MinMax.second;
}

template<class T>
ublas::vector<T> Bezier<T>::operator()(T t) const
{
	return tl::bezier<T>(m_pvecs.get(), m_iN, t);
}



template<class T>
BSpline<T>::BSpline(std::size_t N, const T *px, const T *py, unsigned int iDegree)
	: m_iN(N), m_iDegree(iDegree)
{
	m_pvecs.reset(new ublas::vector<T>[m_iN]);

	for(std::size_t i=0; i<m_iN; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	std::size_t iM = m_iDegree + m_iN + 1;
	m_vecKnots.resize(iM);

	const T eps = std::numeric_limits<T>::epsilon();

	// set knots to uniform, nonperiodic B-Spline
	for(unsigned int i=0; i<m_iDegree+1; ++i)
		m_vecKnots[i] = 0.+i*eps;
	for(unsigned int i=iM-m_iDegree-1; i<iM; ++i)
		m_vecKnots[i] = 1.-i*eps;
	for(unsigned int i=m_iDegree+1; i<iM-m_iDegree-1; ++i)
		m_vecKnots[i] = T(i+1-m_iDegree-1) / T(iM-2*m_iDegree-2 + 1);
}

template<class T>
ublas::vector<T> BSpline<T>::operator()(T t) const
{
	if(m_iN==0)
	{
		ublas::vector<T> vecNull(2);
		vecNull[0] = vecNull[1] = 0.;
		return vecNull;
	}

	ublas::vector<T> vec = tl::bspline<T>(m_pvecs.get(), m_iN, t, m_vecKnots);

	// remove epsilon dependence
	if(t<=0.) vec = m_pvecs[0];
	if(t>=1.) vec = m_pvecs[m_iN-1];

	return vec;
}


template<class T>
LinInterp<T>::LinInterp(std::size_t N, const T *px, const T *py)
	: m_iN(N)
{
	m_pvecs.reset(new ublas::vector<T>[m_iN]);

	for(std::size_t i=0; i<m_iN; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	// ensure that vector is sorted by x values
	std::stable_sort(m_pvecs.get(), m_pvecs.get()+m_iN,
		[](const ublas::vector<T>& vec1, const ublas::vector<T>& vec2) -> bool
		{ return vec1[0] < vec2[0]; });
}

template<class T>
T LinInterp<T>::operator()(T x) const
{
	const auto* iterBegin = m_pvecs.get();
	const auto* iterEnd = m_pvecs.get() + m_iN;

	if(m_iN == 0)
		return T(0);
	if(m_iN == 1)
		return (*iterBegin)[1];

	const auto* iterLower = std::lower_bound(iterBegin, iterEnd, x,
		[](const ublas::vector<T>& vec, const T& x) -> bool
		{ return vec[0] < x; });

	// lower bound at end of range?
	if(iterLower == iterEnd || iterLower == iterEnd-1)
		iterLower = iterEnd - 2;
	const auto* iter2 = iterLower + 1;

	T xrange = (*iter2)[0] - (*iterLower)[0];
	T xpos = (x-(*iterLower)[0]) / xrange;

	return tl::lerp<T,T>((*iterLower)[1], (*iter2)[1], xpos);
}

// ----------------------------------------------------------------------------

template<typename T>
void find_peaks(std::size_t iLen, const T* _px, const T* _py, unsigned int iOrder,
	std::vector<T>& vecMaximaX, std::vector<T>& vecMaximaSize, std::vector<T>& vecMaximaWidth,
	T eps = tl::get_epsilon<T>(), std::size_t iNumSpline = 512)
{
	// allocate memory
	std::unique_ptr<T, std::default_delete<T[]>> uptrMem{
		new T[2*iLen + 4*iNumSpline]};
	T *px = uptrMem.get() + 0*iLen;
	T *py = uptrMem.get() + 1*iLen;
	T *pSplineX = uptrMem.get() + 2*iLen + 0*iNumSpline;
	T *pSplineY = uptrMem.get() + 2*iLen + 1*iNumSpline;
	T *pSplineDiff = uptrMem.get() + 2*iLen + 2*iNumSpline;
	T *pSplineDiff2 = uptrMem.get() + 2*iLen + 3*iNumSpline;


	// sort input values
	std::copy(_px, _px+iLen, px);
	std::copy(_py, _py+iLen, py);
	tl::sort_2<T*>(px, px+iLen, py);


	BSpline<T> spline(iLen, px, py, iOrder);
	const T* pdyMin = std::min_element(py, py+iLen);

	for(std::size_t iSpline=0; iSpline<iNumSpline; ++iSpline)
	{
		const T dT = T(iSpline) / T(iNumSpline-1);
		ublas::vector<T> vec = spline(dT);

		pSplineX[iSpline] = vec[0];
		pSplineY[iSpline] = vec[1];
	}

	tl::diff(iNumSpline, pSplineX, pSplineY, pSplineDiff);
	tl::diff(iNumSpline, pSplineX, pSplineDiff, pSplineDiff2);
	std::vector<std::size_t> vecZeroes = tl::find_zeroes<T>(iNumSpline, pSplineDiff, eps);

	//std::ofstream ofDbg{"0.dat"};
	//for(std::size_t i=0; i<iNumSpline; ++i)
	//	ofDbg << pSplineX[i] << "\t" << pSplineY[i] << "\t" << pSplineDiff[i] << "\t" << pSplineDiff2[i] << std::endl;


	for(std::size_t iZeroIdx=0; iZeroIdx<vecZeroes.size(); ++iZeroIdx)
	{
		const std::size_t iZero = vecZeroes[iZeroIdx];

		// minima / saddle points
		if(pSplineDiff2[iZero] >= 0.)
			continue;

		vecMaximaX.push_back(pSplineX[iZero]);

		int iMinIdxLeft = -1;
		int iMinIdxRight = -1;
		if(iZeroIdx > 0)
			iMinIdxLeft = vecZeroes[iZeroIdx-1];
		if(iZeroIdx+1 < vecZeroes.size())
			iMinIdxRight = vecZeroes[iZeroIdx+1];

		T dHeight = 0.;
		T dWidth = 0.;
		T dDiv = 0.;

		// minimum left of the peak
		if(iMinIdxLeft>=0)
		{
			dHeight += (pSplineY[iZero]-pSplineY[iMinIdxLeft]);
			dWidth += std::abs((pSplineX[iZero]-pSplineX[iMinIdxLeft]));
			dDiv += 1.;
		}

		// minimum right of the peak
		if(iMinIdxRight>=0)
		{
			dHeight += (pSplineY[iZero]-pSplineY[iMinIdxRight]);
			dWidth += std::abs((pSplineX[iZero]-pSplineX[iMinIdxRight]));
			dDiv += 1.;
		}

		// no adjacent minima...
		if(iMinIdxLeft<0 && iMinIdxRight<0)
		{
			dHeight = pSplineY[iZero]- *pdyMin;
			dWidth = (px[iLen-1] - px[0]) / 10.;	// guess something...
			dDiv = 1.;
		}

		if(dDiv != 0.)
		{
			dHeight /= dDiv;
			dWidth /= dDiv;
		}

		vecMaximaSize.push_back(dHeight);
		vecMaximaWidth.push_back(dWidth);
	}

	tl::sort_3<typename std::vector<T>::iterator>(vecMaximaSize.begin(), vecMaximaSize.end(),
		vecMaximaWidth.begin(), vecMaximaX.begin());
	std::reverse(vecMaximaSize.begin(), vecMaximaSize.end());
	std::reverse(vecMaximaWidth.begin(), vecMaximaWidth.end());
	std::reverse(vecMaximaX.begin(), vecMaximaX.end());


	std::ostringstream ostrDbg;
	ostrDbg << "Prefitter found peaks at: ";
	for(std::size_t i=0; i<vecMaximaX.size(); ++i)
	{
		ostrDbg << vecMaximaX[i]
			<< " (int: " << vecMaximaSize[i]
			<< ", width: " << vecMaximaWidth[i] << "), ";
	}
	log_debug(ostrDbg.str());
}

}
#endif
