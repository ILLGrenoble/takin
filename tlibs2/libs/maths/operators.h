/**
 * tlibs2 maths library -- operators
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Additional functions forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Further functions and updates forked on 1-Feb-2021 and 19-Apr-2021 from my privately developed "geo" and "misc" projects (https://github.com/t-weber/geo and https://github.com/t-weber/misc).
 * @note Additional functions forked on 6-Feb-2022 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS2_MATHS_OPS_H__
#define __TLIBS2_MATHS_OPS_H__

#include <cmath>
#include <tuple>
#include <vector>
#include <initializer_list>

#include <boost/algorithm/string.hpp>

#include "../str.h"

#include "decls.h"



namespace tl2_ops {
// ----------------------------------------------------------------------------
// vector operators
// ----------------------------------------------------------------------------

/**
 * unary +
 */
template<class t_vec>
const t_vec& operator+(const t_vec& vec1)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	return vec1;
}


/**
 * unary -
 */
template<class t_vec>
t_vec operator-(const t_vec& vec1)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = -vec1[i];

	return vec;
}


/**
 * binary +
 */
template<class t_vec>
t_vec operator+(const t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	if constexpr(tl2::is_dyn_vec<t_vec>)
		assert((vec1.size() == vec2.size()));

	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] + vec2[i];

	return vec;
}


/**
 * binary -
 */
template<class t_vec>
t_vec operator-(const t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	return vec1 + (-vec2);
}


/**
 * vector * scalar
 */
template<class t_vec>
t_vec operator*(const t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] * d;

	return vec;
}


/**
 * scalar * vector
 */
template<class t_vec>
t_vec operator*(typename t_vec::value_type d, const t_vec& vec)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
	//&& !tl2::is_basic_mat<typename t_vec::value_type>	// hack!
{
	return vec * d;
}


/**
 * vector / scalar
 */
template<class t_vec>
t_vec operator/(const t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	t_vec vec(vec1.size());

	for(std::size_t i=0; i<vec1.size(); ++i)
		vec[i] = vec1[i] / d;

	return vec;
}


/**
 * vector += vector
 */
template<class t_vec>
t_vec& operator+=(t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 + vec2;
	return vec1;
}


/**
 * vector -= vector
 */
template<class t_vec>
t_vec& operator-=(t_vec& vec1, const t_vec& vec2)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 - vec2;
	return vec1;
}


/**
 * vector *= scalar
 */
template<class t_vec>
t_vec& operator*=(t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 * d;
	return vec1;
}

/**
 * vector /= scalar
 */
template<class t_vec>
t_vec& operator/=(t_vec& vec1, typename t_vec::value_type d)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec1 = vec1 / d;
	return vec1;
}


/**
 * operator <<
 */
template<class t_vec>
std::ostream& operator<<(std::ostream& ostr, const t_vec& vec)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	const std::size_t N = vec.size();

	for(std::size_t i=0; i<N; ++i)
	{
		ostr << vec[i];
		if(i < N-1)
			ostr << TL2_COLSEP << " ";
	}

	return ostr;
}


/**
 * operator >>
 */
template<class t_vec>
std::istream& operator>>(std::istream& istr, t_vec& vec)
requires tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	vec.clear();

	std::string str;
	std::getline(istr, str);

	std::vector<std::string> vecstr;
	boost::split(vecstr, str,
		[](auto c) -> bool { return c == TL2_COLSEP; },
		boost::token_compress_on);

	for(auto& tok : vecstr)
	{
		boost::trim(tok);
		typename t_vec::value_type c =
			tl2::stoval<typename t_vec::value_type>(tok);
		vec.emplace_back(std::move(c));
	}

	return istr;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// matrix operators
// ----------------------------------------------------------------------------

/**
 * unary +
 */
template<class t_mat>
const t_mat& operator+(const t_mat& mat1)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return mat1;
}


/**
 * unary -
 */
template<class t_mat>
t_mat operator-(const t_mat& mat1)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	t_mat mat(mat1.size1(), mat1.size2());

	for(std::size_t i = 0; i < mat1.size1(); ++i)
		for(std::size_t j = 0; j < mat1.size2(); ++j)
			mat(i,j) = -mat1(i,j);

	return mat;
}


/**
 * binary +
 */
template<class t_mat>
t_mat operator+(const t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat1.size1() == mat2.size1() && mat1.size2() == mat2.size2()));

	t_mat mat(mat1.size1(), mat1.size2());

	for(std::size_t i = 0; i < mat1.size1(); ++i)
		for(std::size_t j = 0; j < mat1.size2(); ++j)
			mat(i,j) = mat1(i,j) + mat2(i,j);

	return mat;
}


/**
 * binary -
 */
template<class t_mat>
t_mat operator-(const t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return mat1 + (-mat2);
}


/**
 * matrix * scalar
 */
template<class t_mat>
t_mat operator*(const t_mat& mat1, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	t_mat mat(mat1.size1(), mat1.size2());

	for(std::size_t i = 0; i < mat1.size1(); ++i)
		for(std::size_t j = 0; j < mat1.size2(); ++j)
			mat(i,j) = mat1(i,j) * d;

	return mat;
}

/**
 * scalar * matrix
 */
template<class t_mat>
t_mat operator*(typename t_mat::value_type d, const t_mat& mat)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return mat * d;
}


/**
 * matrix / scalar
 */
template<class t_mat>
t_mat operator/(const t_mat& mat, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	using T = typename t_mat::value_type;
	return mat * (T(1)/d);
}


/**
 * matrix-matrix product
 */
template<class t_mat>
t_mat operator*(const t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	return tl2::prod<t_mat>(mat1, mat2);
}


/**
 * matrix *= scalar
 */
template<class t_mat>
t_mat& operator*=(t_mat& mat1, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 * d;
	return mat1;
}


/**
* matrix *= matrix
 */
template<class t_mat>
t_mat& operator*=(t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 * mat2;
	return mat1;
}


/**
 * matrix += matrix
 */
template<class t_mat>
t_mat& operator+=(t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 + mat2;
	return mat1;
}


/**
 * matrix -= matrix
 */
template<class t_mat>
t_mat& operator-=(t_mat& mat1, const t_mat& mat2)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 - mat2;
	return mat1;
}


/**
 * matrix /= scalar
 */
template<class t_mat>
t_mat& operator/=(t_mat& mat1, typename t_mat::value_type d)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	mat1 = mat1 / d;
	return mat1;
}


/**
 * operator <<
 */
template<class t_mat>
std::ostream& operator<<(std::ostream& ostr, const t_mat& mat)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	const std::size_t ROWS = mat.size1();
	const std::size_t COLS = mat.size2();

	for(std::size_t row = 0; row < ROWS; ++row)
	{
		for(std::size_t col = 0; col < COLS; ++col)
		{
			ostr << mat(row, col);
			if(col < COLS-1)
				ostr << TL2_COLSEP << " ";
		}

		if(row < ROWS-1)
			ostr << TL2_ROWSEP << " ";
	}

	return ostr;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// mixed operators
// ----------------------------------------------------------------------------

/**
 * matrix-vector product: c_i = a_ij b_j
 */
template<class t_mat, class t_vec>
t_vec operator*(const t_mat& mat, const t_vec& vec)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
	&& tl2::is_basic_vec<t_vec> && tl2::is_dyn_vec<t_vec>
{
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size2() == vec.size()));
	else
		static_assert(t_mat::size2() == t_vec::size());


	t_vec vecRet(mat.size1());

	for(std::size_t row = 0; row < mat.size1(); ++row)
	{
		vecRet[row] = typename t_vec::value_type{/*0*/};
		for(std::size_t col = 0; col < mat.size2(); ++col)
		{
			auto elem = mat(row, col) * vec[col];
			vecRet[row] = vecRet[row] + elem;
		}
	}

	return vecRet;
}

// ----------------------------------------------------------------------------
}

#endif
