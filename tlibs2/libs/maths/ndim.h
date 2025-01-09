/**
 * tlibs2 maths library -- n-dimensional algorithms
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

#ifndef __TLIBS2_MATHS_NDIM_H__
#define __TLIBS2_MATHS_NDIM_H__

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>
#include <initializer_list>
#include <limits>

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// n-dim algorithms
// ----------------------------------------------------------------------------

/**
 * are two vectors equal within an epsilon range?
 */
template<class t_vec, class t_num = typename t_vec::value_type>
bool equals(const t_vec& vec1, const t_vec& vec2,
	t_num eps = std::numeric_limits<t_num>::epsilon(),
	int _maxSize = -1)
requires is_basic_vec<t_vec>
{
	// size has to be equal
	if(vec1.size() != vec2.size())
		return false;

	std::size_t maxSize = vec1.size();
	if(_maxSize >= 0)
		maxSize = std::min(std::size_t(_maxSize), maxSize);

	// check each element
	for(std::size_t i = 0; i < maxSize; ++i)
	{
		if constexpr(tl2::is_complex<t_num>)
		{
			if(!tl2::equals<t_num>(vec1[i], vec2[i], eps.real()))
				return false;
		}
		else if constexpr(tl2::is_scalar<t_num>)
		{
			if(!tl2::equals<t_num>(vec1[i], vec2[i], eps))
				return false;
		}
	}

	return true;
}


/**
 * are two matrices equal within an epsilon range?
 */
template<class t_mat, class t_real>
bool equals(const t_mat& mat1, const t_mat& mat2,
	t_real eps = std::numeric_limits<t_real>::epsilon(),
	int _maxSize = -1)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;

	if(mat1.size1() != mat2.size1() || mat1.size2() != mat2.size2())
		return false;

	std::size_t maxSize1 = mat1.size1();
	std::size_t maxSize2 = mat1.size2();
	if(_maxSize >= 0)
	{
		maxSize1 = std::min(std::size_t(_maxSize), maxSize1);
		maxSize2 = std::min(std::size_t(_maxSize), maxSize2);
	}

	for(std::size_t i = 0; i < maxSize1; ++i)
	{
		for(std::size_t j = 0; j < maxSize2; ++j)
		{
			if constexpr(is_complex<decltype(eps)>)
			{
				if(!tl2::equals<T>(mat1(i, j), mat2(i, j), eps.real()))
					return false;
			}
			else
			{
				if(!tl2::equals<T>(mat1(i, j), mat2(i, j), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * check if two collections of matrices or vectors are equal
 */
template<class t_obj, template<class...> class t_vec = std::vector>
bool equals_all(const t_vec<t_obj>& vec1, const t_vec<t_obj>& _vec2,
	typename t_obj::value_type eps = std::numeric_limits<typename t_obj::value_type>::epsilon(),
	int maxSize=-1)
{
	auto vec2 = _vec2;
	if(vec1.size() != vec2.size())
		return false;

	for(const auto& obj1 : vec1)
	{
		// find obj1 in vec2
		auto iter = std::find_if(vec2.crbegin(), vec2.crend(),
		[&obj1, eps, maxSize](const t_obj& obj2) -> bool
		{
			return tl2::equals<t_obj>(obj1, obj2, eps, maxSize);
		});

		// not found
		if(iter == vec2.crend())
			return false;

		// remove already checked element
		vec2.erase(iter.base()-1);
	}

	return true;
}


/**
 * remove duplicate vectors or matrices in the container
 */
template<class t_obj, template<class...> class t_cont = std::vector>
t_cont<t_obj> remove_duplicates(const t_cont<t_obj>& objs,
	typename t_obj::value_type eps = std::numeric_limits<typename t_obj::value_type>::epsilon())
{
	t_cont<t_obj> newobjs;
	newobjs.reserve(objs.size());

	for(const auto& elem : objs)
	{
		// find obj in container
		auto iter = std::find_if(newobjs.cbegin(), newobjs.cend(),
		[&elem, eps](const t_obj& elem2) -> bool
		{
			return tl2::equals<t_obj>(elem, elem2, eps);
		});

		// not found
		if(iter == newobjs.cend())
			newobjs.push_back(elem);
	}

	return newobjs;
}


/**
 * set submatrix to unit
 */
template<class t_mat>
void unit(t_mat& mat, std::size_t rows_begin, std::size_t cols_begin,
	std::size_t rows_end, std::size_t cols_end)
requires is_basic_mat<t_mat>
{
	for(std::size_t i=rows_begin; i < rows_end; ++i)
		for(std::size_t j=cols_begin; j < cols_end; ++j)
			mat(i, j) = (i==j ? 1 : 0);
}


/**
 * unit matrix
 */
template<class t_mat>
t_mat unit(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	unit<t_mat>(mat, 0,0, mat.size1(),mat.size2());
	return mat;
}


/**
 * unit matrix
 */
template<class t_mat>
t_mat unit(std::size_t N)
requires is_basic_mat<t_mat>
{
	return unit<t_mat>(N, N);
}


/**
 * is mat a unit matrix?
 */
template<class t_mat, class t_scalar = typename t_mat::value_type>
bool is_unit(const t_mat& mat,
	t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_mat<t_mat>
{
	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		for(std::size_t j = 0; j < mat.size2(); ++j)
		{
			if(i==j && !tl2::equals<t_scalar>(mat(i, j), t_scalar(1), eps))
				return false;
			if(i!=j && !tl2::equals<t_scalar>(mat(i, j), t_scalar(0), eps))
				return false;
		}
	}

	return true;
}


/**
 * zero matrix
 */
template<class t_mat>
t_mat zero(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>
{
	using size_t = decltype(t_mat{}.size1());

	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(size_t i = 0; i < mat.size1(); ++i)
		for(size_t j = 0; j < mat.size2(); ++j)
			mat(i, j) = 0;

	return mat;
}


/**
 * zero matrix
 */
template<class t_mat>
t_mat zero(std::size_t N=0)
requires is_basic_mat<t_mat>
{
	return zero<t_mat>(N, N);
}


/**
 * is mat a zero matrix?
 */
template<class t_mat, class t_scalar = typename t_mat::value_type>
bool is_zero(const t_mat& mat,
	t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_mat<t_mat>
{
	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		for(std::size_t j = 0; j < mat.size2(); ++j)
		{
			if(!tl2::equals<t_scalar>(mat(i, j), t_scalar(0), eps))
				return false;
		}
	}

	return true;
}


/**
 * zero vector
 */
template<class t_vec>
t_vec zero(std::size_t N /* = 0*/)
requires is_basic_vec<t_vec>
{
	using size_t = decltype(t_vec{}.size());

	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(N);

	for(size_t i = 0; i < vec.size(); ++i)
		vec[i] = 0;

	return vec;
}


/**
 * is vec a zero vector?
 */
template<class t_vec, class t_scalar = typename t_vec::value_type>
bool is_zero(const t_vec& vec,
	t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
requires is_vec<t_vec>
{
	for(std::size_t i = 0; i < vec.size(); ++i)
	{
		if(!tl2::equals<t_scalar>(vec[i], t_scalar(0), eps))
			return false;
	}

	return true;
}


/**
 * random scalar
 */
template<typename t_scalar>
t_scalar rand()
requires is_scalar<t_scalar> && (!is_complex<t_scalar>)
{
	return get_rand<t_scalar>(0, 1);
};


/**
 * random complex number
 */
template<typename t_cplx, class t_scalar = typename t_cplx::value_type>
t_cplx rand()
requires is_complex<t_cplx> && is_scalar<t_scalar>
{
	return t_cplx(rand<t_scalar>(), rand<t_scalar>());
};


/**
 * random vector
 */
template<class t_vec>
t_vec rand(std::size_t N)
requires is_basic_vec<t_vec>
{
	using size_t = decltype(t_vec{}.size());
	using t_scalar = typename t_vec::value_type;

	t_vec vec{};
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(N);

	for(size_t i = 0; i < vec.size(); ++i)
		vec[i] = rand<t_scalar>();

	return vec;
}


/**
 * random matrix
 */
template<class t_mat>
t_mat rand(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>
{
	using size_t = decltype(t_mat{}.size1());
	using t_scalar = typename t_mat::value_type;

	t_mat mat{};
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(size_t i = 0; i < mat.size1(); ++i)
		for(size_t j = 0; j < mat.size2(); ++j)
			mat(i, j) = rand<t_scalar>();

	return mat;
}


/**
 * row permutation matrix
 * @see https://en.wikipedia.org/wiki/Permutation_matrix
 */
template<class t_mat>
t_mat perm(std::size_t N1, std::size_t N2, std::size_t from, std::size_t to)
requires is_basic_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	unit<t_mat>(mat, 0,0, mat.size1(),mat.size2());

	mat(from, from) = mat(to, to) = 0;
	mat(from, to) = mat(to, from) = 1;

	return mat;
}


/**
 * create matrix representation of a permutation
 * @see https://en.wikipedia.org/wiki/Permutation_matrix
 */
template<class t_mat, class t_perm=std::vector<std::size_t>>
t_mat perm(const t_perm& perm)
requires is_basic_mat<t_mat>
{
	auto N = perm.size();
	t_mat mat = zero<t_mat>(N, N);

	for(decltype(N) i = 0; i < perm.size(); ++i)
		mat(i, perm[i]) = mat(perm[i], i) = 1;

	return mat;
}


/**
 * diagonal matrix
 */
template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
t_mat diag(const t_vec& vals)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = vals.size();
	t_mat mat = zero<t_mat>(N);

	// static matrix does not necessarily have the required size!
	if constexpr(!tl2::is_dyn_mat<t_mat>)
		assert(mat.size1() == mat.size1() && mat.size1() == N);

	for(std::size_t i = 0; i < std::min(mat.size1(), N); ++i)
		mat(i,i) = vals[i];

	return mat;
}


/**
 * diagonal matrix
 */
template<class t_mat, class t_val = typename t_mat::value_type>
t_mat diag(std::size_t N, const t_val& val)
requires is_basic_mat<t_mat>
{
	t_mat mat = zero<t_mat>(N, N);

	for(std::size_t i = 0; i < N; ++i)
		mat(i, i) = val;

	return mat;
}


/**
 * vector of diagonal matrix elements
 */
template<class t_vec, class t_mat>
t_vec diag_vec(const t_mat& mat)
requires is_vec<t_vec> && is_mat<t_mat>
{
	std::size_t N = std::min(mat.size1(), mat.size2());

	t_vec vec = zero<t_vec>(N);
	for(std::size_t i = 0; i < N; ++i)
		vec[i] = mat(i,i);

	return vec;
}


/**
 * tests for zero scalar
 */
template<class t_real>
bool equals_0(t_real val, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_scalar<t_real>
{
	return tl2::equals<t_real>(val, t_real(0), eps);
}


/**
 * tests for zero complex number
 */
template<class t_cplx, class t_real = typename t_cplx::value_type>
bool equals_0(const t_cplx& val, typename t_cplx::value_type eps =
	std::numeric_limits<typename t_cplx::value_type>::epsilon())
requires is_complex<t_cplx>
{
	return tl2::equals<t_cplx>(val, t_cplx(0, 0), eps);
}


/**
 * tests for zero vector
 */
template<class t_vec>
bool equals_0(const t_vec& vec,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
requires is_basic_vec<t_vec>
{
	return tl2::equals<t_vec>(vec, zero<t_vec>(vec.size()), eps);
}


/**
 * tests for zero matrix
 */
template<class t_mat>
bool equals_0(const t_mat& mat,
	typename t_mat::value_type eps = std::numeric_limits<typename t_mat::value_type>::epsilon())
requires is_mat<t_mat>
{
	return tl2::equals<t_mat>(mat, zero<t_mat>(mat.size1(), mat.size2()), eps);
}


/**
 * tests for symmetric or hermitian matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool is_symm_or_herm(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		for(std::size_t j=i+1; j < mat.size2(); ++j)
		{
			if constexpr(is_complex<t_elem>)
			{
				// not hermitian?
				if(!tl2::equals<t_elem>(mat(i, j), std::conj(mat(j,i)), eps))
					return false;
			}
			else
			{
				// not symmetric?
				if(!tl2::equals<t_elem>(mat(i, j), mat(j,i), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * tests for skew-symmetric or skew-hermitian matrix
 * @see https://en.wikipedia.org/wiki/Skew-Hermitian_matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool is_skew_symm_or_herm(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		for(std::size_t j=i+1; j < mat.size2(); ++j)
		{
			if constexpr(is_complex<t_elem>)
			{
				// not hermitian?
				if(!tl2::equals<t_elem>(mat(i, j), -std::conj(mat(j,i)), eps))
					return false;
			}
			else
			{
				// not skew-symmetric?
				if(!tl2::equals<t_elem>(mat(i, j), -mat(j,i), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * tests for diagonal matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool is_diag(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	using t_elem = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		for(std::size_t j = 0; j < mat.size2(); ++j)
		{
			if(i == j)
				continue;

			if constexpr(is_complex<t_elem>)
			{
				if(!tl2::equals<t_elem>(mat(i, j), t_elem(0, 0), eps))
					return false;
			}
			else
			{
				if(!tl2::equals<t_elem>(mat(i, j), t_elem(0.), eps))
					return false;
			}
		}
	}

	return true;
}


/**
 * transpose matrix
 * WARNING: not possible for static non-square matrix!
 */
template<class t_mat>
t_mat trans(const t_mat& mat)
requires is_mat<t_mat>
{
	using t_idxtype = decltype(mat.size1());

	t_mat mat2;
	if constexpr(is_dyn_mat<t_mat>)
		mat2 = t_mat(mat.size2(), mat.size1());

	for(t_idxtype i = 0; i < mat.size1(); ++i)
		for(t_idxtype j = 0; j < mat.size2(); ++j)
			mat2(j,i) = mat(i, j);

	return mat2;
}


/**
 * split a matrix into a symmetric and a skew-symmetric part
 * @see https://en.wikipedia.org/wiki/Skew-symmetric_matrix
 */
template<class t_mat>
std::pair<t_mat, t_mat> split_symm(const t_mat& mat)
requires is_mat<t_mat>
{
	using namespace tl2_ops;
	const t_mat mat_T = trans<t_mat>(mat);

	return std::make_pair<t_mat, t_mat>(
		0.5 * (mat + mat_T),
		0.5 * (mat - mat_T));
}


// -----------------------------------------------------------------------------
/**
 * set values close to an integer value to that integer
 * scalar version
 */
template<typename t_real>
void set_eps_round(t_real& d, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_scalar<t_real> && (!is_complex<t_real>)
{
	t_real rd = std::round(d);
	if(equals<t_real>(d, rd, eps))
		d = rd;
};


template<typename t_real, typename t_real_2>
void set_eps_round(t_real& d, t_real_2 eps = std::numeric_limits<t_real_2>::epsilon())
requires is_scalar<t_real> && is_scalar<t_real_2> && (!is_complex<t_real> && !is_complex<t_real_2>)
{
	t_real rd = std::round(d);
	if(equals<t_real>(d, rd, eps))
		d = rd;
};


/**
 * set values close to an integer value to that integer
 * complex version
 */
template<typename t_cplx, class t_real = typename t_cplx::value_type>
void set_eps_round(t_cplx& c, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_complex<t_cplx> && is_scalar<t_real>
{
	t_real rd_re = std::round(c.real());
	if(equals<t_real>(c.real(), rd_re, eps))
		c.real(rd_re);

	t_real rd_im = std::round(c.imag());
	if(equals<t_real>(c.imag(), rd_im, eps))
		c.imag(rd_im);
};


/**
 * set values lower than epsilon to zero
 * complex vector version
 */
template<typename t_vec, typename t_val = typename t_vec::value_type>
void set_eps_round(t_vec& vec,
	typename t_val::value_type eps = std::numeric_limits<typename t_val::value_type>::epsilon())
requires is_basic_vec<t_vec> && is_complex<t_val>
{
	for(std::size_t i = 0; i < vec.size(); ++i)
		set_eps_round<t_val, typename t_val::value_type>(vec[i], eps);
};


/**
 * set values close to an integer value to that integer
 * real vector version
 */
template<typename t_vec, typename t_real = typename t_vec::value_type>
void set_eps_round(t_vec& vec, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_basic_vec<t_vec> && (!is_complex<t_real>)
{
	for(std::size_t i = 0; i < vec.size(); ++i)
		set_eps_round<t_real>(vec[i], eps);
};


/**
 * set values close to an integer value to that integer
 * complex matrix version
 */
template<typename t_mat, typename t_val = typename t_mat::value_type>
void set_eps_round(t_mat& mat,
	typename t_val::value_type eps = std::numeric_limits<typename t_val::value_type>::epsilon())
requires is_basic_mat<t_mat> && is_complex<t_val>
{
	for(std::size_t i = 0; i < mat.size1(); ++i)
		for(std::size_t j = 0; j < mat.size2(); ++j)
			set_eps_round<t_val, typename t_val::value_type>(mat(i, j), eps);
};


/**
 * set values close to an integer value to that integer
 * real matrix version
 */
template<typename t_mat, typename t_val = typename t_mat::value_type>
void set_eps_round(t_mat& mat, t_val eps = std::numeric_limits<t_val>::epsilon())
requires is_basic_mat<t_mat> && (!is_complex<t_val>)
{
	for(std::size_t i = 0; i < mat.size1(); ++i)
		for(std::size_t j = 0; j < mat.size2(); ++j)
			set_eps_round<t_val>(mat(i, j), eps);
};


/**
 * set values lower than epsilon to zero
 * scalar version
 */
template<typename t_real>
void set_eps_0(t_real& d, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_scalar<t_real> && (!is_complex<t_real>)
{
	if(std::abs(d) < eps)
		d = t_real(0);
};


template<typename t_real, typename t_real_2>
void set_eps_0(t_real& d, t_real_2 eps = std::numeric_limits<t_real_2>::epsilon())
requires is_scalar<t_real> && is_scalar<t_real_2> && (!is_complex<t_real> && !is_complex<t_real_2>)
{
	if(std::abs(d) < t_real(eps))
		d = t_real(0);
};


/**
 * set values lower than epsilon to zero
 * complex version
 */
template<typename t_cplx, class t_real = typename t_cplx::value_type>
void set_eps_0(t_cplx& c, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_complex<t_cplx> && is_scalar<t_real>
{
	if(std::abs(c.real()) < eps)
		c.real(t_real(0));
	if(std::abs(c.imag()) < eps)
		c.imag(t_real(0));
};


/**
 * set values lower than epsilon to zero
 * complex vector version
 */
template<typename t_vec, typename t_val = typename t_vec::value_type>
void set_eps_0(t_vec& vec,
	typename t_val::value_type eps = std::numeric_limits<typename t_val::value_type>::epsilon())
requires is_basic_vec<t_vec> && is_complex<t_val>
{
	for(std::size_t i = 0; i < vec.size(); ++i)
		set_eps_0<t_val, typename t_val::value_type>(vec[i], eps);
};


/**
 * set values lower than epsilon to zero
 * real vector version
 */
template<typename t_vec, typename t_real = typename t_vec::value_type>
void set_eps_0(t_vec& vec, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_basic_vec<t_vec> && (!is_complex<t_real>)
{
	for(std::size_t i = 0; i < vec.size(); ++i)
		set_eps_0<t_real>(vec[i], eps);
};


/**
 * set values lower than epsilon to zero
 * complex matrix version
 */
template<typename t_mat, typename t_val = typename t_mat::value_type>
void set_eps_0(t_mat& mat,
	typename t_val::value_type eps = std::numeric_limits<typename t_val::value_type>::epsilon())
requires is_basic_mat<t_mat> && is_complex<t_val>
{
	for(std::size_t i = 0; i < mat.size1(); ++i)
		for(std::size_t j = 0; j < mat.size2(); ++j)
			set_eps_0<t_val, typename t_val::value_type>(mat(i, j), eps);
};


/**
 * set values lower than epsilon to zero
 * real matrix version
 */
template<typename t_mat, typename t_val = typename t_mat::value_type>
void set_eps_0(t_mat& mat, t_val eps = std::numeric_limits<t_val>::epsilon())
requires is_basic_mat<t_mat> && (!is_complex<t_val>)
{
	for(std::size_t i = 0; i < mat.size1(); ++i)
		for(std::size_t j = 0; j < mat.size2(); ++j)
			set_eps_0<t_val>(mat(i, j), eps);
};
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
/**
 * prints a matrix in a 'nicely' formatted form
 */
template<class t_mat, class t_elem = typename t_mat::value_type, class t_real = double>
std::ostream& niceprint(std::ostream& ostr, const t_mat& mat,
	t_real eps = 1e-6, unsigned int prec = 6)
requires tl2::is_basic_mat<t_mat> /*&& tl2::is_dyn_mat<t_mat>*/
{
	ostr.precision(prec);

	const std::size_t ROWS = mat.size1();
	const std::size_t COLS = mat.size2();

	for(std::size_t i = 0; i < ROWS; ++i)
	{
		ostr << "(";
		for(std::size_t j = 0; j < COLS; ++j)
		{
			t_elem elem = mat(i, j);
			set_eps_0(elem, eps);

			unsigned int field_width = prec * 1.5 + 4;
			if constexpr(is_complex<t_elem>)
				field_width *= 3;
			ostr << std::setw(field_width) << std::right << elem;
		}
		ostr << ")\n";
	}

	return ostr;
}


/**
 * prints a vector in a 'nicely' formatted form
 */
template<class t_vec, class t_elem = typename t_vec::value_type, class t_real = double>
std::ostream& niceprint(std::ostream& ostr, const t_vec& vec,
	t_real eps = 1e-6, unsigned int prec = 6)
requires tl2::is_basic_vec<t_vec> /*&& tl2::is_dyn_vec<t_vec>*/
{
	ostr.precision(prec);

	const std::size_t ROWS = vec.size();

	for(std::size_t i = 0; i < ROWS; ++i)
	{
		ostr << "(";
		t_elem elem = vec[i];
		set_eps_0(elem, eps);

		unsigned int field_width = prec * 1.5 + 4;
		if constexpr(is_complex<t_elem>)
			field_width *= 3;
		ostr << std::setw(field_width) << std::right << elem;
		ostr << ")\n";
	}

	return ostr;
}
// -----------------------------------------------------------------------------



/**
 * create a vector with given size if it is dynamic
 */
template<class t_vec>
t_vec create(std::size_t size=3)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(size);

	return vec;
}


/**
 * create a matrix with given sizes if it is dynamic
 */
template<class t_mat>
t_mat create(std::size_t size1, std::size_t size2)
requires is_basic_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat{size1, size2};

	return mat;
}


/**
 * create vector from initializer_list
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
t_vec create(const t_cont<typename t_vec::value_type>& lst)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(lst.size());

	auto iterLst = lst.begin();
	auto size = vec.size();
	using local_size_t = std::decay_t<decltype(size)>;
	for(local_size_t i = 0; i < size; ++i)
	{
		if(iterLst != lst.end())
		{
			vec[i] = *iterLst;
			std::advance(iterLst, 1);
		}
		else	// vector larger than given list?
		{
			vec[i] = 0;
		}
	}

	return vec;
}


/**
 * create matrix from nested initializer_lists in columns/rows order
 */
template<class t_mat,
	template<class...> class t_cont_outer = std::initializer_list,
	template<class...> class t_cont = std::initializer_list>
t_mat create_mat(const t_cont_outer<t_cont<typename t_mat::value_type>>& lst)
requires is_mat<t_mat>
{
	const std::size_t iCols = lst.size();
	const std::size_t iRows = lst.begin()->size();

	t_mat mat = unit<t_mat>(iRows, iCols);

	auto iterCol = lst.begin();
	for(std::size_t iCol=0; iCol<iCols; ++iCol)
	{
		auto iterRow = iterCol->begin();
		for(std::size_t iRow=0; iRow<iRows; ++iRow)
		{
			mat(iRow, iCol) = *iterRow;
			std::advance(iterRow, 1);
		}

		std::advance(iterCol, 1);
	}

	return mat;
}


/**
 * create matrix from nested initializer_lists in columns/rows order
 */
template<class t_mat,
	template<class...> class t_cont_outer = std::initializer_list,
	template<class...> class t_cont = std::initializer_list>
t_mat create(const t_cont_outer<t_cont<typename t_mat::value_type>>& lst)
requires is_mat<t_mat>
{
	return create_mat<t_mat, t_cont_outer, t_cont>(lst);
}


/**
 * create matrix from column (or row) vectors
 */
template<class t_mat, class t_vec, template<class...> class t_cont_outer = std::initializer_list>
t_mat create(const t_cont_outer<t_vec>& lst, bool bRow = false)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t iCols = lst.size();
	const std::size_t iRows = lst.begin()->size();

	t_mat mat = unit<t_mat>(iRows, iCols);

	auto iterCol = lst.begin();
	for(std::size_t iCol=0; iCol<iCols; ++iCol)
	{
		for(std::size_t iRow=0; iRow<iRows; ++iRow)
			mat(iRow, iCol) = (*iterCol)[iRow];
		std::advance(iterCol, 1);
	}

	if(bRow) mat = trans<t_mat>(mat);
	return mat;
}


/**
 * create matrix from initializer_list in column/row order
 */
template<class t_mat>
t_mat create_mat(const std::initializer_list<typename t_mat::value_type>& lst)
requires is_mat<t_mat>
{
	const std::size_t N = std::sqrt(lst.size());

	t_mat mat = unit<t_mat>(N, N);

	auto iter = lst.begin();
	for(std::size_t iRow=0; iRow<N; ++iRow)
	{
		for(std::size_t iCol=0; iCol<N; ++iCol)
		{
			mat(iRow, iCol) = *iter;
			std::advance(iter, 1);
		}
	}

	return mat;
}


/**
 * create matrix from initializer_list in column/row order
 */
template<class t_mat>
t_mat create(const std::initializer_list<typename t_mat::value_type>& lst)
requires is_mat<t_mat>
{
	return create_mat<t_mat>(lst);
}


/**
 * linearise a matrix to a vector container
 */
template<class t_vec, class t_mat>
t_vec convert(const t_mat& mat)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	using T_dst = typename t_vec::value_type;
	using t_idx = decltype(mat.size1());

	t_vec vec;
	vec.reserve(mat.size1()*mat.size2());

	for(t_idx iRow=0; iRow<mat.size1(); ++iRow)
		for(t_idx iCol=0; iCol<mat.size2(); ++iCol)
			vec.push_back(T_dst(mat(iRow, iCol)));

	return vec;
}


/**
 * converts matrix containers of different value types
 */
template<class t_mat_dst, class t_mat_src>
t_mat_dst convert(const t_mat_src& mat)
requires is_mat<t_mat_dst> && is_mat<t_mat_src>
{
	using T_dst = typename t_mat_dst::value_type;
	using t_idx = decltype(mat.size1());

	// in case the static size of the destination vector is larger than the source's
	t_idx maxRows = std::max(mat.size1(), t_idx(t_mat_dst{}.size1()));
	t_idx maxCols = std::max(mat.size2(), t_idx(t_mat_dst{}.size2()));

	t_mat_dst matdst = unit<t_mat_dst>(maxRows, maxCols);

	for(t_idx iRow=0; iRow<mat.size1(); ++iRow)
		for(t_idx iCol=0; iCol<mat.size2(); ++iCol)
			matdst(iRow, iCol) = T_dst(mat(iRow, iCol));

	return matdst;
}


/**
 * converts vector containers of different value types
 */
template<class t_vec_dst, class t_vec_src>
t_vec_dst convert(const t_vec_src& vec)
requires is_vec<t_vec_dst> && is_vec<t_vec_src>
{
	using T_dst = typename t_vec_dst::value_type;
	using t_idx = decltype(vec.size());

	t_vec_dst vecdst = create<t_vec_dst>(vec.size());

	for(t_idx i = 0; i < vec.size(); ++i)
		vecdst[i] = T_dst(vec[i]);

	return vecdst;
}


/**
 * converts a container of objects
 */
template<class t_obj_dst, class t_obj_src, template<class...> class t_cont>
t_cont<t_obj_dst> convert(const t_cont<t_obj_src>& src_objs)
requires (is_vec<t_obj_dst> || is_mat<t_obj_dst>) && (is_vec<t_obj_src> || is_mat<t_obj_src>)
{
	t_cont<t_obj_dst> dst_objs;
	dst_objs.reserve(src_objs.size());

	for(const t_obj_src& src_obj : src_objs)
		dst_objs.emplace_back(tl2::convert<t_obj_dst, t_obj_src>(src_obj));

	return dst_objs;
}


/**
 * get a column vector from a matrix
 */
template<class t_mat, class t_vec>
t_vec col(const t_mat& mat, std::size_t col)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(mat.size1());

	for(std::size_t i = 0; i < mat.size1(); ++i)
		vec[i] = mat(i, col);

	return vec;
}


/**
 * get a row vector from a matrix
 */
template<class t_mat, class t_vec>
t_vec row(const t_mat& mat, std::size_t row)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec;
	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(mat.size2());

	auto size = std::min(mat.size2(), vec.size());
	for(std::size_t i = 0; i < std::size_t(size); ++i)
		vec[i] = mat(row, i);

	return vec;
}


/**
 * set a column vector in a matrix
 */
template<class t_mat, class t_vec>
void set_col(t_mat& mat, const t_vec& vec, std::size_t col)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	for(std::size_t i = 0; i < std::min<std::size_t>(mat.size1(), vec.size()); ++i)
		mat(i, col) = vec[i];
}


/**
 * set a row vector in a matrix
 */
template<class t_mat, class t_vec>
void set_row(t_mat& mat, const t_vec& vec, std::size_t row)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	for(std::size_t i = 0; i < std::min<std::size_t>(mat.size1(), vec.size()); ++i)
		mat(row, i) = vec[i];
}


/**
 * reorder matrix columns according to a permutation
 */
template<class t_mat, class t_vec, class t_perm = std::vector<std::size_t>>
t_mat reorder_cols(const t_mat& mat, const t_perm& perm)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_idx = decltype(mat.size1());
	t_mat mat_new = create<t_mat>(mat.size1(), mat.size2());

	for(t_idx col_idx = 0; col_idx < mat.size2(); ++col_idx)
	{
		t_vec colvec = col<t_mat, t_vec>(mat, perm[col_idx]);
		set_col<t_mat, t_vec>(mat_new, colvec, col_idx);
	}

	return mat_new;
}


/**
 * inner product <vec1|vec2>
 */
template<class t_vec>
typename t_vec::value_type inner(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	typename t_vec::value_type val{0};
	auto size = vec1.size();
	using local_size_t = std::decay_t<decltype(size)>;

	for(local_size_t i = 0; i < size; ++i)
	{
		if constexpr(is_complex<typename t_vec::value_type>)
			val += std::conj(vec1[i]) * vec2[i];
		else
			val += vec1[i] * vec2[i];
	}

	return val;
}


/**
 * inner product <vec1|vec2> (without conjugation of complex vector)
 */
template<class t_vec>
typename t_vec::value_type inner_noconj(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	typename t_vec::value_type val{0};
	auto size = vec1.size();
	using local_size_t = std::decay_t<decltype(size)>;

	for(local_size_t i = 0; i < size; ++i)
		val += vec1[i] * vec2[i];

	return val;
}


/**
 * inner product between two vectors of different type
 */
template<class t_vec1, class t_vec2>
typename t_vec1::value_type inner(const t_vec1& vec1, const t_vec2& vec2)
requires is_basic_vec<t_vec1> && is_basic_vec<t_vec2>
{
	if(vec1.size()==0 || vec2.size()==0)
		return typename t_vec1::value_type{};

	// first element
	auto val = vec1[0]*vec2[0];

	// remaining elements
	for(std::size_t i=1; i < std::min(vec1.size(), vec2.size()); ++i)
	{
		if constexpr(is_complex<typename t_vec1::value_type>)
		{
			auto prod = std::conj(vec1[i]) * vec2[i];
			val = val + prod;
		}
		else
		{
			auto prod = vec1[i]*vec2[i];
			val = val + prod;
		}
	}

	return val;
}


/**
 * matrix-matrix product: c_ij = a_ik b_kj
 */
template<class t_mat>
t_mat prod(const t_mat& mat1, const t_mat& mat2, bool assert_sizes/*=true*/)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>
{
	// if not asserting sizes, the inner size will use the minimum of the two matrix sizes
	if(assert_sizes)
	{
		if constexpr(tl2::is_dyn_mat<t_mat>)
			assert((mat1.size2() == mat2.size1()));
		else
			static_assert(t_mat::size2() == t_mat::size1());
	}


	t_mat matRet(mat1.size1(), mat2.size2());
	const std::size_t innersize = std::min(mat1.size2(), mat2.size1());

	for(std::size_t row=0; row<matRet.size1(); ++row)
	{
		for(std::size_t col=0; col<matRet.size2(); ++col)
		{
			matRet(row, col) = 0;
			for(std::size_t i = 0; i < innersize; ++i)
				matRet(row, col) += mat1(row, i) * mat2(i, col);
		}
	}

	return matRet;
}


/**
 * element-wise division
 */
template<class t_mat>
t_mat div_perelem(const t_mat& mat1, const t_mat& mat2, bool assert_sizes=true)
requires tl2::is_basic_mat<t_mat>
{
	if(assert_sizes)
	{
		if constexpr(tl2::is_dyn_mat<t_mat>)
			assert(mat1.size1() == mat2.size1() && mat1.size2() == mat2.size2());
	}


	t_mat matRet = zero<t_mat>(mat1.size1(), mat1.size2());

	for(std::size_t row=0; row<matRet.size1(); ++row)
		for(std::size_t col=0; col<matRet.size2(); ++col)
			matRet(row, col) = mat1(row, col) / mat2(row,col);

	return matRet;
}


/**
 * 2-norm
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_real norm(const t_vec& vec)
requires is_basic_vec<t_vec>
{
	t_real d = static_cast<t_real>(tl2::inner<t_vec>(vec, vec));
	return std::sqrt(d);
}


/**
 * n-norm
 */
template<class t_vec, class t_real = typename t_vec::value_type>
typename t_vec::value_type norm(const t_vec& vec, t_real n)
requires is_basic_vec<t_vec>
{
	t_real d = t_real{0};
	for(std::size_t i = 0; i < vec.size(); ++i)
		d += std::pow(std::abs(vec[i]), n);
	return std::pow(d, t_real(1)/n);
}


/**
 * outer product
 */
template<class t_mat, class t_vec>
t_mat outer(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	const std::size_t N1 = vec1.size();
	const std::size_t N2 = vec2.size();

	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(std::size_t n1=0; n1<N1; ++n1)
	{
		for(std::size_t n2=0; n2<N2; ++n2)
		{
			if constexpr(is_complex<typename t_vec::value_type>)
				mat(n1, n2) = std::conj(vec1[n1]) * vec2[n2];
			else
				mat(n1, n2) = vec1[n1]*vec2[n2];
		}
	}

	return mat;
}


/**
 * outer product (without conjugation of complex vector)
 */
template<class t_mat, class t_vec>
t_mat outer_noconj(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	const std::size_t N1 = vec1.size();
	const std::size_t N2 = vec2.size();

	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(N1, N2);

	for(std::size_t n1=0; n1<N1; ++n1)
		for(std::size_t n2=0; n2<N2; ++n2)
			mat(n1, n2) = vec1[n1]*vec2[n2];

	return mat;
}


/**
 * submatrix removing a column/row from a matrix stored in a vector container
 */
template<class t_vec, class t_matvec = t_vec>
t_vec flat_submat(const t_matvec& mat,
	std::size_t iNumRows, std::size_t iNumCols,
	std::size_t iRemRow, std::size_t iRemCol)
requires is_basic_vec<t_vec>
{
	t_vec vec;
	vec.reserve(mat.size());

	for(std::size_t iRow=0; iRow<iNumRows; ++iRow)
	{
		if(iRow == iRemRow)
			continue;

		for(std::size_t iCol=0; iCol<iNumCols; ++iCol)
		{
			if(iCol == iRemCol)
				continue;
			vec.push_back(mat[iRow*iNumCols + iCol]);
		}
	}

	return vec;
}


/**
 * submatrix removing a column/row from a matrix
 */
template<class t_mat>
t_mat submat(const t_mat& mat, decltype(mat.size1()) iRemRow, decltype(mat.size2()) iRemCol)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());
	t_mat matRet = create<t_mat>(mat.size1()-1, mat.size2()-1);

	size_t iResRow = 0;
	for(size_t iRow=0; iRow<mat.size1(); ++iRow)
	{
		if(iRow == iRemRow)
			continue;

		size_t iResCol = 0;
		for(size_t iCol=0; iCol<mat.size2(); ++iCol)
		{
			if(iCol == iRemCol)
				continue;

			matRet(iResRow, iResCol) = mat(iRow, iCol);
			++iResCol;
		}

		++iResRow;
	}

	return matRet;
}


/**
 * submatrix of a given size starting at given indices
 */
template<class t_mat>
t_mat submat(const t_mat& mat,
	decltype(mat.size1()) row_start, decltype(mat.size2()) col_start,
	decltype(mat.size1()) num_rows, decltype(mat.size2()) num_cols)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());
	t_mat matRet = create<t_mat>(num_rows, num_cols);

	for(size_t row=0; row<num_rows; ++row)
		for(size_t col=0; col<num_cols; ++col)
			matRet(row, col) = mat(row+row_start, col+col_start);

	return matRet;
}


/**
 * set submatrix at given starting indices
 */
template<class t_mat>
void set_submat(t_mat& mat, const t_mat& submat,
	decltype(mat.size1()) row_start, decltype(mat.size2()) col_start,
	decltype(mat.size1()) num_rows=9999, decltype(mat.size2()) num_cols=9999)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());

	num_rows = std::min(num_rows, submat.size1());
	num_cols = std::min(num_cols, submat.size2());

	for(size_t row=0; row<num_rows; ++row)
		for(size_t col=0; col<num_cols; ++col)
			mat(row+row_start, col+col_start) = submat(row, col);
}


/**
 * add submatrix at given starting indices
 */
template<class t_mat>
void add_submat(t_mat& mat, const t_mat& submat,
	decltype(mat.size1()) row_start, decltype(mat.size2()) col_start)
requires is_basic_mat<t_mat> && is_dyn_mat<t_mat>
{
	using size_t = decltype(mat.size1());

	for(size_t row=0; row<submat.size1(); ++row)
	{
		for(size_t col=0; col<submat.size2(); ++col)
			mat(row+row_start, col+col_start) += submat(row, col);
	}
}


/**
 * subvector removing a row from a vector
 */
template<class t_vec>
t_vec subvec(const t_vec& vec, decltype(vec.size()) iRemRow)
requires is_basic_vec<t_vec> && is_dyn_mat<t_vec>
{
	using size_t = decltype(vec.size());
	t_vec vecRet = create<t_vec>(vec.size() - 1);

	size_t iResRow = 0;
	for(size_t iRow = 0; iRow < vec.size(); ++iRow)
	{
		if(iRow == iRemRow)
			continue;

		vecRet[iResRow] = vec[iRow];
		++iResRow;
	}

	return vecRet;
}


/**
 * determinant from a square matrix stored in a vector container
 * @see (Merziger 2006), p. 185
 */
template<class t_vec, class t_matvec = t_vec>
typename t_vec::value_type flat_det(const t_matvec& mat, std::size_t iN)
requires is_basic_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// special cases
	if(iN == 0)
		return 0;
	else if(iN == 1)
		return mat[0];
	else if(iN == 2)
		return mat[0]*mat[3] - mat[1]*mat[2];

	// recursively expand determiant along a row
	T fullDet = T(0);
	std::size_t iRow = 0;

	// get row with maximum number of zeros
	std::size_t iMaxNumZeros = 0;
	for(std::size_t iCurRow = 0; iCurRow < iN; ++iCurRow)
	{
		std::size_t iNumZeros = 0;
		for(std::size_t iCurCol = 0; iCurCol < iN; ++iCurCol)
		{
			if(equals<T>(mat[iCurRow*iN + iCurCol], T(0)))
				++iNumZeros;
		}

		if(iNumZeros > iMaxNumZeros)
		{
			iRow = iCurRow;
			iMaxNumZeros = iNumZeros;
		}
	}

	for(std::size_t iCol = 0; iCol < iN; ++iCol)
	{
		const T elem = mat[iRow*iN + iCol];
		if(equals<T>(elem, 0))
			continue;

		const T sgn = ((iRow+iCol) % 2) == 0 ? T(1) : T(-1);
		const t_vec subMat = tl2::flat_submat<t_vec, t_matvec>(
			mat, iN, iN, iRow, iCol);
		const T subDet = tl2::flat_det<t_vec>(subMat, iN-1) * sgn;

		fullDet += elem * subDet;
	}

	return fullDet;
}


/**
 * trace
 */
template<class t_mat>
typename t_mat::value_type trace(const t_mat& mat)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	T _tr = T(0);

	std::size_t N = std::min(mat.size1(), mat.size2());
	for(std::size_t i = 0; i < N; ++i)
		_tr += mat(i, i);

	return _tr;
}


/**
 * sum of all matrix elements
 */
template<class t_mat>
typename t_mat::value_type sum(const t_mat& mat)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	T s = T(0);

	for(std::size_t i = 0; i < mat.size1(); ++i)
		for(std::size_t j = 0; j < mat.size2(); ++j)
			s += mat(i, j);

	return s;
}


/**
 * gets reciprocal basis vectors |b_i> from real basis vectors |a_i> (and vice versa)
 * c: multiplicative constant (c=2*pi for physical lattices, c=1 for mathematics)
 *
 * Def.: <b_i | a_j> = c * delta(i, j)  =>
 *
 * e.g. 2d case:
 *                   ( a_1x  a_2x )
 *                   ( a_1y  a_2y )
 *
 * ( b_1x  b_1y )    (    1     0 )
 * ( b_2x  b_2y )    (    0     1 )
 *
 * B^t * A = I
 * A = B^(-t)
 */
template<class t_mat, class t_vec,
	template<class...> class t_cont_in = std::initializer_list,
	template<class...> class t_cont_out = std::vector>
t_cont_out<t_vec> recip(const t_cont_in<t_vec>& lstReal, typename t_vec::value_type c=1)
requires is_mat<t_mat> && is_basic_vec<t_vec>
{
	const t_mat basis = create<t_mat, t_vec, t_cont_in>(lstReal);
	auto [basis_inv, bOk] = inv<t_mat>(basis);
	basis_inv *= c;

	t_cont_out<t_vec> lstRecip;
	lstRecip.reserve(basis_inv.size1());

	for(std::size_t currow=0; currow<basis_inv.size1(); ++currow)
	{
		const t_vec rowvec = row<t_mat, t_vec>(basis_inv, currow);
		lstRecip.emplace_back(std::move(rowvec));
	}

	return lstRecip;
}


template<class t_vec>
t_vec cross(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>;


/**
 * general n-dim cross product using determinant definition
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec, template<class...> class t_cont = std::initializer_list>
t_vec cross(const t_cont<t_vec>& vecs)
requires is_basic_vec<t_vec>
{
	using t_size = decltype(t_vec{}.size());
	using T = typename t_vec::value_type;

	// N also has to be equal to the vector size!
	const t_size N = vecs.size()+1;
	t_vec vec = zero<t_vec>(N);

	// 3-dim case
	if(N == 3 && vecs.begin()->size() == 3)
	{
		return tl2::cross<t_vec>(*vecs.begin(), *std::next(vecs.begin(), 1));
	}

	// general case
	else
	{
		for(t_size iComp=0; iComp<N; ++iComp)
		{
			std::vector<T> mat = zero<std::vector<T>>(N*N);
			mat[0*N + iComp] = T(1);

			t_size row_idx = 0;
			for(const t_vec& vec : vecs)
			{
				for(t_size col_idx=0; col_idx<N; ++col_idx)
					mat[(row_idx+1)*N + col_idx] = vec[col_idx];
				++row_idx;
			}

			vec[iComp] = tl2::flat_det<decltype(mat)>(mat, N);
		}
	}

	return vec;
}



// ----------------------------------------------------------------------------
/**
 * QR decomposition of a matrix
 * @returns [ok, Q, R]
 * @see (Scarpino 2011), pp. 269-272
 */
template<class t_mat, class t_vec>
std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>
{
#ifdef __TLIBS2_USE_LAPACK__

	return tl2_la::qr<t_mat, t_vec>(mat);

#else

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();
	const std::size_t N = std::min(cols, rows);

#if __TLIBS2_QR_METHOD__ == 0
	t_mat R = mat;
	t_mat Q = unit<t_mat>(N, N);

	for(std::size_t icol = 0; icol < N-1; ++icol)
	{
		t_vec vecCol = col<t_mat, t_vec>(R, icol);
		t_mat matMirror = ortho_mirror_zero_op<t_mat, t_vec>(vecCol, icol);

		Q = prod(Q, matMirror, false);
		R = prod(matMirror, R);
	}

#elif __TLIBS2_QR_METHOD__ == 1
	std::vector<t_vec> sysM;
	sysM.reserve(mat.size2());
	for(std::size_t i = 0; i < mat.size2(); ++i)
		sysM.push_back(col<t_mat, t_vec>(mat, i));

	std::vector<t_vec> Qsys = orthonorm_sys<t_vec, std::vector, std::vector>(sysM);
	t_mat Q = create<t_mat>(Qsys, false);
	t_mat R = prod(trans(Q), mat);
#endif

	return std::make_tuple(true, Q, R);

#endif  // __TLIBS2_USE_LAPACK__
}


/**
 * inverted matrix
 * @see https://en.wikipedia.org/wiki/Invertible_matrix#In_relation_to_its_adjugate
 * @see https://en.wikipedia.org/wiki/Adjugate_matrix
 */
template<class t_mat>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires is_mat<t_mat>
{
	// fails if matrix is not square, TODO: fix
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

#ifdef __TLIBS2_USE_LAPACK__

	return tl2_la::inv<t_mat>(mat);

#else

	using T = typename t_mat::value_type;
	using t_vec = std::vector<T>;
	const std::size_t N = mat.size1();

	const auto& matFlat = matvec_adapter<t_mat>{mat};
	const T fullDet = tl2::flat_det<t_vec>(matFlat, N);

	// fail if determinant is zero
	if(equals<T>(fullDet, 0))
		return std::make_tuple(t_mat(), false);

	t_mat matInv;
	if constexpr(is_dyn_mat<t_mat>)
		matInv = t_mat(N, N);

	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j < N; ++j)
		{
			const T sgn = ((i+j) % 2) == 0 ? T(1) : T(-1);
			const t_vec subMat = tl2::flat_submat<t_vec>(matFlat, N, N, i, j);
			matInv(j,i) = sgn * tl2::flat_det<t_vec>(subMat, N-1);
		}
	}

	matInv = matInv / fullDet;
	return std::make_tuple(matInv, true);

#endif
}


/**
 * determinant
 *
 * M v = l v
 * M = V^-1 L V
 * det(M) = det(V^-1 L V) = det(V V^-1 L) = det(L) = l_1*...*l_n
 */
template<class t_mat>
typename t_mat::value_type det(const t_mat& mat)
requires is_mat<t_mat>
{
	using T = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return 0;

#ifdef __TLIBS2_USE_LAPACK__

	using t_vec = vec<T>;

	if constexpr(tl2::is_complex<typename t_mat::value_type>)
	{
		const auto [ok, evals, evecs] =
			tl2_la::eigenvec<t_mat, t_vec, T>(mat, true, false, false);

		T detval{1, 0};
		for(const auto& eval : evals)
			detval *= eval;

		return detval;
	}
	else
	{
		const auto [ok, evalsRe, evalsIm, evecsRe, evecsIm] =
			tl2_la::eigenvec<t_mat, t_vec, T>(mat, true, false, false);

		std::complex<T> detval{1, 0};
		for(std::size_t i = 0; i < evalsRe.size(); ++i)
			detval *= std::complex<T>{evalsRe[i], evalsIm[i]};

		return detval.real();
	}

#else

	const auto& matFlat = matvec_adapter<t_mat>{mat};
	return tl2::flat_det<std::vector<T>>(matFlat, mat.size1());

#endif
}


/**
 * matrix power
 */
template<class t_mat>
std::tuple<t_mat, bool> pow(const t_mat& mat, int ipow)
requires is_mat<t_mat>
{
	t_mat themat;
	if(mat.size1() != mat.size2())
		return std::make_tuple(themat, false);

	bool ok = true;
	int ipow_pos = ipow<0 ? -ipow : ipow;

	themat = unit<t_mat>(mat.size1());
	for(int i = 0; i < ipow_pos; ++i)
		themat = themat*mat;

	if(ipow < 0)
		std::tie(themat, ok) = tl2::inv<t_mat>(themat);

	return std::make_tuple(themat, ok);
}


/**
 * least-squares regression, solves for parameters v
 * @see https://en.wikipedia.org/wiki/Least_squares
 * @see (Arens 2015), p. 793
 *
 * exact equation:
 * 	X v = y
 *
 * approx. equation (normal equation):
 * 	X^t X v = X^t y
 * 	v = inv(X^t X) X^t y
 *
 * 	(QR)^t QR v = X^t y
 * 	R^tQ^t QR v = X^t y
 * 	R^tR v = X^t y
 * 	v = inv(R^tR) X^t y
 */
template<class t_vec, class t_mat = mat<typename t_vec::value_type, std::vector>>
std::tuple<t_vec, bool> leastsq(const t_vec& x, const t_vec& y, std::size_t order,
	bool use_qr = true, bool use_pseudoinv = false)
requires is_vec<t_vec> && is_dyn_mat<t_mat>
{
	// check array sizes, TODO: fix
	if constexpr(tl2::is_dyn_vec<t_vec>)
		assert((x.size() == y.size()));


	using namespace tl2_ops;
	using T = typename t_vec::value_type;

	const std::size_t N = x.size();
	t_mat X{N, order+1};

	for(std::size_t j = 0; j <= order; ++j)
		for(std::size_t i = 0; i < N; ++i)
			X(i, j) = std::pow(x[i], static_cast<T>(j));


	t_mat XtX;

	if(use_qr)
	{
		auto [ok_qr, Q, R] = qr<t_mat, t_vec>(X);
		if(!ok_qr)
				return std::make_tuple(t_vec{}, false);

		XtX = trans<t_mat>(R) * R;
	}
	else
	{
		XtX = trans<t_mat>(X) * X;
	}

	t_vec Xty = trans<t_mat>(X) * y;
	t_mat Y;
	bool ok = false;

	if(use_pseudoinv)
	{
#ifdef __TLIBS2_USE_LAPACK__
		std::tie(Y, ok) = tl2_la::pseudoinv<t_mat>(XtX);
#else
		//#pragma message("leastsq: Pseudo-inverse is not available, using standard inverse instead.")
		std::tie(Y, ok) = inv<t_mat>(XtX);
#endif
	}
	else
	{
		std::tie(Y, ok) = inv<t_mat>(XtX);
	}

	t_vec v = Y * Xty;

	return std::make_tuple(v, ok);
}
// ----------------------------------------------------------------------------

}

#endif
