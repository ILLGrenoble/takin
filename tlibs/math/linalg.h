/**
 * basic linalg helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2018
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

#ifndef __TLIBS_LINALG_H__
#define __TLIBS_LINALG_H__

#include "../helper/boost_hacks.h"
#include "../helper/exception.h"
#include "math.h"
#include "../log/log.h"
#include "../log/debug.h"
#include "../helper/traits.h"

#include <initializer_list>
#include <cmath>

#include <boost/algorithm/minmax_element.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/exception.hpp>


namespace tl {

namespace ublas = boost::numeric::ublas;


template<class matrix_type = ublas::matrix<double>>
typename matrix_type::value_type determinant(const matrix_type& mat);


// ----------------------------------------------------------------------------

/**
 * creates a vector
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_lst = std::initializer_list>
t_vec make_vec(t_lst<typename t_vec::value_type>&& lst)
{
	using T = typename t_vec::value_type;
	using t_iter = typename t_lst<T>::const_iterator;

	t_vec vec(lst.size());

	std::size_t i=0;
	for(t_iter iter=lst.begin(); iter!=lst.end(); ++i, ++iter)
		vec[i] = std::move(*iter);

	return vec;
}


/**
 * creates a vector
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_lst = std::initializer_list>
t_vec make_vec(const t_lst<typename t_vec::value_type>& lst)
{
	using T = typename t_vec::value_type;
	using t_iter = typename t_lst<T>::const_iterator;

	t_vec vec(lst.size());

	std::size_t i=0;
	for(t_iter iter=lst.begin(); iter!=lst.end(); ++i, ++iter)
		vec[i] = *iter;

	return vec;
}


/**
 * creates a matrix
 */
template<class t_mat = ublas::matrix<double>,
	template<class...> class t_lst = std::initializer_list>
t_mat make_mat(t_lst<t_lst<typename t_mat::value_type>>&& lst)
{
	using T = typename t_mat::value_type;

	const std::size_t I = lst.size();
	const std::size_t J = lst.begin()->size();

	t_mat mat(I, J);
	typename t_lst<t_lst<T>>::const_iterator iter = lst.begin();

	for(std::size_t i=0; i<I; ++i, ++iter)
	{
		typename t_lst<T>::const_iterator iterinner = iter->begin();
		for(std::size_t j=0; j<J; ++j, ++iterinner)
		{
			mat(i,j) = std::move(*iterinner);
		}
	}

	return mat;
}


/**
 * unit matrix -- general version
 */
template<class t_mat = ublas::matrix<double>,
	typename std::enable_if<!std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat unit_m(std::size_t N)
{
	t_mat mat(N, N);
	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			mat(i,j) = (i==j ? 1 : 0);
	return mat;
}


/**
 * unit matrix -- ublas wrapper
 */
template<class t_mat = ublas::matrix<double>,
typename std::enable_if<std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat unit_m(std::size_t N)
{
	return ublas::identity_matrix<typename t_mat::value_type>(N);
}


/**
 * zero matrix -- general version
 */
template<class t_mat = ublas::matrix<double>,
typename std::enable_if<!std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat zero_m(std::size_t N, std::size_t M)
{
	t_mat mat(N, M);
	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<M; ++j)
			mat(i,j) = 0;
	return mat;
}


/**
 * zero matrix -- ublas wrapper
 */
template<class t_mat = ublas::matrix<double>,
typename std::enable_if<std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat zero_m(std::size_t N, std::size_t M)
{
	return ublas::zero_matrix<typename t_mat::value_type>(N, M);
}


/**
 * zero matrix -- synonym
 */
template<class t_mat = ublas::matrix<double>>
t_mat zero_matrix(std::size_t N, std::size_t M)
{ return zero_m<t_mat>(N, M); }


/**
 * zero vector -- general version
 */
template<class t_vec = ublas::vector<double>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_vec zero_v(std::size_t N)
{
	t_vec vec(N);
	for(std::size_t i=0; i<N; ++i)
		vec(i) = 0;
	return vec;
}


/**
 * zero vector -- ublas wrapper
 */
template<class t_vec = ublas::vector<double>,
	typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_vec zero_v(std::size_t N)
{
	return ublas::zero_vector<typename t_vec::value_type>(N);
}

/**
 * zero vector -- synonym
 */
template<class t_vec = ublas::vector<double>>
t_vec zero_vector(std::size_t N)
{ return zero_v<t_vec>(N); }



/**
 * create a vector of size N filled with value val
 */
template<class vector_type = ublas::vector<double>>
vector_type fill_vector(std::size_t N, typename vector_type::value_type val)
{
	vector_type vec(N);
	for(std::size_t i=0; i<N; ++i)
		vec[i] = val;
	return vec;
}


/**
 * resize matrix, filling up with unity
 */
template<class t_mat = ublas::matrix<double>>
void resize_unity(t_mat& mat, std::size_t N)
{
	const std::size_t iOldSize1 = mat.size1();
	const std::size_t iOldSize2 = mat.size2();

	mat.resize(N,N, true);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
		{
			if(i<iOldSize1 && j<iOldSize2) continue;
			mat(i,j) = (i==j ? 1 : 0);
		}
}


/**
 * converts vector t_vec<t_from> to t_vec<t_to>
 */
template<class t_from, class t_to, template<class...> class t_vec = ublas::vector>
t_vec<t_to> convert_vec(const t_vec<t_from>& vec)
{
	//using t_vec_from = t_vec<t_from>;
	using t_vec_to = t_vec<t_to>;

	t_vec_to vecRet(vec.size());

	for(std::size_t i=0; i<vec.size(); ++i)
		vecRet[i] = t_to(vec[i]);

	return vecRet;
}

/**
 * converts vector t_vec_from<t_from> to t_vec_to<t_to>
 */
template<class t_from, class t_to,
	template<class...> class t_vec_from = std::vector,
	template<class...> class t_vec_to = ublas::vector>
t_vec_to<t_to> convert_vec_full(const t_vec_from<t_from>& vec)
{
	t_vec_to<t_to> vecRet(vec.size());

	for(std::size_t i=0; i<vec.size(); ++i)
		vecRet[i] = t_to(vec[i]);

	return vecRet;
}

// ----------------------------------------------------------------------------


template<class vec_type>
bool vec_equal(const vec_type& vec0, const vec_type& vec1,
	typename tl::_get_epsilon_impl<vec_type>::t_eps eps = get_epsilon<vec_type>())
{
	typedef typename vec_type::value_type T;

	if(vec0.size() != vec1.size())
		return false;

	for(std::size_t i=0; i<vec0.size(); ++i)
		if(!float_equal<T>(vec0[i], vec1[i], eps))
			return false;
	return true;
}


template<class mat_type>
bool mat_equal(const mat_type& mat0, const mat_type& mat1,
	typename tl::_get_epsilon_impl<mat_type>::t_eps eps = get_epsilon<mat_type>())
{
	typedef typename mat_type::value_type T;

	if(mat0.size1() != mat1.size1() || mat0.size2() != mat1.size2())
		return false;

	for(std::size_t i=0; i<mat0.size1(); ++i)
		for(std::size_t j=0; j<mat0.size2(); ++j)
			if(!float_equal<T>(mat0(i,j), mat1(i,j), eps))
				return false;
	return true;
}


// ----------------------------------------------------------------------------


/**
 * transpose -- general version
 */
template<typename t_mat = ublas::matrix<double>,
	typename std::enable_if<!std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat transpose(const t_mat& mat)
{
	t_mat matret(mat.size2(), mat.size1());

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=0; j<mat.size2(); ++j)
			matret(j,i) = mat(i,j);

	return matret;
}


/**
 * transpose -- general version
 */
template<typename t_mat = ublas::matrix<double>,
	typename std::enable_if<std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat transpose(const t_mat& mat)
{
	return ublas::trans(mat);
}


/**
 * conjugate matrix
 */
template<typename t_mat = ublas::matrix<std::complex<double>>>
t_mat conjugate_mat(t_mat mat)
{
	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=0; j<mat.size2(); ++j)
			mat(i,j) = std::conj(mat(i,j));
	return mat;
}


/**
 * hermitian conjugated matrix
 */
template<typename t_mat = ublas::matrix<std::complex<double>>>
t_mat hermitian(const t_mat& mat)
{
	t_mat matret = transpose<t_mat>(mat);
	matret = conjugate_mat<t_mat>(matret);
	return matret;
}


/**
 * conjugate vector
 */
template<typename t_vec = ublas::vector<std::complex<double>>>
t_vec conjugate_vec(t_vec vec)
{
	for(std::size_t i=0; i<vec.size(); ++i)
		vec[i] = std::conj(vec[i]);
	return vec;
}


// ----------------------------------------------------------------------------


/**
 * cross product, c_i = eps_ijk a_j b_k
 */
template<typename vector_type = ublas::vector<double>>
vector_type cross_3(const vector_type& vec0, const vector_type& vec1)
{
	return make_vec<vector_type>
	({
		vec0[1]*vec1[2] - vec1[1]*vec0[2],
		vec0[2]*vec1[0] - vec1[2]*vec0[0],
		vec0[0]*vec1[1] - vec1[0]*vec0[1]
	});
}


/**
 * inner product -- general version
 */
template<typename t_vec = ublas::vector<double>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
typename t_vec::value_type inner(const t_vec& vec0, const t_vec& vec1)
{
	typename t_vec::value_type d(0);
	std::size_t iSize = std::min(vec0.size(), vec1.size());

	for(std::size_t i=0; i<iSize; ++i)
		d += vec0[i]*vec1[i];
	return d;
}


/**
 * complex inner product -- general version
 */
template<typename t_vec = ublas::vector<std::complex<double>>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
typename t_vec::value_type inner_cplx(const t_vec& vec0, const t_vec& vec1)
{
	typename t_vec::value_type d(0);
	std::size_t iSize = std::min(vec0.size(), vec1.size());

	for(std::size_t i=0; i<iSize; ++i)
		d += std::conj(vec0[i])*vec1[i];
	return d;
}


/**
 * inner product -- ublas wrapper
 */
template<typename t_vec = ublas::vector<double>,
	typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
typename t_vec::value_type inner(const t_vec& vec0, const t_vec& vec1)
{
	return ublas::inner_prod(vec0, vec1);
}


/**
 * complex inner product -- ublas wrapper
 */
template<typename t_vec = ublas::vector<std::complex<double>>,
	typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
typename t_vec::value_type inner_cplx(const t_vec& vec0, const t_vec& vec1)
{
	t_vec vec0_c = conjugate_vec<t_vec>(vec0);
	return ublas::inner_prod(vec0_c, vec1);
}


/**
 * outer product -- general version
 */
template<typename t_vec = ublas::vector<double>,
	typename t_mat = ublas::matrix<typename t_vec::value_type>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_mat outer(const t_vec& vec0, const t_vec& vec1)
{
	std::size_t iSize = std::min(vec0.size(), vec1.size());
	t_mat mat(iSize, iSize);

	for(std::size_t i=0; i<iSize; ++i)
		for(std::size_t j=0; j<iSize; ++j)
			mat(i,j) = vec0[i]*vec1[j];

	return mat;
}


/**
 * complex outer product -- general version
 */
template<typename t_vec = ublas::vector<std::complex<double>>,
	typename t_mat = ublas::matrix<typename t_vec::value_type>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_mat outer_cplx(const t_vec& vec0, const t_vec& vec1)
{
	std::size_t iSize = std::min(vec0.size(), vec1.size());
	t_mat mat(iSize, iSize);

	for(std::size_t i=0; i<iSize; ++i)
		for(std::size_t j=0; j<iSize; ++j)
			mat(i,j) = vec0[i]*std::conj(vec1[j]);

	return mat;
}


/**
 * outer product -- ublas wrapper
 */
template<typename t_vec = ublas::vector<double>,
	typename t_mat = ublas::matrix<typename t_vec::value_type>,
	typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_mat outer(const t_vec& vec0, const t_vec& vec1)
{
	return ublas::outer_prod(vec0, vec1);
}


/**
 * complex outer product -- ublas wrapper
 */
template<typename t_vec = ublas::vector<std::complex<double>>,
	typename t_mat = ublas::matrix<typename t_vec::value_type>,
	typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_mat outer_cplx(const t_vec& vec0, const t_vec& vec1)
{
	t_vec vec1_c = conjugate_vec<t_vec>(vec1);
	return ublas::outer_prod(vec0, vec1_c);
}


/**
 * matrix-matrix product -- general version
 * c_ij = a_ik b_kj
 */
template<typename t_mat = ublas::matrix<double>,
typename std::enable_if<!std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat prod_mm(const t_mat& mat0, const t_mat& mat1)
{
	if(mat0.size2() != mat1.size1())
		return t_mat(0,0);

	std::size_t iSize1 = mat0.size1();
	std::size_t iSize2 = mat1.size2();
	std::size_t iSize3 = mat0.size2();
	t_mat mat(iSize1, iSize2);

	for(std::size_t i=0; i<iSize1; ++i)
	{
		for(std::size_t j=0; j<iSize2; ++j)
		{
			mat(i,j) = typename t_mat::value_type(0);
			for(std::size_t k=0; k<iSize3; ++k)
				mat(i,j) += mat0(i,k)*mat1(k,j);
		}
	}

	return mat;
}


/**
 * matrix-matrix product -- ublas wrapper
 */
template<typename t_mat = ublas::matrix<double>,
typename std::enable_if<std::is_convertible<t_mat, ublas::matrix<typename t_mat::value_type>>::value, char>::type=0>
t_mat prod_mm(const t_mat& mat0, const t_mat& mat1)
{
	return ublas::prod(mat0, mat1);
}


/**
 * matrix-vector product -- general version
 * c_i = a_ij b_j
 */
template<typename t_vec = ublas::vector<double>,
typename t_mat = ublas::matrix<typename t_vec::value_type>,
typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_vec prod_mv(const t_mat& mat, const t_vec& vec)
{
	if(mat.size2() != vec.size())
		return t_vec(0);

	std::size_t iSize1 = mat.size1();
	std::size_t iSize2 = mat.size2();
	t_vec vecret(iSize2);

	for(std::size_t i=0; i<iSize1; ++i)
	{
		vecret(i) = typename t_vec::value_type(0);
		for(std::size_t j=0; j<iSize2; ++j)
			vecret(i) += mat(i,j)*vec(j);
	}

	return vecret;
}


/**
 * matrix-vector product -- ublas wrapper
 */
template<typename t_vec = ublas::vector<double>,
typename t_mat = ublas::matrix<typename t_vec::value_type>,
typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_vec prod_mv(const t_mat& mat, const t_vec& vec)
{
	return ublas::prod(mat, vec);
}


/**
 * vector-matrix product
 */
template<typename t_vec = ublas::vector<double>,
typename t_mat = ublas::matrix<typename t_vec::value_type>,
typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
t_vec prod_vm(const t_vec& vec, const t_mat& mat)
{
	return prod_mv(transpose(mat), vec);
}


/**
 * 2-norm -- general real version
 */
template<class t_vec = ublas::vector<double>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value
		&& std::is_floating_point<typename t_vec::value_type>::value, char>::type=0>
typename t_vec::value_type veclen(const t_vec& vec)
{
	using T = typename t_vec::value_type;
	T len(0);

	for(std::size_t i=0; i<vec.size(); ++i)
		len += vec[i]*vec[i];

	return std::sqrt(len);
}


/**
 * 2-norm -- general complex version
 */
template<class t_vec = ublas::vector<double>,
	typename std::enable_if<!std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value
		&& !std::is_floating_point<typename t_vec::value_type>::value, char>::type=0>
typename t_vec::value_type::value_type veclen(const t_vec& vec)
{
        using T = typename t_vec::value_type::value_type;
        T len(0);

        for(std::size_t i=0; i<vec.size(); ++i)
                len += (std::conj(vec[i])*vec[i]).real();

        return std::sqrt(len);
}


/**
 * 2-norm -- ublas wrapper
 */
template<class t_vec = ublas::vector<double>,
	typename std::enable_if<std::is_convertible<t_vec, ublas::vector<typename t_vec::value_type>>::value, char>::type=0>
typename t_vec::value_type veclen(const t_vec& vec)
{
	return ublas::norm_2(vec);
}


/**
 * matrix element <x|M|y>
 */
template<typename t_mat = ublas::matrix<std::complex<double>>,
	typename t_vec = ublas::vector<std::complex<double>>>
typename t_vec::value_type
mat_elem(const t_vec& x, const t_mat& M, const t_vec& y)
{
	t_vec My = prod_mv<t_vec, t_mat>(M, y);
	t_vec x_conj = conjugate_vec<t_vec>(x);

	return inner<t_vec>(x_conj, My);
}


// ----------------------------------------------------------------------------


/**
 * remove an element from a vector
 */
template<class vector_type>
vector_type remove_elem(const vector_type& vec, std::size_t iIdx)
{
	vector_type vecret(vec.size()-1);

	for(std::size_t i=0, j=0; i<vec.size() && j<vecret.size();)
	{
		vecret[j] = vec[i];

		if(i!=iIdx) ++j;
		++i;
	}

	return vecret;
}


/**
 * create a submatrix removing row iRow and column iCol
 */
template<class matrix_type>
matrix_type submatrix(const matrix_type& mat, std::size_t iRow, std::size_t iCol)
{
	matrix_type matret(mat.size1()-1, mat.size2()-1);

	for(std::size_t i=0, i0=0; i<mat.size1() && i0<matret.size1();)
	{
		for(std::size_t j=0, j0=0; j<mat.size2() && j0<matret.size2();)
		{
			matret(i0,j0) = mat(i,j);

			if(j!=iCol) ++j0;
			++j;
		}

		if(i!=iRow) ++i0;
		++i;
	}

	return matret;
}


/**
 * create a submatrix
 */
template<class matrix_type>
matrix_type submatrix_wnd(const matrix_type& mat, std::size_t iSubRows, std::size_t iSubCols,
	std::size_t iBeginRow=0, std::size_t iBeginCol=0)
{
	matrix_type matret(iSubRows, iSubCols);

	for(std::size_t i=0; i<iSubRows; ++i)
		for(std::size_t j=0; j<iSubCols; ++j)
			matret(i, j) = mat(i+iBeginRow, j+iBeginCol);

	return matret;
}


/**
 * unite two matrices copying matrix 1 in the upper left block and matrix 2 in the lower right
 */
template<class matrix_type>
matrix_type block_matrix(const matrix_type& mat1, const matrix_type& mat2)
{
	matrix_type matret(mat1.size1()+mat2.size1(), mat1.size2()+mat2.size2());

	for(std::size_t i=0; i<matret.size1(); ++i)
	{
		for(std::size_t j=0; j<matret.size2(); ++j)
		{
			if(i>=0 && i<mat1.size1() && j>=0 && j<mat1.size2())
				matret(i, j) = mat1(i, j);
			else if(i>=mat1.size1() && j>=mat1.size2())
				matret(i, j) = mat2(i-mat1.size1(), j-mat1.size2());
			else
				matret(i, j) = 0;
		}
	}

	return matret;
}


template<class matrix_type>
matrix_type remove_column(const matrix_type& mat, std::size_t iCol)
{
	matrix_type matret(mat.size1(), mat.size2() - 1);
	for(std::size_t i = 0; i < mat.size1(); ++i)
	{
		for(std::size_t j = 0, j0 = 0; j<mat.size2() && j0 < matret.size2(); ++j)
		{
			matret(i, j0) = mat(i, j);
			if(j != iCol)
				++j0;
		}
	}
	return matret;
}


template<class matrix_type>
matrix_type remove_row(const matrix_type& mat, std::size_t iRow)
{
	matrix_type matret(mat.size1() - 1, mat.size2());
	for(std::size_t j = 0; j < mat.size2(); ++j)
	{
		for(std::size_t i = 0, i0 = 0; i<mat.size1() && i0 < matret.size1(); ++i)
		{
			matret(i0, j) = mat(i, j);
			if(i != iRow)
				++i0;
		}
	}
	return matret;
}


template<class matrix_type>
void submatrix_copy(matrix_type& mat, const matrix_type& sub,
	std::size_t iRowBegin, std::size_t iColBegin)
{
	for(std::size_t i=0; i<sub.size1(); ++i)
		for(std::size_t j=0; j<sub.size2(); ++j)
			mat(iRowBegin+i, iColBegin+j) = sub(i,j);
}


template<class vec_type>
void subvector_copy(vec_type& vec, const vec_type& sub, std::size_t iRowBegin)
{
	for(std::size_t i=0; i<sub.size(); ++i)
		vec[iRowBegin+i] = sub[i];
}


template<class matrix_type>
matrix_type remove_elems(const matrix_type& mat, std::size_t iIdx)
{
	return submatrix(mat, iIdx, iIdx);
}


/**
 * set matrix column
 */
template<class t_vec = ublas::vector<double>,
	class t_mat = ublas::matrix<typename t_vec::value_type>>
void set_column(t_mat& M, std::size_t iCol, const t_vec& vec)
{
	std::size_t s = std::min(vec.size(), M.size1());
	for(std::size_t i=0; i<s; ++i)
		M(i, iCol) = vec[i];
}


/**
 * set matrix row
 */
template<class t_vec = ublas::vector<double>,
	class t_mat = ublas::matrix<typename t_vec::value_type>>
void set_row(t_mat& M, std::size_t iRow, const t_vec& vec)
{
	std::size_t s = std::min(vec.size(), M.size2());
	for(std::size_t i=0; i<s; ++i)
		M(iRow, i) = vec[i];
}


/**
 * get matrix column -- general version
 */
template<class vector_type = ublas::vector<double>,
	class matrix_type = ublas::matrix<typename vector_type::value_type>,
	typename std::enable_if<!std::is_convertible<matrix_type, ublas::matrix<typename vector_type::value_type>>::value, char>::type=0>
vector_type get_column(const matrix_type& mat, std::size_t iCol)
{
	vector_type vecret(mat.size1());

	for(std::size_t i=0; i<mat.size1(); ++i)
		vecret[i] = mat(i, iCol);

	return vecret;
}


/**
 * get matrix column -- ublas wrapper
 */
template<class vector_type = ublas::vector<double>,
	class matrix_type = ublas::matrix<typename vector_type::value_type>,
	typename std::enable_if<std::is_convertible<matrix_type, ublas::matrix<typename vector_type::value_type>>::value, char>::type=0>
vector_type get_column(const matrix_type& mat, std::size_t iRow)
{
	return vector_type(ublas::column(mat, iRow));
}


/**
 * get matrix row -- general version
 */
template<class vector_type = ublas::vector<double>,
	class matrix_type = ublas::matrix<typename vector_type::value_type>,
	typename std::enable_if<!std::is_convertible<matrix_type, ublas::matrix<typename vector_type::value_type>>::value, char>::type=0>
vector_type get_row(const matrix_type& mat, std::size_t iRow)
{
	vector_type vecret(mat.size2());

	for(std::size_t i=0; i<mat.size2(); ++i)
		vecret[i] = mat(iRow, i);

	return vecret;
}


/**
 * get matrix row -- ublas wrapper
 */
template<class vector_type = ublas::vector<double>,
	class matrix_type = ublas::matrix<typename vector_type::value_type>,
	typename std::enable_if<std::is_convertible<matrix_type, ublas::matrix<typename vector_type::value_type>>::value, char>::type=0>
vector_type get_row(const matrix_type& mat, std::size_t iRow)
{
	return vector_type(ublas::row(mat, iRow));
}


template<class vector_type = ublas::vector<double>,
	class matrix_type = ublas::matrix<typename vector_type::value_type>,
	class cont_type = std::vector<vector_type>>
cont_type get_columns(const matrix_type& mat)
{
	cont_type vec;
	vec.reserve(mat.size2());

	for(std::size_t i=0; i<mat.size2(); ++i)
		vec.push_back(get_column(mat, i));

	return vec;
}


// ----------------------------------------------------------------------------


template<class t_mat = ublas::matrix<double>>
t_mat mirror_matrix(std::size_t iSize, std::size_t iComp)
{
	using T = typename t_mat::value_type;

	t_mat mat = unit_m<t_mat>(iSize);
	mat(iComp, iComp) = T(-1);

	return mat;
}


template<class matrix_type = ublas::matrix<double>>
matrix_type rotation_matrix_2d(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;

	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{c, -s},
		{s,  c}});
}


/**
 * generates points in an arc defined by vec1 and vec2 at an angle phi around vec1
 */
template<class t_vec = ublas::vector<double>>
t_vec arc(const t_vec& vec1, const t_vec& vec2, tl::underlying_value_type_t<t_vec> phi)
{
	//using t_real = tl::underlying_value_type_t<t_vec>;
	return std::cos(phi)*vec1 + std::sin(phi)*vec2;
}


/**
 * generates points in a spherical shell
 */
template<class t_vec = ublas::vector<double>>
t_vec sph_shell(const t_vec& vec,
	tl::underlying_value_type_t<t_vec> phi, tl::underlying_value_type_t<t_vec> theta)
{
	using t_real = tl::underlying_value_type_t<t_vec>;

	t_real rho, curphi, curtheta;
	std::tie(rho, curphi, curtheta) = cart_to_sph<t_real>(vec[0], vec[1], vec[2]);

	t_real x,y,z;
	std::tie(x,y,z) = sph_to_cart<t_real>(rho, curphi+phi, curtheta+theta);
	t_vec vecRet = make_vec<t_vec>({x,y,z});
	return vecRet;
}


template<class matrix_type = ublas::matrix<double>>
matrix_type rotation_matrix_3d_x(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{1, 0,  0},
		{0, c, -s},
		{0, s,  c}});
}


template<class matrix_type = ublas::matrix<double>>
matrix_type rotation_matrix_3d_y(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{c,  0, s},
		{0,  1, 0},
		{-s, 0, c}});
}


template<class matrix_type = ublas::matrix<double>>
matrix_type rotation_matrix_3d_z(typename matrix_type::value_type angle)
{
	typedef typename matrix_type::value_type T;

	T s, c;
	if(angle==0.)
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return make_mat<matrix_type>
	({	{c, -s, 0},
		{s,  c, 0},
		{0,  0, 1}});
}


/**
 * cross-product in matrix form
 * @see https://en.wikipedia.org/wiki/Skew-symmetric_matrix
 */
template<class matrix_type = ublas::matrix<double>,
	class vector_type = ublas::vector<typename matrix_type::value_type>>
matrix_type skew(const vector_type& vec)
{
	if(vec.size() == 3)
	{
		return make_mat<matrix_type>
		({	{       0, -vec[2],  vec[1] },
			{  vec[2],       0, -vec[0] },
			{ -vec[1],  vec[0],       0 }});
	}
	else
		throw Err("Skew only defined for three dimensions.");
}


/**
 * diagonal matrix
 */
template<class matrix_type = ublas::matrix<double>,
	class cont_type = std::initializer_list<typename matrix_type::value_type>>
matrix_type diag_matrix(const cont_type& lst)
{
	matrix_type mat = unit_m<matrix_type>(lst.size());

	std::size_t i = 0;
	for(typename cont_type::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter, ++i)
		mat(i,i) = *iter;

	return mat;
}


/**
 * vector of diagonal matrix elements
 */
template<class t_vec = ublas::vector<double>,
	class t_mat = ublas::matrix<double>>
t_vec diag_vec(const t_mat& mat)
{
	std::size_t N = std::min(mat.size1(), mat.size2());

	t_vec vec(N);
	for(std::size_t i=0; i<N; ++i)
		vec[i] = mat(i,i);

	return vec;
}


template<class matrix_type = ublas::matrix<double>,
	class cont_type = std::initializer_list<typename matrix_type::value_type>>
matrix_type scale_matrix(const cont_type& lst)
{
	return diag_matrix<matrix_type, cont_type>(lst);
}


/**
 * translation matrix in homogeneous coords
 */
template<class t_mat = ublas::matrix<double>,
	class t_cont = std::initializer_list<typename t_mat::value_type>>
t_mat translation_matrix(const t_cont& lst)
{
	t_mat mat = unit_m<t_mat>(lst.size()+1);

	const std::size_t iJ = mat.size2();
	std::size_t i = 0;
	for(typename t_cont::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter, ++i)
		mat(i, iJ-1) = *iter;

	return mat;
}


/**
 * is mat a translation matrix in homogeneous coords?
 */
template<class t_mat = ublas::matrix<double>>
bool is_translation_matrix(const t_mat& mat)
{
	using T = typename t_mat::value_type;

	const std::size_t iI = mat.size1();
	const std::size_t iJ = mat.size2();

	for(std::size_t i=0; i<iI-1; ++i)
	{
		if(!float_equal<T>(mat(i, iJ-1), T(0)))
			return true;
	}
	return false;
}


template<class t_mat = ublas::matrix<double>>
bool is_identity_matrix(const t_mat& mat)
{
	using T = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return false;

	const std::size_t iN = mat.size1();

	for(std::size_t i=0; i<iN; ++i)
	{
		for(std::size_t j=0; j<iN; ++j)
		{
			if(i != j)	// off-diagonal elements
			{
				if(!float_equal<T>(mat(i, j), T(0)))
					return false;
			}
			else	// diagonal elements
			{
				if(!float_equal<T>(mat(i, j), T(1)))
					return false;
			}
		}
	}

	return true;
}


/**
 * do the absolute elements form an identity matrix
 * also return the indices of the negative values on the diagonal
 */
template<class t_mat = ublas::matrix<double>>
std::pair<bool, std::vector<std::size_t>> is_abs_identity_matrix(const t_mat& mat)
{
	std::vector<std::size_t> vecMinusses;

	using T = typename t_mat::value_type;
	if(mat.size1() != mat.size2())
		return std::make_pair(false, vecMinusses);

	const std::size_t iN = mat.size1();
	for(std::size_t i=0; i<iN; ++i)
	{
		for(std::size_t j=0; j<iN; ++j)
		{
			if(i != j)	// off-diagonal elements
			{
				if(!float_equal<T>(mat(i, j), T(0)))
					return std::make_pair(false, vecMinusses);
			}
			else	// diagonal elements
			{
				if(!float_equal<T>(std::abs(mat(i, j)), T(1)))
					return std::make_pair(false, vecMinusses);
				if(mat(i,j) < T(0))
					vecMinusses.push_back(i);
			}
		}
	}

	return std::make_pair(true, vecMinusses);
}


template<class t_mat = ublas::matrix<double>>
bool is_inverting_matrix(const t_mat& mat)
{ return is_identity_matrix(-mat); }


/**
 * does the homogeneous matrix mat have a translation component?
 */
template<class t_mat = ublas::matrix<double>>
bool has_translation_components(const t_mat& mat)
{
	using T = typename t_mat::value_type;
	const std::size_t iN = mat.size1();
	if(iN != mat.size2())
		return false;

	// translation?
	for(std::size_t i=0; i<iN-1; ++i)
	{
		if(!float_equal<T>(mat(i, iN-1), T(0)))
			return true;
	}

	return false;
}


/**
 * is mat a centering matrix in homogeneous coords?
 */
template<class t_mat = ublas::matrix<double>>
bool is_centering_matrix(const t_mat& mat)
{
	//if(is_identity_matrix(mat)) return 1;

	const std::size_t iN = mat.size1();
	if(iN != mat.size2())
		return false;

	// left-upper 3x3 unit matrix?
	if(!is_identity_matrix(submatrix(mat, iN-1, iN-1)))
		return false;

	// translation?
	if(has_translation_components<t_mat>(mat))
		return true;
	return false;
}


/**
 * Rodrigues' formula
 * @see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
 * @see (Merziger 1993), p. 208
 * @see (Arens 2015), p. 718 and p. 816
 */
template<class mat_type = ublas::matrix<double>,
	class vec_type = ublas::vector<typename mat_type::value_type>,
	typename T = typename mat_type::value_type>
mat_type rotation_matrix(const vec_type& _vec, T angle)
{
	const vec_type vec = _vec/veclen(_vec);

	T s, c;
	if(angle == T(0))
	{
		s = T(0);
		c = T(1);
	}
	else
	{
		s = std::sin(angle);
		c = std::cos(angle);
	}

	return (T(1) - c) * outer(vec, vec) +
		c * unit_m(vec.size()) +
		s * skew(vec);
}


template<class matrix_type = ublas::matrix<double>>
typename matrix_type::value_type trace(const matrix_type& mat)
{
	typedef typename matrix_type::value_type T;

	if(mat.size1() != mat.size2())
		return T(0);

	T tr = T(0.);
	for(std::size_t i=0; i<mat.size1(); ++i)
		tr += mat(i,i);
	return tr;
}


/**
 * parallel or perspectivic projection matrix
 * @see https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
 * @see https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
 */
template<class t_mat = ublas::matrix<double, ublas::row_major, ublas::bounded_array<double,4*4>>,
	class T = typename t_mat::value_type>
t_mat proj_matrix(T l, T r, T b, T t, T n, T f, bool bParallel)
{
	// scale to [0:2]
	t_mat matScale = scale_matrix<t_mat>({ T(2)/(r-l), T(2)/(t-b), T(1), T(1) });
	// translate to [-1:1]
	t_mat matTrans = translation_matrix<t_mat>({ -T(0.5)*(r+l), -T(0.5)*(t+b), T(0) });
	matScale = prod_mm(matScale, matTrans);

	// project
	t_mat matProj = unit_m<t_mat>(4);

	if(bParallel)	// parallel
	{
		matProj(2,2) = T(2)*f*n / (n-f);
		matProj(2,3) = (n+f) / (n-f);
	}
	else			// perspectivic
	{
		matProj(2,2) = (n+f) / (n-f);
		matProj(2,3) = T(2)*f*n / (n-f);
		matProj(3,2) = T(-1);
		matProj(3,3) = T(0);
	}

	return prod_mm(matScale, matProj);
}


/**
 * parallel projection matrix
 * @see https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
 */
template<class t_mat = ublas::matrix<double, ublas::row_major, ublas::bounded_array<double,4*4>>,
	class T = typename t_mat::value_type>
t_mat ortho_matrix(T l, T r, T b, T t, T n, T f)
{
	return proj_matrix<t_mat, T>(l,r, b,t, n,f, 1);
}


/**
 * perspectivic projection matrix
 * @see https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
 * also see: similar gnomonic projection of spherical coordinates onto a plane
 */
template<class t_mat = ublas::matrix<double, ublas::row_major, ublas::bounded_array<double,4*4>>,
	class T = typename t_mat::value_type>
t_mat perspective_matrix(T yfov, T asp, T n, T f)
{
	const T y = std::tan(T(0.5)*yfov);
	const T x = y*asp;

	return proj_matrix<t_mat, T>(-x,x, -y,y, n,f, 0);
}


// -----------------------------------------------------------------------------
template<typename T, class FKT, const int iDim=get_type_dim<T>::value>
struct is_nan_or_inf_impl {};


template<typename real_type, class FKT>
struct is_nan_or_inf_impl<real_type, FKT, 0>	// scalar impl.
{
	const FKT& m_fkt;
	is_nan_or_inf_impl(const FKT& fkt) : m_fkt(fkt) {}
	bool operator()(real_type d) const { return m_fkt(d); }
};


template<typename vec_type, class FKT>
struct is_nan_or_inf_impl<vec_type, FKT, 1>		// vector impl.
{
	const FKT& m_fkt;
	is_nan_or_inf_impl(const FKT& fkt) : m_fkt(fkt) {}

	bool operator()(const vec_type& vec) const
	{
		for(std::size_t i=0; i<vec.size(); ++i)
			if(m_fkt(vec[i]))
				return true;
		return false;
	}
};


template<typename mat_type, class FKT>
struct is_nan_or_inf_impl<mat_type, FKT, 2>		// matrix impl.
{
	const FKT& m_fkt;
	is_nan_or_inf_impl(const FKT& fkt) : m_fkt(fkt) {}

	bool operator()(const mat_type& mat) const
	{
		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
				if(m_fkt(mat(i,j)))
					return true;
		return false;
	}
};


template<class T = ublas::matrix<double>>
bool isnan(const T& mat)
{
	typedef underlying_value_type_t<T> real_type;

	using fkt = std::function<bool(real_type)>;
	fkt stdisnan = (bool(*)(real_type))std::isnan;
	is_nan_or_inf_impl<T, fkt> _isnan(stdisnan);
	return _isnan(mat);
}


template<class T = ublas::matrix<double>>
bool isinf(const T& mat)
{
	typedef underlying_value_type_t<T> real_type;

	using fkt = std::function<bool(real_type)>;
	fkt stdisinf = (bool(*)(real_type))std::isinf;
	is_nan_or_inf_impl<T, fkt> _isinf(stdisinf);
	return _isinf(mat);
}


template<class T = ublas::matrix<double>>
bool is_nan_or_inf(const T& mat)
{
	typedef underlying_value_type_t<T> real_type;

	using fkt = std::function<bool(real_type)>;
	fkt stdisnaninf = [](real_type d)->bool { return std::isnan(d) || std::isinf(d); };
	is_nan_or_inf_impl<T, fkt> _isnaninf(stdisnaninf);
	return _isnaninf(mat);
}
// -----------------------------------------------------------------------------


/**
 * calculates the matrix inverse
 *
 * @desc code for inverse based on https://github.com/boostorg/ublas/blob/develop/test/test_lu.cpp
 * @desc Boost's test_lu.cpp is (c) 2008 by G. Winkler
 */
template<class mat_type = ublas::matrix<double>>
bool inverse(const mat_type& mat, mat_type& inv)
{
	const typename mat_type::size_type N = mat.size1();

	if(N != mat.size2())
		return false;

	try
	{
		mat_type lu = mat;
		ublas::permutation_matrix<typename mat_type::size_type> perm(N);

		if(ublas::lu_factorize(lu, perm) != 0)
			return false;

		inv = unit_m<mat_type>(N);
		ublas::lu_substitute(lu, perm, inv);
	}
	catch(const std::exception& ex)
	{
		log_err("Matrix inversion failed with exception: ", ex.what(), ".", "\n",
			"Matrix to be inverted was: ", mat, ".");
		return false;
	}

	return true;
}


/**
 * calculates the matrix inverse of a diagonal matrix
 */
template<class mat_type = ublas::matrix<double>>
bool inverse_diag(const mat_type& mat, mat_type& inv)
{
	using T = typename mat_type::value_type;
	const typename mat_type::size_type N = mat.size1();

	if(N != mat.size2())
		return false;

	try
	{
		inv = zero_m<mat_type>(N, N);
		for(std::size_t i=0; i<N; ++i)
			inv(i, i) = T(1) / mat(i, i);
	}
	catch(const std::exception& ex)
	{
		log_err("Diagonal matrix inversion failed with exception: ", ex.what(), ".", "\n",
			"Matrix to be inverted was: ", mat, ".");
		return false;
	}

	return true;
}


/**
 * R = T^(-1) M T
 * bCongr==1: do a congruence trafo
 * bCongr==0: do a similarity trafo
 * @see e.g.: (Merziger 1993), p. 202
 */
template<class mat_type = ublas::matrix<double>>
mat_type transform(const mat_type& mat, const mat_type& matTrafo, bool bCongr = false)
{
	mat_type matTrafoInv;
	if(bCongr)
		matTrafoInv = transpose(matTrafo);
	else
		inverse(matTrafo, matTrafoInv);

	mat_type MT = prod_mm(mat, matTrafo);
	mat_type TinvMT = prod_mm(matTrafoInv, MT);

	return TinvMT;
}


/**
 * R = T M T^(-1)
 * bCongr==1: do a congruence trafo
 * bCongr==0: do a similarity trafo
 * @see e.g.: (Merziger 1993), p. 202
 */
template<class mat_type = ublas::matrix<double>>
mat_type transform_inv(const mat_type& mat, const mat_type& matTrafo, bool bCongr = false)
{
	mat_type matTrafoInv;
	if(bCongr)
		matTrafoInv = transpose(matTrafo);
	else
		inverse(matTrafo, matTrafoInv);

	mat_type MT = prod_mm(mat, matTrafoInv);
	mat_type TinvMT = prod_mm(matTrafo, MT);

	return TinvMT;
}


// -> linalg2.h
template<typename T=double>
bool qr(const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R);


template<typename T>
bool solve_linear_approx(const ublas::matrix<T>& M, const ublas::vector<T>& v,
	ublas::vector<T>& x);


/**
 * solve Mx = v for x
 */
template<typename T = double>
bool solve_linear(const ublas::matrix<T>& M,
	const ublas::vector<T>& v, ublas::vector<T>& x)
{
	if(M.size1() == M.size2())		// determined, TODO: check rank
	{
		try
		{
			const std::size_t N = M.size1();

			ublas::matrix<T> lu = M;
			ublas::permutation_matrix<typename ublas::matrix<T>::size_type> perm(N);

			typename ublas::matrix<T>::size_type sing = ublas::lu_factorize(lu, perm);
			if(sing != 0)
				return false;

			x = v;
			ublas::lu_substitute(lu, perm, x);
		}
		catch(const std::exception& ex)
		{
			log_err("Linear equation solver failed with exception: ", ex.what(), ".");
			return false;
		}
	}
	else if(M.size1() < M.size2())	// underdetermined
	{
		ublas::matrix<T> Q, R;
		if(!qr(M, Q, R))
			return false;
		typedef typename ublas::vector<T>::size_type t_int;

		// M x = v
		// QR x = v
		// R x = Q^T v

		ublas::vector<T> vnew = prod_mv(transpose(Q), v);

		x = zero_v<ublas::vector<T>>(M.size2());
		ublas::vector<T> xnew(R.size1());
		bool bOk = 0;


		// find non-singular right-upper submatrix
		std::vector<t_int> vecDelCols;
		std::size_t iNumToDel = R.size2()-R.size1();
		if(iNumToDel != 1)
		{
			log_err(__func__, " not yet implemented.");
			return false;
		}

		bool bFoundNonSingular = 0;
		ublas::matrix<T> Rsub;
		for(std::ptrdiff_t iCol=std::ptrdiff_t(R.size2()-1); iCol>=0; --iCol)
		{
			Rsub = remove_column(R, (std::size_t)iCol);

			T det = determinant<ublas::matrix<T>>(Rsub);
			if(!float_equal<T>(det, 0.))
			{
				bFoundNonSingular = 1;
				vecDelCols.push_back(iCol);
				break;
			}
		}

		if(!bFoundNonSingular)
		{
			log_err("No non-singluar submatrix found in linear equation solver.");
			return false;
		}

		bOk = solve_linear(Rsub, vnew, xnew);

		for(t_int i=0, i0=0; i<xnew.size() && i0<x.size(); ++i, ++i0)
		{
			while(std::find(vecDelCols.begin(), vecDelCols.end(), i0) != vecDelCols.end())
				++i0;
			x[i0] = xnew[i];
		}

		return bOk;
	}
	else if(M.size1() > M.size2())	// overdetermined
		return solve_linear_approx<T>(M,v,x);
	else
		return false;

	return true;
}


/**
 * solve M^T M x = M^T v for x
 * @see e.g. (Arens 2015), p. 793
 */
template<typename T = double>
bool solve_linear_approx(const ublas::matrix<T>& M,
	const ublas::vector<T>& v, ublas::vector<T>& x)
{
	if(M.size1() <= M.size2())
	{
		return false;
	}

	ublas::matrix<T> Q, R;
	if(!qr(M, Q, R))
		return false;

	// M^T M x = M^T v
	// R^T Q^T Q R x = R^T Q^T v
	// R^T R x = R^T Q^T v

	const ublas::matrix<T> RT = transpose(R);
	const ublas::matrix<T> QT = transpose(Q);
	const ublas::matrix<T> RTR = prod_mm(RT, R);
	const ublas::matrix<T> RTQT = prod_mm(RT, QT);

	const ublas::vector<T> vnew = prod_mv(RTQT, v);
	return solve_linear<T>(RTR, vnew, x);
}


template<class matrix_type = ublas::matrix<double>>
bool is_diag_matrix(const matrix_type& mat)
{
	typedef typename matrix_type::value_type T;

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			if(i==j) continue;

			if(!float_equal(mat(i,j), T(0.)))
				return false;
		}

	return true;
}


// ----------------------------------------------------------------------------


template<class matrix_type = ublas::matrix<double>,
	class vec_type = ublas::vector<typename matrix_type::value_type>,
	class container_type = std::initializer_list<vec_type>, const bool bRowMat>
inline matrix_type row_col_matrix(const container_type& vecs)
{
	if(vecs.size() == 0)
		return matrix_type(0,0);

	const std::size_t N = vecs.size();
	const std::size_t M = vecs.begin()->size();

	matrix_type mat(bRowMat?N:M, bRowMat?M:N);
	std::size_t j=0;
	for(typename container_type::const_iterator iter=vecs.begin(); iter!=vecs.end(); ++iter)
	{
		const vec_type& vec = *iter;

		for(std::size_t i=0; i<vec.size(); ++i)
		{
			if(bRowMat)
				mat(j,i) = vec[i];
			else
				mat(i,j) = vec[i];
		}

		++j;
	}

	return mat;
}


/**
 * vectors form rows of matrix
 */
template<class matrix_type = ublas::matrix<double>,
	class vec_type = ublas::vector<typename matrix_type::value_type>,
	class container_type = std::initializer_list<vec_type>>
matrix_type row_matrix(const container_type& vecs)
{
	return row_col_matrix<matrix_type, vec_type, container_type, true>(vecs);
}


/**
 * vectors form columns of matrix
 */
template<class matrix_type = ublas::matrix<double>,
	class vec_type = ublas::vector<typename matrix_type::value_type>,
	class container_type = std::initializer_list<vec_type>>
matrix_type column_matrix(const container_type& vecs)
{
	return row_col_matrix<matrix_type, vec_type, container_type, false>(vecs);
}


// ----------------------------------------------------------------------------


/**
 * determinant
 * @see e.g.: (Merziger 1993), p. 185
 * @see e.g.: https://en.wikipedia.org/wiki/Determinant
 */
template<class t_mat/*=ublas::matrix<double>*/>
typename t_mat::value_type determinant(const t_mat& mat)
{
	typedef typename t_mat::value_type T;
	typedef typename t_mat::size_type t_size;

	if(mat.size1() != mat.size2())
		return T(0);

	if(mat.size1()==0)
		return T(0);
	else if(mat.size1()==1)
		return mat(0,0);
	else if(mat.size1()==2)
		return mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
	else if(mat.size1()==3)
	{
		ublas::vector<T> vec0 = get_column(mat, 0);
		ublas::vector<T> vec1 = get_column(mat, 1);
		ublas::vector<T> vec2 = get_column(mat, 2);

		ublas::vector<T> vecCross = cross_3<ublas::vector<T> >(vec1, vec2);
		return inner(vec0, vecCross);
	}
	else if(mat.size1()>3 && mat.size1()<6)		// recursive expansion, complexity: O(n!)
	{
		t_size i = 0;
		t_size iZeros = 0;

		// count zeros
		for(t_size _i=0; _i<mat.size1(); ++_i)
		{
			ublas::vector<T> vecRow = get_row<ublas::vector<T>, t_mat>(mat, _i);
			t_size iNewZeros = std::count_if(vecRow.begin(), vecRow.end(),
				[](T d) -> bool { return float_equal<T>(d, 0.); });
			if(iNewZeros > iZeros)
			{
				i = _i;
				iZeros = iNewZeros;
			}
		}


		T val = T(0);

		for(t_size j=0; j<mat.size2(); ++j)
		{
			if(float_equal<T>(mat(i,j), 0.))
				continue;

			T dSign = 1.;
			if(is_odd<std::size_t>(i+j))
				dSign = -1.;

			t_mat matSub = submatrix(mat, i, j);
			val += dSign * mat(i,j) * determinant<t_mat>(matSub);
		}

		return val;
	}
	else if(mat.size1()>=6)				// LU decomposition, complexity: O(n^3)
	{
		t_mat lu = mat;
		t_size N = mat.size1();
		ublas::permutation_matrix<typename t_mat::size_type> perm(N);

		ublas::lu_factorize(lu, perm);

		t_mat L = ublas::triangular_adaptor<t_mat, ublas::unit_lower>(lu);
		t_mat U = ublas::triangular_adaptor<t_mat, ublas::upper>(lu);

		T dDet = T(1.);
		for(t_size i=0; i<mat.size1(); ++i)
			dDet *= L(i,i)*U(i,i);

		std::size_t iNumSwaps=0;
		for(t_size iSwap=0; iSwap<perm.size(); ++iSwap)
			if(iSwap != perm(iSwap))
				++iNumSwaps;

		if(is_odd<std::size_t>(iNumSwaps))
			dDet *= T(-1.);

		return dDet;
	}

	return T(0);
}


/**
 * minor determinant
 * @see e.g.: https://en.wikipedia.org/wiki/Minor_(linear_algebra)
 */
template<class t_mat = ublas::matrix<double>>
typename t_mat::value_type minor_det(const t_mat& mat, std::size_t iRow, std::size_t iCol)
{
	t_mat M = submatrix(mat, iRow, iCol);
	return determinant<t_mat>(M);
}


/**
 * cofactor
 * @see e.g.: https://en.wikipedia.org/wiki/Minor_(linear_algebra)
 */
template<class t_mat = ublas::matrix<double>>
typename t_mat::value_type cofactor(const t_mat& mat, std::size_t iRow, std::size_t iCol)
{
	using T = typename t_mat::value_type;

	T m = minor_det(mat, iRow, iCol);
	T s = std::pow(T(-1), T(iRow+1 + iCol+1));

	return m*s;
}


/**
 * adjugate matrix
 * @see e.g.: https://en.wikipedia.org/wiki/Adjugate_matrix
 */
template<class t_mat = ublas::matrix<double>>
t_mat adjugate(const t_mat& mat, bool bTranspose=1)
{
	using T = typename t_mat::value_type;

	t_mat matRet(mat.size1(), mat.size2());

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=0; j<mat.size2(); ++j)
		{
			T c = cofactor<t_mat>(mat, i, j);
			matRet(i,j) = c;
		}

	if(bTranspose)
		matRet = transpose(matRet);
	return matRet;
}


template<class matrix_type = ublas::matrix<double>>
typename matrix_type::value_type get_volume(const matrix_type& mat)
{
	//typedef typename matrix_type::value_type T;
	return std::abs(determinant<matrix_type>(mat));
}


template<class matrix_type = ublas::matrix<double>>
typename matrix_type::value_type get_ellipsoid_volume(const matrix_type& mat)
{
	typedef typename matrix_type::value_type T;
	T tDet = std::abs(determinant<matrix_type>(mat));

	return T(4./3.) * get_pi<T>() * std::sqrt(T(1)/tDet);
}


// ----------------------------------------------------------------------------


/**
 * calculate fractional coordinate basis vectors from angles
 * @see http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_75.html
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 * for the reciprocal lattice this is equal to the B matrix from Acta Cryst. (1967), 22, 457
 */
template<class t_vec>
bool fractional_basis_from_angles(typename t_vec::value_type a,
	typename t_vec::value_type b,
	typename t_vec::value_type c,
	typename t_vec::value_type alpha,
	typename t_vec::value_type beta,
	typename t_vec::value_type gamma,
	t_vec& veca, t_vec& vecb, t_vec& vecc)
{
	typedef typename t_vec::value_type T;

	const T dSG = std::sin(gamma), dCG = std::cos(gamma);
	const T dCA = std::cos(alpha), dCB = std::cos(beta);
	const T dCA2 = dCA*dCA, dCB2 = dCB*dCB, dCG2 = dCG*dCG;

	const T dVol = a*b*c *
		std::sqrt(T(1) - dCA2 - dCB2 - dCG2 + T(2)*dCA*dCB*dCG);
	if(std::isinf(dVol) || std::isnan(dVol))
		return false;

	if(veca.size() != 3) veca.resize(3);
	if(vecb.size() != 3) vecb.resize(3);
	if(vecc.size() != 3) vecc.resize(3);

	veca[0] = a;
	veca[1] = T(0);
	veca[2] = T(0);

	vecb[0] = b*dCG;
	vecb[1] = b*dSG;
	vecb[2] = T(0);

	vecc[0] = c*dCB;
	vecc[1] = c*(dCA - dCB*dCG) / dSG;
	vecc[2] = dVol / (a*b*dSG);

	return true;
}

// ----------------------------------------------------------------------------


/**
 * signed angle wrt basis
 */
template<typename vec_type>
typename vec_type::value_type vec_angle(const vec_type& vec)
{
	if(vec.size() == 2)
		return std::atan2(vec[1], vec[0]);

	throw Err("vec_angle not yet implemented for size != 2.");
}


// -----------------------------------------------------------------------------
/**
 * set values lower than epsilon to zero
 */
template<typename T> void set_eps_0(T& d, underlying_value_type_t<T> eps=-1.);
template<typename T, LinalgType ty=get_linalg_type<T>::value> struct set_eps_0_impl {};


/**
 * set values lower than epsilon to zero
 * scalar version
 */
template<typename real_type>
struct set_eps_0_impl<real_type, LinalgType::REAL>
{
	real_type eps = get_epsilon<real_type>();

	void operator()(real_type& d) const
	{
		if(std::abs(d) < eps)
			d = real_type(0);
	}
};


/**
 * set values lower than epsilon to zero
 * vector version
 */
template<typename vec_type>
struct set_eps_0_impl<vec_type, LinalgType::VECTOR>
{
	using real_type = typename vec_type::value_type;
	real_type eps = get_epsilon<real_type>();

	void operator()(vec_type& vec) const
	{
		for(real_type& d : vec)
			set_eps_0<real_type>(d, eps);
	}
};


/**
 * set values lower than epsilon to zero
 * matrix version
 */
template<typename mat_type>
struct set_eps_0_impl<mat_type, LinalgType::MATRIX>
{
	using real_type = typename mat_type::value_type;
	real_type eps = get_epsilon<real_type>();

	void operator()(mat_type& mat) const
	{
		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
				set_eps_0<real_type>(mat(i,j), eps);
	}
};


template<typename T>
void set_eps_0(T& d, underlying_value_type_t<T> eps)
{
	set_eps_0_impl<T, get_linalg_type<T>::value> op;
	if(eps >= underlying_value_type_t<T>(0))
		op.eps = eps;
	op(d);
}
// -----------------------------------------------------------------------------


template<typename t_vec, typename T = typename t_vec::value_type>
bool vec_is_collinear(const t_vec& _vec1, const t_vec& _vec2, T eps = get_epsilon<T>())
{
	const t_vec vec1 = _vec1 / veclen(_vec1);
	const t_vec vec2 = _vec2 / veclen(_vec2);

	T tdot = std::abs(inner(vec1, vec2));
	return float_equal<T>(tdot, 1, eps);
}


/**
 * signed angle between two vectors
 */
template<typename vec_type>
typename vec_type::value_type vec_angle(const vec_type& vec0,
	const vec_type& vec1, const vec_type* pvec_norm=nullptr)
{
	typedef typename vec_type::value_type real_type;

	if(vec0.size() != vec1.size())
		throw Err("In vec_angle: Vector sizes do not match.");

	if(vec0.size() == 2)
	{
		return vec_angle<vec_type>(vec0) - vec_angle<vec_type>(vec1);
	}
	if(vec0.size() == 3)
	{
		real_type dC = inner(vec0, vec1);
		vec_type veccross = cross_3<vec_type>(vec0, vec1);
		real_type dS = veclen(veccross);

		real_type dAngle = std::atan2(dS, dC);

		// get signed angle
		if(pvec_norm)
		{
			if(inner(veccross, *pvec_norm) < real_type(0))
				dAngle = -dAngle;
		}

		return dAngle;
	}

	throw Err("vec_angle only implemented for size == 2 and size == 3.");
}


template<class T, LinalgType ty=get_linalg_type<T>::value>
struct vec_angle_unsigned_impl {};


/**
 * unsigned angle between two vectors
 */
template<class T>
struct vec_angle_unsigned_impl<T, LinalgType::VECTOR>
{
	typename T::value_type operator()(const T& q1, const T& q2) const
	{
		typedef typename T::value_type REAL;

		if(q1.size() != q2.size())
			return REAL();

		REAL dot = REAL();
		REAL len1 = REAL();
		REAL len2 = REAL();
		for(std::size_t i=0; i<q1.size(); ++i)
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
};


template<class T>
typename T::value_type vec_angle_unsigned(const T& q1, const T& q2)
{
	return vec_angle_unsigned_impl<T>()(q1, q2);
}

// -----------------------------------------------------------------------------


/**
 * @see K. Shoemake, "Animating rotation with quaternion curves", http://dx.doi.org/10.1145/325334.325242
 * @see (Desktop Bronstein 2008), formula 4.207
 * @see (Bronstein 2008), p. 306, formula 4.155
 */
template<class T>
T slerp(const T& q1, const T& q2, typename T::value_type t)
{
	typedef typename T::value_type REAL;

	REAL angle = vec_angle_unsigned<T>(q1, q2);

	T q = std::sin((REAL(1)-t)*angle)/std::sin(angle) * q1 +
		std::sin(t*angle)/std::sin(angle) * q2;

	return q;
}



// --------------------------------------------------------------------------------

}


#include <boost/version.hpp>

#if BOOST_VERSION >= 106600
	#include <boost/integer/common_factor_rt.hpp>
	namespace integer = boost::integer;
#else
	#include <boost/math/common_factor_rt.hpp>
	namespace integer = boost::math;
#endif


namespace tl{

template<class t_vec=ublas::vector<int>>
t_vec get_gcd_vec(const t_vec& vec)
{
	if(vec.size() <= 1)
		return vec;

	typedef typename t_vec::value_type t_int;

	t_int igcd_total = 1;
	for(std::size_t i=0; i<vec.size()-1; ++i)
	{
		t_int i0 = vec[i];
		t_int i1 = vec[i+1];

		t_int igcd = integer::gcd<t_int>(i0, i1);

		if(i==0)
			igcd_total = igcd;
		else
			igcd_total = integer::gcd<t_int>(igcd, igcd_total);
	}

	if(igcd_total == 0)
		return vec;

	return vec/igcd_total;
}


// --------------------------------------------------------------------------------

/**
 * Householder reflection matrix
 * @see (Scarpino 2011), p. 268
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
	t_mat reflection_matrix(const t_vec& vecNorm)
{
	// projection of "pt" onto (normalised) vector "norm":
	// proj = (norm^t * pt) * norm
	// proj = (norm * norm^t) * pt
	// "Lotfusspunkt" = pt - proj
	// mirror_point = pt - 2*proj
	t_mat mat = -T(2) * outer(vecNorm, vecNorm);
	mat /= inner(vecNorm, vecNorm);

	for(std::size_t i=0; i<vecNorm.size(); ++i)
		mat(i,i) += T(1);

	return mat;
}


/**
 * Householder reflection
 * @see (Scarpino 2011), p. 268
 */
template<class t_vec = ublas::vector<double>,
	class t_mat = ublas::matrix<typename t_vec::value_type>,
	typename T = typename t_mat::value_type>
	t_vec reflection(const t_vec& vec, const t_vec& vecNorm)
{
	t_mat mat = reflection_matrix<t_mat, t_vec, T>(vecNorm);
	return prod_mv(mat, vec);
}


/**
 * add a nxn unit matrix to the upper left of a matrix
 */
template<class t_mat = ublas::matrix<double>,
	typename T = typename t_mat::value_type>
	t_mat insert_unity(const t_mat& M, std::size_t n)
{
	if(M.size1()!=M.size2())
		throw Err("Non-square matrix not yet supported.");

	std::size_t m = M.size1();
	t_mat M2 = t_mat(m+n, m+n);

	for(std::size_t iR=0; iR<m+n; ++iR)
	{
		for(std::size_t jR=0; jR<m+n; ++jR)
		{
			if(iR<n || jR<n)
				M2(iR, jR) = (iR==jR ? 1 : 0);
			else
				M2(iR, jR) = M(iR-n, jR-n);
		}
	}

	return M2;
}


/**
 * QR decomposition via householder reflections
 * @see (Scarpino 2011), pp. 269-272
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool qr_decomp(const t_mat& M, t_mat& Q, t_mat& R)
{
	std::size_t m = M.size1();
	std::size_t n = M.size2();

	t_mat A = M;
	std::vector<t_mat> vecRefls;

	for(std::size_t i=0; i<std::min(m-1,n); ++i)
	{
		t_vec vec0 = get_column(A, 0);

		// vector of form [123.4 0 0 0] ?
		t_vec vec0_rest = ublas::subrange(vec0, 1, vec0.size());
		if(vec_equal<t_vec>(vec0_rest, zero_v<t_vec>(vec0_rest.size())))
		{
			t_mat matReflM = unit_m(m);
			vecRefls.push_back(matReflM);
			continue;
		}

		t_vec vecE0 = zero_v<t_vec>(vec0.size());
		vecE0[0] = veclen(vec0);

		t_vec vecReflNorm = vec0 - vecE0;
		t_mat matRefl = reflection_matrix(vecReflNorm);

		A = prod_mm(matRefl, A);
		A = submatrix(A,0,0);

		t_mat matReflM = insert_unity(matRefl, m-matRefl.size1());
		vecRefls.push_back(matReflM);
	}

	if(vecRefls.size() == 0)
		return false;

	Q = unit_m(m);
	for(const t_mat& matRefl : vecRefls)
	{
		t_mat matReflT = transpose(matRefl);
		Q = prod_mm(Q, matReflT);
	}

	t_mat QT = transpose(Q);
	R = prod_mm(QT, M);

	return true;
}


template<typename t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type>
std::vector<t_vec> gram_schmidt(const std::vector<t_vec>& vecs, bool bNorm = true);


/**
 * QR decomposition via gram-schmidt orthogonalisation
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool qr_decomp_gs(const t_mat& M, t_mat& Q, t_mat& R)
{
	Q = column_matrix(gram_schmidt(get_columns(M), 1));

	// M = QR  =>  Q^T M = R
	R = prod_mm(transpose(Q), M);
	return 1;
}


template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
t_mat norm_col_vecs(const t_mat& M)
{
	t_mat N(M.size1(), M.size2());

	for(std::size_t i=0; i<M.size2(); ++i)
	{
		t_vec vec0 = get_column(M, i);
		vec0 /= veclen(vec0);

		set_column(N, i, vec0);
	}

	return N;
}


template<class t_mat=ublas::matrix<double>, class t_real=underlying_value_type_t<t_mat>>
bool is_symmetric(const t_mat& mat, t_real eps = get_epsilon<t_real>())
{
	if(mat.size1() != mat.size2())
		return false;

	for(std::size_t i=0; i<mat.size1(); ++i)
		for(std::size_t j=i+1; j<mat.size2(); ++j)
			if(!float_equal(mat(i,j), mat(j,i), eps))
				return false;

	return true;
}


// -----------------------------------------------------------------------------
template<class T, LinalgType ty=get_linalg_type<T>::value> struct apply_fkt_impl {};

template<class T>
struct apply_fkt_impl<T, LinalgType::REAL>
{
	T operator()(T t, const std::function<T(T)>& fkt) const
	{
		return fkt(t);
	}
};


template<class t_vec>
struct apply_fkt_impl<t_vec, LinalgType::VECTOR>
{
	using value_type = underlying_value_type_t<t_vec>;

	t_vec operator()(const t_vec& vec, const std::function<value_type(value_type)>& fkt) const
	{
		t_vec v;
		v.resize(vec.size());

		for(std::size_t i=0; i<vec.size(); ++i)
			v[i] = fkt(vec[i]);

		return v;
	}
};


template<class t_mat>
struct apply_fkt_impl<t_mat, LinalgType::MATRIX>
{
	using value_type = underlying_value_type_t<t_mat>;

	t_mat operator()(const t_mat& mat, const std::function<value_type(value_type)>& fkt) const
	{
		t_mat m;
		m.resize(mat.size1(), mat.size2());

		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
				m(i,j) = fkt(mat(i,j));

		return m;
	}
};


template<class T, class t_val=underlying_value_type_t<T>>
T apply_fkt(const T& t, const std::function<t_val(t_val)>& fkt)
{
	apply_fkt_impl<T> impl;
	return impl(t, fkt);
}


template<class T, class t_val=underlying_value_type_t<T>>
inline T apply_fkt(const T& t, t_val(*pfkt)(t_val))
{
	std::function<t_val(t_val)> fkt(pfkt);
	return apply_fkt<T, t_val>(t, fkt);
}
// -----------------------------------------------------------------------------


template<class T, LinalgType ty=get_linalg_type<T>::value>
struct get_minmax_impl {};


template<class T>
struct get_minmax_impl<T, LinalgType::REAL>
{
	std::pair<T, T>
	operator()(T t) const
	{
		return std::pair<T,T>(t,t);
	}
};


template<class t_vec>
struct get_minmax_impl<t_vec, LinalgType::VECTOR>
{
	using t_val = underlying_value_type_t<t_vec>;

	std::pair<t_val, t_val>
	operator()(const t_vec& vec) const
	{
		t_val tmin = std::numeric_limits<t_val>::max();
		t_val tmax = -tmin;

		for(std::size_t i=0; i<vec.size(); ++i)
		{
			if(vec[i] < tmin) tmin = vec[i];
			if(vec[i] > tmax) tmax = vec[i];
		}

		return std::pair<t_val, t_val>(tmin, tmax);
	}
};


template<class t_mat>
struct get_minmax_impl<t_mat, LinalgType::MATRIX>
{
	using t_val = underlying_value_type_t<t_mat>;

	std::pair<t_val, t_val>
	operator()(const t_mat& mat) const
	{
		t_val tmin = std::numeric_limits<t_val>::max();
		t_val tmax = -tmin;

		for(std::size_t i=0; i<mat.size1(); ++i)
			for(std::size_t j=0; j<mat.size2(); ++j)
			{
				if(mat(i,j) < tmin) tmin = mat(i,j);
				if(mat(i,j) > tmax) tmax = mat(i,j);
			}

		return std::pair<t_val, t_val>(tmin, tmax);
	}
};


template<class T>
std::pair<underlying_value_type_t<T>, underlying_value_type_t<T>>
get_minmax(const T& t)
{
	get_minmax_impl<T> impl;
	return impl(t);
}

// -----------------------------------------------------------------------------


/**
 * calculates the dominant eigenvector/eigenvalue for symmetric matrices
 * @see (Desktop Bronstein 2008), equs. (4.148)-(4.151)
 * @see (Bronstein 2008), p. 324
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool eigenvec_dominant_sym(const t_mat& mat, t_vec& evec, T& eval,
	t_vec vecInit = tl::make_vec<t_vec>({1,0,0}),
	std::size_t iMaxIter = 50)
{
	if(mat.size1() != mat.size2())
	{
		log_err("Matrix ", mat, " is not square.");
		return false;
	}

#ifndef NDEBUG
	t_mat matAbs = apply_fkt(mat, std::function<T(T)>((T(*)(T))std::abs));
	T _dEps = get_minmax(matAbs).second / 100.;	// 1% accuracy
	if(!tl::is_symmetric(mat, _dEps)) log_warn("Matrix ", mat, " is not symmetric.");
#endif

	t_vec vecPrev;
	for(std::size_t iIter=0; iIter<iMaxIter; ++iIter)
	{
		if(iIter == iMaxIter-1)
			vecPrev = vecInit;
		vecInit = prod_mv(mat, vecInit);
	}

	const T normInit = veclen(vecInit);
	const T normPrev = veclen(vecPrev);

	eval = normInit / normPrev;
	evec = vecInit / normInit;
	return true;
}


/**
 * calculates the least dominant eigenvector/eigenvalue for symmetric matrices
 * @see (Desktop Bronstein 2008), equs. (4.148)-(4.151)
 * @see (Bronstein 2008), p. 324
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool eigenvec_least_dominant_sym(const t_mat& mat, t_vec& evec, T& eval,
	t_vec vecInit = tl::make_vec<t_vec>({1,0,0}),
	std::size_t iMaxIter = 50)
{
	t_mat M;
	if(!tl::inverse(mat, M))
		return false;

	if(!eigenvec_dominant_sym(M, evec, eval, vecInit, iMaxIter))
		return false;

	eval = T(1)/eval;
	return true;
}


/**
 * calculates the eigenvectors/eigenvalues for symmetric matrices
 * using the qr algorithm
 * ! for large matrices use eigenvec_sym  !
 * @see https://en.wikipedia.org/wiki/QR_algorithm
 */
template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool eigenvec_sym_simple(const t_mat& mat, std::vector<t_vec>& evecs, std::vector<T>& evals,
	std::size_t MAX_ITER = 1024, T tEps = std::cbrt(get_epsilon<T>()))
{
	if(mat.size1() != mat.size2())
	{
		log_err("Matrix ", mat, " is not square.");
		return false;
	}

#ifndef NDEBUG
	t_mat matAbs = apply_fkt(mat, std::function<T(T)>((T(*)(T))std::abs));
	T _dEps = get_minmax(matAbs).second / 100.;	// 1% accuracy
	if(!tl::is_symmetric(mat, _dEps))
		log_warn("Matrix ", mat, " is not symmetric.");
#endif

	const std::size_t n = mat.size1();
	t_mat I = unit_m<t_mat>(n);
	t_mat M = mat;

	bool bConverged = false;
	std::size_t iIter = 0;
	for(iIter = 0; iIter < MAX_ITER; ++iIter)
	{
		t_mat Q, R;
		if(!qr_decomp(M, Q, R))
		{
			log_err("QR decomposition failed for matrix ", M);
			return false;
		}

		t_mat Mlast = M;
		M = prod_mm(R, Q);
		I = prod_mm(I, Q);

		bConverged = true;
		for(std::size_t iVal = 0; iVal < n; ++iVal)
		{
			if(std::abs(M(iVal,iVal) - Mlast(iVal,iVal)) > tEps)
			{
				bConverged = false;
				break;
			}
		}
		if(bConverged)
			break;
	}

	evals.resize(n);
	evecs.resize(n);

	for(std::size_t iVal=0; iVal<n; ++iVal)
	{
		evals[iVal] = M(iVal, iVal);
		evecs[iVal] = get_column(I, iVal);
	}

	return bConverged;
}


template<class t_mat = ublas::matrix<double>,
	class t_vec = ublas::vector<typename t_mat::value_type>,
	typename T = typename t_mat::value_type>
bool eigenvec_approxsym_simple(const t_mat& mat, std::vector<t_vec>& evecs, std::vector<T>& evals,
	std::size_t MAX_ITER = 1024, T tEps = std::cbrt(get_epsilon<T>()))
{
	t_mat MtM = prod_mm(transpose(mat), mat);
	bool bOk = eigenvec_sym_simple(MtM, evecs, evals, MAX_ITER, tEps);

	for(T& eval : evals)
		eval = std::sqrt(std::abs(eval));
	return bOk;
}


template<typename T=double>
void sort_eigenvecs(std::vector<ublas::vector<T>>& evecs,
	std::vector<T>& evals, bool bOrder=0, T (*pEvalFkt)(T)=0,
	ublas::vector<T>* pUserVec = nullptr)
{
	if(evecs.size() != evals.size())
		return;

	struct Evec
	{
		ublas::vector<T> vec;
		T val;

		T userval = 0;
	};

	std::vector<Evec> myevecs;
	myevecs.reserve(evecs.size());

	for(std::size_t i=0; i<evecs.size(); ++i)
	{
		Evec ev;
		ev.vec = evecs[i];
		ev.val = evals[i];
		if(pUserVec) ev.userval = (*pUserVec)[i];

		myevecs.push_back(ev);
	}


	std::sort(myevecs.begin(), myevecs.end(),
		[&](const Evec& evec1, const Evec& evec2) -> bool
		{
			bool b;
			if(pEvalFkt)
				b = pEvalFkt(evec1.val) < pEvalFkt(evec2.val);
			else
				b = evec1.val < evec2.val;

			if(bOrder) b = !b;
			return b;
		});


	for(std::size_t i=0; i<evecs.size(); ++i)
	{
		evecs[i] = myevecs[i].vec;
		evals[i] = myevecs[i].val;
		if(pUserVec) (*pUserVec)[i] = myevecs[i].userval;
	}
}


// --------------------------------------------------------------------------------

/**
 * project vec1 onto vec2
 * proj_op = |vec2><vec2| / len(vec2)^2,  len(vec2) = sqrt(<vec2|vec2>)
 * proj = proj_op * vec1 = |vec2> * <vec2|vec1> / <vec2|vec2>
 */
template<typename t_vec = ublas::vector<double>>
t_vec proj_vec(const t_vec& vec1, const t_vec& vec2)
{
	using T = typename t_vec::value_type;

	T tnum = inner(vec1, vec2);
	T tden = inner(vec2, vec2);

	t_vec vecProj = tnum/tden * vec2;
	return vecProj;
}


/**
 * Gram-Schmidt orthogonalisation of basis vectors
 * @see e.g. https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
 * @see e.g. (Arens 2015), p. 744
 */
template<typename t_vec /*= ublas::vector<double>*/,
	typename T /*= typename t_vec::value_type*/ >
std::vector<t_vec> gram_schmidt(const std::vector<t_vec>& vecs, bool bNorm)
{
	std::vector<t_vec> vecsOut;
	if(vecs.size() == 0)
		return vecsOut;

	vecsOut.resize(vecs.size());

	// iterate through all basis vectors i
	for(std::size_t i=0; i<vecs.size(); ++i)
	{
		vecsOut[i] = vecs[i];

		// iterate through all previous basis vectors j<i
		// and remove projected contributions from i vectors
		for(std::size_t j=0; j<i; ++j)
			vecsOut[i] -= proj_vec<t_vec>(vecs[i], vecsOut[j]);
	}

	// normalise basis?
	if(bNorm)
		for(t_vec& vec : vecsOut)
			vec /= veclen(vec);

	return vecsOut;
}


template<typename t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type>
std::vector<t_vec> get_ortho_rhs(const std::vector<t_vec>& vecs)
{
	assert(vecs.size() == 2);

	std::vector<t_vec> vecOrtho = gram_schmidt(vecs, true);
	t_vec vecUp = cross_3(vecOrtho[0], vecOrtho[1]);
	vecOrtho.push_back(vecUp);

	return vecOrtho;
}

}

#endif
