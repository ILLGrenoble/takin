/**
 * tlibs2 maths library -- lapack(e) interface
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

#ifndef __TLIBS2_MATHS_LAPACKE__H__
#define __TLIBS2_MATHS_LAPACKE__H__

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

#include "decls.h"
#include "complex.h"
#include "ndim.h"



// ----------------------------------------------------------------------------
// lapack wrappers
// ----------------------------------------------------------------------------
#ifdef __TLIBS2_USE_LAPACK__

namespace tl2_la {

/**
 * LU decomposition of a matrix, mat = P * L * U, returning raw results
 * @returns [ok, LU, perm]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgetrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgetrf.html
 */
template<class t_mat, template<class...> class t_vec>
std::tuple<bool, t_vec<typename t_mat::value_type>, t_vec<lapack_int>> _lu_raw(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();
	const std::size_t minor = std::min(rows, cols);


	t_vec<t_scalar> outmat(rows*cols);
	t_vec<lapack_int> outpivots(minor);

	for(std::size_t i = 0; i < rows; ++i)
		for(std::size_t j = 0; j < cols; ++j)
			outmat[i*cols + j] = mat(i, j);

	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outpivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}

	return std::make_tuple(err == 0, outmat, outpivots);
}


/**
 * LU decomposition of a matrix, mat = P * L * U
 * @returns [ok, P, L, U]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgetrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgetrf.html
 */
template<class t_mat, template<class...> class t_vec>
std::tuple<bool, t_mat, t_mat, t_mat> lu(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();


	auto [ ok, lumat, pivots ] = _lu_raw<t_mat, t_vec>(mat);

	t_mat P = tl2::unit<t_mat>(rows, cols);
	t_mat L = tl2::unit<t_mat>(rows, cols);
	t_mat U = tl2::unit<t_mat>(rows, cols);

	// L and U
	for(std::size_t i = 0; i < rows; ++i)
	{
		for(std::size_t j = 0; j < cols; ++j)
		{
			if(j>=i)
				U(i, j) = lumat[i*cols + j];
			else
				L(i, j) = lumat[i*cols + j];
		}
	}

	// permutation matrix P
	for(std::size_t i = 0; i < pivots.size(); ++i)
		P = tl2::prod<t_mat>(P, tl2::perm<t_mat>(rows, cols, i, pivots[i]-1));

	return std::make_tuple(ok, P, L, U);
}


/**
 * inverted matrix
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgetri.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgetri.html
 */
template<class t_mat>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	// fail if matrix is not square
	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();

	t_mat I = tl2::unit<t_mat>(rows, cols);


	// lu factorisation
	auto [ ok, lumat, pivots ] = _lu_raw<t_mat, std::vector>(mat);
	if(!ok)
		return std::make_tuple(I, false);


	// inversion
	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
				rows, lumat.data(), rows,
				pivots.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}


	for(std::size_t i = 0; i < rows; ++i)
		for(std::size_t j = 0; j < cols; ++j)
			I(i, j) = lumat[i*cols + j];

	return std::make_tuple(I, err == 0);
}


/**
 * QR decomposition of a matrix, mat = QR
 * @returns [ok, Q, R]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgeqrf.html
 */
template<class t_mat, class t_vec>
std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();
	const std::size_t minor = std::min(rows, cols);

	const t_mat I = tl2::unit<t_mat>(minor);
	t_mat Q = I, R = mat;


	t_vec outmat(rows*cols), outvec(minor);

	for(std::size_t i = 0; i < rows; ++i)
		for(std::size_t j = 0; j < cols; ++j)
			outmat[i*cols + j] = mat(i, j);

	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,
				rows, cols, outmat.data(), cols,
				outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}

	for(std::size_t i = 0; i < rows; ++i)
		for(std::size_t j = 0; j < cols; ++j)
			R(i, j) = (j>=i ? outmat[i*cols + j] : t_real{0});


	t_vec v = tl2::zero<t_vec>(minor);

	for(std::size_t k=1; k<= minor; ++k)
	{
		for(std::size_t i=1; i <= k-1; ++i)
			v[i-1] = t_real{0};
		v[k-1] = t_real{1};

		for(std::size_t i=k+1; i <= minor; ++i)
			v[i-1] = outmat[(i-1)*cols + (k-1)];

		Q = Q * (I - outvec[k-1]*tl2::outer<t_mat, t_vec>(v, v));
	}

	return std::make_tuple(err == 0, Q, R);
}


/**
 * upper cholesky decomposition of a hermitian, positive-definite matrix, mat = C^H C
 * @returns [ok, C]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zpotrf.html
 * @see http://www.netlib.org/utk/papers/factor/node9.html
 */
template<class t_mat, class t_vec>
std::tuple<bool, t_mat> chol(const t_mat& mat, bool blocked)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	const std::size_t N = mat.size1();
	int err = -1;

	std::vector<t_scalar> outmat(N*N);

	for(std::size_t i = 0; i < N; ++i)
		for(std::size_t j = 0; j < N; ++j)
			outmat[i*N + j] = (j >= i ? mat(i, j) : 0.);

	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			if(blocked)
				err = LAPACKE_cpotrf(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
			else
				err = LAPACKE_cpotrf2(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			if(blocked)
				err = LAPACKE_zpotrf(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
			else
				err = LAPACKE_zpotrf2(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			if(blocked)
				err = LAPACKE_spotrf(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
			else
				err = LAPACKE_spotrf2(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			if(blocked)
				err = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
			else
				err = LAPACKE_dpotrf2(LAPACK_ROW_MAJOR,
					'U', N, outmat.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}

	t_mat C = tl2::create<t_mat>(N, N);

	for(std::size_t i = 0; i < N; ++i)
		for(std::size_t j = 0; j < N; ++j)
			C(i, j) = outmat[i*N + j];

	return std::make_tuple(err == 0, C);
}


/**
 * upper cholesky decomposition of a hermitian, positive-semidefinite matrix, mat = C^H C
 * @returns [ok, C, P]
 * @see https://netlib.org/lapack/explore-html/db/dba/zpstrf_8f_source.html
 */
template<class t_mat, class t_vec, class t_real>
std::tuple<bool, t_mat, t_mat> chol_semi(const t_mat& mat, t_real eps)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;

	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	const std::size_t N = mat.size1();
	int err = -1;

	std::vector<t_scalar> outmat(N*N);
	std::vector<int> pivot(N, 0);
	int steps = 0;

	for(std::size_t i = 0; i < N; ++i)
		for(std::size_t j = 0; j < N; ++j)
			outmat[i*N + j] = (j >= i ? mat(i, j) : 0.);

	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cpstrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N, pivot.data(), &steps, eps);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zpstrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N, pivot.data(), &steps, eps);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_spstrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N, pivot.data(), &steps, eps);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dpstrf(LAPACK_ROW_MAJOR,
				'U', N, outmat.data(), N, pivot.data(), &steps, eps);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}

	t_mat C = tl2::create<t_mat>(N, N);
	t_mat P = tl2::zero<t_mat>(N, N);    // pivot matrix

	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j < N; ++j)
			C(i, j) = outmat[i*N + j];

		P(pivot[i]-1, i) = 1.;
	}

	return std::make_tuple(err == 0, C, P);
}


/**
 * cholesky decomposition of a hermitian matrix, mat = C D C^H
 * @returns [ok, C, D]
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zhptrf.html
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dsptrf.html
 */
template<class t_mat, class t_vec>
std::tuple<bool, t_mat, t_mat> chol_herm(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using namespace tl2_ops;
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	if constexpr(tl2::is_dyn_mat<t_mat>)
		assert((mat.size1() == mat.size2()));
	else
		static_assert(t_mat::size1() == t_mat::size2());

	const std::size_t N = mat.size1();
	int err = -1;

	std::vector<t_scalar> outmat(N*(N+1)/2);
	std::vector<int> outvec(N);

	std::size_t lin_idx = 0;
	for(std::size_t i = 0; i < N; ++i)
		for(std::size_t j = 0; j <= i; ++j)
			outmat[lin_idx++] = mat(i, j);

	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_chptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zhptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_ssptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dsptrf(LAPACK_ROW_MAJOR,
				'L', N, outmat.data(), outvec.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}

	t_mat C = tl2::unit<t_mat>(N, N);
	t_mat D = tl2::unit<t_mat>(N);

	lin_idx = 0;
	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j <= i; ++j)
		{
			if(i == j)
				D(i, j) = outmat[lin_idx++];
			else
				C(i, j) = outmat[lin_idx++];
		}
	}

	// TODO: only correct for unit permutation so far...
	for(int& i : outvec)
		i = std::abs(i)-1;
	t_mat perm = tl2::perm<t_mat>(outvec);
	C = perm * C;

	return std::make_tuple(err == 0, C, D);
}


/**
 * eigenvectors and -values of a complex matrix
 * @returns [ok, evals, evecs]
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zheev.html
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgeev.html
 */
template<class t_mat_cplx, class t_vec_cplx, class t_cplx, class t_real>
std::tuple<bool, std::vector<t_cplx>, std::vector<t_vec_cplx>>
eigenvec(const t_mat_cplx& mat,
	bool only_evals, bool is_hermitian, bool normalise,
	t_real mineval, t_real maxeval, t_real eps)
	requires tl2::is_mat<t_mat_cplx> && tl2::is_vec<t_vec_cplx> && tl2::is_complex<t_cplx>
{
	bool only_selected_evals = (mineval <= maxeval);
	bool use_selective_func = only_selected_evals;
	//use_selective_func = true;

	std::vector<t_cplx> evals;
	std::vector<t_vec_cplx> evecs;

	if(mat.size1() != mat.size2() || mat.size1() == 0)
		return std::make_tuple(0, evals, evecs);

	const std::size_t N = mat.size1();
	evals.resize(N);

	if(!only_evals)
		evecs.resize(N, tl2::zero<t_vec_cplx>(N));

	std::vector<t_cplx> inmat(N*N, t_cplx{0,0}),
		outevecs(only_evals ? 0 : N*N, t_cplx{0,0});

	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j < N; ++j)
		{
			if(is_hermitian)
				inmat[i*N + j] = (j>=i ? mat(j,i) : t_real{0});
			else
				inmat[i*N + j] = mat(j,i);
		}
	}


	int err = -1;

	if(is_hermitian)
	{
		// evals of hermitian matrix are purely real
		std::vector<t_real> outevals_real(N, t_real{0});

		// all eigenvalues
		if(!use_selective_func)
		{
			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_cheev(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					'L', N, inmat.data(), N,
					outevals_real.data());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_zheev(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					'L', N, inmat.data(), N,
					outevals_real.data());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			}
		}

		// only selected eigenvalues
		else
		{
			int minidx = 1, maxidx = N;
			int iNumFound = 0;

			std::unique_ptr<int, std::default_delete<int[]>>
			uptrIdxArr(new int[2*N]);

			// use maximum precision if none given
			if(eps < t_real{0})
			{
				if constexpr(std::is_same_v<t_real, float>)
				{
					eps = LAPACKE_slamch('S');
				}
				else if constexpr(std::is_same_v<t_real, double>)
				{
					eps = LAPACKE_dlamch('S');
				}
				else
				{
					static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
				}
			}

			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_cheevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound,
					outevals_real.data(), outevecs.data(),
					N, uptrIdxArr.get());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_zheevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound,
					outevals_real.data(), outevecs.data(),
					N, uptrIdxArr.get());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			}

			// resize to actual number of eigenvalues and -vectors
			if(std::size_t(iNumFound) != N)
			{
				evals.resize(iNumFound, t_real{0});
				evecs.resize(iNumFound, tl2::zero<t_vec_cplx>(N));
			}
		}

		// copy to complex output vector
		for(std::size_t i = 0; i < evals.size(); ++i)
			evals[i] = outevals_real[i];
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N, evals.data(), nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N, evals.data(), nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}


	if(!only_evals)
	{
		for(std::size_t i = 0; i < evecs.size(); ++i)
		{
			// hermitian algo overwrites original matrix!
			for(std::size_t j = 0; j < N; ++j)
				evecs[i][j] = (is_hermitian && !use_selective_func)
					? inmat[i*N + j] : outevecs[i*N + j];

			if(normalise && (err == 0))
			{
				t_cplx n = tl2::norm(evecs[i]);
				if(!tl2::equals<t_cplx>(n, t_cplx{0,0}))
					evecs[i] /= n;
			}
		}
	}

	return std::make_tuple(err == 0, evals, evecs);
}


/**
 * eigenvectors and -values of a real matrix
 * @returns [ok, evals_re, evals_im, evecs_re, evecs_im]
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dsyev.html
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgeev.html
 */
template<class t_mat, class t_vec, class t_real>
std::tuple<bool, std::vector<t_real>, std::vector<t_real>, std::vector<t_vec>, std::vector<t_vec>>
eigenvec(const t_mat& mat, bool only_evals, bool is_symmetric, bool normalise,
	t_real mineval, t_real maxeval, t_real eps)
	requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && !tl2::is_complex<t_real>)
{
	bool only_selected_evals = (mineval <= maxeval);
	bool use_selective_func = only_selected_evals;
	//use_selective_func = true;

	std::vector<t_real> evals_re, evals_im;
	std::vector<t_vec> evecs_re, evecs_im;

	if(mat.size1() != mat.size2() || mat.size1() == 0)
		return std::make_tuple(false, evals_re, evals_im, evecs_re, evecs_im);

	const std::size_t N = mat.size1();
	evals_re.resize(N, t_real{0});
	evals_im.resize(N, t_real{0});

	if(!only_evals)
	{
		evecs_re.resize(N, tl2::zero<t_vec>(N));
		evecs_im.resize(N, tl2::zero<t_vec>(N));
	}

	std::vector<t_real> inmat(N*N, t_real{0}),
		outevecs(only_evals ? 0 : N*N, t_real{0});

	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j < N; ++j)
		{
			if(is_symmetric)
				inmat[i*N + j] = (j>=i ? mat(j,i) : t_real{0});
			else
				inmat[i*N + j] = mat(j,i);
		}
	}

	int err = -1;

	if(is_symmetric)
	{
		// all eigenvalues
		if(!use_selective_func)
		{
			// evals of symmetric matrix are purely real
			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_ssyev(LAPACK_COL_MAJOR,
					(only_evals ? 'N' : 'V'),
					'L', N, inmat.data(), N,
					evals_re.data());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_dsyev(LAPACK_COL_MAJOR,
					(only_evals ? 'N' : 'V'),
					'L', N, inmat.data(), N,
					evals_re.data());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			}
		}

		// only selected eigenvalues
		else
		{
			int minidx = 1, maxidx = N;
			int iNumFound = 0;

			std::unique_ptr<int, std::default_delete<int[]>>
			uptrIdxArr(new int[2*N]);

			// use maximum precision if none given
			if(eps < t_real{0})
			{
				if constexpr(std::is_same_v<t_real, float>)
				{
					eps = LAPACKE_slamch('S');
				}
				else if constexpr(std::is_same_v<t_real, double>)
				{
					eps = LAPACKE_dlamch('S');
				}
				else
				{
					static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
				}
			}

			if constexpr(std::is_same_v<t_real, float>)
			{
				err = LAPACKE_ssyevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound, evals_re.data(),
					outevecs.data(), N, uptrIdxArr.get());
			}
			else if constexpr(std::is_same_v<t_real, double>)
			{
				err = LAPACKE_dsyevr(LAPACK_COL_MAJOR,
					only_evals ? 'N' : 'V',
					only_selected_evals ? 'V' : 'A',
					'L', N, inmat.data(), N,
					mineval, maxeval, minidx, maxidx,
					eps, &iNumFound, evals_re.data(),
					outevecs.data(), N, uptrIdxArr.get());
			}
			else
			{
				static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
			}

			// resize to actual number of eigenvalues and -vectors
			if(std::size_t(iNumFound) != N)
			{
				evals_re.resize(iNumFound, t_real{0});
				evals_im.resize(iNumFound, t_real{0});
				evecs_re.resize(iNumFound, tl2::zero<t_vec>(N));
				evecs_im.resize(iNumFound, tl2::zero<t_vec>(N));
			}
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N,
				evals_re.data(), evals_im.data(),
				nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgeev(LAPACK_COL_MAJOR,
				'N', only_evals ? 'N' : 'V', N,
				inmat.data(), N,
				evals_re.data(), evals_im.data(),
				nullptr, N,
				only_evals ? nullptr : outevecs.data(), N);
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}


	// evecs
	if(!only_evals)
	{
		if((is_symmetric && !use_selective_func))
		{
			for(std::size_t i = 0; i < evals_re.size(); ++i)
			{
				// symmetric algo overwrites original matrix!
				for(std::size_t j = 0; j < N; ++j)
					evecs_re[i][j] = inmat[i*N + j];
			}
		}
		else
		{
			for(std::size_t i = 0; i < evals_re.size(); ++i)
			{
				if(tl2::equals<t_real>(evals_im[i], 0))
				{
					for(std::size_t j = 0; j < N; ++j)
					{
						evecs_re[i][j] = outevecs[i*N + j];
						evecs_im[i][j] = 0;
					}
				}
				else
				{
					for(std::size_t j = 0; j < N; ++j)
					{
						evecs_re[i][j] = outevecs[i*N + j];
						evecs_im[i][j] = outevecs[(i+1)*N + j]; // imag part of evec follows next in array

						evecs_re[i+1][j] = evecs_re[i][j];      // next evec is the conjugated one
						evecs_im[i+1][j] = -evecs_im[i][j];
					}
					++i;  // already used two values in array
				}
			}
		}

		if(normalise && (err == 0))
		{
			for(std::size_t i = 0; i < evecs_re.size(); ++i)
			{
				t_real sum{0};
				for(std::size_t j = 0; j < N; ++j)
					sum += std::norm(
						std::complex(evecs_re[i][j],
							evecs_im[i][j]));
				sum = std::sqrt(sum);

				if(!tl2::equals<t_real>(sum, 0))
				{
					evecs_re[i] /= sum;
					evecs_im[i] /= sum;
				}
			}
		}
	}

	return std::make_tuple(err == 0, evals_re, evals_im, evecs_re, evecs_im);
}


/**
 * singular values of a real or complex matrix mat = U * diag{vals} * V^h
 * @returns [ ok, U, Vh, vals ]
 * @see http://www.math.utah.edu/software/lapack/lapack-z/zgesvd.html
 * @see http://www.math.utah.edu/software/lapack/lapack-d/dgesvd.html
 */
template<class t_mat, class t_scalar, class t_real>
std::tuple<bool, t_mat, t_mat, std::vector<t_real>>
singval(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	const std::size_t rows = mat.size1();
	const std::size_t cols = mat.size2();

	const auto [Nmin, Nmax] = std::minmax(rows, cols);
	std::vector<t_scalar> inmat(rows*cols), outU(rows*rows), outVh(cols*cols);
	std::vector<t_real> vals(Nmin);
	std::vector<t_real> _tmp(Nmax * Nmax * 2);

	for(std::size_t i = 0; i < rows; ++i)
		for(std::size_t j = 0; j < cols; ++j)
			inmat[i*cols + j] = mat(i,j);

	int err = -1;
	if constexpr(tl2::is_complex<t_scalar>)
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_cgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols,
				inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(),
				cols, _tmp.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_zgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols,
				inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(),
				cols, _tmp.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}
	else
	{
		if constexpr(std::is_same_v<t_real, float>)
		{
			err = LAPACKE_sgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols, inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(), cols, _tmp.data());
		}
		else if constexpr(std::is_same_v<t_real, double>)
		{
			err = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,
				'A', 'A', rows, cols, inmat.data(), cols,
				vals.data(), outU.data(),
				rows, outVh.data(), cols, _tmp.data());
		}
		else
		{
			static_assert(tl2::bool_value<0, t_real>, "Invalid element type");
		}
	}


	t_mat U = tl2::unit<t_mat>(rows);
	t_mat Vh = tl2::unit<t_mat>(cols);

	for(std::size_t i = 0; i < Nmax; ++i)
	{
		for(std::size_t j = 0; j < Nmax; ++j)
		{
			if(i < U.size1() && j < U.size2()) U(i,j) = outU[i*cols + j];
			if(i < Vh.size1() && j < Vh.size2()) Vh(i,j) = outVh[i*cols + j];
		}
	}

	return std::make_tuple(err==0, U, Vh, vals);
}


/**
 * pseudoinverse M+ of a matrix
 * @see https://de.wikipedia.org/wiki/Pseudoinverse#Berechnung
 * @see (Arens 2015), pp. 788-792
 *
 * M  = U D (V*)^h
 * M+ = V D+ (U*)^h
 */
template<class t_mat>
std::tuple<t_mat, bool> pseudoinv(const t_mat& mat)
requires tl2::is_mat<t_mat>
{
	using t_scalar = typename t_mat::value_type;
	using t_real = tl2::underlying_value_type<t_scalar>;

	auto [ ok, U, Vh, vals ] = singval<t_mat>(mat);

	auto V = tl2::herm(Vh);
	auto Uh = tl2::herm(U);

	for(t_real& d : vals)
	{
		if(!tl2::equals<t_real>(d, t_real(0)))
			d = t_real(1)/d;
	}

	auto diag = tl2::diag<t_mat>(vals);
	return std::make_tuple(tl2::prod<t_mat>(V, tl2::prod(diag, Uh)), ok);
}



// ----------------------------------------------------------------------------
// equation solvers
// ----------------------------------------------------------------------------

/**
 * system of ODEs with constant coefficients C
 * @see (Arens 2015), pp. 1049-1051
 *
 * f'(x) = C f(x) and f(x0) = f0
 * => f(x) = f0 * exp(C(x-x0)) = sum_i norm_i * evec_i * exp(eval_i * (x-x0))
 *    norm = (evec_i)^(-1) * f0
 */
template<class t_mat, class t_vec, class t_val>
std::tuple<bool, t_vec>
odesys_const(const t_mat& C, const t_val& x, const t_val& x0, const t_vec& f0)
requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec>)
{
	const std::size_t rows = C.size1();
	const std::size_t cols = C.size2();
	if(rows != cols)
		return std::make_tuple(false, t_vec{});

	const auto [ok, evals, evecs] = eigenvec<t_mat, t_vec, t_val>(C);
	if(!ok)
		return std::make_tuple(false, t_vec{});

	const t_mat basis = tl2::create<t_mat, t_vec, std::vector>(evecs, false);
	const auto [basis_inv, invok] = tl2::inv<t_mat>(basis);
	if(!invok)
		return std::make_tuple(false, t_vec{});
	const t_vec norm = basis_inv * f0;

	t_vec f = tl2::zero<t_vec>(cols);
	for(std::size_t i = 0; i < cols; ++i)
		f += norm[i] * evecs[i] * std::exp(evals[i]*(x-x0));

	return std::make_tuple(true, f);
}


/**
 * system of difference equations with constant coefficients C
 */
template<class t_mat, class t_vec, class t_val>
std::tuple<bool, t_vec>
diffsys_const(const t_mat& C, const t_val& n, const t_val& n0, const t_vec& f0)
requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec>)
{
	const std::size_t rows = C.size1();
	const std::size_t cols = C.size2();
	if(rows != cols)
		return std::make_tuple(false, t_vec{});

	const auto [ok, evals, evecs] = eigenvec<t_mat, t_vec, t_val>(C);
	if(!ok)
		return std::make_tuple(false, t_vec{});

	const t_mat basis = tl2::create<t_mat, t_vec, std::vector>(evecs, false);
	const auto [basis_inv, invok] = tl2::inv<t_mat>(basis);
	if(!invok)
		return std::make_tuple(false, t_vec{});
	const t_vec norm = basis_inv * f0;

	t_vec f = tl2::zero<t_vec>(cols);
	for(std::size_t i = 0; i < cols; ++i)
		f += norm[i] * evecs[i] * std::pow(evals[i], n-n0);

	return std::make_tuple(true, f);
}
// ----------------------------------------------------------------------------

}

#endif	// __TLIBS2_USE_LAPACK__


#endif
