/**
 * tlibs2 maths library -- (forward) declarations
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

#ifndef __TLIBS2_MATHS_DECLS_H__
#define __TLIBS2_MATHS_DECLS_H__

#include <tuple>
#include <vector>

#include "../traits.h"



// numbers
#if __has_include(<numbers>)
	#define __TLIBS2_USE_NUMBERS__
#endif



// lapack(e)
#if __has_include(<lapacke.h>) && USE_LAPACK
	#define __TLIBS2_USE_LAPACK__

	// ----------------------------------------------------------------------------
	// forward declarations
	// ----------------------------------------------------------------------------
	#define lapack_complex_double std::complex<double>
	#define lapack_complex_double_real(z) (z.real())
	#define lapack_complex_double_imag(z) (z.imag())

	#define lapack_complex_float std::complex<float>
	#define lapack_complex_float_real(z) (z.real())
	#define lapack_complex_float_imag(z) (z.imag())

	#include <lapacke.h>

namespace tl2_la {

	template<class t_mat, template<class...> class t_vec = std::vector>
	std::tuple<bool, t_vec<typename t_mat::value_type>, t_vec<lapack_int>> _lu_raw(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat, template<class...> class t_vec = std::vector>
	std::tuple<bool, t_mat, t_mat, t_mat> lu(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat>
	std::tuple<t_mat, bool> inv(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
	std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
	std::tuple<bool, t_mat> chol(const t_mat& mat, bool blocked = true)
	requires tl2::is_mat<t_mat>;

	template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>,
		class t_real = tl2::underlying_value_type<typename t_mat::value_type>>
	std::tuple<bool, t_mat, t_mat> chol_semi(const t_mat& mat, t_real eps)
	requires tl2::is_mat<t_mat>;

	template<class t_mat, class t_vec = std::vector<typename t_mat::value_type>>
	std::tuple<bool, t_mat, t_mat> chol_herm(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat_cplx, class t_vec_cplx, class t_cplx = typename t_mat_cplx::value_type,
		class t_real = typename t_cplx::value_type>
	std::tuple<bool, std::vector<t_cplx>, std::vector<t_vec_cplx>>
	eigenvec(const t_mat_cplx& mat,
		bool only_evals = false, bool is_hermitian = false, bool normalise = false,
		t_real mineval = -1, t_real maxeval = -2, t_real eps = -1)
		requires tl2::is_mat<t_mat_cplx> && tl2::is_vec<t_vec_cplx> && tl2::is_complex<t_cplx>;

	template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
	std::tuple<bool, std::vector<t_real>, std::vector<t_real>, std::vector<t_vec>, std::vector<t_vec>>
	eigenvec(const t_mat& mat, bool only_evals = false, bool is_symmetric = false, bool normalise = false,
		t_real mineval = -1, t_real maxeval = -2, t_real eps = -1)
		requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && !tl2::is_complex<t_real>);

	template<class t_mat, class t_scalar = typename t_mat::value_type,
		class t_real = tl2::underlying_value_type<t_scalar>>
	std::tuple<bool, t_mat, t_mat, std::vector<t_real>>
	singval(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat>
	std::tuple<t_mat, bool> pseudoinv(const t_mat& mat)
	requires tl2::is_mat<t_mat>;

	template<class t_mat, class t_vec, class t_val = typename t_vec::value_type>
	std::tuple<bool, t_vec>
	odesys_const(const t_mat& C, const t_val& x, const t_val& x0, const t_vec& f0)
	requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec>);

	template<class t_mat, class t_vec, class t_val = typename t_vec::value_type>
	std::tuple<bool, t_vec>
	diffsys_const(const t_mat& C, const t_val& n, const t_val& n0, const t_vec& f0)
	requires (tl2::is_mat<t_mat> && tl2::is_vec<t_vec>);
}

#else
	//#pragma message("tlibs2: Disabling Lapack(e) library (not found).")
#endif



// fadeeva
#if __has_include(<Faddeeva.hh>)
	#define __TLIBS2_USE_FADDEEVA__
#else
	//#pragma message("tlibs2: Disabling Faddeeva library (not found).")
#endif



// qhull
#if __has_include(<Qhull.h>)
	#define __TLIBS2_USE_QHULL__
#else
	//#pragma message("tlibs2: Disabling QHull library (not found).")
#endif



// separator tokens
#define TL2_COLSEP ';'
#define TL2_ROWSEP '|'



namespace tl2 {
// ----------------------------------------------------------------------------
// forward declarations
// ----------------------------------------------------------------------------
template<class t_mat>
t_mat prod(const t_mat& mat1, const t_mat& mat2, bool assert_sizes = true)
requires tl2::is_basic_mat<t_mat> && tl2::is_dyn_mat<t_mat>;

template<class t_mat>
t_mat zero(std::size_t N1, std::size_t N2)
requires is_basic_mat<t_mat>;

template<class t_vec>
t_vec zero(std::size_t N = 0)
requires is_basic_vec<t_vec>;

template<class t_mat>
t_mat unit(std::size_t N = 0)
requires is_basic_mat<t_mat>;

template<class t_mat>
t_mat unit(std::size_t N, std::size_t M)
requires is_basic_mat<t_mat>;

template<class t_mat, class t_vec>
std::tuple<bool, t_mat, t_mat> qr(const t_mat& mat)
requires is_mat<t_mat> && is_vec<t_vec>;

template<class t_mat>
std::tuple<t_mat, bool> inv(const t_mat& mat)
requires is_mat<t_mat>;

template<class t_mat>
typename t_mat::value_type det(const t_mat& mat)
requires is_mat<t_mat>;

template<class t_elem, template<class...> class t_cont = std::vector>
t_elem mean(const t_cont<t_elem>& vec)
requires is_basic_vec<t_cont<t_elem>>;
// ----------------------------------------------------------------------------


}

#endif
