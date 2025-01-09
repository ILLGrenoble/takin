/**
 * tlibs2 maths library -- operations with metric
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

#ifndef __TLIBS2_MATHS_TENSOR_H__
#define __TLIBS2_MATHS_TENSOR_H__

#include <cmath>
#include <vector>
#include <initializer_list>

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// operations with metric
// ----------------------------------------------------------------------------

/**
 * covariant metric tensor: g_{i,j} = e_i * e_j
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec, template<class...> class t_cont=std::initializer_list>
t_mat metric(const t_cont<t_vec>& basis_co)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = basis_co.size();

	t_mat g_co;
	if constexpr(is_dyn_mat<t_mat>)
		g_co = t_mat(N, N);

	auto iter_i = basis_co.begin();
	for(std::size_t i=0; i<N; ++i)
	{
		auto iter_j = basis_co.begin();
		for(std::size_t j=0; j<N; ++j)
		{
			g_co(i,j) = inner<t_vec>(*iter_i, *iter_j);
			std::advance(iter_j, 1);
		}
		std::advance(iter_i, 1);
	}

	return g_co;
}


/**
 * covariant metric tensor: g_{i,j} = e_i * e_j
 * @see (Arens 2015), p. 815
 */
template<class t_mat>
t_mat metric(const t_mat& basis_co)
requires is_basic_mat<t_mat>
{
	t_mat basis_trans = trans(basis_co);
	return basis_trans * basis_co;
}


/**
 * get levi-civita symbol in fractional coordinates
 * @see (Arens 2015), p. 815
 */
template<class t_mat>
typename t_mat::value_type levi(const t_mat& basis_co,
	const std::initializer_list<std::size_t>& indices)
requires is_basic_mat<t_mat>
{
	using size_t = decltype(basis_co.size1());
	using t_vec = vec<typename t_mat::value_type>;
	t_mat mat = create<t_mat>(basis_co.size1(), basis_co.size2());

	auto iter = indices.begin();
	size_t maxcols = std::min(indices.size(), basis_co.size2());
	for(size_t i=0; i<maxcols; ++i)
	{
		set_col<t_mat, t_vec>(
			mat, col<t_mat, t_vec>(basis_co, *iter), i);
		std::advance(iter, 1);
	}

	return det<t_mat>(mat);
}


/**
 * cross product in fractional coordinates: c^l = eps_ijk g^li a^j b^k
 * @see (Arens 2015), p. 815
 */
template<class t_mat, class t_vec>
t_vec cross(const t_mat& B, const t_vec& a, const t_vec& b)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using size_t = decltype(B.size1());
	using t_real = typename t_mat::value_type;

	// metric and its inverse
	auto G = metric<t_mat>(B);
	auto [G_inv, ok] = inv(G);

	// maximum indices
	const size_t i_max = G_inv.size1();
	const size_t j_max = a.size();
	const size_t k_max = b.size();
	const size_t l_max = G_inv.size2();

	// cross product result vector
	t_vec c = zero<t_vec>(l_max);

	for(size_t i=0; i<i_max; ++i)
	{
		for(size_t j=0; j<j_max; ++j)
		{
			for(size_t k=0; k<k_max; ++k)
			{
				t_real eps = levi<t_mat>(B, {i,j,k});
				for(size_t l=0; l<l_max; ++l)
					c[l] += eps * G_inv(l,i) * a[j] * b[k];
			}
		}
	}

	return c;
}


/**
 * lower index using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
t_vec lower_index(const t_mat& metric_co, const t_vec& vec_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = vec_contra.size();
	t_vec vec_co = zero<t_vec>(N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			vec_co[i] += metric_co(i,j) * vec_contra[j];

	return vec_co;
}


/**
 * raise index using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
t_vec raise_index(const t_mat& metric_contra, const t_vec& vec_co)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	const std::size_t N = vec_co.size();
	t_vec vec_contra = zero<t_vec>(N);

	for(std::size_t i=0; i<N; ++i)
		for(std::size_t j=0; j<N; ++j)
			vec_contra[i] += metric_contra(i,j) * vec_co[j];

	return vec_contra;
}


/**
 * inner product using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type inner(const t_mat& metric_co,
	const t_vec& vec1_contra, const t_vec& vec2_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	t_vec vec2_co = lower_index<t_mat, t_vec>(metric_co, vec2_contra);
	return inner<t_vec>(vec1_contra, vec2_co);
}


/**
 * 2-norm using metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type norm(const t_mat& metric_co, const t_vec& vec_contra)
requires is_basic_vec<t_vec>
{
	return std::sqrt(inner<t_mat, t_vec>(metric_co, vec_contra, vec_contra));
}


/**
 * angle between vectors under a given metric
 * @see (Arens 2015), p. 808
 */
template<class t_mat, class t_vec>
typename t_vec::value_type angle(const t_mat& metric_co,
	const t_vec& vec1_contra, const t_vec& vec2_contra)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	using t_real = typename t_mat::value_type;

	t_real len1 = norm<t_mat, t_vec>(metric_co, vec1_contra);
	t_real len2 = norm<t_mat, t_vec>(metric_co, vec2_contra);

	t_real c = inner<t_mat, t_vec>(metric_co, vec1_contra, vec2_contra);
	c /= len1 * len2;
	c = std::clamp<t_real>(c, -1, 1);

	return std::acos(c);
}

// ----------------------------------------------------------------------------

}

#endif
