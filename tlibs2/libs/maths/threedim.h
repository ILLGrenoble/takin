/**
 * tlibs2 -- maths library -- 3-dim algos
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

#ifndef __TLIBS2_MATHS_3D_H__
#define __TLIBS2_MATHS_3D_H__

#include <cmath>
#include <tuple>
#include <vector>
#include <limits>

#include "decls.h"
#include "constants.h"
#include "projectors.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// 3-dim algos
// ----------------------------------------------------------------------------

/**
 * 3-dim cross product
 * @see https://en.wikipedia.org/wiki/Cross_product
 */
template<class t_vec>
t_vec cross(const t_vec& vec1, const t_vec& vec2)
requires is_basic_vec<t_vec>
{
	t_vec vec;

	// only valid for 3-vectors -> use first three components
	if(vec1.size() < 3 || vec2.size() < 3)
		return vec;

	if constexpr(is_dyn_vec<t_vec>)
		vec = t_vec(3);

	for(int i=0; i<3; ++i)
		vec[i] = vec1[(i+1)%3]*vec2[(i+2)%3] - vec1[(i+2)%3]*vec2[(i+1)%3];

	return vec;
}


/**
 * cross product matrix (3x3)
 * @see https://en.wikipedia.org/wiki/Skew-symmetric_matrix
 */
template<class t_mat, class t_vec>
t_mat skewsymmetric(const t_vec& vec)
requires is_basic_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat;
	if constexpr(is_dyn_mat<t_mat>)
		mat = t_mat(3,3);

	// if static matrix is larger than 3x3 (e.g. for homogeneous coordinates), initialise as identity
	if(mat.size1() > 3 || mat.size2() > 3)
		mat = unit<t_mat>(mat.size1(), mat.size2());

	mat(0,0) = 0;       mat(0,1) = -vec[2]; mat(0,2) = vec[1];
	mat(1,0) = vec[2];  mat(1,1) = 0;       mat(1,2) = -vec[0];
	mat(2,0) = -vec[1]; mat(2,1) = vec[0];  mat(2,2) = 0;

	return mat;
}


/**
 * givens rotation matrix
 * @see https://en.wikipedia.org/wiki/Givens_rotation
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat givens(std::size_t N, std::size_t i, std::size_t j, t_real angle)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(N, N);

	const t_real s = std::sin(angle);
	const t_real c = std::cos(angle);

	mat(i, j) = -s;
	mat(j, i) = +s;
	mat(i, i) = mat(j, j) = c;

	return mat;
}


/**
 * givens rotation matrix with precalculated sin and cos constants
 * @see https://en.wikipedia.org/wiki/Givens_rotation
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat givens(std::size_t N, std::size_t i, std::size_t j, t_real s, t_real c)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(N, N);

	mat(i, j) = -s;
	mat(j, i) = +s;
	mat(i, i) = mat(j, j) = c;

	return mat;
}


/**
 * SO(2) rotation matrix
 * @see https://en.wikipedia.org/wiki/Rotation_matrix
 */
template<class t_mat>
t_mat rotation_2d(const typename t_mat::value_type angle)
requires tl2::is_mat<t_mat>
{
	return givens<t_mat>(2, 0, 1, angle);
}


/**
 * SO(3) matrix to rotate around an axis (Rodrigues' formula)
 * @see (Arens 2015), p. 718 and p. 816
 * @see (Merziger 2006), p. 208
 * @see https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
 */
template<class t_mat, class t_vec>
t_mat rotation(const t_vec& axis, const typename t_vec::value_type angle,
	bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_real = typename t_vec::value_type;

	const t_real c = std::cos(angle);
	const t_real s = std::sin(angle);

	t_real len = 1;
	if(!is_normalised)
		len = tl2::norm<t_vec>(axis);

	// ----------------------------------------------------
	// special cases: rotations around [100], [010], [001]
	if(tl2::equals(axis, create<t_vec>({ len, 0, 0 })))
		return tl2::create<t_mat>({{1,0,0}, {0,c,s}, {0,-s,c}});
	else if(tl2::equals(axis, create<t_vec>({ 0, len, 0 })))
		return tl2::create<t_mat>({{c,0,-s}, {0,1,0}, {s,0,c}});
	else if(tl2::equals(axis, create<t_vec>({ 0, 0, len })))
		return tl2::create<t_mat>({{c,s,0}, {-s,c,0}, {0,0,1}});

	// ----------------------------------------------------
	// general case
	// project along rotation axis
	t_mat matProj1 = tl2::projector<t_mat, t_vec>(axis, is_normalised);

	// project along axis 2 in plane perpendicular to rotation axis
	t_mat matProj2 = tl2::ortho_projector<t_mat, t_vec>(axis, is_normalised) * c;

	// project along axis 3 in plane perpendicular to rotation axis and axis 2
	t_mat matProj3 = tl2::skewsymmetric<t_mat, t_vec>(axis/len) * s;

	t_mat matProj = matProj1 + matProj2 + matProj3;

	// if matrix is larger than 3x3 (e.g. for homogeneous cooridnates), fill up with identity
	tl2::unit<t_mat>(matProj, 3,3, matProj.size1(), matProj.size2());
	return matProj;
}


/**
 * get a vector perpendicular to vec
 */
template<class t_vec, class t_scalar = typename t_vec::value_type>
requires is_vec<t_vec>
t_vec perp(const t_vec& vec, t_scalar eps = std::numeric_limits<t_scalar>::epsilon())
{
	if(vec.size() == 2)
	{
		t_vec vecret = create<t_vec>({vec[1], vec[0]});
		if(std::abs(vecret[0]) > std::abs(vecret[1]))
			vecret[0] = -vecret[0];
	}

	else if(vec.size() == 3 || vec.size() == 4)
	{
		while(1)
		{
			t_vec rand = tl2::rand<t_vec>(3);
			t_vec perp = tl2::cross<t_vec>({ vec, rand });
			t_scalar dot = tl2::inner<t_vec>(perp, perp);

			if(dot > eps)
			{
				perp /= std::sqrt(dot);
				return perp;
			}
		}
	}

	return t_vec{};
}


/**
 * matrix to rotate a vector into [1, 0, 0, ...]
 * @see (Zhelezov 2017) O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
t_mat rotation_x0(const t_vec& _vec, t_real eps = 1e-6)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_size = decltype(_vec.size());
	const t_size N = _vec.size();

	t_vec vec = _vec;
	t_mat mat = unit<t_mat>(N);

	// equations (6) and (7) from (Zhelezov 2017)
	for(t_size n=1; n<N; ++n)
	{
		t_size i = N - n - 1;
		t_size j = i + 1;

		t_real len = vec[i]*vec[i] + vec[j]*vec[j];
		if(!equals_0<t_real>(len, eps))
		{
			len = std::sqrt(len);
			t_real s = -vec[j] / len;
			t_real c = vec[i] / len;

			t_mat rot_ij = givens<t_mat, t_real>(N, i, j, s, c);
			vec = rot_ij * vec;
			mat = rot_ij * mat;
		}
	}

	return mat;
}


/**
 * signed angle between two vectors
 * @see (Zhelezov 2017) O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<typename t_vec>
typename t_vec::value_type angle(const t_vec& vec0, const t_vec& vec1,
	const t_vec* pvec_norm = nullptr, t_vec* rotaxis = nullptr,
	bool force_3dim = true)
requires is_vec<t_vec>
{
	using namespace tl2_ops;
	using t_real = typename t_vec::value_type;

	assert(vec0.size() == vec1.size());

	// 2-dim case
	if(vec0.size() == 2)
	{
		// signed angles wrt basis
		t_real angle0 = std::atan2(vec0[1], vec0[0]);
		t_real angle1 = std::atan2(vec1[1], vec1[0]);

		return angle1 - angle0;
	}

	// 3-dim and 4-dim homogeneous case
	else if(vec0.size() == 3 || (vec0.size() == 4 && force_3dim))
	{
		// cross product gives sine and rotation axis
		t_vec axis = cross<t_vec>({ vec0, vec1 });
		t_real s = tl2::norm(axis);
		if(rotaxis)
			*rotaxis = axis / s;

		// dot product gives cosine
		t_real c = inner(vec0, vec1);
		t_real angle = std::atan2(s, c);

		// get signed angle
		if(pvec_norm)
		{
			// see if the cross product points along the direction
			// of the given normal
			if(inner(axis, *pvec_norm) < t_real{0})
				angle = -angle;
		}

		return angle;
	}

	// only implemented for two- or three-dimensional vectors
	assert(false);
	return t_real(0);
}


/**
 * matrix to rotate vector vec1 into vec2
 * (the vectors do not need to be normalised)
 * @see (Zhelezov 2017) O. I. Zhelezov, American Journal of Computational and Applied Mathematics 7(2), pp. 51-57 (2017), doi: 10.5923/j.ajcam.20170702.04
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
t_mat rotation(const t_vec& vec1, const t_vec& vec2,
	const t_vec *ortho_vec = nullptr, t_real eps = 1e-6, bool force_3dim = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	assert(vec1.size() == vec2.size());

	using t_size = decltype(vec1.size());
	const t_size dim = vec1.size();

	// get rotation axis and angle
	t_vec axis{};
	t_real angle{};
	if(dim == 2 || dim == 3 || (dim == 4 && force_3dim))
		angle = tl2::angle<t_vec>(vec1, vec2, nullptr, &axis, force_3dim);

	// 2-dim case
	if(dim == 2)
	{
		return rotation_2d<t_mat>(angle);
	}

	// 3-dim and 4-dim homogeneous case
	else if(dim == 3 || (dim == 4 && force_3dim))
	{
		// collinear vectors?
		if(equals<t_real>(angle, 0, eps))
			return tl2::unit<t_mat>(vec1.size());

		// antiparallel vectors?
		if(equals<t_real>(std::abs(angle), pi<t_real>, eps))
		{
			if(ortho_vec)
			{
				axis = *ortho_vec;
				axis /= tl2::norm<t_vec>(axis);
			}
			else
			{
				axis = tl2::perp<t_vec>(vec1, eps);
			}

			angle = pi<t_real>;
		}

		return tl2::rotation<t_mat, t_vec>(axis, angle, true);
	}

	// general case, equation (8) from (Zhelezov 2017)
	else
	{
		// matrix to rotate vec1 to [1, 0, 0, ...]
		t_mat mat1_x0 = tl2::rotation_x0<t_mat, t_vec, t_real>(vec1, eps);

		// matrix to rotate vec2 to [1, 0, 0, ...]
		t_mat mat2_x0 = tl2::rotation_x0<t_mat, t_vec, t_real>(vec2, eps);

		// matrix to rotate [1, 0, 0, ...] to vec2
		auto [mat2_x0_inv, mat2_ok] = tl2::inv<t_mat>(mat2_x0);

		// rotate vec1 to [1, 0, 0, ...], then rotate [1, 0, 0, ...] to vec2
		return mat2_x0_inv * mat1_x0;
	}
}


/**
 * transforms vertices and normals using a matrix
 */
template<class t_mat, class t_vec, template<class...> class t_cont = std::vector>
void transform_obj(t_cont<t_vec>& verts, t_cont<t_vec>& norms, const t_mat& mat, bool is_3dhom = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using size_t = decltype(mat.size1());

	// make sure a 3-vector and a 4-matrix are handled correctly in homogeneous coordinates
	if(is_3dhom && mat.size1() == 4)
	{
		t_mat mat3 = mat;
		for(size_t i = 0; i < 3; ++i)
		{
			mat3(3,i) = 0;
			mat3(i,3) = 0;
		}
		mat3(3,3) = 1;

		for(auto& vert : verts)
		{
			vert = mat3 * vert;

			// add translation and normalise
			for(size_t i = 0; i < 3; ++i)
			{
				vert[i] += mat(i,3);
				vert[i] /= mat(3,3);
			}
		}

		for(auto& norm : norms)
			norm = mat3 * norm;
	}

	// standard case: just multiply
	else
	{
		for(auto& vert : verts)
			vert = mat * vert;

		for(auto& norm : norms)
			norm = mat * norm;
	}
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// 3-dim algos in homogeneous coordinates
// ----------------------------------------------------------------------------

/**
 * does the homogeneous matrix mat have a translation component?
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool hom_has_translation_components(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	const std::size_t N = mat.size1();
	if(N != mat.size2())
		return false;

	// translation?
	for(std::size_t i=0; i<N-1; ++i)
	{
		if(!equals<t_real>(mat(i, N-1), t_real(0), eps))
			return true;
	}

	return false;
}


/**
 * is mat a centering matrix in homogeneous coords?
 */
template<class t_mat, class t_real = typename t_mat::value_type>
bool hom_is_centring(const t_mat& mat,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_mat<t_mat>
{
	const std::size_t N = mat.size1();
	if(N != mat.size2())
		return false;

	// is the left-upper 3x3 a rotation matrix (and no unit matrix)?
	if(!is_unit(submat<t_mat>(mat, 0, 0, 3, 3), eps))
		return false;

	// translation?
	if(hom_has_translation_components<t_mat, t_real>(mat, eps))
		return true;

	return false;
}


/**
 * project a homogeneous vector to screen coordinates
 * @returns [vecPersp, vecScreen]
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluProject.xml
 */
template<class t_mat, class t_vec>
std::tuple<t_vec, t_vec> hom_to_screen_coords(const t_vec& vec4,
	const t_mat& matModelView, const t_mat& matProj, const t_mat& matViewport,
	bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	// perspective trafo and divide
	t_vec vecPersp = matProj * matModelView * vec4;
	vecPersp /= vecPersp[3];

	// viewport trafo
	t_vec vec = matViewport * vecPersp;

	// flip y coordinate
	if(bFlipY) vec[1] = matViewport(1,1)*2 - vec[1];
	// flip x coordinate
	if(bFlipX) vec[0] = matViewport(0,0)*2 - vec[0];

	return std::make_tuple(vecPersp, vec);
}


/**
 * calculate world coordinates from screen coordinates
 * (vary zPlane to get the points of the z-line at constant (x,y))
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluUnProject.xml
 */
template<class t_mat, class t_vec>
t_vec hom_from_screen_coords(
	typename t_vec::value_type xScreen, typename t_vec::value_type yScreen, typename t_vec::value_type zPlane,
	const t_mat& matModelView_inv, const t_mat& matProj_inv, const t_mat& matViewport_inv,
	const t_mat* pmatViewport = nullptr, bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_vec vecScreen = create<t_vec>({xScreen, yScreen, zPlane, 1.});

	// flip y coordinate
	if(pmatViewport && bFlipY) vecScreen[1] = (*pmatViewport)(1,1)*2 - vecScreen[1];
	// flip x coordinate
	if(pmatViewport && bFlipX) vecScreen[0] = (*pmatViewport)(0,0)*2 - vecScreen[0];

	t_vec vecWorld = matModelView_inv * matProj_inv * matViewport_inv * vecScreen;

	vecWorld /= vecWorld[3];
	return vecWorld;
}


/**
 * calculate line from screen coordinates
 * @returns [pos, dir]
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluUnProject.xml
 */
template<class t_mat, class t_vec>
std::tuple<t_vec, t_vec> hom_line_from_screen_coords(
	typename t_vec::value_type xScreen, typename t_vec::value_type yScreen,
	typename t_vec::value_type z1, typename t_vec::value_type z2,
	const t_mat& matModelView_inv, const t_mat& matProj_inv, const t_mat& matViewport_inv,
	const t_mat* pmatViewport = nullptr, bool bFlipY = false, bool bFlipX = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	const t_vec lineOrg = tl2::hom_from_screen_coords<t_mat, t_vec>(
		xScreen, yScreen, z1, matModelView_inv, matProj_inv,
		matViewport_inv, pmatViewport, bFlipY, bFlipX);
	const t_vec linePos2 = tl2::hom_from_screen_coords<t_mat, t_vec>(
		xScreen, yScreen, z2, matModelView_inv, matProj_inv,
		matViewport_inv, pmatViewport, bFlipY, bFlipX);

	t_vec lineDir = linePos2 - lineOrg;
	lineDir /= tl2::norm<t_vec>(lineDir);

	return std::make_tuple(lineOrg, lineDir);
}


/**
 * perspective matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_perspective(
	t_scalar n = 0.01, t_scalar f = 100.,
	t_scalar fov = 0.5*pi<typename t_mat::value_type>, t_scalar ratio = 3./4.,
	bool bLHS = true, bool bZ01 = false)
requires is_mat<t_mat>
{
	const t_scalar c = 1./std::tan(0.5 * fov);
	const t_scalar n0 = bZ01 ? t_scalar{0} : n;
	const t_scalar sc = bZ01 ? t_scalar{1} : t_scalar{2};
	const t_scalar zs = bLHS ? t_scalar{-1} : t_scalar{1};
	const t_scalar range_nf = std::abs(f-n);

	//         ( x*c*r                           )      ( -x*c*r/z                         )
	//         ( y*c                             )      ( -y*c/z                           )
	// P * x = ( z*(n0+f)/(n-f) + w*sc*n*f/(n-f) )  =>  ( -(n0+f)/(n-f) - w/z*sc*n*f/(n-f) )
	//         ( -z                              )      ( 1                                )
	return create<t_mat>({
		c*ratio,  0.,  0.,                   0.,
		0,        c,   0.,                   0.,
		0.,       0.,  zs*(n0+f)/range_nf,   -sc*n*f/range_nf,
		0.,       0.,  zs,                   0.
	});
}


/**
 * orthographic projection matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
 * @see https://en.wikipedia.org/wiki/Orthographic_projection
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_ortho_sym(
	t_scalar n = 0.01, t_scalar f = 100.,
	t_scalar range_lr = 2., t_scalar range_bt = 2.,
	bool bLHS = true, bool bMap05 = false)
requires is_mat<t_mat>
{
	// map ranges into [-0.5, 0.5] or [-1, 1] else
	const t_scalar sc = bMap05 ? t_scalar{1} : t_scalar{2};
	const t_scalar zs = bLHS ? t_scalar{-1} : t_scalar{1};

	const t_scalar range_nf = std::abs(f) + std::abs(n);

	// centring
	const t_scalar tr_z = sc*t_scalar{0.5} * (n+f) / range_nf;

	//         ( sc_x*x        )
	//         ( sc_y*y        )
	// P * x = ( sc_z*z - tr_z )
	//         ( 1             )
	return create<t_mat>({
		sc/range_lr, 0.,          0.,              0.,
		0.,          sc/range_bt, 0.,              0.,
		0.,          0.,          zs*sc/range_nf,  zs*tr_z,
		0.,          0.,          0.,              1.
	});
}


/**
 * orthographic projection matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml
 * @see https://en.wikipedia.org/wiki/Orthographic_projection
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_ortho(
	t_scalar n = 0.01, t_scalar f = 100.,
	t_scalar l = -1., t_scalar r = 1.,
	t_scalar b = -1., t_scalar t = 1.,
	bool bLHS = true, bool bMap05 = false)
requires is_mat<t_mat>
{
	// map ranges into [-0.5, 0.5] or [-1, 1] else
	const t_scalar sc = bMap05 ? t_scalar{1} : t_scalar{2};
	const t_scalar zs = bLHS ? t_scalar{-1} : t_scalar{1};

	const t_scalar range_lr = r - l;
	const t_scalar range_bt = t - b;
	const t_scalar range_nf = f - n;

	// centring
	const t_scalar tr_x = sc*t_scalar{0.5} * (l+r) / range_lr;
	const t_scalar tr_y = sc*t_scalar{0.5} * (b+t) / range_bt;
	const t_scalar tr_z = sc*t_scalar{0.5} * (n+f) / range_nf;

	//         ( sc_x*x - tr_x )
	//         ( sc_y*y - tr_y )
	// P * x = ( sc_z*z - tr_z )
	//         ( 1             )
	return create<t_mat>({
		sc/range_lr, 0.,          0.,             -tr_x,
		0.,          sc/range_bt, 0.,             -tr_y,
		0.,          0.,          zs*sc/range_nf, zs*tr_z,
		0.,          0.,          0.,             1.
	});
}


/**
 * viewport matrix (homogeneous 4x4)
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glViewport.xml
 */
template<class t_mat, typename t_scalar=typename t_mat::value_type>
t_mat hom_viewport(t_scalar w, t_scalar h, t_scalar n = 0, t_scalar f = 1)
requires is_mat<t_mat>
{
	t_scalar sc{0.5};

	return create<t_mat>({
		sc*w,  0.,    0.,        sc*w,
		0,     sc*h,  0.,        sc*h,
		0.,    0.,    sc*(f-n),  sc*(f+n),
		0.,    0.,    0.,        1.
	});
}


/**
 * translation matrix in homogeneous coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat hom_translation(t_real x, t_real y, t_real z)
requires is_mat<t_mat>
{
	return create<t_mat>({
		1.,  0.,  0.,  x,
		0.,  1.,  0.,  y,
		0.,  0.,  1.,  z,
		0.,  0.,  0.,  1.
	});
}


/**
 * scaling matrix in homogeneous coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat hom_scaling(t_real x, t_real y, t_real z)
requires is_mat<t_mat>
{
	return create<t_mat>({
		x,   0.,  0.,  0.,
		0.,  y,   0.,  0.,
		0.,  0.,  z,   0.,
		0.,  0.,  0.,  1.
	});
}


/**
 * mirror matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_mirror(const t_vec& axis, bool is_normalised=1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat = ortho_mirror_op<t_mat, t_vec>(axis, is_normalised);

	t_mat mat_hom = unit<t_mat>(4,4);
	tl2::convert<t_mat, t_mat>(mat_hom, mat);

	// in case the matrix is statically sized and was already larger than 4x4
	mat_hom(3, 3) = 1.;
	return mat_hom;
}


/**
 * mirror matrix in homogeneous coordinates with translation
 */
template<class t_mat, class t_vec>
t_mat hom_mirror(const t_vec& axis, const t_vec& pos, bool is_normalised=1)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mirr = hom_mirror<t_mat, t_vec>(axis, is_normalised);

	t_mat offs = hom_translation<t_mat, t_vec>(pos);
	t_mat offs_inv = hom_translation<t_mat, t_vec>(-pos);

	return offs * mirr * offs_inv;
}


/**
 * "look at" matrix in homogeneous coordinates
 * @see (Sellers 2014), pp. 78-79
 * @see https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluLookAt.xml
 */
template<class t_mat, class t_vec>
t_mat hom_lookat(const t_vec& pos, const t_vec& target, const t_vec& _up)
requires is_vec<t_vec> && is_mat<t_mat>
{
	using t_real = typename t_mat::value_type;

	// create orthonormal system
	t_vec dir = -(target - pos);
	dir = dir / tl2::norm<t_vec>(dir);

	t_vec side = tl2::cross<t_vec>({_up, dir});
	side = side / tl2::norm<t_vec>(side);

	t_vec up = cross<t_vec>({dir, side});
	//up = up / tl2::norm<t_vec>(up);

	// inverted/transposed rotation matrix
	t_mat rot_inv = unit<t_mat>(4);
	tl2::set_row<t_mat, t_vec>(rot_inv, side, 0);
	tl2::set_row<t_mat, t_vec>(rot_inv, up, 1);
	tl2::set_row<t_mat, t_vec>(rot_inv, dir, 2);

	// inverted translation matrix
	t_mat trans_inv = hom_translation<t_mat, t_real>(
		-pos[0], -pos[1], -pos[2]);

	return rot_inv * trans_inv;
}


/**
 * shear matrix
 * @see https://en.wikipedia.org/wiki/Shear_matrix
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat shear(std::size_t ROWS, std::size_t COLS, std::size_t i, std::size_t j, t_real s)
requires is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(ROWS, COLS);
	mat(i,j) = s;
	return mat;
}


/**
 * rotation matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_rotation(const t_vec& vec1, const t_vec& vec2,
	const t_vec *perp_vec = nullptr)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat rot = tl2::rotation<t_mat, t_vec>(vec1, vec2, perp_vec);

	return tl2::create<t_mat>({
		rot(0,0), rot(0,1), rot(0,2), 0.,
		rot(1,0), rot(1,1), rot(1,2), 0.,
		rot(2,0), rot(2,1), rot(2,2), 0.,
		0.,       0.,       0.,       1.
	});
}


/**
 * rotation matrix in homogeneous coordinates
 */
template<class t_mat, class t_vec>
t_mat hom_rotation(const t_vec& axis, typename t_vec::value_type angle,
	bool is_normalised = true)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat rot = tl2::rotation<t_mat, t_vec>(axis, angle, is_normalised);

	return tl2::create<t_mat>({
		rot(0,0), rot(0,1), rot(0,2), 0.,
		rot(1,0), rot(1,1), rot(1,2), 0.,
		rot(2,0), rot(2,1), rot(2,2), 0.,
		0.,       0.,       0.,       1.
	});
}

// ----------------------------------------------------------------------------



/**
 * arrow matrix
 */
template<class t_vec, class t_mat, class t_real = typename t_vec::value_type>
t_mat get_arrow_matrix(
	const t_vec& vecTo,
	t_real postscale = 1, const t_vec& vecPostTrans = create<t_vec>({0,0,0.5}),
	const t_vec& vecFrom = create<t_vec>({0,0,1}),
	t_real prescale = 1, const t_vec& vecPreTrans = create<t_vec>({0,0,0}),
	const t_vec* perp_vec = nullptr)
requires is_vec<t_vec> && is_mat<t_mat>
{
	t_mat mat = unit<t_mat>(4);

	mat *= hom_translation<t_mat>(vecPreTrans[0], vecPreTrans[1], vecPreTrans[2]);
	mat *= hom_scaling<t_mat>(prescale, prescale, prescale);

	mat *= hom_rotation<t_mat, t_vec>(vecFrom, vecTo, perp_vec);

	mat *= hom_scaling<t_mat>(postscale, postscale, postscale);
	mat *= hom_translation<t_mat>(vecPostTrans[0], vecPostTrans[1], vecPostTrans[2]);

	return mat;
}



/**
 * real crystallographic A matrix
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 *
 * derivation from dot products:
 *   (1) <a|b> = ab cos(_cc), (2) <a|c> = ac cos(_bb), (3) <b|c> = bc cos(_aa),
 *   (4) <a|a> = a^2,         (5) <b|b> = b^2,         (6) <c|c> = c^2
 *
 *   (1) => b_0 = b cos(_cc)
 *   (5) => b_1 = b sin(_cc)

 *   (2) => c_0 = c cos(_bb)
 *   (3) => b_0*c_0 + b_1*c_1 = bc cos(_aa)
 *       => bc cos(_cc)*cos(_bb) + b sin(_cc) c_1 = bc cos(_aa)
 *       => c_1 = (c cos(_aa) - c cos(_cc)*cos(_bb)) / sin(_cc)
 *   (6) => c_2^2 = c^2 - c_0^2 - c_1^2
 *       => c_2^2 = c^2 (sin^2(_bb) - (cos(_aa) - cos(_cc)*cos(_bb))^2 / sin^2(_cc))
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat A_matrix(t_real a, t_real b, t_real c, t_real _aa, t_real _bb, t_real _cc)
requires is_mat<t_mat>
{
	const t_real ca = std::cos(_aa);
	const t_real cb = std::cos(_bb);
	const t_real cc = std::cos(_cc);
	const t_real sc = std::sin(_cc);
	const t_real sb = std::sin(_bb);

	return create<t_mat>({
		a,         b*cc,      c*cb,
		t_real{0}, b*sc,      c*(ca - cc*cb)/sc,
		t_real{0}, t_real{0}, c*std::sqrt(sb*sb - std::pow((ca - cc*cb)/sc, t_real{2}))
	});
}


/**
 * reciprocal crystallographic B matrix, B = 2pi * A^(-T)
 * Q [1/A] = B * Q [rlu]
 * @see https://en.wikipedia.org/wiki/Fractional_coordinates
 */
template<class t_mat, class t_real = typename t_mat::value_type>
t_mat B_matrix(t_real a, t_real b, t_real c, t_real _aa, t_real _bb, t_real _cc)
requires is_mat<t_mat>
{
	const t_real sc = std::sin(_cc);
	const t_real ca = std::cos(_aa);
	const t_real cb = std::cos(_bb);
	const t_real cc = std::cos(_cc);
	const t_real rr = std::sqrt(
		t_real{1} + t_real{2}*ca*cb*cc - (ca*ca + cb*cb + cc*cc));

	return t_real{2}*pi<t_real> * create<t_mat>({
		t_real{1}/a,             t_real{0},                   t_real{0},
		-t_real{1}/a * cc/sc,    t_real{1}/b * t_real{1}/sc,  t_real{0},
		(cc*ca - cb)/(a*sc*rr),  (cb*cc - ca)/(b*sc*rr),      sc/(c*rr)
	});
}


/**
 * UB orientation matrix, Q_exp [1/A] = U*B * Q [rlu]
 * @see https://dx.doi.org/10.1107/S0021889805004875
 */
template<class t_mat, class t_vec>
t_mat UB_matrix(const t_mat& B,
	const t_vec& vec1_rlu, const t_vec& vec2_rlu, const t_vec& vec3_rlu)
requires is_mat<t_mat> && is_vec<t_vec>
{
	t_vec vec1_lab = B * vec1_rlu;
	t_vec vec2_lab = B * vec2_rlu;
	t_vec vec3_lab = B * vec3_rlu;

	vec1_lab /= norm<t_vec>(vec1_lab);
	vec2_lab /= norm<t_vec>(vec2_lab);
	vec3_lab /= norm<t_vec>(vec3_lab);

	t_mat U_lab = unit<t_mat>(B.size1(), B.size2());
	set_row<t_mat, t_vec>(U_lab, vec1_lab, 0);
	set_row<t_mat, t_vec>(U_lab, vec2_lab, 1);
	set_row<t_mat, t_vec>(U_lab, vec3_lab, 2);

	return U_lab * B;
}


/**
 * get correct distance in unit cell, considering wrapping-around
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
t_real get_dist_uc(const t_mat& matA, const t_vec& vec1, const t_vec& vec2)
requires is_mat<t_mat> && is_vec<t_vec>
{
	t_vec vec1A = matA * vec1;

	// all supercell position to try
	std::vector<t_vec> vecSCs
	{{
		create<t_vec>({0,  0,  0}),

		create<t_vec>({+1,  0,  0}), create<t_vec>({ 0, +1,  0}), create<t_vec>({ 0,  0, +1}),

		create<t_vec>({ 0, +1, +1}), create<t_vec>({ 0, +1, -1}),
		create<t_vec>({+1,  0, +1}), create<t_vec>({+1,  0, -1}),
		create<t_vec>({+1, +1,  0}), create<t_vec>({+1, -1,  0}),

		create<t_vec>({+1, +1, +1}), create<t_vec>({+1, +1, -1}),
		create<t_vec>({+1, -1, +1}), create<t_vec>({+1, -1, -1}),
	}};

	t_real thedist = std::numeric_limits<t_real>::max();

	for(const t_vec& _vecSC : vecSCs)
	{
		for(const t_vec& vecSC : { _vecSC, -_vecSC, t_real{2}*_vecSC, t_real{-2}*_vecSC })
		{
			t_vec vec2A = matA * (vec2 + vecSC);
			t_real dist = tl2::norm(vec1A - vec2A);

			thedist = std::min(thedist, dist);
		}
	}

	return thedist;
}


/**
 * wrap atom positions back to unit cell
 */
template<class t_vec, class t_real = typename t_vec::value_type>
t_vec keep_atom_in_uc(const t_vec& _atom, std::size_t dim = 3)
requires is_vec<t_vec>
{
	t_vec newatom = _atom;

	for(std::size_t i=0; i</*newatom.size()*/ dim; ++i)
	{
		newatom[i] = std::fmod(newatom[i], t_real{1});
		if(newatom[i] <= t_real{-0.5})
			newatom[i] += std::abs(std::floor(newatom[i]));
		if(newatom[i] > t_real{0.5})
			newatom[i] -= std::abs(std::ceil(newatom[i]));
	}

	return newatom;
}



/**
 * wrap collection of atom positions back to unit cell
 */
template<class t_vec, class t_real = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
t_cont<t_vec> keep_atoms_in_uc(const t_cont<t_vec>& _atoms, std::size_t dim = 3)
requires is_vec<t_vec>
{
	t_cont<t_vec> newatoms;
	newatoms.reserve(_atoms.size());

	for(const auto& _atom : _atoms)
		newatoms.emplace_back(keep_atom_in_uc<t_vec, t_real>(_atom, dim));

	return newatoms;
}



/**
 * get the supercell vector and the unitcell atom index from a supercell atom position
 */
template<class t_vec, class t_real = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
std::tuple<bool, std::size_t, t_vec>
get_supercell(const t_vec& sc_pos, const t_cont<t_vec>& _uc_sites,
	std::size_t dim = 3,
	t_real eps = std::numeric_limits<t_real>::eps())
requires is_vec<t_vec>
{
	t_cont<t_vec> uc_sites = keep_atoms_in_uc<t_vec, t_real, t_cont>(_uc_sites, dim);

	// unit cell position corresponding to supercell position
	t_vec uc_pos = keep_atom_in_uc<t_vec, t_real>(sc_pos, dim);

	// get corresponding unit cell index
	bool found = false;
	std::size_t uc_site_idx = 0;
	t_vec sc_vec = sc_pos - uc_pos;

	if(auto iter = std::find_if(uc_sites.begin(), uc_sites.end(),
		[&uc_pos, &eps](const t_vec& vec) -> bool
	{
		return equals<t_vec>(uc_pos, vec, eps);
	}); iter != uc_sites.end())
	{
		found = true;
		uc_site_idx = iter - uc_sites.begin();
	}

	return std::make_tuple(found, uc_site_idx, sc_vec);
}



/**
 * transform vectors, e.g. atomic positions, using the given symmetry operations
 */
template<class t_vec, class t_mat, class t_real = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
t_cont<t_vec> apply_ops_hom(const t_vec& _vec, const t_cont<t_mat>& ops,
	t_real eps = std::numeric_limits<t_real>::epsilon(),
	bool keepInUnitCell = true, bool ignore_occupied = false,
	bool ret_hom = false, bool is_pseudovector = false)
requires is_vec<t_vec> && is_mat<t_mat>
{
	// convert vector to homogeneous coordinates
	t_vec vec = _vec;
	if(vec.size() == 3)
		vec = create<t_vec>({ vec[0], vec[1], vec[2], 1. });

	t_cont<t_vec> newvecs;
	newvecs.reserve(ops.size());

	for(const auto& op : ops)
	{
		// apply symmetry operator
		auto newvec = op * vec;

		if(is_pseudovector)
		{
			t_real oldw = newvec[3];
			newvec *= det<t_mat>(submat<t_mat>(op, 3, 3));
			newvec[3] = oldw;
		}

		// return a homogeneous 4-vector or a 3-vector
		if(!ret_hom)
			newvec = create<t_vec>({ newvec[0], newvec[1], newvec[2] });

		if(keepInUnitCell)
			newvec = keep_atom_in_uc<t_vec>(newvec, 3);

		// position already occupied?
		if(ignore_occupied ||
			std::find_if(newvecs.begin(), newvecs.end(),
				[&newvec, eps](const t_vec& vec) -> bool
		{
			return tl2::equals<t_vec>(vec, newvec, eps);
		}) == newvecs.end())
		{
			newvecs.emplace_back(std::move(newvec));
		}
	}

	return newvecs;
}



/**
 * transform matrices using the given symmetry operations
 */
template<class t_mat, class t_real = typename t_mat::value_type,
template<class...> class t_cont = std::vector>
t_cont<t_mat> apply_ops_hom(const t_mat& _mat, const t_cont<t_mat>& ops,
	bool ret_hom = false)
requires is_mat<t_mat>
{
	// convert matrix to homogeneous coordinates
	t_mat mat = _mat;
	if(mat.size1() == 3)
	{
		mat = create<t_mat>({
			mat(0, 0), mat(0, 1), mat(0, 2),  0,
			mat(1, 0), mat(1, 1), mat(1, 2),  0,
			mat(2, 0), mat(2, 1), mat(2, 2),  0,
			        0,         0,         0,  1
		});
	}

	t_cont<t_mat> newmatrices;
	newmatrices.reserve(ops.size());

	for(const auto& op : ops)
	{
		// apply symmetry operator
		auto newmatrix = trans(op) * mat * op;

		// return a homogeneous 4-vector or a 3-vector
		if(!ret_hom)
		{
			newmatrix = create<t_mat>({
				newmatrix(0, 0), newmatrix(0, 1), newmatrix(0, 2),
				newmatrix(1, 0), newmatrix(1, 1), newmatrix(1, 2),
				newmatrix(2, 0), newmatrix(2, 1), newmatrix(2, 2)
			});
		}

		newmatrices.emplace_back(std::move(newmatrix));
	}

	return newmatrices;
}

// ----------------------------------------------------------------------------

}

#endif
