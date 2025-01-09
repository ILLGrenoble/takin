/**
 * tlibs2 maths library -- quaternion algorithms
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

#ifndef __TLIBS2_MATHS_QUAT_H__
#define __TLIBS2_MATHS_QUAT_H__

#include <cmath>

#include "decls.h"
#include "constants.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// quaternion algorithms
// @see (Kuipers 2002) for infos
// ----------------------------------------------------------------------------

/**
 * unsigned angle between two quaternions
 * <q1|q2> / (|q1| |q2|) = cos(alpha)
 */
template<class t_quat>
typename t_quat::value_type angle_unsigned(const t_quat& q1, const t_quat& q2)
requires is_quat<t_quat>
{
	using t_real = typename t_quat::value_type;

	t_real dot = q1.R_component_1() * q2.R_component_1() +
		q1.R_component_2() * q2.R_component_2() +
		q1.R_component_3() * q2.R_component_3() +
		q1.R_component_4() * q2.R_component_4();

	t_real len1 = q1.R_component_1() * q1.R_component_1() +
		q1.R_component_2() * q1.R_component_2() +
		q1.R_component_3() * q1.R_component_3() +
		q1.R_component_4() * q1.R_component_4();

	t_real len2 = q2.R_component_1() * q2.R_component_1() +
		q2.R_component_2() * q2.R_component_2() +
		q2.R_component_3() * q2.R_component_3() +
		q2.R_component_4() * q2.R_component_4();

	len1 = std::sqrt(len1);
	len2 = std::sqrt(len2);

	dot /= len1;
	dot /= len2;

	return std::acos(dot);
}


/**
 * slerp
 * @see K. Shoemake, "Animating rotation with quaternion curves", http://dx.doi.org/10.1145/325334.325242
 * @see (Desktop Bronstein 2008), formula 4.207
 * @see (Bronstein 2008), p. 306, formula 4.155
 */
template<class T>
T slerp(const T& q1, const T& q2, typename T::value_type t)
{
	using t_real = typename T::value_type;
	t_real angle = angle_unsigned<T>(q1, q2);

	T q = std::sin((t_real(1)-t)*angle)/std::sin(angle) * q1 +
	std::sin(t*angle)/std::sin(angle) * q2;

	return q;
}


/**
 * set values close to an integer value to that integer
 * quaternion version
 */
template<typename t_quat, typename t_real = typename t_quat::value_type>
void set_eps_round(t_quat& quat, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_quat<t_quat>
{
	t_real re = quat.R_component_1();
	t_real im1 = quat.R_component_2();
	t_real im2 = quat.R_component_3();
	t_real im3 = quat.R_component_4();

	set_eps_round<t_real>(re, eps);
	set_eps_round<t_real>(im1, eps);
	set_eps_round<t_real>(im2, eps);
	set_eps_round<t_real>(im3, eps);

	quat = t_quat(re, im1, im2, im3);
};


/**
 * set values lower than epsilon to zero
 * quaternion version
 */
template<typename t_quat, typename t_real = typename t_quat::value_type>
void set_eps_0(t_quat& quat, t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_quat<t_quat>
{
	t_real re = quat.R_component_1();
	t_real im1 = quat.R_component_2();
	t_real im2 = quat.R_component_3();
	t_real im3 = quat.R_component_4();

	set_eps_0<t_real>(re, eps);
	set_eps_0<t_real>(im1, eps);
	set_eps_0<t_real>(im2, eps);
	set_eps_0<t_real>(im3, eps);

	quat = t_quat(re, im1, im2, im3);
};


template<class t_quat> t_quat unit_quat() requires is_quat<t_quat>
{
	return t_quat(1, 0,0,0);
}


/**
 * calculates the quaternion inverse
 * @see (Bronstein 2008), Ch. 4
 * @see (Kuipers 2002), p. 112
 * @see https://en.wikipedia.org/wiki/Quaternion#Conjugation,_the_norm,_and_reciprocal
 */
template<class t_quat> t_quat inv(const t_quat& q) requires is_quat<t_quat>
{
	t_quat qc{q.R_component_1(), -q.R_component_2(), -q.R_component_3(), -q.R_component_4()};
	return qc / (q*qc);
}


/**
 * quaternion product
 * @see (Kuipers 2002), pp. 106-110
 * @see https://en.wikipedia.org/wiki/Quaternion#Scalar_and_vector_parts
 */
template<class t_quat, class t_vec>
t_quat prod(const t_quat& q1, const t_quat& q2)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_quat::value_type;

	T r1 = q1.R_component_1();
	T r2 = q2.R_component_1();

	t_vec vec1 = create<t_vec>({q1.R_component_2(), q1.R_component_3(), q1.R_component_4()});
	t_vec vec2 = create<t_vec>({q2.R_component_2(), q2.R_component_3(), q2.R_component_4()});

	T r = r1*r2 - inner<t_vec>(vec1, vec2);
	t_vec vec = r1*vec2 + r2*vec1 + cross<t_vec>({vec1, vec2});;

	return t_quat(r, vec[0], vec[1], vec[2]);
}


/**
 * 3x3 matrix -> quat
 * @desc algo from: http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q55
 */
template<class t_mat, class t_quat>
t_quat rot3_to_quat(const t_mat& rot)
requires is_quat<t_quat> && is_mat<t_mat>
{
	using T = typename t_quat::value_type;
	const T tr = trace<t_mat>(rot);
	T v[3]{}, w{};

	if(tr > T(0))	// scalar component is largest
	{
		w = T(0.5) * std::sqrt(tr+T(1));
		v[0] = (rot(2,1) - rot(1,2)) / (T(4)*w);
		v[1] = (rot(0,2) - rot(2,0)) / (T(4)*w);
		v[2] = (rot(1,0) - rot(0,1)) / (T(4)*w);
	}
	else
	{
		for(std::size_t iComp = 0; iComp < 3; ++iComp)	// find largest vector component
		{
			const std::size_t iM = iComp;		// major comp.
			const std::size_t im1 = (iComp+1)%3;	// minor comp. 1
			const std::size_t im2 = (iComp+2)%3;	// minor comp. 2

			if(rot(iM,iM) >= rot(im1,im1) && rot(iM,iM) >= rot(im2,im2))
			{
				v[iM] = T(0.5) * std::sqrt(T(1) + rot(iM,iM) - rot(im1,im1) - rot(im2,im2));
				v[im1] = (rot(im1, iM) + rot(iM, im1)) / (v[iM]*T(4));
				v[im2] = (rot(iM, im2) + rot(im2, iM)) / (v[iM]*T(4));
				w = (rot(im2,im1) - rot(im1,im2)) / (v[iM]*T(4));

				break;
			}

			if(iComp >= 2)
				assert(false);  // should not get here
		}
	}

	t_quat quatRet{w, v[0], v[1], v[2]};
	T norm_eucl = std::sqrt(w*w + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	return quatRet / norm_eucl;
}


/**
 * quat -> 3x3 matrix
 * @see (Desktop Bronstein 2008), formulas (4.162a/b)
 * @see (Bronstein 2008), p. 296, formulas (4.109a/b)
 */
template<class t_quat, class t_mat>
t_mat quat_to_rot3(const t_quat& quat)
requires is_quat<t_quat> && is_mat<t_mat>
{
	t_quat qc{quat.R_component_1(), -quat.R_component_2(), -quat.R_component_3(), -quat.R_component_4()};
	const t_quat i{0,1,0,0}, j{0,0,1,0}, k{0,0,0,1};
	const t_quat cols[] = { quat*i*qc, quat*j*qc, quat*k*qc };

	t_mat mat = unit<t_mat>(3);
	for(std::size_t icol = 0; icol < 3; ++icol)
	{
		mat(0, icol) = cols[icol].R_component_2();
		mat(1, icol) = cols[icol].R_component_3();
		mat(2, icol) = cols[icol].R_component_4();
	}

	return mat;
}


/**
 * vector -> quat
 * @see (Kuipers 2002), p. 114
 */
template<class t_vec, class t_quat>
t_quat vec3_to_quat(const t_vec& vec)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	return t_quat{T(0), vec[0], vec[1], vec[2]};
}


/**
 * quat, vector product
 * @see (Kuipers 2002), p. 127
 */
template<class t_quat, class t_vec>
t_vec quat_vec_prod(const t_quat& q, const t_vec& v)
requires is_quat<t_quat> && is_vec<t_vec>
{
	t_quat qv = vec3_to_quat<t_vec, t_quat>(v);
	t_quat qc{q.R_component_1(), -q.R_component_2(), -q.R_component_3(), -q.R_component_4()};
	t_quat qvq =  q * qv * qc;

	return create<t_vec>({ qvq.R_component_2(), qvq.R_component_3(), qvq.R_component_4() });
}


/**
 * quat -> complex 2x2 matrix
 * @see (Scherer 2010), p.173
 * @see (Desktop Bronstein 2008), ch. 4, equations (4.163a) and (4.163b)
 * @see (Bronstein 2008), ch. 4, p. 296, equation (4.110a) and (4.110b)
 */
template<class t_mat, class t_mat_cplx, class t_quat>
t_mat_cplx quat_to_cmat(const t_quat& quat)
requires is_quat<t_quat> && is_mat<t_mat> && is_mat<t_mat_cplx>
{
	using t_cplx = typename t_mat_cplx::value_type;

	const auto matI = unit<t_mat_cplx>(2);
	t_mat_cplx mat =
		t_cplx(quat.R_component_1()) * matI +
		t_cplx(quat.R_component_2()) * su2_matrix<t_mat_cplx>(0) +
		t_cplx(quat.R_component_3()) * su2_matrix<t_mat_cplx>(1) +
		t_cplx(quat.R_component_4()) * su2_matrix<t_mat_cplx>(2);

	return mat;
}


/**
 * rotation angle
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
 */
template<class t_quat>
typename t_quat::value_type rotation_angle(const t_quat& quat)
requires is_quat<t_quat>
{
	using t_real = typename t_quat::value_type;
	return t_real{2}*std::acos(quat.R_component_1());
}


/**
 * quat -> rotation axis
 * @see (Bronstein 2008), Ch. 4, pp. 301-302
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 */
template<class t_quat, class t_vec>
std::pair<t_vec, typename t_vec::value_type> rotation_axis(const t_quat& quat)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec = create<t_vec>({
		quat.R_component_2(),
		quat.R_component_3(),
		quat.R_component_4()
	});

	T angle = rotation_angle(quat);
	vec /= std::sin(T{0.5}*angle);

	return std::make_pair(vec, angle);
}


/**
 * rotation axis -> quat
 * @see (Desktop Bronstein 2008), formula (4.193)
 * @see (Bronstein 2008), ch. 4, pp. 301-302
 * @see https://en.wikipedia.org/wiki/Quaternion#Exponential,_logarithm,_and_power_functions
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 */
template<class t_vec, class t_quat>
t_quat rotation_quat(const t_vec& vec, typename t_vec::value_type angle)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	const T s = std::sin(T(0.5)*angle);
	const T c = std::cos(T(0.5)*angle);
	const T n = norm<t_vec>(vec);

	const T x = s * vec[0] / n;
	const T y = s * vec[1] / n;
	const T z = s * vec[2] / n;
	const T r = c;

	return t_quat{r, x,y,z};
}


/**
 * quaternion to rotate vec0 into vec1
 */
template<class t_vec, class t_quat>
t_quat rotation_quat(const t_vec& _vec0, const t_vec& _vec1)
requires is_quat<t_quat> && is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec0 = _vec0 / norm<t_vec>(_vec0);
	t_vec vec1 = _vec1 / norm<t_vec>(_vec1);

	// parallel vectors -> do nothing
	if(equals<t_vec>(vec0, vec1))
	{
		return unit_quat<t_quat>();
	}

	// antiparallel vectors -> rotate about any perpendicular axis
	else if(equals<t_vec>(vec0, -vec1))
	{
		t_vec vecPerp = perp<t_vec>(vec0);
		return rotation_quat<t_vec, t_quat>(vecPerp, pi<T>);
	}

	// rotation axis from cross product
	t_vec vecaxis = cross<t_vec>({vec0, vec1});

	T dC = inner<t_vec>(vec0, vec1);
	T dS = norm<t_vec>(vecaxis);

	// rotation angle
	T dAngle = std::atan2(dS, dC);

	return rotation_quat<t_vec, t_quat>(vecaxis, dAngle);
}


template<class t_quat>
t_quat rotation_quat_x(typename t_quat::value_type angle)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return t_quat{std::cos(T(0.5)*angle), std::sin(T(0.5)*angle), T(0), T(0)};
}


template<class t_quat>
t_quat rotation_quat_y(typename t_quat::value_type angle)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return t_quat{std::cos(T(0.5)*angle), T(0), std::sin(T(0.5)*angle), T(0)};
}


template<class t_quat>
t_quat rotation_quat_z(typename t_quat::value_type angle)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return t_quat{std::cos(T(0.5)*angle), T(0), T(0), std::sin(T(0.5)*angle)};
}



/**
 * XYZ euler angles -> quat
 * @see (Kuipers 2002), pp. 166-167
 */
template<class t_quat>
t_quat euler_to_quat_xyz(
	typename t_quat::value_type phi, typename t_quat::value_type theta, typename t_quat::value_type psi)
requires is_quat<t_quat>
{
	t_quat q1 = rotation_quat_x<t_quat>(phi);
	t_quat q2 = rotation_quat_y<t_quat>(theta);
	t_quat q3 = rotation_quat_z<t_quat>(psi);

	return q3 * q2 * q1;
}


/**
 * ZXZ euler angles -> quat
 * @see (Kuipers 2002), pp. 166-167
 */
template<class t_quat>
t_quat euler_to_quat_zxz(
	typename t_quat::value_type phi, typename t_quat::value_type theta, typename t_quat::value_type psi)
requires is_quat<t_quat>
{
	t_quat q1 = rotation_quat_z<t_quat>(phi);
	t_quat q2 = rotation_quat_x<t_quat>(theta);
	t_quat q3 = rotation_quat_z<t_quat>(psi);

	return q3 * q2 * q1;
}


/**
 * quat -> XYZ euler angles
 * @see http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
 */
template<class t_quat>
std::vector<typename t_quat::value_type>
quat_to_euler_xyz(const t_quat& quat)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	T q[] = { quat.R_component_1(), quat.R_component_2(),
		quat.R_component_3(), quat.R_component_4() };

	// formulas from:
	// http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	T phi = std::atan2(T(2)*(q[0]*q[1] + q[2]*q[3]), T(1)-T(2)*(q[1]*q[1] + q[2]*q[2]));
	T theta = std::asin(T(2)*(q[0]*q[2] - q[3]*q[1]));
	T psi = std::atan2(T(2)*(q[0]*q[3] + q[1]*q[2]), T(1)-T(2)*(q[2]*q[2] + q[3]*q[3]));

	return std::vector<T>({ phi, theta, psi });
}


/**
 * @see (Desktop Bronstein 2008), formula (4.217)
 * @see (Bronstein 2008), p. 308, formula (4.165)
 */
template<class t_quat>
t_quat stereo_proj(const t_quat& quat)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return (T{1}+quat) / (T{1}-quat);
}


/**
 * @see (Desktop Bronstein 2008), formula (4.217)
 * @see (Bronstein 2008), p. 308, formula (4.165)
 */
template<class t_quat>
t_quat stereo_proj_inv(const t_quat& quat)
requires is_quat<t_quat>
{
	using T = typename t_quat::value_type;
	return (T{1}-quat) / (T{1}+quat);
}

// ----------------------------------------------------------------------------


}

#endif
