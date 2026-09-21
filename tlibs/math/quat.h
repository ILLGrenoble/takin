/**
 * basic quaternion helpers
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 2013-2016
 * @license GPLv2 or GPLv3
 * @see (Kuipers 2002) for a very good reference on quaternions
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

#ifndef __TLIBS_QUAT_H__
#define __TLIBS_QUAT_H__

#include <boost/math/quaternion.hpp>
#include <limits>
#include "linalg.h"
#include "linalg_ops.h"
#include "rand.h"
#include "../phys/spin.h"

namespace tl {

namespace math = boost::math;


// ------------------------------------------------------------------------------------------------
// ops

template<class t_quat = math::quaternion<double>>
t_quat unit_quat()
{
	return t_quat(1, 0,0,0);
}


/**
 * calculates the quaternion inverse
 * @see e.g.: (Bronstein 2008), Ch. 4
 */
template<class t_quat = math::quaternion<double>>
t_quat quat_inverse(const t_quat& q)
{
	t_quat qc = math::conj(q);
	return qc / (q*qc);
}


/**
 * quaternion product
 * @see: (Kuipers 2002), p. 110
 */
template<class t_quat = math::quaternion<double>>
t_quat quat_prod(const t_quat& q1, const t_quat& q2)
{
	using T = typename t_quat::value_type;
	using t_vec = ublas::vector<T>;

	T r1 = q1.R_component_1();
	T r2 = q2.R_component_1();

	t_vec vec1 = make_vec<t_vec>(
		{q1.R_component_2(), q1.R_component_3(), q1.R_component_4()});
	t_vec vec2 = make_vec<t_vec>(
		{q2.R_component_2(), q2.R_component_3(), q2.R_component_4()});

	T r = r1*r2 - mult<t_vec, t_vec>(vec1, vec2);
	t_vec vec = r1*vec2 + r2*vec1 + cross_3(vec1, vec2);

	return t_quat(r, vec[0], vec[1], vec[2]);
}


// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
// SO(3)

/**
 * 3x3 matrix -> quat
 * @see algo from: http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q55
 */
template<class mat_type = ublas::matrix<double>,
	class quat_type = math::quaternion<typename mat_type::value_type>>
quat_type rot3_to_quat(const mat_type& rot)
{
	using T = typename quat_type::value_type;
	const T tr = trace(rot);
	T v[3], w;

	if(tr > T(0))								// scalar component is largest
	{
		w = T(0.5) * std::sqrt(tr+T(1));
		v[0] = (rot(2,1) - rot(1,2)) / (T(4)*w);
		v[1] = (rot(0,2) - rot(2,0)) / (T(4)*w);
		v[2] = (rot(1,0) - rot(0,1)) / (T(4)*w);
	}
	else
	{
		for(std::size_t iComp=0; iComp<3; ++iComp)	// find largest vector component
		{
			const std::size_t iM = iComp;			// major comp.
			const std::size_t im1 = (iComp+1)%3;	// minor comp. 1
			const std::size_t im2 = (iComp+2)%3;	// minor comp. 2

			if(rot(iM,iM)>=rot(im1,im1) && rot(iM,iM)>=rot(im2,im2))
			{
				v[iM] = T(0.5) * std::sqrt(T(1) + rot(iM,iM) - rot(im1,im1) - rot(im2,im2));
				v[im1] = (rot(im1, iM) + rot(iM, im1)) / (v[iM]*T(4));
				v[im2] = (rot(iM, im2) + rot(im2, iM)) / (v[iM]*T(4));
				w = (rot(im2,im1) - rot(im1,im2)) / (v[iM]*T(4));

				break;
			}

			if(iComp>=2) throw Err("rot3_to_quat: Invalid condition.");
		}
	}

	quat_type quatRet(w, v[0],v[1],v[2]);
	T norm_eucl = math::abs(quatRet);
	return quatRet / norm_eucl;
}


/**
 * quat -> 3x3 matrix
 * @see (Desktop Bronstein 2008), formulas (4.162a/b)
 * @see (Bronstein 2008), p. 296, formulas (4.109a/b)
 */
template<class mat_type=ublas::matrix<double>,
	class quat_type=math::quaternion<typename mat_type::value_type>>
mat_type quat_to_rot3(const quat_type& quat)
{
	const quat_type cquat = math::conj(quat);
	const quat_type i(0,1,0,0), j(0,0,1,0), k(0,0,0,1);

	const quat_type cols[] =
	{
		quat * i * cquat,
		quat * j * cquat,
		quat * k * cquat
	};

	mat_type mat(3,3);
	for(std::size_t icol=0; icol<3; ++icol)
	{
		mat(0, icol) = cols[icol].R_component_2();
		mat(1, icol) = cols[icol].R_component_3();
		mat(2, icol) = cols[icol].R_component_4();
	}

	return mat;
}

// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
// vector ops

/**
 * vector -> quat
 * @see (Kuipers 2002), p. 114
 */
template<class t_vec = ublas::vector<double>,
	class t_quat = math::quaternion<typename t_vec::value_type>>
t_quat vec3_to_quat(const t_vec& vec)
{
	using T = typename t_vec::value_type;
	return t_quat(T(0), vec[0], vec[1], vec[2]);
}


/**
 * quat, vector product
 * @see (Kuipers 2002), p. 127
 */
template<class t_vec = ublas::vector<double>,
	class t_quat = math::quaternion<typename t_vec::value_type>>
t_vec quat_vec_prod(const t_quat& q, const t_vec& v)
{
	t_quat qv = vec3_to_quat<t_vec, t_quat>(v);
	t_quat qvq =  q * qv * math::conj(q);

	t_vec vec(3);
	vec[0] = qvq.R_component_2();
	vec[1] = qvq.R_component_3();
	vec[2] = qvq.R_component_4();
	return vec;
}

// ------------------------------------------------------------------------------------------------


/**
 * quat -> complex 2x2 matrix
 * @see e.g. (Scherer 2010), p.173
 */
template<template<class...> class t_mat = ublas::matrix,
	class t_real = double,
	class t_quat = math::quaternion<t_real>>
t_mat<std::complex<t_real>> quat_to_cmat(const t_quat& quat)
{
	const auto vecS = get_spin_matrices<t_mat, ublas::vector, t_real>();
	const auto matI = unit_m<t_mat<std::complex<t_real>>>(2);

	t_mat<std::complex<t_real>> mat =
		std::complex<t_real>(quat.R_component_1()) * matI +
		std::complex<t_real>(quat.R_component_2()) * vecS[0] +
		std::complex<t_real>(quat.R_component_3()) * vecS[1] +
		std::complex<t_real>(quat.R_component_4()) * vecS[2];
	return mat;
}


// ------------------------------------------------------------------------------------------------
// rotation axis

template<class quat_type = math::quaternion<double>,
	typename T = typename quat_type::value_type>
std::vector<T> quat_to_euler(const quat_type& quat);


template<typename T=double, class... Args>
std::vector<T> rotation_angle(const ublas::matrix<T, Args...>& rot)
{
	std::vector<T> vecResult;

	if(rot.size1()!=rot.size2())
		return vecResult;
	if(rot.size1()<2)
		return vecResult;

	if(rot.size2()==2)
	{
		// rot = ( c -s )
		//       ( s  c )
		T angle = std::atan2(rot(1,0), rot(0,0));
		vecResult.push_back(angle);
	}
	else if(rot.size2()==3)
	{
		math::quaternion<T> quat =
			rot3_to_quat<decltype(rot), math::quaternion<T>>(rot);
		vecResult = quat_to_euler(quat);
	}

	return vecResult;
}


/**
 * rotation angle
 * @see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
 */
template<typename T=double>
T rotation_angle(const math::quaternion<T>& quat)
{
	//return 2.*std::asin(math::abs(math::unreal(quat)));
	return T(2)*std::acos(quat.R_component_1());
}


/**
 * quat -> rotation axis
 * @see (Bronstein 2008), Ch. 4, pp. 301-302
 */
template<class t_vec = ublas::vector<double>>
t_vec rotation_axis(const math::quaternion<typename t_vec::value_type>& quat)
{
	using T = typename t_vec::value_type;

	t_vec vec(3);
	vec[0] = quat.R_component_2();
	vec[1] = quat.R_component_3();
	vec[2] = quat.R_component_4();

	typename t_vec::value_type angle = rotation_angle(quat);
	vec /= std::sin(T(0.5)*angle);

	return vec;
}


/**
 * rotation axis -> quat
 * @see (Desktop Bronstein 2008), formula (4.193)
 * @see (Bronstein 2008), p. 302
 */
template<class quat_type = math::quaternion<double>,
	class vec_type = ublas::vector<typename quat_type::value_type>,
	typename T = typename quat_type::value_type>
quat_type rotation_quat(const vec_type& vec, const T angle)
{
	const T s = std::sin(T(0.5)*angle);
	const T c = std::cos(T(0.5)*angle);
	const T n = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

	const T x = s * vec[0] / n;
	const T y = s * vec[1] / n;
	const T z = s * vec[2] / n;
	const T r = c;

	return quat_type(r, x,y,z);
}


/**
 * quaternion to rotate vec0 into vec1
 */
template<class t_quat = math::quaternion<double>,
	class t_vec = ublas::vector<typename t_quat::value_type>,
	typename T = typename t_quat::value_type>
t_quat rotation_quat(const t_vec& _vec0, const t_vec& _vec1)
{
	t_vec vec0 = _vec0 / veclen(_vec0);
	t_vec vec1 = _vec1 / veclen(_vec1);

	if(vec_equal(vec0, vec1))
	{ // parallel vectors -> do nothing
		return unit_quat<t_quat>();
	}
	else if(vec_equal(vec0, t_vec(-vec1)))
	{ // antiparallel vectors -> rotate about any perpendicular axis
		t_vec vecRnd(3), vecPerp(3);
		T lenPerp2 = 0.;
		const T eps = std::numeric_limits<T>::epsilon();

		while(lenPerp2 <= eps)
		{
			vecRnd[0] = 1. - 2.*tl::rand01<T>(); // vec0[2];
			vecRnd[1] = 1. - 2.*tl::rand01<T>(); // 0;
			vecRnd[2] = 1. - 2.*tl::rand01<T>(); // -vec0[0];

			vecPerp = cross_3<t_vec>(vec0, vecRnd);
			lenPerp2 = vecPerp[0]*vecPerp[0] + vecPerp[1]*vecPerp[1] + vecPerp[2]*vecPerp[2];
		}

		return rotation_quat<t_quat, t_vec, T>(vecPerp, get_pi<T>());
	}

	t_vec veccross = cross_3<t_vec>(vec0, vec1);

	T dC = mult<t_vec, t_vec>(vec0, vec1);
	T dS = veclen(veccross);
	T dAngle = std::atan2(dS, dC);

	return rotation_quat<t_quat, t_vec, T>(veccross, dAngle);
}


template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type rotation_quat_x(typename quat_type::value_type angle)
{
	return quat_type(std::cos(T(0.5)*angle),
		std::sin(T(0.5)*angle), T(0), T(0));
}


template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type rotation_quat_y(typename quat_type::value_type angle)
{
	return quat_type(std::cos(T(0.5)*angle),
		T(0), std::sin(T(0.5)*angle), T(0));
}


template<class quat_type=math::quaternion<double>, typename T = typename quat_type::value_type>
quat_type rotation_quat_z(typename quat_type::value_type angle)
{
	return quat_type(std::cos(T(0.5)*angle),
		T(0), T(0), std::sin(T(0.5)*angle));
}

// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
// Euler angles

/**
 * XYZ euler angles -> quat
 * @see (Kuipers 2002), pp. 166, 167
 */
template<class t_quat = math::quaternion<double>,
	typename T = typename t_quat::value_type>
t_quat euler_to_quat_xyz(T phi, T theta, T psi)
{
	t_quat q1 = rotation_quat_x<t_quat>(phi);
	t_quat q2 = rotation_quat_y<t_quat>(theta);
	t_quat q3 = rotation_quat_z<t_quat>(psi);

	return q3 * q2 * q1;
}


/**
 * ZXZ euler angles -> quat
 * @see (Kuipers 2002), pp. 166, 167
 */
template<class t_quat = math::quaternion<double>,
	typename T = typename t_quat::value_type>
t_quat euler_to_quat_zxz(T phi, T theta, T psi)
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
template<class quat_type, typename T>
std::vector<T> quat_to_euler_xyz(const quat_type& quat)
{
	T q[] = { quat.R_component_1(), quat.R_component_2(),
		quat.R_component_3(), quat.R_component_4() };

	// formulas from:
	// http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
	T phi = std::atan2(T(2)*(q[0]*q[1] + q[2]*q[3]), T(1)-T(2)*(q[1]*q[1] + q[2]*q[2]));
	T theta = std::asin(T(2)*(q[0]*q[2] - q[3]*q[1]));
	T psi = std::atan2(T(2)*(q[0]*q[3] + q[1]*q[2]), T(1)-T(2)*(q[2]*q[2] + q[3]*q[3]));

	return std::vector<T>({ phi, theta, psi });
}


// use XYZ version by default
template<class t_quat = math::quaternion<double>,
	typename T = typename t_quat::value_type>
t_quat euler_to_quat(T phi, T theta, T psi)
{
	return euler_to_quat_xyz<t_quat, T>(phi, theta, psi);
}


template<class quat_type, typename T>
std::vector<T> quat_to_euler(const quat_type& quat)
{
	return quat_to_euler_xyz<quat_type, T>(quat);
}

// ------------------------------------------------------------------------------------------------



/**
 * @see (Desktop Bronstein 2008), formula (4.217)
 * @see (Bronstein 2008), p. 308, formula (4.165)
 */
template<class quat_type=math::quaternion<double>,
	typename T = typename quat_type::value_type>
quat_type stereo_proj(const quat_type& quat)
{
	return (T(1)+quat) / (T(1)-quat);
}


/**
 * @see (Desktop Bronstein 2008), formula (4.217)
 * @see (Bronstein 2008), p. 308, formula (4.165)
 */
template<class quat_type=math::quaternion<double>,
	typename T = typename quat_type::value_type>
quat_type stereo_proj_inv(const quat_type& quat)
{
	return (T(1)-quat) / (T(1)+quat);
}


template<class QUAT>
struct vec_angle_unsigned_impl<QUAT, LinalgType::QUATERNION>
{
	typename QUAT::value_type operator()(const QUAT& q1, const QUAT& q2) const
	{
		typedef typename QUAT::value_type REAL;
		REAL dot = q1.R_component_1() * q2.R_component_1() +
			q1.R_component_2() * q2.R_component_2() +
			q1.R_component_3() * q2.R_component_3() +
			q1.R_component_4() * q2.R_component_4();

		dot /= math::abs(q1);
		dot /= math::abs(q2);

		return std::acos(dot);
	}
};

}
#endif
