/**
 * Bravais Lattice Calculations
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014-2016
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

#ifndef __TLIBS_LATTICE_H__
#define __TLIBS_LATTICE_H__


#include "../math/linalg.h"
#include "../math/quat.h"
#include "../math/math.h"
#include "../math/tensor.h"
#include "../math/geo.h"

#include "neutrons.h"

#include <ostream>


namespace tl {


/**
 * reciprocal matrix to A: B = 2*pi * A^(-T)
 * @see e.g.: https://en.wikipedia.org/wiki/Reciprocal_lattice
 */
template<typename T = double>
bool reciprocal(const ublas::matrix<T>& matReal, ublas::matrix<T>& matRecip)
{
	ublas::matrix<T> matInv;
	if(!inverse<ublas::matrix<T>>(transpose(matReal), matInv))
		return false;

	matRecip = T(2)*get_pi<T>()*matInv;
	return true;
}


template<typename T = double>
class Lattice
{
	public:
		using t_vec = ublas::vector<T>;
		using t_mat = ublas::matrix<T>;

	protected:
		t_vec m_vecs[3];

	public:
		Lattice(T a, T b, T c, T alpha, T beta, T gamma);
		Lattice(const t_vec& vec0, const t_vec& vec1, const t_vec& vec2);
		Lattice(const Lattice<T>& lattice);
		const Lattice<T>& operator=(const Lattice<T>& lattice);

		Lattice() = default;
		~Lattice() = default;

		bool IsInited() const { return m_vecs[0].size()!=0; }

		/**
		 * Euler ZXZ rotation
		 */
		void RotateEuler(T dPhi, T dTheta, T dPsi);

		/**
		 * Euler vecRecipZ vecRecipX vecRecipZ rotation
		 */
		void RotateEulerRecip(const t_vec& vecRecipX,
			const t_vec& vecRecipY, const t_vec& vecRecipZ,
			T dPhi, T dTheta, T dPsi);

		Lattice GetRecip() const;
		Lattice GetAligned() const;

		t_vec GetPos(T h, T k, T l) const;
		t_vec GetHKL(const t_vec& vec) const;

		T GetAlpha() const;
		T GetBeta() const;
		T GetGamma() const;

		T GetA() const;
		T GetB() const;
		T GetC() const;

		T GetVol() const;

		const t_vec& GetVec(std::size_t i) const { return m_vecs[i]; }

		/**
		 * covariant / contravariant base matrix
		 */
		t_mat GetBaseMatrixCov() const;
		t_mat GetBaseMatrixCont() const;

		/**
		 * covariant / contravariant metric
		 */
		t_mat GetMetricCov() const;
		t_mat GetMetricCont() const;

		bool IsCubic() const;
};


template<typename T>
Lattice<T>::Lattice(T a, T b, T c, T alpha, T beta, T gamma)
{
	m_vecs[0].resize(3,0);
	m_vecs[1].resize(3,0);
	m_vecs[2].resize(3,0);

	fractional_basis_from_angles(a, b, c,
		alpha, beta, gamma,
		m_vecs[0], m_vecs[1], m_vecs[2]);
}


template<typename T>
Lattice<T>::Lattice(const t_vec& vec0, const t_vec& vec1, const t_vec& vec2)
{
	this->m_vecs[0] = vec0;
	this->m_vecs[1] = vec1;
	this->m_vecs[2] = vec2;
}


template<typename T>
Lattice<T>::Lattice(const Lattice<T>& lattice)
{
	this->operator=(lattice);
}


template<typename T>
const Lattice<T>& Lattice<T>::operator=(const Lattice<T>& lattice)
{
	this->m_vecs[0] = lattice.m_vecs[0];
	this->m_vecs[1] = lattice.m_vecs[1];
	this->m_vecs[2] = lattice.m_vecs[2];
	return *this;
}


template<typename T>
void Lattice<T>::RotateEuler(T dPhi, T dTheta, T dPsi)
{
	t_mat mat1 = rotation_matrix_3d_z(dPhi);
	t_mat mat2 = rotation_matrix_3d_x(dTheta);
	t_mat mat3 = rotation_matrix_3d_z(dPsi);

	t_mat mat21 = prod_mm(mat2,mat1);
	t_mat mat = prod_mm(mat3, mat21);

	for(std::size_t i=0; i<3; ++i)
		m_vecs[i] = prod_mv(mat, m_vecs[i]);
}


template<typename T>
void Lattice<T>::RotateEulerRecip(const t_vec& vecRecipX,
	const t_vec& vecRecipY, const t_vec& vecRecipZ,
	T dPhi, T dTheta, T dPsi)
{
	// get real vectors
	const std::size_t iDim=3;
	t_mat matReal = column_matrix({ vecRecipX, vecRecipY, vecRecipZ });
	if(matReal.size1()!=matReal.size2() || matReal.size1()!=iDim)
		throw Err("Invalid real matrix.");

	t_mat matRecip;
	if(!reciprocal(matReal, matRecip))
		throw Err("Reciprocal matrix could not be calculated.");

	t_vec vecX = get_column(matRecip,0);
	t_vec vecY = get_column(matRecip,1);
	t_vec vecZ = get_column(matRecip,2);

	T dLenX = veclen(vecX);
	T dLenY = veclen(vecY);
	T dLenZ = veclen(vecZ);

	if(float_equal<T>(dLenX, 0.) || float_equal<T>(dLenY, 0.) || float_equal<T>(dLenZ, 0.)
		|| std::isnan(dLenX) || std::isnan(dLenY) || std::isnan(dLenZ))
	{
		throw Err("Invalid reciprocal matrix.");
	}

	vecX /= dLenX;
	vecY /= dLenY;
	vecZ /= dLenZ;


	// rotate around real vectors
	t_mat mat1 = rotation_matrix(vecZ, dPhi);
	t_mat mat2 = rotation_matrix(vecX, dTheta);
	t_mat mat3 = rotation_matrix(vecZ, dPsi);

	t_mat mat21 = prod_mm(mat2, mat1);
	t_mat mat = prod_mm(mat3, mat21);

	for(std::size_t i=0; i<3; ++i)
		m_vecs[i] = prod_mv(mat, m_vecs[i]);
}


template<typename T> T Lattice<T>::GetAlpha() const
{ return std::acos(inner(m_vecs[1]/GetB(), m_vecs[2]/GetC())); }
template<typename T> T Lattice<T>::GetBeta() const
{ return std::acos(inner(m_vecs[0]/GetA(), m_vecs[2]/GetC())); }
template<typename T> T Lattice<T>::GetGamma() const
{ return std::acos(inner(m_vecs[0]/GetA(), m_vecs[1]/GetB())); }


template<typename T> T Lattice<T>::GetA() const { return veclen(m_vecs[0]); }
template<typename T> T Lattice<T>::GetB() const { return veclen(m_vecs[1]); }
template<typename T> T Lattice<T>::GetC() const { return veclen(m_vecs[2]); }


template<typename T>
T Lattice<T>::GetVol() const
{
	return get_volume(column_matrix({ m_vecs[0], m_vecs[1], m_vecs[2] }));
}


/**
 (x)   (v0_x v1_x v2_x) (h)
 (y) = (v0_y v1_y v2_y) (k)
 (z)   (v0_z v1_z v2_z) (l)
 */
template<typename T>
typename Lattice<T>::t_vec Lattice<T>::GetPos(T h, T k, T l) const
{
	return h*m_vecs[0] + k*m_vecs[1] + l*m_vecs[2];
}


/**
 (h)   (v0_x v1_x v2_x)^(-1) (x)
 (k) = (v0_y v1_y v2_y)      (y)
 (l)   (v0_z v1_z v2_z)      (z)
 */
template<typename T>
typename Lattice<T>::t_vec Lattice<T>::GetHKL(const t_vec& vec) const
{
	t_mat mat = column_matrix({ m_vecs[0], m_vecs[1], m_vecs[2] });

	t_mat matInv;
	if(!inverse(mat, matInv))
		throw Err("Miller indices could not be calculated.");

	return prod_mv(matInv, vec);
}


template<typename T>
Lattice<T> Lattice<T>::GetRecip() const
{
	const std::size_t iDim = 3;
	t_mat matReal = column_matrix({ m_vecs[0], m_vecs[1], m_vecs[2] });
	if(matReal.size1()!=matReal.size2() || matReal.size1()!=iDim)
		throw Err("Invalid real lattice matrix.");

	t_mat matRecip;
	if(!reciprocal(matReal, matRecip))
		throw Err("Reciprocal lattice could not be calculated.");

	// warning: first axis does not (necessarily) coincide with assumed first
	// orientation vector [0,0,1] anymore!
	return Lattice<T>(get_column(matRecip,0), get_column(matRecip,1),
		get_column(matRecip,2));
}


template<typename T>
Lattice<T> Lattice<T>::GetAligned() const
{
	// construct new, correctly oriented reciprocal lattice with first axis along
	// [0,0,1]
	return Lattice<T>(GetA(), GetB(), GetC(), GetAlpha(), GetBeta(), GetGamma());
}


template<typename T>
typename Lattice<T>::t_mat Lattice<T>::GetBaseMatrixCov() const
{
	t_mat matBase = column_matrix({ m_vecs[0], m_vecs[1], m_vecs[2] });
	set_eps_0(matBase);
	return matBase;
}


template<typename T>
typename Lattice<T>::t_mat Lattice<T>::GetBaseMatrixCont() const
{
	t_mat matCov = GetBaseMatrixCov();
	t_mat matGCont = GetMetricCont() /** T(2)*get_pi<T>()*/;

	t_mat matBase = prod_mm(matCov, matGCont);

	set_eps_0(matBase);
	return matBase;
}


template<typename T>
typename Lattice<T>::t_mat Lattice<T>::GetMetricCov() const
{
	t_mat matG = make_metric_cov({ m_vecs[0], m_vecs[1], m_vecs[2] });
	set_eps_0(matG);
	return matG;
}


template<typename T>
typename Lattice<T>::t_mat Lattice<T>::GetMetricCont() const
{
	t_mat matGCov = GetMetricCov();

	t_mat matGCont;
	inverse(matGCov, matGCont);

	set_eps_0(matGCont);
	return matGCont;
}



template<typename T>
bool Lattice<T>::IsCubic() const
{
	bool bEqualLen = float_equal<T>(GetA(), GetB()) && float_equal<T>(GetB(), GetC());
	if(!bEqualLen) return false;

	static const T pihalf = get_pi<T>()*T(0.5);
	bool b90Deg = float_equal<T>(GetAlpha(), pihalf) &&
		float_equal<T>(GetBeta(), pihalf) &&
		float_equal<T>(GetGamma(), pihalf);

	return b90Deg;
}


// -----------------------------------------------------------------------------


/**
 * B matrix converts rlu to 1/A
 * @see e.g.: https://doi.org/10.1107/S0021889805004875
 */
template<typename T = double>
ublas::matrix<T> get_B(const Lattice<T>& lattice, bool bIsRealLattice = true)
{
	using t_mat = ublas::matrix<T>;

	t_mat matB;
	if(bIsRealLattice)
		matB = lattice.GetRecip().GetBaseMatrixCov();
	else
		matB = lattice.GetBaseMatrixCov();

	return matB;
}


/**
 * U matrix expresses the coordinates in the basis of the scattering plane
 * @see e.g.: https://doi.org/10.1107/S0021889805004875
 */
template<typename T = double>
ublas::matrix<T> get_U(const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2,
	const ublas::matrix<T>* pmatB = nullptr)
{
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

	t_vec vec1, vec2;

	if(pmatB)
	{
		// in 1/A
		vec1 = prod_mv(*pmatB, _vec1);
		vec2 = prod_mv(*pmatB, _vec2);
	}
	else
	{
		// in rlu
		vec1 = _vec1;
		vec2 = _vec2;
	}

	// U: scattering plane coordinate system
	t_mat matU = row_matrix(get_ortho_rhs({ vec1, vec2 }));
	return matU;
}


/**
 * UB matrix converts rlu to 1/A and expresses it in the scattering plane coords:
 * Q = U*B*hkl
 * @see e.g.: https://doi.org/10.1107/S0021889805004875
 */
template<typename T = double>
ublas::matrix<T> get_UB(const Lattice<T>& lattice_real,
	const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2)
{
	using t_mat = ublas::matrix<T>;

	t_mat matB = get_B(lattice_real, true);		// rlu to 1/A
	t_mat matU = get_U(_vec1, _vec2, &matB);	// scattering in 1/A

	t_mat matUB = prod_mm(matU, matB);
	return matUB;
}


/**
 * get the scattering plane from a UB matrix
 */
template<typename t_real = double>
std::array<ublas::vector<t_real>, 3> get_scattering_plane(
	const ublas::matrix<t_real>& UB,
	const tl::Lattice<t_real>& latt_real,
	t_real gu = 0., t_real gl = 0.)
{
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;

	t_mat B = tl::get_B(latt_real, true);
	t_mat B_inv;
	tl::inverse(B, B_inv);

	t_mat U = tl::prod_mm(UB, B_inv);

	std::array<ublas::vector<t_real>, 3> plane;
	for(int i = 0; i < 3; ++i)
	{
		t_vec vec = tl::make_vec<t_vec>({ U(i, 0), U(i, 1), U(i, 2) });
		vec = prod_mv(B_inv, vec);
		vec /= tl::veclen(vec);

		plane[i] = std::move(vec);
	}

	return plane;
}


/**
 * get the scattering plane from a UB matrix
 */
template<typename t_real = double>
std::array<ublas::vector<t_real>, 3> get_scattering_plane(
	const ublas::matrix<t_real>& UB,
	const std::array<t_real, 3>& lattice,
	const std::array<t_real, 3>& angles,
	t_real gu = 0., t_real gl = 0.)
{
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;

	tl::Lattice<t_real> latt(
		lattice[0], lattice[1], lattice[2],
		angles[0], angles[1], angles[2]);
	latt.RotateEuler(-gl, -gu, 0.);

	return get_scattering_plane<t_real>(UB, latt);
}
// -----------------------------------------------------------------------------


/**
 * hklE -> TAS angles
 * @see e.g.: https://doi.org/10.1107/S0021889805004875
 */
template<typename T = double>
void get_tas_angles(const Lattice<T>& lattice_real,
	const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2,
	T dKi, T dKf,
	T dh, T dk, T dl,
	bool bSense,
	T *pTheta, T *pTwoTheta,
	ublas::vector<T>* pVecQ = nullptr)
{
	// distance for point to be considered inside scattering plane
	static const T dDelta = std::cbrt(get_epsilon<T>());
	static const auto angs = get_one_angstrom<T>();
	static const auto rad = get_one_radian<T>();
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

	t_mat matUB = get_UB(lattice_real, _vec1, _vec2);

	t_vec vechkl = make_vec({ dh, dk, dl });
	t_vec vecQ = prod_mv(matUB, vechkl);
	if(pVecQ) *pVecQ = vecQ;

	if(std::fabs(vecQ[2]) > dDelta)
	{
		std::ostringstream ostrErr;
		ostrErr << "Position ("
			<< vechkl[0] << vechkl[1] << vechkl[2]
			<< ") is not in scattering plane.";
		throw Err(ostrErr.str());
	}

	T dQ = veclen(vecQ);
	*pTwoTheta = get_sample_twotheta(dKi/angs, dKf/angs, dQ/angs, bSense) / rad;
	T dKiQ = get_angle_ki_Q(dKi/angs, dKf/angs, dQ/angs, bSense) / rad;
	vecQ.resize(2, true);

	// sample rotation = angle between ki and first orientation reflex
	// (plus an arbitrary, but fixed constant)
	T dAngleKiOrient1 = dKiQ + vec_angle(vecQ);
	*pTheta = -dAngleKiOrient1 + get_pi<T>()/T(2);	// a3 convention would be: kiorient1 - pi
}


/**
 * TAS angles -> hklE
 * @see e.g.: https://doi.org/10.1107/S0021889805004875
 */
template<typename T = double>
void get_hkl_from_tas_angles(const Lattice<T>& lattice_real,
	const ublas::vector<T>& _vec1, const ublas::vector<T>& _vec2,
	T dm, T da, T th_m, T th_a, T _th_s, T _tt_s,
	bool bSense_m, bool bSense_a, bool bSense_s,
	T* h, T* k, T* l,
	T* pki=0, T* pkf=0, T* pE=0, T* pQ=0,
	ublas::vector<T>* pVecQ = nullptr)
{
	static const auto angs = get_one_angstrom<T>();
	static const auto rad = get_one_radian<T>();
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;

	T th_s = _th_s;
	T tt_s = _tt_s;
	/*if(!bSense_s)
	{
		th_s = -th_s;
		tt_s = -tt_s;
	}*/

	T ki = get_mono_k(th_m*rad, dm*angs, bSense_m)*angs;
	T kf = get_mono_k(th_a*rad, da*angs, bSense_a)*angs;
	T E = get_energy_transfer(ki/angs, kf/angs) / get_one_meV<T>();
	T Q = get_sample_Q(ki/angs, kf/angs, tt_s*rad)*angs;
	T kiQ = get_angle_ki_Q(ki/angs, kf/angs, Q/angs, bSense_s) / rad;
	T Qvec1 = get_pi<T>()/T(2) - th_s - kiQ;	// a3 convention

	t_mat matUB = get_UB(lattice_real, _vec1, _vec2);
	t_mat matUBinv;
	if(!inverse(matUB, matUBinv))
		throw Err("Cannot invert UB.");

	t_mat rot = rotation_matrix_3d_z(Qvec1);
	t_vec vecQ = prod_mv(rot, make_vec({ Q, 0., 0. }));
	t_vec vechkl = prod_mv(matUBinv, vecQ);

	if(pVecQ) *pVecQ = vecQ;

	if(vechkl.size() != 3)
		throw Err("Cannot determine hkl.");

	*h = vechkl[0];
	*k = vechkl[1];
	*l = vechkl[2];

	if(pki) *pki = ki;
	if(pkf) *pkf = kf;
	if(pE) *pE = E;
	if(pQ) *pQ = Q;
}


// -----------------------------------------------------------------------------



/**
 * hklE -> quaternion of crystal orientation
 * also rotate towards up vector if given
 */
template<typename T = double>
math::quaternion<T> get_hkl_orient(const Lattice<T>& lattice_real,
	T dh, T dk, T dl,				// original Bragg peak
	T dh_new, T dk_new, T dl_new,	// Bragg peak to rotate into
	T up_h = 0, T up_k = 0, T up_l = 1,
	ublas::vector<T>* pvecG = nullptr,
	const ublas::vector<T>& vecQx_rlu = make_vec<ublas::vector<T>>({ T(1), T(0), T(0) }),	// front
	const ublas::vector<T>& vecQy_rlu = make_vec<ublas::vector<T>>({ T(0), T(1), T(0) }))	// side
{
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;
	using t_quat = math::quaternion<T>;

	// scattering plane normal
	ublas::vector<T> vecQz_rlu = cross_3(vecQx_rlu, vecQy_rlu);

	t_mat matUB = get_UB(lattice_real, vecQx_rlu, vecQy_rlu);

	t_vec vecG_rlu = make_vec({ dh, dk, dl });
	t_vec vecG = prod_mv(matUB, vecG_rlu);

	t_vec vecGnew_rlu = make_vec({ dh_new, dk_new, dl_new });
	t_vec vecGnew = prod_mv(matUB, vecGnew_rlu);

	// scattering plane
	t_vec vecQx = prod_mv(matUB, vecQx_rlu);
	t_vec vecQz = prod_mv(matUB, vecQz_rlu);


	if(pvecG) *pvecG = vecG;

	// quaternion to rotate G into Gnew
	t_quat quatRot = rotation_quat<t_quat, t_vec, T>(vecG, vecGnew);


	// ------------------------------------------------------------------------
	t_vec vechklUp = make_vec({ up_h, up_k, up_l });
	t_vec vecUp = prod_mv(matUB, vechklUp);

	// get angle between vecUp and vecQz
	t_vec vec0 = make_vec<t_vec>({ T(0), T(0), T(0) });
	Plane<T> plane(vec0, vecQx, vecQz);
	T angle = plane.GetAngle(vecUp);

	// rotate around G towards up vector
	t_quat quatRotAroundG = rotation_quat<t_quat, t_vec, T>(vecG, angle);
	quatRot *= quatRotAroundG;
	// ------------------------------------------------------------------------


	return quatRot;
}



/**
 * hklE -> euler angles
 */
template<typename T = double>
math::quaternion<T> get_euler_angles(const Lattice<T>& lattice_real,
	T dKi, T dh, T dk, T dl,		// original Bragg peak
	T dh_new, T dk_new, T dl_new,	// Bragg peak to rotate into
	T *pRelTheta, T *pThetaX, T *pTwoTheta, T *pChi, T *pPsi,
	T up_h = 0, T up_k = 0, T up_l = 1,
	const ublas::vector<T>& vecQx_rlu = make_vec<ublas::vector<T>>({ T(1), T(0), T(0) }),	// front
	const ublas::vector<T>& vecQy_rlu = make_vec<ublas::vector<T>>({ T(0), T(1), T(0) }))	// side
{
	static const auto angs = get_one_angstrom<T>();
	static const auto rad = get_one_radian<T>();
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;
	using t_quat = math::quaternion<T>;

	t_vec vecG;
	// matrix to rotate G into Q
	t_quat quatRot = get_hkl_orient<T>(lattice_real,
		dh, dk, dl,
		dh_new, dk_new, dl_new,
		up_h, up_k, up_l,
		&vecG,
		vecQx_rlu, vecQy_rlu);

	// two theta
	T dG = veclen(vecG);
	*pTwoTheta = get_sample_twotheta(dKi/angs, dKi/angs, dG/angs, 1) / rad;


	// theta angle of vecQx
	bool bSense = true;
	t_mat matUB = get_UB(lattice_real, vecQx_rlu, vecQy_rlu);
	t_vec vecQx = prod_mv(matUB, vecQx_rlu);
	T dKiQ = get_angle_ki_Q(dKi/angs, dKi/angs, dG/angs, bSense) / rad;

	vecQx.resize(2, true);
	T dAngleKiOrient1 = dKiQ; //+ vec_angle(vecQx);
	*pThetaX = dAngleKiOrient1 - get_pi<T>()/T(2);
	if(bSense) *pThetaX = -*pThetaX;


	// rel. euler angles of normalized hkl direction towards vecQx
	std::vector<T> vecEuler = quat_to_euler(quatRot);
	*pRelTheta = vecEuler[2];
	*pChi = vecEuler[1];
	*pPsi = vecEuler[0];

	return quatRot;
}



// -----------------------------------------------------------------------------


template<typename T = double>
std::ostream& operator<<(std::ostream& ostr, const Lattice<T>& lat)
{
	ostr << "a = " << lat.GetA() << ", ";
	ostr << "b = " << lat.GetB() << ", ";
	ostr << "c = " << lat.GetC() << "; ";

	ostr << "alpha = " << tl::r2d(lat.GetAlpha()) << " deg, ";
	ostr << "beta = " << tl::r2d(lat.GetBeta()) << " deg, ";
	ostr << "gamma = " << tl::r2d(lat.GetGamma()) << " deg; ";

	return ostr;
}

}

#endif
