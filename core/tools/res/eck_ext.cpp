/**
 * implementation of Mechthild's extended version of the Eckold-Sobolev method (see [end25])
 * WARNING: WORK IN PROGRESS
 *
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2025
 * @license GPLv2
 *
 * @desc for extended algorithm: [end25] M. Enderle, personal communication (7/apr/2025)
 * @desc for original algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
 * @desc for vertical scattering modification: [eck20] G. Eckold, personal communication, 2020.
 * @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include "eck_ext.h"
#include "ellipse.h"

#include "tlibs/math/linalg.h"
#include "tlibs/math/math.h"
#include "ellipse.h"

#include <tuple>
#include <future>
#include <string>
#include <iostream>


typedef t_real_reso t_real;
typedef ublas::matrix<t_real> t_mat;
typedef ublas::vector<t_real> t_vec;

using angle = tl::t_angle_si<t_real>;
using wavenumber = tl::t_wavenumber_si<t_real>;
using energy = tl::t_energy_si<t_real>;
using length = tl::t_length_si<t_real>;
using volume = tl::t_volume_si<t_real>;
using inv_length = tl::t_length_inverse_si<t_real>;

static const auto angs = tl::get_one_angstrom<t_real>();
static const auto rads = tl::get_one_radian<t_real>();
static const auto meV = tl::get_one_meV<t_real>();
static const t_real pi = tl::get_pi<t_real>();
static const t_real sig2fwhm = tl::get_SIGMA2FWHM<t_real>();
static const t_real fwhm2sig = tl::get_FWHM2SIGMA<t_real>();


enum EckQE : std::size_t
{
	ECK_Q_X = 0, ECK_Q_Y, ECK_Q_Z,
	ECK_OM, ECK_K_Y, ECK_K_Z,

	ECK_NUM_QE
};


enum EckKiKfIdx : std::size_t
{
	ECK_KI_X = 0, ECK_KI_Y, ECK_KI_Z,
	ECK_KF_X, ECK_KF_Y, ECK_KF_Z,

	ECK_NUM_KIKF
};


static std::tuple<t_mat, t_mat, t_mat, t_real>
get_mono_vals(const length& src_w, const length& src_h,
	const length& mono_w, const length& mono_h,
	const length& dist_vsrc_mono, const length& dist_hsrc_mono,
	const length& dist_mono_sample,
	const wavenumber& ki, const angle& thetam,
	const angle& coll_h_pre_mono, const angle& coll_h_pre_sample,
	const angle& coll_v_pre_mono, const angle& coll_v_pre_sample,
	const angle& mono_mosaic, const angle& mono_mosaic_v,
	const inv_length& inv_mono_curvh, const inv_length& inv_mono_curvv,
	t_real dRefl)
{
	const t_real s_th_m = units::abs(units::sin(thetam));
	const t_real t_th_m = units::tan(thetam);

	// A matrix: formula 26 in [eck14]
	t_mat A = ublas::identity_matrix<t_real>(3);
	{
		const auto A_t0 = t_real(1) / mono_mosaic;
		const auto A_tx = inv_mono_curvh*dist_mono_sample / s_th_m;
		const auto A_t1 = A_t0*A_tx;

		A(0, 0) = t_real(0.5) / (ki*angs*ki*angs) * t_th_m*t_th_m *
		(
/*a*/			+ units::pow<2>(t_real(2)/coll_h_pre_mono) *rads*rads
/*b*/			+ units::pow<2>(t_real(2)*dist_hsrc_mono/src_w)
/*c*/			+ A_t0*A_t0 * rads*rads
		);
		A(0, 1) = A(1, 0) = t_real(0.5) / (ki*angs*ki*angs) * t_th_m *
		(
/*w*/			+ t_real(2)*tl::my_units_pow2(t_real(1)/coll_h_pre_mono) *rads*rads
/*x*/			+ t_real(2)*dist_hsrc_mono*(dist_hsrc_mono-dist_mono_sample)/(src_w*src_w)
/*y*/			+ A_t0*A_t0 * rads*rads
/*z*/			- A_t0*A_t1 * rads*rads
		);
		A(1, 1) = t_real(0.5) / (ki*angs*ki*angs) *
		(
/*1*/			+ units::pow<2>(t_real(1)/coll_h_pre_mono) *rads*rads
/*2*/			+ units::pow<2>(t_real(1)/coll_h_pre_sample) *rads*rads
/*3*/			+ units::pow<2>((dist_hsrc_mono-dist_mono_sample)/src_w)
/*4*/			+ units::pow<2>(dist_mono_sample/(mono_w*s_th_m))

/*5*/			+ A_t0*A_t0 * rads*rads
/*6*/			- t_real(2)*A_t0*A_t1 * rads*rads
/*7*/			+ A_t1*A_t1 * rads*rads
		);
	}

	// Av matrix: formula 38 in [eck14]
	// some typos in paper leading to the (false) result of a better Qz resolution when focusing
	// => trying to match terms in Av with corresponding terms in A
	// corresponding pre-mono terms commented out in Av, as they are not considered there
	t_mat Av(2, 2);
	{
		const auto Av_t0 = t_real(0.5) / (mono_mosaic_v*s_th_m);
		const auto Av_t1 = inv_mono_curvv*dist_mono_sample / mono_mosaic_v;

		Av(0, 0) = t_real(0.5) / (ki*angs*ki*angs) *
		(
/*1*/		//	+ units::pow<2>(t_real(1) / coll_v_pre_mono) *rads*rads	// missing in paper?
/*2*/			+ units::pow<2>(t_real(1) / coll_v_pre_sample) *rads*rads
/*~3*/			+ units::pow<2>(dist_mono_sample / src_h)
/*4*/			+ units::pow<2>(dist_mono_sample / mono_h)

/*5*/			+ Av_t0*Av_t0 * rads*rads
/*6*/			- t_real(2)*Av_t0*Av_t1 * rads*rads     // typo in paper?
/*7*/			+ Av_t1*Av_t1 * rads*rads               // missing in paper?
		);
		Av(0, 1) = Av(1, 0) = t_real(0.5) / (ki*angs*ki*angs) *
		(
/*w*/		//	- units::pow<2>(1./coll_v_pre_mono) *rads*rads   // missing in paper?
/*~x*/			+ dist_vsrc_mono*dist_mono_sample/(src_h*src_h)
/*y*/			- Av_t0*Av_t0 * rads*rads
/*z*/			+ Av_t0*Av_t1 * rads*rads
		);
		Av(1, 1) = t_real(0.5) / (ki*angs*ki*angs) *
		(
/*a*/			+ units::pow<2>(t_real(1)/coll_v_pre_mono) *rads*rads
/*b*/			+ units::pow<2>(dist_vsrc_mono/src_h)
/*c*/			+ Av_t0*Av_t0 *rads*rads
		);
	}

	// B matrix from equ. 2.3 in [end25] which corresponds to the B vector from equ. 27 in [eck14]
	t_mat B = ublas::zero_matrix<t_real>(3, 3);
	{
		const auto B_t0 = inv_mono_curvh / (mono_mosaic*mono_mosaic*s_th_m);

		B(0, 1) = 1. / ki * t_th_m *
		(
/*i*/			+ t_real(2)*dist_hsrc_mono / (src_w*src_w)
/*j*/			+ B_t0 *rads*rads
		);
		B(1, 1) = 1. / ki *
		(
/*r*/			- dist_mono_sample / (units::pow<2>(mono_w*s_th_m))
/*s*/			+ B_t0 * rads*rads
/*t*/			- B_t0 * rads*rads * inv_mono_curvh*dist_mono_sample / s_th_m
/*u*/			+ (dist_hsrc_mono-dist_mono_sample) / (src_w*src_w)
		);
	}

	// Bzz component from equ. 2.3 in [end25] which corresponds to the Bv vector from equ. 39 in [eck14]
	t_vec Bv(2);
	{
		const auto Bv_t0 = inv_mono_curvv / (mono_mosaic_v*mono_mosaic_v);

		Bv(0) = 1. / ki * t_real(-1.) *
		(
/*r*/			+ dist_mono_sample / (mono_h*mono_h)    // typo in paper?
/*~s*/			- t_real(0.5)*Bv_t0 *rads*rads / s_th_m
/*~t*/			+ Bv_t0 * rads*rads * inv_mono_curvv*dist_mono_sample
/*~u*/			+ dist_mono_sample / (src_h*src_h)      // typo in paper?
		);
		Bv(1) = 1. / ki * t_real(-1.) *
		(
/*i*/			+ dist_vsrc_mono / (src_h*src_h)        // typo in paper?
/*j*/			+ t_real(0.5)*Bv_t0/s_th_m * rads*rads
		);
	}


	// C matrix from equ. 2.12 in [end25] which corresponds to the C scalar from equ. 28 in [eck14]
	t_mat C = ublas::zero_matrix<t_real>(3, 3);
	C(1, 1) = t_real(0.5) * angs*angs*
	(
		t_real(1)/(src_w*src_w) +
		units::pow<2>(t_real(1)/(mono_w*s_th_m)) +
		units::pow<2>(inv_mono_curvh/(mono_mosaic * s_th_m)) * rads*rads
	);

	// Czz component from equ. 2.14 in [end25] which corresponds to the Cv scalar from equ. 40 in [eck14]
	t_real Cv = t_real(0.5) * angs*angs*
	(
		t_real(1)/(src_h*src_h) +
		t_real(1)/(mono_h*mono_h) +
		units::pow<2>(inv_mono_curvv/mono_mosaic_v) * rads*rads
	);


	// Bzz component from equ. 2.3 in [end25] which corresponds to the Bv vector from equ. 42 in [eck14]
	A(2, 2) = Av(0, 0) - Av(0, 1)*Av(0, 1)/Av(1, 1);
	B(2, 2) = Bv[0] - Bv[1]*Av(0, 1)/Av(1, 1);
	// Czz component from equ. 2.14 in [end25] which corresponds to the Cv scalar from equ. 40 in [eck14]
	C(2, 2) = Cv - t_real(0.25)*Bv[1]*Bv[1]/Av(1, 1);  // typo in paper? (thanks to F. Bourdarot for pointing this out)


	// [eck14], equ. 54, in th paper the sqrt factor is missing in some other equations
	t_real refl = dRefl * std::sqrt(pi / (Av(1, 1)));

	return std::make_tuple(A, B, C, refl);
}


ResoResults calc_eck_ext(const EckParams& eck)
{
	// if the user moved the scattering angle to the other
	// side of the scattering sense indicated by the flag
	t_real manually_changed_sense = t_real(1);
	if(eck.twotheta/rads < t_real(0))
		manually_changed_sense = t_real(-1);

	angle twotheta = eck.twotheta * eck.dsample_sense;
	angle thetaa = eck.thetaa * eck.dana_sense;
	angle thetam = eck.thetam * eck.dmono_sense;
	angle ki_Q = eck.angle_ki_Q * eck.dsample_sense * manually_changed_sense;
	angle kf_Q = eck.angle_kf_Q * eck.dsample_sense * manually_changed_sense;
	//kf_Q = ki_Q + twotheta;


	// --------------------------------------------------------------------
	// mono/ana focus
	// --------------------------------------------------------------------
	length mono_curvh = eck.mono_curvh, mono_curvv = eck.mono_curvv;
	length ana_curvh = eck.ana_curvh, ana_curvv = eck.ana_curvv;

	if(eck.bMonoIsOptimallyCurvedH)
		mono_curvh = tl::foc_curv(eck.dist_hsrc_mono, eck.dist_mono_sample, units::abs(t_real(2)*thetam), false);
	if(eck.bMonoIsOptimallyCurvedV)
		mono_curvv = tl::foc_curv(eck.dist_vsrc_mono, eck.dist_mono_sample, units::abs(t_real(2)*thetam), true);
	if(eck.bAnaIsOptimallyCurvedH)
		ana_curvh = tl::foc_curv(eck.dist_sample_ana, eck.dist_ana_det, units::abs(t_real(2)*thetaa), false);
	if(eck.bAnaIsOptimallyCurvedV)
		ana_curvv = tl::foc_curv(eck.dist_sample_ana, eck.dist_ana_det, units::abs(t_real(2)*thetaa), true);

	//mono_curvh *= eck.dmono_sense; mono_curvv *= eck.dmono_sense;
	//ana_curvh *= eck.dana_sense; ana_curvv *= eck.dana_sense;

	inv_length inv_mono_curvh = t_real(0)/angs, inv_mono_curvv = t_real(0)/angs;
	inv_length inv_ana_curvh = t_real(0)/angs, inv_ana_curvv = t_real(0)/angs;

	if(eck.bMonoIsCurvedH)
		inv_mono_curvh = t_real(1)/mono_curvh;
	if(eck.bMonoIsCurvedV)
		inv_mono_curvv = t_real(1)/mono_curvv;
	if(eck.bAnaIsCurvedH)
		inv_ana_curvh = t_real(1)/ana_curvh;
	if(eck.bAnaIsCurvedV)
		inv_ana_curvv = t_real(1)/ana_curvv;
	// --------------------------------------------------------------------


	angle coll_h_pre_mono = eck.coll_h_pre_mono;
	angle coll_v_pre_mono = eck.coll_v_pre_mono;

	if(eck.bGuide)
	{
		const length lam = tl::k2lam(eck.ki);

		coll_h_pre_mono = lam*(eck.guide_div_h/angs);
		coll_v_pre_mono = lam*(eck.guide_div_v/angs);
	}


	// -------------------------------------------------------------------------

	// - if the instruments works in kf=const mode and the scans are counted for
	//   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
	// - if the instrument works in ki=const mode the kf^3 factor is needed.
	const auto tupScFact = get_scatter_factors(eck.flags, eck.thetam, eck.ki, eck.thetaa, eck.kf);

	t_real dmono_refl = eck.dmono_refl * std::get<0>(tupScFact);
	t_real dana_effic = eck.dana_effic * std::get<1>(tupScFact);
	if(eck.mono_refl_curve)
		dmono_refl *= (*eck.mono_refl_curve)(eck.ki);
	if(eck.ana_effic_curve)
		dana_effic *= (*eck.ana_effic_curve)(eck.kf);
	t_real dxsec = std::get<2>(tupScFact);
	t_real dmonitor = std::get<3>(tupScFact);


	// if no vertical mosaic is given, use the horizontal one
	angle mono_mosaic_v = eck.mono_mosaic_v;
	angle ana_mosaic_v = eck.ana_mosaic_v;
	angle sample_mosaic_v = eck.sample_mosaic_v;
	if(tl::float_equal<t_real>(mono_mosaic_v/rads, 0.), 0.)
		mono_mosaic_v = eck.mono_mosaic;
	if(tl::float_equal<t_real>(ana_mosaic_v/rads, 0.), 0.)
		ana_mosaic_v = eck.ana_mosaic;
	if(tl::float_equal<t_real>(sample_mosaic_v/rads, 0.), 0.)
		sample_mosaic_v = eck.sample_mosaic;


	// sample position
	t_vec sample_pos(3);
	sample_pos[0] = eck.pos_x / angs;
	sample_pos[1] = eck.pos_y / angs;
	sample_pos[2] = eck.pos_z / angs;


	//--------------------------------------------------------------------------
	// mono & ana calculations, equ. 43 in [eck14]
	//--------------------------------------------------------------------------
	std::launch lpol = std::launch::async;
	std::future<std::tuple<t_mat, t_mat, t_mat, t_real>> futMono
		= std::async(lpol, get_mono_vals,
			fwhm2sig * eck.src_w, fwhm2sig * eck.src_h,
			fwhm2sig * eck.mono_w, fwhm2sig * eck.mono_h,
			eck.dist_vsrc_mono,  eck.dist_hsrc_mono,
			eck.dist_mono_sample,
			eck.ki, thetam,
			fwhm2sig * coll_h_pre_mono, fwhm2sig * eck.coll_h_pre_sample,
			fwhm2sig * coll_v_pre_mono, fwhm2sig * eck.coll_v_pre_sample,
			fwhm2sig * eck.mono_mosaic, fwhm2sig * mono_mosaic_v,
			inv_mono_curvh, inv_mono_curvv,
			dmono_refl);

	std::future<std::tuple<t_mat, t_mat, t_mat, t_real>> futAna
		= std::async(lpol, get_mono_vals,
			fwhm2sig * eck.det_w, fwhm2sig * eck.det_h,
			fwhm2sig * eck.ana_w, fwhm2sig * eck.ana_h,
			eck.dist_ana_det, eck.dist_ana_det,
			eck.dist_sample_ana,
			eck.kf, -thetaa,
			fwhm2sig * eck.coll_h_post_ana, fwhm2sig * eck.coll_h_post_sample,
			fwhm2sig * eck.coll_v_post_ana, fwhm2sig * eck.coll_v_post_sample,
			fwhm2sig * eck.ana_mosaic, fwhm2sig * ana_mosaic_v,
			inv_ana_curvh, inv_ana_curvv,
			dana_effic);
	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	// get mono & ana results
	//--------------------------------------------------------------------------
	std::tuple<t_mat, t_mat, t_mat,t_real> tupMono = futMono.get();
	const t_mat& A = std::get<0>(tupMono);
	const t_mat& B = std::get<1>(tupMono);
	const t_mat& C = std::get<2>(tupMono);
	const t_real& dReflM = std::get<3>(tupMono);

	std::tuple<t_mat, t_mat, t_mat, t_real> tupAna = futAna.get();
	t_mat& E = std::get<0>(tupAna);
	t_mat& F = std::get<1>(tupAna);
	t_mat& G = std::get<2>(tupAna);
	const t_real& dReflA = std::get<3>(tupAna);

	// vertical scattering in kf axis, formula from [eck20]
	bool bKfVertical = (eck.angle_kf/rads > t_real(0));
	if(bKfVertical)
	{
		t_mat matTvert = tl::rotation_matrix_3d_x(-eck.angle_kf/rads);

		// T_vert has to be applied at the same positions in the formulas as Dtwotheta, see eck.cpp
		E = tl::transform(E, matTvert, true);
		F = tl::transform(F, matTvert, true);
		G = tl::transform(G, matTvert, true);

		//sample_pos = ublas::prod(matTvert, sample_pos);
	}
	//--------------------------------------------------------------------------


	ResoResults res;

	res.Q_avg.resize(4);
	res.Q_avg[0] = eck.Q*angs;
	res.Q_avg[1] = 0.;
	res.Q_avg[2] = 0.;
	res.Q_avg[3] = eck.E/meV;


	// equ. 4 & equ. 53 in [eck14]
	const t_real dE = (eck.ki*eck.ki - eck.kf*eck.kf) / (t_real(2)*eck.Q*eck.Q);
	const wavenumber kipara = eck.Q * (t_real(0.5) + dE);
	const wavenumber kfpara = eck.Q - kipara;
	wavenumber kperp = tl::my_units_sqrt<wavenumber>(units::abs(kipara*kipara - eck.ki*eck.ki));
	kperp *= eck.dsample_sense * manually_changed_sense;

	const t_real ksq2E = tl::get_KSQ2E<t_real>();

	// trafo, equ. 52 in [eck14]
	t_mat T = ublas::zero_matrix<t_real>(ECK_NUM_QE, ECK_NUM_KIKF);
	T(ECK_Q_X, ECK_KI_X) = T(ECK_Q_Y, ECK_KI_Y) = T(ECK_Q_Z, ECK_KI_Z) = +1.;
	T(ECK_Q_X, ECK_KF_X) = T(ECK_Q_Y, ECK_KF_Y) = T(ECK_Q_Z, ECK_KF_Z) = -1.;
	T(ECK_OM, ECK_KI_X) = t_real(2)*ksq2E * kipara * angs;
	T(ECK_OM, ECK_KF_X) = t_real(2)*ksq2E * kfpara * angs;
	T(ECK_OM, ECK_KI_Y) = t_real(2)*ksq2E * kperp * angs;
	T(ECK_OM, ECK_KF_Y) = t_real(-2)*ksq2E * kperp * angs;
	T(ECK_K_Y, ECK_KI_Y) = T(ECK_K_Z, ECK_KI_Z) = (0.5 - dE);
	T(ECK_K_Y, ECK_KF_Y) = T(ECK_K_Z, ECK_KF_Z) = (0.5 + dE);
	t_mat Tinv;
	if(!tl::inverse(T, Tinv))
	{
		res.bOk = false;
		res.strErr = "Matrix T cannot be inverted.";
		return res;
	}

	// equ. 54 in [eck14]
	t_mat Dalph_i = tl::rotation_matrix_3d_z(-ki_Q/rads);
	t_mat Dalph_f = tl::rotation_matrix_3d_z(-kf_Q/rads);
	t_mat Dtwotheta = tl::rotation_matrix_3d_z(-twotheta/rads);
	t_mat Arot = tl::transform(A, Dalph_i, true);
	t_mat Erot = tl::transform(E, Dalph_f, true);

	t_mat matAE = ublas::zero_matrix<t_real>(Arot.size1() + Erot.size1(), Arot.size2() + Erot.size2());
	tl::submatrix_copy(matAE, Arot, 0, 0);
	tl::submatrix_copy(matAE, Erot, Arot.size1(), Arot.size2());

	// U1 matrix
	t_mat U1 = tl::transform(matAE, Tinv, true);  // typo in paper in quadric trafo in equ. 54 (top)?

	// V matrix from equ. 2.9 [end25], corresponds to V1 vector in [eck14]
	t_mat matBrot = ublas::prod(ublas::trans(Dalph_i), B);
	// the transpose of Dtwotheta is part of Dalph_f^T
	t_mat matFrot0 = ublas::prod(ublas::trans(Dalph_f), F);
	t_mat matFrot = ublas::prod(matFrot0, Dtwotheta);
	t_mat matBF = ublas::zero_matrix<t_real>(matBrot.size1() + matFrot.size1(), matBrot.size2());
	tl::submatrix_copy(matBF, matBrot, 0, 0);
	tl::submatrix_copy(matBF, matFrot, matBrot.size1(), 0);
	t_mat matV = ublas::prod(ublas::trans(Tinv), matBF);


	//--------------------------------------------------------------------------
	// integrate last 2 vars -> equs 57 & 58 in [eck14]
	//--------------------------------------------------------------------------
	t_mat U2 = quadric_proj(U1, ECK_K_Z);
	// careful: factor -0.5*... missing in U matrix compared to normal gaussian!
	t_mat U = t_real(2) * quadric_proj(U2, ECK_K_Y);

	// P matrix from equ. 2.21 in [end25]
	// quadric_proj_mat() gives the same as equ. 2.21 in [end25]
	t_mat matV2 = quadric_proj_mat(matV, U1, ECK_K_Z);
	t_mat matP = quadric_proj_mat(matV2, U2, ECK_K_Y);

	// K matrix from equ. 2.11 in [end25]
	t_mat matGrot = tl::transform(G, Dtwotheta, true);
	t_mat matK = C + matGrot;

	// equ. 2.19 in [end25], corresponds to equ. 57 & 58 in [eck14], "W -= ..." in eck.cpp
	// squares in Vs missing in paper? (thanks to F. Bourdarot for pointing this out)
	for(std::size_t i = 0; i < 3; ++i)
		for(std::size_t j = 0; j < 3; ++j)
			matK(i, j) -= t_real(0.25) *
				(matV(ECK_K_Z, i) * matV(ECK_K_Z, j) / U1(ECK_K_Z, ECK_K_Z)
				+ matV2(ECK_K_Y, i) * matV2(ECK_K_Y, j) / U2(ECK_K_Y, ECK_K_Y));

	// C_all,0 in [end25], equ. 1.1, 2.1
	t_real Z0 = std::sqrt(pi/std::abs(U1(ECK_K_Z, ECK_K_Z)))
		* std::sqrt(pi/std::abs(U2(ECK_K_Y, ECK_K_Y)));
	t_real Z = dReflM * dReflA * Z0;
	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	// include sample mosaic, see also cn.cpp
	//--------------------------------------------------------------------------
	const t_real mos_Q_sq[2] =
	{
		(fwhm2sig*eck.sample_mosaic/rads * eck.Q*angs) * (fwhm2sig*eck.sample_mosaic/rads * eck.Q*angs),
		(fwhm2sig*sample_mosaic_v/rads * eck.Q*angs) * (fwhm2sig*sample_mosaic_v/rads * eck.Q*angs)
	};

	// sample mosaic, gives the same as equs. 3.3 and 7.7 in [end25]
	t_mat mosaic = tl::submatrix_wnd(U, 2, 2, 1, 1);
	mosaic(0, 0) += 1. / mos_Q_sq[0];  // horizontal
	mosaic(1, 1) += 1. / mos_Q_sq[1];  // vertical

	// equs. 3.4 and 7.14 in [end25]
	Z *= pi*pi / std::sqrt(tl::determinant(mosaic));
	Z *= fwhm2sig / std::sqrt(2. * pi * mos_Q_sq[0] * mos_Q_sq[1]);

	for(int comp = 0; comp < 2; ++comp)  // horizontal and vertical components
	{
		t_vec Pvec = tl::get_row<t_vec>(matP, comp + 1);
		t_vec Uvec = tl::get_column<t_vec>(U, comp + 1);
		t_real Mnorm = 1./mos_Q_sq[comp] + U(comp + 1, comp + 1);

		// gives the same as equs. 3.5 and 7.8 in [end25]
		matK -= 0.25 * ublas::outer_prod(Pvec, Pvec) / Mnorm;
		// gives the same as equs. 3.7 and 7.10 in [end25]
		matP -= ublas::outer_prod(Uvec, Pvec) / Mnorm;
		// gives the same as equs. 3.6 and 7.9 in [end25]
		U -= ublas::outer_prod(Uvec, Uvec) / Mnorm;
	}
	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	// integrate over sample shape
	//--------------------------------------------------------------------------
	// TODO: sample rotation
	t_mat T_E = ublas::identity_matrix<t_real>(3);

	// cuboid sample integration, equ. 6.6 in [end25]
	t_real sample_var[3] = { pi, pi, pi };
	volume V_sample = eck.sample_w_perpq * eck.sample_w_q * eck.sample_h;
	if(!eck.bSampleCub)
	{
		// cylindrical sample integration, equ. 6.3 in [end25]
		sample_var[0] = sample_var[1] = 4.;
		V_sample = pi * 0.5*eck.sample_w_perpq * 0.5*eck.sample_w_q * eck.sample_h;
	}

	// equs. 4.4, 6.3, and 6.6 in [end25]
	t_mat matN = matK;
	matN += tl::transform<t_mat>(tl::diag_matrix<t_mat>({
		sample_var[0]/(eck.sample_w_perpq*eck.sample_w_perpq /angs/angs),
		sample_var[1]/(eck.sample_w_q*eck.sample_w_q /angs/angs),
		sample_var[2]/(eck.sample_h*eck.sample_h /angs/angs) }),
		ublas::trans(T_E), true);

	t_real detN = tl::determinant(matN);
	t_mat Nadj = tl::adjugate(matN, true);

	// page 9 and equs. 7.11, 7.12 and 7.13 in [end25]
	U -= 0.25 / detN * tl::transform<t_mat>(Nadj, ublas::trans(matP), true);
	t_mat NadjK = ublas::prod(Nadj, matK);
	matP -= 1. / detN * ublas::prod(matP, NadjK);
	matK -= 1. / detN * ublas::prod(matK, NadjK);

	// page 9 in [end25]
	Z *= pi*pi*pi / detN;

	// normalise R0 to sample volume
	if(eck.flags & NORM_TO_SAMPLE)
		Z /= V_sample*V_sample / units::pow<6>(angs);
	//--------------------------------------------------------------------------


	// quadratic part of quadric (matrix U)
	res.reso = U;
	// linear and constant part of quadric (V and W in [eck14], equ. 2.2 in [end25])
	res.reso_v = ublas::prod(matP, sample_pos);
	res.reso_s = ublas::inner_prod(sample_pos, ublas::prod(matK, sample_pos));


	// prefactor and volume
	res.dResVol = tl::get_ellipsoid_volume(res.reso);
	bool use_monitor = (eck.flags & CALC_MON) != 0;

	res.dR0 = Z;
	if(use_monitor)
		res.dR0 /= dReflM;

	if(!(eck.flags & NORM_TO_RESVOL))
	{
		// missing volume prefactor to normalise gaussian,
		// cf. equ. 56 in [eck14] to  equ. 1 in [pop75] and equ. A.57 in [mit84]
		// res.dR0 /= std::sqrt(std::abs(tl::determinant(res.reso))) / (2.*pi*2.*pi);
		res.dR0 *= res.dResVol * pi * t_real(3.);  // TODO: check
	}

	res.dR0 *= std::exp(-res.reso_s);
	res.dR0 *= dxsec * dmonitor;
	res.dR0 = std::abs(res.dR0);

	// Bragg widths
	const std::vector<t_real> vecFwhms = calc_bragg_fwhms(res.reso);
	std::copy(vecFwhms.begin(), vecFwhms.end(), res.dBraggFWHMs);

	if(tl::is_nan_or_inf(res.dR0) || tl::is_nan_or_inf(res.reso))
	{
		res.strErr = "Invalid result.";
		res.bOk = false;
		return res;
	}

	res.bOk = true;
	return res;
}
