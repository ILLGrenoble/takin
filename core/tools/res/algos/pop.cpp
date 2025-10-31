/**
 * popovici calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2
 *
 * @desc This is a reimplementation in C++ of the file rc_popma.m of the
 *		rescal5 package by Zinkin, McMorrow, Tennant, Farhi, and Wildes (ca. 1995-2007):
 *		http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/
 * @desc see: - [pop75] M. Popovici, Acta Cryst. A 31, 507 (1975), doi: 10.1107/S0567739475001088
 *            - [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
 *            - [zhe07] A. Zheludev, ResLib 3.4 manual (2007), https://ethz.ch/content/dam/ethz/special-interest/phys/solid-state-physics/neutron-scattering-and-magnetism-dam/images/research/manual.pdf
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "pop.h"
#include "../ellipse.h"

#include "tlibs/math/linalg.h"
#include "tlibs/math/math.h"

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
static const auto cm = tl::get_one_centimeter<t_real>();
static const t_real sig2fwhm = tl::get_SIGMA2FWHM<t_real>();
static const t_real pi = tl::get_pi<t_real>();


/**
 * position indices, see V vector definition in [pop75], p. 509
 */
enum PopPosIdx : std::size_t
{
	POP_SRC_Y = 0, POP_SRC_Z,                 // w, h
	POP_MONO_X, POP_MONO_Y, POP_MONO_Z,       // d, w, h
	POP_SAMPLE_X, POP_SAMPLE_Y, POP_SAMPLE_Z, // perp Q, para Q, h
	POP_ANA_X, POP_ANA_Y, POP_ANA_Z,          // d, w, h
	POP_DET_Y, POP_DET_Z,                     // w, h

	POP_NUM_POS
};


enum PopCoordIdx : std::size_t
{
	POP_MONO_H = 0, POP_MONO_V,
	POP_ANA_H, POP_ANA_V,

	POP_NUM_COORDS
};


enum PopCompIdx : std::size_t
{
	POP_PREMONO_H = 0, POP_PREMONO_V,
	POP_PRESAMPLE_H, POP_PRESAMPLE_V,
	POP_POSTSAMPLE_H, POP_POSTSAMPLE_V,
	POP_POSTANA_H, POP_POSTANA_V,

	POP_NUM_COMPS
};


enum PopKiKfIdx : std::size_t
{
	POP_KI_X = 0, POP_KI_Y, POP_KI_Z,
	POP_KF_X, POP_KF_Y, POP_KF_Z,

	POP_NUM_KIKF
};


/**
 * get corresponding ki position index from a kf index
 */
static inline PopPosIdx get_ki_pos(PopPosIdx posidx)
{
	switch(posidx)
	{
		case POP_ANA_X: return POP_MONO_X;
		case POP_ANA_Y: return POP_MONO_Y;
		case POP_ANA_Z: return POP_MONO_Z;
		case POP_DET_Y: return POP_SRC_Y;
		case POP_DET_Z: return POP_SRC_Z;
		default: return posidx;
	}
}


ResoResults calc_pop(const PopParams& pop)
{
	ResoResults res;

	res.Q_avg.resize(4);
	res.Q_avg[0] = pop.Q * angs;
	res.Q_avg[1] = 0.;
	res.Q_avg[2] = 0.;
	res.Q_avg[3] = pop.E / meV;

	// if the user moved the scattering angle to the other
	// side of the scattering sense indicated by the flag
	t_real manually_changed_sense = t_real(1);
	if(pop.twotheta/rads < t_real(0))
		manually_changed_sense = t_real(-1);

	angle twotheta = pop.twotheta * pop.dsample_sense;
	angle thetaa = pop.thetaa * pop.dana_sense;
	angle thetam = pop.thetam * pop.dmono_sense;
	angle ki_Q = pop.angle_ki_Q * pop.dsample_sense * manually_changed_sense;
	angle kf_Q = pop.angle_kf_Q * pop.dsample_sense * manually_changed_sense;

	// B matrix, [pop75], Appendix 1 -> U matrix in CN
	t_mat B_trafo_QE = get_trafo_dkidkf_dQdE(ki_Q, kf_Q, pop.ki, pop.kf);
	B_trafo_QE.resize(4, POP_NUM_KIKF, true);

	angle coll_h_pre_mono = pop.coll_h_pre_mono;
	angle coll_v_pre_mono = pop.coll_v_pre_mono;

	if(pop.bGuide)
	{
		length lam = tl::k2lam(pop.ki);

		coll_h_pre_mono = lam*(pop.guide_div_h/angs);
		coll_v_pre_mono = lam*(pop.guide_div_v/angs);
	}

	// if no vertical mosaic is given, use the horizontal one
	angle mono_mosaic_v = pop.mono_mosaic_v;
	angle sample_mosaic_v = pop.sample_mosaic_v;
	angle ana_mosaic_v = pop.ana_mosaic_v;
	if(tl::float_equal<t_real>(mono_mosaic_v/rads, 0.), 0.)
		mono_mosaic_v = pop.mono_mosaic;
	if(tl::float_equal<t_real>(sample_mosaic_v/rads, 0.), 0.)
		sample_mosaic_v = pop.sample_mosaic;
	if(tl::float_equal<t_real>(ana_mosaic_v/rads, 0.), 0.)
		ana_mosaic_v = pop.ana_mosaic;


	// --------------------------------------------------------------------
	// instrument property matrices
	// --------------------------------------------------------------------
	// collimator covariance matrix G, [pop75], Appendix 1
	t_mat G_collis = tl::zero_matrix(POP_NUM_COMPS, POP_NUM_COMPS);

	G_collis(POP_PREMONO_H, POP_PREMONO_H) =
		t_real(1) / (coll_h_pre_mono*coll_h_pre_mono /rads/rads);
	G_collis(POP_PRESAMPLE_H, POP_PRESAMPLE_H) =
		t_real(1) / (pop.coll_h_pre_sample*pop.coll_h_pre_sample /rads/rads);

	G_collis(POP_POSTSAMPLE_H, POP_POSTSAMPLE_H) =
		t_real(1) / (pop.coll_h_post_sample*pop.coll_h_post_sample /rads/rads);
	G_collis(POP_POSTANA_H, POP_POSTANA_H) =
		t_real(1) / (pop.coll_h_post_ana*pop.coll_h_post_ana /rads/rads);

	G_collis(POP_PREMONO_V, POP_PREMONO_V) =
		t_real(1) / (coll_v_pre_mono*coll_v_pre_mono /rads/rads);
	G_collis(POP_PRESAMPLE_V, POP_PRESAMPLE_V) =
		t_real(1) / (pop.coll_v_pre_sample*pop.coll_v_pre_sample /rads/rads);

	G_collis(POP_POSTSAMPLE_V, POP_POSTSAMPLE_V) =
		t_real(1) / (pop.coll_v_post_sample*pop.coll_v_post_sample /rads/rads);
	G_collis(POP_POSTANA_V, POP_POSTANA_V) =
		t_real(1) / (pop.coll_v_post_ana*pop.coll_v_post_ana /rads/rads);

	// crystal mosaic covariance matrix F, [pop75], Appendix 1
	t_mat F_mosaics = tl::zero_matrix(POP_NUM_COORDS, POP_NUM_COORDS);
	F_mosaics(POP_MONO_H, POP_MONO_H) =
		t_real(1)/(pop.mono_mosaic*pop.mono_mosaic /rads/rads);
	F_mosaics(POP_MONO_V, POP_MONO_V) =
		t_real(1)/(mono_mosaic_v*mono_mosaic_v /rads/rads);
	F_mosaics(POP_ANA_H, POP_ANA_H) =
		t_real(1)/(pop.ana_mosaic*pop.ana_mosaic /rads/rads);
	F_mosaics(POP_ANA_V, POP_ANA_V) =
		t_real(1)/(ana_mosaic_v*ana_mosaic_v /rads/rads);

	const t_real s_th_m = units::sin(thetam);
	const t_real c_th_m = units::cos(thetam);
	const t_real s_th_a = units::sin(thetaa);
	const t_real c_th_a = units::cos(thetaa);
	const t_real cot_th_m = c_th_m / s_th_m;
	const t_real cot_th_a = c_th_a / s_th_a;
	const t_real s_th_s = units::sin(t_real(0.5)*twotheta);
	const t_real c_th_s = units::cos(t_real(0.5)*twotheta);
	const t_real sign_z = -1.; // to compare with the calculations by F. Bourdarot -> should cancel out

	// A matrix, [pop75], Appendix 1
	auto get_div_monosample = [](t_real ki, t_real cot_th_m) -> std::tuple<t_real, t_real>
	{
		t_real div_mono = t_real(0.5) * ki * cot_th_m;
		t_real div_sample = t_real(-0.5) * ki * cot_th_m;
		return std::make_tuple(div_mono, div_sample);
	};

	t_mat A_div_kikf_trafo = ublas::zero_matrix<t_real>(POP_NUM_KIKF, POP_NUM_COMPS);

	auto div_ki = get_div_monosample(pop.ki*angs, cot_th_m);
	A_div_kikf_trafo(POP_KI_X, POP_PREMONO_H) = std::get<0>(div_ki);
	A_div_kikf_trafo(POP_KI_X, POP_PRESAMPLE_H) = std::get<1>(div_ki);
	A_div_kikf_trafo(POP_KI_Y, POP_PRESAMPLE_H) = pop.ki * angs;
	A_div_kikf_trafo(POP_KI_Z, POP_PRESAMPLE_V) = sign_z * pop.ki * angs;

	auto div_kf = get_div_monosample(pop.kf*angs, -cot_th_a);
	A_div_kikf_trafo(POP_KF_X, POP_POSTANA_H) = std::get<0>(div_kf);
	A_div_kikf_trafo(POP_KF_X, POP_POSTSAMPLE_H) = std::get<1>(div_kf);
	A_div_kikf_trafo(POP_KF_Y, POP_POSTSAMPLE_H) = pop.kf * angs;
	A_div_kikf_trafo(POP_KF_Z, POP_POSTSAMPLE_V) = sign_z * pop.kf * angs;

	// covariance matrix of component geometries, S, [pop75], Appendices 2 and 3
	// uniform distribution variance = L^2/12, see: https://en.wikipedia.org/wiki/Continuous_uniform_distribution
	// circular distribution variance = R^2/4 = D^2/16, see: https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
	// these factors are used to approximate a uniform or circular variance by a gaussian variance
	t_real var_uniform = t_real(1. / 12.);
	t_real var_circular = t_real(1. / 16.);
	t_real dMultSrc = pop.bSrcRect ? var_uniform : var_circular;
	t_real dMultMono = var_uniform;
	t_real dMultSample = pop.bSampleCub ? var_uniform : var_circular;
	t_real dMultAna = var_uniform;
	t_real dMultDet = pop.bDetRect ? var_uniform : var_circular;
	t_real dMultMonitor = pop.bMonitorRect ? var_uniform : var_circular;

	t_mat SI_geo = tl::zero_matrix(POP_NUM_POS, POP_NUM_POS);
	SI_geo(POP_SRC_Y, POP_SRC_Y) = dMultSrc * pop.src_w*pop.src_w /cm/cm;
	SI_geo(POP_SRC_Z, POP_SRC_Z) = dMultSrc * pop.src_h*pop.src_h /cm/cm;

	SI_geo(POP_MONO_X, POP_MONO_X) = dMultMono * pop.mono_thick*pop.mono_thick /cm/cm;
	SI_geo(POP_MONO_Y, POP_MONO_Y) = dMultMono * pop.mono_w*pop.mono_w /cm/cm;
	SI_geo(POP_MONO_Z, POP_MONO_Z) = dMultMono * pop.mono_h*pop.mono_h /cm/cm;

	SI_geo(POP_SAMPLE_X, POP_SAMPLE_X) = dMultSample * pop.sample_w_perpq*pop.sample_w_perpq /cm/cm;
	SI_geo(POP_SAMPLE_Y, POP_SAMPLE_Y) = dMultSample * pop.sample_w_q*pop.sample_w_q /cm/cm;
	SI_geo(POP_SAMPLE_Z, POP_SAMPLE_Z) = var_uniform * pop.sample_h*pop.sample_h /cm/cm;

	SI_geo(POP_ANA_X, POP_ANA_X) = dMultAna * pop.ana_thick*pop.ana_thick /cm/cm;
	SI_geo(POP_ANA_Y, POP_ANA_Y) = dMultAna * pop.ana_w*pop.ana_w /cm/cm;
	SI_geo(POP_ANA_Z, POP_ANA_Z) = dMultAna * pop.ana_h*pop.ana_h /cm/cm;

	SI_geo(POP_DET_Y, POP_DET_Y) = dMultDet * pop.det_w*pop.det_w /cm/cm;
	SI_geo(POP_DET_Z, POP_DET_Z) = dMultDet * pop.det_h*pop.det_h /cm/cm;

	SI_geo *= sig2fwhm*sig2fwhm;  // convert gaussian variance to gaussian fwhm

	t_mat S_geo;
	if(!tl::inverse_diag(SI_geo, S_geo))
	{
		res.bOk = false;
		res.strErr = "S matrix cannot be inverted.";
		return res;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// mono/ana focus
	// --------------------------------------------------------------------
	length mono_curvh = pop.mono_curvh, mono_curvv = pop.mono_curvv;
	length ana_curvh = pop.ana_curvh, ana_curvv = pop.ana_curvv;

	if(pop.bMonoIsOptimallyCurvedH)
		mono_curvh = tl::foc_curv(pop.dist_hsrc_mono, pop.dist_mono_sample, units::abs(t_real(2)*thetam), false);
	if(pop.bAnaIsOptimallyCurvedH)
		ana_curvh = tl::foc_curv(pop.dist_sample_ana, pop.dist_ana_det, units::abs(t_real(2)*thetaa), false);
	if(pop.bMonoIsOptimallyCurvedV)
		mono_curvv = tl::foc_curv(pop.dist_vsrc_mono, pop.dist_mono_sample, units::abs(t_real(2)*thetam), true);
	if(pop.bAnaIsOptimallyCurvedV)
		ana_curvv = tl::foc_curv(pop.dist_sample_ana, pop.dist_ana_det, units::abs(t_real(2)*thetaa), true);

	mono_curvh *= pop.dmono_sense;
	ana_curvh *= pop.dana_sense;
	mono_curvv *= pop.dmono_sense;
	ana_curvv *= pop.dana_sense;

	inv_length inv_mono_curvh = pop.bMonoIsCurvedH ? t_real(1)/mono_curvh : t_real(0)/cm;
	inv_length inv_ana_curvh = pop.bAnaIsCurvedH ? t_real(1)/ana_curvh : t_real(0)/cm;
	inv_length inv_mono_curvv = pop.bMonoIsCurvedV ? t_real(1)/mono_curvv : t_real(0)/cm;
	inv_length inv_ana_curvv = pop.bAnaIsCurvedV ? t_real(1)/ana_curvv : t_real(0)/cm;

	const auto tupScFact = get_scatter_factors(pop.flags, pop.thetam, pop.ki, pop.thetaa, pop.kf);

	// cf. [zhe07], p. 10, equ. 8
	t_real dmono_refl = pop.dmono_refl * std::get<0>(tupScFact) * pi / (t_real(2.) * s_th_m);
	t_real dana_effic = pop.dana_effic * std::get<1>(tupScFact) * pi / (t_real(2.) * s_th_a);
	if(pop.mono_refl_curve)
		dmono_refl *= (*pop.mono_refl_curve)(pop.ki);
	if(pop.ana_effic_curve)
		dana_effic *= (*pop.ana_effic_curve)(pop.kf);
	t_real dxsec = std::get<2>(tupScFact);
	t_real dmonitor = std::get<3>(tupScFact);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// trafo matrices
	// --------------------------------------------------------------------
	// T matrix to transform the mosaic cov. matrix, [pop75], Appendix 2
	auto get_mosaic_trafo = [](const length& dist_vsrc_mono,
		const length& dist_hsrc_mono, const length& dist_mono_sample,
		t_real s_th_m, t_real c_th_m, t_real s_th_s, t_real c_th_s,
		const inv_length& inv_curvh, const inv_length& inv_curvv)
		-> std::array<t_real, 8>
	{
		std::array<t_real, 8> arr;

		// horizontal
		arr[POP_SRC_Y] = t_real(-0.5) / dist_hsrc_mono * cm;

		arr[POP_MONO_X] = c_th_m * (t_real(0.5) / dist_mono_sample * cm
			- t_real(0.5) / dist_hsrc_mono * cm),
		arr[POP_MONO_Y] = s_th_m * (t_real(0.5) / dist_hsrc_mono * cm
			+ t_real(0.5) / dist_mono_sample * cm)
			- inv_curvh * cm;

		arr[POP_SAMPLE_X] = t_real(0.5) * s_th_s / dist_mono_sample * cm;
		arr[POP_SAMPLE_Y] = t_real(0.5) * c_th_s / dist_mono_sample * cm;

		// vertical
		arr[POP_SRC_Z] = t_real(-0.5) / (dist_vsrc_mono/cm * s_th_m);
		arr[POP_MONO_Z] = ( + t_real(0.5) / dist_vsrc_mono * cm
			+ t_real(0.5) / dist_mono_sample * cm ) / s_th_m
			- inv_curvv * cm;
		arr[POP_SAMPLE_Z] = t_real(-0.5) / (dist_mono_sample/cm * s_th_m);

		return arr;
	};

	t_mat T_mosaic_trafo = ublas::zero_matrix<t_real>(POP_NUM_COORDS, POP_NUM_POS);

	auto mono_mosaic_trafo = get_mosaic_trafo(
		pop.dist_vsrc_mono, pop.dist_hsrc_mono, pop.dist_mono_sample,
		s_th_m, c_th_m, s_th_s, c_th_s,
		inv_mono_curvh, inv_mono_curvv);
	T_mosaic_trafo(POP_MONO_H, POP_SRC_Y) = mono_mosaic_trafo[POP_SRC_Y]; // sign forgotten in paper
	T_mosaic_trafo(POP_MONO_H, POP_MONO_X) = mono_mosaic_trafo[POP_MONO_X];
	T_mosaic_trafo(POP_MONO_H, POP_MONO_Y) = mono_mosaic_trafo[POP_MONO_Y];
	T_mosaic_trafo(POP_MONO_H, POP_SAMPLE_X) = mono_mosaic_trafo[POP_SAMPLE_X];
	T_mosaic_trafo(POP_MONO_H, POP_SAMPLE_Y) = mono_mosaic_trafo[POP_SAMPLE_Y];
	T_mosaic_trafo(POP_MONO_V, POP_SRC_Z) = sign_z*mono_mosaic_trafo[POP_SRC_Z];
	T_mosaic_trafo(POP_MONO_V, POP_MONO_Z) = sign_z*mono_mosaic_trafo[POP_MONO_Z];
	T_mosaic_trafo(POP_MONO_V, POP_SAMPLE_Z) = sign_z*mono_mosaic_trafo[POP_SAMPLE_Z];

	auto ana_mosaic_trafo = get_mosaic_trafo(
		pop.dist_ana_det, pop.dist_ana_det, pop.dist_sample_ana,
		s_th_a, -c_th_a, s_th_s, -c_th_s,
		inv_ana_curvh, inv_ana_curvv);
	T_mosaic_trafo(POP_ANA_H, POP_DET_Y) = -ana_mosaic_trafo[get_ki_pos(POP_DET_Y)];
	T_mosaic_trafo(POP_ANA_H, POP_ANA_X) = ana_mosaic_trafo[get_ki_pos(POP_ANA_X)];
	T_mosaic_trafo(POP_ANA_H, POP_ANA_Y) = ana_mosaic_trafo[get_ki_pos(POP_ANA_Y)];
	T_mosaic_trafo(POP_ANA_H, POP_SAMPLE_X) = ana_mosaic_trafo[get_ki_pos(POP_SAMPLE_X)];
	T_mosaic_trafo(POP_ANA_H, POP_SAMPLE_Y) = ana_mosaic_trafo[get_ki_pos(POP_SAMPLE_Y)];
	T_mosaic_trafo(POP_ANA_V, POP_DET_Z) = sign_z*ana_mosaic_trafo[get_ki_pos(POP_DET_Z)];
	T_mosaic_trafo(POP_ANA_V, POP_ANA_Z) = sign_z*ana_mosaic_trafo[get_ki_pos(POP_ANA_Z)];
	T_mosaic_trafo(POP_ANA_V, POP_SAMPLE_Z) = sign_z*ana_mosaic_trafo[get_ki_pos(POP_SAMPLE_Z)];

	// D matrix to transform spatial to divergence variables, [pop75], Appendix 2
	auto get_geo_trafo = [](const length& dist_vsrc_mono,
		const length& dist_hsrc_mono, const length& dist_mono_sample,
		t_real s_th_m, t_real c_th_m, t_real s_th_s, t_real c_th_s)
		-> std::array<t_real, 11>
	{
		std::array<t_real, 11> arr;

		// POP_PREMONO_H
		arr[0] = t_real(-1) / dist_hsrc_mono * cm;   // POP_SRC_Y
		arr[1] = -c_th_m / dist_hsrc_mono * cm;      // POP_MONO_X
		arr[2] = s_th_m / dist_hsrc_mono * cm;       // POP_MONO_Y

		// POP_PRESAMPLE_H
		arr[3] = c_th_m / dist_mono_sample * cm;     // POP_MONO_X
		arr[4] = s_th_m / dist_mono_sample * cm;     // POP_MONO_Y
		arr[5] = s_th_s / dist_mono_sample * cm;     // POP_SAMPLE_X
		arr[6] = c_th_s / dist_mono_sample * cm;     // POP_SAMPLE_Y

		// POP_PREMONO_V
		arr[7] = t_real(-1) / dist_vsrc_mono * cm;   // POP_SRC_Z
		arr[8] = t_real(1) / dist_vsrc_mono * cm;    // POP_MONO_Z

		// POP_PRESAMPLE_V
		arr[9] = t_real(-1) / dist_mono_sample * cm; // POP_MONO_Z
		arr[10] = t_real(1) / dist_mono_sample * cm; // POP_SAMPLE_Z

		return arr;
	};

	t_mat D_geo_div_trafo = ublas::zero_matrix<t_real>(POP_NUM_COMPS, POP_NUM_POS);

	auto mono_geo_trafo = get_geo_trafo(pop.dist_vsrc_mono,
		pop.dist_hsrc_mono, pop.dist_mono_sample,
		s_th_m, c_th_m, s_th_s, c_th_s);
	D_geo_div_trafo(POP_PREMONO_H, POP_SRC_Y) = mono_geo_trafo[0];
	D_geo_div_trafo(POP_PREMONO_H, POP_MONO_X) = mono_geo_trafo[1];
	D_geo_div_trafo(POP_PREMONO_H, POP_MONO_Y) = mono_geo_trafo[2];
	D_geo_div_trafo(POP_PRESAMPLE_H, POP_MONO_X) = mono_geo_trafo[3];
	D_geo_div_trafo(POP_PRESAMPLE_H, POP_MONO_Y) = mono_geo_trafo[4];
	D_geo_div_trafo(POP_PRESAMPLE_H, POP_SAMPLE_X) = mono_geo_trafo[5];
	D_geo_div_trafo(POP_PRESAMPLE_H, POP_SAMPLE_Y) = mono_geo_trafo[6];
	D_geo_div_trafo(POP_PREMONO_V, POP_SRC_Z) = sign_z*mono_geo_trafo[7];
	D_geo_div_trafo(POP_PREMONO_V, POP_MONO_Z) = sign_z*mono_geo_trafo[8];
	D_geo_div_trafo(POP_PRESAMPLE_V, POP_MONO_Z) = sign_z*mono_geo_trafo[9];
	D_geo_div_trafo(POP_PRESAMPLE_V, POP_SAMPLE_Z) = sign_z*mono_geo_trafo[10];

	auto ana_geo_trafo = get_geo_trafo(pop.dist_ana_det,
		pop.dist_ana_det, pop.dist_sample_ana,
		s_th_a, -c_th_a, s_th_s, -c_th_s);
	D_geo_div_trafo(POP_POSTANA_H, POP_DET_Y) = -ana_geo_trafo[0];
	D_geo_div_trafo(POP_POSTANA_H, POP_ANA_X) = ana_geo_trafo[1];
	D_geo_div_trafo(POP_POSTANA_H, POP_ANA_Y) = ana_geo_trafo[2];
	D_geo_div_trafo(POP_POSTSAMPLE_H, POP_ANA_X) = ana_geo_trafo[3];
	D_geo_div_trafo(POP_POSTSAMPLE_H, POP_ANA_Y) = ana_geo_trafo[4];
	D_geo_div_trafo(POP_POSTSAMPLE_H, POP_SAMPLE_X) = ana_geo_trafo[5];
	D_geo_div_trafo(POP_POSTSAMPLE_H, POP_SAMPLE_Y) = ana_geo_trafo[6]; // forgotten in paper
	D_geo_div_trafo(POP_POSTANA_V, POP_DET_Z) = -sign_z*ana_geo_trafo[7];
	D_geo_div_trafo(POP_POSTANA_V, POP_ANA_Z) = -sign_z*ana_geo_trafo[8];
	D_geo_div_trafo(POP_POSTSAMPLE_V, POP_ANA_Z) = -sign_z*ana_geo_trafo[9];
	D_geo_div_trafo(POP_POSTSAMPLE_V, POP_SAMPLE_Z) = -sign_z*ana_geo_trafo[10];
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// covariance matrix calculation
	// --------------------------------------------------------------------
	// [pop75], equ. 20
	// [T] = 1/cm, [F] = 1/rad^2, [pop75], equ. 15
	t_mat K_geo = S_geo + tl::transform(F_mosaics, T_mosaic_trafo, true);
	t_mat Ki_geo;
	if(!tl::inverse(K_geo, Ki_geo))
	{
		res.bOk = false;
		res.strErr = "Matrix K cannot be inverted.";
		return res;
	}


	// [pop75], equ. 17
	t_mat Hi_div = tl::transform_inv(Ki_geo, D_geo_div_trafo, true);
	t_mat H_div;
	if(!tl::inverse(Hi_div, H_div))
	{
		res.bOk = false;
		res.strErr = "Matrix H^(-1) cannot be inverted.";
		return res;
	}

	// [pop75], equ. 18
	t_mat H_G_div = H_div + G_collis;
	t_mat H_Gi_div;
	if(!tl::inverse(H_G_div, H_Gi_div))
	{
		res.bOk = false;
		res.strErr = "Matrix H+G cannot be inverted.";
		return res;
	}

	//
	// [pop75], equ. 20, resolution matrix:
	// R = ((BA) ((                  H                      ) + G_collis)^(-1) (BA)^T)^(-1)
	// R = ((BA) ((D (        K              )^(-1) D^T)^(-1) + G_collis)^(-1) (BA)^T)^(-1)
	// R = ((BA) ((D (T^T F_mosaics T + S_geo)^(-1) D^T)^(-1) + G_collis)^(-1) (BA)^T)^(-1)
	//
	// compare to cn resolution matrix, [pop75], equ. 11:
	// R = ((BA) (                      H_cn                            )^(-1) (BA)^T)^(-1)
	// R = ((BA) (    C^T F_mosaics C                         + G_collis)^(-1) (BA)^T)^(-1)
	//  => C = T D^(-1) for zero instrument lengths
	//
	// [pop75], equ. 13a & 16, R0 intensity factor:
	// R0 ~ sqrt(|F_mosaics| |S_geo| / (|K| |(D S_geo^(-1) D^T)^(-1) + G_collis|))
	//
	// compare to cn R0 factor, [pop75], equ. 9 and equ. 5:
	// R0 ~ sqrt(|F_mosaics| / |H_cn|)
	//
	t_mat BA = ublas::prod(B_trafo_QE, A_div_kikf_trafo);
	t_mat cov = tl::transform_inv(H_Gi_div, BA, true);

	// include sample mosaic, see [zhe07], equs. 12-14
	t_real mos_h = pop.Q*pop.Q*angs*angs * pop.sample_mosaic*pop.sample_mosaic /rads/rads;
	t_real mos_v = pop.Q*pop.Q*angs*angs * sample_mosaic_v*sample_mosaic_v /rads/rads;
	//t_mat cov_nomosaic = cov;
	cov(1, 1) += mos_h;
	cov(2, 2) += mos_v;

	if(!tl::inverse(cov, res.reso))
	{
		res.bOk = false;
		res.strErr = "Covariance matrix cannot be inverted.";
		return res;
	}
	// -------------------------------------------------------------------------


	// --------------------------------------------------------------------
	// r0 intensity scaling factor and resolution volume calculation
	// --------------------------------------------------------------------
	res.reso *= sig2fwhm*sig2fwhm;  // convert back to sigmas
	res.reso_v = ublas::zero_vector<t_real>(4);
	res.reso_s = 0.;

	/*if(pop.dsample_sense < 0.)
	{
		// mirror Q_perp
		t_mat matMirror = tl::mirror_matrix<t_mat>(res.reso.size1(), 1);
		res.reso = tl::transform(res.reso, matMirror, true);
		res.reso_v[1] = -res.reso_v[1];
	}*/

	res.dResVol = tl::get_ellipsoid_volume(res.reso);
	res.dR0 = dmono_refl * dana_effic * dxsec * dmonitor;

	// include sample mosaic, see [zhe07], equs. 12-14 and cf. equs. 15 & 16
	// in [zhe07], equ. 16 the (det(...)/det(...))^(1/2) factor should be
	// (det(...)/(det(...)))^(-1/2), cf. volume normalisation in [zhe07], equ. 8
	// this just divides out the old factor (res_nomosaic) and multiplies the new (res)
	// and is already taken care of by the total normalisation further below
	//res.dR0 /= std::sqrt(std::abs(tl::determinant(cov) / tl::determinant(cov_nomosaic)));

	// --------------------------------------------------------------------
	// mono parts of the matrices until the monitor, see: [zhe07], p. 10
	// --------------------------------------------------------------------
	// mono - monitor distance and monitor dimensions
	length dist_mono_monitor = pop.dist_mono_monitor;
	length monitor_w = pop.monitor_w;
	length monitor_h = pop.monitor_h;
	length monitor_thick = pop.monitor_thick;

	// if the given monitor values are invalid, use the mono dimensions as fallback values
	if(tl::float_equal<t_real>(dist_mono_monitor / cm, 0.))
		dist_mono_monitor = pop.dist_mono_sample;
	if(tl::float_equal<t_real>(monitor_w / cm, 0.))
		monitor_w = pop.mono_w;
	if(tl::float_equal<t_real>(monitor_h / cm, 0.))
		monitor_h = pop.mono_h;
	if(tl::float_equal<t_real>(monitor_thick / cm, 0.))
		monitor_thick = pop.mono_thick;

	t_mat G_mono_collis = G_collis;
	G_mono_collis.resize(POP_PRESAMPLE_V+1, POP_PRESAMPLE_V+1, true);

	t_mat F_mono_mosaics = F_mosaics;
	F_mono_mosaics.resize(POP_MONO_V+1, POP_MONO_V+1, true);

	// these factors are used to approximate a uniform or circular variance by a gaussian variance
	t_mat SI_mono_geo = tl::zero_matrix(POP_SAMPLE_Z+1, POP_SAMPLE_Z+1);
	SI_mono_geo(POP_SRC_Y, POP_SRC_Y) = dMultSrc * pop.src_w*pop.src_w /cm/cm;
	SI_mono_geo(POP_SRC_Z, POP_SRC_Z) = dMultSrc * pop.src_h*pop.src_h /cm/cm;

	SI_mono_geo(POP_MONO_X, POP_MONO_X) = dMultMono * pop.mono_thick*pop.mono_thick /cm/cm;
	SI_mono_geo(POP_MONO_Y, POP_MONO_Y) = dMultMono * pop.mono_w*pop.mono_w /cm/cm;
	SI_mono_geo(POP_MONO_Z, POP_MONO_Z) = dMultMono * pop.mono_h*pop.mono_h /cm/cm;

	// monitor in sample position of calculation
	SI_mono_geo(POP_SAMPLE_X, POP_SAMPLE_X) = var_uniform * monitor_thick*monitor_thick /cm/cm;;
	SI_mono_geo(POP_SAMPLE_Y, POP_SAMPLE_Y) = dMultMonitor * monitor_w*monitor_w /cm/cm;
	SI_mono_geo(POP_SAMPLE_Z, POP_SAMPLE_Z) = dMultMonitor * monitor_h*monitor_h /cm/cm;

	SI_mono_geo *= sig2fwhm*sig2fwhm;  // convert gaussian variance to gaussian fwhm

	t_mat S_mono_geo;
	if(!tl::inverse_diag(SI_mono_geo, S_mono_geo))
	{
		res.bOk = false;
		res.strErr = "S_mono matrix cannot be inverted.";
		return res;
	}

	t_real sign_z_mon = 1.; //sign_z;
	t_mat T_mono_mosaic_trafo = ublas::zero_matrix<t_real>(POP_MONO_V+1, POP_SAMPLE_Z+1);
	auto mon_mosaic_trafo = get_mosaic_trafo(
		pop.dist_vsrc_mono, pop.dist_hsrc_mono,
		dist_mono_monitor,  // monitor in sample position of calculation
		s_th_m, c_th_m, 0., 1.,
		inv_mono_curvh, inv_mono_curvv);
	T_mono_mosaic_trafo(POP_MONO_H, POP_SRC_Y) = mon_mosaic_trafo[POP_SRC_Y]; // sign forgotten in paper
	T_mono_mosaic_trafo(POP_MONO_H, POP_MONO_X) = mon_mosaic_trafo[POP_MONO_X];
	T_mono_mosaic_trafo(POP_MONO_H, POP_MONO_Y) = mon_mosaic_trafo[POP_MONO_Y];
	T_mono_mosaic_trafo(POP_MONO_H, POP_SAMPLE_X) = mon_mosaic_trafo[POP_SAMPLE_X];
	T_mono_mosaic_trafo(POP_MONO_H, POP_SAMPLE_Y) = mon_mosaic_trafo[POP_SAMPLE_Y];
	T_mono_mosaic_trafo(POP_MONO_V, POP_SRC_Z) = sign_z_mon*mon_mosaic_trafo[POP_SRC_Z];
	T_mono_mosaic_trafo(POP_MONO_V, POP_MONO_Z) = sign_z_mon*mon_mosaic_trafo[POP_MONO_Z];
	T_mono_mosaic_trafo(POP_MONO_V, POP_SAMPLE_Z) = sign_z_mon*mon_mosaic_trafo[POP_SAMPLE_Z];

	t_mat D_mono_geo_div_trafo = ublas::zero_matrix<t_real>(POP_PRESAMPLE_V+1, POP_SAMPLE_Z+1);
	auto mon_geo_trafo = get_geo_trafo(pop.dist_vsrc_mono,
		pop.dist_hsrc_mono, dist_mono_monitor,  // monitor in sample position of calculation with "scattering" angle 0
		s_th_m, c_th_m, 0., 1.);
	D_mono_geo_div_trafo(POP_PREMONO_H, POP_SRC_Y) = mon_geo_trafo[0];
	D_mono_geo_div_trafo(POP_PREMONO_H, POP_MONO_X) = mon_geo_trafo[1];
	D_mono_geo_div_trafo(POP_PREMONO_H, POP_MONO_Y) = mon_geo_trafo[2];
	D_mono_geo_div_trafo(POP_PRESAMPLE_H, POP_MONO_X) = mon_geo_trafo[3];
	D_mono_geo_div_trafo(POP_PRESAMPLE_H, POP_MONO_Y) = mon_geo_trafo[4];
	D_mono_geo_div_trafo(POP_PRESAMPLE_H, POP_SAMPLE_X) = mon_geo_trafo[5];
	D_mono_geo_div_trafo(POP_PRESAMPLE_H, POP_SAMPLE_Y) = mon_geo_trafo[6];
	D_mono_geo_div_trafo(POP_PREMONO_V, POP_SRC_Z) = sign_z_mon*mon_geo_trafo[7];
	D_mono_geo_div_trafo(POP_PREMONO_V, POP_MONO_Z) = sign_z_mon*mon_geo_trafo[8];
	D_mono_geo_div_trafo(POP_PRESAMPLE_V, POP_MONO_Z) = sign_z_mon*mon_geo_trafo[9];
	D_mono_geo_div_trafo(POP_PRESAMPLE_V, POP_SAMPLE_Z) = sign_z_mon*mon_geo_trafo[10];

	t_mat K_mono_geo = S_mono_geo + tl::transform(F_mono_mosaics, T_mono_mosaic_trafo, true);
	t_mat Ki_mono_geo;
	if(!tl::inverse(K_mono_geo, Ki_mono_geo))
	{
		res.bOk = false;
		res.strErr = "Matrix K_mono cannot be inverted.";
		return res;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// R0 calculation methods
	// --------------------------------------------------------------------
	t_real dDetF = tl::determinant(F_mosaics);

	if(pop.flags & CALC_ALT_R0)
	{
		// resolution volume, [pop75], equ. 13a & 16
		// [D] = 1/cm, [SI] = cm^2
		t_mat DSiDt = tl::transform_inv(SI_geo, D_geo_div_trafo, true);
		t_mat DSiDti;
		if(!tl::inverse(DSiDt, DSiDti))
		{
			res.bOk = false;
			res.strErr = "R0 factor cannot be calculated.";
			return res;
		}
		DSiDti += G_collis;

		t_real dDetS = tl::determinant(S_geo);
		t_real dDetK = tl::determinant(K_geo);
		t_real dDetDSiDti = tl::determinant(DSiDti);

		// [pop75], equs. 13a & 16
		res.dR0 *= std::sqrt(dDetS*dDetF / (dDetK*dDetDSiDti));

		if(pop.flags & CALC_MON)
		{
			// mono part
			t_mat DSiDt_mono = tl::transform_inv(SI_mono_geo, D_mono_geo_div_trafo, true);
			t_mat DSiDti_mono;
			if(!tl::inverse(DSiDt_mono, DSiDti_mono))
			{
				res.bOk = false;
				res.strErr = "R0_mono factor cannot be calculated.";
				return res;
			}
			DSiDti_mono += G_mono_collis;

			t_real dDetS_mono = tl::determinant(S_mono_geo);
			t_real dDetF_mono = tl::determinant(F_mono_mosaics);
			t_real dDetK_mono = tl::determinant(K_mono_geo);
			t_real dDetDSiDti_mono = tl::determinant(DSiDti_mono);

			res.dR0 /= std::sqrt(std::abs(dDetS_mono*dDetF_mono / (dDetK_mono*dDetDSiDti_mono)));
			res.dR0 /= dmono_refl;  // removes all monochromator intensity factors
		}
	}
	else
	{
		// more general calculation from [zhe07], p. 10, equ. 8
		t_real dDetHG = tl::determinant(H_G_div);
		res.dR0 *= std::sqrt(dDetF / dDetHG);

		if(pop.flags & CALC_MON)
		{
			t_mat Hi_mono_div = tl::transform_inv(Ki_mono_geo, D_mono_geo_div_trafo, true);
			t_mat HG_mono_div;
			if(!tl::inverse(Hi_mono_div, HG_mono_div))
			{
				res.bOk = false;
				res.strErr = "Matrix H_mono^(-1) cannot be inverted.";
				return res;
			}
			HG_mono_div += G_mono_collis;

			// mono part, [zhe07], p. 10, equ. 10
			res.dR0 /= std::sqrt(std::abs(tl::determinant(F_mono_mosaics) / tl::determinant(HG_mono_div)));
			res.dR0 /= dmono_refl;  // removes all monochromator intensity factors
		}
	}

	res.dR0 = std::abs(res.dR0);

	if(pop.flags & NORM_TO_SAMPLE)
	{
		// normalise R0 to sample volume
		volume V_sample = pop.bSampleCub
			? pop.sample_w_perpq * pop.sample_w_q * pop.sample_h
			: pi * 0.5*pop.sample_w_perpq * 0.5*pop.sample_w_q * pop.sample_h;
		res.dR0 /= V_sample*V_sample / units::pow<6>(cm);
	}


	if(pop.flags & NORM_TO_RESVOL)
	{
		// rest of the prefactors, equ. 1 in [pop75], together with the mono and and ana reflectivities
		// (defining the resolution volume), these give the same correction as in [mit84] equ. A.57
		// NOTE: these factors are not needed if the normalisation of the 4d gaussian distribution
		// is already taken care of in the MC step by the employed std::normal_distribution function
		// res.dR0 *= std::sqrt(std::abs(tl::determinant(res.reso))) / (2.*pi*2.*pi);
		// except for the (unimportant) prefactors this is the same as dividing by the resolution volume
		res.dR0 /= res.dResVol * pi * t_real(3.);
	}
	// --------------------------------------------------------------------

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
