/**
 * popovici calculation of the cooper-nathans resolution
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2022
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
#include "ellipse.h"

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
using inv_length = tl::t_length_inverse_si<t_real>;

static const auto angs = tl::get_one_angstrom<t_real>();
static const auto rads = tl::get_one_radian<t_real>();
static const auto meV = tl::get_one_meV<t_real>();
static const auto cm = tl::get_one_centimeter<t_real>();
static const t_real sig2fwhm = tl::get_SIGMA2FWHM<t_real>();


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


ResoResults calc_pop_cn(const CNParams& pop)
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

	angle thetaa = pop.thetaa * pop.dana_sense;
	angle thetam = pop.thetam * pop.dmono_sense;
	angle ki_Q = pop.angle_ki_Q * pop.dsample_sense * manually_changed_sense;
	angle kf_Q = pop.angle_kf_Q * pop.dsample_sense * manually_changed_sense;

	// B matrix, [pop75], Appendix 1 -> U matrix in CN
	t_mat B_trafo_QE = get_trafo_dkidkf_dQdE(ki_Q, kf_Q, pop.ki, pop.kf);
	B_trafo_QE.resize(4, POP_NUM_KIKF, true);

	angle coll_h_pre_mono = pop.coll_h_pre_mono;
	angle coll_v_pre_mono = pop.coll_v_pre_mono;

	/*if(pop.bGuide)
	{
		length lam = tl::k2lam(pop.ki);
		coll_h_pre_mono = lam*(pop.guide_div_h/angs);
		coll_v_pre_mono = lam*(pop.guide_div_v/angs);
	}*/

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


	const auto tupScFact = get_scatter_factors(pop.flags, pop.thetam, pop.ki, pop.thetaa, pop.kf);

	t_real dmono_refl = pop.dmono_refl * std::get<0>(tupScFact);
	t_real dana_effic = pop.dana_effic * std::get<1>(tupScFact);
	if(pop.mono_refl_curve) dmono_refl *= (*pop.mono_refl_curve)(pop.ki);
	if(pop.ana_effic_curve) dana_effic *= (*pop.ana_effic_curve)(pop.kf);
	t_real dxsec = std::get<2>(tupScFact);
	t_real dmonitor = std::get<3>(tupScFact);


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
	const t_real sign_z = -1.;  // to compare with the calculations by F. Bourdarot

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
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// trafo matrix
	// --------------------------------------------------------------------
	// C matrix, [pop75], Appendix 1
	t_mat C_trafo = ublas::zero_matrix<t_real>(POP_NUM_COORDS, POP_NUM_COMPS);
	C_trafo(POP_MONO_H, POP_PRESAMPLE_H) = 0.5;
	C_trafo(POP_MONO_H, POP_PREMONO_H) = 0.5;
	C_trafo(POP_MONO_V, POP_PREMONO_V) = t_real(0.5)/s_th_m;
	C_trafo(POP_MONO_V, POP_PRESAMPLE_V) = t_real(-0.5)/s_th_m;  // typo in paper
	C_trafo(POP_ANA_H, POP_POSTANA_H) = 0.5;
	C_trafo(POP_ANA_H, POP_POSTSAMPLE_H) = 0.5;
	C_trafo(POP_ANA_V, POP_POSTSAMPLE_V) = t_real(0.5)/s_th_a;
	C_trafo(POP_ANA_V, POP_POSTANA_V) = t_real(-0.5)/s_th_a;
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// covariance matrix calculation
	// --------------------------------------------------------------------
	// [pop75], equ. 8
	// [T] = 1/cm, [F] = 1/rad^2, [pop75], equ. 8
	t_mat H = G_collis + tl::transform(F_mosaics, C_trafo, true);
	t_mat H_inv;
	if(!tl::inverse(H, H_inv))
	{
		res.bOk = false;
		res.strErr = "Matrix H cannot be inverted.";
		return res;
	}

	//
	// [pop75], equ. 11, resolution matrix:
	// R = ((BA) (              H           )^(-1) (BA)^T)^(-1)
	// R = ((BA) (C^T F_mosaics C + G_collis)^(-1) (BA)^T)^(-1)
	//
	// cn R0 factor, [pop75], equ. 9 and equ. 5:
	// R0 ~ sqrt(|F_mosaics| / |H|)
	//
	t_mat BA = ublas::prod(B_trafo_QE, A_div_kikf_trafo);
	t_mat cov = tl::transform_inv(H_inv, BA, true);

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
	res.reso *= sig2fwhm*sig2fwhm;
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
	// this just corrects the resolution volume, which is already normalised in the MC step
	//res.dR0 /= std::sqrt(std::abs(tl::determinant(cov) / tl::determinant(cov_nomosaic)));


	// --------------------------------------------------------------------
	// mono parts of the matrices, see: [zhe07], p. 10
	// --------------------------------------------------------------------
	t_mat G_mono_collis = G_collis;
	G_mono_collis.resize(POP_PRESAMPLE_V+1, POP_PRESAMPLE_V+1, true);

	t_mat F_mono_mosaics = F_mosaics;
	F_mono_mosaics.resize(POP_MONO_V+1, POP_MONO_V+1, true);

	t_mat C_mono_trafo = C_trafo;
	C_mono_trafo.resize(POP_MONO_V+1, POP_PRESAMPLE_V+1, true);

	t_mat H_mono = G_mono_collis + tl::transform(F_mono_mosaics, C_mono_trafo, true);
	//t_mat H_mono = H;
	//H_mono.resize(POP_PRESAMPLE_V+1, POP_PRESAMPLE_V+1, true);
	// --------------------------------------------------------------------


	const t_real pi = tl::get_pi<t_real>();

	// R0 calculation methods
	// cancels out in [pop75] equ. 5 and equ. 9
	//t_real dDetG = tl::determinant(G_collis);
	t_real dDetF = tl::determinant(F_mosaics);
	t_real dDetH = tl::determinant(H);

	// [pop75], equ. 9 and equ. 5 (and [zhe07], equ. 7)
	res.dR0 *= t_real((2.*pi)*(2.*pi)*(2.*pi)*(2.*pi));
	res.dR0 *= std::sqrt(dDetF / dDetH);
	res.dR0 /= t_real(8.*pi*8.*pi) * s_th_m * s_th_a;

	if(pop.flags & CALC_MON)
	{
		// cancels out in [pop75] equ. 5 and equ. 9
		//t_real dDetG_mono = tl::determinant(G_mono_collis);

		// [zhe07], equ. 9
		t_real dDetF_mono = tl::determinant(F_mono_mosaics);
		t_real dDetH_mono = tl::determinant(H_mono);

		res.dR0 /= std::sqrt(dDetF_mono / dDetH_mono);
		res.dR0 *= t_real(2.)/pi * s_th_m / dmono_refl;
	}

	res.dR0 = std::abs(res.dR0);

	// rest of the prefactors, equ. 1 in [pop75], together with the mono and and ana reflectivities
	// (defining the resolution volume) these give the same correction as in [mit84] equ. A.57
	// NOTE: these factors are not needed, because the normalisation of the 4d gaussian distribution
	// is already taken care of in the MC step by the employed std::normal_distribution function
	//res.dR0 *= std::sqrt(std::abs(tl::determinant(res.reso))) / (2.*pi*2.*pi);
	// except for the (unimportant) prefactors this is the same as dividing by the resolution volume
	//res.dR0 /= res.dResVol * pi * t_real(3.);
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
