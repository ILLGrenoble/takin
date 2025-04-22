/**
 * cooper-nathans calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2
 *
 * @desc This is a reimplementation in C++ of the file rc_cnmat.m of the
 *		rescal5 package by Zinkin, McMorrow, Tennant, Farhi, and Wildes (ca. 1995-2007):
 *		http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/
 * @desc see:
 *		[cn67] M. J. Cooper and R. Nathans, Acta Cryst. 23, 357 (1967), doi: 10.1107/S0365110X67002816
 *		[ch73] N. J. Chesser and J. D. Axe, Acta Cryst. A 29, 160 (1973), doi: 10.1107/S0567739473000422
 *		[mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
 *		[pop75] M. Popovici, Acta Cryst. A 31, 507 (1975), doi: 10.1107/S0567739475001088
 *		[zhe07] A. Zheludev, ResLib 3.4 manual (2007), https://ethz.ch/content/dam/ethz/special-interest/phys/solid-state-physics/neutron-scattering-and-magnetism-dam/images/research/manual.pdf
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

#include "cn.h"
#include "ellipse.h"

#include "tlibs/math/geo.h"
#include "tlibs/math/math.h"
#include "tlibs/log/log.h"

#include <string>
#include <future>
#include <iostream>


typedef t_real_reso t_real;
typedef ublas::matrix<t_real> t_mat;
typedef ublas::vector<t_real> t_vec;

using angle = tl::t_angle_si<t_real>;
using wavenumber = tl::t_wavenumber_si<t_real>;
using energy = tl::t_energy_si<t_real>;
using length = tl::t_length_si<t_real>;

static const auto angs = tl::get_one_angstrom<t_real>();
static const auto rads = tl::get_one_radian<t_real>();
static const auto meV = tl::get_one_meV<t_real>();
static const auto sec = tl::get_one_second<t_real>();
static const auto mn = tl::get_m_n<t_real>();
static const auto hbar = tl::get_hbar<t_real>();
static const t_real pi = tl::get_pi<t_real>();
static const t_real sig2fwhm = tl::get_SIGMA2FWHM<t_real>();


// -----------------------------------------------------------------------------
/**
 * scattering factors
 */
std::tuple<t_real, t_real, t_real, t_real> get_scatter_factors(
	std::size_t flags,
	const angle& thetam, const wavenumber& ki,
	const angle& thetaa, const wavenumber& kf)
{
	t_real dmono = t_real(1);
	t_real dana = t_real(1);
	t_real dSqwToXSec = t_real(1);
	t_real dmonitor = t_real(1);

	if(flags & CALC_KI3)
		dmono *= tl::ana_effic_factor(ki, units::abs(thetam));
	if(flags & CALC_KF3)
		dana *= tl::ana_effic_factor(kf, units::abs(thetaa));
	if(flags & CALC_KFKI)
		dSqwToXSec *= kf/ki;  // kf/ki factor, see Shirane, equ. (2.7)
	if(flags & CALC_MONKI)
		dmonitor *= ki*angs;  // monitor 1/ki factor, see [zhe07], p. 10

	return std::make_tuple(dmono, dana, dSqwToXSec, dmonitor);
}


/**
 * transformation matrix -> [mit84], equ. A.15 and [pop75], Appendix 1
 *        dki part                dkf part
 * (  Ti11   Ti12      0  |   Tf11   Tf12      0 )   ( dki_x )   ( dQ_x  )
 * (  Ti12   Ti22      0  |   Tf12   Tf22      0 )   ( dki_y )   ( dQ_y  )
 * (     0      0      1  |      0      0     -1 ) * ( dki_z ) = ( dQ_z  )
 * ( 2ki*c      0      0  | -2kf*c      0      0 )   ( dkf_x )   ( dE    )
 * (     1      0      0  |      0      0      0 )   ( dkf_y )   ( dki_x )
 * (     0      0      1  |      0      0      0 )   ( dkf_z )   ( dki_z )
 */
t_mat get_trafo_dkidkf_dQdE(const angle& ki_Q, const angle& kf_Q,
	const wavenumber& ki, const wavenumber& kf)
{
	t_mat Ti = tl::rotation_matrix_2d(ki_Q/rads);
	t_mat Tf = -tl::rotation_matrix_2d(kf_Q/rads);

	// dQ_{x,y} = dki_{x,y} - dkf_{x,y}
	t_mat U = ublas::zero_matrix<t_real>(6, 6);
	tl::submatrix_copy(U, Ti, 0, 0);
	tl::submatrix_copy(U, Tf, 0, 3);

	// dQ_z = dki_z - dkf_z
	U(2 /*dQ_z*/, 2 /*dki_z*/) = 1.;
	U(2 /*dQ_z*/, 5 /*dkf_z*/) = -1.;

	//  E ~ ki^2 - kf^2
	// dE ~ 2ki*dki - 2kf*dkf
	U(3 /*dE*/, 0 /*dki_x*/) = +t_real(2)*ki * tl::get_KSQ2E<t_real>() * angs;
	U(3 /*dE*/, 3 /*dkf_x*/) = -t_real(2)*kf * tl::get_KSQ2E<t_real>() * angs;

	// simply copy the same variables
	U(4 /*dki_x*/, 0 /*dki_x*/) = 1.;
	U(5 /*dki_z*/, 2 /*dki_z*/) = 1.;

	return U;
}
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
/**
 * R0 factor from formula (2) in [ch73]
 */
t_real chess_R0(bool norm_to_ki_vol,
	wavenumber ki, wavenumber kf,
	const t_mat& M, const t_mat& V,
	angle theta_m, angle theta_a, angle twotheta_s,
	angle mos_m, angle mos_a, angle coll_pre_mono_v, angle coll_post_ana_v,
	t_real refl_m, t_real refl_a)
{
	auto R0_P = [](angle theta, angle coll, angle mosaic) -> t_real
	{
		t_real tS = units::sin(theta);
		return std::sqrt(t_real(2)*pi) / rads *
			tl::my_units_sqrt<angle>(t_real(1) / (
				t_real(1)/(coll*coll) + t_real(1)/(t_real(4)*mosaic*mosaic*tS*tS)));
	};

	auto R0_N = [](angle theta, angle mosaic, t_real refl) -> t_real
	{
		t_real tS = units::sin(theta);
		return (refl / (t_real(2)*mosaic * tS)) / std::sqrt(t_real(2)*pi) * rads;
	};

	t_real s_tt = units::sin(twotheta_s);
	t_real R0 = mn/hbar / (ki*ki * kf*kf*kf * s_tt) / angs/angs/angs/sec;

	// kf volumne
	R0 *= R0_P(theta_a, coll_post_ana_v, mos_a) * R0_N(theta_a, mos_a, refl_a);

	// ki volume
	if(!norm_to_ki_vol)
		R0 *= R0_P(theta_m, coll_pre_mono_v, mos_m) * R0_N(theta_m, mos_m, refl_m);

	// R0x see: [mit84], equ. A.43
	R0 /= std::sqrt(
		  M(1,1)*V(1,4)*V(1,4)    + M(3,3)*V(3,4)*V(3,4)
		+ M(4,4)*V(4,4)*V(4,4)    + M(0,0)*V(0,4)*V(0,4)
		+ M(0,1)*V(0,4)*V(1,4)*2. + M(3,4)*V(4,4)*V(3,4)*2.);
	// TODO: R0z
	R0 *= t_real(2)*pi;

	return R0;
}
// -----------------------------------------------------------------------------



ResoResults calc_cn(const CNParams& cn)
{
	ResoResults res;

	res.Q_avg.resize(4);
	res.Q_avg[0] = cn.Q * angs;
	res.Q_avg[1] = 0.;
	res.Q_avg[2] = 0.;
	res.Q_avg[3] = cn.E / meV;

	// if the user moved the scattering angle to the other
	// side of the scattering sense indicated by the flag
	t_real manually_changed_sense = t_real(1);
	if(cn.twotheta/rads < t_real(0))
		manually_changed_sense = t_real(-1);

	angle thetaa = cn.thetaa * cn.dana_sense;
	angle thetam = cn.thetam * cn.dmono_sense;
	angle ki_Q = cn.angle_ki_Q * cn.dsample_sense * manually_changed_sense;
	angle kf_Q = cn.angle_kf_Q * cn.dsample_sense * manually_changed_sense;

	// transformation matrix U and its inverse V
	t_mat U_trafo_QE = get_trafo_dkidkf_dQdE(ki_Q, kf_Q, cn.ki, cn.kf);

	// V matrix -> [mit84], equ. A.16
	t_mat V_inv_trafo_QE(6, 6);
	if(!tl::inverse(U_trafo_QE, V_inv_trafo_QE))
	{
		res.bOk = false;
		res.strErr = "Transformation matrix cannot be inverted.";
		return res;
	}
	// -------------------------------------------------------------------------


	const auto tupScFact = get_scatter_factors(cn.flags, cn.thetam, cn.ki, cn.thetaa, cn.kf);

	t_real dmono_refl = cn.dmono_refl * std::get<0>(tupScFact);
	t_real dana_effic = cn.dana_effic * std::get<1>(tupScFact);
	if(cn.mono_refl_curve) dmono_refl *= (*cn.mono_refl_curve)(cn.ki);
	if(cn.ana_effic_curve) dana_effic *= (*cn.ana_effic_curve)(cn.kf);
	t_real dxsec = std::get<2>(tupScFact);
	t_real dmonitor = std::get<3>(tupScFact);


	// if no vertical mosaic is given, use the horizontal one
	angle mono_mosaic_v = cn.mono_mosaic_v;
	angle ana_mosaic_v = cn.ana_mosaic_v;
	angle sample_mosaic_v = cn.sample_mosaic_v;
	if(tl::float_equal<t_real>(mono_mosaic_v/rads, 0.), 0.)
		mono_mosaic_v = cn.mono_mosaic;
	if(tl::float_equal<t_real>(ana_mosaic_v/rads, 0.), 0.)
		ana_mosaic_v = cn.ana_mosaic;
	if(tl::float_equal<t_real>(sample_mosaic_v/rads, 0.), 0.)
		sample_mosaic_v = cn.sample_mosaic;


	// -------------------------------------------------------------------------
	// resolution matrix, [mit84], equ. A.5
	auto calc_mono_ana_res =
		[](angle theta, wavenumber k,
		angle mosaic, angle mosaic_v,
		angle coll1, angle coll2,
		angle coll1_v, angle coll2_v) -> t_mat
	{
		t_real t_th = units::tan(theta);

		// horizontal part
		t_vec vecMos(2);
		vecMos[0] = t_th;
		vecMos[1] = 1.;
		vecMos /= k*angs * mosaic/rads;

		t_vec vecColl1(2);
		vecColl1[0] = 2.*t_th;
		vecColl1[1] = 1.;
		vecColl1 /= (k*angs * coll1/rads);

		t_vec vecColl2(2);
		vecColl2[0] = 0;
		vecColl2[1] = 1.;
		vecColl2 /= (k*angs * coll2/rads);

		t_mat matHori = ublas::outer_prod(vecMos, vecMos) +
			ublas::outer_prod(vecColl1, vecColl1) +
			ublas::outer_prod(vecColl2, vecColl2);

		t_real s_th = units::sin(theta);

		// vertical part, [mit84], equ. A.9 & A.13
		t_real dVert = t_real(1)/(k*k * angs*angs) * rads*rads *
		(
			t_real(1) / (coll2_v * coll2_v) +
			t_real(1) / (
				(t_real(2)*s_th * mosaic_v) *
				(t_real(2)*s_th * mosaic_v) +
				coll1_v * coll1_v)
		);

		t_mat matFull = ublas::zero_matrix<t_real>(3, 3);
		tl::submatrix_copy(matFull, matHori, 0, 0);
		matFull(2, 2) = dVert;

		return matFull;
	};

	std::launch lpol = /*std::launch::deferred |*/ std::launch::async;
	std::future<t_mat> futMono
		= std::async(lpol, calc_mono_ana_res,
			thetam, cn.ki,
			cn.mono_mosaic, mono_mosaic_v,
			cn.coll_h_pre_mono, cn.coll_h_pre_sample,
			cn.coll_v_pre_mono, cn.coll_v_pre_sample);
	std::future<t_mat> futAna
		= std::async(lpol, calc_mono_ana_res,
			-thetaa, cn.kf,
			cn.ana_mosaic, ana_mosaic_v,
			cn.coll_h_post_ana, cn.coll_h_post_sample,
			cn.coll_v_post_ana, cn.coll_v_post_sample);

	t_mat matMono = futMono.get();
	t_mat matAna = futAna.get();

	t_mat M = ublas::zero_matrix<t_real>(6, 6);
	tl::submatrix_copy(M, matMono, 0, 0);
	tl::submatrix_copy(M, matAna, 3, 3);
	// -------------------------------------------------------------------------


	// transform reso matrix, see [mit84], p. 158
	t_mat M_trafo = tl::transform(M, V_inv_trafo_QE, 1);

	// integrate components, see [mit84], p. 159
	M_trafo = quadric_proj(M_trafo, 5);
	M_trafo = quadric_proj(M_trafo, 4);

	// add horizontal sample mosaic
	const t_real mos_Q_sq =
		(cn.sample_mosaic/rads * cn.Q*angs) *
		(cn.sample_mosaic/rads * cn.Q*angs);
	t_vec vec1 = tl::get_column<t_vec>(M_trafo, 1);
	//const t_real M_vert = M_trafo(2, 2);
	M_trafo -= ublas::outer_prod(vec1, vec1) / (1./mos_Q_sq + M_trafo(1, 1));
	//M_trafo(2, 2) = M_vert;

	// add vertical sample mosaic
	const t_real mos_v_Q_sq =
		(sample_mosaic_v/rads * cn.Q*angs) *
		(sample_mosaic_v/rads * cn.Q*angs);
	t_vec vec2 = tl::get_column<t_vec>(M_trafo, 2);
	M_trafo -= ublas::outer_prod(vec2, vec2) / (1./mos_v_Q_sq + M_trafo(2, 2));

	res.reso = M_trafo;
	res.reso *= sig2fwhm*sig2fwhm;
	res.reso_v = ublas::zero_vector<t_real>(4);
	res.reso_s = 0.;

	/*if(cn.dsample_sense < 0.)
	{
		// mirror Q_perp
		t_mat matMirror = tl::mirror_matrix<t_mat>(res.reso.size1(), 1);
		res.reso = tl::transform(res.reso, matMirror, true);
		res.reso_v[1] = -res.reso_v[1];
	}*/

	// -------------------------------------------------------------------------

	bool use_monitor = (cn.flags & CALC_MON) != 0;

	res.dResVol = tl::get_ellipsoid_volume(res.reso);
	res.dR0 = dana_effic * dxsec * dmonitor;
	if(!use_monitor)
		res.dR0 *= dmono_refl;

	t_real s_th_m = units::sin(thetam);
	t_real s_th_a = units::sin(thetaa);
	t_real t_th_m = units::tan(thetam);

	if(use_monitor)
	{
		// [mit84], equ. B.10/B.11 (similar to equ. B.13 with coll. idx 5 -> idx 1) and equ. B.3
		res.dR0 *=
			  t_real(1)/std::sqrt((M(1, 1) + M(4, 4)) * (M(0, 0) + M(3, 3)) - std::pow(M(0, 1) + M(3, 4), t_real(2)))
			* t_real(1)/std::sqrt(M(2, 2) + M(5, 5))
			* t_real(1)/std::sqrt(t_real(1) + std::pow(t_real(2)*s_th_a*ana_mosaic_v/cn.coll_v_post_ana, t_real(2)))
			* rads * units::sqrt(
				  t_real(1)/(cn.coll_v_pre_sample*cn.coll_v_pre_sample)
				+ t_real(1)/(cn.coll_v_pre_mono*cn.coll_v_pre_mono + t_real(4)*s_th_m*s_th_m*mono_mosaic_v*mono_mosaic_v))
			* rads*rads * units::sqrt(
				  t_real(1)/(cn.coll_h_pre_mono*cn.coll_h_pre_mono*cn.coll_h_pre_sample*cn.coll_h_pre_sample)
				+ t_real(1)/(t_real(4)*cn.coll_h_pre_mono*cn.coll_h_pre_mono*cn.mono_mosaic*cn.mono_mosaic)
				+ t_real(1)/(t_real(4)*cn.coll_h_pre_sample*cn.coll_h_pre_sample*cn.mono_mosaic*cn.mono_mosaic));
		res.dR0 /= cn.ki*cn.ki*cn.ki*cn.ki * angs*angs*angs*angs * t_th_m;
		//res.dR0 *= hbar/mn /angs/angs*sec / std::sqrt(2.*pi);
	}
	else
	{
		// [mit84], equ. B.10
		res.dR0 *=
			  t_real(1)/std::sqrt((M(1, 1) + M(4, 4)) * (M(0, 0) + M(3, 3)) - std::pow(M(0, 1) + M(3, 4), 2.))
			* t_real(1)/std::sqrt(M(2, 2) + M(5, 5))
			* t_real(1)/std::sqrt(t_real(1) + std::pow(t_real(2)*s_th_a*ana_mosaic_v/cn.coll_v_post_ana, t_real(2)))
			* t_real(1)/std::sqrt(t_real(1) + std::pow(t_real(2)*s_th_m*mono_mosaic_v/cn.coll_v_pre_mono, t_real(2)));
		res.dR0 /= cn.ki*cn.ki*cn.ki * angs*angs*angs;
		//res.dR0 *= 2.*pi*hbar/mn /angs/angs*sec;
	}

	/*res.dR0 *= chess_R0(use_monitor,
		cn.ki, cn.kf, M, V_inv_trafo_QE,
		thetam, thetaa, cn.twotheta,
		cn.mono_mosaic, cn.ana_mosaic,
		cn.coll_v_pre_mono, cn.coll_v_post_ana,
		dmono_refl, dana_effic);*/

	if(cn.flags & NORM_TO_RESVOL)
		res.dR0 /= res.dResVol * pi * t_real(3.);

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
