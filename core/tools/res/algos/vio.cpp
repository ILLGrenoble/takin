/**
 * implementation of Violini's TOF reso algorithm
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2016
 * @license GPLv2
 *
 * @desc for used algo, see: [vio14] N. Violini et al., NIM A 736 (2014) pp. 31-39, doi: 10.1016/j.nima.2013.10.042
 * @desc results checked with: [ehl11] G. Ehlers et al., http://arxiv.org/abs/1109.1482 (2011)
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

#include "vio.h"
#include "../ellipse.h"

#include "tlibs/math/linalg.h"
#include "tlibs/math/geo.h"
#include "tlibs/math/math.h"
#include "tlibs/log/log.h"

#include <string>
#include <iostream>
#include <numeric>


typedef t_real_reso t_real;
typedef ublas::matrix<t_real> t_mat;
typedef ublas::vector<t_real> t_vec;

using angle = tl::t_angle_si<t_real>;
using wavenumber = tl::t_wavenumber_si<t_real>;
using velocity = tl::t_velocity_si<t_real>;
using t_time = tl::t_time_si<t_real>;
using energy = tl::t_energy_si<t_real>;
using length = tl::t_length_si<t_real>;
using mass = tl::t_mass_si<t_real>;

static const auto rads = tl::get_one_radian<t_real>();
static const length angs = tl::get_one_angstrom<t_real>();
static const energy meV = tl::get_one_meV<t_real>();
static const t_time sec = tl::get_one_second<t_real>();
static const length meter = tl::get_one_meter<t_real>();


ResoResults calc_vio(const VioParams& params)
{
	ResoResults res;
	res.Q_avg.resize(4);

	const angle& tt = params.twotheta;
	const wavenumber &ki = params.ki, &kf = params.kf;
	const energy E = tl::get_energy_transfer(ki, kf);
	const wavenumber Q = tl::get_sample_Q(ki, kf, tt);

	res.Q_avg[0] = Q * angs;
	res.Q_avg[1] = 0.;
	res.Q_avg[2] = 0.;
	res.Q_avg[3] = E / meV;

	const velocity vi = tl::k2v(ki);
	const velocity vf = tl::k2v(kf);

	const length& lp = params.len_pulse_mono,
		&lm = params.len_mono_sample,
		&ls = params.len_sample_det;
	const length& slp = params.sig_len_pulse_mono,
		&slm = params.sig_len_mono_sample,
		&sls = params.sig_len_sample_det;
	const angle &tt_i = params.twotheta_i,
		&ph_i = params.angle_outplane_i,
		&ph_f = params.angle_outplane_f;

	const t_time &sp = params.sig_pulse,
		&sm = params.sig_mono,
		&sd = params.sig_det;
	const angle &s2ti = params.sig_twotheta_i,
		&s2tf = params.sig_twotheta_f,
		&sphi = params.sig_outplane_i,
		&sphf = params.sig_outplane_f;
	const t_time ti = lp / vi,
		tf = ls / vf;

#ifndef NDEBUG
	tl::log_debug("ki = ", ki, ", kf = ", kf);
	tl::log_debug("vi = ", vi, ", vf = ", vf);
	tl::log_debug("ti = ", ti, ", tf = ", tf);
	tl::log_debug("Q = ", Q, ", E = ", E);
#endif

	const mass mn = tl::get_m_n<t_real>();
	const auto mn_hbar = mn / tl::get_hbar<t_real>();

	// formulas 20 & 21 in [vio14]
	const t_time st = tl::my_units_sqrt<t_time>(sp*sp + sm*sm);
	const t_time stm = tl::my_units_sqrt<t_time>(sd*sd + sm*sm);

	std::vector<t_real> vecsigs = { st/sec, stm/sec, slp/meter, slm/meter, sls/meter,
		s2ti/rads, sphi/rads, s2tf/rads, sphf/rads };



	// --------------------------------------------------------------------
	// E formulas 14-18 in [vio14]

	std::vector<std::function<t_real()>> vecEderivs =
	{
		/*1*/ [&]()->t_real { return (-mn*lp*lp/(ti*ti*ti) -mn*ls*ls/(tf*tf*tf)*lm/lp) /meV*sec; },
		/*2*/ [&]()->t_real { return mn*ls*ls/(tf*tf*tf) /meV*sec; },
		/*3*/ [&]()->t_real { return (mn*lp/(ti*ti) + mn*(ls*ls)/(tf*tf*tf)*ti/(lp*lp)*lm) /meV*meter; },
		/*4*/ [&]()->t_real { return -mn*(ls*ls)/(tf*tf*tf) * ti/lp /meV*meter; },
		/*5*/ [&]()->t_real { return -mn*ls/(tf*tf) /meV*meter; },
	};


#ifndef NDEBUG
	// formula 19 in [vio14]
	t_real sigE = std::sqrt(
		std::inner_product(vecEderivs.begin(), vecEderivs.end(), vecsigs.begin(), t_real(0),
		[](t_real r1, t_real r2)->t_real { return r1 + r2; },
		[](const std::function<t_real()>& f1, t_real r2)->t_real { return f1()*f1()*r2*r2; } ));
	tl::log_debug("dE (Vanadium fwhm) = ", tl::get_SIGMA2FWHM<t_real>()*sigE);


	// --------------------------------------------------------------------
	// checking results against the E resolution formula from [ehl11], p. 6
	auto vi3 = vi*vi*vi/lp; auto vf3 = vf*vf*vf/ls;
	ublas::vector<energy> vecE2 = tl::make_vec<ublas::vector<energy>>(
	{
		mn * sp * (vi3 + vf3*lm/lp),
		mn * sm * (vi3 + vf3*(lp+lm)/lp),
		mn * sd * vf3
	});
	t_real sigE2 = tl::my_units_norm2<energy>(vecE2) / meV;
	tl::log_debug("dE (Vanadium fwhm, check) = ", tl::get_SIGMA2FWHM<t_real>()*sigE2);
	// --------------------------------------------------------------------

#endif
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// spherical: Q formulas in appendices A.1 and p. 34 of [vio14]
	// cylindrical: Q formulas in appendices A.3 and p. 35 of [vio14]

	t_real ctt_i = std::cos(tt_i/rads), stt_i = std::sin(tt_i/rads);
	t_real ctt_f = std::cos(tt/rads), stt_f = std::sin(tt/rads);
	t_real cph_i = std::cos(ph_i/rads), sph_i = std::sin(ph_i/rads);
	t_real cph_f = std::cos(ph_f/rads), sph_f = std::sin(ph_f/rads);
	t_real tph_f = std::tan(ph_f/rads), dtph_f = t_real(1) + std::pow(std::tan(ph_f/rads), t_real(2));

	t_mat R, R_tt_f, R_ph_f;
	t_mat R_tt_i = tl::make_mat(	// R derivs w.r.t the angles
		{{ -stt_i*cph_i, t_real(0) },
		{ ctt_i*cph_i, 	t_real(0) },
		{ t_real(0), 	t_real(0) }});
	t_mat R_ph_i = tl::make_mat(
		{{ -ctt_i*sph_i, t_real(0) },
		{ -stt_i*sph_i,  t_real(0) },
		{ cph_i, 	 t_real(0) }});

	if(params.det_shape == TofDetShape::SPH)
	{
		R = tl::make_mat<t_mat>(	// R: spherical coordinates
			{{ ctt_i*cph_i,	-ctt_f*cph_f },
			{ stt_i*cph_i, 	-stt_f*cph_f },
			{ sph_i, 	-sph_f }});
		R_tt_f = tl::make_mat(		// R derivs w.r.t the angles
			{{ t_real(0),	stt_f*cph_f },
			{  t_real(0),	-ctt_f*cph_f },
			{  t_real(0),	t_real(0) }});
		R_ph_f = tl::make_mat(
			{{ t_real(0),	ctt_f*sph_f },
			{  t_real(0), 	stt_f*sph_f },
			{  t_real(0), 	-cph_f }});
	}
	else if(params.det_shape == TofDetShape::CYL)
	{	// modified original formulas; TODO: check if still correct
		R = tl::make_mat<t_mat>(	// R: cylindrical coordinates
			{{ ctt_i*cph_i,	-ctt_f },
			{ stt_i*cph_i, 	-stt_f },
			{ sph_i, 	-tph_f }});
		R_tt_f = tl::make_mat(		// R derivs w.r.t the angles
			{{ t_real(0),	stt_f },
			{  t_real(0),	-ctt_f },
			{  t_real(0),	t_real(0) }});
		R_ph_f = tl::make_mat(
			{{ t_real(0),	t_real(0) },
			{  t_real(0), 	t_real(0) },
			{  t_real(0), 	-dtph_f }});
	}
	else
	{
		res.bOk = false;
		res.strErr = "Unknown detector shape.";
		return res;
	}
	//	tl::log_debug("R = ", R);

	t_vec vecViVf = tl::make_vec<t_vec>({ mn_hbar*vi *angs, mn_hbar*vf *angs });
	std::vector<std::function<t_vec()>> vecQderivs =
	{
		/*1*/ [&]()->t_vec { return ublas::prod(R, tl::make_vec<t_vec>(
			{ -mn_hbar*vi/ti *angs*sec, mn_hbar*vf/tf * lm/lp *angs*sec })); },
		/*2*/ [&]()->t_vec { return ublas::prod(R, tl::make_vec<t_vec>(
			{ t_real(0), -mn_hbar*vf/tf *angs*sec })); },
		/*3*/ [&]()->t_vec { return ublas::prod(R, tl::make_vec<t_vec>(
			{ mn_hbar/ti *angs*meter, -mn_hbar*vf/tf * lm/(vi*lp) *angs*meter })); },
		/*4*/ [&]()->t_vec { return ublas::prod(R, tl::make_vec<t_vec>(
			{ t_real(0), mn_hbar*vf/tf / vi *angs*meter })); },
		/*5*/ [&]()->t_vec { return ublas::prod(R, tl::make_vec<t_vec>(
			{ t_real(0), mn_hbar/tf *angs*meter })); },

		/*6*/ [&]()->t_vec { return ublas::prod(R_tt_i, vecViVf); },
		/*7*/ [&]()->t_vec { return ublas::prod(R_ph_i, vecViVf); },
		/*8*/ [&]()->t_vec { return ublas::prod(R_tt_f, vecViVf); },
		/*9*/ [&]()->t_vec { return ublas::prod(R_ph_f, vecViVf); },
	};


#ifndef NDEBUG
	t_vec vecQsq = std::inner_product(vecQderivs.begin(), vecQderivs.end(),
		vecsigs.begin(), tl::make_vec<t_vec>({0,0,0}),
		[](const t_vec& vec1, const t_vec& vec2) -> t_vec { return vec1 + vec2; },
		[](const std::function<t_vec()>& f1, const t_real r2) -> t_vec
		{
			const t_vec vec1 = f1();
			return ublas::element_prod(vec1, vec1) * r2*r2;
		});

	t_vec sigQ = tl::apply_fkt(vecQsq, static_cast<t_real(*)(t_real)>(std::sqrt));
	tl::log_debug("dQ (Vanadium fwhm) = ", tl::get_SIGMA2FWHM<t_real>()*sigQ);


	// --------------------------------------------------------------------
	// checking results against the Q resolution formulas from [ehl11], pp. 6-7
	auto vi2 = vi*vi/lp; auto vf2 = vf*vf/ls;
	ublas::vector<wavenumber> vecQ2[] = {
	tl::make_vec<ublas::vector<wavenumber>>(
	{
		mn_hbar * sp * (vi2 + vf2*lm/lp * ctt_f),
		mn_hbar * sm * (vi2 + vf2*(lp+lm)/lp * ctt_f),
		mn_hbar * sd * (vf2 * ctt_f),
		mn_hbar * s2tf/rads * (vf * stt_f)
	}),
	tl::make_vec<ublas::vector<wavenumber>>(
	{
		mn_hbar * sp * (vf2*lm/lp * stt_f),
		mn_hbar * sm * (vf2*(lp+lm)/lp * stt_f),
		mn_hbar * sd * (vf2 * stt_f),
		mn_hbar * s2tf/rads * (vf * ctt_f)
	}) };

	t_vec sigQ2 = tl::make_vec<t_vec>({
		tl::my_units_norm2<wavenumber>(vecQ2[0]) * angs,
		tl::my_units_norm2<wavenumber>(vecQ2[1]) * angs });
	tl::log_debug("dQ (Vanadium fwhm, check) = ", tl::get_SIGMA2FWHM<t_real>()*sigQ2);
	// --------------------------------------------------------------------

#endif
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// formulas 10 & 11 in [vio14]
	t_mat matSigSq = tl::diag_matrix({
		st*st /sec/sec, stm*stm /sec/sec,
		slp*slp /meter/meter, slm*slm /meter/meter, sls*sls /meter/meter,
		s2ti*s2ti /rads/rads, sphi*sphi /rads/rads,
		s2tf*s2tf /rads/rads, sphf*sphf /rads/rads });

	std::size_t N = matSigSq.size1();
	t_mat matJacobiInstr(4, N, t_real(0));
	for(std::size_t iDeriv = 0; iDeriv < vecQderivs.size(); ++iDeriv)
		tl::set_column(matJacobiInstr, iDeriv, vecQderivs[iDeriv]());
	for(std::size_t iDeriv = 0; iDeriv < vecEderivs.size(); ++iDeriv)
		matJacobiInstr(3, iDeriv) = vecEderivs[iDeriv]();

	t_mat matSigQE = tl::transform_inv(matSigSq, matJacobiInstr, true);
	if(!tl::inverse(matSigQE, res.reso))
	{
		res.bOk = false;
		res.strErr = "Jacobi matrix cannot be inverted.";
		return res;
	}

#ifndef NDEBUG
	tl::log_debug("J_instr = ", matJacobiInstr);
	tl::log_debug("J_QE = ", matSigQE);
	tl::log_debug("Reso = ", res.reso);
#endif
	// --------------------------------------------------------------------

	// transform from  (ki, ki_perp, Q_z)  to  (Q_para, Q_perp, Q_z)  system
	t_mat matKiQ = tl::rotation_matrix_2d(-params.angle_ki_Q / rads);
	matKiQ.resize(4, 4, true);
	matKiQ(2, 2) = matKiQ(3, 3) = 1.;
	matKiQ(2, 0) = matKiQ(2, 1) = matKiQ(2, 3) = matKiQ(3, 0) = matKiQ(3, 1) =
	matKiQ(3, 2) = matKiQ(0, 2) = matKiQ(0, 3) = matKiQ(1, 2) = matKiQ(1, 3) = 0.;

	res.reso = tl::transform(res.reso, matKiQ, true);
	//res.reso *= tl::get_SIGMA2FWHM<t_real>()*tl::get_SIGMA2FWHM<t_real>();
	res.dResVol = tl::get_ellipsoid_volume(res.reso);
	res.dR0 = 1.;   // TODO

	// Bragg widths
	const std::vector<t_real> vecFwhms = calc_bragg_fwhms(res.reso);
	std::copy(vecFwhms.begin(), vecFwhms.end(), res.dBraggFWHMs);

	res.reso_v = ublas::zero_vector<t_real>(4);
	res.reso_s = 0.;

	res.bOk = true;
	return res;
}
