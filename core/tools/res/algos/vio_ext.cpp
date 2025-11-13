/**
 * re-implementation of V. Mecoli's extension of Violini's TOF reso algorithm [mec25, mec25b]
 * @author Tobias Weber <tobias.weber@tum.de>
 * @modif_by Victor Mecoli <mecoli@ill.fr> - nov-2025
 * @date oct-2025
 * @license GPLv2
 *
 * @desc for extended method, see: [mec25] V. Mecoli, PhD thesis in preparation, https://github.com/ILLGrenoble/takin-pytools/blob/main/resolution/algos/vio_cov_ext2.py
 * @desc for extended method, see: [mec25b] V. Mecoli, PhD thesis in preparation, https://github.com/ILLGrenoble/takin-pytools/blob/main/resolution/algos/vio_ext2.py
 * @desc for original method, see: [vio14] N. Violini et al., NIM A 736 (2014) pp. 31-39, doi: 10.1016/j.nima.2013.10.042
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

#include "vio_ext.h"
#include "../ellipse.h"

#include "tlibs/math/linalg.h"
#include "tlibs/math/geo.h"
#include "tlibs/math/math.h"
#include "tlibs/math/rand.h"
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
static const auto degs = tl::get_one_deg<t_real>();
static const length angs = tl::get_one_angstrom<t_real>();
static const energy meV = tl::get_one_meV<t_real>();
static const t_time sec = tl::get_one_second<t_real>();
static const length meter = tl::get_one_meter<t_real>();
static const mass kg = tl::get_one_kg<t_real>();
static const t_real pi = tl::get_pi<t_real>();


/**
 * re-implementation of the function "length" from [mec25], l. 52ff
 */
std::tuple<t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real>
static mc_length(t_real radS, t_real heiS,
	t_real L_PE, t_real L_ME, t_real L_ES,
	t_real wEy, t_real wEz, t_real moyPx,
	t_real sigPx, t_real moyMx, t_real sigMx,
	unsigned int mc_points = 100000)
{
	t_real LPS = 0., LPSx = 0., LPSy = 0., LPSz = 0.;
	t_real LMS = 0., LMSx = 0., LMSy = 0., LMSz = 0.;

	// TODO: thread pool
	for(unsigned int pt_idx = 0; pt_idx < mc_points; ++pt_idx)
	{
		t_real Sr = radS*std::sqrt(tl::rand01<t_real>());
		t_real theta = 2.*pi*tl::rand01<t_real>();
		t_real Sz = tl::rand_minmax<t_real>(-heiS/2., heiS/2.);
		t_real Sx = Sr*std::cos(theta);
		t_real Sy = Sr*std::sin(theta);

		const t_real Pfact = (L_PE + L_ES + Sx) / (L_ES + Sx);
		t_real Px = tl::rand_norm<t_real>(moyPx, sigPx);
		t_real Py = tl::rand_minmax<t_real>(Sy + Pfact * (-wEy - Sy), Sy + Pfact * (+wEy - Sy));
		t_real Pz = tl::rand_minmax<t_real>(Sz + Pfact * (-wEz - Sz), Sz + Pfact * (+wEz - Sz));

		const t_real Mfact = (L_ME + L_ES + Sx) / (L_ES + Sx);
		t_real Mx = tl::rand_norm<t_real>(moyMx, sigMx);
		t_real My = tl::rand_minmax<t_real>(Sy + Mfact * (-wEy - Sy), Sy + Mfact * (+wEy - Sy));
		t_real Mz = tl::rand_minmax<t_real>(Sz + Mfact * (-wEz - Sz), Sz + Mfact * (+wEz - Sz));

		LPS += std::sqrt(std::pow(Sx - Px, 2.) + std::pow(Sy - Py, 2.) + std::pow(Sz - Pz, 2.));
		LPSx += Sx - Px;
		LPSy += Sy - Py;
		LPSz += Sz - Pz;

		LMS += std::sqrt(std::pow(Sx - Mx, 2.) + std::pow(Sy - My, 2.) + std::pow(Sz - Mz, 2.));
		LMSx += Sx - Mx;
		LMSy += Sy - My;
		LMSz += Sz - Mz;
	}

	// normalise
	for(t_real *val : {&LPS, &LPSx, &LPSy, &LPSz, &LMS, &LMSx, &LMSy, &LMSz})
		*val /= t_real(mc_points);

	return std::make_tuple(LPS - LMS, LPSx - LMSx, LPSy - LMSy, LPSz - LMSz, LMS, LMSx, LMSy, LMSz);
}


/**
 * re-implementation of the function "jacobTerms" from [mec25], l. 96ff
 */
static t_mat jacobian(t_real v_i, t_real v_f,
	t_real L_PM, t_real L_PMx, t_real L_PMy, t_real L_PMz,
	t_real L_MS, t_real L_MSx, t_real L_MSy, t_real L_MSz,
	t_real L_SD, t_real L_SDx, t_real L_SDy, t_real L_SDz)
{
	const t_real meV2J = meV * sec*sec / kg / meter / meter;
	const t_real m_n = tl::get_m_n<t_real>() / kg;
	const t_real hbar = tl::get_hbar<t_real>() * sec / kg / meter / meter;;

	const t_real mn_div_lpm = m_n / L_PM / meV2J;  // includes meV->J conversion factor
	const t_real mn_div_lsd = m_n / L_SD / meV2J;  // includes meV->J conversion factor
	const t_real mn_div_hbarlpm = m_n / hbar / L_PM;
	const t_real mn_div_hbarlsd = m_n / hbar / L_SD;

	const t_real vi_sq = v_i * v_i;
	const t_real vf_sq = v_f * v_f;
	const t_real vi_cb = vi_sq * v_i;
	const t_real vf_cb = vf_sq * v_f;
	const t_real vf3_div_vi = vf_cb / v_i;
	const t_real vf2_div_vi = vf_sq / v_i;

	const t_real lms_div_lsd = L_MS / L_SD;
	const t_real lpm_div_lsd = L_PM / L_SD;

	const t_real lpmx_div_lrm = L_PMx / L_PM;
	const t_real lpmy_div_lrm = L_PMy / L_PM;
	const t_real lpmz_div_lrm = L_PMz / L_PM;

	const t_real lsdx_div_lsd = L_SDx / L_SD;
	const t_real lsdy_div_lsd = L_SDy / L_SD;
	const t_real lsdz_div_lsd = L_SDz / L_SD;

	const t_real lmsx_div_lms = L_MSx / L_MS;
	const t_real lmsy_div_lms = L_MSy / L_MS;
	const t_real lmsz_div_lms = L_MSz / L_MS;

	// Qx
	t_real dQxdPx = -mn_div_hbarlpm * (v_i + vf2_div_vi*lms_div_lsd*lsdx_div_lsd*lpmx_div_lrm);
	t_real dQxdPy = -mn_div_hbarlpm * (vf2_div_vi*lms_div_lsd*lsdx_div_lsd*lpmy_div_lrm);
	t_real dQxdPz = -mn_div_hbarlpm * (vf2_div_vi*lms_div_lsd*lsdx_div_lsd*lpmz_div_lrm);
	t_real dQxdMx = mn_div_hbarlpm * (v_i + vf2_div_vi*(lpm_div_lsd*lmsx_div_lms + lms_div_lsd*lpmx_div_lrm)*lsdx_div_lsd);
	t_real dQxdMy = mn_div_hbarlpm * (vf2_div_vi*(lpm_div_lsd*lmsy_div_lms + lms_div_lsd*lpmy_div_lrm)*lsdx_div_lsd);
	t_real dQxdMz = mn_div_hbarlpm * (vf2_div_vi*(lpm_div_lsd*lmsz_div_lms + lms_div_lsd*lpmz_div_lrm)*lsdx_div_lsd);
	t_real dQxdSx = mn_div_hbarlsd * (v_f - vf2_div_vi*lsdx_div_lsd*lmsx_div_lms);
	t_real dQxdSy = -mn_div_hbarlsd * (vf2_div_vi*lsdx_div_lsd*lmsy_div_lms);
	t_real dQxdSz = -mn_div_hbarlsd * (vf2_div_vi*lsdx_div_lsd*lmsz_div_lms);
	t_real dQxdDx = -mn_div_hbarlsd * v_f;
	t_real dQxdDy = 0.;
	t_real dQxdDz = 0.;
	t_real dQxdtp = mn_div_hbarlpm * (vi_sq*lpmx_div_lrm + vf_sq*lms_div_lsd*lsdx_div_lsd);
	t_real dQxdtm = -mn_div_hbarlpm * (vi_sq*lpmx_div_lrm + vf_sq*(lpm_div_lsd + lms_div_lsd)*lsdx_div_lsd);
	t_real dQxdtd = mn_div_hbarlsd * vf_sq*lsdx_div_lsd;

	// Qy
	t_real dQydPx = -mn_div_hbarlpm * (vf2_div_vi*lms_div_lsd*lsdy_div_lsd*lpmx_div_lrm);
	t_real dQydPy = -mn_div_hbarlpm * (v_i + vf2_div_vi*lms_div_lsd*lsdy_div_lsd*lpmy_div_lrm);
	t_real dQydPz = -mn_div_hbarlpm * (vf2_div_vi*lms_div_lsd*lsdy_div_lsd*lpmz_div_lrm);
	t_real dQydMx = mn_div_hbarlpm * (vf2_div_vi*(lpm_div_lsd*lmsx_div_lms + lms_div_lsd*lpmx_div_lrm)*lsdy_div_lsd);
	t_real dQydMy = mn_div_hbarlpm * (v_i + vf2_div_vi*(lpm_div_lsd*lmsy_div_lms + lms_div_lsd*lpmy_div_lrm)*lsdy_div_lsd);
	t_real dQydMz = mn_div_hbarlpm * (vf2_div_vi*(lpm_div_lsd*lmsz_div_lms + lms_div_lsd*lpmz_div_lrm)*lsdy_div_lsd);
	t_real dQydSx = -mn_div_hbarlsd * (vf2_div_vi*lsdy_div_lsd*lmsx_div_lms);
	t_real dQydSy = mn_div_hbarlsd * (v_f - vf2_div_vi*lsdy_div_lsd*lmsy_div_lms);
	t_real dQydSz = -mn_div_hbarlsd * (vf2_div_vi*lsdy_div_lsd*lmsz_div_lms);
	t_real dQydDx = 0.;
	t_real dQydDy = -mn_div_hbarlsd * v_f;
	t_real dQydDz = 0.;
	t_real dQydtp = mn_div_hbarlpm * (vi_sq*lpmy_div_lrm + vf_sq*lms_div_lsd*lsdy_div_lsd);
	t_real dQydtm = -mn_div_hbarlpm * (vi_sq*lpmy_div_lrm + vf_sq*(lpm_div_lsd + lms_div_lsd)*lsdy_div_lsd);
	t_real dQydtd = mn_div_hbarlsd * vf_sq*lsdy_div_lsd;

	// Qz
	t_real dQzdPx = -mn_div_hbarlpm * (vf2_div_vi*lms_div_lsd*lsdz_div_lsd*lpmx_div_lrm);
	t_real dQzdPy = -mn_div_hbarlpm * (vf2_div_vi*lms_div_lsd*lsdz_div_lsd*lpmy_div_lrm);
	t_real dQzdPz = -mn_div_hbarlpm * (v_i + vf2_div_vi*lms_div_lsd*lsdz_div_lsd*lpmz_div_lrm);
	t_real dQzdMx = mn_div_hbarlpm * (vf2_div_vi*(lpm_div_lsd*lmsx_div_lms + lms_div_lsd*lpmx_div_lrm)*lsdz_div_lsd);
	t_real dQzdMy = mn_div_hbarlpm * (vf2_div_vi*(lpm_div_lsd*lmsy_div_lms + lms_div_lsd*lpmy_div_lrm)*lsdz_div_lsd);
	t_real dQzdMz = mn_div_hbarlpm * (v_i + vf2_div_vi*(lpm_div_lsd*lmsz_div_lms + lms_div_lsd*lpmz_div_lrm)*lsdz_div_lsd);
	t_real dQzdSx = -mn_div_hbarlsd * (vf2_div_vi*lsdz_div_lsd*lmsx_div_lms);
	t_real dQzdSy = -mn_div_hbarlsd * (vf2_div_vi*lsdz_div_lsd*lmsy_div_lms);
	t_real dQzdSz = mn_div_hbarlsd * (v_f - vf2_div_vi*lsdz_div_lsd*lmsz_div_lms);
	t_real dQzdDx = 0.;
	t_real dQzdDy = 0.;
	t_real dQzdDz = -mn_div_hbarlsd * v_f;
	t_real dQzdtp = mn_div_hbarlpm * (vi_sq*lpmz_div_lrm + vf_sq*lms_div_lsd*lsdz_div_lsd);
	t_real dQzdtm = -mn_div_hbarlpm * (vi_sq*lpmz_div_lrm + vf_sq*(lpm_div_lsd + lms_div_lsd)*lsdz_div_lsd);
	t_real dQzdtd = mn_div_hbarlsd * vf_sq*lsdz_div_lsd;

	// E
	const t_real dEdPxyz = -mn_div_lpm * (vi_sq + vf3_div_vi*lms_div_lsd);
	t_real dEdPx = dEdPxyz * lpmx_div_lrm;
	t_real dEdPy = dEdPxyz * lpmy_div_lrm;
	t_real dEdPz = dEdPxyz * lpmz_div_lrm;
	t_real dEdMx = mn_div_lpm * ((vi_sq + vf3_div_vi*lms_div_lsd) * lpmx_div_lrm + vf3_div_vi*lpm_div_lsd*lmsx_div_lms);
	t_real dEdMy = mn_div_lpm * ((vi_sq + vf3_div_vi*lms_div_lsd) * lpmy_div_lrm + vf3_div_vi*lpm_div_lsd*lmsy_div_lms);
	t_real dEdMz = mn_div_lpm * ((vi_sq + vf3_div_vi*lms_div_lsd) * lpmz_div_lrm + vf3_div_vi*lpm_div_lsd*lmsz_div_lms);
	t_real dEdSx = mn_div_lsd * (vf_sq*lsdx_div_lsd - vf3_div_vi*lmsx_div_lms);
	t_real dEdSy = mn_div_lsd * (vf_sq*lsdy_div_lsd - vf3_div_vi*lmsy_div_lms);
	t_real dEdSz = mn_div_lsd * (vf_sq*lsdz_div_lsd - vf3_div_vi*lmsz_div_lms);
	t_real dEdDx = -mn_div_lsd * vf_sq*lsdx_div_lsd;
	t_real dEdDy = -mn_div_lsd * vf_sq*lsdy_div_lsd;
	t_real dEdDz = -mn_div_lsd * vf_sq*lsdz_div_lsd;
	t_real dEdtp = mn_div_lpm * (vi_cb + vf_cb*lms_div_lsd);
	t_real dEdtm = -mn_div_lpm * (vi_cb + vf_cb*(lpm_div_lsd + lms_div_lsd));
	t_real dEdtd = mn_div_lsd * vf_cb;

	return tl::make_mat<t_mat>(
	{
		{ dQxdPx, dQxdPy, dQxdPz, dQxdMx, dQxdMy, dQxdMz, dQxdSx, dQxdSy, dQxdSz, dQxdDx, dQxdDy, dQxdDz, dQxdtp, dQxdtm, dQxdtd },
		{ dQydPx, dQydPy, dQydPz, dQydMx, dQydMy, dQydMz, dQydSx, dQydSy, dQydSz, dQydDx, dQydDy, dQydDz, dQydtp, dQydtm, dQydtd },
		{ dQzdPx, dQzdPy, dQzdPz, dQzdMx, dQzdMy, dQzdMz, dQzdSx, dQzdSy, dQzdSz, dQzdDx, dQzdDy, dQzdDz, dQzdtp, dQzdtm, dQzdtd },
		{ dEdPx,   dEdPy,  dEdPz,  dEdMx,  dEdMy,  dEdMz,  dEdSx,  dEdSy,  dEdSz,  dEdDx,  dEdDy,  dEdDz,  dEdtp,  dEdtm,  dEdtd }
	}) / 1e10 / 1e10 /* angstrom conversion factor */;
}


ResoResults calc_vio_ext(const VioExtParams& params)
{
	// beginning as in vio.cpp
	ResoResults res;
	res.Q_avg.resize(4);

	const t_real &M = params.M_coating, &n_b = params.n_b;

	const angle& tt = params.twotheta;
	const wavenumber &ki = params.ki, &kf = params.kf;
	const energy E = tl::get_energy_transfer(ki, kf);
	const wavenumber Q = tl::get_sample_Q(ki, kf, tt);

	res.Q_avg[0] = Q * angs;
	res.Q_avg[1] = 0.;
	res.Q_avg[2] = 0.;
	res.Q_avg[3] = E / meV;

	const t_real vi = tl::k2v(ki) * sec / angs;
	const t_real vf = tl::k2v(kf) * sec / angs;
	const t_real c_tt = std::cos(tt / rads);
	const t_real s_tt = std::sin(tt / rads);


	// --------------------------------------------------------------------------------
	// [mec25b], l. 48ff
	// --------------------------------------------------------------------------------
	// sample
	t_real sample_rad = params.sample_width / 2. / angs;
	t_real sample_height = params.sample_height / angs;

	// pulse-shaping chopper
	t_real ch_pulse_rpm_mult = params.ch_pulse_counterrot ? 2. : 1.;
	t_real chopperP_wnd_angle = params.ch_pulse_angle_win / degs;
	t_real chopperP_beam_angle = params.ch_pulse_angle_beam / degs;
	t_real chopperP_width = params.ch_pulse_width / angs;
	t_real chopperP_rpm = ch_pulse_rpm_mult * params.ch_pulse_rpm;

	// monochromatising chopper
	t_real ch_mono_rpm_mult = params.ch_mono_counterrot ? 2. : 1.;
	t_real chopperM_wnd_angle = params.ch_mono_angle_win / degs;
	t_real chopperM_beam_angle = params.ch_mono_angle_beam / degs;
	t_real chopperM_width = params.ch_mono_width / angs;
	t_real chopperM_rpm = ch_mono_rpm_mult * params.ch_mono_rpm;

	// guide dimensions
	t_real endguide_ywidth = params.endguide_width / angs / 2.;
	t_real endguide_zheight = params.endguide_height / angs / 2.;

	// detector geometry
	t_real det_height = params.det_height / angs;
	t_real det_tube_w = params.det_tube_width / angs;
	t_real det_z = params.det_z / angs;

	// distances
	t_real dist_chP_endguide = params.len_ch_pulse_guide / angs;
	t_real dist_chM_endguide = params.len_ch_mono_guide / angs;
	t_real dist_endguide_sample = params.len_guide_sample / angs;
	t_real dist_sample_det = params.len_sample_det2 / angs;  // TODO: unite with len_sample_det

	const t_real thetacrit = M*std::arcsin((2.*pi/ki)/10.*np.sqrt(n_b/pi));
	const t_real Hcrit = endguide_zheight - dist_endguide_sample*std::tan(thetacrit);
	const t_real Hmax = endguide_zheight + dist_endguide_sample*std::tan(thetacrit);

	const t_real Dr_sq = dist_sample_det*dist_sample_det;
	const t_real Sr_sq = sample_rad*sample_rad;
	const t_real Sh_sq = sample_height*sample_height;
	const t_real Lpe_sq = dist_chP_endguide*dist_chP_endguide;
	const t_real Lme_sq = dist_chM_endguide*dist_chM_endguide;
	const t_real Les_sq = dist_endguide_sample*dist_endguide_sample;
	const t_real Eyh_sq = endguide_ywidth*endguide_ywidth;
	const t_real Ezh_sq = endguide_zheight*endguide_zheight;

	const t_real Hc_sq = Hcrit*Hcrit

	const t_real c0 = std::sqrt(1. - Sr_sq / Les_sq);
	const t_real c1 = 1. - c0;
	const t_real c2 = 1./c0 - 1.;

	t_real VarDr = det_tube_w*det_tube_w / 12.;
	t_real Vartd = VarDr / (vf*vf);

	t_real VarDtheta = std::pow(2.*det_tube_w*(dist_sample_det - std::sqrt(Dr_sq - Sr_sq)) / Sr_sq, 2.);
	t_real Vartp = (chopperP_wnd_angle*chopperP_wnd_angle + chopperP_beam_angle*chopperP_beam_angle) / (12.*std::pow(6.*chopperP_rpm, 2.));
	t_real Vartm = (chopperM_wnd_angle*chopperM_wnd_angle + chopperM_beam_angle*chopperM_beam_angle) / (12.*std::pow(6.*chopperM_rpm, 2.));

	t_real VarPx = std::pow(vi*std::sqrt(Vartp) - chopperP_width, 2.);
	t_real VarPy = 1./12.*( 3.*std::pow(Sr, 2.) + (4.*std::pow(dist_chP_endguide + dist_endguide_sample, 2.) + std::pow(Sr, 2.))*std::tan(thetacrit));
	t_real VarPz = 1
	if (sample_height <= Hcrit) {
		VarPz = 1./12. * ( Sh_sq + (4.*std::pow(dist_chP_endguide+dist_endguide_sample, 2.) + Sr_sq)*std::pow(std::tan(thetacrit), 2.) )
	} else if (sample_height <= Hmax) {
		VarPz = 1./(36.*sample_height) * ( 1./(Sr_sq*std::sqrt(Les_sq-Sr_sq))*(sample_height-2.*Hcrit)*(
        (3.*Hc_sq + std::pow(Hcrit + sample_height, 2.))*(2.*dist_chP_endguide*(Les_sq - Sr_sq + dist_chP_endguide*dist_endguide_sample) - std::sqrt(Les_sq - Sr_sq)*(2.*Lpe_sq - Sr_sq + 2.*dist_chP_endguide*dist_endguide_sample))
        + 3.*endguide_zheight*(sample_height + 2.*Hcrit)*(2.*dist_chP_endguide*(Les_sq - Sr_sq - 2.*dist_chP_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(4.*Lpe_sq + Sr_sq - 2.*dist_chP_endguide*dist_endguide_sample))
        + 12.*Ezh_sq*(2.*dist_chP_endguide*(2.*Sr_sq - 2.*Les_sq + dist_chP_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(Sr_sq - 2.*Lpe_sq + 4*dist_chP_endguide*dist_endguide_sample))
        - 3.*std::tan(thetacrit)*(sample_height + 2.*Hcrit)*(2.*Lpe_sq*(Les_sq - Sr_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*dist_endguide_sample*Sr_sq + dist_chP_endguide*Sr_sq - 2.*Lpe_sq*dist_endguide_sample))
        - 12.*std::tan(thetacrit)*endguide_zheight*(2.*Lpe_sq*(Sr_sq-Les_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*Lpe_sq*dist_endguide_sample + dist_endguide_sample*Sr_sq + 2.*dist_chP_endguide*Sr_sq))
        + 3.*np.square(Sr)*np.sqrt(np.square(Les) - np.square(Sr))*(np.square(2*Lpe + 2*Les)+np.square(Sr))*np.square(np.tan(thetacrit)))
    + 6.*(4*std::pow(Hcrit, 3) + Hcrit*(4*std::pow(dist_chP_endguide + dist_endguide_sample, 2.) + Sr_sq)*std::pow(std::tan(thetacrit), 2.)))
	} else {
		VarPz = 1./(18.*sample_height) * ( 1./(Sr_sq*std::sqrt(Les_sq-Sr_sq))*(Hmax-Hcrit)*(
        (3.*Hc_sq + std::pow(Hcrit + 2.*Hmax, 2.))*(2.*dist_chP_endguide*(Les_sq - Sr_sq + dist_chP_endguide*dist_endguide_sample) - std::sqrt(Les_sq - Sr_sq)*(2.*Lpe_sq - Sr_sq + 2.*dist_chP_endguide*dist_endguide_sample))
        + 3.*endguide_zheight*(2.*Hmax + 2.*Hcrit)*(2.*dist_chP_endguide*(Les_sq - Sr_sq - 2.*dist_chP_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(4.*Lpe_sq + Sr_sq - 2.*dist_chP_endguide*dist_endguide_sample))
        + 12.*Ezh_sq*(2.*dist_chP_endguide*(2.*Sr_sq - 2.*Les_sq + dist_chP_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(Sr_sq - 2.*Lpe_sq + 4*dist_chP_endguide*dist_endguide_sample))
        - 3.*std::tan(thetacrit)*(2.*Hmax + 2.*Hcrit)*(2.*Lpe_sq*(Les_sq - Sr_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*dist_endguide_sample*Sr_sq + dist_chP_endguide*Sr_sq - 2.*Lpe_sq*dist_endguide_sample))
        - 12.*std::tan(thetacrit)*endguide_zheight*(2.*Lpe_sq*(Sr_sq-Les_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*Lpe_sq*dist_endguide_sample + dist_endguide_sample*Sr_sq + 2.*dist_chP_endguide*Sr_sq))
        + 3.*np.square(Sr)*np.sqrt(np.square(Les) - np.square(Sr))*(np.square(2*Lpe + 2*Les)+np.square(Sr))*np.square(np.tan(thetacrit)))
    + 3.*(4*std::pow(Hcrit, 3) + Hcrit*(4*std::pow(dist_chP_endguide + dist_endguide_sample, 2.) + Sr_sq)*std::pow(std::tan(thetacrit), 2.)))
	}

	t_real VarMx = std::pow(vi*std::sqrt(Vartm) - chopperM_width, 2.);
	t_real VarMy = 1./12.*( 3.*std::pow(Sr, 2.) + (4.*std::pow(dist_chM_endguide + dist_endguide_sample, 2.) + std::pow(Sr, 2.))*std::tan(thetacrit));
	t_real VarMz = 1
	if (sample_height <= Hcrit) {
		VarMz = 1./12. * ( Sh_sq + (4.*std::pow(dist_chM_endguide+dist_endguide_sample, 2.) + Sr_sq)*std::pow(std::tan(thetacrit), 2.) )
	} else if (sample_height <= Hmax) {
		VarMz = 1./(36.*sample_height) * ( 1./(Sr_sq*std::sqrt(Les_sq-Sr_sq))*(sample_height-2.*Hcrit)*(
        (3.*Hc_sq + std::pow(Hcrit + sample_height, 2.))*(2.*dist_chM_endguide*(Les_sq - Sr_sq + dist_chM_endguide*dist_endguide_sample) - std::sqrt(Les_sq - Sr_sq)*(2.*Lme_sq - Sr_sq + 2.*dist_chM_endguide*dist_endguide_sample))
        + 3.*endguide_zheight*(sample_height + 2.*Hcrit)*(2.*dist_chM_endguide*(Les_sq - Sr_sq - 2.*dist_chM_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(4.*Lme_sq + Sr_sq - 2.*dist_chM_endguide*dist_endguide_sample))
        + 12.*Ezh_sq*(2.*dist_chM_endguide*(2.*Sr_sq - 2.*Les_sq + dist_chM_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(Sr_sq - 2.*Lme_sq + 4*dist_chM_endguide*dist_endguide_sample))
        - 3.*std::tan(thetacrit)*(sample_height + 2.*Hcrit)*(2.*Lme_sq*(Les_sq - Sr_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*dist_endguide_sample*Sr_sq + dist_chM_endguide*Sr_sq - 2.*Lme_sq*dist_endguide_sample))
        - 12.*std::tan(thetacrit)*endguide_zheight*(2.*Lme_sq*(Sr_sq-Les_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*Lme_sq*dist_endguide_sample + dist_endguide_sample*Sr_sq + 2.*dist_chM_endguide*Sr_sq))
        + 3.*np.square(Sr)*np.sqrt(np.square(Les) - np.square(Sr))*(np.square(2*Lpe + 2*Les)+np.square(Sr))*np.square(np.tan(thetacrit)))
    + 6.*(4*std::pow(Hcrit, 3) + Hcrit*(4*std::pow(dist_chM_endguide + dist_endguide_sample, 2.) + Sr_sq)*std::pow(std::tan(thetacrit), 2.)))
	} else {
		VarMz = 1./(18.*sample_height) * ( 1./(Sr_sq*std::sqrt(Les_sq-Sr_sq))*(Hmax-Hcrit)*(
        (3.*Hc_sq + std::pow(Hcrit + 2.*Hmax, 2.))*(2.*dist_chM_endguide*(Les_sq - Sr_sq + dist_chM_endguide*dist_endguide_sample) - std::sqrt(Les_sq - Sr_sq)*(2.*Lme_sq - Sr_sq + 2.*dist_chM_endguide*dist_endguide_sample))
        + 3.*endguide_zheight*(2.*Hmax + 2.*Hcrit)*(2.*dist_chM_endguide*(Les_sq - Sr_sq - 2.*dist_chM_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(4.*Lme_sq + Sr_sq - 2.*dist_chM_endguide*dist_endguide_sample))
        + 12.*Ezh_sq*(2.*dist_chM_endguide*(2.*Sr_sq - 2.*Les_sq + dist_chM_endguide*dist_endguide_sample) + std::sqrt(Les_sq - Sr_sq)*(Sr_sq - 2.*Lme_sq + 4*dist_chM_endguide*dist_endguide_sample))
        - 3.*std::tan(thetacrit)*(2.*Hmax + 2.*Hcrit)*(2.*Lme_sq*(Les_sq - Sr_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*dist_endguide_sample*Sr_sq + dist_chM_endguide*Sr_sq - 2.*Lme_sq*dist_endguide_sample))
        - 12.*std::tan(thetacrit)*endguide_zheight*(2.*Lme_sq*(Sr_sq-Les_sq) + std::sqrt(Les_sq - Sr_sq)*(2.*Lme_sq*dist_endguide_sample + dist_endguide_sample*Sr_sq + 2.*dist_chM_endguide*Sr_sq))
        + 3.*np.square(Sr)*np.sqrt(np.square(Les) - np.square(Sr))*(np.square(2*Lpe + 2*Les)+np.square(Sr))*np.square(np.tan(thetacrit)))
    + 3.*(4*std::pow(Hcrit, 3) + Hcrit*(4*std::pow(dist_chM_endguide + dist_endguide_sample, 2.) + Sr_sq)*std::pow(std::tan(thetacrit), 2.)))
	}
	t_real VarMy = Eyh_sq/3. + Lme_sq*(2.*Les_sq/Sr_sq*c1 - 1.) + 4.*dist_chM_endguide*dist_endguide_sample/(3.*Sr_sq)*Eyh_sq*c1 + 2.*Lme_sq/(3.*Sr_sq)*Eyh_sq*c2;
	t_real VarMz = Ezh_sq/3. + Lme_sq/(3.*Sr_sq)*(Sh_sq/2. + 2.*Ezh_sq)*c2 + 4.*dist_chM_endguide*dist_endguide_sample/(3.*Sr_sq)*Ezh_sq*c1;

	t_real VarSx = Sr_sq/4.;
	t_real VarSy = Sr_sq/4.;
	t_real VarSz = Sh_sq/12.;

	t_real VarDx = c_tt*c_tt*VarDr + Dr_sq*s_tt*s_tt*VarDtheta;
	t_real VarDy = s_tt*s_tt*VarDr + Dr_sq*c_tt*c_tt*VarDtheta;
	t_real VarDz = std::pow(det_height / 100., 2.);
	t_real CovDxDy = c_tt*s_tt*VarDr - Dr_sq*c_tt*s_tt*VarDtheta;

	t_mat cov = tl::diag_matrix(
	{
		VarPx, VarPy, VarPz,
		VarMx, VarMy, VarMz,
		VarSx, VarSy, VarSz,
		VarDx, VarDy, VarDz,
		Vartp, Vartm, Vartd
	});
	cov(9, 10) = cov(10, 9) = CovDxDy;

	t_real LPM = 0., LPMx = 0., LPMy = 0., LPMz = 0., LMS = 0, LMSx = 0., LMSy = 0., LMSz = 0.;
	std::tie(LPM, LPMx, LPMy, LPMz, LMS, LMSx, LMSy, LMSz) = mc_length(
		sample_rad, sample_height,
		dist_chP_endguide,  dist_chM_endguide, dist_endguide_sample,
		endguide_ywidth, endguide_zheight, -(dist_chP_endguide + dist_endguide_sample),
		std::sqrt(VarPx), -dist_chM_endguide - dist_endguide_sample, std::sqrt(VarMx),
		params.mc_lengths);
	t_real LSDz = det_z;
	t_real LSDx = dist_sample_det*c_tt;
	t_real LSDy = dist_sample_det*s_tt;

	t_mat J = jacobian(vi, vf,
		LPM, LPMx, LPMy, LPMz,
		LMS, LMSx, LMSy, LMSz,
		dist_sample_det, LSDx, LSDy, LSDz);
	// --------------------------------------------------------------------------------


	// rest as in vio.cpp
	t_mat covQhw = tl::transform_inv(cov, J, true);
	if(!tl::inverse(covQhw, res.reso))
	{
		res.bOk = false;
		res.strErr = "Jacobi matrix cannot be inverted.";
		return res;
	}

	// transform from  (ki, ki_perp, Q_z)  to  (Q_para, Q_perp, Q_z)  system
	t_mat matKiQ = tl::rotation_matrix_2d(-params.angle_ki_Q / rads);
	matKiQ.resize(4, 4, true);
	matKiQ(2, 2) = matKiQ(3, 3) = 1.;
	matKiQ(2, 0) = matKiQ(2, 1) = matKiQ(2, 3) = matKiQ(3, 0) = matKiQ(3, 1) =
	matKiQ(3, 2) = matKiQ(0, 2) = matKiQ(0, 3) = matKiQ(1, 2) = matKiQ(1, 3) = 0.;

	res.reso = tl::transform(res.reso, matKiQ, true);
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
