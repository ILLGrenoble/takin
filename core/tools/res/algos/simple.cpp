/**
 * simple resolution calculation including only ki and kf errors
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2016
 * @license GPLv2
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

#include "simple.h"
#include "../ellipse.h"

#include "tlibs/math/linalg.h"
#include "tlibs/math/geo.h"
#include "tlibs/math/math.h"
#include "tlibs/log/log.h"

#include <string>
#include <iostream>


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


ResoResults calc_simplereso(const SimpleResoParams& params)
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

#ifndef NDEBUG
	tl::log_debug("ki = ", ki, ", kf = ", kf);
	tl::log_debug("Q = ", Q, ", E = ", E);
#endif

	const mass mn = tl::get_m_n<t_real>();
	const auto hbar = tl::get_hbar<t_real>();


	const t_real ctt = units::cos(tt);
	const t_real stt = units::sin(tt);

	// --------------------------------------------------------------------
	// sigmas
	const wavenumber& sig_kix = params.sig_ki;
	const wavenumber sig_kiy = params.sig_ki_perp;
	const wavenumber sig_kiz = params.sig_ki_z;

	const wavenumber sig_kfx = params.sig_kf * ctt - params.sig_kf_perp * stt;
	const wavenumber sig_kfy = params.sig_kf * stt + params.sig_kf_perp * ctt;
	const wavenumber sig_kfz = params.sig_kf_z;

	t_mat matSigSq = tl::diag_matrix({
		sig_kix*sig_kix * angs*angs,
		sig_kiy*sig_kiy * angs*angs,
		sig_kiz*sig_kiz * angs*angs,
		sig_kfx*sig_kfx * angs*angs,
		sig_kfy*sig_kfy * angs*angs,
		sig_kfz*sig_kfz * angs*angs });


	// --------------------------------------------------------------------
	// values
	const wavenumber& kix = ki;
	const wavenumber kiy = t_real(0)/angs;
	const wavenumber kiz = t_real(0)/angs;

	const wavenumber kfx = kf * ctt;
	const wavenumber kfy = kf * stt;
	const wavenumber kfz = t_real(0)/angs;

	// Matrix of derivatives of:
	// Q = ki - kf and E = hbar^2/(2mn)*(ki^2 - kf^2)
	t_mat matJacobiInstr(4, matSigSq.size1(), t_real(0));
	matJacobiInstr(0,0) = matJacobiInstr(1,1) = matJacobiInstr(2,2) = t_real(1);
	matJacobiInstr(0,3) = matJacobiInstr(1,4) = matJacobiInstr(2,5) = t_real(-1);
	matJacobiInstr(3,0) = hbar*hbar/mn * kix /angs/meV;
	matJacobiInstr(3,1) = hbar*hbar/mn * kiy /angs/meV;
	matJacobiInstr(3,2) = hbar*hbar/mn * kiz /angs/meV;
	matJacobiInstr(3,3) = -hbar*hbar/mn * kfx /angs/meV;
	matJacobiInstr(3,4) = -hbar*hbar/mn * kfy /angs/meV;
	matJacobiInstr(3,5) = -hbar*hbar/mn * kfz /angs/meV;

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

	// transform from  (ki, ki_perp, Qz)  to  (Q_perp, Q_para, Q_z)  system
	t_mat matKiQ = tl::rotation_matrix_2d(-params.angle_ki_Q / rads);
	matKiQ.resize(4,4, true);
	matKiQ(2,2) = matKiQ(3,3) = 1.;
	matKiQ(2,0) = matKiQ(2,1) = matKiQ(2,3) = matKiQ(3,0) = matKiQ(3,1) =
	matKiQ(3,2) = matKiQ(0,2) = matKiQ(0,3) = matKiQ(1,2) = matKiQ(1,3) = 0.;

	res.reso = tl::transform(res.reso, matKiQ, true);
	// uncomment to have parameter deviations in FWHM
	//res.reso *= tl::get_SIGMA2FWHM<t_real>()*tl::get_SIGMA2FWHM<t_real>();

	res.dResVol = tl::get_ellipsoid_volume(res.reso);
	res.dR0 = 0.;   // TODO

	// Bragg widths
	const std::vector<t_real> vecFwhms = calc_bragg_fwhms(res.reso);
	std::copy(vecFwhms.begin(), vecFwhms.end(), res.dBraggFWHMs);

	res.reso_v = ublas::zero_vector<t_real>(4);
	res.reso_s = 0.;

	res.bOk = true;
	return res;
}
