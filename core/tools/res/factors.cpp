/**
 * common factors for the resolution algorithms
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
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "defs.h"

#include "tlibs/math/math.h"
#include "tlibs/phys/neutrons.h"

namespace units = boost::units;


typedef t_real_reso t_real;
typedef ublas::matrix<t_real> t_mat;

using angle = tl::t_angle_si<t_real>;
using wavenumber = tl::t_wavenumber_si<t_real>;

static const auto angs = tl::get_one_angstrom<t_real>();
static const auto rads = tl::get_one_radian<t_real>();


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

	if(flags & CALC_KI3)         // monochromator reflectivity factor
		dmono *= tl::ana_effic_factor(ki, units::abs(thetam));
	if(flags & CALC_KF3)         // analyser reflectivity factor
		dana *= tl::ana_effic_factor(kf, units::abs(thetaa));
	if(flags & CALC_KF)
		dSqwToXSec *= kf * angs;   // kf part of the kf/ki factor, see Shirane, equ. (2.7)
	if(flags & CALC_KI)
		dSqwToXSec *= t_real(1)/ki / angs; // 1/ki part of the kf/ki factor, see Shirane, equ. (2.7)
//	if(flags & CALC_KFKI)
//		dSqwToXSec *= kf/ki;     // kf/ki factor, see Shirane, equ. (2.7)
	if(flags & CALC_MONKI)
		dmonitor *= ki*angs;       // monitor 1/ki factor, see [zhe07], p. 10

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
