/**
 * tlibs2 -- physics library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2012 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Additional functions were forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Further functions and updates forked on 1-Feb-2021 and 19-Apr-2021 from my privately developed "geo" and "misc" projects (https://github.com/t-weber/geo and https://github.com/t-weber/misc).
 *
 * @note for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", and "misc" projects
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

#ifndef __TLIBS2_PHYS__
#define __TLIBS2_PHYS__

#include "units.h"
#include "maths.h"

#include <stdexcept>
#include <optional>
#include <boost/units/pow.hpp>


namespace tl2 {


// --------------------------------------------------------------------------------
// constants
// --------------------------------------------------------------------------------
// import scipy.constants as co
// E2KSQ = 2.*co.neutron_mass/(co.Planck/co.elementary_charge*1000./2./co.pi)**2. / co.elementary_charge*1000. * 1e-20
template<class T = double> constexpr T KSQ2E = T(0.5) * hbar<T>/angstrom<T>/m_n<T> * hbar<T>/angstrom<T>/meV<T>;
template<class T = double> constexpr T E2KSQ = T(1) / KSQ2E<T>;
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * differentiated Bragg equation:
 * n lam = 2d sin(th)				| diff
 * n dlam = 2dd sin(th) + 2d cos(th) dth	| / Bragg equ
 * dlam/lam = dd/d + cos(th)/sin(th) dth
 *
 * n G = 2k sin(th)
 * n dG = 2dk sin(th) + 2k cos(th) dth
 * dG/G = dk/k + cos(th)/sin(th) dth
 */
template<class Sys, class T = double>
T bragg_diff(T dDoverD, const t_angle<Sys, T>& theta, T dTheta)
{
	T dLamOverLam = dDoverD + units::cos(theta)/units::sin(theta) * dTheta;
	return dLamOverLam;
}


/**
 * kinematic plane
 * @see (ILL Neutron Data Booklet), sec. 2.6-2
 *
 * Q_vec = ki_vec - kf_vec
 * Q^2 = ki^2 + kf^2 - 2ki kf cos 2th	| * hbar^2 / (2 mn)
 *
 * using: Ei = hbar^2 ki^2 / (2 mn)
 * Q^2 * hbar^2 / (2 mn) = Ei + Ef - 2 ki kf cos(2th) * hbar^2 / (2 mn)
 *
 * using: ki^2 = 2 mn Ei / hbar^2
 * Q^2 = [Ei + Ef - 2 sqrt(Ei) sqrt(Ef) cos 2th] * 2 mn / hbar^2
 *
 * using: dE = Ei - Ef, Ef = Ei - dE
 * Q^2 = [2 Ei - dE - 2 sqrt(Ei (Ei - dE)) cos 2th] * 2 mn / hbar^2
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T> kinematic_plane(bool bFixedKi,
	const t_energy<Sys, T>& EiEf, const t_energy<Sys, T>& DeltaE,
	const t_angle<Sys, T>& twotheta)
{
	t_energy<Sys, T> dE = DeltaE;
	if(bFixedKi)
		dE = -dE;

	auto c = T(2.)*m_n<T> / (hbar<T>*hbar<T>);
	t_wavenumber<Sys, T> Q =
		my_units_sqrt<t_wavenumber<Sys, T>>(c *
			(T(2.)*EiEf + dE - T(2.)*units::cos(twotheta) *
			my_units_sqrt<t_wavenumber<Sys, T>>(EiEf*(EiEf + dE))));

	return Q;
}


/**
 * kinematic plane
 * @see (ILL Neutron Data Booklet), sec. 2.6-2
 *
 * solving the above equation for dE using sage:
 *   Q, Ei, dE, ctt, c = var("Q, Ei, dE, ctt, c")
 *   equ = (Q^2 -2*Ei*c + dE*c)^2 == 4*Ei*(Ei-dE)*c^2*ctt^2
 *   equ.solve(dE)
 */
template<class Sys, class T = double>
t_energy<Sys, T> kinematic_plane(bool bFixedKi, bool bBranch,
	const t_energy<Sys, T>& EiEf, const t_wavenumber<Sys, T>& Q,
	const t_angle<Sys, T>& twotheta)
{
	using t_cE = units::quantity<units::unit<typename units::derived_dimension<
		units::length_base_dimension, -2>::type,
		Sys>, T>;

	auto c = T(2.)*m_n<T> / (hbar<T>*hbar<T>);
	auto c2 = c*c;
	auto EiEf2 = EiEf*EiEf;

	T ctt = units::cos(twotheta);
	T ctt2 = ctt*ctt;

	T dSign = bBranch ? T(1.) : T(-1.);
	T dSignFixedKf = bFixedKi ? T(-1.) : T(1.);

	t_energy<Sys, T> dE =
		dSignFixedKf*T(2.) * EiEf * ctt2
		- dSignFixedKf*T(2.) * EiEf
		+ dSignFixedKf * Q*Q / c
		+ dSign*T(2.) * ctt/c * my_units_sqrt<t_cE>(
			c2*ctt2*EiEf2 - c2*EiEf2 + c*EiEf*Q*Q);

	return dE;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// tas calculations
// @see M. D. Lumsden, J. L. Robertson, and M. Yethiraj, J. Appl. Crystallogr. 38(3), pp. 405–411 (2005), doi: 10.1107/S0021889805004875.
// @see (Shirane 2002), Ch. 1.3
// --------------------------------------------------------------------------------

/**
 * angle between ki and kf in the scattering triangle (a4)
 * @returns nullopt if the angle can't be reached
 *
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * 2*<ki|kf> = ki^2 + kf^2 - Q^2
 * cos phi = (ki^2 + kf^2 - Q^2) / (2 ki kf)
 */
template<typename t_real = double>
std::optional<t_real> calc_tas_angle_ki_kf(
	t_real ki, t_real kf, t_real Q, t_real sense = 1)
{
	t_real c = (ki*ki + kf*kf - Q*Q) / (t_real(2)*ki*kf);
	if(std::abs(c) > t_real(1))
		return std::nullopt;
	return sense*std::acos(c);
}


/**
 * angle between ki and kf in the scattering triangle (a4)
 * (version with units)
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_angle_ki_kf(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf, const t_wavenumber<Sys, T>& Q,
	bool bPosSense = true)
{
	t_dimensionless<Sys, T> ttCos = (ki*ki + kf*kf - Q*Q)/(T(2.)*ki*kf);
	if(units::abs(ttCos) > T(1.))
		throw std::runtime_error("Scattering triangle not closed.");

	t_angle<Sys, T> tt = units::acos(ttCos);

	if(!bPosSense)
		tt = -tt;
	return tt;
}


/**
 * angle between ki and Q in the scattering triangle
 * @returns nullopt if the angle can't be reached
 *
 * |Q> = |ki> - |kf>
 * |kf> = |ki> + |Q>
 * kf^2 = ki^2 + Q^2 - 2*<ki|Q>
 * 2*<ki|Q> = ki^2 + Q^2 - kf^2
 * cos phi = (ki^2 + Q^2 - kf^2) / (2 ki*Q)
 */
template<typename t_real = double>
std::optional<t_real> calc_tas_angle_ki_Q(
	t_real ki, t_real kf, t_real Q, t_real sense = 1)
{
	t_real c = (ki*ki + Q*Q - kf*kf) / (t_real(2)*ki*Q);
	if(std::abs(c) > t_real(1))
		return std::nullopt;
	return sense*std::acos(c);
}


/**
 * angle between ki and Q in the scattering triangle
 * (version with units)
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_angle_ki_Q(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf,
	const t_wavenumber<Sys, T>& Q,
	bool bPosSense = true,
	bool bAngleOutsideTriag = false)
{
	t_angle<Sys, T> angle;

	if(Q*angstrom<T> == T(0.))
	{
		angle = pi<T>/T(2) * radians<T>;
	}
	else
	{
		auto c = (ki*ki - kf*kf + Q*Q) / (T(2.)*ki*Q);
		if(units::abs(c) > T(1.))
			throw std::runtime_error("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(bAngleOutsideTriag)
		angle = pi<T>*radians<T> - angle;
	if(!bPosSense)
		angle = -angle;

	return angle;
}


/**
 * angle between kf and Q in the scattering triangle
 * (version with units)
 *
 * Q_vec = ki_vec - kf_vec
 * ki_vec = Q_vec + kf_vec
 * ki^2 = Q^2 + kf^2 + 2Q kf cos th
 * cos th = (ki^2 - Q^2 - kf^2) / (2Q kf)
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_angle_kf_Q(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf,
	const t_wavenumber<Sys, T>& Q,
	bool bPosSense = true,
	bool bAngleOutsideTriag = true)
{
	t_angle<Sys, T> angle;

	if(Q*angstrom<T> == T(0.))
		angle = pi<T>/T(2) * radians<T>;
	else
	{
		auto c = (ki*ki - kf*kf - Q*Q) / (T(2.)*kf*Q);
		if(units::abs(c) > T(1.))
			throw std::runtime_error("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(!bAngleOutsideTriag)
		angle = pi<T>*radians<T> - angle;
	if(!bPosSense)
		angle = -angle;

	return angle;
}


/**
 * get length of Q
 * |Q> = |ki> - |kf>
 * Q^2 = ki^2 + kf^2 - 2*<ki|kf>
 * Q^2 = ki^2 + kf^2 - 2*ki*kf*cos(a4)
 */
template<typename t_real = double>
t_real calc_tas_Q_len(t_real ki, t_real kf, t_real a4)
{
	t_real Qsq = ki*ki + kf*kf - t_real(2)*ki*kf*std::cos(a4);
	return std::sqrt(Qsq);
}


/**
 * get length of Q (version with units)
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T>
calc_tas_Q_len(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf, const t_angle<Sys, T>& tt)
{
	t_dimensionless<Sys, T> ctt = units::cos(tt);
	decltype(ki*ki) Qsq = ki*ki + kf*kf - T(2.)*ki*kf*ctt;

	if(T(Qsq*angstrom<T>*angstrom<T>) < T(0.))
	{
		// TODO
		Qsq = -Qsq;
	}

	t_wavenumber<Sys, T> Q = my_units_sqrt<t_wavenumber<Sys, T>>(Qsq);
	return Q;
}


/**
 * get tas a3 and a4 angles
 * @return [a3, a4, distance of Q to the scattering plane]
 * @see M. D. Lumsden, et al., doi: 10.1107/S0021889805004875.
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::tuple<bool, t_real, t_real, t_real> calc_tas_a3a4(
	const t_mat& B, t_real ki_lab, t_real kf_lab,
	const t_vec& Q_rlu, const t_vec& orient_rlu, const t_vec& orient_up_rlu,
	t_real sample_sense = 1, t_real a3_offs = pi<t_real>)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	// metric from crystal B matrix
	t_mat G = metric<t_mat>(B);

	// length of Q vector
	t_real Q_len_lab = norm<t_mat, t_vec>(G, Q_rlu);

	// angle xi between Q and orientation reflex
	t_real xi = angle<t_mat, t_vec>(G, Q_rlu, orient_rlu);

	// sign/direction of xi
	t_vec xivec = cross<t_mat, t_vec>(G, orient_rlu, Q_rlu);
	t_real xidir = inner<t_mat, t_vec>(G, xivec, orient_up_rlu);
	if(xidir < t_real(0))
		xi = -xi;

	// angle psi between ki and Q
	std::optional<t_real> psi =
		calc_tas_angle_ki_Q<t_real>(ki_lab, kf_lab, Q_len_lab, sample_sense);
	if(!psi)
		return std::make_tuple(false, 0, 0, 0);

	// crystal and scattering angle
	t_real a3 = - *psi - xi + a3_offs;
	std::optional<t_real> a4 =
		calc_tas_angle_ki_kf<t_real>(ki_lab, kf_lab, Q_len_lab);
	if(!a4)
		return std::make_tuple(false, a3, 0, 0);
	*a4 *= sample_sense;

	// distance of Q to the scattering plane
	t_real dist_Q_plane = inner<t_mat, t_vec>(G, Q_rlu, orient_up_rlu);
	dist_Q_plane /= norm<t_mat, t_vec>(G, orient_up_rlu);

	return std::make_tuple(true, a3, *a4, dist_Q_plane);
}


/**
 * get hkl position of a tas
 * @return Q_rlu
 * @see M. D. Lumsden, et al., doi: 10.1107/S0021889805004875.
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::optional<t_vec> calc_tas_hkl(
	const t_mat& B, t_real ki_lab, t_real kf_lab, t_real Q_len_lab, t_real a3,
	const t_vec& orient_rlu, const t_vec& orient_up_rlu,
	t_real sample_sense = 1, t_real a3_offs = pi<t_real>)
requires is_basic_mat<t_mat> && is_basic_vec<t_vec>
{
	auto [Binv, ok] = inv<t_mat>(B);
	if(!ok)
		return std::nullopt;

	// angle psi between ki and Q
	std::optional<t_real> psi =
		calc_tas_angle_ki_Q<t_real>(ki_lab, kf_lab, Q_len_lab, sample_sense);
	if(!psi)
		return std::nullopt;

	// angle xi between Q and orientation reflex
	t_real xi = a3_offs - a3 - *psi;

	t_vec rotaxis_lab = B * orient_up_rlu;
	t_mat rotmat = rotation<t_mat, t_vec>(rotaxis_lab, xi, false);

	t_vec orient_lab = B * orient_rlu;
	t_vec Q_lab = rotmat * orient_lab;
	Q_lab /= norm<t_vec>(Q_lab);
	Q_lab *= Q_len_lab;

	t_vec Q_rlu = Binv * Q_lab;
	return Q_rlu;
}


/**
 * get a1 or a5 angle
 * @returns nullopt of the angle can't be reached
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 *
 * Bragg: n lam = 2d sin(theta)
 * n 2pi / k = 2d sin(theta)
 * n pi / k = d sin(theta)
 * theta = asin(n pi / (k d))
 */
template<class t_real = double>
std::optional<t_real> calc_tas_a1(t_real k, t_real d)
{
	t_real sintheta = pi<t_real> / (k*d);
	if(std::abs(sintheta) > t_real(1))
		return std::nullopt;
	return std::asin(sintheta);
}


/**
 * get a2 or a6 angle
 * (version with units)
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 */
template<class Sys, class T = double>
t_angle<Sys, T> calc_tas_a1(const t_wavenumber<Sys, T>& k,
	const t_length<Sys, T>& d, bool bPosSense = true)
{
	const T order = T(1.);
	t_length<Sys, T> lam = T(2.)*pi<T> / k;
	auto dS = order*lam/(T(2.)*d);
	if(std::abs(T(dS)) > T(1))
		throw std::runtime_error("Invalid twotheta angle.");

	t_angle<Sys, T> theta = units::asin(dS);
	if(!bPosSense)
		theta = -theta;
	return theta;
}


/**
 * get k from crystal angle
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 *
 * k = n pi / (d sin(theta))
 */
template<class t_real = double>
t_real calc_tas_k(t_real theta, t_real d)
{
	t_real sintheta = std::abs(std::sin(theta));
	return pi<t_real> / (d * sintheta);
}


/**
 * get k from crystal angle
 * (version with units)
 * @see https://en.wikipedia.org/wiki/Bragg's_law
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T> calc_tas_k(const t_angle<Sys, T>& _theta,
	const t_length<Sys, T>& d, bool bPosSense = true)
{
	t_angle<Sys, T> theta = _theta;
	if(!bPosSense)
		theta = -theta;

	const T order = T(1.);

	// https://en.wikipedia.org/wiki/Bragg%27s_law
	t_length<Sys, T> lam = T(2.)*d/order * units::sin(theta);
	t_wavenumber<Sys, T> k = T(2.)*pi<T> / lam;

	return k;
}


/**
 * get ki from kf and energy transfer
 */
template<class t_real = double>
t_real calc_tas_ki(t_real kf, t_real E)
{
	return std::sqrt(kf*kf + E2KSQ<t_real>*E);
}


/**
 * get kf from ki and energy transfer
 */
template<class t_real = double>
t_real calc_tas_kf(t_real ki, t_real E)
{
	return std::sqrt(ki*ki - E2KSQ<t_real>*E);
}


/**
 * get energy transfer from ki and kf
 */
template<class t_real>
t_real calc_tas_E(t_real ki, t_real kf)
{
	return (ki*ki - kf*kf) / E2KSQ<t_real>;
}


template<class Sys, class T = double>
t_energy<Sys, T> k2E(const t_wavenumber<Sys, T>& k)
{
	T dk = k*angstrom<T>;
	T dE = KSQ2E<T> * dk*dk;
	return dE * meV<T>;
}


template<class Sys, class T = double>
t_wavenumber<Sys, T> E2k(const t_energy<Sys, T>& _E, bool &bImag)
{
	bImag = (_E < T(0.)*meV<T>);
	t_energy<Sys, T> E = bImag ? -_E : _E;
	const T dE = E / meV<T>;
	const T dk = std::sqrt(E2KSQ<T> * dE);
	return dk / angstrom<T>;
}


/**
 * get energy transfer from ki and kf
 * (version with units)
 */
template<class Sys, class T = double>
t_energy<Sys, T> get_energy_transfer(const t_wavenumber<Sys, T>& ki,
	const t_wavenumber<Sys, T>& kf)
{
	return k2E<Sys, T>(ki) - k2E<Sys, T>(kf);
}


/**
 * (hbar*ki)^2 / (2*mn)  -  (hbar*kf)^2 / (2mn)  =  E
 * 1) ki^2  =  +E * 2*mn / hbar^2  +  kf^2
 * 2) kf^2  =  -E * 2*mn / hbar^2  +  ki^2
 */
template<class Sys, class T = double>
t_wavenumber<Sys, T> get_other_k(const t_energy<Sys, T>& E,
	const t_wavenumber<Sys, T>& kfix, bool bFixedKi)
{
	auto kE_sq = E*T(2.)*(m_n<T>/hbar<T>)/hbar<T>;
	if(bFixedKi)
		kE_sq = -kE_sq;

	auto k_sq = kE_sq + kfix*kfix;
	if(k_sq*angstrom<T>*angstrom<T> < T(0.))
		throw std::runtime_error("Scattering triangle not closed.");

	return my_units_sqrt<t_wavenumber<Sys, T>>(k_sq);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * kf^3 mono/ana reflectivity factor
 * @see (Shirane 2002) p. 125
 */
template<class Sys, class T = double>
T ana_effic_factor(const t_wavenumber<Sys, T>& kf, const t_angle<Sys, T>& theta)
{
	return kf*kf*kf / units::tan(theta) * angstrom<T>*angstrom<T>*angstrom<T>;
}


/**
 * kf^3 mono/ana reflectivity factor,
 * @see (Shirane 2002) p. 125
 */
template<class Sys, class T = double>
T ana_effic_factor(const t_wavenumber<Sys, T>& kf, const t_length<Sys, T>& d)
{
	t_angle<Sys, T> theta = units::abs(calc_tas_a1<Sys, T>(kf, d, true));
	return ana_effic_factor<Sys, T>(kf, theta);
}

// --------------------------------------------------------------------------------



/**
 * Bose distribution (occupation number including detailed balance)
 * @see (Shirane 2002), p. 28
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 *
 * bose(+E, T) / bose(-E, T) = [ 1/(exp(E/kT) - 1) + 1 ] / [ 1/(exp(E/kT) - 1) ]
 *                           = 1 + 1/1/(exp(E/kT) - 1)
 *                           = exp(E/kT)
 * which is the detailed balance, S(+Q, +E) / S(-Q, -E), see (Shirane 2002), p. 26.
 */
template<class t_real = double>
t_real bose(t_real E, t_real T)
{
	const t_real _kB = kB<t_real> * kelvin<t_real>/meV<t_real>;

	t_real n = t_real(1)/(std::exp(std::abs(E)/(_kB*T)) - t_real(1));
	if(E >= t_real(0))
		n += t_real(1);

	return n;
}


/**
 * Bose factor with a lower cutoff energy
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class t_real = double>
t_real bose_cutoff(t_real E, t_real T, t_real E_cutoff=t_real(0.02))
{
	t_real dB;

	E_cutoff = std::abs(E_cutoff);
	if(std::abs(E) < E_cutoff)
		dB = bose<t_real>(sign(E)*E_cutoff, T);
	else
		dB = bose<t_real>(E, T);

	return dB;
}


/**
 * Bose factor
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class Sys, class T = double>
T bose(const t_energy<Sys, T>& E, const t_temperature<Sys, T>& temp,
	t_energy<Sys, T> E_cutoff = -meV<T>)
{
	if(E_cutoff < T(0)*meV<T>)
		return bose<T>(T(E/meV<T>), T(temp/kelvin<T>));
	else
		return bose_cutoff<T>(T(E/meV<T>), T(temp/kelvin<T>),
			T(E_cutoff/meV<T>));
}


/**
 * DHO
 * @see B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108, doi: https://doi.org/10.1016/S0921-4526(97)00121-X
 */
template<class t_real = double>
t_real DHO_model(t_real E, t_real T, t_real E0, t_real hwhm, t_real amp, t_real offs)
{
	return std::abs(bose<t_real>(E, T)*amp/(E0*pi<t_real>) *
		(hwhm/((E-E0)*(E-E0) + hwhm*hwhm) - hwhm/((E+E0)*(E+E0) + hwhm*hwhm)))
		+ offs;
}


// --------------------------------------------------------------------------------

/**
 * Fermi distribution
 * @see https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics
 */
template<class t_real=double>
t_real fermi(t_real E, t_real mu, t_real T)
{
	const t_real _kB = kB<t_real> * kelvin<t_real>/meV<t_real>;
	t_real n = t_real(1)/(std::exp((E-mu)/(_kB*T)) + t_real(1));
	return n;
}


/**
 * Fermi distribution
 * @see https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics
 */
template<class Sys, class T = double>
T fermi(const t_energy<Sys, T>& E, const t_energy<Sys, T>& mu,
	const t_temperature<Sys, T>& temp)
{
	return fermi<T>(T(E/meV<T>), T(mu/meV<T>), T(temp/kelvin<T>));
}

// --------------------------------------------------------------------------------


/**
 * get macroscopic from microscopic cross-section
 */
template<class Sys, class T = double>
t_length_inverse<Sys, T> macro_xsect(const t_area<Sys, T>& xsect,
	unsigned int iNumAtoms, const t_volume<Sys, T>& volUC)
{
	return xsect * T(iNumAtoms) / volUC;
}



// --------------------------------------------------------------------------------

/**
 * thin lens equation: 1/f = 1/lenB + 1/lenA
 * @see https://en.wikipedia.org/wiki/Thin_lens
 */
template<class Sys, class T = double>
t_length<Sys, T> focal_len(const t_length<Sys, T>& lenBefore, const t_length<Sys, T>& lenAfter)
{
	const t_length_inverse<Sys, T> f_inv = T(1)/lenBefore + T(1)/lenAfter;
	return T(1) / f_inv;
}


/**
 * optimal mono/ana curvature,
 * @see (Shirane 2002) p. 66
 * @see NICOS: https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/nicos/nicos-core.git/plain/nicos/devices/tas/mono.py
 * @see McStas: https://github.com/McStasMcXtrace/McCode/blob/master/mcstas-comps/optics/Monochromator_curved.comp
 */
template<class Sys, class T = double>
t_length<Sys, T> foc_curv(const t_length<Sys, T>& lenBefore, const t_length<Sys, T>& lenAfter,
	const t_angle<Sys, T>& tt, bool vert_focus)
{
	const t_length<Sys, T> f = focal_len<Sys, T>(lenBefore, lenAfter);
	const T s = T(units::abs(units::sin(T(0.5)*tt)));

	const t_length<Sys, T> curv = vert_focus ? T(2)*f*s : T(2)*f/s;
	return curv;
}


/**
 * optimal mono/ana curvature,
 * @see e.g. (Shirane 2002) p. 66
 * @see e.g. McStas: https://github.com/McStasMcXtrace/McCode/blob/master/mcstas-comps/optics/Monochromator_curved.comp
 * @see e.g. Nicos: https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/nicos/nicos.git/tree/nicos/devices/tas/mono.py
 * @see e.g. [eck14], equs. 59-61
 */
template<class Sys, class T = double>
t_length<Sys, T> foc_curv(const t_length<Sys, T>& lenBefore, const t_length<Sys, T>& lenAfter,
	const t_wavenumber<Sys, T>& k, const t_length<Sys, T>& d, bool vert_focus)
{
	const t_length<Sys, T> f = focal_len<Sys, T>(lenBefore, lenAfter);
	const T s = T(units::abs(pi<T> / d / k));

	const t_length<Sys, T> curv = vert_focus ? T(2)*f*s : T(2)*f/s;
	return curv;
}

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
/**
 * @brief disc chopper burst time
 * @param r chopper radius
 * @param L chopper window length
 * @param om chopper frequency
 * @param bCounterRot single disc or two counter-rotating discs?
 * @param bSigma burst time in sigma or fwhm?
 * @return burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class T = double>
t_time<Sys, T> burst_time(const t_length<Sys, T>& r,
	const t_length<Sys, T>& L, const t_freq<Sys, T>& om, bool bCounterRot,
	bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	T tScale = bCounterRot ? T(2) : T(1);
	return L / (r * om * tScale) * tSig;
}


/**
 * @brief disc chopper burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class T = double>
t_length<Sys, T> burst_time_L(const t_length<Sys, T>& r,
	const t_time<Sys, T>& dt, const t_freq<Sys, T>& om, bool bCounterRot,
	bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	T tScale = bCounterRot ? T(2) : T(1);
	return dt * r * om * tScale / tSig;
}


/**
 * @brief disc chopper burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class T = double>
t_length<Sys, T> burst_time_r(const t_time<Sys, T>& dt,
	const t_length<Sys, T>& L, const t_freq<Sys, T>& om, bool bCounterRot,
	bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	T tScale = bCounterRot ? T(2) : T(1);
	return L / (dt * om * tScale) * tSig;
}


/**
 * @brief disc chopper burst time
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class T = double>
t_freq<Sys, T> burst_time_om(const t_length<Sys, T>& r,
	const t_length<Sys, T>& L, const t_time<Sys, T>& dt, bool bCounterRot,
	bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	T tScale = bCounterRot ? T(2) : T(1);
	return L / (r * dt * tScale) * tSig;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * @brief collimation
 * @param L length of collimator
 * @param w distance between blade
 * @param bSigma calculate sigma or fwhm?
 * @return angular divergence
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class T = double>
t_angle<Sys, T> colli_div(const t_length<Sys, T>& L, const t_length<Sys, T>& w, bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	return units::atan(w/L) * tSig;
}


/**
 * @brief collimation
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class T = double>
t_length<Sys, T> colli_div_L(const t_angle<Sys, T>& ang, const t_length<Sys, T>& w, bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	return w/units::tan(ang/tSig);
}


/**
 * @brief collimation
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class T = double>
t_length<Sys, T> colli_div_w(const t_length<Sys, T>& L, const t_angle<Sys, T>& ang, bool bSigma = true)
{
	const T tSig = bSigma ? FWHM2SIGMA<T> : T(1);
	return units::tan(ang/tSig) * L;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
/**
 * @brief velocity selector
 * @return selector angular frequency
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class T = double>
t_freq<Sys, T> vsel_freq(const t_length<Sys, T>& lam,
	const t_length<Sys, T>& len, const t_angle<Sys, T>& twist)
{
	// https://en.wikiversity.org/wiki/De_Broglie_wavelength
	t_wavenumber<Sys, T> k = T(2.)*pi<T> / lam;
	t_velocity<Sys, T> v_n = hbar<T>*k / m_n<T>;

	return v_n*twist / (len * radian<T>);
}


/**
 * @brief velocity selector
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class T = double>
t_length<Sys, T> vsel_len(const t_length<Sys, T>& lam,
	const t_freq<Sys, T>& om, const t_angle<Sys, T>& twist)
{
	// https://en.wikiversity.org/wiki/De_Broglie_wavelength
	t_wavenumber<Sys, T> k = T(2.)*pi<T> / lam;
	t_velocity<Sys, T> v_n = hbar<T>*k / m_n<T>;

	return v_n*twist / (om * radian<T>);
}


/**
 * @brief velocity selector
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class T = double>
t_angle<Sys, T> vsel_twist(const t_length<Sys, T>& lam,
	const t_freq<Sys, T>& om, const t_length<Sys, T>& len)
{
	// https://en.wikiversity.org/wiki/De_Broglie_wavelength
	t_wavenumber<Sys, T> k = T(2.)*pi<T> / lam;
	t_velocity<Sys, T> v_n = hbar<T>*k / m_n<T>;

	return  (len * om * radian<T>) / v_n;
}


/**
 * @brief velocity selector
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class T = double>
t_length<Sys, T> vsel_lam(const t_angle<Sys, T>& twist,
	const t_freq<Sys, T>& om, const t_length<Sys, T>& len)
{
	t_velocity<Sys, T> v_n = (len * om * radian<T>) / twist;
	t_wavenumber<Sys, T> k = m_n<T>*v_n/hbar<T>;
	t_length<Sys, T> lam = T(2.)*pi<T> / k;

	return lam;
}

// --------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Larmor precession
//------------------------------------------------------------------------------

/**
 * gamma*B = omega
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class T = double>
t_freq<Sys, T> larmor_om(const t_flux<Sys, T>& B)
{
	return co::gamma_n * B;
}


/**
 * B = omega/gamma
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class T = double>
t_flux<Sys, T> larmor_B(const t_freq<Sys, T>& om)
{
	return om/co::gamma_n;
}


/* omega = -gamma*B
 * omega*t = -gamma*B*t
 * phi = - gamma * B * l/v
 * B = -phi*v / (gamma*l)
 * phi = -pi  =>  B = pi*v / (gamma*l)
 *
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class T = double>
t_flux<Sys, T> larmor_field(const t_length<Sys, T>& lam,
	const t_length<Sys, T>& len,
	const t_angle<Sys, T>& phi)
{
	t_velocity<Sys, T> v = h<T> / lam / co::m_n;
	t_freq<Sys, T> om = -T(phi/radians<T>)*v/len;
	return om/co::gamma_n;
}

//------------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// polarisation
// ----------------------------------------------------------------------------

/**
 * polarisation density matrix
 *   (based on a proof from a lecture by P. J. Brown, 2006)
 *
 * eigenvector expansion of a state: |psi> = a_i |xi_i>
 * mean value of operator with mixed states:
 * <A> = p_i * <a_i|A|a_i>
 * <A> = tr( A * p_i * |a_i><a_i| )
 * <A> = tr( A * rho )
 * polarisation density matrix: rho = 0.5 * (1 + <P|sigma>)
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9
 * @see (Desktop Bronstein 2008), Ch. 21 (Zusatzkapitel.pdf), pp. 11-12 and p. 24
 */
template<class t_vec, class t_mat>
t_mat pol_density_mat(const t_vec& P, typename t_vec::value_type c=0.5)
requires is_vec<t_vec> && is_mat<t_mat>
{
	return (unit<t_mat>(2,2) + proj_su2<t_vec, t_mat>(P, true)) * c;
}


/**
 * Blume-Maleev equation
 * calculate equation indirectly with density matrix
 *   (based on a proof from a lecture by P. J. Brown, 2006)
 *
 * V   = N*1 + <Mperp|sigma>
 * I   = tr( <V|V> rho )
 * P_f = tr( <V|sigma|V> rho ) / I
 *
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_mat, class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec>
blume_maleev_indir(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_mat<t_mat> && is_vec<t_vec>
{
	// spin-1/2
	constexpr t_cplx c = 0.5;

	// vector of pauli matrices
	const auto sigma = su2_matrices<std::vector<t_mat>>(false);

	// density matrix
	const auto density = pol_density_mat<t_vec, t_mat>(P_i, c);

	// potential
	const auto V_mag = proj_su2<t_vec, t_mat>(Mperp, true);
	const auto V_nuc = N * unit<t_mat>(2);
	const auto V = V_nuc + V_mag;
	const auto VConj = herm(V);

	// scattering intensity
	t_cplx I = trace(VConj*V * density);

	// ------------------------------------------------------------------------
	// scattered polarisation vector
	const auto m0 = (VConj * sigma[0]) * V * density;
	const auto m1 = (VConj * sigma[1]) * V * density;
	const auto m2 = (VConj * sigma[2]) * V * density;

	t_vec P_f = create<t_vec>({ trace(m0), trace(m1), trace(m2) });
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * Blume-Maleev equation
 * @returns scattering intensity and final polarisation vector
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_vec>
blume_maleev(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_vec<t_vec>
{
	using t_real = typename t_cplx::value_type;

	const t_vec MperpConj = conj(Mperp);
	const t_cplx NConj = std::conj(N);
	constexpr t_cplx imag(0, 1);

	t_cplx N2 = N * NConj;
	t_cplx M2 = inner<t_vec>(Mperp, Mperp);
	t_vec Mx2 = cross<t_vec>({ MperpConj, Mperp });

	// ------------------------------------------------------------------------
	// intensity
	// nuclear and magnetic non-chiral
	t_cplx I = N2 + M2;

	// magnetic chiral
	I += imag * inner<t_vec>(P_i, Mx2);

	// nuclear-magnetic
	t_cplx I_nm = N * inner<t_vec>(Mperp, P_i);
	I += I_nm + std::conj(I_nm);
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// polarisation vector
	// nuclear
	t_vec P_f = N2 * P_i;                           // rotates P

	// magnetic non-chiral
	t_vec rot_ch = Mperp * inner<t_vec>(Mperp, P_i);
	P_f += rot_ch + tl2::conj(rot_ch);              // rotates P
	P_f -= M2 * P_i;                                // rotates P

	// magnetic chiral
	P_f -= imag * Mx2;                              // creates P

	// nuclear-magnetic
	t_vec rot_nm = imag * NConj * cross<t_vec>({ Mperp, P_i });
	t_vec create_nm = NConj * Mperp;
	P_f += rot_nm + tl2::conj(rot_nm);              // rotates P
	P_f += create_nm + tl2::conj(create_nm);        // creates P
	// ------------------------------------------------------------------------

	return std::make_tuple(I, P_f/I);
}


/**
 * Blume-Maleev in tensor form
 *   (based on a lecture by P. J. Brown, 2006, 2009)
 *
 * @see https://doi.org/10.1016/B978-044451050-1/50006-9 - p. 225-226
 */
template<class t_mat, class t_vec, typename t_cplx = typename t_vec::value_type>
std::tuple<t_cplx, t_mat, t_vec, t_vec>
blume_maleev_tensor(const t_vec& P_i, const t_vec& Mperp, const t_cplx& N)
requires is_mat<t_mat> && is_vec<t_vec>
{
	using t_real = typename t_cplx::value_type;

	const t_vec MperpConj = conj(Mperp);
	const t_cplx NConj = std::conj(N);

	t_cplx N2 = N * NConj;
	t_cplx M2 = inner<t_vec>(Mperp, Mperp);
	t_mat Mo = t_real(2) * outer<t_mat, t_vec>(Mperp, Mperp);
	t_vec NM = t_real(2) * N * MperpConj;

	auto [ Mor, Moi ] = split_cplx<t_mat, t_mat>(Mo);
	auto [ NMr, NMi ] = split_cplx<t_vec, t_vec>(NM);

	// cross product vector of imaginary component of Mo
	t_vec Moivec = create<t_vec>({ Moi(2, 1), Moi(0, 2), Moi(1, 0) });

	// ------------------------------------------------------------------------
	// intensity
	// nuclear-magnetic and magnetic chiral intensity
	t_cplx I = inner<t_vec>(P_i, NMr + Moivec);

	// nuclear and magnetic non-chiral intensity
	I += N2 + M2;
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// rotates polarisation
	// nuclear and magnetic non-chiral components
	t_mat Prot = diag<t_mat>({ N2 - M2, N2 - M2, N2 - M2 });

	// magnetic non-chiral component
	// Mor * P_i corresponds to rot_ch + tl2::conj(rot_ch) in the case above
	Prot += Mor;

	// nuclear-magnetic component
	Prot += skewsymmetric<t_mat>(NMi);

	Prot /= I;
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// creates polarisation (nuclear-magnetic and magnetic chiral components)
	// Moivec correspond to imag * Mx2 in the case above
	t_vec Pcreate = NMr - Moivec;
	Pcreate /= I;
	// ------------------------------------------------------------------------

	t_vec P_f = Prot * P_i + Pcreate;
	return std::make_tuple(I, Prot, Pcreate, P_f);
}


/**
 * general structure factor calculation
 * e.g. type T as vector (complex number) for magnetic (nuclear) structure factor
 * Ms_or_bs:
	- nuclear scattering lengths for nuclear neutron scattering or
	- atomic form factors for x-ray scattering
	- magnetisation (* magnetic form factor) for magnetic neutron scattering
 * Rs: atomic positions
 * Q: scattering vector G for nuclear scattering or G+k for magnetic scattering with propagation vector k
 * fs: optional magnetic form factors
 *
 * @see (Shirane 2002), p. 25, equ. 2.26 for nuclear structure factor
 * @see (Shirane 2002), p. 40, equ. 2.81 for magnetic structure factor
 * @see https://doi.org/10.1016/B978-044451050-1/50002-1
 */
template<class t_vec, class T = t_vec, template<class...> class t_cont = std::vector,
	class t_cplx = std::complex<typename t_vec::value_type>>
T structure_factor(const t_cont<T>& Ms_or_bs, const t_cont<t_vec>& Rs,
	const t_vec& Q, const t_vec* fs = nullptr)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_cplx::value_type;
	constexpr t_cplx cI{0,1};
	constexpr t_real twopi = pi<t_real> * t_real{2};
	constexpr t_real expsign = -1;

	T F{};
	if(Rs.size() == 0)
		return F;
	if constexpr(is_vec<T>)
		F = zero<T>(Rs.begin()->size());	// always 3 dims...
	else if constexpr(is_complex<T>)
		F = T(0);

	auto iterM_or_b = Ms_or_bs.begin();
	auto iterR = Rs.begin();
	typename t_vec::const_iterator iterf;
	if(fs) iterf = fs->begin();

	while(iterM_or_b != Ms_or_bs.end() && iterR != Rs.end())
	{
		// if form factors are given, use them, otherwise set to 1
		t_real f = t_real(1);
		if(fs)
		{
			auto fval = *iterf;
			if constexpr(is_complex<decltype(fval)>)
				f = fval.real();
			else
				f = fval;
		}

		// structure factor
		F += (*iterM_or_b) * f * std::exp(expsign * cI * twopi * inner<t_vec>(Q, *iterR));

		// next M or b if available (otherwise keep current)
		auto iterM_or_b_next = std::next(iterM_or_b, 1);
		if(iterM_or_b_next != Ms_or_bs.end())
			iterM_or_b = iterM_or_b_next;

		if(fs)
		{
			// next form factor if available (otherwise keep current)
			auto iterf_next = std::next(iterf, 1);
			if(iterf_next != fs->end())
				iterf = iterf_next;
		}

		// next atomic position
		std::advance(iterR, 1);
	}

	return F;
}
// ----------------------------------------------------------------------------

}
#endif
