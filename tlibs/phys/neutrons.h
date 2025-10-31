/**
 * neutron formulas
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2012-2016
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

#ifndef __TLIBS_NEUTRONS__
#define __TLIBS_NEUTRONS__

#include "units.h"
#include "../math/math.h"
#include "../math/linalg.h"
#include "../helper/exception.h"

#include <boost/units/pow.hpp>
#include <cmath>


namespace tl {


// --------------------------------------------------------------------------------

template<typename T=double> T get_KSQ2E()
{
	const auto _A = get_one_angstrom<T>();
	const auto _meV = get_one_meV<T>();
	const auto _mn = get_m_n<T>();
	const auto _hbar = get_hbar<T>();

	// factor order important for small value types like "float"!
	return T(0.5) * _hbar/_A/_mn * _hbar/_A/_meV;
}


template<typename T=double> T get_E2KSQ()
{
	return T(1)/get_KSQ2E<T>();
}


template<class T=double> T t_KSQ2E = get_KSQ2E<T>();
template<class T=double> T t_E2KSQ = T(1)/t_KSQ2E<T>;


// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// de Broglie: lam = h/p
// @see https://en.wikiversity.org/wiki/De_Broglie_wavelength
// --------------------------------------------------------------------------------

template<class Sys, class Y>
t_momentum<Sys,Y> lam2p(const t_length<Sys,Y>& lam)
{
	return get_h<Y>() / lam;
}


template<class Sys, class Y>
t_length<Sys,Y> p2lam(const t_momentum<Sys,Y>& p)
{
	return get_h<Y>() / p;
}


// lam = 2pi/k
template<class Sys, class Y>
t_length<Sys,Y> k2lam(const t_wavenumber<Sys,Y>& k)
{
	return Y(2.)*get_pi<Y>() / k;
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> lam2k(const t_length<Sys,Y>& lam)
{
	return Y(2.)*get_pi<Y>() / lam;
}


template<class Sys, class Y>
t_momentum<Sys,Y> k2p(const t_wavenumber<Sys,Y>& k)
{
	return get_hbar<Y>()*k;
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> p2k(const t_momentum<Sys,Y>& p)
{
	return p/get_hbar<Y>();
}


template<class Sys, class Y>
t_velocity<Sys,Y> k2v(const t_wavenumber<Sys,Y>& k)
{
	return k2p(k) / get_m_n<Y>();
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> v2k(const t_velocity<Sys,Y>& v)
{
	return get_m_n<Y>()*v/get_hbar<Y>();
}
// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------
// E = hbar*omega
// --------------------------------------------------------------------------------

template<class Sys, class Y>
t_energy<Sys,Y> omega2E(const t_freq<Sys,Y>& omega)
{
	return get_hbar<Y>() * omega;
}


template<class Sys, class Y>
t_freq<Sys,Y> E2omega(const t_energy<Sys,Y>& en)
{
	return en / get_hbar<Y>();
}


template<class Sys, class Y>
t_energy<Sys,Y> k2E_direct(const t_wavenumber<Sys,Y>& k)
{
	t_momentum<Sys,Y> p = get_hbar<Y>()*k;
	t_energy<Sys,Y> E = p*p / (Y(2.)*get_m_n<Y>());
	return E;
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> E2k_direct(const t_energy<Sys,Y>& _E, bool &bImag)
{
	bImag = (_E < Y(0.)*get_one_meV<Y>());
	t_energy<Sys,Y> E = bImag ? -_E : _E;

	auto pp = Y(2.) * get_m_n<Y>() * E;
	t_momentum<Sys,Y> p = my_units_sqrt<t_momentum<Sys,Y>>(pp);
	t_wavenumber<Sys,Y> k = p / get_hbar<Y>();
	return k;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// indirect calculations using conversion factors for numerical stability

template<class Sys, class Y>
t_energy<Sys,Y> k2E(const t_wavenumber<Sys,Y>& k)
{
	Y dk = k*get_one_angstrom<Y>();
	Y dE = get_KSQ2E<Y>() * dk*dk;
	return dE * get_one_meV<Y>();
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> E2k(const t_energy<Sys,Y>& _E, bool &bImag)
{
	bImag = (_E < Y(0.)*get_one_meV<Y>());
	t_energy<Sys,Y> E = bImag ? -_E : _E;
	const Y dE = E / get_one_meV<Y>();
	const Y dk = std::sqrt(get_E2KSQ<Y>() * dE);
	return dk / get_one_angstrom<Y>();
}

// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------
/**
 * Bragg equation
 * real: n * lam = 2d * sin(twotheta/2)
 * @see https://en.wikipedia.org/wiki/Bragg%27s_law
 */
template<class Sys, class Y>
t_length<Sys,Y> bragg_real_lam(const t_length<Sys,Y>& d,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return Y(2.)*d/n * units::sin(twotheta/Y(2.));
}


template<class Sys, class Y>
t_length<Sys,Y> bragg_real_d(const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return n * lam / (Y(2.)*units::sin(twotheta/Y(2.)));
}


template<class Sys, class Y>
t_angle<Sys,Y> bragg_real_twotheta(const t_length<Sys,Y>& d,
	const t_length<Sys,Y>& lam, Y n = Y(1))
{
	auto dS = n*lam/(Y(2.)*d);
	if(std::abs(Y(dS)) > Y(1))
		throw Err("Invalid twotheta angle.");
	return units::asin(dS) * Y(2.);
}


/**
 * reciprocal Bragg equation: G * lam = 4pi * sin(twotheta/2)
 */
template<class Sys, class Y>
t_angle<Sys,Y> bragg_recip_twotheta(const t_wavenumber<Sys,Y>& G,
	const t_length<Sys,Y>& lam, Y n = Y(1))
{
	auto dS = G*n*lam/(Y(4)*get_pi<Y>());
	if(std::abs(Y(dS)) > Y(1))
		throw Err("Invalid twotheta angle.");
	return units::asin(dS) * Y(2);
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> bragg_recip_G(const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return Y(4)*get_pi<Y>() / (n*lam) * units::sin(twotheta/Y(2));
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> bragg_recip_Q(const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return bragg_recip_G<Sys,Y>(lam,twotheta,n);
}


template<class Sys, class Y>
t_length<Sys,Y> bragg_recip_lam(const t_wavenumber<Sys,Y>& G,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return Y(4)*get_pi<Y>() / G * units::sin(twotheta/Y(2)) / n;
}


/**
 * reciprocal Bragg equation [2]: n * G = 2*k * sin(twotheta/2)
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y> bragg_recip_G(const t_wavenumber<Sys,Y>& k,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return Y(2)*k / n * units::sin(twotheta/Y(2));
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> bragg_recip_k(const t_wavenumber<Sys,Y>& G,
	const t_angle<Sys,Y>& twotheta, Y n = Y(1))
{
	return n*G / (Y(2) * units::sin(twotheta/Y(2)));
}


template<class Sys, class Y>
t_angle<Sys,Y> bragg_recip_twotheta(const t_wavenumber<Sys,Y>& G,
	const t_wavenumber<Sys,Y>& k, Y n = Y(1))
{
	auto dS = n * G / (Y(2) * k);
	if(std::abs(Y(dS)) > Y(1))
		throw Err("Invalid twotheta angle.");
	return units::asin(dS) * Y(2);
}


/**
 * lattice vector
 * G = 2pi / d
 */
template<class Sys, class Y>
t_length<Sys,Y> G2d(const t_wavenumber<Sys,Y>& G)
{
	return Y(2.)*get_pi<Y>() / G;
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> d2G(const t_length<Sys,Y>& d)
{
	return Y(2.)*get_pi<Y>() / d;
}

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
template<class Sys, class Y = double>
Y bragg_diff(Y dDoverD, const t_angle<Sys,Y>& theta, Y dTheta)
{
	Y dLamOverLam = dDoverD + units::cos(theta)/units::sin(theta) * dTheta;
	return dLamOverLam;
}

// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------

/**
 * kinematic plane
 * @see e.g. (ILL Neutron Data Booklet), sec. 2.6-2
 *
 * Q_vec = ki_vec - kf_vec
 * Q^2 = ki^2 + kf^2 - 2ki kf cos 2th	| * hbar^2 / (2 mn)
 *
 * using:
 * Ei = hbar^2 ki^2 / (2 mn)
 *
 * Q^2 * hbar^2 / (2 mn) = Ei + Ef - 2 ki kf cos(2th) * hbar^2 / (2 mn)
 *
 * using:
 * ki^2 = 2 mn Ei / hbar^2
 *
 * Q^2 = [Ei + Ef - 2 sqrt(Ei) sqrt(Ef) cos 2th] * 2 mn / hbar^2
 *
 * using:
 * dE = Ei - Ef
 * Ef = Ei - dE
 *
 * Q^2 = [2 Ei - dE - 2 sqrt(Ei (Ei - dE)) cos 2th] * 2 mn / hbar^2
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y> kinematic_plane(bool bFixedKi,
	const t_energy<Sys,Y>& EiEf, const t_energy<Sys,Y>& DeltaE,
	const t_angle<Sys,Y>& twotheta)
{
	t_energy<Sys,Y> dE = DeltaE;
	if(bFixedKi)
		dE = -dE;

	auto c = Y(2.)*get_m_n<Y>() / (get_hbar<Y>()*get_hbar<Y>());
	t_wavenumber<Sys,Y> Q =
		my_units_sqrt<t_wavenumber<Sys,Y>>(c *
			(Y(2.)*EiEf + dE - Y(2.)*units::cos(twotheta) *
			my_units_sqrt<t_wavenumber<Sys,Y>>(EiEf*(EiEf + dE))));

	return Q;
}


/**
 * kinematic plane
 * @see e.g. (ILL Neutron Data Booklet), sec. 2.6-2
 *
 * solving the above equation for dE using sage:
 *   Q, Ei, dE, ctt, c = var("Q, Ei, dE, ctt, c")
 *   equ = (Q^2 -2*Ei*c + dE*c)^2 == 4*Ei*(Ei-dE)*c^2*ctt^2
 *   equ.solve(dE)
 */
template<class Sys, class Y>
t_energy<Sys,Y> kinematic_plane(bool bFixedKi, bool bBranch,
	const t_energy<Sys,Y>& EiEf, const t_wavenumber<Sys,Y>& Q,
	const t_angle<Sys,Y>& twotheta)
{
	using t_cE = units::quantity<units::unit<typename units::derived_dimension<
		units::length_base_dimension, -2>::type,
		Sys>, Y>;

	auto c = Y(2.)*get_m_n<Y>() / (get_hbar<Y>()*get_hbar<Y>());
	auto c2 = c*c;
	auto EiEf2 = EiEf*EiEf;

	Y ctt = units::cos(twotheta);
	Y ctt2 = ctt*ctt;

	Y dSign = bBranch ? Y(1.) : Y(-1.);
	Y dSignFixedKf = bFixedKi ? Y(-1.) : Y(1.);

	t_energy<Sys,Y> dE =
		dSignFixedKf*Y(2.) * EiEf * ctt2
		- dSignFixedKf*Y(2.) * EiEf
		+ dSignFixedKf * Q*Q / c
		+ dSign*Y(2.) * ctt/c * my_units_sqrt<t_cE>(
			c2*ctt2*EiEf2 - c2*EiEf2 + c*EiEf*Q*Q);

	return dE;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// Debye-Waller factor
// @see e.g. (Shirane 2002) p. 24, (Squires 2012) p. 34-35
// --------------------------------------------------------------------------------

template<class Sys, class Y>
Y debye_waller_high_T(const t_temperature<Sys,Y>& T_D,
	const t_temperature<Sys,Y>& T, const t_mass<Sys,Y>& M,
	const t_wavenumber<Sys,Y>& Q, t_length_square<Sys,Y>* pZeta_sq = nullptr)
{
	t_length_square<Sys,Y> zeta_sq;
	zeta_sq = Y(9.)*get_hbar<Y>()/get_kB<Y>() / (T_D * M) * T/T_D * get_hbar<Y>();
	Y dwf = units::exp(Y(-1./3.) * Q*Q * zeta_sq);

	if(pZeta_sq)
		*pZeta_sq = zeta_sq;
	return dwf;
}


template<class Sys, class Y>
Y debye_waller_low_T(const t_temperature<Sys,Y>& T_D,
	const t_temperature<Sys,Y>& T, const t_mass<Sys,Y>& M,
	const t_wavenumber<Sys,Y>& Q, t_length_square<Sys,Y>* pZeta_sq = nullptr)
{
	t_length_square<Sys,Y> zeta_sq;
	zeta_sq = Y(9.)*get_hbar<Y>()/get_kB<Y>() / (Y(4.)*T_D*M) * get_hbar<Y>() *
		(Y(1.) + Y(2./3.) * get_pi<Y>()*get_pi<Y>() * (T/T_D)*(T/T_D));
	Y dwf = units::exp(Y(-1./3.) * Q*Q * zeta_sq);

	if(pZeta_sq)
		*pZeta_sq = zeta_sq;
	return dwf;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// scattering triangle / TAS stuff
// @see (Shirane 2002), Ch. 1.3
// --------------------------------------------------------------------------------

/**
 * Q_vec = ki_vec - kf_vec
 * kf_vec = ki_vec - Q_vec
 * kf^2 = ki^2 + Q^2 - 2ki Q cos th
 * cos th = (-kf^2 + ki^2 + Q^2) / (2kiQ)
 */
template<class Sys, class Y>
t_angle<Sys,Y> get_angle_ki_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_wavenumber<Sys,Y>& Q,
	bool bPosSense = true, bool bAngleOutsideTriag = false)
{
	t_angle<Sys,Y> angle;

	if(Q*get_one_angstrom<Y>() == Y(0.))
		angle = get_pi<Y>()/Y(2) * get_one_radian<Y>();
	else
	{
		auto c = (ki*ki - kf*kf + Q*Q) / (Y(2.)*ki*Q);
		if(units::abs(c) > Y(1.))
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(bAngleOutsideTriag) angle = get_pi<Y>()*get_one_radian<Y>() - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}


/**
 * Q_vec = ki_vec - kf_vec
 * ki_vec = Q_vec + kf_vec
 * ki^2 = Q^2 + kf^2 + 2Q kf cos th
 * cos th = (ki^2 - Q^2 - kf^2) / (2Q kf)
 */
template<class Sys, class Y>
t_angle<Sys,Y> get_angle_kf_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_wavenumber<Sys,Y>& Q,
	bool bPosSense = true, bool bAngleOutsideTriag = true)
{
	t_angle<Sys,Y> angle;

	if(Q*get_one_angstrom<Y>() == Y(0.))
		angle = get_pi<Y>()/Y(2) * get_one_radian<Y>();
	else
	{
		auto c = (ki*ki - kf*kf - Q*Q) / (Y(2.)*kf*Q);
		if(units::abs(c) > Y(1.))
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(!bAngleOutsideTriag) angle = get_pi<Y>()*get_one_radian<Y>() - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}


template<class Sys, class Y>
t_angle<Sys,Y> get_mono_twotheta(const t_wavenumber<Sys,Y>& k,
	const t_length<Sys,Y>& d, bool bPosSense = true)
{
	const Y dOrder = Y(1.);
	t_angle<Sys,Y> tt = bragg_real_twotheta(d, k2lam(k), dOrder);
	if(!bPosSense)
		tt = -tt;
	return tt;
}


template<class Sys, class Y>
t_wavenumber<Sys,Y> get_mono_k(const t_angle<Sys,Y>& _theta,
	const t_length<Sys,Y>& d, bool bPosSense = true)
{
	t_angle<Sys,Y> theta = _theta;
	if(!bPosSense)
		theta = -theta;

	const Y dOrder = Y(1.);
	return lam2k(bragg_real_lam(d, Y(2.)*theta, dOrder));
}


/**
 * Q_vec = ki_vec - kf_vec
 * Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
 * cos 2th = (-Q^2 + ki^2 + kf^2) / (2ki kf)
 */
template<class Sys, class Y>
t_angle<Sys,Y> get_sample_twotheta(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_wavenumber<Sys,Y>& Q,
	bool bPosSense = true)
{
	t_dimensionless<Sys,Y> ttCos = (ki*ki + kf*kf - Q*Q)/(Y(2.)*ki*kf);
	if(units::abs(ttCos) > Y(1.))
		throw Err("Scattering triangle not closed.");

	t_angle<Sys,Y> tt;
	tt = units::acos(ttCos);

	if(!bPosSense)
		tt = -tt;
	return tt;
}


/**
 * again cos theorem:
 * Q_vec = ki_vec - kf_vec
 * Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
 * Q = sqrt(ki^2 + kf^2 - 2ki kf cos 2th)
 */
template<class Sys, class Y>
const t_wavenumber<Sys,Y>
get_sample_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_angle<Sys,Y>& tt)
{
	t_dimensionless<Sys,Y> ctt = units::cos(tt);
	decltype(ki*ki) Qsq = ki*ki + kf*kf - Y(2.)*ki*kf*ctt;
	if(Y(Qsq*get_one_angstrom<Y>()*get_one_angstrom<Y>()) < Y(0.))
	{
		// TODO

		Qsq = -Qsq;
	}

	t_wavenumber<Sys,Y> Q = my_units_sqrt<t_wavenumber<Sys,Y>>(Qsq);
	return Q;
}



template<class Sys, class Y>
t_energy<Sys,Y> get_energy_transfer(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf)
{
	return k2E<Sys,Y>(ki) - k2E<Sys,Y>(kf);
}


/**
 * (hbar*ki)^2 / (2*mn)  -  (hbar*kf)^2 / (2mn)  =  E
 * 1) ki^2  =  +E * 2*mn / hbar^2  +  kf^2
 * 2) kf^2  =  -E * 2*mn / hbar^2  +  ki^2
 */
template<class Sys, class Y>
t_wavenumber<Sys,Y> get_other_k(const t_energy<Sys,Y>& E,
	const t_wavenumber<Sys,Y>& kfix, bool bFixedKi)
{
	auto kE_sq = E*Y(2.)*(get_m_n<Y>()/get_hbar<Y>())/get_hbar<Y>();
	if(bFixedKi) kE_sq = -kE_sq;

	auto k_sq = kE_sq + kfix*kfix;
	if(k_sq*get_one_angstrom<Y>()*get_one_angstrom<Y>() < Y(0.))
		throw Err("Scattering triangle not closed.");

	return my_units_sqrt<t_wavenumber<Sys,Y>>(k_sq);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------

/**
 * kf^3 mono/ana reflectivity factor
 * @see e.g. (Shirane 2002) p. 125
 */
template<class Sys, class Y>
Y ana_effic_factor(const t_wavenumber<Sys, Y>& kf, const t_angle<Sys, Y>& theta)
{
	return kf*kf*kf / units::tan(theta) *
		get_one_angstrom<Y>()*get_one_angstrom<Y>()*get_one_angstrom<Y>();
}


/**
 * kf^3 mono/ana reflectivity factor
 * @see e.g. (Shirane 2002) p. 125
 */
template<class Sys, class Y>
Y ana_effic_factor(const t_wavenumber<Sys, Y>& kf, const t_length<Sys, Y>& d)
{
	t_angle<Sys, Y> theta = Y(0.5)*units::abs(get_mono_twotheta<Sys, Y>(kf, d, true));
	return ana_effic_factor<Sys, Y>(kf, theta);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// spurions

/**
 * Bragg tail
 * @see (Shirane 2002) p. 152
 */
template<class Sys, class Y>
t_energy<Sys,Y> get_bragg_tail(t_wavenumber<Sys,Y> k,
	t_wavenumber<Sys,Y> q, bool bConstEi=0)
{
	Y t = q / (Y(2.)*k);
	if(!bConstEi)
		t = -t;
	t_energy<Sys,Y> E = get_hbar<Y>()/get_m_n<Y>() * k*q*(Y(1.)+t) * get_hbar<Y>();
	return E;
}


/**
 * higher-order inelastic spurions
 * @see (Shirane 2002) pp. 146-148
 */
template<class Sys, class Y>
t_energy<Sys,Y> get_inelastic_spurion(bool bConstEi, t_energy<Sys,Y> E,
	unsigned int iOrderMono, unsigned int iOrderAna)
{
	const Y dOrderMonoSq = Y(iOrderMono)*Y(iOrderMono);
	const Y dOrderAnaSq = Y(iOrderAna)*Y(iOrderAna);

	t_energy<Sys,Y> E_sp;

	// formulas from (Shirane 2002), p. 147
	if(bConstEi)
		E_sp = (Y(1.) - dOrderMonoSq/dOrderAnaSq) * E;
	else
		E_sp = (dOrderAnaSq/dOrderMonoSq - Y(1.)) * E;

	return E_sp;
}


template<class Y=double>
struct InelasticSpurion
{
	Y dE_meV = Y(0.);
	unsigned int iOrderMono = 1;
	unsigned int iOrderAna = 1;
};


template<class Sys, class Y>
std::vector<InelasticSpurion<Y>> check_inelastic_spurions(bool bConstEi,
	t_energy<Sys,Y> Ei, t_energy<Sys,Y> Ef,
	t_energy<Sys,Y> E, unsigned int iMaxOrder = 5)
{
	const Y dESensitivity = Y(0.25);	// meV

	std::vector<InelasticSpurion<Y>> vecSpuris;

	for(unsigned int iOrder=1; iOrder<=iMaxOrder; ++iOrder)
	{
		InelasticSpurion<Y> spuri;
		t_energy<Sys,Y> EiEf;

		if(bConstEi)
		{
			spuri.iOrderAna = iOrder;
			EiEf = Ei;
		}
		else
		{
			spuri.iOrderMono = iOrder;
			EiEf = Ef;
		}

		spuri.dE_meV = get_inelastic_spurion(bConstEi, EiEf,
			spuri.iOrderMono, spuri.iOrderAna) / get_one_meV<Y>();

		if(spuri.dE_meV!=Y(0.) && float_equal<Y>(spuri.dE_meV, Y(E/get_one_meV<Y>()), dESensitivity))
			vecSpuris.push_back(spuri);
	}

	return vecSpuris;
}


struct ElasticSpurion
{
	bool bAType = false;
	bool bMType = false;

	bool bAKfSmallerKi = false;
	bool bMKfSmallerKi = false;
};


/**
 * accidental elastic (currat-axe) spurions
 * @see (Shirane 2002) pp. 150-155 (esp. fig. 6.2)
 */
template<typename T=double>
ElasticSpurion check_elastic_spurion(const ublas::vector<T>& ki,
	const ublas::vector<T>& kf, const ublas::vector<T>& q)
{
	const T dKi = veclen(ki);
	const T dKf = veclen(kf);
	const T dq = veclen(q);

	const T dAngleSensitivity = T(2.);
	const T dQSensitivity = std::max(dKi, dKf) / T(50.);


	ElasticSpurion result;

	ublas::vector<T> ki_norm = ki;	ki_norm /= dKi;
	ublas::vector<T> kf_norm = kf;	kf_norm /= dKf;

	// Q, q and G point in the opposite direction in Shirane!
	// Shirane: Q = kf - ki, E = Ei - Ef
	// here: Q = ki - kf, E = Ei - Ef
	ublas::vector<T> q_norm = -q;	q_norm /= dq;

	T dAngleKfq = std::acos(inner(kf_norm, q_norm));
	T dAngleKiq = std::acos(inner(ki_norm, q_norm));

	bool bKiqParallel = false, bkiqAntiParallel = false;
	bool bKfqParallel = false, bKfqAntiParallel = false;

	if(float_equal<T>(dAngleKiq, 0., d2r(dAngleSensitivity)))
		bKiqParallel = true;
	else if(float_equal<T>(dAngleKiq, get_pi<T>(), d2r(dAngleSensitivity)))
		bkiqAntiParallel = true;
	if(float_equal<T>(dAngleKfq, 0., d2r(dAngleSensitivity)))
		bKfqParallel = true;
	else if(float_equal<T>(dAngleKfq, get_pi<T>(), d2r(dAngleSensitivity)))
		bKfqAntiParallel = true;

	// type A: q || kf, kf > ki
	if(bKfqParallel)
	{
		T dApparentKf = dKf - dq;

		if(float_equal<T>(dApparentKf, dKi, dQSensitivity))
		{
			result.bAType = true;
			result.bAKfSmallerKi = false;
		}
	}
	// type A: q || kf, kf < ki
	else if(bKfqAntiParallel)
	{
		T dApparentKf = dKf + dq;

		if(float_equal<T>(dApparentKf, dKi, dQSensitivity))
		{
			result.bAType = true;
			result.bAKfSmallerKi = true;
		}
	}

	// type M: q || ki, kf > ki
	if(bKiqParallel)
	{
		T dApparentKi = dKi + dq;

		if(float_equal<T>(dApparentKi, dKf, dQSensitivity))
		{
			result.bMType = true;
			result.bMKfSmallerKi = false;
		}
	}
	// type M: q || ki, kf < ki
	else if(bkiqAntiParallel)
	{
		T dApparentKi = dKi - dq;

		if(float_equal<T>(dApparentKi, dKf, dQSensitivity))
		{
			result.bMType = true;
			result.bMKfSmallerKi = true;
		}
	}

	return result;
}


// --------------------------------------------------------------------------------

/**
 * Bose distribution
 * see e.g.: (Shirane 2002), p. 28
 * @see https://en.wikipedia.org/wiki/Bose%E2%80%93Einstein_statistics
 */
template<class t_real=double>
t_real bose(t_real E, t_real T)
{
	const t_real kB = get_kB<t_real>() * get_one_kelvin<t_real>()/get_one_meV<t_real>();

	t_real n = t_real(1)/(std::exp(std::abs(E)/(kB*T)) - t_real(1));
	if(E >= t_real(0))
		n += t_real(1);

	return n;
}


/**
 * Bose factor with a lower cutoff energy
 */
template<class t_real=double>
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


template<class Sys, class Y>
Y bose(const t_energy<Sys,Y>& E, const t_temperature<Sys,Y>& T,
	t_energy<Sys,Y> E_cutoff = -get_one_meV<Y>())
{
	if(E_cutoff < Y(0)*get_one_meV<Y>())
		return bose<Y>(Y(E/get_one_meV<Y>()), Y(T/kelvin));
	else
		return bose_cutoff<Y>(Y(E/get_one_meV<Y>()), Y(T/kelvin),
			Y(E_cutoff/get_one_meV<Y>()));
}


/**
 * @see B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108, doi: https://doi.org/10.1016/S0921-4526(97)00121-X
 */
template<class t_real=double>
t_real DHO_model(t_real E, t_real T, t_real E0, t_real hwhm, t_real amp, t_real offs)
{
	return std::abs(bose<t_real>(E, T)*amp/(E0*get_pi<t_real>()) *
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
	const t_real kB = get_kB<t_real>() * get_one_kelvin<t_real>()/get_one_meV<t_real>();
	t_real n = t_real(1)/(std::exp((E-mu)/(kB*T)) + t_real(1));
	return n;
}


template<class Sys, class Y>
Y fermi(const t_energy<Sys,Y>& E, const t_energy<Sys,Y>& mu,
	const t_temperature<Sys,Y>& T)
{
	return fermi<Y>(Y(E/get_one_meV<Y>()), Y(mu/get_one_meV<Y>()),
		Y(T/kelvin));
}

// --------------------------------------------------------------------------------


/**
 * get macroscopic from microscopic cross-section
 */
template<class Sys, class Y=double>
t_length_inverse<Sys, Y> macro_xsect(const t_area<Sys, Y>& xsect,
	unsigned int iNumAtoms, const t_volume<Sys, Y>& volUC)
{
	return xsect * Y(iNumAtoms) / volUC;
}



// --------------------------------------------------------------------------------

/**
 * thin lens equation: 1/f = 1/lenB + 1/lenA
 * @see https://en.wikipedia.org/wiki/Thin_lens
 */
template<class Sys, class Y=double>
t_length<Sys, Y> focal_len(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter)
{
	const t_length_inverse<Sys, Y> f_inv = Y(1)/lenBefore + Y(1)/lenAfter;
	return Y(1) / f_inv;
}


/**
 * optimal mono/ana curvature,
 * @see e.g. (Shirane 2002) p. 66
 * @see e.g. McStas: https://github.com/McStasMcXtrace/McCode/blob/master/mcstas-comps/optics/Monochromator_curved.comp
 * @see e.g. Nicos: https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/nicos/nicos.git/tree/nicos/devices/tas/mono.py
 * @see e.g. [eck14], equs. 59-61
 */
template<class Sys, class Y=double>
t_length<Sys, Y> foc_curv(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter,
	const t_angle<Sys, Y>& tt, bool bVert)
{
	const t_length<Sys, Y> f = focal_len<Sys, Y>(lenBefore, lenAfter);
	const Y s = Y(units::abs(units::sin(Y(0.5)*tt)));

	const t_length<Sys, Y> curv = bVert ? Y(2)*f*s : Y(2)*f/s;
	return curv;
}


/**
 * optimal mono/ana curvature,
 * @see e.g. (Shirane 2002) p. 66
 * @see e.g. McStas: https://github.com/McStasMcXtrace/McCode/blob/master/mcstas-comps/optics/Monochromator_curved.comp
 * @see e.g. Nicos: https://forge.frm2.tum.de/cgit/cgit.cgi/frm2/nicos/nicos.git/tree/nicos/devices/tas/mono.py
 * @see e.g. [eck14], equs. 59-61
 */
template<class Sys, class Y=double>
t_length<Sys, Y> foc_curv(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter,
	const t_wavenumber<Sys,Y>& k, const t_length<Sys,Y>& d, bool bVert)
{
	const t_length<Sys, Y> f = focal_len<Sys, Y>(lenBefore, lenAfter);
	const Y s = Y(units::abs(get_pi<Y>() / d / k));

	const t_length<Sys, Y> curv = bVert ? Y(2)*f*s : Y(2)*f/s;
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
 *
 * @see NIMA 492, pp. 97-104 (2002), doi: https://doi.org/10.1016/S0168-9002(02)01285-8
 */
template<class Sys, class Y=double>
t_time<Sys,Y> burst_time(const t_length<Sys,Y>& r,
	const t_length<Sys,Y>& L, const t_freq<Sys,Y>& om, bool bCounterRot,
	bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return L / (r * om * tScale) * tSig;
}


template<class Sys, class Y=double>
t_length<Sys,Y> burst_time_L(const t_length<Sys,Y>& r,
	const t_time<Sys,Y>& dt, const t_freq<Sys,Y>& om, bool bCounterRot,
	bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return dt * r * om * tScale / tSig;
}


template<class Sys, class Y=double>
t_length<Sys,Y> burst_time_r(const t_time<Sys,Y>& dt,
	const t_length<Sys,Y>& L, const t_freq<Sys,Y>& om, bool bCounterRot,
	bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
	return L / (dt * om * tScale) * tSig;
}


template<class Sys, class Y=double>
t_freq<Sys,Y> burst_time_om(const t_length<Sys,Y>& r,
	const t_length<Sys,Y>& L, const t_time<Sys,Y>& dt, bool bCounterRot,
	bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	Y tScale = bCounterRot ? Y(2) : Y(1);
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
 *
 * @see (Shirane 2002), Ch. 3.3
 */
template<class Sys, class Y=double>
t_angle<Sys,Y> colli_div(const t_length<Sys,Y>& L, const t_length<Sys,Y>& w, bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	return units::atan(w/L) * tSig;
}

template<class Sys, class Y=double>
t_length<Sys,Y> colli_div_L(const t_angle<Sys,Y>& ang, const t_length<Sys,Y>& w, bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	return w/units::tan(ang/tSig);
}

template<class Sys, class Y=double>
t_length<Sys,Y> colli_div_w(const t_length<Sys,Y>& L, const t_angle<Sys,Y>& ang, bool bSigma = true)
{
	const Y tSig = bSigma ? get_FWHM2SIGMA<Y>() : Y(1);
	return units::tan(ang/tSig) * L;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
/**
 * @brief velocity selector
 * @return selector angular frequency
 * @see https://doi.org/10.1016/0921-4526(95)00336-8
 */
template<class Sys, class Y=double>
t_freq<Sys, Y> vsel_freq(const t_length<Sys,Y>& lam,
	const t_length<Sys,Y>& len, const t_angle<Sys,Y>& twist)
{
	t_velocity<Sys,Y> v_n = k2v<Sys,Y>(lam2k<Sys,Y>(lam));
	return v_n*twist / (len * get_one_radian<Y>());
}


template<class Sys, class Y=double>
t_length<Sys,Y> vsel_len(const t_length<Sys,Y>& lam,
	const t_freq<Sys, Y>& om, const t_angle<Sys,Y>& twist)
{
	t_velocity<Sys,Y> v_n = k2v<Sys,Y>(lam2k<Sys,Y>(lam));
	return v_n*twist / (om * get_one_radian<Y>());
}


template<class Sys, class Y=double>
t_angle<Sys,Y> vsel_twist(const t_length<Sys,Y>& lam,
	const t_freq<Sys, Y>& om, const t_length<Sys,Y>& len)
{
	t_velocity<Sys,Y> v_n = k2v<Sys,Y>(lam2k<Sys,Y>(lam));
	return  (len * om * get_one_radian<Y>()) / v_n;
}


template<class Sys, class Y=double>
t_length<Sys,Y> vsel_lam(const t_angle<Sys,Y>& twist,
	const t_freq<Sys, Y>& om, const t_length<Sys,Y>& len)
{
	t_velocity<Sys,Y> v_n = (len * om * get_one_radian<Y>()) / twist;
	return k2lam<Sys,Y>(v2k<Sys,Y>(v_n));
}

// --------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Larmor precession
//------------------------------------------------------------------------------

/**
 * gamma*B = omega
 * @see https://en.wikipedia.org/wiki/Larmor_precession
 */
template<class Sys, class Y=double>
t_freq<Sys,Y> larmor_om(const t_flux<Sys,Y>& B)
{
	return co::gamma_n * B;
}


template<class Sys, class Y=double>
t_flux<Sys,Y> larmor_B(const t_freq<Sys,Y>& om)
{
	return om/co::gamma_n;
}


/* omega = -gamma*B
 * omega*t = -gamma*B*t
 * phi = - gamma * B * l/v
 * B = -phi*v / (gamma*l)
 * phi = -pi  =>  B = pi*v / (gamma*l)
 */
template<class Sys, class Y=double>
t_flux<Sys,Y> larmor_field(const t_length<Sys,Y>& lam,
	const t_length<Sys,Y>& len,
	const t_angle<Sys,Y>& phi)
{
	t_velocity<Sys,Y> v = lam2p(lam) / co::m_n;
	t_freq<Sys,Y> om = -Y(phi/radians)*v/len;
	return om/co::gamma_n;
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
/**
 * normalise detector (y) by monitor (m) counts
 * y_new = y / m
 * dy_new = 1/m dy - y/m^2 dm
 */
template<class t_real = double>
std::pair<t_real, t_real> norm_cnts_to_mon(t_real y, t_real dy, t_real m, t_real dm)
{
	t_real val = y / m;
	t_real err = std::sqrt(std::pow(dy/m, 2.) + std::pow(dm*y/(m*m), 2.));

	return std::make_pair(val, err);
}
//------------------------------------------------------------------------------

}
#endif
