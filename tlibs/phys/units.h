/**
 * wrapper for boost.units
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS_UNITS__
#define __TLIBS_UNITS__

#include "../helper/boost_hacks.h"

#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimensionless_quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/io.hpp>

#include <boost/units/systems/si.hpp>
#include <boost/units/systems/angle/degrees.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/neutron_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>

#include <boost/numeric/ublas/vector.hpp>

#include "../math/math.h"


namespace tl {

namespace units = boost::units;
namespace co = boost::units::si::constants::codata;


// general quantities
template<class Sys, class T = double> using t_length =
	units::quantity<units::unit<units::length_dimension, Sys>, T>;
template<class Sys, class T = double> using t_momentum =
	units::quantity<units::unit<units::momentum_dimension, Sys>, T>;
template<class Sys, class T = double> using t_wavenumber =
	units::quantity<units::unit<units::wavenumber_dimension, Sys>, T>;
template<class Sys, class T = double> using t_velocity =
	units::quantity<units::unit<units::velocity_dimension, Sys>, T>;
template<class Sys, class T = double> using t_frequency =
	units::quantity<units::unit<units::frequency_dimension, Sys>, T>;
template<class Sys, class T = double> using t_energy =
	units::quantity<units::unit<units::energy_dimension, Sys>, T>;
template<class Sys, class T = double> using t_angle =
	units::quantity<units::unit<units::plane_angle_dimension, Sys>, T>;
template<class Sys, class T = double> using t_temperature =
	units::quantity<units::unit<units::temperature_dimension, Sys>, T>;
template<class Sys, class T = double> using t_mass =
	units::quantity<units::unit<units::mass_dimension, Sys>, T>;
template<class Sys, class T = double> using t_time =
	units::quantity<units::unit<units::time_dimension, Sys>, T>;
template<class Sys, class T = double> using t_flux =
	units::quantity<units::unit<units::magnetic_flux_density_dimension, Sys>, T>;
template<class Sys, class T = double> using t_inductance =
	units::quantity<units::unit<units::inductance_dimension, Sys>, T>;
template<class Sys, class T = double> using t_area =
	units::quantity<units::unit<units::area_dimension, Sys>, T>;
template<class Sys, class T = double> using t_volume =
	units::quantity<units::unit<units::volume_dimension, Sys>, T>;

template<class Sys, class T = double> using t_length_inverse =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, -1>::type, Sys>, T>;
template<class Sys, class T = double> using t_length_square =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, 2>::type, Sys>, T>;
template<class Sys, class T = double> using t_momentum_square =
	units::quantity<units::unit<units::derived_dimension<units::momentum_dimension, 2>::type, Sys>, T>;
template<class Sys, class T = double> using t_action =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2, units::time_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T = double> using t_energy_per_temperature =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2,
	units::time_base_dimension,-2, units::temperature_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T = double> using t_energy_per_field =
	units::quantity<units::unit<typename units::derived_dimension
	<units::current_base_dimension,1, units::length_base_dimension,2>::type, Sys>, T>;
template<class Sys, class T = double> using t_inductance_per_length =
	units::quantity<units::unit<typename units::derived_dimension
	<units::current_base_dimension,-2, units::length_base_dimension,1,
	units::time_base_dimension,-2, units::mass_base_dimension,1>::type, Sys>, T>;
template<class Sys, class T = double> using t_inv_flux_time =
	units::quantity<units::unit<typename units::derived_dimension
	<units::current_base_dimension,1, units::time_base_dimension,1, units::mass_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T = double> using t_dimensionless =
	units::quantity<units::unit<units::dimensionless_type, Sys>, T>;


// synonyms
template<class Sys, class T = double> using t_freq = t_frequency<Sys, T>;
template<class Sys, class T = double> using t_temp = t_temperature<Sys, T>;


// si quantities -- partial specialisations
template<class T = double> using t_length_si = t_length<units::si::system, T>;
template<class T = double> using t_length_inverse_si = t_length_inverse<units::si::system, T>;
template<class T = double> using t_volume_si = t_volume<units::si::system, T>;
template<class T = double> using t_momentum_si = t_momentum<units::si::system, T>;
template<class T = double> using t_wavenumber_si = t_wavenumber<units::si::system, T>;
template<class T = double> using t_velocity_si = t_velocity<units::si::system, T>;
template<class T = double> using t_frequency_si = t_frequency<units::si::system, T>;
template<class T = double> using t_energy_si = t_energy<units::si::system, T>;
template<class T = double> using t_angle_si = t_angle<units::si::system, T>;
template<class T = double> using t_temperature_si = t_temperature<units::si::system, T>;
template<class T = double> using t_mass_si = t_mass<units::si::system, T>;
template<class T = double> using t_time_si = t_time<units::si::system, T>;
template<class T = double> using t_flux_si = t_flux<units::si::system, T>;
template<class T = double> using t_inductance_si = t_inductance<units::si::system, T>;
template<class T = double> using t_area_si = t_area<units::si::system, T>;
template<class T = double> using t_action_si = t_action<units::si::system, T>;
template<class T = double> using t_energy_per_temperature_si = t_energy_per_temperature<units::si::system, T>;
template<class T = double> using t_inv_flux_time_si = t_inv_flux_time<units::si::system, T>;


// si quantities -- full specialisations
using length = t_length_si<>;
using inv_length = decltype(1./length());
using momentum = t_momentum_si<>;
using wavenumber = t_wavenumber_si<>;
using velocity = t_velocity_si<>;
using frequency = t_frequency_si<>;
using energy = t_energy_si<>;
using angle = t_angle_si<>;
using temperature = t_temperature_si<>;
using mass = t_mass_si<>;
using time = t_time_si<>;
using flux = t_flux_si<>;
using area = t_area_si<>;
using action = t_action_si<>;


// synonyms
typedef frequency freq;
typedef temperature temp;


// constants
template<class T = double> t_energy<units::si::system, T> get_one_meV()
{ return T(1e-3) * T(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class T = double> t_energy<units::si::system, T> get_one_eV()
{ return T(1) * T(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class T = double> t_energy<units::si::system, T> get_one_MeV()
{ return T(1e6) * T(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class T = double> t_length<units::si::system, T> get_one_angstrom()
{ return T(1e-10) * units::si::meters; }
template<class T = double> t_length<units::si::system, T> get_one_meter()
{ return T(1) * units::si::meters; }
template<class T = double> t_length<units::si::system, T> get_one_femtometer()
{ return T(1e-15) * units::si::meters; }
template<class T = double> t_mass<units::si::system, T> get_one_kg()
{ return T(1) * units::si::kilograms; }
template<class T = double> t_area<units::si::system, T> get_one_barn()
{ return T(1e-28) * units::si::meters*units::si::meters; }
template<class T = double> t_temperature<units::si::system, T> get_one_kelvin()
{ return T(1) * units::si::kelvins; }
template<class T = double> t_length<units::si::system, T> get_one_centimeter()
{ return T(1e-2) * units::si::meters; }
template<class T = double> t_time<units::si::system, T> get_one_second()
{ return T(1) * units::si::seconds; }
template<class T = double> t_time<units::si::system, T> get_one_picosecond()
{ return T(1e-12) * units::si::seconds; }
template<class T = double> t_angle<units::si::system, T> get_one_radian()
{ return T(1) * units::si::radians; }
template<class T = double> t_angle<units::si::system, T> get_one_deg()
{ return T(get_pi<T>()/T(180.)) * units::si::radians; }
template<class T = double> t_flux<units::si::system, T> get_one_tesla()
{ return T(1) * units::si::teslas; }
template<class T = double> t_flux<units::si::system, T> get_one_kilogauss()
{ return T(0.1) * units::si::teslas; }
template<class T = double> t_inductance<units::si::system, T> get_one_henry()
{ return T(1) * units::si::henry; }

template<class T = double> t_mass<units::si::system, T> get_m_n()
{ return T(co::m_n/units::si::kilograms)*units::si::kilograms; }
template<class T = double> t_mass<units::si::system, T> get_m_e()
{ return T(co::m_e/units::si::kilograms)*units::si::kilograms; }
template<class T = double> t_mass<units::si::system, T> get_amu()
{ return T(co::m_u/units::si::kilograms)*units::si::kilograms; }
template<class T = double> t_action<units::si::system, T> get_hbar()
{ return T(co::hbar/units::si::joules/units::si::seconds)*units::si::joules*units::si::seconds; }
template<class T = double> t_action<units::si::system, T> get_h()
{ return get_hbar<T>() * T(2)*get_pi<T>(); }
template<class T = double> t_velocity<units::si::system, T> get_c()
{ return T(co::c/units::si::meters*units::si::seconds)*units::si::meters/units::si::seconds; }
template<class T = double> t_energy_per_temperature<units::si::system, T> get_kB()
{ return T(co::k_B*units::si::kelvin/units::si::joules)/units::si::kelvin*units::si::joules; }
template<class T = double> t_energy_per_field<units::si::system, T> get_mu_B()
{ return T(co::mu_B/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class T = double> t_energy_per_field<units::si::system, T> get_mu_n()
{ return T(co::mu_n/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class T = double> t_energy_per_field<units::si::system, T> get_mu_N()
{ return T(co::mu_N/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class T = double> t_energy_per_field<units::si::system, T> get_mu_e()
{ return T(co::mu_e/units::si::joules*units::si::tesla)*units::si::joules/units::si::tesla; }
template<class T = double> t_inductance_per_length<units::si::system, T> get_mu_0()
{ return T(co::mu_0/units::si::henry*units::si::meter)*units::si::henry/units::si::meter; }
template<class T = double> T get_g_n() { return T(co::g_n.value()); }
template<class T = double> T get_g_e() { return T(co::g_e.value()); }
template<class T = double> t_inv_flux_time<units::si::system, T> get_gamma_n()
{ return T(co::gamma_n*units::si::tesla*units::si::seconds)/units::si::tesla/units::si::seconds; }
template<class T = double> t_inv_flux_time<units::si::system, T> get_gamma_e()
{ return T(co::gamma_e*units::si::tesla*units::si::seconds)/units::si::tesla/units::si::seconds; }
template<class T = double> t_length<units::si::system, T> get_r_e()
{ return T(co::r_e/units::si::meters)*units::si::meters; }


// template constants
template<class T = double> const t_length_si<T> t_meters = get_one_meter<T>();
template<class T = double> const t_flux_si<T> t_teslas = T(1)*units::si::teslas;
template<class T = double> const t_time_si<T> t_seconds = T(1)*units::si::seconds;
template<class T = double> const t_temperature_si<T> t_kelvins = T(1)*units::si::kelvins;
template<class T = double> const t_area_si<T> t_barns = T(1e-28)*units::si::meters*units::si::meters;
template<class T = double> const t_angle_si<T> t_radians = T(1)*units::si::radians;
template<class T = double> const t_angle_si<T> t_degrees = get_pi<T>()/T(180)*units::si::radians;

template<class T = double> const t_energy_si<T> t_meV = get_one_meV<T>();
template<class T = double> const t_length_si<T> t_angstrom = get_one_angstrom<T>();


// constants
static const length meters = 1.*units::si::meters;
static const flux teslas = 1.*units::si::teslas;
static const time seconds = 1.*units::si::seconds;
static const angle radians = 1.*units::si::radians;
static const temp kelvins = 1.*units::si::kelvins;
static const mass amu = co::m_u;
static const area barns = 1e-28 * meters*meters;

static const energy one_meV = 1e-3 * co::e * units::si::volts;
static const energy one_eV = co::e * units::si::volts;
static const length angstrom = 1e-10 * meters;
static const length cm = meters/100.;
static const time ps = 1e-12 * seconds;


// synonyms
static const temp kelvin = kelvins;
static const length meter = meters;
static const time second = seconds;
static const energy meV = one_meV;
static const flux tesla = teslas;
static const area barn = barns;


// helper functions
template<class t_quant>
t_quant my_units_sqrt(const decltype(t_quant() * t_quant())& val)
{
	using t_quant_sq = decltype(t_quant() * t_quant());
	using T = typename t_quant::value_type;

	t_quant one_quant = t_quant::from_value(T(1));
	t_quant_sq one_quant_sq = t_quant_sq::from_value(T(1));

	T valsq = T(val / one_quant_sq);
	return std::sqrt(valsq) * one_quant;
}


template<class t_quant>
decltype(t_quant()*t_quant()) my_units_pow2(const t_quant& val)
{
	return val*val;
}


template<class t_elem, template<class...> class t_vec=boost::numeric::ublas::vector>
t_elem my_units_norm2(const t_vec<t_elem>& vec)
{
	using t_elem_sq = decltype(t_elem()*t_elem());
	t_elem_sq tRet = t_elem_sq();

	for(std::size_t i=0; i<vec.size(); ++i)
		tRet += vec[i]*vec[i];

	return tl::my_units_sqrt<t_elem>(tRet);
}

}
#endif
