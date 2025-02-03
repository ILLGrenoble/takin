/**
 * tlibs2
 * type traits and concepts library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2014-2021
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Further functions forked on 19-Apr-2021 from my privately developed "geo" project (https://github.com/t-weber/geo).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * magtools, geo, and "misc" projects
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

#ifndef __TLIBS2_TRAITS_H__
#define __TLIBS2_TRAITS_H__

#include <type_traits>
#include <vector>
#include <array>
#include <list>
#include <initializer_list>
#include <utility>
#include <concepts>
#include <complex>


namespace tl2 {
// ----------------------------------------------------------------------------
// concepts
// ----------------------------------------------------------------------------

/**
 * requirements for having a value_type
 */
template<class T>
concept has_value_type = requires { typename T::value_type; };


/**
 * requirements for a scalar type
 */
template<class T>
concept is_scalar =
	std::is_floating_point_v<T> ||
	std::is_integral_v<T> /*|| std::is_arithmetic_v<T>*/;


/**
 * requirements for a basic vector container like std::vector
 */
template<class T>
concept is_basic_vec = requires(const T& a)
{
	a.size();		// must have a size() member function
	a.operator[](1);	// must have an operator[]
} && has_value_type<T>;


/**
 * requirements of a vector type with a dynamic size
 */
template<class T>
concept is_dyn_vec = requires(const T& a)
{
	T(3);			// constructor
};


/**
 * requirements for a vector container
 */
template<class T>
concept is_vec = requires(const T& a)
{
	a + a;			// operator+
	a - a;			// operator-
	a[0] * a;		// operator*
	a * a[0];
	a / a[0];		// operator/
} && is_basic_vec<T>;


/**
 * requirements for a quaternion container
 * (modelled after Boost's quaternion interface)
 */
template<class T>
concept is_quat = requires(const T& a)
{
	a + a;			// operator+
	a - a;			// operator-
	a * a;			// operator*
	a / a;			// operator/

	a + a.R_component_1();
	a.R_component_1() + a;

	a.R_component_1() * a;	// operator*
	a * a.R_component_1();
	a / a.R_component_1();	// operator/

	a.R_component_1();
	a.R_component_2();
	a.R_component_3();
	a.R_component_4();
};


/**
 * requirements for a basic matrix container
 */
template<class T>
concept is_basic_mat = requires(const T& a)
{
	a.size1();		// must have a size1() member function
	a.size2();		// must have a size2() member function
	a.operator()(1,1);	// must have an operator()
} && has_value_type<T>;


/**
 * requirements of a matrix type with a dynamic size
 */
template<class T>
concept is_dyn_mat = requires(const T& a)
{
	T(3,3);			// constructor
};


/**
 * requirements for a matrix container
 */
template<class T>
concept is_mat = requires(const T& a)
{
	a + a;			// operator+
	a - a;			// operator-
	a(0,0) * a;		// operator*
	a * a(0,0);
	a / a(0,0);		// operator/
} && is_basic_mat<T>;



/**
 * requirements for a complex number
 */
template<class T>
concept is_complex = requires(const T& a)
{
	T(0, 0);		// constructor

	std::conj(a);
	a.real();		// must have a real() member function
	a.imag();		// must have an imag() member function

	a + a;
	a - a;
	a * a;
	a / a;
} && has_value_type<T>;


/**
 * requirements for an iterable container
 */
template<class T>
concept is_iterable = requires(const T& a)
{
	a.begin();
	a.end();

	{ a.begin() == a.end() } -> std::same_as<bool>;
	{ a.begin() != a.end() } -> std::same_as<bool>;
};
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
/**
 * if it exists, get T's underlying value_type, otherwise use the type T directly
 */
template<class T, bool b = has_value_type<T>> struct _underlying_value_type {};
template<class T> struct _underlying_value_type<T, false> { using ty = T; };
template<class T> struct _underlying_value_type<T, true> { using ty = typename T::value_type; };
template<class T> using underlying_value_type = typename _underlying_value_type<T>::ty;
// ----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
template<class T>
struct remove_constref
{
	typedef typename std::remove_const<
		typename std::remove_reference<T>::type
			>::type type;
};

template<class T>
using remove_constref_t = typename remove_constref<T>::type;
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
/**
 * function call implementation
 */
template<class t_func,
	class t_arg = double, template<class ...> class t_cont = std::vector,
	std::size_t... idx>
t_arg _call_impl(t_func func, const t_cont<t_arg>& args,
	const std::integer_sequence<std::size_t, idx...>&)
{
	return func(args[idx]...);
}

/**
 * function call implementation (specialisation for std::array)
 */
template<class t_func, class t_arg, std::size_t... idx>
t_arg _call_impl(t_func func, const std::array<t_arg, sizeof...(idx)>& args,
	const std::integer_sequence<std::size_t, idx...>&)
{
	return func(args[idx]...);
}


/**
 * call a function with the args from an STL container
 */
template<std::size_t iNumArgs, class t_func,
	class t_arg = double, template<class ...> class t_cont = std::vector>
t_arg call(t_func func, const t_cont<t_arg>& args)
{
	using t_seq = std::make_integer_sequence<std::size_t, iNumArgs>;
	return _call_impl<t_func, t_arg, t_cont>(func, args, t_seq());
}

/**
 * call a function with the args from a std::array
 */
template<std::size_t iNumArgs, class t_func, class t_arg = double>
t_arg call(t_func func, const std::array<t_arg, iNumArgs>& args)
{
	using t_seq = std::make_integer_sequence<std::size_t, iNumArgs>;
	return _call_impl<t_func, t_arg>(func, args, t_seq());
}


// -----------------------------------------------------------------------------


template<typename t_arg, std::size_t ...idx>
using _t_fkt_vararg_impl = t_arg(*)(
	typename std::remove_reference<
		decltype(std::declval<t_arg*>()[idx])
	>::type...);

template<typename t_arg, std::size_t ...idx>
static _t_fkt_vararg_impl<t_arg, idx...>
_tstfkt_vararg(const std::integer_sequence<std::size_t, idx...>&)
{ return nullptr; /* not interested in return value, only its type */ }


/**
 * constructs a function type with 'iNumArgs' arguments: t_arg (*) (t_arg, t_arg, ...)
 */
template<typename t_arg, std::size_t iNumArgs>
using t_fkt_vararg = decltype(
	_tstfkt_vararg<t_arg>(
		std::make_integer_sequence<std::size_t, iNumArgs>()));
// -----------------------------------------------------------------------------


// boolean value
template<bool value, class=void> constexpr bool bool_value = value;

}
#endif
