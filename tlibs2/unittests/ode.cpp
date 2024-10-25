/**
 * test ode calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date 28-oct-20
 * @license GPLv3, see 'LICENSE' file
 * @desc forked from https://github.com/t-weber/misc/blob/master/boost/ode3.cpp
 *
 * g++ -std=c++20 -DUSE_LAPACK -I.. -I/usr/include/lapacke -Iext/lapacke/include -Lext/lapacke/lib -o ode ode.cpp -llapacke
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#define BOOST_TEST_MODULE Ode Test
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"
using namespace tl2_ops;

#define DEBUG_OUTPUT


// ----------------------------------------------------------------------------
// numeric integration
#include <boost/numeric/odeint.hpp>
namespace odeint = boost::numeric::odeint;

// mark custom vector as resizeable
template<> struct odeint::is_resizeable<tl2::vec<float, std::vector>> : std::true_type {};
template<> struct odeint::is_resizeable<tl2::vec<double, std::vector>> : std::true_type {};

template<class t_mat, class t_vec, class t_val=typename t_vec::value_type>
t_vec odesys(const t_mat& C, const t_vec& y0, t_val x_start, t_val x_end, t_val x_step = 0.01)
{
	t_vec y = y0;
	odeint::integrate_adaptive(odeint::runge_kutta4<t_vec>{},
		[&C](const t_vec& y, t_vec& y_diff, t_val x) -> void
		{
			y_diff = C*y;
		}, y, x_start, x_end, x_step);

	return y;
}
// ----------------------------------------------------------------------------


using t_types = std::tuple<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_ode, t_real, t_types)
{
	const t_real eps = 1e-3;

	using t_cplx = std::complex<t_real>;
	using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
	using t_mat_cplx = tl2::mat<t_cplx, std::vector>;
	using t_vec = tl2::vec<t_real, std::vector>;
	using t_mat = tl2::mat<t_real, std::vector>;


	t_mat_cplx coeff = tl2::create<t_mat_cplx>({
		0., 1.,
		1., 1. });

	t_vec_cplx f0 = tl2::create<t_vec_cplx>({1., 1.});
	t_cplx x0 = 0.;
	t_cplx n0 = 0.;

	const auto [coeff_re, coeff_im] = tl2::split_cplx<t_mat_cplx, t_mat>(coeff);
	const auto [f0_re, f0_im] = tl2::split_cplx<t_vec_cplx, t_vec>(f0);


#ifdef DEBUG_OUTPUT
	std::cout << "coeff = " << coeff << std::endl;
	std::cout << "x0 = " << x0 << ", f0 = " << f0 << std::endl;
	std::cout << std::endl;
#endif

	for(t_cplx x = 0; x.real() < 10; x += 1)
	{
		auto [ok, f] = tl2_la::odesys_const<t_mat_cplx, t_vec_cplx, t_cplx>(coeff, x, x0, f0);
		BOOST_TEST(ok);

		// compare with numerical result
		t_vec num_val = odesys<t_mat, t_vec>(coeff_re, f0_re, t_real{0.}, x.real(), t_real{0.01});
		BOOST_TEST(f[0].real() == num_val[0], testtools::tolerance(eps));
		BOOST_TEST(f[1].real() == num_val[1], testtools::tolerance(eps));
		BOOST_TEST(f[0].imag() == t_real{0}, testtools::tolerance(eps));
		BOOST_TEST(f[1].imag() == t_real{0}, testtools::tolerance(eps));

#ifdef DEBUG_OUTPUT
		std::cout << "x = " << x << std::endl;
		std::cout << "ok = " << std::boolalpha << ok << std::endl;

		for(std::size_t i = 0; i < f.size(); ++i)
			std::cout << "f_" << i << " = " << f[i] << std::endl;
		std::cout << std::endl;
#endif
	}


	std::size_t idx=0;
	const t_cplx fibo[] = { 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89 };
	for(t_cplx n = 0; n.real() < 10; n += 1)
	{
		auto [ok, f] = tl2_la::diffsys_const<t_mat_cplx, t_vec_cplx, t_cplx>(coeff, n, n0, f0);
		BOOST_TEST(ok);
		BOOST_TEST(tl2::equals(f[0], fibo[idx], eps));
		BOOST_TEST(tl2::equals(f[1], fibo[idx+1], eps));

#ifdef DEBUG_OUTPUT
		std::cout << "n = " << n << std::endl;
		std::cout << "ok = " << std::boolalpha << ok << std::endl;

		for(std::size_t i=0; i<f.size(); ++i)
			std::cout << "f_" << i << " = " << f[i] << std::endl;
		std::cout << std::endl;
#endif
		++idx;
	}
}
