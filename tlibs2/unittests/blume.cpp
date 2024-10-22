/**
 * phys lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 18-oct-2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// g++ -std=c++20 -I.. -o blume blume.cpp

#define BOOST_TEST_MODULE Blume0
#include <boost/test/included/unit_test.hpp>
namespace test = boost::unit_test;
namespace testtools = boost::test_tools;

#include <iostream>
#include <vector>

#include "libs/maths.h"
#include "libs/phys.h"


using t_types = std::tuple<double, float>;
BOOST_AUTO_TEST_CASE_TEMPLATE(test_blume0, t_real, t_types)
{
	using namespace tl2_ops;

	using t_cplx = std::complex<t_real>;
	using t_vec = tl2::vec<t_cplx, std::vector>;
	using t_mat = tl2::mat<t_cplx, std::vector>;

	constexpr const t_cplx imag(0, 1);

	t_vec Pi = tl2::create<t_vec>({
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1) });
	Pi /= std::sqrt((Pi[0]*Pi[0] + Pi[1]*Pi[1] + Pi[2]*Pi[2]).real());

	t_vec Mperp = tl2::create<t_vec>({
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1) }) +
		imag * tl2::create<t_vec>({
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1),
		tl2::get_rand<t_real>(-1, 1) });
	t_cplx N = tl2::get_rand<t_real>(0, 1) + imag*tl2::get_rand<t_real>(0, 1);

	//Mperp = ortho_project(Mperp, Pi);

	auto [ I, Pf ] = tl2::blume_maleev<t_vec, t_cplx>(Pi, Mperp, N);	
	auto [ I2, Pf2 ] = tl2::blume_maleev_indir<t_mat, t_vec, t_cplx>(Pi, Mperp, N);	
	auto [ I3, Prot3, Pcreate3, Pf3 ] = tl2::blume_maleev_tensor<t_mat, t_vec, t_cplx>(Pi, Mperp, N);	

	t_real lenPf = std::sqrt((Pf[0]*Pf[0] + Pf[1]*Pf[1] + Pf[2]*Pf[2]).real());
	t_real lenPf2 = std::sqrt((Pf2[0]*Pf2[0] + Pf2[1]*Pf2[1] + Pf2[2]*Pf2[2]).real());
	t_real lenPf3 = std::sqrt((Pf3[0]*Pf3[0] + Pf3[1]*Pf3[1] + Pf3[2]*Pf3[2]).real());
	//Pf3 /= t_cplx(lenPf3);

	t_cplx eps = std::pow(std::numeric_limits<t_real>::epsilon(), 1./2.);
	std::cout << "eps = " << eps.real() << std::endl;

	std::cout << "Mperp = " << Mperp << std::endl;
	std::cout << "N = " << N << std::endl;
	std::cout << "Pi = " << Pi << std::endl;
	std::cout << "Pf1 = " << Pf << std::endl;
	std::cout << "Pf2 = " << Pf2 << std::endl;
	std::cout << "Pf3 = " << Pf3 << std::endl;
	std::cout << "I1 = " << I << std::endl;
	std::cout << "I2 = " << I2 << std::endl;
	std::cout << "I3 = " << I3 << std::endl;
	std::cout << "|Pf1| = " << lenPf << std::endl;
	std::cout << "|Pf2| = " << lenPf2 << std::endl;
	std::cout << "|Pf3| = " << lenPf3 << std::endl;

	BOOST_TEST(tl2::equals(I, I2, eps));
	BOOST_TEST(tl2::equals(Pf, Pf2, eps));
	BOOST_TEST(tl2::equals(lenPf, lenPf2, eps.real()));
	BOOST_TEST(tl2::equals(I, I3, eps));
	BOOST_TEST(tl2::equals(Pf, Pf3, eps));
	BOOST_TEST(tl2::equals(lenPf, lenPf3, eps.real()));
}
