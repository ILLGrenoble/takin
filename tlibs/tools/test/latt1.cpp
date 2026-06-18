 /**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 9-may-17
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

// g++ -I.. -o latt1 latt1.cpp ../../log/log.cpp -std=c++17


#include "../phys/lattice.h"
#include <iostream>

using T = double;
using t_mat = tl::ublas::matrix<T>;
using t_vec = tl::ublas::vector<T>;


int main()
{
	T alpha = 90., beta = 60., gamma = 45.;
	tl::Lattice<T> latt(4., 5., 6., alpha/180.*M_PI, beta/180.*M_PI, gamma/180.*M_PI);
	tl::Lattice<T> recip = latt.GetRecip();
	//std::cout << recip.GetPos(1., 2., 3.) << std::endl;

	t_mat matBaseCont = latt.GetBaseMatrixCont() * 2.*M_PI;
	t_mat matBaseCov = latt.GetBaseMatrixCov();
	t_mat matRecipCov = recip.GetBaseMatrixCov();

	// test reciprocal lattice covariant vector and real lattice contravariant vector equivalence
	std::cout << "latt cont = " << matBaseCont << std::endl;
	//std::cout << "latt cov = " << matBaseCov << std::endl;
	std::cout << "reci cov = " << matRecipCov << std::endl;
	//std::cout << "reci cov aligned = " << recip.GetAligned().GetBaseMatrixCov() << std::endl;

	t_mat matTst = tl::ublas::prod(tl::ublas::trans(matBaseCont), matBaseCov);
	tl::set_eps_0(matTst);
	std::cout << "base_cont * base_cov = " << matTst << std::endl;
	std::cout << std::endl;

	t_vec vecCol0 = tl::get_column(matBaseCont, 0); vecCol0 /= tl::ublas::norm_2(vecCol0);
	t_vec vecCol1 = tl::get_column(matBaseCont, 1); vecCol1 /= tl::ublas::norm_2(vecCol1);
	t_vec vecCol2 = tl::get_column(matBaseCont, 2); vecCol2 /= tl::ublas::norm_2(vecCol2);

	t_vec vecRec0 = tl::get_column(matRecipCov, 0); vecRec0 /= tl::ublas::norm_2(vecRec0);
	t_vec vecRec1 = tl::get_column(matRecipCov, 1); vecRec1 /= tl::ublas::norm_2(vecRec1);
	t_vec vecRec2 = tl::get_column(matRecipCov, 2); vecRec2 /= tl::ublas::norm_2(vecRec2);

	std::cout << "angle 01 = " << std::acos(tl::ublas::inner_prod(vecCol0, vecCol1)) / M_PI*180. << std::endl;
	std::cout << "angle 02 = " << std::acos(tl::ublas::inner_prod(vecCol0, vecCol2)) / M_PI*180. << std::endl;
	std::cout << "angle 12 = " << std::acos(tl::ublas::inner_prod(vecCol1, vecCol2)) / M_PI*180. << std::endl;

	std::cout << "angle 01 = " << std::acos(tl::ublas::inner_prod(vecRec0, vecRec1)) / M_PI*180. << std::endl;
	std::cout << "angle 02 = " << std::acos(tl::ublas::inner_prod(vecRec0, vecRec2)) / M_PI*180. << std::endl;
	std::cout << "angle 12 = " << std::acos(tl::ublas::inner_prod(vecRec1, vecRec2)) / M_PI*180. << std::endl;

	return 0;
}
