 /**
 * tlibs test file
 * @author Tobias Weber <tweber@ill.fr>
 * @date 18-june-2026
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

// g++ -I.. -o latt2 latt2.cpp ../../log/log.cpp -std=c++17


#include "../phys/lattice.h"
#include <iostream>

using T = double;
using t_mat = tl::ublas::matrix<T>;
using t_vec = tl::ublas::vector<T>;


int main()
{
	t_vec vec1 = tl::make_vec<t_vec>({ 1., -1., 0. });
	t_vec vec2 = tl::make_vec<t_vec>({ 0.,  0., 1. });

	T alpha = 90., beta = 90., gamma = 120.;
	tl::Lattice<T> latt(4., 5., 6., alpha/180.*M_PI, beta/180.*M_PI, gamma/180.*M_PI);

	t_mat UB = tl::get_UB(latt, vec1, vec2);
	std::cout << "UB = " << UB << std::endl;

	std::array<t_vec, 3> plane = tl::get_scattering_plane(UB, latt);
	std::cout << "vec1 = " << plane[0] << std::endl;
	std::cout << "vec2 = " << plane[1] << std::endl;
	std::cout << "vec3 = " << plane[2] << std::endl;

	return 0;
}
