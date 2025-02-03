/**
 * polarisation test
 * @author Tobias Weber <tweber@ill.fr>
 * @date aug-18
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 *
 * g++ -std=c++20 -I../.. -o pol pol.cpp
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "magtools" project
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
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

#include <iostream>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/phys.h"
#include "libs/defs.h"


using t_vec = tl2::vec<t_cplx, std::vector>;
using t_mat = tl2::mat<t_cplx, std::vector>;
using t_matvec = std::vector<t_mat>;


std::istringstream get_istr(std::istream& istr)
{
	std::string _str;
	std::getline(istr, _str);
	return std::istringstream(_str);
}


int main()
{
	using namespace tl2_ops;

	constexpr bool check = true;
	t_real eps = 1e-5;

	t_cplx N(0,0);
	t_vec Mperp = tl2::create<t_vec>({ 0, 0, t_cplx(1,2) });
	t_vec P = tl2::create<t_vec>({ 0, 1, 1 });


	std::cout << "N = "; get_istr(std::cin) >> N;
	std::cout << "Mperp = "; std::cin >> Mperp; Mperp.resize(3);
	std::cout << "P = "; std::cin >> P; P.resize(3);

	std::cout << "\n";
	std::cout << "N = " << N << "\n";
	std::cout << "Mperp = " << Mperp << "\n";
	std::cout << "P = " << P << "\n";
	std::cout << std::endl;

	std::cout << "density matrix = " << tl2::pol_density_mat<t_vec, t_mat>(P) << "\n";
	std::cout << std::endl;

	auto [ I, P_f ] = tl2::blume_maleev_indir<t_mat, t_vec, t_cplx>(P, Mperp, N);
	std::cout << "I = " << I << "\nP_f = " << P_f << std::endl;

	// double-check results
	if constexpr(check)
	{
		auto [ I2, P_f2 ] = tl2::blume_maleev<t_vec, t_cplx>(P, Mperp, N);
		if(!tl2::equals<t_cplx>(I, I2, eps) || !equals<t_vec, t_cplx>(P_f, P_f2, eps))
		{
			std::cerr << "Mismatch between blume_maleev() and blume_maleev_indir()!" << std::endl;
		}
	}

	return 0;
}
