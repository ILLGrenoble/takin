/**
 * math lib test
 * @author Tobias Weber <tweber@ill.fr>
 * @date 16-feb-2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include <iostream>

#include "libs/maths.h"


int main()
{
	using t_int = int;

	std::cout
		<< std::setw(15) << std::left << "n"
		<< std::setw(15) << std::left << "k"
		<< std::setw(15) << std::left << "var_rep"
		<< std::setw(15) << std::left << "var_norep"
		<< std::setw(15) << std::left << "comb_rep"
		<< std::setw(15) << std::left << "comb_norep"
		<< std::endl;

	for(t_int n = 1; n < 10; ++n)
	for(t_int k = 0; k < 10; ++k)
	{
		t_int var_rep = tl2::combinatorics<t_int>(n, k, true, true);
		t_int var_norep = tl2::combinatorics<t_int>(n, k, true, false);
		t_int comb_rep = tl2::combinatorics<t_int>(n, k, false, true);
		t_int comb_norep = tl2::combinatorics<t_int>(n, k, false, false);

		std::cout
			<< std::setw(15) << std::left << n
			<< std::setw(15) << std::left << k
			<< std::setw(15) << std::left << var_rep
			<< std::setw(15) << std::left << var_norep
			<< std::setw(15) << std::left << comb_rep
			<< std::setw(15) << std::left << comb_norep
			<< std::endl;
	}

	return 0;
}
