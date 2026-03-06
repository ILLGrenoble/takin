/**
 * tlibs test file
 * @author Tobias Weber <tweber@ill.fr>
 * @date 6-mar-2026
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// clang++ -std=c++17 -DNO_IOSTR -I.. -o loadfrmold ../../log/log.cpp ../../file/loadinstr.cpp loadfrmold.cpp -lboost_iostreams

#include <iostream>
#include "../file/loadinstr.h"


int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cerr << "Please give a file name." << std::endl;
		return -1;
	}

	tl::FileFrmOld<double> dat;

	if(!dat.Load(argv[1]))
	{
		std::cerr << "Could not load \"" << argv[1] << "\"." << std::endl;
		return -2;
	}

	const auto& col = dat.GetCol("E");
	for(auto val : col)
		std::cout << val << " ";
	std::cout << std::endl;

	return 0;
}
