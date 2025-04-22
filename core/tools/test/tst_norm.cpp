/**
 * test normalisation of random normal distribution
 * @author Tobias Weber <tobias.weber@tum.de>
 * @data may-2019
 * @license GPLv2
 *
 * g++ -I../../ -o tst_norm tst_norm.cpp ../../tlibs/math/rand.cpp ../../tlibs/log/log.cpp
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/histogram.hpp>
#include "tlibs/math/linalg.h"
#include "tlibs/math/rand.h"

using t_real = double;

int main()
{
	constexpr const std::size_t ITERS = 100000;
	constexpr const std::size_t BINS = 64;

	t_real sigma = 2.5;
	t_real mu = 5.;

	auto histo_axis = std::vector<boost::histogram::axis::regular<t_real>>
	{
		boost::histogram::axis::regular<t_real>{ BINS, mu-4.*sigma /*min*/, mu+4.*sigma /*max*/ }
	};
	auto histo = boost::histogram::make_histogram(histo_axis);

	for(std::size_t i = 0; i < ITERS; ++i)
	{
		t_real val = tl::rand_norm<t_real>(mu, sigma);
		histo(val);
	}

	std::ostream& ostr = std::cout;
	ostr.precision(5);

	ostr << std::left
		<< std::setw(12) << "# x" << " "
		<< std::setw(12) << "y_mc" << " "
		<< std::setw(12) << "y_model" << " "
		<< std::setw(12) << "y_model/y_mc"
		<< std::endl;

	for(const auto& val : boost::histogram::indexed(histo))
	{
		t_real x = val.bin().lower() + 0.5*(val.bin().upper() - val.bin().lower());
		t_real yMC = *val / t_real{ITERS} * t_real{BINS} / M_PI;
		t_real yModel = tl::gauss_model_amp<t_real>(x, mu, sigma, 1., 0.);

		ostr << std::left << std::setw(12) << x << " " <<
			std::left << std::setw(12) << yMC << " " <<
			std::left << std::setw(12) << yModel << " " <<
			std::left << std::setw(12) << yModel/yMC << "\n";
	}

	return 0;
}
