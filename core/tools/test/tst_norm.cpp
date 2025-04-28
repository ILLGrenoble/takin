/**
 * test normalisation of random normal distribution
 * @author Tobias Weber <tobias.weber@tum.de>
 * @data may-2019
 * @license GPLv2
 *
 * g++ -std=c++14 -I../../ -o tst_norm tst_norm.cpp ../../tlibs/math/rand.cpp ../../tlibs/log/log.cpp
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/histogram.hpp>
#include "tlibs/math/linalg.h"
#include "tlibs/math/rand.h"


/**
 * test to ensure that the histogram doesn't perform an internal normalisation
 */
template<class t_real>
void tst_histo()
{
	auto histo_axis = std::vector<boost::histogram::axis::regular<t_real>>
	{
		boost::histogram::axis::regular<t_real>{
			5 /* bins*/, 0 /* min */, 5 /* max */ }
	};

	auto histo = boost::histogram::make_histogram(histo_axis);

	histo(0);
	histo(0);
	histo(1);
	histo(1);
	histo(1);
	histo(3);

	std::ostream& ostr = std::cout;
	ostr.precision(5);

	ostr << std::left
		<< std::setw(12) << "lower" << " "
		<< std::setw(12) << "upper" << " "
		<< std::setw(12) << "count" << " "
		<< std::endl;

	for(const auto& val : boost::histogram::indexed(histo))
	{
		t_real x = val.bin().lower() + 0.5*(val.bin().upper() - val.bin().lower());

		ostr << std::left
			<< std::setw(12) << val.bin().lower() << " "
			<< std::setw(12) << val.bin().upper() << " "
			<< std::setw(12) << *val << "\n";
	}

	ostr << std::endl;
}


/**
 * tests normalisation of random gauss distribution
 */
template<class t_real>
void tst_norm(t_real sigma, t_real mu, unsigned int iters, unsigned int bins, std::ostream& ostr = std::cout)
{
	t_real range = 6. * sigma;
	//range = 1.;
	t_real bin_density = static_cast<t_real>(bins) / range;

	auto histo_axis = std::vector<boost::histogram::axis::regular<t_real>>
	{
		boost::histogram::axis::regular<t_real>
		{
			bins,
			mu - range / 2. /* min */,
			mu + range / 2. /* max */
		}
	};

	auto histo = boost::histogram::make_histogram(histo_axis);

	for(std::size_t i = 0; i < iters; ++i)
	{
		t_real val = tl::rand_norm<t_real>(mu, sigma);
		histo(val);
	}

	ostr.precision(5);
	ostr << std::left
		<< std::setw(12) << "# x" << " "
		<< std::setw(12) << "y_mc" << " "
		<< std::setw(12) << "y_amp" << " "
		<< std::setw(12) << "y_area" << " "
		<< std::setw(12) << "y_amp/y_mc" << " "
		<< std::setw(12) << "y_area/y_mc"
		<< std::endl;

	for(const auto& val : boost::histogram::indexed(histo))
	{
		t_real x = val.bin().lower() + 0.5*(val.bin().upper() - val.bin().lower());
		t_real yMC = *val / static_cast<t_real>(iters);
		//yMC *= bin_density * 8. * sigma / M_PI;  // amp normalisation
		//yMC *= bin_density * 8. / M_PI / std::sqrt(2. * M_PI);  // area normalisation
		yMC *= bin_density * std::sqrt(32. / M_PI/M_PI/M_PI);  // area normalisation
		t_real yModel_amp = tl::gauss_model_amp<t_real>(x, mu, sigma, 1., 0.);
		t_real yModel_area = tl::gauss_model<t_real>(x, mu, sigma, 1., 0.);

		ostr << std::left
			<< std::setw(12) << x << " "
			<< std::setw(12) << yMC << " "
			<< std::setw(12) << yModel_amp << " "
			<< std::setw(12) << yModel_area << " "
			<< std::setw(12) << yModel_amp/yMC << " "
			<< std::setw(12) << yModel_area/yMC << "\n";
	}

	ostr << std::endl;
}


int main(int argc, char** argv)
{
	using t_real = double;

	if(argc >= 6)
	{
		// run gaussian test with given parameters
		t_real sig = 1.;
		t_real mu = 0.;
		unsigned int iters = 100000;
		unsigned int bins = 64;
		std::string filename = "gauss.dat";

		std::istringstream{argv[1]} >> sig;
		std::istringstream{argv[2]} >> mu;
		std::istringstream{argv[3]} >> iters;
		std::istringstream{argv[4]} >> bins;
		std::istringstream{argv[5]} >> filename;

		std::ofstream ofstr{filename};
		tst_norm<t_real>(sig, mu, iters, bins, ofstr);

		std::cerr << "Plot with \"plot \"" << filename << "\" u 1:2 w l lw 2, "
			<< "\"" << filename << "\" u 1:4 w l lw 2\"."
			<< std::endl;
	}
	else
	{
		std::cerr << "Running standard tests. To use specific parameters, use the format:\n"
			<< "\t" << argv[0] << " <sigma> <mu> <iters> <bins> <filename>\n"
			<< std::endl;

		// run some tests
		tst_histo<t_real>();
		tst_norm<t_real>(2.5, 0., 250000, 50);
		tst_norm<t_real>(3.5, 0., 250000, 100);
	}

	return 0;
}
