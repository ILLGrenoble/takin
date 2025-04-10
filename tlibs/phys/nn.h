/**
 * nearest neighbours
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date may-2016
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

#ifndef __TLIBS_NN_H__
#define __TLIBS_NN_H__

#include <initializer_list>
#include <vector>
#include <tuple>
#include <cmath>
#include <complex>

#include "../math/linalg.h"
#include "../helper/misc.h"


namespace tl {


template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector>
t_cont<t_vec>
get_atoms_by_idx(const t_cont<t_vec>& vecAtoms, const t_cont<std::size_t>& vecIdx)
{
	t_cont<t_vec> vecRet;
	vecRet.reserve(vecIdx.size());

	for(std::size_t iIdx : vecIdx)
		vecRet.push_back(vecAtoms[iIdx]);

	return vecRet;
}



/**
 * get nearest neighbours
 * @return vector of neighbours
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
t_cont<t_cont<std::size_t>>
get_neighbours(const t_cont<t_vec>& vecAtoms, const t_vec& vecCentre,
	t_real epsShell = t_real(1e-3))
{
	t_cont<t_cont<std::size_t>> vecN;

	// generate lengths
	t_cont<t_real> vecLens;
	vecLens.reserve(vecAtoms.size());
	for(const t_vec& vec : vecAtoms)
		vecLens.push_back(veclen(vec - vecCentre));


	// sort by lengths
	auto fktSort = [](t_real len0, t_real len1) -> bool
	{
		return len0 < len1;
	};

	t_cont<std::size_t> vecIdx = sorted_idx<std::vector, t_real, decltype(fktSort)>
		(vecLens, fktSort);

	if(vecIdx.size() == 0) return vecN;

	// sort by nearest neighbour, next-nearest neighbour, etc.
	t_real distLast = vecLens[vecIdx[0]];
	t_cont<std::size_t> vecNearest;

	for(typename decltype(vecIdx)::const_iterator iter=vecIdx.begin(); iter!=vecIdx.end(); ++iter)
	{
		std::size_t iIdx = *iter;
		//std::cout << vecLens[iIdx] << ": " << vecAtoms[iIdx] << std::endl;

		t_real distCur = vecLens[iIdx];
		if(float_equal(distCur, distLast, epsShell))
		{
			vecNearest.push_back(iIdx);
		}
		else
		{
			vecN.push_back(std::move(vecNearest));

			vecNearest.clear();
			vecNearest.push_back(iIdx);
			distLast = distCur;
		}

		if(iter+1 == vecIdx.end() && vecNearest.size())
			vecN.push_back(vecNearest);
	}
	return vecN;
}



// ----------------------------------------------------------------------------
// cubic systems

// atom positions
enum class UCType { SIMPLE, FCC, BCC, };

/**
 * Next neighbours
 * iDist == 0: nearest neighbours
 * iDist == 1: next-nearest neighbours
 */
template<typename T=double>
std::vector<ublas::vector<T>> get_neighbour_atoms(UCType crys, int iDist=0, T a=1.)
{
	std::vector<ublas::vector<T>> vecAtoms;

	// generated with Takin: one atom at (000) and sc space group e.g. P-43m,
	//                       "Real Space" -> "Information..." -> "Unit Cell"
	if(crys == UCType::SIMPLE)
	{
		if(iDist == 0)
			vecAtoms = {
				make_vec({1., 0., 0.}),
				make_vec({0., 1., 0.}),
				make_vec({0., 0., 1.}),
				make_vec({-1., 0., 0.}),
				make_vec({0., -1., 0.}),
				make_vec({0., 0., -1.}) };
		else if(iDist == 1)
			vecAtoms = {
				make_vec({1., 1., 0.}),
				make_vec({1., -1., 0.}),
				make_vec({-1., 1., 0.}),
				make_vec({-1., -1., 0.}),
				make_vec({1., 0., 1.}),
				make_vec({1., 0., -1.}),
				make_vec({-1., 0., 1.}),
				make_vec({-1., 0., -1.}),
				make_vec({0., 1., 1.}),
				make_vec({0., 1., -1.}),
				make_vec({0., -1., 1.}),
				make_vec({0., -1., -1.}) };
	}
	// generated with Takin: one atom at (000) and fcc space group e.g. F-43m
	//                       "Real Space" -> "Information..." -> "Unit Cell"
	else if(crys == UCType::FCC)
	{
		if(iDist == 0)
			vecAtoms = {
				make_vec({0.5, 0.5, 0.}),
				make_vec({0.5, -0.5, 0.}),
				make_vec({-0.5, 0.5, 0.}),
				make_vec({-0.5, -0.5, 0.}),
				make_vec({0.5, 0., 0.5}),
				make_vec({0.5, 0., -0.5}),
				make_vec({-0.5, 0., 0.5}),
				make_vec({-0.5, 0., -0.5}),
				make_vec({0., 0.5, 0.5}),
				make_vec({0., 0.5, -0.5}),
				make_vec({0., -0.5, 0.5}),
				make_vec({0., -0.5, -0.5}) };
		else if(iDist == 1)
			vecAtoms = {
				make_vec({1., 0., 0.}),
				make_vec({0., 1., 0.}),
				make_vec({0., 0., 1.}),
				make_vec({-1., 0., 0.}),
				make_vec({0., -1., 0.}),
				make_vec({0., 0., -1.}) };
	}
	// generated with Takin: one atom at (000) and bcc space group e.g. I-43m
	//                       "Real Space" -> "Information..." -> "Unit Cell"
	else if(crys == UCType::BCC)
	{
		if(iDist == 0)
			vecAtoms = {
				make_vec({0.5, 0.5, 0.5}),
				make_vec({0.5, 0.5, -0.5}),
				make_vec({0.5, -0.5, 0.5}),
				make_vec({0.5, -0.5, -0.5}),
				make_vec({-0.5, 0.5, 0.5}),
				make_vec({-0.5, 0.5, -0.5}),
				make_vec({-0.5, -0.5, 0.5}),
				make_vec({-0.5, -0.5, -0.5}) };
		else if(iDist == 1)
			vecAtoms = {
				make_vec({1., 0., 0.}),
				make_vec({0., 1., 0.}),
				make_vec({0., 0., 1.}),
				make_vec({-1., 0., 0.}),
				make_vec({0., -1., 0.}),
				make_vec({0., 0., -1.}) };
	}

	if(!float_equal<T>(a, T(1)))
	{
		for(ublas::vector<T>& vec : vecAtoms)
			vec *= a;
	}

	return vecAtoms;
}
// ----------------------------------------------------------------------------

}
#endif
