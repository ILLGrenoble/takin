/**
 * tlibs2 maths library -- qhull interface
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_MATHS_QHULL_H__
#define __TLIBS2_MATHS_QHULL_H__

#include <cmath>
#include <tuple>
#include <vector>

#include "decls.h"
#include "operators.h"


#ifdef __TLIBS2_USE_QHULL__
	#include <Qhull.h>
	#include <QhullFacetList.h>
	#include <QhullVertexSet.h>
#else
	//#pragma message("tlibs2: Disabling QHull library (not found).")
#endif



#ifdef __TLIBS2_USE_QHULL__

namespace tl2_qh {

/**
 * calculates the convex hull
 * @see https://github.com/t-weber/misc/blob/master/geo/qhulltst.cpp
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_cont<t_vec>>, t_cont<t_vec>, t_cont<typename t_vec::value_type>>
get_convexhull(const t_cont<t_vec>& vecVerts)
requires tl2::is_vec<t_vec>
{
	using namespace tl2_ops;

	using t_real = typename t_vec::value_type;
	using t_real_qh = double;
	using t_facetlist_iter = typename orgQhull::QhullLinkedList<orgQhull::QhullFacet>::iterator;
	using t_vertexset_iter = typename orgQhull::QhullSet<orgQhull::QhullVertex>::iterator;


	// copy vertices
	int dim = vecVerts[0].size();
	std::size_t len = vecVerts.size()*dim;
	std::unique_ptr<t_real_qh[]> mem{new t_real_qh[len]};

	std::size_t i=0;
	for(const t_vec& vert : vecVerts)
	{
		for(int d=0; d<dim; ++d)
		{
			mem[i] = t_real_qh(vert[d]);
			++i;
		}
	}


	orgQhull::Qhull qhull{"tlibs2", dim, int(vecVerts.size()), mem.get(), "Qt"};
	orgQhull::QhullFacetList facets = qhull.facetList();

	t_cont<t_cont<t_vec>> vecPolys;
	t_cont<t_vec> vecNormals;
	t_cont<t_real> vecDists;
	vecPolys.reserve(facets.size());
	vecNormals.reserve(facets.size());
	vecDists.reserve(facets.size());
	//const t_vec vecCentre = tl2::mean(vecVerts);

	for(t_facetlist_iter iter=facets.begin(); iter!=facets.end(); ++iter)
	{
		// triangulable?
		if(iter->isUpperDelaunay())
			continue;

		orgQhull::QhullVertexSet vertices = iter->vertices();

		t_cont<t_vec> vecPoly;
		vecPoly.reserve(vertices.size());

		for(t_vertexset_iter iterVertex=vertices.begin(); iterVertex!=vertices.end(); ++iterVertex)
		{
			orgQhull::QhullPoint point = (*iterVertex).point();
			t_vec vecPoint(dim);
			for(int i=0; i<dim; ++i)
				vecPoint[i] = t_real(point[i]);

			vecPoly.emplace_back(std::move(vecPoint));
		}


		t_vec vecNormal; //= tl2::sort_poly_verts<t_vec, t_cont>(vecPoly, vecCentre, true);
		vecNormal.resize(dim);

		orgQhull::QhullHyperplane plane = iter->hyperplane();
		const t_real_qh* planenorm = plane.coordinates();
		const t_real_qh planedist = plane.offset();
		for(int i = 0; i < dim; ++i)
			vecNormal[i] = t_real(planenorm[i]);

		vecPolys.emplace_back(std::move(vecPoly));
		vecNormals.emplace_back(std::move(vecNormal));
		vecDists.emplace_back(planedist);
	}

	// too few polygons => remove polyhedron
	if(vecPolys.size() < 3)
		vecPolys = decltype(vecPolys){};

	return std::make_tuple(vecPolys, vecNormals, vecDists);
}


}
#endif  // __TLIBS2_USE_QHULL__


#endif
