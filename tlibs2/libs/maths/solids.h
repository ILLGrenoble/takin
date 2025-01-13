/**
 * tlibs2 maths library -- 3-dim solids
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2025
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

#ifndef __TLIBS2_MATHS_SOLIDS_H__
#define __TLIBS2_MATHS_SOLIDS_H__

#include <cmath>
#include <vector>
#include <limits>

#include "decls.h"
#include "constants.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// 3-dim solids
// ----------------------------------------------------------------------------

/**
 * create the faces of a sphere
 * input: [triangle vertices, normals, uvs] (like subdivide_triangles)
 * @returns [triangles, face normals, vertex uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>
spherify(const std::tuple<t_cont<t_vec>, t_cont<t_vec>, t_cont<t_vec>>& tup,
	typename t_vec::value_type rad = 1)
requires is_vec<t_vec>
{
	const t_cont<t_vec>& vertices = std::get<0>(tup);
	//const t_cont<t_vec>& normals = std::get<1>(tup);
	const t_cont<t_vec>& uvs = std::get<2>(tup);

	t_cont<t_vec> vertices_new;
	t_cont<t_vec> normals_new;
	vertices_new.reserve(vertices.size());
	normals_new.reserve(vertices.size());


	// vertices
	for(t_vec vec : vertices)
	{
		vec /= tl2::norm<t_vec>(vec);
		vec *= rad;
		vertices_new.emplace_back(std::move(vec));
	}


	// face normals
	auto itervert = vertices.begin();
	// iterate over triplets forming triangles
	while(itervert != vertices.end())
	{
		const t_vec& vec1 = *itervert;
		std::advance(itervert, 1);
		if(itervert == vertices.end())
			break;
		const t_vec& vec2 = *itervert;
		std::advance(itervert, 1);
		if(itervert == vertices.end())
			break;
		const t_vec& vec3 = *itervert;
		std::advance(itervert, 1);

		t_vec vecmid = mean<t_vec>({ vec1, vec2, vec3 });
		vecmid /= tl2::norm<t_vec>(vecmid);
		normals_new.emplace_back(std::move(vecmid));
	}

	return std::make_tuple(vertices_new, normals_new, uvs);
}


/**
 * create a plane
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_mat, class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_plane(const t_vec& norm, typename t_vec::value_type lx = 1, typename t_vec::value_type ly = 1)
requires is_vec<t_vec>
{
	t_vec norm_old = create<t_vec>({ 0, 0, -1 });
	t_vec rot_vec = create<t_vec>({ 1, 0, 0 });
	t_mat rot = rotation<t_mat, t_vec>(norm_old, norm, &rot_vec);

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ -lx, -ly, 0. }),	// vertex 0
		create<t_vec>({ -lx, +ly, 0. }),	// vertex 1
		create<t_vec>({ +lx, +ly, 0. }),	// vertex 2
		create<t_vec>({ +lx, -ly, 0. }),	// vertex 3
	};

	// rotate according to given normal
	for(t_vec& vec : vertices)
		vec = rot * vec;

	t_cont<t_cont<std::size_t>> faces = { { 0, 1, 2, 3 } };
	t_cont<t_vec> normals = { norm };

	t_cont<t_cont<t_vec>> uvs =
	{{
		create<t_vec>({0, 0}),
		create<t_vec>({0, 1}),
		create<t_vec>({1, 1}),
		create<t_vec>({1, 0}),
	}};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a patch, z = f(x, y)
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_func, class t_mat, class t_vec,
	template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_patch(const t_func& func,
	t_real width = 1., t_real height = 1.,
	std::size_t num_points_x = 16, std::size_t num_points_y = 16)
requires is_vec<t_vec>
{
	// create 2d grid in (x, y) for patch
	t_cont<t_vec> vertices;
	t_cont<bool> valid_vertices;
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	vertices.reserve(num_points_x * num_points_y);
	valid_vertices.reserve(num_points_x * num_points_y);
	faces.reserve((num_points_x - 1) * (num_points_y - 1));
	normals.reserve((num_points_x - 1) * (num_points_y - 1));
	uvs.reserve((num_points_x - 1) * (num_points_y - 1));

	for(std::size_t j = 0; j < num_points_y; ++j)
	{
		t_real y = -height*0.5 + height *
			static_cast<t_real>(j)/static_cast<t_real>(num_points_y - 1);

		t_real v0 = static_cast<t_real>(j - 1) / static_cast<t_real>(num_points_x - 1);
		t_real v1 = static_cast<t_real>(j) / static_cast<t_real>(num_points_y - 1);

		for(std::size_t i = 0; i < num_points_x; ++i)
		{
			// create vertices
			t_real x = -width*0.5 + width *
				static_cast<t_real>(i)/static_cast<t_real>(num_points_x - 1);

			auto [z, valid] = func(x, y, i, j);
			vertices.emplace_back(tl2::create<t_vec>({ x, y, z }));
			valid_vertices.emplace_back(valid);

			if(i == 0 || j == 0)
				continue;

			// create faces
			std::size_t idx_ij = j*num_points_x + i;
			std::size_t idx_im1j = j*num_points_x + i - 1;
			std::size_t idx_i1jm1 = (j - 1)*num_points_x + i;
			std::size_t idx_im1jm1 = (j - 1)*num_points_x + i - 1;

			if(!valid_vertices[idx_ij] || !valid_vertices[idx_im1j] ||
				!valid_vertices[idx_i1jm1] || !valid_vertices[idx_im1jm1])
				continue;

			faces.emplace_back(t_cont<std::size_t>{{
				idx_im1jm1, idx_i1jm1, idx_ij, idx_im1j
			}});

			// create normals
			t_vec n = cross<t_vec>({
				vertices[idx_i1jm1] - vertices[idx_im1jm1],
				vertices[idx_ij] - vertices[idx_im1jm1]
			});
			n /= norm<t_vec>(n);
			normals.emplace_back(n);

			// create uv coordinates
			t_real u0 = static_cast<t_real>(i - 1) / static_cast<t_real>(num_points_x - 1);
			t_real u1 = static_cast<t_real>(i) / static_cast<t_real>(num_points_x - 1);

			uvs.emplace_back(t_cont<t_vec>{{
				create<t_vec>({ u0, v0 }),  // face vertex 0
				create<t_vec>({ u1, v0 }),  // face vertex 1
				create<t_vec>({ u1, v1 }),  // face vertex 2
				create<t_vec>({ u0, v1 }),  // face vertex 3
			}});
		}
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a disk
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_disk(typename t_vec::value_type r = 1, std::size_t num_points = 32)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// vertices and uvs
	t_cont<t_vec> vertices;
	t_cont<t_vec> all_uvs;
	vertices.reserve(num_points);
	all_uvs.reserve(num_points);

	// vertices
	for(std::size_t pt = 0; pt < num_points; ++pt)
	{
		const t_real phi = t_real(pt)/t_real(num_points) * t_real(2)*pi<t_real>;
		const t_real c = std::cos(phi);
		const t_real s = std::sin(phi);

		// outer vertices
		t_vec vert = create<t_vec>({ r*c, r*s, 0 });
		vertices.emplace_back(std::move(vert));

		t_vec uv = create<t_vec>({
			(t_real(1)+c) / t_real(2),
			(t_real(1)+s) / t_real(2) });
		all_uvs.emplace_back(std::move(uv));
	}

	// faces, normals & uvs
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	t_cont<std::size_t> face(num_points);
	std::iota(face.begin(), face.end(), 0);
	faces.emplace_back(std::move(face));

	normals.emplace_back(create<t_vec>({0,0,1}));

	t_cont<t_vec> uv;
	uv.reserve(num_points);
	for(std::size_t vert_idx : faces[0])
		uv.push_back(all_uvs[vert_idx]);
	uvs.emplace_back(std::move(uv));

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a cone
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cone(typename t_vec::value_type r = 1, typename t_vec::value_type h = 1,
	bool bWithCap = true, std::size_t num_points = 32)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// vertices
	t_cont<t_vec> vertices;
	vertices.reserve(num_points + 1);

	// inner vertex
	vertices.emplace_back(create<t_vec>({ 0, 0, h }));

	for(std::size_t pt = 0; pt < num_points; ++pt)
	{
		const t_real phi = t_real(pt)/t_real(num_points) * t_real(2)*pi<t_real>;
		const t_real c = std::cos(phi);
		const t_real s = std::sin(phi);

		// outer vertices
		t_vec vert = create<t_vec>({ r*c, r*s, 0 });
		vertices.emplace_back(std::move(vert));
	}

	// faces, normals & uvs
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;	// TODO

	faces.reserve(num_points);
	normals.reserve(num_points);
	uvs.reserve(num_points);

	for(std::size_t face = 0; face < num_points; ++face)
	{
		std::size_t idx0 = face + 1;	// outer 1
		std::size_t idx1 = (face == num_points-1 ? 1 : face + 2);	// outer 2
		std::size_t idx2 = 0;	// inner

		faces.push_back({ idx0, idx1, idx2 });

		t_vec n = tl2::cross<t_vec>({
			vertices[idx2]-vertices[idx0],
			vertices[idx1]-vertices[idx0] });
		n /= tl2::norm<t_vec>(n);

		normals.emplace_back(std::move(n));
	}

	if(bWithCap)
	{
		const auto [disk_vertices, disk_faces, disk_normals, disk_uvs]
			= create_disk<t_vec, t_cont>(r, num_points);

		// vertex indices have to be adapted for merging
		const std::size_t vert_start_idx = vertices.size();
		vertices.insert(vertices.end(), disk_vertices.begin(), disk_vertices.end());

		auto disk_faces_bottom = disk_faces;
		for(auto& disk_face : disk_faces_bottom)
		{
			for(auto& disk_face_idx : disk_face)
				disk_face_idx += vert_start_idx;
			std::reverse(disk_face.begin(), disk_face.end());
		}
		faces.insert(faces.end(), disk_faces_bottom.begin(), disk_faces_bottom.end());

		for(const auto& normal : disk_normals)
			normals.push_back(-normal);

		uvs.insert(uvs.end(), disk_uvs.begin(), disk_uvs.end());
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create a cylinder
 * cyltype: 0 (no caps), 1 (with caps), 2 (arrow)
 * @returns [vertices, face vertex indices, face normals, face uvs]
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cylinder(typename t_vec::value_type r = 1, typename t_vec::value_type h = 1,
	int cyltype = 0, std::size_t num_points = 32,
	typename t_vec::value_type arrow_r = 1.5, typename t_vec::value_type arrow_h = 0.5)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	// vertices
	t_cont<t_vec> vertices;
	t_cont<t_real> vertices_u;

	vertices.reserve(num_points*2);
	vertices_u.reserve(num_points*2);

	for(std::size_t pt = 0; pt < num_points; ++pt)
	{
		const t_real u = t_real(pt)/t_real(num_points);
		const t_real phi = u * t_real(2)*pi<t_real>;
		const t_real c = std::cos(phi);
		const t_real s = std::sin(phi);

		t_vec top = tl2::create<t_vec>({ r*c, r*s, h*t_real(0.5) });
		t_vec bottom = tl2::create<t_vec>({ r*c, r*s, -h*t_real(0.5) });

		vertices.emplace_back(std::move(top));
		vertices.emplace_back(std::move(bottom));

		vertices_u.push_back(u);
	}

	// faces, normals & uvs
	t_cont<t_cont<std::size_t>> faces;
	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	faces.reserve(num_points);
	normals.reserve(num_points);
	uvs.reserve(num_points);

	for(std::size_t face = 0; face < num_points; ++face)
	{
		std::size_t idx0 = face*2 + 0;	// top 1
		std::size_t idx1 = face*2 + 1;	// bottom 1
		std::size_t idx2 = (face >= num_points-1 ? 1 : face*2 + 3);	// bottom 2
		std::size_t idx3 = (face >= num_points-1 ? 0 : face*2 + 2);	// top 2

		t_vec n = tl2::cross<t_vec>({
			vertices[idx1]-vertices[idx0],
			vertices[idx3]-vertices[idx0] });
		n /= tl2::norm<t_vec>(n);

		faces.push_back({ idx0, idx1, idx2, idx3 });
		normals.emplace_back(std::move(n));

		t_real u1 = vertices_u[idx0];
		t_real u2 = (face >= num_points-1 ? 1 : vertices_u[idx3]);
		uvs.push_back({
			tl2::create<t_vec>({u1, 1}), tl2::create<t_vec>({u1, 0}),
			tl2::create<t_vec>({u2, 0}), tl2::create<t_vec>({u2, 1})
		});
	}

	if(cyltype > 0)
	{
		const auto [disk_vertices, disk_faces, disk_normals, disk_uvs]
			= create_disk<t_vec, t_cont>(r, num_points);

		// bottom lid
		// vertex indices have to be adapted for merging
		std::size_t vert_start_idx = vertices.size();
		const t_vec top = tl2::create<t_vec>({ 0, 0, h*t_real(0.5) });

		for(const auto& disk_vert : disk_vertices)
			vertices.push_back(disk_vert - top);

		auto disk_faces_bottom = disk_faces;
		for(auto& disk_face : disk_faces_bottom)
		{
			for(auto& disk_face_idx : disk_face)
				disk_face_idx += vert_start_idx;
			std::reverse(disk_face.begin(), disk_face.end());
		}
		faces.insert(faces.end(), disk_faces_bottom.begin(), disk_faces_bottom.end());

		for(const auto& normal : disk_normals)
			normals.push_back(-normal);

		uvs.insert(uvs.end(), disk_uvs.begin(), disk_uvs.end());


		vert_start_idx = vertices.size();

		if(cyltype == 1)	// top lid
		{
			for(const auto& disk_vert : disk_vertices)
				vertices.push_back(disk_vert + top);

			auto disk_faces_top = disk_faces;
			for(auto& disk_face : disk_faces_top)
				for(auto& disk_face_idx : disk_face)
					disk_face_idx += vert_start_idx;
			faces.insert(faces.end(), disk_faces_top.begin(), disk_faces_top.end());

			for(const auto& normal : disk_normals)
				normals.push_back(normal);

			uvs.insert(uvs.end(), disk_uvs.begin(), disk_uvs.end());
		}
		else if(cyltype == 2)	// arrow top
		{
			// no need to cap the arrow if the radii are equal
			bool bConeCap = !equals<t_real>(r, arrow_r);

			const auto [cone_vertices, cone_faces, cone_normals, cone_uvs] =
				tl2::create_cone<t_vec, t_cont>(arrow_r, arrow_h, bConeCap, num_points);

			for(const auto& cone_vert : cone_vertices)
				vertices.push_back(cone_vert + top);

			auto cone_faces_top = cone_faces;
			for(auto& cone_face : cone_faces_top)
				for(auto& cone_face_idx : cone_face)
					cone_face_idx += vert_start_idx;
			faces.insert(faces.end(), cone_faces_top.begin(), cone_faces_top.end());

			for(const auto& normal : cone_normals)
				normals.push_back(normal);

			uvs.insert(uvs.end(), cone_uvs.begin(), cone_uvs.end());
		}
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a cuboid
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cuboid(typename t_vec::value_type lx = 1, typename t_vec::value_type ly = 1,
	typename t_vec::value_type lz = 1)
requires is_vec<t_vec>
{
	t_cont<t_vec> vertices =
	{
		create<t_vec>({ +lx, -ly, -lz }),	// vertex 0
		create<t_vec>({ -lx, -ly, -lz }),	// vertex 1
		create<t_vec>({ -lx, +ly, -lz }),	// vertex 2
		create<t_vec>({ +lx, +ly, -lz }),	// vertex 3

		create<t_vec>({ -lx, -ly, +lz }),	// vertex 4
		create<t_vec>({ +lx, -ly, +lz }),	// vertex 5
		create<t_vec>({ +lx, +ly, +lz }),	// vertex 6
		create<t_vec>({ -lx, +ly, +lz }),	// vertex 7
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 0, 1, 2, 3 },	// -z face
		{ 4, 5, 6, 7 },	// +z face
		{ 1, 0, 5, 4 }, // -y face
		{ 7, 6, 3, 2 },	// +y face
		{ 1, 4, 7, 2 },	// -x face
		{ 5, 0, 3, 6 },	// +x face
	};

	t_cont<t_vec> normals =
	{
		create<t_vec>({ 0, 0, -1 }),	// -z face
		create<t_vec>({ 0, 0, +1 }),	// +z face
		create<t_vec>({ 0, -1, 0 }),	// -y face
		create<t_vec>({ 0, +1, 0 }),	// +y face
		create<t_vec>({ -1, 0, 0 }),	// -x face
		create<t_vec>({ +1, 0, 0 }),	// +x face
	};

	t_cont<t_cont<t_vec>> uvs =
	{
		// -z face
		{
			create<t_vec>({0, 0}),
			create<t_vec>({1, 0}),
			create<t_vec>({1, 1}),
			create<t_vec>({0, 1})
		},
		// +z face
		{
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0})
		},
		// -y face
		{
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0}),
		},
		// +y face
		{
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1})
		},
		// -x face
		{
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 1})
		},
		// +x face
		{
			create<t_vec>({0, 0}),
			create<t_vec>({0, 1}),
			create<t_vec>({1, 1}),
			create<t_vec>({1, 0})
		},
	};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a cube
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_cube(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	return create_cuboid<t_vec, t_cont>(l, l, l);
}


/**
 * create the faces of a icosahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_icosahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	const T g = golden<T>;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ 0, -l, -g*l }), create<t_vec>({ 0, -l, +g*l }),
		create<t_vec>({ 0, +l, -g*l }), create<t_vec>({ 0, +l, +g*l }),

		create<t_vec>({ -g*l, 0, -l }), create<t_vec>({ -g*l, 0, +l }),
		create<t_vec>({ +g*l, 0, -l }), create<t_vec>({ +g*l, 0, +l }),

		create<t_vec>({ -l, -g*l, 0 }), create<t_vec>({ -l, +g*l, 0 }),
		create<t_vec>({ +l, -g*l, 0 }), create<t_vec>({ +l, +g*l, 0 }),
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 4, 2, 0 }, { 0, 6, 10 }, { 10, 7, 1 }, { 1, 3, 5 }, { 5, 9, 4 },
		{ 7, 10, 6 }, { 6, 0, 2 }, { 2, 4, 9 }, { 9, 5, 3 }, { 3, 1, 7 },
		{ 0, 10, 8 }, { 10, 1, 8 }, { 1, 5, 8 }, { 5, 4, 8 }, { 4, 0, 8 },
		{ 3, 7, 11 }, { 7, 6, 11 }, { 6, 2, 11 }, { 2, 9, 11 }, { 9, 3, 11 },
	};


	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	normals.reserve(faces.size());
	uvs.reserve(faces.size());

	for(const auto& face : faces)
	{
		auto iter = face.begin();
		const t_vec& vec1 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec2 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec3 = *(vertices.begin() + *iter);

		const t_vec vec12 = vec2 - vec1;
		const t_vec vec13 = vec3 - vec1;

		// TODO
		uvs.emplace_back(t_cont<t_vec>{{
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}) }});

		t_vec n = tl2::cross<t_vec>({vec12, vec13});
		n /= tl2::norm<t_vec>(n);
		normals.emplace_back(std::move(n));
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a dodecahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_dodecahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;
	const T g = golden<T>;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ l, l, l }), create<t_vec>({ l, l, -l }),
		create<t_vec>({ l, -l, l }), create<t_vec>({ l, -l, -l }),

		create<t_vec>({ -l, l, l }), create<t_vec>({ -l, l, -l }),
		create<t_vec>({ -l, -l, l }), create<t_vec>({ -l, -l, -l }),

		create<t_vec>({ 0, T{l}/g, g*l }), create<t_vec>({ 0, T{l}/g, -g*l }),
		create<t_vec>({ 0, -T{l}/g, g*l }), create<t_vec>({ 0, -T{l}/g, -g*l }),

		create<t_vec>({ g*l, 0, T{l}/g }), create<t_vec>({ g*l, 0, -T{l}/g }),
		create<t_vec>({ -g*l, 0, T{l}/g }), create<t_vec>({ -g*l, 0, -T{l}/g }),

		create<t_vec>({ T{l}/g, g*l, 0 }), create<t_vec>({ T{l}/g, -g*l, 0 }),
		create<t_vec>({ -T{l}/g, g*l, 0 }), create<t_vec>({ -T{l}/g, -g*l, 0 }),
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 0, 16, 18, 4, 8 }, { 0, 8, 10, 2, 12 }, { 0, 12, 13, 1, 16 },
		{ 1, 9, 5, 18, 16 }, { 1, 13, 3, 11, 9 }, { 2, 17, 3, 13, 12 },
		{ 3, 17, 19, 7, 11 }, { 2, 10, 6, 19, 17 }, { 4, 14, 6, 10, 8 },
		{ 4, 18, 5, 15, 14 }, { 5, 9, 11, 7, 15 }, { 6, 14, 15, 7, 19 },
	};


	t_cont<t_vec> normals;
	t_cont<t_cont<t_vec>> uvs;

	normals.reserve(faces.size());
	uvs.reserve(faces.size());

	for(const auto& face : faces)
	{
		auto iter = face.begin();
		const t_vec& vec1 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec2 = *(vertices.begin() + *iter); std::advance(iter,1);
		const t_vec& vec3 = *(vertices.begin() + *iter);

		const t_vec vec12 = vec2 - vec1;
		const t_vec vec13 = vec3 - vec1;

		// TODO
		uvs.emplace_back(t_cont<t_vec>{{
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}),
			create<t_vec>({0, 0}) }});

		t_vec n = tl2::cross<t_vec>({vec12, vec13});
		n /= tl2::norm<t_vec>(n);
		normals.emplace_back(std::move(n));
	}

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a octahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_octahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ +l, 0, 0 }),	// vertex 0
		create<t_vec>({ 0, +l, 0 }),	// vertex 1
		create<t_vec>({ 0, 0, +l }),	// vertex 2

		create<t_vec>({ -l, 0, 0 }),	// vertex 3
		create<t_vec>({ 0, -l, 0 }),	// vertex 4
		create<t_vec>({ 0, 0, -l }),	// vertex 5
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 2, 0, 1 }, { 0, 5, 1 }, { 5, 3, 1 }, { 3, 2, 1 },	// upper half
		{ 0, 2, 4 }, { 5, 0, 4 }, { 3, 5, 4 }, { 2, 3, 4 },	// lower half
	};


	const T len = std::sqrt(3);

	t_cont<t_vec> normals =
	{
		create<t_vec>({ +1/len, +1/len, +1/len }),
		create<t_vec>({ +1/len, +1/len, -1/len }),
		create<t_vec>({ -1/len, +1/len, -1/len }),
		create<t_vec>({ -1/len, +1/len, +1/len }),

		create<t_vec>({ +1/len, -1/len, +1/len }),
		create<t_vec>({ +1/len, -1/len, -1/len }),
		create<t_vec>({ -1/len, -1/len, -1/len }),
		create<t_vec>({ -1/len, -1/len, +1/len }),
	};

	t_cont<t_cont<t_vec>> uvs =
	{
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },

		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
	};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * create the faces of a tetrahedron
 * @returns [vertices, face vertex indices, face normals, face uvs]
 * @see https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_cont<t_vec>, t_cont<t_cont<std::size_t>>, t_cont<t_vec>, t_cont<t_cont<t_vec>>>
create_tetrahedron(typename t_vec::value_type l = 1)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_cont<t_vec> vertices =
	{
		create<t_vec>({ -l, -l, +l }),	// vertex 0
		create<t_vec>({ +l, +l, +l }),	// vertex 1
		create<t_vec>({ -l, +l, -l }),	// vertex 2
		create<t_vec>({ +l, -l, -l }),	// vertex 3
	};

	t_cont<t_cont<std::size_t>> faces =
	{
		{ 1, 2, 0 }, { 2, 1, 3 },	// connected to upper edge
		{ 0, 3, 1 }, { 3, 0, 2 },	// connected to lower edge
	};


	const T len = std::sqrt(3);

	t_cont<t_vec> normals =
	{
		create<t_vec>({ -1/len, +1/len, +1/len }),
		create<t_vec>({ +1/len, +1/len, -1/len }),
		create<t_vec>({ +1/len, -1/len, +1/len }),
		create<t_vec>({ -1/len, -1/len, -1/len }),
	};

	t_cont<t_cont<t_vec>> uvs =
	{
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
		{ create<t_vec>({0,0}), create<t_vec>({1,0}), create<t_vec>({0.5,1}) },
	};

	return std::make_tuple(vertices, faces, normals, uvs);
}


/**
 * calculates the bounding sphere of a collection of vertices
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, typename t_vec::value_type> bounding_sphere(
	const t_cont<t_vec>& verts)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	t_real rad{};
	t_vec center = tl2::mean<t_vec, t_cont>(verts);

	for(const t_vec& vec : verts)
	{
		t_vec vecCur = vec-center;

		t_real dot = tl2::inner<t_vec>(vecCur, vecCur);
		rad = std::max(rad, dot);
	}

	rad = std::sqrt(rad);
	return std::make_tuple(center, rad);
}


/**
 * calculates the bounding box of a collection of vertices
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, t_vec> bounding_box(
	const t_cont<t_vec>& verts,
	const t_vec* oldmin = nullptr, const t_vec* oldmax = nullptr)
requires is_vec<t_vec>
{
	if(verts.size() == 0)
		return std::make_tuple(t_vec{}, t_vec{});

	using t_real = typename t_vec::value_type;
	const std::size_t dim = verts[0].size();
	if(dim == 0)
		return std::make_tuple(t_vec{}, t_vec{});

	t_vec min, max;
	if(oldmin && oldmax)
	{
		min = *oldmin;
		max = *oldmax;
	}
	else
	{
		min = create<t_vec>(dim);
		max = create<t_vec>(dim);

		for(std::size_t i = 0; i < dim; ++i)
		{
			min[i] = std::numeric_limits<t_real>::max();
			max[i] = std::numeric_limits<t_real>::lowest();
		}
	}

	for(const t_vec& vert : verts)
	{
		for(std::size_t i = 0; i < dim; ++i)
		{
			min[i] = std::min(min[i], vert[i]);
			max[i] = std::max(max[i], vert[i]);
		}
	}

	return std::make_tuple(min, max);
}


/**
 * calculates the bounding box of a collection of vertices
 * @see https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
 */
template<class t_vec, template<class...> class t_cont = std::vector>
std::tuple<t_vec, t_vec> bounding_box(
	const t_cont<t_cont<t_vec>>& allverts, std::size_t dim = 3)
requires is_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	t_vec min = create<t_vec>(dim);
	t_vec max = create<t_vec>(dim);

	for(std::size_t i = 0; i < dim; ++i)
	{
		min[i] = std::numeric_limits<t_real>::max();
		max[i] = std::numeric_limits<t_real>::lowest();
	}

	for(const t_cont<t_vec>& verts : allverts)
		std::tie(min, max) = bounding_box<t_vec, t_cont>(verts, &min, &max);

	return std::make_tuple(min, max);
}


/**
 * calculates the bounding box of a sphere
 */
template<class t_vec, class t_real=typename t_vec::value_type> requires is_vec<t_vec>
std::tuple<t_vec, t_vec> sphere_bounding_box(const t_vec& pos, t_real rad)
{
	const std::size_t dim = pos.size();
	if(dim == 0)
		return std::make_tuple(t_vec{}, t_vec{});

	t_vec min = create<t_vec>(dim);
	t_vec max = create<t_vec>(dim);

	for(std::size_t i = 0; i < dim; ++i)
	{
		min[i] = pos[i] - rad;
		max[i] = pos[i] + rad;
	}

	return std::make_tuple(min, max);
}


/**
 * calculates the bounding box of a collection of spheres
 */
template<class t_vec, template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec> sphere_bounding_box(
	const t_cont<std::tuple<t_vec, t_real>>& spheres, std::size_t dim = 3)
requires is_vec<t_vec>
{
	t_vec min = create<t_vec>(dim);
	for(std::size_t i = 0; i < dim; ++i)
		min[i] = std::numeric_limits<t_real>::max();
	t_vec max = -min;

	for(const auto& sphere : spheres)
	{
		auto [spheremin, spheremax] = sphere_bounding_box<t_vec, t_real>(
			std::get<0>(sphere), std::get<1>(sphere));

		for(std::size_t i = 0; i < dim; ++i)
		{
			min[i] = std::min(min[i], spheremin[i]);
			max[i] = std::max(max[i], spheremax[i]);
		}
	}

	return std::make_tuple(min, max);
}


/**
 * checks if a point is inside a bounding box
 * @see https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
 */
template<class t_vec> requires is_vec<t_vec>
bool in_bounding_box(
	const t_vec& pt, const std::tuple<t_vec, t_vec>& bb)
{
	const std::size_t dim = pt.size();

	const t_vec& min = std::get<0>(bb);
	const t_vec& max = std::get<1>(bb);

	for(std::size_t i = 0; i < dim; ++i)
	{
		if(pt[i] < min[i] || pt[i] > max[i])
			return false;
	}

	return true;
}


/**
 * checks for bounding box intersection
 * @see https://developer.mozilla.org/en-US/docs/Games/Techniques/3D_collision_detection
 */
template<class t_vec> requires is_vec<t_vec>
bool collide_bounding_boxes(
	const std::tuple<t_vec, t_vec>& bb1,
	const std::tuple<t_vec, t_vec>& bb2)
{
	// invalid bounding boxes?
	if(std::get<0>(bb1).size() == 0 || std::get<1>(bb1).size() == 0 ||
		std::get<0>(bb2).size() == 0 || std::get<1>(bb2).size() == 0)
		return false;

	const std::size_t dim = std::min(std::get<0>(bb1).size(), std::get<0>(bb2).size());
	for(std::size_t i = 0; i < dim; ++i)
	{
		if(std::get<0>(bb1)[i] > std::get<1>(bb2)[i])
			return false;
		if(std::get<0>(bb2)[i] > std::get<1>(bb1)[i])
			return false;
	}
	return true;
}
// ----------------------------------------------------------------------------

}

#endif
