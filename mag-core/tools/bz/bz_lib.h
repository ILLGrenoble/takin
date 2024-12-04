/**
 * brillouin zone calculation library
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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


#ifndef __TAKIN_BZLIB_H__
#define __TAKIN_BZLIB_H__

#include <vector>
#include <optional>

#include "tlibs2/libs/maths.h"
#include "libs/symops.h"
#include "libs/loadcif.h"

#if __has_include("pathslib/libs/voronoi.h")
	#include "pathslib/libs/voronoi.h"
#else
	#include "libs/voronoi.h"
#endif


/**
 * brillouin zone calculation
 */
template<class t_mat, class t_vec, class t_real = typename t_vec::value_type>
//requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
class BZCalc
{
public:
	BZCalc() = default;
	~BZCalc() = default;


	/**
	 * clear old BZ results
	 */
	void ClearBZ()
	{
		m_vertices.clear();

		m_triags.clear();
		m_triags_idx.clear();

		m_all_triags.clear();
		m_all_triags_idx.clear();

		m_face_polygons.clear();
		m_face_norms.clear();
		m_face_dists.clear();
	}


	/**
	 * clear old symops
	 */
	void ClearSymOps()
	{
		m_symops.clear();
	}


	/**
	 * clear old peaks
	 */
	void ClearPeaks()
	{
		m_peaks.clear();
		m_peaks_invA.clear();
	}


	// --------------------------------------------------------------------------------
	// getter and setter
	// --------------------------------------------------------------------------------
	static std::size_t GetErrIdx()
	{
		return s_erridx;
	}


	void SetEps(t_real eps)
	{
		m_eps = eps;
	}


	void SetCrystalB(const t_mat& B)
	{
		m_crystB = B;
		m_crystB_ortho = orthonorm_sys<t_mat, t_vec>(m_crystB);
	}

	/**
	 * sets up a crystal lattice and angles
	 */
	void SetCrystal(t_real a, t_real b, t_real c,
		t_real alpha = 90., t_real beta = 90., t_real gamma = 90.)
	{
		t_mat crystB = tl2::B_matrix<t_mat>(a, b, c,
			tl2::d2r<t_real>(alpha), tl2::d2r<t_real>(beta), tl2::d2r<t_real>(gamma));

		SetCrystalB(crystB);
	}


	void SetPeaks(const std::vector<t_vec>& peaks)
	{
		m_peaks = peaks;
	}


	const std::vector<t_vec>& GetPeaks() const
	{
		return m_peaks;
	}


	void SetPeaksInvA(const std::vector<t_vec>& peaks)
	{
		m_peaks_invA = peaks;
	}


	const std::vector<t_vec>& GetPeaksInvA() const
	{
		return m_peaks_invA;
	}


	const std::vector<t_vec>& GetVertices() const
	{
		return m_vertices;
	}


	const std::vector<std::vector<t_vec>>& GetTriangles() const
	{
		return m_triags;
	}


	const std::vector<std::vector<std::size_t>>& GetTrianglesIndices() const
	{
		return m_triags_idx;
	}


	const std::vector<std::vector<std::size_t>>& GetFacesIndices() const
	{
		return m_face_polygons;
	}


	const std::vector<t_vec>& GetFaceNormals() const
	{
		return m_face_norms;
	}


	const std::vector<t_real>& GetFaceDistances() const
	{
		return m_face_dists;
	}


	const std::vector<t_vec>& GetAllTriangles() const
	{
		return m_all_triags;
	}


	const std::vector<std::size_t>& GetAllTrianglesIndices() const
	{
		return m_all_triags_idx;
	}


	const t_mat& GetCrystalB(bool ortho = false) const
	{
		return ortho ? m_crystB_ortho : m_crystB;
	}


	/**
	 * get the index of the (000) bragg peak
	 */
	std::size_t Get000Peak() const
	{
		if(m_idx000)
			return *m_idx000;
		return s_erridx;
	}


	void AddSymOp(
		t_real r00, t_real r01, t_real r02,
		t_real r10, t_real r11, t_real r12,
		t_real r20, t_real r21, t_real r22,
		t_real t0, t_real t1, t_real t2)
	{
		t_mat symop = tl2::unit<t_mat>(4);
		// rotation
		symop(0, 0) = r00; symop(0, 1) = r01; symop(0, 2) = r02;
		symop(1, 0) = r10; symop(1, 1) = r11; symop(1, 2) = r12;
		symop(2, 0) = r20; symop(2, 1) = r21; symop(2, 2) = r22;
		// translation
		symop(0, 3) = t0;  symop(1, 3) = t1;  symop(2, 3) = t2;

		m_symops.emplace_back(std::move(symop));
	}


	/**
	 * set up a list of symmetry operations (given by the space group)
	 * @returns number of actually set symops
	 */
	std::size_t SetSymOps(const std::vector<t_mat>& ops, bool are_centring = false)
	{
		if(are_centring)
		{
			// symops are already purely centring, add all
			m_symops = ops;
		}
		else
		{
			// only add centring ops
			m_symops.clear();
			m_symops.reserve(ops.size());

			for(const t_mat& op : ops)
			{
				if(!tl2::hom_is_centring<t_mat>(op, m_eps))
					continue;

				m_symops.push_back(op);
			}
		}

		return m_symops.size();
	}


	std::size_t SetSymOpsFromSpaceGroup(const std::string& sgname)
	{
		std::vector<t_mat> ops = get_sg_ops<t_mat, t_real>(sgname);
		return SetSymOps(ops, false);
	}
	// --------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------
	// calculations
	// --------------------------------------------------------------------------------
	/**
	 * calculate the nuclear bragg peaks in lab coordinates
	 * @returns number of created peaks
	 */
	std::size_t CalcPeaksInvA()
	{
		// calculate the peaks in lab coordinates
		m_peaks_invA.clear();
		m_peaks_invA.reserve(m_peaks.size());

		for(const t_vec& Q : m_peaks)
		{
			if(!is_reflection_allowed<t_mat, t_vec, t_real>(Q, m_symops, m_eps).first)
				continue;

			// also get the index of the (000) peak
			if(tl2::equals_0(Q, m_eps))
				m_idx000 = m_peaks_invA.size();

			m_peaks_invA.emplace_back(m_crystB * Q);
		}

		return m_peaks_invA.size();
	}


	/**
	 * create nuclear bragg peaks up to the given order
	 * @returns number of created peaks
	 */
	std::size_t CalcPeaks(int order, bool cleate_invA = false)
	{
		m_peaks.clear();
		m_peaks.reserve((2*order+1)*(2*order+1)*(2*order+1));

		for(int h = -order; h <= order; ++h)
		for(int k = -order; k <= order; ++k)
		for(int l = -order; l <= order; ++l)
		{
			m_peaks.emplace_back(
				tl2::create<t_vec>(
					{ t_real(h), t_real(k), t_real(l) }));
		}

		if(cleate_invA)
			CalcPeaksInvA();

		return m_peaks.size();
	}


	/**
	 * calculate the index of the nuclear (000) peak
	 */
	void Calc000Peak()
	{
		for(const t_vec& Q : m_peaks)
		{
			if(!tl2::equals_0(Q, m_eps))
				continue;
			m_idx000 = m_peaks_invA.size();
			break;
		}
	}


	/**
	 * calculate the brillouin zone
	 */
	bool CalcBZ()
	{
		ClearBZ();

		if(!m_idx000)
			Calc000Peak();

		// calculate the voronoi diagram's vertices
		std::tie(m_vertices, std::ignore, std::ignore) =
			geo::calc_delaunay(3, m_peaks_invA, false, false, m_idx000);
		m_vertices = tl2::remove_duplicates(m_vertices, m_eps);
		if(!m_vertices.size())
			return false;

		for(t_vec& vertex : m_vertices)
			tl2::set_eps_0(vertex, m_eps);

		// calculate the faces of the BZ
		std::tie(std::ignore, m_triags, std::ignore) =
			geo::calc_delaunay(3, m_vertices, true, false);
		if(!m_triags.size())
			return false;

		// calculate all BZ triangles
		for(std::size_t triag_idx = 0; triag_idx < m_triags.size(); ++triag_idx)
		{
			std::vector<t_vec>& bz_triag = m_triags[triag_idx];
			if(bz_triag.size() < 3)
				continue;

			// calculate face plane
			t_vec norm = tl2::cross<t_vec>(
				{ bz_triag[1] - bz_triag[0], bz_triag[2] - bz_triag[0] });
			t_real norm_len = tl2::norm<t_vec>(norm);
			if(!tl2::equals_0<t_real>(norm_len, m_eps))
				norm /= tl2::norm<t_vec>(norm);

			t_real dist = tl2::inner<t_vec>(bz_triag[0], norm);
			if(dist < 0.)
			{
				norm = -norm;
				dist = -dist;

				std::reverse(bz_triag.begin(), bz_triag.end());
			}

			// find out if we've already got this face
			bool face_found = false;
			std::size_t face_idx = 0;
			for(face_idx=0; face_idx<m_face_polygons.size(); ++face_idx)
			{
				if(tl2::equals<t_vec>(m_face_norms[face_idx], norm, m_eps)
					&& tl2::equals<t_real>(m_face_dists[face_idx], dist, m_eps))
				{
					face_found = true;
					break;
				}
			}

			// add the triangle to the face
			if(face_found)
			{
				m_face_polygons[face_idx].push_back(triag_idx);
			}
			else
			{
				m_face_polygons.emplace_back(std::vector<std::size_t>({triag_idx}));

				tl2::set_eps_0(norm, m_eps);
				tl2::set_eps_0(dist, m_eps);

				m_face_norms.emplace_back(std::move(norm));
				m_face_dists.push_back(dist);
			}

			std::vector<std::size_t> triagindices;
			for(t_vec& vert : bz_triag)
			{
				tl2::set_eps_0(vert, m_eps);

				// find index of vert among voronoi vertices
				std::ptrdiff_t voroidx = -1;
				if(auto voro_iter = std::find_if(m_vertices.begin(), m_vertices.end(),
					[&vert, this](const t_vec& vec) -> bool
					{
						return tl2::equals<t_vec>(vec, vert, m_eps);
					}); voro_iter != m_vertices.end())
				{
					voroidx = voro_iter - m_vertices.begin();
				}

				std::size_t idx = (voroidx >= 0 ? voroidx : s_erridx);
				triagindices.push_back(idx);

				m_all_triags.push_back(vert);
				m_all_triags_idx.push_back(idx);
			}  // vertices

			m_triags_idx.emplace_back(std::move(triagindices));
		}  // triangles

		return true;
	}
	// --------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------
	// output
	// --------------------------------------------------------------------------------
	/**
	 * print a description of the bz
	 */
	std::string Print(int prec = 6) const
	{
		using namespace tl2_ops;
		t_mat BorthoT = tl2::trans(GetCrystalB(true));

		std::ostringstream ostr;
		ostr.precision(prec);

#ifdef DEBUG
		ostr << "# centring symmetry operations\n";
		for(const t_mat& op : m_symops)
			ostr << op << "\n";
#endif

		// voronoi vertices forming the vertices of the bz
		// (rotated by orthonormal part of B^(-1))
		const std::vector<t_vec>& voronoiverts = GetVertices();
		ostr << "# Brillouin zone vertices (Å⁻¹, B⁻¹-rotated)\n";
		for(std::size_t idx=0; idx<voronoiverts.size(); ++idx)
		{
			t_vec voro = BorthoT * voronoiverts[idx];
			tl2::set_eps_0(voro, m_eps);

			ostr << "vertex " << idx << ": (" << voro << ")\n";
		}

		ostr << "\n# Brillouin zone vertices (Å⁻¹, raw)\n";
		for(std::size_t idx=0; idx<voronoiverts.size(); ++idx)
		{
			const t_vec& voro = voronoiverts[idx];
			ostr << "raw vertex " << idx << ": (" << voro << ")\n";
		}

		// voronoi bisectors
		const auto& bz_polys = GetTriangles();
		const auto& bz_polys_idx = GetTrianglesIndices();

		ostr << "\n# Brillouin zone polygons (Å⁻¹)\n";
		for(std::size_t idx_triag=0; idx_triag<bz_polys.size(); ++idx_triag)
		{
			const auto& triag = bz_polys[idx_triag];
			const auto& triag_idx = bz_polys_idx[idx_triag];

			ostr << "polygon " << idx_triag << ": \n";
			for(std::size_t idx_vert=0; idx_vert<triag.size(); ++idx_vert)
			{
				const t_vec& vert = triag[idx_vert];
				std::size_t voroidx = triag_idx[idx_vert];

				ostr << "\tvertex " << voroidx << ": (" << vert << ")\n";
			}
		}

		return ostr.str();
	}


	/**
	 * export a description of the bz in json format
	 */
	std::string PrintJSON(int prec = 6) const
	{
		const t_mat& B = GetCrystalB(false);
		const t_mat& Bortho = GetCrystalB(true);
		t_mat BorthoT = tl2::trans(Bortho);

		std::ostringstream ostr;
		ostr.precision(prec);

		ostr << "{\n";

		// voronoi vertices forming the vertices of the bz
		// (rotated by orthonormal part of B^(-1))
		const std::vector<t_vec>& voronoiverts = GetVertices();
		ostr << "\"vertices\" : [\n";
		for(std::size_t idx=0; idx<voronoiverts.size(); ++idx)
		{
			t_vec voro = BorthoT * voronoiverts[idx];
			tl2::set_eps_0(voro, m_eps);

			ostr << "\t[ " << voro[0] << ", " << voro[1] << ", " << voro[2] << " ]";
			if(idx < voronoiverts.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// raw voronoi vertices
		ostr << "\"vertices_nonrot\" : [\n";
		for(std::size_t idx=0; idx<voronoiverts.size(); ++idx)
		{
			const t_vec& voro = voronoiverts[idx];
			ostr << "\t[ " << voro[0] << ", " << voro[1] << ", " << voro[2] << " ]";
			if(idx < voronoiverts.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// voronoi bisectors
		const auto& bz_polys_idx = GetTrianglesIndices();
		ostr << "\"polygons\" : [\n";
		for(std::size_t idx_triag=0; idx_triag<bz_polys_idx.size(); ++idx_triag)
		{
			const auto& triag_idx = bz_polys_idx[idx_triag];

			ostr << "\t[ ";
			for(std::size_t idx_vert=0; idx_vert<triag_idx.size(); ++idx_vert)
			{
				ostr << triag_idx[idx_vert];
				if(idx_vert < triag_idx.size() - 1)
					ostr << ", ";
			}
			ostr << " ]";
			if(idx_triag < bz_polys_idx.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// faces
		const auto& bz_faces_idx = GetFacesIndices();
		ostr << "\"faces\" : [\n";
		for(std::size_t idx_triag=0; idx_triag<bz_faces_idx.size(); ++idx_triag)
		{
			const auto& triag_idx = bz_faces_idx[idx_triag];

			ostr << "\t[ ";
			for(std::size_t idx_vert=0; idx_vert<triag_idx.size(); ++idx_vert)
			{
				ostr << triag_idx[idx_vert];
				if(idx_vert < triag_idx.size() - 1)
					ostr << ", ";
			}
			ostr << " ]";
			if(idx_triag < bz_faces_idx.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// face planes
		const std::vector<t_vec>& norms = GetFaceNormals();
		const std::vector<t_real>& dists = GetFaceDistances();
		ostr << "\"face_planes\" : [\n";
		for(std::size_t idx=0; idx<norms.size(); ++idx)
		{
			const t_vec& norm = norms[idx];
			t_real dist = dists[idx];

			ostr << "\t[ " << norm[0] << ", " << norm[1] << ", " << norm[2] << ", " << dist << " ]";
			if(idx < norms.size() - 1)
				ostr << ",";
			ostr << "\n";
		}
		ostr << "],\n\n";

		// crystal B matrix
		ostr << "\"crystal_B\" : [ ";
		for(std::size_t i=0; i<B.size1(); ++i)
		{
			for(std::size_t j=0; j<B.size2(); ++j)
			{
				t_real elem = B(i, j);
				tl2::set_eps_0(elem, m_eps);
				ostr << elem;
				if(i < B.size1() - 1 || j < B.size2() - 1)
					ostr << ",";
				ostr << " ";
			}
		}
		ostr << "],\n\n";

		// orthonormal (rotation) part of crystal B matrix
		ostr << "\"crystal_B_ortho\" : [ ";
		for(std::size_t i=0; i<Bortho.size1(); ++i)
		{
			for(std::size_t j=0; j<Bortho.size2(); ++j)
			{
				t_real elem = Bortho(i, j);
				tl2::set_eps_0(elem, m_eps);
				ostr << elem;
				if(i < Bortho.size1() - 1 || j < Bortho.size2() - 1)
					ostr << ",";
				ostr << " ";
			}
		}
		ostr << "]\n";

		ostr << "}\n";
		return ostr.str();
	}
	// --------------------------------------------------------------------------------


private:
	t_real m_eps{ 1e-7 };                          // calculation epsilon

	t_mat m_crystB{tl2::unit<t_mat>(3)};           // crystal B matrix
	t_mat m_crystB_ortho{tl2::unit<t_mat>(3)};     // orthonormal part of crystal B matrix

	std::vector<t_mat> m_symops{ };                // space group centring symmetry operations
	std::vector<t_vec> m_peaks{ };                 // nuclear bragg peaks
	std::vector<t_vec> m_peaks_invA { };           // nuclear bragg peaks in lab coordinates
	std::optional<std::size_t> m_idx000{};         // index of the (000) peak

	std::vector<t_vec> m_vertices{};               // voronoi vertices

	std::vector<std::vector<t_vec>> m_triags{};    // bz triangles
	std::vector<std::vector<std::size_t>> m_triags_idx{}; // ... and the version with voronoi vertex indices

	std::vector<t_vec> m_all_triags {};            // all brillouin zone triangles
	std::vector<std::size_t> m_all_triags_idx {};  // ... and the version with voronoi vertex indices

	std::vector<std::vector<std::size_t>> m_face_polygons{};
	std::vector<t_vec> m_face_norms{};
	std::vector<t_real> m_face_dists{};

	static const std::size_t s_erridx{0xffffffff}; // index for reporting errors
};


#endif
