/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc The magdyn library implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS2_MAGDYN_STRUCTS_H__
#define __TLIBS2_MAGDYN_STRUCTS_H__

#include <vector>
#include <array>
#include <string>

#include "../maths.h"



namespace tl2_mag {

// ----------------------------------------------------------------------------
// input- and output struct templates
// ----------------------------------------------------------------------------

using t_strarr3 = std::array<std::string, 3>;
using t_strarr33 = std::array<std::array<std::string, 3>, 3>;



/**
 * magnetic sites
 */
template<class t_mat, class t_vec, class t_vec_real,
	class t_size, class t_real = typename t_vec_real::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
struct t_MagneticSite
{
	// ------------------------------------------------------------------------
	// input properties
	std::string name{};          // identifier
	t_size sym_idx{};            // groups positions belonging to the same symmetry group (0: none)

	t_strarr3 pos{};             // magnetic site position

	t_strarr3 spin_dir{};        // spin direction
	t_strarr3 spin_ortho{};      // spin orthogonal vector

	std::string spin_mag{};      // spin magnitude
	t_mat g_e{};                 // electron g factor
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// calculated properties
	t_vec_real pos_calc{};       // magnetic site position

	t_vec_real spin_dir_calc{};  // spin vector
	t_vec trafo_z_calc{};        // trafo z vector (3rd column in trafo matrix)
	t_vec trafo_plane_calc{};    // trafo orthogonal plane (1st and 2nd coumns)
	t_vec trafo_plane_conj_calc{};

	t_vec ge_trafo_z_calc{};     // g_e * trafo z vector
	t_vec ge_trafo_plane_calc{}; // g_e * trafo orthogonal plane
	t_vec ge_trafo_plane_conj_calc{};

	t_real spin_mag_calc{};      // spin magnitude
	// ------------------------------------------------------------------------
};



/**
 * couplings between magnetic sites
 */
template<class t_mat, class t_vec, class t_vec_real,
	class t_size, class t_cplx = typename t_mat::value_type,
	class t_real = typename t_vec_real::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
struct t_ExchangeTerm
{
	// ------------------------------------------------------------------------
	// input properties (parsable expressions)
	std::string name{};          // identifier
	t_size sym_idx{};            // groups couplings belonging to the same symmetry group (0: none)

	std::string site1{};         // name of first magnetic site
	std::string site2{};         // name of second magnetic site
	t_strarr3 dist{};            // distance between unit cells

	std::string J{};             // Heisenberg interaction
	t_strarr3 dmi{};             // Dzyaloshinskij-Moriya interaction
	t_strarr33 Jgen{};           // general exchange interaction
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// calculated properties
	t_size site1_calc{};         // index of first magnetic site
	t_size site2_calc{};         // index of second magnetic site
	t_vec_real dist_calc{};      // distance between unit cells
	t_real length_calc{};        // length of the coupling (in lab units)

	t_cplx J_calc{};             // Heisenberg interaction
	t_vec dmi_calc{};            // Dzyaloshinskij-Moriya interaction
	t_mat Jgen_calc{};           // general exchange interaction
	// ------------------------------------------------------------------------
};



/**
 * terms related to an external magnetic field
 */
template<class t_vec_real, class t_real>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_vec<t_vec_real>
#endif
struct t_ExternalField
{
	bool align_spins{};          // align spins along external field

	t_vec_real dir{};            // field direction
	t_real mag{};                // field magnitude
};



/**
 * eigenenergies and spin-spin correlation matrix
 */
template<class t_mat, class t_real, class t_cplx = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
struct t_EnergyAndWeight
{
	t_real E{};

	// full dynamical structure factor
	t_mat S{};
	t_cplx S_sum{};
	t_real weight_full{};

	// projected dynamical structure factor for neutron scattering
	t_mat S_perp{};
	t_cplx S_perp_sum{};
	t_real weight{};
};



template<class t_mat, class t_real, class t_cplx = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
struct t_SofQE
{
	t_real h, k, l;
	std::vector<t_EnergyAndWeight<t_mat, t_real, t_cplx>> E_and_S;
};



/**
 * variables for the expression parser
 */
template<class t_cplx>
struct t_Variable
{
	std::string name{};
	t_cplx value{};
};
// ----------------------------------------------------------------------------

}
#endif
