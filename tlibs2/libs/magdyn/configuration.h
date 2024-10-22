/**
 * tlibs2 -- magnetic dynamics -- loading and saving
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
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

#ifndef __TLIBS2_MAGDYN_CFG_FILE_H__
#define __TLIBS2_MAGDYN_CFG_FILE_H__

#include <unordered_set>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "../algos.h"
#include "../maths.h"
#include "../str.h"

#include "magdyn.h"



// --------------------------------------------------------------------
// loading and saving of the configuration
// --------------------------------------------------------------------
/**
 * load a configuration from a file
 */
MAGDYN_TEMPL
bool MAGDYN_INST::Load(const std::string& filename)
{
	try
	{
		// properties tree
		boost::property_tree::ptree node;

		// read xml file
		std::ifstream ifstr{filename};
		boost::property_tree::read_xml(ifstr, node);

		const auto &magdyn = node.get_child("magdyn");
		return Load(magdyn);
	}
	catch(const std::exception& ex)
	{
		CERR_OPT << "Magdyn error: Could not load \""
			<< filename << "\"."
			<< " Reason: " << ex.what()
			<< std::endl;

		return false;
	}
}



/**
 * save the configuration to a file
 */
MAGDYN_TEMPL
bool MAGDYN_INST::Save(const std::string& filename) const
{
	// properties tree
	boost::property_tree::ptree node;

	if(!Save(node))
		return false;

	boost::property_tree::ptree root_node;
	root_node.put_child("magdyn", node);

	// write xml file
	std::ofstream ofstr{filename};
	if(!ofstr)
		return false;

	ofstr.precision(m_prec);
	boost::property_tree::write_xml(ofstr, root_node,
		boost::property_tree::xml_writer_make_settings(
			'\t', 1, std::string{"utf-8"}));
	return true;
}



/**
 * load a configuration from a property tree
 */
MAGDYN_TEMPL
bool MAGDYN_INST::Load(const boost::property_tree::ptree& node)
{
	// check signature
	if(auto optInfo = node.get_optional<std::string>("meta.info");
		!optInfo || !(*optInfo == std::string{"magdyn_tool"}))
	{
		return false;
	}

	Clear();

	// variables
	if(auto vars = node.get_child_optional("variables"); vars)
	{
		m_variables.reserve(vars->size());

		for(const auto &var : *vars)
		{
			auto name = var.second.get_optional<std::string>("name");
			if(!name)
				continue;

			Variable variable;
			variable.name = *name;
			variable.value = var.second.get<t_cplx>("value", 0.);

			AddVariable(std::move(variable));
		}
	}

	// magnetic sites
	if(auto sites = node.get_child_optional("atom_sites"); sites)
	{
		m_sites.reserve(sites->size());
		std::unordered_set<std::string> seen_names;
		t_size unique_name_counter = 1;

		for(const auto &site : *sites)
		{
			MagneticSite magnetic_site;

			magnetic_site.name = site.second.get<std::string>("name", "");
			if(magnetic_site.name == "")
				magnetic_site.name = "site_" + tl2::var_to_str(GetMagneticSitesCount());

			if(seen_names.find(magnetic_site.name) != seen_names.end())
			{
				// try to create a unique name
				magnetic_site.name += "_" + tl2::var_to_str(unique_name_counter, m_prec);
				++unique_name_counter;
			}
			else
			{
				seen_names.insert(magnetic_site.name);
			}

			magnetic_site.pos_calc = tl2::create<t_vec_real>(
			{
				site.second.get<t_real>("position_x", 0.),
				site.second.get<t_real>("position_y", 0.),
				site.second.get<t_real>("position_z", 0.),
			});

			magnetic_site.sym_idx = site.second.get<t_size>("symmetry_index", 0);

			magnetic_site.pos[0] = site.second.get<std::string>("position_x", "0");
			magnetic_site.pos[1] = site.second.get<std::string>("position_y", "0");
			magnetic_site.pos[2] = site.second.get<std::string>("position_z", "0");

			magnetic_site.spin_dir[0] = site.second.get<std::string>("spin_x", "0");
			magnetic_site.spin_dir[1] = site.second.get<std::string>("spin_y", "0");
			magnetic_site.spin_dir[2] = site.second.get<std::string>("spin_z", "1");

			magnetic_site.spin_ortho[0] = site.second.get<std::string>("spin_ortho_x", "");
			magnetic_site.spin_ortho[1] = site.second.get<std::string>("spin_ortho_y", "");
			magnetic_site.spin_ortho[2] = site.second.get<std::string>("spin_ortho_z", "");

			magnetic_site.spin_mag = site.second.get<std::string>("spin_magnitude", "1");

			if(magnetic_site.g_e.size1() == 0 || magnetic_site.g_e.size2() == 0)
				magnetic_site.g_e = tl2::g_e<t_real> * tl2::unit<t_mat>(3);

			for(std::uint8_t i = 0; i < 3; ++i)
			for(std::uint8_t j = 0; j < 3; ++j)
			{
				std::string g_name = std::string{"gfactor_"} +
					std::string{g_comp_names[i]} +
					std::string{g_comp_names[j]};

				if(auto g_comp = site.second.get_optional<t_cplx>(g_name); g_comp)
					magnetic_site.g_e(i, j) = *g_comp;
			}

			AddMagneticSite(std::move(magnetic_site));
		}
	}

	// exchange terms / couplings
	if(auto terms = node.get_child_optional("exchange_terms"); terms)
	{
		m_exchange_terms.reserve(terms->size());
		std::unordered_set<std::string> seen_names;
		t_size unique_name_counter = 1;

		for(const auto &term : *terms)
		{
			ExchangeTerm exchange_term;

			exchange_term.name = term.second.get<std::string>("name", "");
			if(exchange_term.name == "")
				exchange_term.name = "coupling_" + tl2::var_to_str(GetExchangeTermsCount());

			if(seen_names.find(exchange_term.name) != seen_names.end())
			{
				// try to create a unique name
				exchange_term.name += "_" + tl2::var_to_str(unique_name_counter, m_prec);
				++unique_name_counter;
			}
			else
			{
				seen_names.insert(exchange_term.name);
			}

			exchange_term.site1_calc = term.second.get<t_size>("atom_1_index", 0);
			exchange_term.site2_calc = term.second.get<t_size>("atom_2_index", 0);

			if(auto name1 = term.second.get_optional<std::string>("atom_1_name"); name1)
			{
				// get the magnetic site index via the name
				if(const MagneticSite* sites1 = FindMagneticSite(*name1); sites1)
				{
					exchange_term.site1 = sites1->name;
				}
				else
				{
					CERR_OPT << "Magdyn error: Site 1 name \"" << *name1 << "\" "
						<< "was not found in coupling \"" << exchange_term.name
						<< "\"." << std::endl;
				}
			}
			else
			{
				// get the magnetic site name via the index
				exchange_term.site1 = GetMagneticSite(exchange_term.site1_calc).name;
			}

			if(auto name2 = term.second.get_optional<std::string>("atom_2_name"); name2)
			{
				if(const MagneticSite* sites2 = FindMagneticSite(*name2); sites2)
				{
					exchange_term.site2 = sites2->name;
				}
				else
				{
					CERR_OPT << "Magdyn error: Site 2 name \"" << *name2 << "\" "
						<< "was not found in coupling \"" << exchange_term.name
						<< "\"." << std::endl;
				}
			}
			else
			{
				// get the magnetic site name via the index
				exchange_term.site2 = GetMagneticSite(exchange_term.site2_calc).name;
			}

			exchange_term.dist_calc = tl2::create<t_vec_real>(
			{
				term.second.get<t_real>("distance_x", 0.),
				term.second.get<t_real>("distance_y", 0.),
				term.second.get<t_real>("distance_z", 0.),
			});

			exchange_term.dist[0] = term.second.get<std::string>("distance_x", "0");
			exchange_term.dist[1] = term.second.get<std::string>("distance_y", "0");
			exchange_term.dist[2] = term.second.get<std::string>("distance_z", "0");

			exchange_term.sym_idx = term.second.get<t_size>("symmetry_index", 0);

			exchange_term.J = term.second.get<std::string>("interaction", "0");

			for(std::uint8_t i = 0; i < 3; ++i)
			{
				exchange_term.dmi[i] = term.second.get<std::string>(
					std::string{"dmi_"} + std::string{g_comp_names[i]}, "0");

				for(std::uint8_t j = 0; j < 3; ++j)
				{
					exchange_term.Jgen[i][j] = term.second.get<std::string>(
						std::string{"gen_"} +
						std::string{g_comp_names[i]} +
						std::string{g_comp_names[j]}, "0");
				}
			}

			AddExchangeTerm(std::move(exchange_term));
		}
	}

	// external field
	if(auto field = node.get_child_optional("field"); field)
	{
		ExternalField thefield;

		thefield.mag = 0.;
		thefield.align_spins = false;

		thefield.dir = tl2::create<t_vec_real>(
		{
			field->get<t_real>("direction_h", 0.),
			field->get<t_real>("direction_k", 0.),
			field->get<t_real>("direction_l", 1.),
		});

		if(auto optVal = field->get_optional<t_real>("magnitude"))
			thefield.mag = *optVal;

		if(auto optVal = field->get_optional<bool>("align_spins"))
			thefield.align_spins = *optVal;

		SetExternalField(thefield);
	}

	// temperature
	m_temperature = node.get<t_real>("temperature", -1.);

	// form factor
	SetMagneticFormFactor(node.get<std::string>("magnetic_form_factor", ""));

	// ordering vector
	if(auto ordering = node.get_child_optional("ordering"); ordering)
	{
		t_vec_real ordering_vec = tl2::create<t_vec_real>(
		{
			ordering->get<t_real>("h", 0.),
			ordering->get<t_real>("k", 0.),
			ordering->get<t_real>("l", 0.),
		});

		SetOrderingWavevector(ordering_vec);
	}

	// rotation axis
	if(auto axis = node.get_child_optional("rotation_axis"); axis)
	{
		t_vec_real rotaxis = tl2::create<t_vec_real>(
		{
			axis->get<t_real>("h", 1.),
			axis->get<t_real>("k", 0.),
			axis->get<t_real>("l", 0.),
		});

		SetRotationAxis(rotaxis);
	}

	// crystal lattice
	t_real a = node.get<t_real>("xtal.a", 5.);
	t_real b = node.get<t_real>("xtal.b", 5.);
	t_real c = node.get<t_real>("xtal.c", 5.);
	t_real alpha = tl2::d2r<t_real>(node.get<t_real>("xtal.alpha", 90.));
	t_real beta = tl2::d2r<t_real>(node.get<t_real>("xtal.beta", 90.));
	t_real gamma = tl2::d2r<t_real>(node.get<t_real>("xtal.gamma", 90.));
	SetCrystalLattice(a, b, c, alpha, beta, gamma);

	// scattering plane
	t_real ah = node.get<t_real>("xtal.plane_ah", 1.);
	t_real ak = node.get<t_real>("xtal.plane_ak", 0.);
	t_real al = node.get<t_real>("xtal.plane_al", 0.);
	t_real bh = node.get<t_real>("xtal.plane_bh", 0.);
	t_real bk = node.get<t_real>("xtal.plane_bk", 1.);
	t_real bl = node.get<t_real>("xtal.plane_bl", 0.);
	SetScatteringPlane(ah, ak, al, bh, bk, bl);

	CalcExternalField();
	CalcMagneticSites();
	CalcExchangeTerms();

	return true;
}



/**
 * save the configuration to a property tree
 */
MAGDYN_TEMPL
bool MAGDYN_INST::Save(boost::property_tree::ptree& node) const
{
	// write signature
	node.put<std::string>("meta.info", "magdyn_tool");
	node.put<std::string>("meta.date", tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()));
	node.put<std::string>("meta.doi_tlibs", "https://doi.org/10.5281/zenodo.5717779");

	// external field
	if(m_field.dir.size() == 3)
	{
		node.put<t_real>("field.direction_h", m_field.dir[0]);
		node.put<t_real>("field.direction_k", m_field.dir[1]);
		node.put<t_real>("field.direction_l", m_field.dir[2]);
	}
	node.put<t_real>("field.magnitude", m_field.mag);
	node.put<bool>("field.align_spins", m_field.align_spins);

	// ordering vector
	if(m_ordering.size() == 3)
	{
		node.put<t_real>("ordering.h", m_ordering[0]);
		node.put<t_real>("ordering.k", m_ordering[1]);
		node.put<t_real>("ordering.l", m_ordering[2]);
	}

	// rotation axis
	if(m_rotaxis.size() == 3)
	{
		node.put<t_real>("rotation_axis.h", m_rotaxis[0]);
		node.put<t_real>("rotation_axis.k", m_rotaxis[1]);
		node.put<t_real>("rotation_axis.l", m_rotaxis[2]);
	}

	// temperature
	node.put<t_real>("temperature", m_temperature);

	// form factor
	node.put<std::string>("magnetic_form_factor", GetMagneticFormFactor());

	// variables
	for(const Variable& var : GetVariables())
	{
		boost::property_tree::ptree itemNode;
		itemNode.put<std::string>("name", var.name);
		itemNode.put<t_cplx>("value", var.value);

		node.add_child("variables.variable", itemNode);
	}

	// magnetic sites
	for(const MagneticSite& site : GetMagneticSites())
	{
		boost::property_tree::ptree itemNode;
		itemNode.put<std::string>("name", site.name);

		itemNode.put<std::string>("position_x", site.pos[0]);
		itemNode.put<std::string>("position_y", site.pos[1]);
		itemNode.put<std::string>("position_z", site.pos[2]);

		itemNode.put<t_size>("symmetry_index", site.sym_idx);

		itemNode.put<std::string>("spin_x", site.spin_dir[0]);
		itemNode.put<std::string>("spin_y", site.spin_dir[1]);
		itemNode.put<std::string>("spin_z", site.spin_dir[2]);

		itemNode.put<std::string>("spin_ortho_x", site.spin_ortho[0]);
		itemNode.put<std::string>("spin_ortho_y", site.spin_ortho[1]);
		itemNode.put<std::string>("spin_ortho_z", site.spin_ortho[2]);

		itemNode.put<std::string>("spin_magnitude", site.spin_mag);

		for(std::uint8_t i = 0; i < site.g_e.size1(); ++i)
		for(std::uint8_t j = 0; j < site.g_e.size2(); ++j)
		{
			itemNode.put<t_cplx>(std::string{"gfactor_"} +
				std::string{g_comp_names[i]} +
				std::string{g_comp_names[j]}, site.g_e(i, j));
		}

		node.add_child("atom_sites.site", itemNode);
	}

	// exchange terms
	for(const ExchangeTerm& term : GetExchangeTerms())
	{
		boost::property_tree::ptree itemNode;
		itemNode.put<std::string>("name", term.name);

		// save the magnetic site names and indices
		itemNode.put<t_size>("atom_1_index", term.site1_calc);
		itemNode.put<t_size>("atom_2_index", term.site2_calc);
		itemNode.put<std::string>("atom_1_name", term.site1);
		itemNode.put<std::string>("atom_2_name", term.site2);

		itemNode.put<std::string>("distance_x", term.dist[0]);
		itemNode.put<std::string>("distance_y", term.dist[1]);
		itemNode.put<std::string>("distance_z", term.dist[2]);

		itemNode.put<t_size>("symmetry_index", term.sym_idx);

		itemNode.put<std::string>("interaction", term.J);

		for(std::uint8_t i = 0; i < 3; ++i)
		{
			itemNode.put<std::string>(std::string{"dmi_"} +
				std::string{g_comp_names[i]}, term.dmi[i]);

			for(std::uint8_t j = 0; j < 3; ++j)
			{
				itemNode.put<std::string>(std::string{"gen_"} +
					std::string{g_comp_names[i]} +
					std::string{g_comp_names[j]},
					term.Jgen[i][j]);
			}
		}

		node.add_child("exchange_terms.term", itemNode);
	}

	// crystal lattice
	node.put<t_real>("xtal.a", m_xtallattice[0]);
	node.put<t_real>("xtal.b", m_xtallattice[1]);
	node.put<t_real>("xtal.c", m_xtallattice[2]);
	node.put<t_real>("xtal.alpha", tl2::r2d<t_real>(m_xtalangles[0]));
	node.put<t_real>("xtal.beta", tl2::r2d<t_real>(m_xtalangles[1]));
	node.put<t_real>("xtal.gamma", tl2::r2d<t_real>(m_xtalangles[2]));

	// scattering plane
	// x vector
	node.put<t_real>("xtal.plane_ah", m_scatteringplane[0][0]);
	node.put<t_real>("xtal.plane_ak", m_scatteringplane[0][1]);
	node.put<t_real>("xtal.plane_al", m_scatteringplane[0][2]);
	// y vector
	node.put<t_real>("xtal.plane_bh", m_scatteringplane[1][0]);
	node.put<t_real>("xtal.plane_bk", m_scatteringplane[1][1]);
	node.put<t_real>("xtal.plane_bl", m_scatteringplane[1][2]);
	// up vector (saving is optional)
	node.put<t_real>("xtal.plane_ch", m_scatteringplane[2][0]);
	node.put<t_real>("xtal.plane_ck", m_scatteringplane[2][1]);
	node.put<t_real>("xtal.plane_cl", m_scatteringplane[2][2]);

	return true;
}
// --------------------------------------------------------------------

#endif
