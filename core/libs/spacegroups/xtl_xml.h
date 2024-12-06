/**
 * libcrystal xml loader
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 8-may-2019
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __LIB_XTL_XML_H__
#define __LIB_XTL_XML_H__

#include "tlibs/log/log.h"
#include "tlibs/file/prop.h"
#include "libs/spacegroups/spacegroup.h"

#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <unordered_set>


namespace xtl {


template<class t_real, class t_vec>
t_vec get_vec(const std::string& str)
{
	std::vector<t_real> _vec;
	tl::get_tokens<t_real, std::string>(str, ",;", _vec);

	t_vec vec(_vec.size());
	for(std::size_t i=0; i<_vec.size(); ++i)
		vec[i] = _vec[i];

	return vec;
}


template<class t_real, class t_mat, class t_vec>
t_mat get_mat(const std::string& str)
{
	std::vector<std::string> _mat;
	tl::get_tokens<std::string, std::string>(str, ";", _mat);

	t_mat mat(_mat.size(), _mat.size());
	for(std::size_t i=0; i<mat.size1(); ++i)
	{
		t_vec vec = get_vec<t_real, t_vec>(_mat[i]);
		for(std::size_t j=0; j<std::min(vec.size(), mat.size2()); ++j)
			mat(i,j) = vec[j];
	}

	return mat;
}


/**
 * try to identify the space group from the given symmetry operations
 */
template<class t_real, class t_mat>
std::string find_sg_from_ops(const std::vector<t_mat>& ops)
{
	std::string strsg;
	std::unordered_set<std::string> setSG;
	for(const t_mat& mat : ops)
		setSG.emplace(tl::var_to_str(mat, g_iPrec));


	std::shared_ptr<const xtl::SpaceGroups<t_real>> _sgs = xtl::SpaceGroups<t_real>::GetInstance();
	const typename xtl::SpaceGroups<t_real>::t_vecSpaceGroups* sgs = _sgs->get_space_groups_vec();
	for(const xtl::SpaceGroup<t_real>* sg : *sgs)
	{
		std::unordered_set<std::string> setTableSG;

		const std::vector<typename xtl::SpaceGroup<t_real>::t_mat>& vecTrafos = sg->GetTrafos();
		for(const typename xtl::SpaceGroup<t_real>::t_mat& mat : vecTrafos)
			setTableSG.emplace(tl::var_to_str(mat, g_iPrec));

		if(setSG == setTableSG)
		{
			//std::cout << strsg << std::endl;
			strsg = sg->GetName();
			break;
		}
	}

	return strsg;
}


/**
 * load crystal definitions from an xml file
 */
template<class t_real, class t_vec, class t_mat>
std::tuple<
	bool,				// ok?
	t_real, t_real, t_real, 	// lattice constants
	t_real, t_real, t_real,	 	// lattice angles
	std::vector<std::string>,	// nuclei names
	std::vector<t_vec>,		// nuclei basic positions
	std::vector<std::vector<t_vec>>, // all symmetry-equivalent nuclei positions
	std::string			// space group name, if available
	>
load_xml(std::istream& istr)
{
	std::vector<std::string> vecNames;
	std::vector<t_vec> vecPos;
	std::vector<std::vector<t_vec>> vecAllPos;
	std::vector<t_mat> vecOps;
	std::string strsg;


	tl::Prop<std::string> xml;
	if(!xml.Load(istr, tl::PropType::XML))
	{
		tl::log_err("Could not open crystal XML stream.\n");
		return std::make_tuple(false,
			t_real{0},t_real{0},t_real{0}, t_real{0},t_real{0},t_real{0},
			vecNames, vecPos, vecAllPos, strsg);
	}


	const std::string strXmlRoot("xtal/");


	// lattice
	bool bOkA=0, bOkB=0, bOkC=0, bOkAlpha=0, bOkBeta=0, bOkGamma=0;
	t_real a = xml.Query<t_real>(strXmlRoot + "lattice/a", 0, &bOkA);
	t_real b = xml.Query<t_real>(strXmlRoot + "lattice/b", 0, &bOkB);
	t_real c = xml.Query<t_real>(strXmlRoot + "lattice/c", 0, &bOkC);
	t_real alpha = tl::d2r(xml.Query<t_real>(strXmlRoot + "lattice/alpha", 0, &bOkAlpha));
	t_real beta = tl::d2r(xml.Query<t_real>(strXmlRoot + "lattice/beta", 0, &bOkBeta));
	t_real gamma = tl::d2r(xml.Query<t_real>(strXmlRoot + "lattice/gamma", 0, &bOkGamma));

	if(!bOkA || !bOkB || !bOkC || !bOkAlpha || !bOkBeta || !bOkGamma)
	{
		tl::log_err("Crystal XML stream has no valid lattice definition.");
		return std::make_tuple(false, a,b,c, alpha,beta,gamma, vecNames, vecPos, vecAllPos, strsg);
	}


	// nuclei names and basic positions
	for(std::size_t iNucl = 0; ; ++iNucl)
	{
		std::string strIdx = tl::var_to_str<std::size_t>(iNucl);

		bool bOkPos = 0, bOkName = 0;
		std::string pos = xml.Query<std::string>(strXmlRoot + "atoms/pos" + strIdx, "", &bOkPos);
		std::string name = xml.Query<std::string>(strXmlRoot + "atoms/name" + strIdx, "", &bOkName);

		// end of list?
		if(!bOkPos || pos == "")
			break;

		if(!bOkName || name == "")
			tl::log_warn("Nuclei ", iNucl, " has no name defined.");


		vecPos.emplace_back(get_vec<t_real, t_vec>(pos));
		vecNames.emplace_back(std::move(name));
	}


	// symmetry operations
	for(std::size_t iOp = 0; ; ++iOp)
	{
		std::string strIdx = tl::var_to_str<std::size_t>(iOp);

		bool bOk = 0;
		std::string op = xml.Query<std::string>(strXmlRoot + "symops/op" + strIdx, "", &bOk);

		// end of list?
		if(!bOk || op == "")
			break;

		vecOps.emplace_back(get_mat<t_real, t_mat, t_vec>(op));
	}

	// try to identify space group from symmetry operations
	strsg = find_sg_from_ops<t_real, t_mat>(vecOps);


	// all symmetry-equivalent nuclei positions
	for(std::size_t iNucl = 0; ; ++iNucl)
	{
		std::string strIdxNucl = tl::var_to_str<std::size_t>(iNucl);

		// end of list ?
		if(!xml.PathExists(strXmlRoot + "generated_atoms/pos" + strIdxNucl))
			break;

		std::vector<t_vec> vecEquivPos;
		for(std::size_t iNuclEquiv = 0; ; ++iNuclEquiv)
		{
			std::string strIdxNuclEquiv = tl::var_to_str<std::size_t>(iNuclEquiv);
			bool bOkPos = 0;
			std::string pos = xml.Query<std::string>(
				strXmlRoot + "generated_atoms/pos" + strIdxNucl + "/gen" + strIdxNuclEquiv, "", &bOkPos);

			// end of list?
			if(!bOkPos || pos == "")
				break;

			vecEquivPos.emplace_back(get_vec<t_real, t_vec>(pos));
		}

		vecAllPos.emplace_back(std::move(vecEquivPos));
	}


	return std::make_tuple(true, a,b,c, alpha,beta,gamma, vecNames, vecPos, vecAllPos, strsg);
}


}

#endif
