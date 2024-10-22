/**
 * Loads tabulated spacegroups
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2016
 * @license GPLv2
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

#ifndef __SG_TAB_IMPL_H__
#define __SG_TAB_IMPL_H__

#include <sstream>
#include "libs/globals.h"	// find_resource


namespace xtl {

template<class t_real>
SpaceGroups<t_real>::SpaceGroups(const std::string& strFile, const std::string& strXmlRoot)
{
	// load general space group list
	m_bOk = LoadSpaceGroups(strFile, 1, strXmlRoot);

	// load custom-defined space groups
	LoadSpaceGroups("res/data/sg_user.xml", 0 /*, strXmlRoot*/);
	LoadSpaceGroups("data/sg_user.xml", 0 /*, strXmlRoot*/);
}


template<class t_real>
SpaceGroups<t_real>::~SpaceGroups() {}


template<class t_real>
const SpaceGroup<t_real>* SpaceGroups<t_real>::Find(const std::string& strSG) const
{
	if(!IsOk()) return nullptr;

	// search space groups by name
	typename t_mapSpaceGroups::const_iterator iterSG = g_mapSpaceGroups.find(strSG);
	if(iterSG != g_mapSpaceGroups.end())
		return &iterSG->second;

	// search space groups by number
	if(tl::str_is_digits(strSG))
	{
		unsigned int iSGNr = tl::str_to_var<unsigned int>(strSG);

		const t_vecSpaceGroups* pVec = get_space_groups_vec();
		for(const SpaceGroup<t_real>* pSG : *pVec)
		{
			if(pSG->GetNr() == iSGNr)
				return pSG;
		}
	}

	// nothing found
	return nullptr;
}


template<class t_real>
bool SpaceGroups<t_real>::LoadSpaceGroups(const std::string& strFile, bool bMandatory, const std::string& strXmlRoot)
{
	using t_mat = typename SpaceGroup<t_real>::t_mat;
	//using t_vec = typename SpaceGroup<t_real>::t_vec;

	std::string strTabFile = find_resource(strFile, bMandatory);
	if(strTabFile == "")
		return false;
	tl::log_debug("Loading space groups from file \"", strTabFile, "\".");

	tl::Prop<std::string> xml;
	if(!xml.Load(strTabFile.c_str(), tl::PropType::XML))
		return false;

	//unsigned int iNumSGs = xml.Query<unsigned int>(strXmlRoot + "/sgroups/num_groups", 230);
	typedef typename t_mapSpaceGroups::value_type t_val;

	unsigned int iSg = 0;
	while(1)
	{
		std::ostringstream ostrGroup;
		ostrGroup << "sgroups/group_" << iSg;
		std::string strGroup = ostrGroup.str();

		if(!xml.Exists(strGroup.c_str()))
			break;

		unsigned int iSgNr = xml.Query<unsigned int>(strXmlRoot + "/" + strGroup + "/number");
		std::string strName = tl::trimmed(xml.Query<std::string>(strXmlRoot + "/" + strGroup + "/name"));
		std::string strLaue = tl::trimmed(xml.Query<std::string>(strXmlRoot + "/" + strGroup + "/lauegroup"));
		unsigned int iNumTrafos = xml.Query<unsigned int>(strXmlRoot + "/" + strGroup + "/num_trafos", 0);

		std::vector<t_mat> vecTrafos;
		std::vector<unsigned int> vecInvTrafos, vecPrimTrafos, vecCenterTrafos, vecTrans;
		vecTrafos.reserve(iNumTrafos);

		unsigned int iTrafo = 0;
		while(true)
		{
			std::ostringstream ostrTrafo;
			ostrTrafo << strGroup << "/trafo_" << iTrafo;
			std::string strTrafo = ostrTrafo.str();

			if(!xml.Exists(strTrafo.c_str()))
				break;

			std::string strTrafoVal = xml.Query<std::string>(strXmlRoot + "/" + strTrafo);
			std::pair<std::string, std::string> pairSg = tl::split_first(strTrafoVal, std::string(";"), 1);

			std::istringstream istrMat(pairSg.first);
			t_mat mat;
			istrMat >> mat;

			for(typename std::string::value_type c : pairSg.second)
			{
				if(std::tolower(c)=='p') vecPrimTrafos.push_back(iTrafo);
				if(std::tolower(c)=='i') vecInvTrafos.push_back(iTrafo);
				if(std::tolower(c)=='c') vecCenterTrafos.push_back(iTrafo);
				if(std::tolower(c)=='t') vecTrans.push_back(iTrafo);
			}

			vecTrafos.push_back(std::move(mat));
			++iTrafo;
		}

		if(iNumTrafos != 0 && iTrafo != iNumTrafos)
		{
			tl::log_warn(iTrafo, " transformations were given, but ", iNumTrafos,
				" were expected in space group ", iSgNr);
		}

		SpaceGroup<t_real> sg;
		sg.SetNr(iSgNr);
		sg.SetName(strName);
		sg.SetLaueGroup(strLaue);
		sg.SetTrafos(std::move(vecTrafos));
		sg.SetInvTrafos(std::move(vecInvTrafos));
		sg.SetPrimTrafos(std::move(vecPrimTrafos));
		sg.SetCenterTrafos(std::move(vecCenterTrafos));
		sg.SetTransTrafos(std::move(vecTrans));

		auto pairSG = g_mapSpaceGroups.insert(t_val(sg.GetName(), std::move(sg)));
		g_vecSpaceGroups.push_back(&pairSG.first->second);
		++iSg;
	}

	// sort by sg number
	std::sort(g_vecSpaceGroups.begin(), g_vecSpaceGroups.end(),
		[](const SpaceGroup<t_real>* sg1, const SpaceGroup<t_real>* sg2) -> bool
		{ return sg1->GetNr() < sg2->GetNr(); });

	if(s_strSrc == "")
		s_strSrc = xml.Query<std::string>(strXmlRoot + "/sgroups/source", "");
	if(s_strUrl == "")
		s_strUrl = xml.Query<std::string>(strXmlRoot + "/sgroups/source_url", "");


	if(g_vecSpaceGroups.size() < 230)
		tl::log_warn("Less than 230 space groups are defined!");

	return true;
}


template<class t_real>
std::shared_ptr<const SpaceGroups<t_real>> SpaceGroups<t_real>::GetInstance(const char *pcFile)
{
	std::lock_guard<std::mutex> _guard(s_mutex);

	if(!s_inst)
	{
		s_inst = std::shared_ptr<SpaceGroups>(new SpaceGroups(
			pcFile ? pcFile : "res/data/sgroups.xml", ""));
	}

	return s_inst;
}


template<class t_real>
const typename SpaceGroups<t_real>::t_mapSpaceGroups* SpaceGroups<t_real>::get_space_groups() const
{
	if(g_mapSpaceGroups.empty())
		tl::log_warn("Space Groups have not been initialised properly.");

	return &g_mapSpaceGroups;
}


template<class t_real>
const typename SpaceGroups<t_real>::t_vecSpaceGroups* SpaceGroups<t_real>::get_space_groups_vec() const
{
	if(g_vecSpaceGroups.empty())
		tl::log_warn("Space Groups have not been initialised properly.");

	return &g_vecSpaceGroups;
}


template<class t_real>
const std::string& SpaceGroups<t_real>::get_sgsource(bool bUrl) const
{
	return bUrl ? s_strUrl : s_strSrc ;
}


// statics
template<class t_real> std::shared_ptr<SpaceGroups<t_real>> SpaceGroups<t_real>::s_inst = nullptr;
template<class t_real> std::mutex SpaceGroups<t_real>::s_mutex;


}
#endif
