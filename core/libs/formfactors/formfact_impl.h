/**
 * Form factor and scattering length tables
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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

#ifndef __XTL_FFACT_IMPL_H__
#define __XTL_FFACT_IMPL_H__

#include "formfact.h"
#include "tlibs/math/math.h"
#include "tlibs/math/stat.h"
#include "tlibs/log/log.h"
#include "tlibs/file/prop.h"
#include "tlibs/string/string.h"
#include "libs/globals.h"


namespace xtl {

template<typename T> std::shared_ptr<PeriodicSystem<T>> PeriodicSystem<T>::s_inst = nullptr;
template<typename T> std::mutex PeriodicSystem<T>::s_mutex;

template<typename T>
PeriodicSystem<T>::PeriodicSystem(const std::string& strFile, const std::string& strXmlRoot)
{
	std::string strTabFile = find_resource(strFile);
	tl::log_debug("Loading periodic table from file \"", strTabFile, "\".");

	tl::Prop<std::string> xml;
	if(!xml.Load(strTabFile.c_str(), tl::PropType::XML))
		return;

	const std::size_t iNumDat = xml.Query<std::size_t>(strXmlRoot + "/pte/num_elems", 0);
	if(!iNumDat)
	{
		tl::log_err("No data in periodic table of elements.");
		return;
	}

	for(std::size_t iElem=0; iElem<iNumDat; ++iElem)
	{
		elem_type elem;
		std::string strAtom = "pte/elem_" + tl::var_to_str(iElem);

		elem.strAtom = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/name", "");
		if(elem.strAtom == "")
			continue;

		elem.iNr = xml.Query<int>(strXmlRoot + "/" + strAtom + "/num", -1.);
		elem.iPeriod = xml.Query<int>(strXmlRoot + "/" + strAtom + "/period", -1.);
		elem.iGroup = xml.Query<int>(strXmlRoot + "/" + strAtom + "/group", -1.);

		elem.strOrbitals = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/orbitals", "");
		elem.strBlock = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/block", "");

		elem.dMass = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/m", -1.);

		elem.dRadCov = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/r_cov", -1.);
		elem.dRadVdW = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/r_vdW", -1.);

		elem.dEIon = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/E_ion", -1.);
		elem.dEAffin = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/E_affin", -1.);

		elem.dTMelt = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/T_melt", -1.);
		elem.dTBoil = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/T_boil", -1.);

		s_vecAtoms.push_back(std::move(elem));
	}

	s_strSrc = xml.Query<std::string>(strXmlRoot + "/pte/source", "");
	s_strSrcUrl = xml.Query<std::string>(strXmlRoot + "/pte/source_url", "");
}

template<typename T> PeriodicSystem<T>::~PeriodicSystem() {}

template<typename T>
std::shared_ptr<const PeriodicSystem<T>> PeriodicSystem<T>::GetInstance(const char* pcFile)
{
	std::lock_guard<std::mutex> _guard(s_mutex);

	if(!s_inst)
	{
		//s_inst = std::make_shared<PeriodicSystem<T>>("");
		s_inst = std::shared_ptr<PeriodicSystem<T>>(new PeriodicSystem<T>(
			pcFile ? pcFile : "res/data/elements.xml", ""));
	}

	return s_inst;
}

template<typename T>
const typename PeriodicSystem<T>::elem_type* PeriodicSystem<T>::Find(const std::string& strElem) const
{
	typedef typename decltype(s_vecAtoms)::const_iterator t_iter;

	// elements
	t_iter iter = std::find_if(s_vecAtoms.begin(), s_vecAtoms.end(),
		[&strElem](const elem_type& elem)->bool
		{
			//std::cout << elem.GetAtomIdent() << std::endl;
			return elem.GetAtomIdent() == strElem;
		});
	if(iter != s_vecAtoms.end())
		return &*iter;

	return nullptr;
}



// =============================================================================



template<typename T>
std::shared_ptr<FormfactList<T>> FormfactList<T>::s_inst = nullptr;

template<typename T>
std::mutex FormfactList<T>::s_mutex;


template<typename T>
FormfactList<T>::FormfactList(const std::string& strFile, const std::string& strXmlRoot)
{
	std::string strTabFile = find_resource(strFile);
	tl::log_debug("Loading atomic form factors from file \"", strTabFile, "\".");

	tl::Prop<std::string> xml;
	if(!xml.Load(strTabFile.c_str(), tl::PropType::XML))
		return;

	const std::size_t iNumDat = xml.Query<std::size_t>(strXmlRoot + "/ffacts/num_atoms", 0);
	if(!iNumDat)
	{
		tl::log_err("No data in atomic form factor list.");
		return;
	}

	bool bIonStart = false;
	for(std::size_t iSf=0; iSf<iNumDat; ++iSf)
	{
		elem_type ffact;
		std::string strAtom = "ffacts/atom_" + tl::var_to_str(iSf);

		ffact.strAtom = xml.Query<std::string>(strXmlRoot + "/" +strAtom + "/name", "");
		tl::get_tokens<value_type, std::string, std::vector<value_type>>
			(xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/a", ""), " \t", ffact.a);
		tl::get_tokens<value_type, std::string, std::vector<value_type>>
			(xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/b", ""), " \t", ffact.b);
		ffact.c = xml.Query<value_type>(strXmlRoot + "/" + strAtom + "/c", 0.);

		if(!bIonStart && ffact.strAtom.find_first_of("+-") != std::string::npos)
			bIonStart = true;

		if(!bIonStart)
			s_vecAtoms.push_back(std::move(ffact));
		else
			s_vecIons.push_back(std::move(ffact));
	}

	s_strSrc = xml.Query<std::string>(strXmlRoot + "/ffacts/source", "");
	s_strSrcUrl = xml.Query<std::string>(strXmlRoot + "/ffacts/source_url", "");
}

template<typename T>
FormfactList<T>::~FormfactList()
{}

template<typename T>
std::shared_ptr<const FormfactList<T>> FormfactList<T>::GetInstance(const char* pcFile)
{
	std::lock_guard<std::mutex> _guard(s_mutex);

	if(!s_inst)
	{
		//s_inst = std::make_shared<FormfactList<T>>("");
		s_inst = std::shared_ptr<FormfactList<T>>(new FormfactList<T>(
			pcFile ? pcFile : "res/data/ffacts.xml", ""));
	}

	return s_inst;
}

template<typename T>
const typename FormfactList<T>::elem_type* FormfactList<T>::Find(const std::string& strElem) const
{
	typedef typename decltype(s_vecAtoms)::const_iterator t_iter;

	// atoms
	t_iter iter = std::find_if(s_vecAtoms.begin(), s_vecAtoms.end(),
		[&strElem](const elem_type& elem)->bool
		{
			//std::cout << elem.GetAtomIdent() << std::endl;
			return elem.GetAtomIdent() == strElem;
		});
	if(iter != s_vecAtoms.end())
		return &*iter;

	// ions
	iter = std::find_if(s_vecIons.begin(), s_vecIons.end(),
		[&strElem](const elem_type& elem)->bool
		{
			//std::cout << elem.GetAtomIdent() << std::endl;
			return elem.GetAtomIdent() == strElem;
		});
	if(iter != s_vecIons.end())
		return &*iter;

	return nullptr;
}



// =============================================================================



template<typename T>
std::shared_ptr<MagFormfactList<T>> MagFormfactList<T>::s_inst = nullptr;

template<typename T>
std::mutex MagFormfactList<T>::s_mutex;


template<typename T>
MagFormfactList<T>::MagFormfactList(const std::string& strFile, const std::string& strXmlRoot)
{
	std::string strTabFile = find_resource(strFile);
	tl::log_debug("Loading magnetic form factors from file \"", strTabFile, "\".");

	tl::Prop<std::string, true> xml;
	if(!xml.Load(strTabFile.c_str(), tl::PropType::XML))
		return;

	const std::size_t iNumDat = xml.Query<std::size_t>(strXmlRoot + "/magffacts/num_atoms", 0);
	if(!iNumDat)
	{
		tl::log_err("No data in magnetic form factor list.");
		return;
	}

	for(std::size_t iSf=0; iSf<iNumDat; ++iSf)
	{
		elem_type ffact;
		std::string strAtom = "magffacts/j0/atom_" + tl::var_to_str(iSf);

		std::string strvecA = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/A");
		std::string strveca = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/a");

		tl::get_tokens<value_type>(strvecA, std::string(";"), ffact.A0);
		tl::get_tokens<value_type>(strveca, std::string(";"), ffact.a0);

		ffact.strAtom = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/name", "");

		s_vecAtoms.push_back(std::move(ffact));
	}

	for(std::size_t iSf=0; iSf<iNumDat; ++iSf)
	{
		std::string strAtom = "magffacts/j2/atom_" + tl::var_to_str(iSf);
		std::string strAtomName = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/name", "");

		MagFormfactList<T>::elem_type* pElem =
			const_cast<MagFormfactList<T>::elem_type*>(Find(strAtomName));
		if(!pElem)
		{
			tl::log_err("Mismatch in j0 and j2 form factor tables.");
			continue;
		}

		std::string strvecA = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/A");
		std::string strveca = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/a");

		tl::get_tokens<value_type>(strvecA, std::string(";"), pElem->A2);
		tl::get_tokens<value_type>(strveca, std::string(";"), pElem->a2);
	}

	s_strSrc = xml.Query<std::string>(strXmlRoot + "/magffacts/source", "");
	s_strSrcUrl = xml.Query<std::string>(strXmlRoot + "/magffacts/source_url", "");
}

template<typename T>
MagFormfactList<T>::~MagFormfactList()
{}

template<typename T>
std::shared_ptr<const MagFormfactList<T>> MagFormfactList<T>::GetInstance(const char* pcFile)
{
	std::lock_guard<std::mutex> _guard(s_mutex);

	if(!s_inst)
	{
		//s_inst = std::make_shared<MagFormfactList<T>>("");
		s_inst = std::shared_ptr<MagFormfactList<T>>(new MagFormfactList<T>(
			pcFile ? pcFile : "res/data/magffacts.xml", ""));
	}

	return s_inst;
}

template<typename T>
const typename MagFormfactList<T>::elem_type* MagFormfactList<T>::Find(const std::string& strElem) const
{
	typedef typename decltype(s_vecAtoms)::const_iterator t_iter;

	t_iter iter = std::find_if(s_vecAtoms.begin(), s_vecAtoms.end(),
		[&strElem](const elem_type& elem) -> bool
		{
			return elem.GetAtomIdent() == strElem;
		});
	if(iter != s_vecAtoms.end())
		return &*iter;

	return nullptr;
}



// =============================================================================


template<typename T>
std::shared_ptr<ScatlenList<T>> ScatlenList<T>::s_inst = nullptr;

template<typename T>
std::mutex ScatlenList<T>::s_mutex;


template<typename T>
ScatlenList<T>::ScatlenList(const std::string& strFile, const std::string& strXmlRoot)
{
	std::string strTabFile = find_resource(strFile);
	tl::log_debug("Loading neutron scattering lengths from file \"", strTabFile, "\".");

	tl::Prop<std::string> xml;
	if(!xml.Load(strTabFile.c_str(), tl::PropType::XML))
		return;

	const std::size_t iNumDat = xml.Query<std::size_t>(strXmlRoot + "/scatlens/num_atoms", 0);
	if(!iNumDat)
	{
		tl::log_err("No data in scattering length list.");
		return;
	}

	for(std::size_t iSl=0; iSl<iNumDat; ++iSl)
	{
		ScatlenList<T>::elem_type slen;
		std::string strAtom = "scatlens/atom_" + tl::var_to_str(iSl);

		slen.strAtom = xml.Query<std::string>(strXmlRoot + "/" + strAtom + "/name", "");
		slen.coh = xml.Query<ScatlenList<T>::value_type>(strXmlRoot + "/" + strAtom + "/coh", 0.);
		slen.incoh = xml.Query<ScatlenList<T>::value_type>(strXmlRoot + "/" + strAtom + "/incoh", 0.);

		if(xml.Exists((strAtom + "/xsec_coh").c_str()))
			slen.xsec_coh = xml.Query<ScatlenList<T>::value_type>(strXmlRoot + "/" + strAtom + "/xsec_coh", 0.);
		else
			slen.xsec_coh = (slen.coh*std::conj(slen.coh)).real()*T(4)*tl::get_pi<T>();

		if(xml.Exists((strAtom + "/xsec_incoh").c_str()))
			slen.xsec_incoh = xml.Query<ScatlenList<T>::value_type>(strXmlRoot + "/" + strAtom + "/xsec_incoh", 0.);
		else
			slen.xsec_incoh = (slen.incoh*std::conj(slen.incoh)).real()*T(4)*tl::get_pi<T>();

		if(xml.Exists((strAtom + "/xsec_scat").c_str()))
		{
			slen.xsec_scat = xml.Query<ScatlenList<T>::value_type>(strXmlRoot + "/" + strAtom + "/xsec_scat", 0.);
			//tl::log_debug("Total scattering xsec exists for: ", slen.strAtom, ".");
		}
		else
		{
			slen.xsec_scat = slen.xsec_coh + slen.xsec_incoh;
		}

		slen.xsec_abs = xml.Query<ScatlenList<T>::value_type>(strXmlRoot + "/" + strAtom + "/xsec_absorp", 0.);


		slen.abund = xml.QueryOpt<ScatlenList<T>::real_type>(strXmlRoot + "/" + strAtom + "/abund");
		slen.hl = xml.QueryOpt<ScatlenList<T>::real_type>(strXmlRoot + "/" + strAtom + "/hl");

		if(std::isdigit(slen.strAtom[0]))
			s_vecIsotopes.push_back(std::move(slen));	// pure isotopes
		else
			s_vecElems.push_back(std::move(slen));		// isotope mixtures
	}

	// link pure isotopes to isotope mixtures
	for(const auto& isotope : s_vecIsotopes)
	{
		const std::string& strIso = isotope.GetAtomIdent();
		std::string strElem = tl::remove_chars(strIso, std::string("0123456789"));
		tl::trim(strElem);

		auto iterElem = std::find_if(s_vecElems.begin(), s_vecElems.end(),
			[&strElem](const elem_type& elem) { return (elem.GetAtomIdent() == strElem); });
		if(iterElem == s_vecElems.end())
		{
			tl::log_err("Mixture for isotope \"", strElem, "\" was not found in scattering lengths list.");
			continue;
		}

		iterElem->m_vecIsotopes.push_back(&isotope);
	}

	s_strSrc = xml.Query<std::string>(strXmlRoot + "/scatlens/source", "");
	s_strSrcUrl = xml.Query<std::string>(strXmlRoot + "/scatlens/source_url", "");

#ifndef NDEBUG
	// testing scattering lengths
	for(const auto& elem : s_vecElems)
	{
		const auto& vecIsotopes = elem.GetIsotopes();
		if(!vecIsotopes.size()) continue;

		std::vector<ScatlenList<T>::value_type> vecbcoh, vecbincoh;
		std::vector<ScatlenList<T>::real_type> vecAbund;

		for(const auto* isotope : vecIsotopes)
		{
			const auto& abund = isotope->GetAbundance();
			if(abund)
			{
				vecAbund.push_back(*abund);
				vecbcoh.push_back(isotope->GetCoherent());
				vecbincoh.push_back(isotope->GetIncoherent());
			}
		}

		auto themean = tl::mean_value<decltype(vecAbund), decltype(vecbcoh)>
			(vecAbund, vecbcoh);
		auto thedev = tl::std_dev<decltype(vecAbund), decltype(vecbcoh)>
			(vecAbund, vecbcoh);

		tl::log_debug(elem.GetAtomIdent(), ": mean b: ", themean,
			": stddev b: ", thedev,
			", orig coh: ", elem.GetCoherent(),
			", orig inc: ", elem.GetIncoherent());
	}
#endif
}


template<typename T>
ScatlenList<T>::~ScatlenList()
{}

template<typename T>
std::shared_ptr<const ScatlenList<T>> ScatlenList<T>::GetInstance(const char* pcFile)
{
	std::lock_guard<std::mutex> _guard(s_mutex);

	if(!s_inst)
	{
		//s_inst = std::make_shared<ScatlenList<T>>("");
		s_inst = std::shared_ptr<ScatlenList<T>>(new ScatlenList<T>(
			pcFile ? pcFile : "res/data/scatlens.xml", ""));
	}

	return s_inst;
}

template<typename T>
const typename ScatlenList<T>::elem_type* ScatlenList<T>::Find(const std::string& strElem) const
{
	typedef typename decltype(s_vecElems)::const_iterator t_iter;

	// elements
	t_iter iter = std::find_if(s_vecElems.begin(), s_vecElems.end(),
		[&strElem](const elem_type& elem)->bool
		{
			//std::cout << elem.GetAtomIdent() << std::endl;
			return elem.GetAtomIdent() == strElem;
		});
	if(iter != s_vecElems.end())
		return &*iter;

	// isotopes
	iter = std::find_if(s_vecIsotopes.begin(), s_vecIsotopes.end(),
		[&strElem](const elem_type& elem)->bool
		{
			//std::cout << elem.GetAtomIdent() << std::endl;
			return elem.GetAtomIdent() == strElem;
		});
	if(iter != s_vecIsotopes.end())
		return &*iter;

	return nullptr;
}

}
#endif
