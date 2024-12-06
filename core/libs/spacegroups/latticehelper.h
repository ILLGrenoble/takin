/**
 * Lattice Helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015-2016
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

#ifndef __LATTICE_HELPER_H__
#define __LATTICE_HELPER_H__

#include "tlibs/math/linalg.h"
#include "tlibs/math/geo2.h"
#include "tlibs/phys/lattice.h"

#include "spacegroup.h"
#include "../formfactors/formfact.h"


namespace xtl {

namespace ublas = tl::ublas;


template<class t_real = double>
struct AtomPos
{
	std::string strAtomName;

	// position in fractional units
	ublas::vector<t_real> vecPos;

	// spin
	ublas::vector<t_real> vecSpin;
};


/**
 * additional infos for the atom
 */
template<class t_real = double>
struct AtomPosAux
{
	// position in angstrom units
	ublas::vector<t_real> vecPosAA;

	// (supercell) indices of neighbour atoms
	std::vector<std::vector<std::size_t>> vecIdxNN;

	// coordination polyhedron
	std::vector<std::vector<ublas::vector<t_real>>> vecPolys;
};


/**
 * common lattice calculations needed for several program modules
 */
template<class t_real = double>
struct LatticeCommon
{
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;

	tl::Lattice<t_real> lattice, recip;
	tl::Plane<t_real> planeRLU, planeFrac;

	const SpaceGroup<t_real>* pSpaceGroup = nullptr;
	const std::vector<AtomPos<t_real>>* pvecAtomPos = nullptr;
	std::vector<AtomPosAux<t_real>> vecAtomPosAux;

	t_mat matPlane, matPlaneRLU, matPlane_inv;
	t_mat matPlaneReal, matPlaneReal_inv;

	std::vector<std::string> vecAllNames;
	std::vector<t_vec> vecAllAtoms, vecAllAtomsFrac;
	std::vector<std::size_t> vecAllAtomTypes;
	std::vector<std::complex<t_real>> vecScatlens;
	std::vector<AtomPosAux<t_real>> vecAllAtomPosAux;

	std::vector<t_vec> vecAtomsSC;
	std::vector<std::size_t> vecIdxSC;
	std::vector<std::string> vecNamesSC;

	t_vec dir0RLU, dir1RLU;
	t_mat matA, matB;

	tl::Plane<t_real> plane, planeReal;

	// real unit cell volume in A^3
	t_real dVol = 0;

	// microscopic cross-sections in fm^2
	t_real dsigCoh=0, dsigInc=0, dsigScat=0, dsigAbs=0;
	// macroscopic cross-sections in 1/cm
	t_real dSigScat=0, dSigAbs=0;

	LatticeCommon() = default;
	~LatticeCommon() = default;


	/**
	 * structure factors
	 */
	void CalcStructFacts()
	{
		if(!pSpaceGroup || !pvecAtomPos)
			return;

		std::shared_ptr<const ScatlenList<t_real>> lstsl
			= ScatlenList<t_real>::GetInstance();

		const std::vector<t_mat>* pvecSymTrafos = nullptr;
		if(pSpaceGroup)
			pvecSymTrafos = &pSpaceGroup->GetTrafos();

		if(pvecSymTrafos && pvecSymTrafos->size() && /*g_bHasFormfacts &&*/
			g_bHasScatlens && pvecAtomPos && pvecAtomPos->size())
		{
			std::vector<t_vec> vecAtoms;
			std::vector<std::string> vecNames;
			for(std::size_t iAtom=0; iAtom<pvecAtomPos->size(); ++iAtom)
			{
				vecAtoms.push_back((*pvecAtomPos)[iAtom].vecPos);
				vecNames.push_back((*pvecAtomPos)[iAtom].strAtomName);
			}

			std::tie(vecAllNames, vecAllAtoms, vecAllAtomsFrac, vecAllAtomTypes) =
			tl::generate_all_atoms<t_mat, t_vec, std::vector>
				(*pvecSymTrafos, vecAtoms, &vecNames, matA,
				g_dEps, &matB);

			for(std::size_t iAtom=0; iAtom<vecAllAtoms.size(); ++iAtom)
			{
				tl::set_eps_0(vecAllAtoms[iAtom]);
				tl::set_eps_0(vecAllAtomsFrac[iAtom]);
			}

			for(const std::string& strElem : vecAllNames)
			{
				const typename ScatlenList<t_real>::elem_type* pElem = lstsl->Find(strElem);
				vecScatlens.push_back(pElem ? pElem->GetCoherent() : std::complex<t_real>(0.,0.));
				if(!pElem)
					tl::log_err("Element \"", strElem, "\" not found in scattering length table.",
						" Using b = 0.");
			}
		}
	}


	/**
	 * super cell & environments
	 */
	void CalcSC()
	{
		const unsigned iSC = 3;
		std::vector<std::complex<t_real>> vecDummy;
		std::tie(vecAtomsSC, std::ignore, vecIdxSC) =
			tl::generate_supercell<t_vec, std::vector, t_real>
				(lattice, vecAllAtoms, vecDummy, iSC);
		for(std::size_t iIdxSC : vecIdxSC)
			vecNamesSC.push_back((*pvecAtomPos)[vecAllAtomTypes[iIdxSC]].strAtomName);


		auto fktgetNN = [this](const t_vec& vecCentre) -> AtomPosAux<t_real>
		{
			AtomPosAux<t_real> atomaux;

			atomaux.vecPosAA = vecCentre;
			if(tl::is_nan_or_inf(vecCentre))
				tl::log_err("Invalid atomic position: ", vecCentre, ".");
			tl::set_eps_0(atomaux.vecPosAA, g_dEps);


			// neighbours
			const t_real dEpsShell = 0.01;
			atomaux.vecIdxNN =
				tl::get_neighbours<t_vec, std::vector, t_real>
					(vecAtomsSC, atomaux.vecPosAA, dEpsShell);


			// nearest neighbours
			if(atomaux.vecIdxNN.size() > 1)
			{
				std::vector<t_vec> vecNN;
				for(std::size_t iIdxNN : atomaux.vecIdxNN[1])
				{
					t_vec vecThisAA = vecAtomsSC[iIdxNN] - atomaux.vecPosAA;
					tl::set_eps_0(vecThisAA, g_dEps);
					//const std::string& strThisAtom = vecNamesSC[iIdxNN];

					vecNN.emplace_back(std::move(vecThisAA));
				}


				// calculate coordination polyhedron if enough next neighbours are in list
				if(vecNN.size() >= 4)
					atomaux.vecPolys = tl::get_convexhull<t_vec, std::vector, t_real>(vecNN, g_dEps);
			}

			return atomaux;
		};

		// neighbours of all atoms in primitive UC
		for(const AtomPos<t_real>& atomPos : *pvecAtomPos)
		{
			t_vec vecPosAA = tl::mult<t_mat, t_vec>(matA, atomPos.vecPos);
			vecAtomPosAux.emplace_back(std::move(fktgetNN(vecPosAA)));
		}

		// neighbours of all atoms in conventional UC
		for(const t_vec& vecCentre : vecAllAtoms)
			vecAllAtomPosAux.emplace_back(std::move(fktgetNN(vecCentre)));
	}


	/**
	 * calculate scattering cross sections, etc.
	 */
	void CalcCrysInfos()
	{
		dsigCoh = dsigInc = dsigScat = dsigAbs = 0;
		dSigScat = dSigAbs = 0;

		dVol = lattice.GetVol();

		std::shared_ptr<const ScatlenList<t_real>> lstsl
			= ScatlenList<t_real>::GetInstance();

		for(const std::string& strElem : vecAllNames)
		{
   			const typename ScatlenList<t_real>::elem_type* pElem = lstsl->Find(strElem);
			if(!pElem)
			{
				tl::log_err("Element \"", strElem, "\" not found in scattering length table. Ignoring in cross-sections!");
				continue;
			}

			// microscopic cross-sections
			dsigCoh += pElem->GetXSecCoherent().real();
			dsigInc += pElem->GetXSecIncoherent().real();
			dsigScat += pElem->GetXSecScatter().real();
			dsigAbs += pElem->GetXSecAbsorption().real();

			// macroscopic cross-sections
			dSigScat += tl::macro_xsect(pElem->GetXSecScatter().real()*tl::get_one_femtometer<t_real>()*tl::get_one_femtometer<t_real>(),
				1, dVol*tl::get_one_angstrom<t_real>()*tl::get_one_angstrom<t_real>()*tl::get_one_angstrom<t_real>()) * tl::get_one_centimeter<t_real>();

			dSigAbs += tl::macro_xsect(pElem->GetXSecAbsorption().real()*tl::get_one_femtometer<t_real>()*tl::get_one_femtometer<t_real>(),
				1, dVol*tl::get_one_angstrom<t_real>()*tl::get_one_angstrom<t_real>()*tl::get_one_angstrom<t_real>()) * tl::get_one_centimeter<t_real>();
		}
	}


	/**
	 * calculate everything
	 */
	bool Calc(const tl::Lattice<t_real>& lat, const tl::Lattice<t_real>& rec,
		const tl::Plane<t_real>& plRLU, const tl::Plane<t_real>& plFrac,
		const SpaceGroup<t_real> *pSG, const std::vector<AtomPos<t_real>>* pvecAt)
	{
		lattice = lat;
		recip = rec;

		planeRLU = plRLU;	// for recip lattice
		planeFrac = plFrac;	// for real lattice

		pSpaceGroup = pSG;
		pvecAtomPos = pvecAt;

		{	// reciprocal lattice
			t_vec vecX0 = ublas::zero_vector<t_real>(3);
			t_vec vecPlaneX = planeRLU.GetDir0()[0] * recip.GetVec(0) +
				planeRLU.GetDir0()[1] * recip.GetVec(1) +
				planeRLU.GetDir0()[2] * recip.GetVec(2);
			t_vec vecPlaneY = planeRLU.GetDir1()[0] * recip.GetVec(0) +
				planeRLU.GetDir1()[1] * recip.GetVec(1) +
				planeRLU.GetDir1()[2] * recip.GetVec(2);

			std::vector<t_vec> vecOrthRlu =
				tl::gram_schmidt<t_vec>(
					{ planeRLU.GetDir0(), planeRLU.GetDir1(), planeRLU.GetNorm() },
					true);
			matPlaneRLU = tl::column_matrix(vecOrthRlu);

			plane = tl::Plane<t_real>(vecX0, vecPlaneX, vecPlaneY);
			std::vector<t_vec> vecOrth =
				tl::gram_schmidt<t_vec>(
					{ plane.GetDir0(), plane.GetDir1(), plane.GetNorm() },
					true);
			matPlane = tl::column_matrix(vecOrth);
			t_real dDir0Len = ublas::norm_2(vecOrth[0]), dDir1Len = ublas::norm_2(vecOrth[1]);

			if(tl::float_equal<t_real>(dDir0Len, 0.) || tl::float_equal<t_real>(dDir1Len, 0.)
				|| tl::is_nan_or_inf<t_real>(dDir0Len) || tl::is_nan_or_inf<t_real>(dDir1Len))
			{
				tl::log_err("Invalid scattering plane.");
				return false;
			}

			bool bInv = tl::inverse(matPlane, matPlane_inv);
			if(!bInv)
			{
				tl::log_err("Cannot invert scattering plane metric.");
				return false;
			}
		}

		{	// real lattice
			t_vec vecX0 = ublas::zero_vector<t_real>(3);
			t_vec vecPlaneX = planeFrac.GetDir0()[0] * lattice.GetVec(0) +
				planeFrac.GetDir0()[1] * lattice.GetVec(1) +
				planeFrac.GetDir0()[2] * lattice.GetVec(2);
			t_vec vecPlaneY = planeFrac.GetDir1()[0] * lattice.GetVec(0) +
				planeFrac.GetDir1()[1] * lattice.GetVec(1) +
				planeFrac.GetDir1()[2] * lattice.GetVec(2);
			planeReal = tl::Plane<t_real>(vecX0, vecPlaneX, vecPlaneY);

			std::vector<t_vec> vecOrth =
				tl::gram_schmidt<t_vec>(
					{planeReal.GetDir0(), planeReal.GetDir1(), planeReal.GetNorm()},
					true);
			matPlaneReal = tl::column_matrix(vecOrth);
			t_real dDir0Len = ublas::norm_2(vecOrth[0]), dDir1Len = ublas::norm_2(vecOrth[1]);

			if(tl::float_equal<t_real>(dDir0Len, 0.) || tl::float_equal<t_real>(dDir1Len, 0.)
				|| tl::is_nan_or_inf<t_real>(dDir0Len) || tl::is_nan_or_inf<t_real>(dDir1Len))
			{
				tl::log_err("Invalid plane in real lattice.");
				return false;
			}

			bool bInv = tl::inverse(matPlaneReal, matPlaneReal_inv);
			if(!bInv)
			{
				tl::log_err("Cannot invert plane metric in real lattice.");
				return false;
			}
		}


		matA = lattice.GetBaseMatrixCov();
		matB = recip.GetBaseMatrixCov();

		dir0RLU = planeRLU.GetDir0();
		dir1RLU = planeRLU.GetDir1();


		CalcStructFacts();
		CalcSC();
		CalcCrysInfos();

		return true;
	}


	bool CanCalcStructFact() const { return vecScatlens.size() != 0; }


	std::tuple<std::complex<t_real>, t_real, t_real> GetStructFact(const t_vec& vecPeak) const
	{
		std::complex<t_real> cF =
			tl::structfact<t_real, std::complex<t_real>, t_vec, std::vector>
				(vecAllAtoms, vecPeak, vecScatlens);
		t_real dFsq = (std::conj(cF)*cF).real();
		t_real dF = std::sqrt(dFsq);

		return std::make_tuple(cF, dF, dFsq);
	}
};

}
#endif
