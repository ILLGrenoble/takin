/**
 * atoms and structural calculations
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015-2016
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

#ifndef __TLIBS_ATOMS_H__
#define __TLIBS_ATOMS_H__


#include "../math/linalg.h"
#include "../math/linalg_ops.h"
#include "../math/rt.h"
#include "lattice.h"
#include <tuple>


namespace tl{

/**
 * Maps fractional atom position back to [tMin, tMax]
 */
template<class t_vec>
void restrict_to_uc(t_vec& vec,
	typename t_vec::value_type tMin=0, typename t_vec::value_type tMax=1)
{
	using T = typename t_vec::value_type;
	const T dSize = std::abs(tMax-tMin);

	for(std::size_t i=0; i<vec.size(); ++i)
	{
		vec[i] = std::fmod(vec[i], dSize);

		while(vec[i] < tMin) vec[i] += dSize;
		while(vec[i] >= tMax) vec[i] -= dSize;
	}
}


/**
 * Maps atom position back to units cell using an rt index tree (in angs coordinates!)
 */
template<class t_vec, class t_mat>
void restrict_to_uc(const Rt<typename t_vec::value_type, 3, 8, 0>& rtPeaks,
	const t_mat& matA, const t_mat& matB, t_vec& vecAtom,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
{
	using T = typename t_vec::value_type;

	t_vec vecAtomAA = mult<t_mat, t_vec>(matA, vecAtom);
	std::vector<T> nodeNearestAA =
		rtPeaks.GetNearestNode(std::vector<T>{{ vecAtomAA[0], vecAtomAA[1], vecAtomAA[2] }});

	// if atom is not in the (000) zone, subtract a lattice vector to shift it there
	if(!float_equal<T>(nodeNearestAA[0], 0, eps) ||
		!float_equal<T>(nodeNearestAA[1], 0, eps) ||
		!float_equal<T>(nodeNearestAA[2], 0, eps))
	{
		t_mat matAinv = transpose(matB) / (tl::get_pi<T>()*T(2));

		const t_vec vecNearestAA = tl::make_vec<t_vec>({nodeNearestAA[0], nodeNearestAA[1], nodeNearestAA[2]});
		const t_vec vecNearestFrac = mult<t_mat, t_vec>(matAinv, vecNearestAA);
		//log_debug("nearest lattice vector: ", vecNearestFrac, ", ", vecNearestAA);


		// don't shift atoms on the border of the next cell
		const T dAtomTo000 = veclen<t_vec>(vecAtomAA);
		const T dAtomToNearest = veclen<t_vec>(vecAtomAA-vecNearestAA);

		if(dAtomToNearest < dAtomTo000-eps)
			vecAtom -= vecNearestFrac;
	}
}


/**
 * test if position (in fractional units) and equivalent
 * positions are already occupied
 */
template<class t_vec, template<class...> class t_cont>
bool position_occupied_frac(const t_cont<t_vec>& vecvecPos, const t_vec& vecPos,
	typename t_vec::value_type eps)
{
	static const t_cont<t_vec> vecEquiv
	{{
		make_vec<t_vec>( {0., 0., 0.} ), make_vec<t_vec>( {0., 0., 1.} ),
		make_vec<t_vec>( {0., 0., -1.} ), make_vec<t_vec>( {0., 1., 0.} ),
		make_vec<t_vec>( {0., -1., 0.} ), make_vec<t_vec>( {1., 0., 0.} ),
		make_vec<t_vec>( {-1., 0., 0.} ),

		make_vec<t_vec>( {0., 1., 1.} ), make_vec<t_vec>( {0., 1., -1.} ),
		make_vec<t_vec>( {0., -1., 1.} ), make_vec<t_vec>( {0., -1., -1.} ),
		make_vec<t_vec>( {1., 1., 0.} ), make_vec<t_vec>( {1., -1., 0.} ),
		make_vec<t_vec>( {-1., 1., 0.} ), make_vec<t_vec>( {-1., -1., 0.} ),
		make_vec<t_vec>( {1., 0., 1.} ), make_vec<t_vec>( {1., 0., -1.} ),
		make_vec<t_vec>( {-1., 0., 1.} ), make_vec<t_vec>( {-1., 0., -1.} ),

		make_vec<t_vec>( {1., 1., 1.} ), make_vec<t_vec>( {1., 1., -1.} ),
		make_vec<t_vec>( {1., -1., 1.} ), make_vec<t_vec>( {-1., 1., 1.} ),
		make_vec<t_vec>( {1., -1., -1.} ), make_vec<t_vec>( {-1., 1., -1.} ),
		make_vec<t_vec>( {-1., -1., 1.} ), make_vec<t_vec>( {-1., -1., -1.} ),
	}};

	for(const t_vec& vecOld : vecvecPos)
		for(const t_vec& vecShift : vecEquiv)
			if(vec_equal<t_vec>(vecPos+vecShift, vecOld, eps))
				return true;

	return false;
}


/**
 * Generates atom positions using trafo matrices
 */
template<class t_mat, class t_vec, template<class...> class t_cont>
t_cont<t_vec> generate_atoms(const t_cont<t_mat>& trafos, const t_vec& _vecAtom,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon(),
	const t_mat* pmatA = nullptr, const t_mat* pmatB = nullptr)
{
	using T = typename t_vec::value_type;


	// calculate lattice in angs units to calculate unit cell infos
	constexpr const T max_peak = 8;
	Rt<T, 3, 8, 0> rtPeaks;
	if(pmatA && pmatB)
	{
		// generate spatial index tree using angs coords
		for(T h=-max_peak; h<=max_peak; h+=1)
			for(T k=-max_peak; k<=max_peak; k+=1)
				for(T l=-max_peak; l<=max_peak; l+=1)
				{
					t_vec vecPeakFrac = make_vec<t_vec>({h, k, l});
					t_vec vecPeakAA = mult<t_mat, t_vec>(*pmatA, vecPeakFrac);
					//log_debug("lattice vector: ", vecPeakFrac, ", ", vecPeakAA);

					rtPeaks.InsertPoint(
						std::vector<T>{{vecPeakAA[0], vecPeakAA[1], vecPeakAA[2]}}, 0);
				}
	}


	// trafos are in homogeneous coordinates
	t_vec vecAtom = _vecAtom;
	vecAtom.resize(4,1); vecAtom[3] = T(1);

	t_cont<t_vec> vecvecRes;

	for(const t_mat& mat : trafos)
	{
		t_vec vecRes = mult<t_mat, t_vec>(mat, vecAtom);
		// back to non-homogeneous coordinates
		vecRes.resize(3,1);

		// restrict atoms to fractional coordinates of [-0.5, 0.5]
		restrict_to_uc<t_vec>(vecRes, T(-0.5), T(0.5));
		// if crystal matrix is available, restrict atoms to unit cell
		if(pmatA && pmatB)
			restrict_to_uc<t_vec, t_mat>(rtPeaks, *pmatA, *pmatB, vecRes, eps);

		// skip already occupied positions
		if(!position_occupied_frac<t_vec, t_cont>(vecvecRes, vecRes, eps))
			vecvecRes.push_back(std::move(vecRes));
	}

	return vecvecRes;
}


/**
 * Generates atom positions using trafo matrices for all atoms in unit cell
 * @return tuple of (names, positions, positions in rlu, atom types)
 */
template<class t_mat, class t_vec, template<class...> class t_cont,
	class t_str=std::string, class t_real = typename t_mat::value_type>
std::tuple<t_cont<t_str>, t_cont<t_vec>, t_cont<t_vec>, t_cont<std::size_t>>
generate_all_atoms(const t_cont<t_mat>& trafos,
	const t_cont<t_vec>& vecAtoms, const t_cont<t_str>* pvecNames,
	const t_mat& matA, t_real eps = std::numeric_limits<t_real>::epsilon(),
	const t_mat* pmatB = nullptr)
{
	t_cont<t_vec> vecAllAtoms, vecAllAtomsFrac;
	t_cont<t_str> vecAllNames;
	t_cont<std::size_t> vecAllAtomTypes;

	for(std::size_t iAtom=0; iAtom<vecAtoms.size(); ++iAtom)
	{
		const t_vec& vecAtom = vecAtoms[iAtom];
		t_str strNone;

		const t_str& strElem = pvecNames ? (*pvecNames)[iAtom] : strNone;

		t_cont<t_vec> vecOtherAtoms = vecAtoms;
		vecOtherAtoms.erase(vecOtherAtoms.begin() + iAtom);

		t_cont<t_vec> vecSymPos =
			tl::generate_atoms<t_mat, t_vec, t_cont>
				(trafos, vecAtom, eps, &matA, pmatB);


		//std::size_t iGeneratedAtoms = 0;
		for(t_vec vecThisAtom : vecSymPos)
		{
			// is the atom position in the unit cell still free?
			if(std::find_if(vecAllAtomsFrac.begin(), vecAllAtomsFrac.end(),
				[&vecThisAtom, eps](const t_vec& _v) -> bool
				{ return tl::vec_equal(_v, vecThisAtom, eps); }) == vecAllAtomsFrac.end()
				&& // and is it not at a given initial atom position?
				std::find_if(vecOtherAtoms.begin(), vecOtherAtoms.end(),
				[&vecThisAtom, eps](const t_vec& _v) -> bool
				{ return tl::vec_equal(_v, vecThisAtom, eps); }) == vecOtherAtoms.end())
			{
				vecAllAtomsFrac.push_back(vecThisAtom);

				// converts from fractional coordinates
				vecThisAtom = mult<t_mat, t_vec>(matA, vecThisAtom);
				vecAllAtoms.push_back(std::move(vecThisAtom));
				vecAllNames.push_back(strElem);
				vecAllAtomTypes.push_back(iAtom);

				//++iGeneratedAtoms;
			}
			else
			{
				tl::log_warn("Position ", vecThisAtom, " is already occupied,",
					" skipping current ", strElem, " atom.");
			}
		}
		//tl::log_info("Unit cell has ", iGeneratedAtoms, " ", strElem, " atom(s).");
	}

	return std::make_tuple(vecAllNames,
		vecAllAtoms, vecAllAtomsFrac, vecAllAtomTypes);
}


/**
 * Generates supercell
 * @return tuple of positions, (user-defined) factors and indices
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class t_real = typename t_vec::value_type>
std::tuple<t_cont<t_vec>, t_cont<std::complex<t_real>>, t_cont<std::size_t>>
generate_supercell(const Lattice<t_real>& latt,
	const t_cont<t_vec>& vecAtomsUC,
	const t_cont<std::complex<t_real>>& vecFactsUC,
	std::ptrdiff_t N)
{
	using t_cplx = std::complex<t_real>;

	t_cont<t_vec> vecAllAtoms;
	t_cont<t_cplx> vecAllFacts;
	t_cont<std::size_t> vecAllIdx;

	for(std::ptrdiff_t h=-N+1; h<N; ++h)
	for(std::ptrdiff_t k=-N+1; k<N; ++k)
	for(std::ptrdiff_t l=-N+1; l<N; ++l)
	{
		t_vec vecPos = latt.GetPos(h,k,l);

		for(std::size_t iAtom=0; iAtom<vecAtomsUC.size(); ++iAtom)
		{
			vecAllIdx.push_back(iAtom);
			const t_vec& vecAtom = vecAtomsUC[iAtom];
			t_cplx cFact;
			if(vecFactsUC.size() == vecAtomsUC.size())
				cFact = vecFactsUC[iAtom];
			else if(vecFactsUC.size() == 1)		// use the same for all atoms
				cFact = vecFactsUC[0];
			else if(vecFactsUC.size() == 0)		// use 0 for all atoms
				cFact = t_real(0);

			vecAllAtoms.push_back(vecPos + vecAtom);
			if(vecFactsUC.size() != 0)
				vecAllFacts.push_back(cFact);
		}
	}

	return std::make_tuple(vecAllAtoms, vecAllFacts, vecAllIdx);
}
// ----------------------------------------------------------------------------


/**
 * calculates atomic form factors
 * @param G Length of lattice vector
 * @param vecA "a" coefficients
 * @param vecB "b" coefficients
 * @param c "c" coefficient
 * @return form factor
 *
 * @see Waasmaier and Kirfel, Acta Cryst. A51, 416-431 (1995), doi: https://doi.org/10.1107/S0108767394013292
 */
template<class T=double, template<class...> class t_cont>
T formfact(T G, const t_cont<T>& vecA, const t_cont<T>& vecB, T c=0)
{
	assert(vecA.size() == vecB.size() /*|| vecA.size() == vecB.size()+1*/);

	T ff = T(0);
	T s = G / (T(4)*get_pi<T>());

	// 'c' coefficient directly stored in 'vecA'
	//if(vecA.size() == vecB.size()+1)
	//	c = vecA.rbegin();

	typename t_cont<T>::const_iterator iterA = vecA.begin();
	typename t_cont<T>::const_iterator iterB = vecB.begin();

	for(; iterA!=vecA.end() && iterB!=vecB.end(); ++iterA, ++iterB)
		ff += (*iterA)*std::exp(-(*iterB)*s*s);
	ff += c;

	return ff;
}


/**
 * calculates the structure factor F
 * @param lstAtoms List of atom positions
 * @param vecG Lattice vector
 * @param lstf G-dependent Atomic form factors (x-rays) or coherent scattering length (neutrons)
 * @param pF0 optional total form factor.
 * @param dVuc optionally normalise by the unit cell volume
 * @return structure factor
 *
 * @see (Shirane 2002), p. 25, equ. 2.26
 */
template<typename T = double, typename t_ff = std::complex<T>,
	class t_vec = ublas::vector<T>,
	template<class ...> class t_cont = std::initializer_list>
std::complex<T> structfact(const t_cont<t_vec>& lstAtoms, const t_vec& vecG,
	const t_cont<t_ff>& lstf = t_cont<t_ff>(),
	t_ff *pF0 = nullptr, T dVuc = T(-1))
{
	constexpr std::complex<T> i(0., 1.);
	std::complex<T> F(0., 0.);

	using t_iter_atoms = typename t_cont<t_vec>::const_iterator;
	using t_iter_ffact = typename t_cont<t_ff>::const_iterator;

	t_iter_atoms iterAtom = lstAtoms.begin();
	t_iter_ffact iterFFact = lstf.begin();

	if(pF0) *pF0 = t_ff(0);

	for(; iterAtom != lstAtoms.end(); ++iterAtom)
	{
		// only use form factors or scattering lengths when available
		t_ff tFF = T(1);
		if(iterFFact != lstf.end())
			tFF = *iterFFact;

		F += tFF * std::exp(i * (mult<t_vec, t_vec>(vecG, *iterAtom)));
		if(pF0) *pF0 += tFF;

		// if there is only one form factor in the list, use it for all positions
		if(iterFFact!=lstf.end() && std::next(iterFFact)!=lstf.end())
			++iterFFact;
	}

	// if unit cell volume is given, normalise by it
	if(dVuc > T(0)) F /= dVuc;
	return F;
}


/**
 * Lorentz factor
 * @param twotheta Scattering angle in rad
 * @see (Shirane 2002), pp. 170-172
 * @see e.g.: https://dictionary.iucr.org/Lorentz%E2%80%93polarization_correction
 */
template<typename T=double>
T lorentz_factor(T twotheta)
{
	T theta = 0.5*twotheta;
	return T(0.25) / (std::sin(theta)*std::sin(theta) * std::cos(theta));
}


/**
 * Lorentz polarisation factor (only for x-rays)
 * @param twotheta Scattering angle in rad
 * @see (Shirane 2002), pp. 170-172
 * @see e.g.: https://dictionary.iucr.org/Lorentz%E2%80%93polarization_correction
 */
template<typename T=double>
T lorentz_pol_factor(T twotheta)
{
	return T(0.5) + T(0.5)*std::cos(twotheta)*std::cos(twotheta);
}

}
#endif
