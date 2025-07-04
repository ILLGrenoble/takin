/**
 * Convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
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

#ifndef __CONVOFIT_SCAN_H__
#define __CONVOFIT_SCAN_H__

#include <vector>
#include <string>
#include <boost/optional.hpp>

#include "tlibs/math/math.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/file/loadinstr.h"

#include "../res/defs.h"
using t_real_sc = t_real_reso;


struct ScanPoint
{
	t_real_sc h, k, l;
	tl::t_wavenumber_si<t_real_sc> ki, kf;
	tl::t_energy_si<t_real_sc> Ei, Ef, E;
};

struct Sample
{
	t_real_sc a, b, c;
	t_real_sc alpha, beta, gamma;
};

struct Plane
{
	t_real_sc vec1[3];
	t_real_sc vec2[3];
};


struct Filter
{
	boost::optional<t_real_sc> dLower;
	boost::optional<t_real_sc> dUpper;
	boost::optional<std::pair<std::string, std::string>> colEquals;
};


struct Scan
{
	Sample sample;
	Plane plane;
	bool bKiFixed = false;
	t_real_sc dKFix = 2.662;

	std::string strTempCol = "TT";
	t_real_sc dTemp = 100., dTempErr = 0.;

	std::string strFieldCol = "";
	t_real_sc dField = 0., dFieldErr = 0.;

	std::string strCntCol = "";
	std::string strMonCol = "";
	std::string strCntErrCol = "";
	std::string strMonErrCol = "";
	std::vector<ScanPoint> vecPoints;

	std::vector<t_real_sc> vechklE[4];
	std::vector<t_real_sc> vecX;
	std::vector<t_real_sc> vecCts, vecMon;
	std::vector<t_real_sc> vecCtsErr, vecMonErr;

	t_real_sc vecScanOrigin[4];      // scan origin
	t_real_sc vecScanDir[4];         // actual scan axis
	t_real_sc vecMainScanDir[4];     // principal scan axis
	unsigned int m_iScIdx = 0;       // index of principal scan axis


	ScanPoint InterpPoint(std::size_t i, std::size_t N) const
	{
		const ScanPoint& ptBegin = *vecPoints.cbegin();
		const ScanPoint& ptEnd = *vecPoints.crbegin();

		ScanPoint pt;

		pt.h = tl::lerp(ptBegin.h, ptEnd.h, t_real_sc(i)/t_real_sc(N-1));
		pt.k = tl::lerp(ptBegin.k, ptEnd.k, t_real_sc(i)/t_real_sc(N-1));
		pt.l = tl::lerp(ptBegin.l, ptEnd.l, t_real_sc(i)/t_real_sc(N-1));
		pt.E = tl::lerp(ptBegin.E, ptEnd.E, t_real_sc(i)/t_real_sc(N-1));
		pt.Ei = tl::lerp(ptBegin.Ei, ptEnd.Ei, t_real_sc(i)/t_real_sc(N-1));
		pt.Ef = tl::lerp(ptBegin.Ef, ptEnd.Ef, t_real_sc(i)/t_real_sc(N-1));
		bool bImag=0;
		pt.ki = tl::E2k(pt.Ei, bImag);
		pt.kf = tl::E2k(pt.Ef, bImag);

		return pt;
	}
};


extern bool load_file(const std::vector<std::string>& vecFiles, Scan& scan,
	bool bNormToMon = true, const Filter& filter = Filter(),
	bool bFlipCoords = false, bool bAllowScanMerging = false,
	bool bUseFirstAndLastPoints = false,
	unsigned iScanAxis = 0, bool bVerbose = true);

extern bool load_file(const char* pcFile, Scan& scan,
	bool bNormToMon = true, const Filter& filter = Filter(),
	bool bFlipCoords = false, bool bAllowScanMerging = false,
	bool bUseFirstAndLastPoints = false);

extern bool save_file(const char* pcFile, const Scan& sc, const char* comment = nullptr);


#endif
