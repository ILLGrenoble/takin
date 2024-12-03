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

#include "scan.h"
#include "tlibs/log/log.h"
#include "tlibs/math/stat.h"

#include <fstream>

#ifndef USE_BOOST_REX
	#include <regex>
	namespace rex = ::std;
#else
	#include <boost/tr1/regex.hpp>
	namespace rex = ::boost;
#endif



/**
 * saving a scan file
 */
bool save_file(const char* pcFile, const Scan& sc)
{
	std::ofstream ofstr(pcFile);
	if(!ofstr)
		return false;

	ofstr.precision(16);

	ofstr << "# scan_origin = "
		<< sc.vecScanOrigin[0] << " "
		<< sc.vecScanOrigin[1] << " "
		<< sc.vecScanOrigin[2] << " "
		<< sc.vecScanOrigin[3] << "\n";
	ofstr << "# scan_dir = "
		<< sc.vecScanDir[0] << " "
		<< sc.vecScanDir[1] << " "
		<< sc.vecScanDir[2] << " "
		<< sc.vecScanDir[3] << "\n";
	ofstr << "# T = " << sc.dTemp << " +- " << sc.dTempErr << "\n";
	ofstr << "# B = " << sc.dField << " +- " << sc.dFieldErr << "\n";

	ofstr << "#\n";

	ofstr << std::left << std::setw(21) << "# x"
		<< std::left << std::setw(21) << "counter"
		<< std::left << std::setw(21) << "counter errors"
		<< std::left << std::setw(21) << "monitor"
		<< std::left << std::setw(21) << "monitor errors" << "\n";

	const std::size_t iNum = sc.vecX.size();
	for(std::size_t i=0; i<iNum; ++i)
	{
		ofstr << std::left << std::setw(20) << sc.vecX[i] << " "
			<< std::left << std::setw(20) << sc.vecCts[i] << " "
			<< std::left << std::setw(20) << sc.vecCtsErr[i] << " "
			<< std::left << std::setw(20) << sc.vecMon[i] << " "
			<< std::left << std::setw(20) << sc.vecMonErr[i] << "\n";
	}

	return true;
}



/**
 * loading multiple scan files
 */
bool load_file(const std::vector<std::string>& vecFiles, Scan& scan, bool bNormToMon,
	const Filter& filter, bool bFlipCoords, bool bAllowScanMerging,
	bool bUseFirstAndLastPoints, unsigned iScanAxis, bool bVerbose)
{
	if(!vecFiles.size())
		return 0;
	tl::log_info("Loading \"", vecFiles[0], "\".");

	std::unique_ptr<tl::FileInstrBase<t_real_sc>>
		pInstr(tl::FileInstrBase<t_real_sc>::LoadInstr(vecFiles[0].c_str()));
	if(!pInstr)
	{
		tl::log_err("Cannot load \"", vecFiles[0], "\".");
		return false;
	}

	for(std::size_t iFile=1; iFile<vecFiles.size(); ++iFile)
	{
		tl::log_info("Loading \"", vecFiles[iFile], "\" for merging.");
		std::unique_ptr<tl::FileInstrBase<t_real_sc>>
			pInstrM(tl::FileInstrBase<t_real_sc>::LoadInstr(vecFiles[iFile].c_str()));
		if(!pInstrM)
		{
			tl::log_err("Cannot load \"", vecFiles[iFile], "\".");
			continue;
		}

		pInstr->MergeWith(pInstrM.get(), bAllowScanMerging);
	}


	// defaults
	std::string strCountVar = pInstr->GetCountVar();
	std::string strMonVar = pInstr->GetMonVar();
	std::string strCountErr = pInstr->GetCountErr();
	std::string strMonErr = pInstr->GetMonErr();

	// overrides
	if(scan.strCntCol != "")
	{
		if(pInstr->HasCol(scan.strCntCol))
			strCountVar = scan.strCntCol;
		else
			tl::log_err("Invalid counter column \"", scan.strCntCol, "\", using default: \"", strCountVar, "\".");
	}

	if(scan.strMonCol != "")
	{
		if(pInstr->HasCol(scan.strMonCol))
			strMonVar = scan.strMonCol;
		else
			tl::log_err("Invalid monitor column \"", scan.strMonCol, "\", using default: \"", strMonVar, "\".");
	}

	if(scan.strCntErrCol != "")
	{
		if(pInstr->HasCol(scan.strCntErrCol))
			strCountErr = scan.strCntErrCol;
		else
			tl::log_err("Invalid counter error column \"", scan.strCntErrCol, "\", using default: \"", strCountErr, "\".");
	}

	if(scan.strMonErrCol != "")
	{
		if(pInstr->HasCol(scan.strMonErrCol))
			strMonErr = scan.strMonErrCol;
		else
			tl::log_err("Invalid monitor error column \"", scan.strMonErrCol, "\", using default: \"", strMonErr, "\".");
	}

	scan.vecCts = pInstr->GetCol(strCountVar);
	scan.vecMon = pInstr->GetCol(strMonVar);

	tl::log_info("Counter column: ", strCountVar, ", monitor column: ", strMonVar, ".");

	// error
	std::function<t_real_sc(t_real_sc)> funcErr = [](t_real_sc d) -> t_real_sc
	{
		if(tl::float_equal<t_real_sc>(d, 0.))
			return t_real_sc(1);
		return std::sqrt(std::abs(d));
	};

	if(strCountErr != "")
	{
		// use the given counter error column
		scan.vecCtsErr = pInstr->GetCol(strCountErr);
		tl::log_info("Counter error column: ", strCountErr, ".");
	}
	else
	{
		// calculate the error from the counter
		scan.vecCtsErr = tl::apply_fkt(scan.vecCts, funcErr);
	}

	if(strMonErr != "")
	{
		// use the given monitor error column
		scan.vecMonErr = pInstr->GetCol(strMonErr);
		tl::log_info("Monitor error column: ", strMonErr, ".");
	}
	else
	{
		// calculate the error from the monitor
		scan.vecMonErr = tl::apply_fkt(scan.vecMon, funcErr);
	}


	if(bNormToMon)
	{
		// normalise couter to monitor
		for(std::size_t iPos = 0; iPos < scan.vecCts.size(); ++iPos)
		{
			t_real_sc y = scan.vecCts[iPos];
			t_real_sc dy = scan.vecCtsErr[iPos];
			t_real_sc m = scan.vecMon[iPos];
			t_real_sc dm  = scan.vecMonErr[iPos];

			scan.vecCts[iPos] /= m;
			scan.vecCtsErr[iPos] = std::sqrt(dy/m * dy/m  +  y*dm/(m*m) * y*dm/(m*m));
		}
	}

	const std::array<t_real_sc, 3> latt = pInstr->GetSampleLattice();
	const std::array<t_real_sc, 3> ang = pInstr->GetSampleAngles();

	scan.sample.a = latt[0]; scan.sample.b = latt[1]; scan.sample.c = latt[2];
	scan.sample.alpha = ang[0]; scan.sample.beta = ang[1]; scan.sample.gamma = ang[2];

	tl::log_info("Sample lattice: ", scan.sample.a, " ", scan.sample.b, " ", scan.sample.c, ".");
	tl::log_info("Sample angles: ", tl::r2d(scan.sample.alpha), " ", tl::r2d(scan.sample.beta), " ", tl::r2d(scan.sample.gamma), ".");


	const std::array<t_real_sc, 3> vec1 = pInstr->GetScatterPlane0();
	const std::array<t_real_sc, 3> vec2 = pInstr->GetScatterPlane1();
	t_real_sc dFlip = bFlipCoords ? t_real_sc(-1) : t_real_sc(1);
	scan.plane.vec1[0] = vec1[0]; scan.plane.vec1[1] = vec1[1]; scan.plane.vec1[2] = vec1[2];
	scan.plane.vec2[0] = dFlip*vec2[0]; scan.plane.vec2[1] = dFlip*vec2[1]; scan.plane.vec2[2] = dFlip*vec2[2];

	tl::log_info("Scattering plane: [", vec1[0], vec1[1], vec1[2], "], "
		"[", vec2[0], vec2[1], vec2[2], "].");
	if(bFlipCoords)
		tl::log_info("Flipped RHS <-> LHS coordinate system.");


	scan.bKiFixed = pInstr->IsKiFixed();
	scan.dKFix = pInstr->GetKFix();
	if(scan.bKiFixed)
		tl::log_info("ki = ", scan.dKFix, " = const.");
	else
		tl::log_info("kf = ", scan.dKFix, " = const.");


	if(scan.strTempCol != "")
	{
		const tl::FileInstrBase<t_real_sc>::t_vecVals& vecTemp = pInstr->GetCol(scan.strTempCol);

		if(vecTemp.size() == 0)
		{
			tl::log_warn("Sample temperature column \"", scan.strTempCol, "\" not found.");
		}
		else
		{
			scan.dTemp = tl::mean_value(vecTemp);
			scan.dTempErr = tl::std_dev(vecTemp);
			tl::log_info("Sample temperature: ", scan.dTemp, " +- ", scan.dTempErr, ".");
		}
	}

	if(scan.strFieldCol != "")
	{
		const tl::FileInstrBase<t_real_sc>::t_vecVals& vecField = pInstr->GetCol(scan.strFieldCol);

		if(vecField.size() == 0)
		{
			tl::log_warn("Sample field column \"", scan.strFieldCol, "\" not found.");
		}
		else
		{
			scan.dField = tl::mean_value(vecField);
			scan.dFieldErr = tl::std_dev(vecField);
			tl::log_info("Sample field: ", scan.dField, " +- ", scan.dFieldErr, ".");
		}
	}


	// component-wise minima and maxima
	ScanPoint ptMin, ptMax;
	ptMin.h = ptMin.k = ptMin.l = std::numeric_limits<t_real_sc>::max();
	ptMax.h = ptMax.k = ptMax.l = -std::numeric_limits<t_real_sc>::max();
	ptMin.E = std::numeric_limits<t_real_sc>::max() * tl::get_one_meV<t_real_sc>();
	ptMax.E = -std::numeric_limits<t_real_sc>::max() * tl::get_one_meV<t_real_sc>();

	const std::size_t iNumPts = pInstr->GetScanCount();
	for(std::size_t iPt = 0; iPt < iNumPts; ++iPt)
	{
		const std::array<t_real_sc, 5> sc = pInstr->GetScanHKLKiKf(iPt);

		ScanPoint pt;
		pt.h = sc[0]; pt.k = sc[1]; pt.l = sc[2];
		pt.ki = sc[3] / tl::get_one_angstrom<t_real_sc>();
		pt.kf = sc[4] / tl::get_one_angstrom<t_real_sc>();
		pt.Ei = tl::k2E(pt.ki);
		pt.Ef = tl::k2E(pt.kf);
		pt.E = pt.Ei - pt.Ef;

		// component-wise minima and maxima
		if(pt.h > ptMax.h) ptMax.h = pt.h;
		if(pt.k > ptMax.k) ptMax.k = pt.k;
		if(pt.l > ptMax.l) ptMax.l = pt.l;
		if(pt.E > ptMax.E) ptMax.E = pt.E;

		if(pt.h < ptMin.h) ptMin.h = pt.h;
		if(pt.k < ptMin.k) ptMin.k = pt.k;
		if(pt.l < ptMin.l) ptMin.l = pt.l;
		if(pt.E < ptMin.E) ptMin.E = pt.E;

		if(bVerbose)
		{
			tl::log_info("Point ", iPt+1, ": ", "h=", pt.h, ", k=", pt.k, ", l=", pt.l,
				", ki=", t_real_sc(pt.ki * tl::get_one_angstrom<t_real_sc>()),
				", kf=", t_real_sc(pt.kf * tl::get_one_angstrom<t_real_sc>()),
				", E=", t_real_sc(pt.E / tl::get_one_meV<t_real_sc>())/*, ", Q=", pt.Q*tl::angstrom*/,
				", Cts=", scan.vecCts[iPt]/*, "+-", scan.vecCtsErr[iPt]*/,
				", Mon=", scan.vecMon[iPt]/*, "+-", scan.vecMonErr[iPt]*/, ".");
		}

		scan.vecPoints.emplace_back(std::move(pt));
	}


	const ScanPoint* ptBegin = &ptMin;
	const ScanPoint* ptEnd = &ptMax;

	// old behaviour: first and last points instead of component-wise min and max
	if(bUseFirstAndLastPoints)
	{
		ptBegin = &*scan.vecPoints.cbegin();
		ptEnd = &*scan.vecPoints.crbegin();
	}

	scan.vecScanOrigin[0] = ptBegin->h;
	scan.vecScanOrigin[1] = ptBegin->k;
	scan.vecScanOrigin[2] = ptBegin->l;
	scan.vecScanOrigin[3] = ptBegin->E / tl::get_one_meV<t_real_sc>();

	scan.vecScanDir[0] = ptEnd->h - ptBegin->h;
	scan.vecScanDir[1] = ptEnd->k - ptBegin->k;
	scan.vecScanDir[2] = ptEnd->l - ptBegin->l;
	scan.vecScanDir[3] = (ptEnd->E - ptBegin->E) / tl::get_one_meV<t_real_sc>();

	// determine principal scan axis
	const t_real_sc dEps = 0.01;

	for(unsigned iAx=0; iAx<4; ++iAx)
		scan.vecMainScanDir[iAx] = 0.;

	scan.m_iScIdx = 0;

	if(iScanAxis >= 1 && iScanAxis <= 4)
	{
		scan.m_iScIdx = iScanAxis-1;
		scan.vecMainScanDir[scan.m_iScIdx] = 1.;
	}
	else	// automatic determination
	{
		for(unsigned int i=0; i<4; ++i)
		{
			if(!tl::float_equal<t_real_sc>(scan.vecScanDir[i], 0., dEps))
			{
				scan.vecMainScanDir[i] = 1.;
				scan.m_iScIdx = i;
				break;
			}
		}
	}

	tl::log_info("Scan origin: (", scan.vecScanOrigin[0], " ", scan.vecScanOrigin[1], " ", scan.vecScanOrigin[2], " ", scan.vecScanOrigin[3], ").");
	tl::log_info("Scan dir: [", scan.vecScanDir[0], " ", scan.vecScanDir[1], " ", scan.vecScanDir[2], " ", scan.vecScanDir[3], "].");
	tl::log_info("Principal scan axis: [", scan.vecMainScanDir[0], " ", scan.vecMainScanDir[1], " ", scan.vecMainScanDir[2], " ", scan.vecMainScanDir[3], "].");


	if(scan.m_iScIdx >= 4)
	{
		tl::log_err("No scan variable found!");
		return false;
	}

	for(std::size_t iPt = 0; iPt < iNumPts; ++iPt)
	{
		const ScanPoint& pt = scan.vecPoints[iPt];

		t_real_sc dPos[] = { pt.h, pt.k, pt.l, pt.E/tl::get_one_meV<t_real_sc>() };
		scan.vecX.push_back(dPos[scan.m_iScIdx]);

		for(int ihklE=0; ihklE<4; ++ihklE)
			scan.vechklE[ihklE].push_back(dPos[ihklE]);
	}


	// filter
	decltype(scan.vecPoints) vecPointsNew;
	decltype(scan.vecX) vecXNew;
	decltype(scan.vecCts) vecCtsNew, vecMonNew, vecCtsErrNew, vecMonErrNew;
	decltype(scan.vechklE) vechklENew{};

	for(std::size_t i = 0; i < iNumPts; ++i)
	{
		// below lower limit?
		if(filter.dLower && scan.vecX[i] <= *filter.dLower)
			continue;
		// above upper limit?
		if(filter.dUpper && scan.vecX[i] >= *filter.dUpper)
			continue;

		// keep row if the value in the column equals the given one
		if(filter.colEquals && filter.colEquals->first!="" && filter.colEquals->second!="")
		{
			bool skip_pt = false;

			const std::string& colName = filter.colEquals->first;
			rex::regex rxCol(colName, rex::regex::ECMAScript | rex::regex_constants::icase);
			rex::smatch matchCol;

			const tl::FileInstrBase<t_real_sc>::t_vecColNames& colNames
				= pInstr->GetColNames();
			for(std::size_t col = 0; col < colNames.size(); ++col)
			{
				// get first matching column
				if(rex::regex_match(colNames[col], matchCol, rxCol))
				{
					const std::vector<t_real_sc>& filterCol = pInstr->GetCol(colNames[col]);

					const std::string& rowVal = filter.colEquals->second;
					rex::regex rxVal(rowVal, rex::regex::ECMAScript | rex::regex_constants::icase);
					rex::smatch matchVal;

					std::string curRowVal = tl::var_to_str(filterCol[i]);
					skip_pt = !rex::regex_match(curRowVal, matchVal, rxVal);

					break;
				}
			}

			if(skip_pt)
				continue;
		}

		vecPointsNew.push_back(scan.vecPoints[i]);
		vecXNew.push_back(scan.vecX[i]);
		vecCtsNew.push_back(scan.vecCts[i]);
		vecMonNew.push_back(scan.vecMon[i]);
		vecCtsErrNew.push_back(scan.vecCtsErr[i]);
		vecMonErrNew.push_back(scan.vecMonErr[i]);

		for(int ihklE=0; ihklE<4; ++ihklE)
			vechklENew[ihklE].push_back(scan.vechklE[ihklE][i]);
	}

	scan.vecPoints = std::move(vecPointsNew);
	scan.vecX = std::move(vecXNew);
	scan.vecCts = std::move(vecCtsNew);
	scan.vecMon = std::move(vecMonNew);
	scan.vecCtsErr = std::move(vecCtsErrNew);
	scan.vecMonErr = std::move(vecMonErrNew);

	for(int ihklE = 0; ihklE < 4; ++ihklE)
		scan.vechklE[ihklE] = std::move(vechklENew[ihklE]);

	return true;
}



/**
 * loading a scan file
 */
bool load_file(const char* pcFile, Scan& scan, bool bNormToMon, const Filter& filter,
	bool bFlip, bool bAllowScanMerging, bool bUseFirstAndLastPoints)
{
	std::vector<std::string> vec{pcFile};
	return load_file(vec, scan, bNormToMon, filter, bFlip, bAllowScanMerging, bUseFirstAndLastPoints);
}
