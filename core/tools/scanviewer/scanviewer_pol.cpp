/**
 * scan viewer -- polarisation matrix calculation
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2023
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

#include "scanviewer.h"

#include <iostream>
#include <unordered_set>
#include <string>

#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"

using t_real = t_real_glob;


/**
 * calculates the polarisation matrix (of the first selected file)
 */
void ScanViewerDlg::CalcPol()
{
	editPolMat->clear();

	if(!m_instrs.size())
		return;

	tl::FileInstrBase<t_real_glob> *instr = m_instrs[0];
	if(!instr)
		return;

	const std::string strPolVec1 = editPolVec1->text().toStdString();
	const std::string strPolVec2 = editPolVec2->text().toStdString();
	const std::string strPolCur1 = editPolCur1->text().toStdString();
	const std::string strPolCur2 = editPolCur2->text().toStdString();
	const std::string strFlip1 = editFlip1->text().toStdString();
	const std::string strFlip2 = editFlip2->text().toStdString();
	const std::string strXYZ = editXYZ->text().toStdString();

	instr->SetPolNames(strPolVec1.c_str(), strPolVec2.c_str(), strPolCur1.c_str(), strPolCur2.c_str());
	instr->SetLinPolNames(strFlip1.c_str(), strFlip2.c_str(), strXYZ.c_str());
	instr->ParsePolData();

	const std::vector<std::array<t_real, 6>>& vecPolStates = instr->GetPolStates();
	const std::size_t iNumPolStates = vecPolStates.size();
	if(iNumPolStates == 0)
	{
		editPolMat->setHtml("<html><body><font size=\"5\" color=\"#ff0000\">No polarisation data.</font></body></html>");
		return;
	}


	// get the SF state to a given NSF state
	auto find_spinflip_state_idx = [&vecPolStates, iNumPolStates]
		(const std::array<t_real, 6>& state) -> std::size_t
	{
		t_real dPix1 = state[0], dPiy1 = state[1], dPiz1 = state[2];
		t_real dPfx1 = state[3], dPfy1 = state[4], dPfz1 = state[5];

		for(std::size_t iPol=0; iPol<iNumPolStates; ++iPol)
		{
			t_real dPix2 = vecPolStates[iPol][0];
			t_real dPiy2 = vecPolStates[iPol][1];
			t_real dPiz2 = vecPolStates[iPol][2];
			t_real dPfx2 = vecPolStates[iPol][3];
			t_real dPfy2 = vecPolStates[iPol][4];
			t_real dPfz2 = vecPolStates[iPol][5];

			if(tl::float_equal(dPix1, dPix2, g_dEps) &&
				tl::float_equal(dPiy1, dPiy2, g_dEps) &&
				tl::float_equal(dPiz1, dPiz2, g_dEps) &&
				tl::float_equal(dPfx1, -dPfx2, g_dEps) &&
				tl::float_equal(dPfy1, -dPfy2, g_dEps) &&
				tl::float_equal(dPfz1, -dPfz2, g_dEps))
			{
				return iPol;
			}
		}

		// none found
		return iNumPolStates;
	};


	auto propagate_pol_err = [](t_real x, t_real y, t_real dx, t_real dy) -> t_real
	{
		// d((x-y)/(x+y)) = dx * 2*y/(x+y)^2 - dy * 2*x/(x+y)^2
		return std::sqrt(dx*2.*y/((x+y)*(x+y))*dx*2.*y/((x+y)*(x+y))
			+ dy*2.*x/((x+y)*(x+y))*dy*2.*x/((x+y)*(x+y)));
	};


	// convert polarisation vector to string representation
	auto polvec_str = [](t_real x, t_real y, t_real z) -> std::string
	{
		std::ostringstream ostr;
		ostr.precision(g_iPrec);

		if(tl::float_equal<t_real>(x, 1., g_dEps) &&
			tl::float_equal<t_real>(y, 0., g_dEps) &&
			tl::float_equal<t_real>(z, 0., g_dEps))
			ostr << "x";
		else if(tl::float_equal<t_real>(x, -1., g_dEps) &&
			tl::float_equal<t_real>(y, 0., g_dEps) &&
			tl::float_equal<t_real>(z, 0., g_dEps))
			ostr << "-x";
		else if(tl::float_equal<t_real>(x, 0., g_dEps) &&
			tl::float_equal<t_real>(y, 1., g_dEps) &&
			tl::float_equal<t_real>(z, 0., g_dEps))
			ostr << "y";
		else if(tl::float_equal<t_real>(x, 0., g_dEps) &&
			tl::float_equal<t_real>(y, -1., g_dEps) &&
			tl::float_equal<t_real>(z, 0., g_dEps))
			ostr << "-y";
		else if(tl::float_equal<t_real>(x, 0., g_dEps) &&
			tl::float_equal<t_real>(y, 0., g_dEps) &&
			tl::float_equal<t_real>(z, 1., g_dEps))
			ostr << "z";
		else if(tl::float_equal<t_real>(x, 0., g_dEps) &&
			tl::float_equal<t_real>(y, 0., g_dEps) &&
			tl::float_equal<t_real>(z, -1., g_dEps))
			ostr << "-z";
		else
			ostr << "[" << x << " " << y << " " << z << "]";

		return ostr.str();
	};


	std::vector<bool> vecHasSFPartner;
	// indices to spin-flipped states
	std::vector<std::size_t> vecSFIdx;

	for(std::size_t iPol=0; iPol<iNumPolStates; ++iPol)
	{
		const std::array<t_real, 6>& state = vecPolStates[iPol];
		std::size_t iIdx = find_spinflip_state_idx(state);

		vecHasSFPartner.push_back(iIdx < iNumPolStates);
		vecSFIdx.push_back(iIdx);
	}


	const std::vector<std::string> vecScanVars = instr->GetScannedVars();
	std::string strX;
	if(vecScanVars.size())
		strX = vecScanVars[0];
	const std::vector<t_real>& vecX = instr->GetCol(strX.c_str());
	const std::vector<t_real>& vecCnts = instr->GetCol(instr->GetCountVar().c_str());
	const std::vector<t_real>& vecMons = instr->GetCol(instr->GetMonVar().c_str());
	bool has_mon = (vecCnts.size() == vecMons.size());

	// raw counts per polarisation channel
	std::ostringstream ostrCnts;
	ostrCnts.precision(g_iPrec);
	ostrCnts << "<p><h2>Counts in Polarisation Channels</h2>";

	// iterate over scan points
	for(std::size_t iPt=0; iPt<vecCnts.size();)
	{
		ostrCnts << "<p><b>Scan Point " << (iPt/iNumPolStates+1);
		if(strX != "" && iPt < vecX.size())
			ostrCnts << " (" << strX << " = " << vecX[iPt + 0] << ")";
		ostrCnts << "</b>";
		ostrCnts << "<table border=\"1\" cellpadding=\"0\">";
		ostrCnts << "<tr><th>Init. Pol. Vec.</th>";
		ostrCnts << "<th>Fin. Pol. Vec.</th>";
		ostrCnts << "<th>Counts</th>";
		ostrCnts << "<th>Error</th>";
		if(has_mon)
		{
			ostrCnts << "<th>Monitor</th>";
			ostrCnts << "<th>Norm. Counts</th>";
			ostrCnts << "<th>Norm. Error</th>";
		}
		ostrCnts << "</tr>";

		// iterate over polarisation states
		for(std::size_t iPol=0; iPol<iNumPolStates; ++iPol, ++iPt)
		{
			t_real dPix = vecPolStates[iPol][0];
			t_real dPiy = vecPolStates[iPol][1];
			t_real dPiz = vecPolStates[iPol][2];

			t_real dPfx = vecPolStates[iPol][3];
			t_real dPfy = vecPolStates[iPol][4];
			t_real dPfz = vecPolStates[iPol][5];

			std::size_t iCnts = 0, iMon = 0;
			if(iPt < vecCnts.size())
				iCnts = std::size_t(vecCnts[iPt]);
			if(iPt < vecMons.size())
				iMon = std::size_t(vecMons[iPt]);
			t_real dErr = (iCnts==0 ? 1 : std::sqrt(t_real(iCnts)));
			t_real dMonErr = (iMon==0 ? 1 : std::sqrt(t_real(iMon)));

			ostrCnts << "<tr><td>" << polvec_str(dPix, dPiy, dPiz) << "</td>";
			ostrCnts << "<td>" << polvec_str(dPfx, dPfy, dPfz) << "</td>";
			ostrCnts << "<td><b>" << iCnts << "</b></td>";
			ostrCnts << "<td><b>" << dErr << "</b></td>";
			if(has_mon)
			{
				t_real norm_cts = 0., norm_err = 0.;
				std::tie(norm_cts, norm_err) = norm_cnts_to_mon(
					t_real(iCnts), dErr, t_real(iMon), dMonErr);

				ostrCnts << "<td><b>" << iMon << "</b></td>";
				ostrCnts << "<td><b>" << norm_cts << "</b></td>";
				ostrCnts << "<td><b>" << norm_err << "</b></td>";
			}
			ostrCnts << "</tr>";
		}
		ostrCnts << "</table></p>";
	}
	ostrCnts << "</p><br><hr><br>";


	// polarisation matrix elements
	std::ostringstream ostrPol;
	ostrPol.precision(g_iPrec);
	ostrPol << "<p><h2>Point-wise Polarisation Matrix Elements</h2>";
	bool bHasAnyData = false;

	// iterate over scan points
	for(std::size_t iPt=0; iPt<vecCnts.size()/iNumPolStates; ++iPt)
	{
		ostrPol << "<p><b>Scan Point " << (iPt+1);
		if(strX != "" && iPt*iNumPolStates < vecX.size())
			ostrPol << " (" << strX << " = " << vecX[iPt*iNumPolStates + 0] << ")";
		ostrPol << "</b>";
		ostrPol << "<table border=\"1\" cellpadding=\"0\">";
		ostrPol << "<tr><th>Index 1</th>";
		ostrPol << "<th>Index 2</th>";
		ostrPol << "<th>Polarisation</th>";
		ostrPol << "<th>Error</th></tr>";

		// iterate over all polarisation states which have a SF partner
		std::unordered_set<std::size_t> setPolAlreadySeen;
		for(std::size_t iPol=0; iPol<iNumPolStates; ++iPol)
		{
			if(!vecHasSFPartner[iPol]) continue;
			if(setPolAlreadySeen.find(iPol) != setPolAlreadySeen.end())
				continue;

			const std::size_t iSF = vecSFIdx[iPol];
			const std::array<t_real, 6>& state = vecPolStates[iPol];

			setPolAlreadySeen.insert(iPol);
			setPolAlreadySeen.insert(iSF);

			t_real dCntsNSF = t_real(0);
			t_real dCntsSF = t_real(0);

			if(iPt*iNumPolStates + iPol < vecCnts.size())
				dCntsNSF = vecCnts[iPt*iNumPolStates + iPol];
			if(iPt*iNumPolStates + iSF < vecCnts.size())
				dCntsSF = vecCnts[iPt*iNumPolStates + iSF];
			t_real dNSFErr = std::sqrt(dCntsNSF);
			t_real dSFErr = std::sqrt(dCntsSF);
			if(tl::float_equal(dCntsNSF, t_real(0), g_dEps))
				dNSFErr = 1.;
			if(tl::float_equal(dCntsSF, t_real(0), g_dEps))
				dSFErr = 1.;

			// normalise to monitor if available
			if(has_mon)
			{
				t_real dMonNSF = t_real(0);
				t_real dMonSF = t_real(0);

				if(iPt*iNumPolStates + iPol < vecMons.size())
					dMonNSF = vecMons[iPt*iNumPolStates + iPol];
				if(iPt*iNumPolStates + iSF < vecMons.size())
					dMonSF = vecMons[iPt*iNumPolStates + iSF];

				t_real dMonNSFErr = std::sqrt(dMonNSF);
				t_real dMonSFErr = std::sqrt(dMonSF);

				if(tl::float_equal(dMonNSF, t_real(0), g_dEps))
					dMonNSFErr = 1.;
				if(tl::float_equal(dMonSF, t_real(0), g_dEps))
					dMonSFErr = 1.;

				std::tie(dCntsNSF, dNSFErr) = norm_cnts_to_mon(
					dCntsNSF, dNSFErr, dMonNSF, dMonNSFErr);
				std::tie(dCntsSF, dSFErr) = norm_cnts_to_mon(
					dCntsSF, dSFErr, dMonSF, dMonSFErr);
			}

			bool bInvalid = tl::float_equal(dCntsNSF+dCntsSF, t_real(0), g_dEps);
			t_real dPolElem = 0., dPolErr = 1.;
			if(!bInvalid)
			{
				dPolElem = /*std::abs*/((dCntsSF-dCntsNSF) / (dCntsSF+dCntsNSF));
				dPolErr = propagate_pol_err(dCntsNSF, dCntsSF, dNSFErr, dSFErr);
			}

			// polarisation matrix elements, e.g. <[100] | P | [010]> = <x|P|y>
			ostrPol << "<tr><td>" << polvec_str(state[0], state[1], state[2]) << "</td>"
				<< "<td>" << polvec_str(state[3], state[4], state[5]) << "</td>"
				<< "<td><b>" << (bInvalid ? "--- ": tl::var_to_str(dPolElem, g_iPrec)) << "</b></td>"
				<< "<td><b>" << (bInvalid ? "--- ": tl::var_to_str(dPolErr, g_iPrec)) << "</b></td>"
				<< "</tr>";

			bHasAnyData = true;
		}
		ostrPol << "</table></p>";
	}
	ostrPol << "</p>";

	if(!bHasAnyData)
		ostrPol << "<font size=\"5\" color=\"#ff0000\">Insufficient Data.</font>";

	ostrPol << "<br><hr><br>";


	// maybe a peak-background pair
	std::ostringstream ostrPeakBkgrd;

	if(bHasAnyData && vecCnts.size()/iNumPolStates == 2)
	{
		ostrPeakBkgrd.precision(g_iPrec);
		ostrPeakBkgrd << "<p><h2>Polarisation Matrix Elements for Peak-Background</h2>";

		std::size_t iPtFg=0, iPtBg=1;
		if(vecCnts[0*iNumPolStates + 0] >= vecCnts[1*iNumPolStates + 0])
		{
			iPtFg = 0;
			iPtBg = 1;
		}
		else
		{
			iPtFg = 1;
			iPtBg = 0;
		}

		ostrPeakBkgrd << "<p><b>Foreground Scan Point: " << (iPtFg+1) << ", Background: " << (iPtBg+1);
		if(strX != "" && iPtFg*iNumPolStates < vecX.size())
			ostrPeakBkgrd << " (" << strX << " = " << vecX[iPtFg*iNumPolStates + 0] << ")";
		ostrPeakBkgrd << "</b>";
		ostrPeakBkgrd << "<table border=\"1\" cellpadding=\"0\">";
		ostrPeakBkgrd << "<tr><th>Index 1</th>";
		ostrPeakBkgrd << "<th>Index 2</th>";
		ostrPeakBkgrd << "<th>Polarisation</th>";
		ostrPeakBkgrd << "<th>Error</th></tr>";

		// iterate over all polarisation states which have a SF partner
		std::unordered_set<std::size_t> setPolAlreadySeen;
		for(std::size_t iPol=0; iPol<iNumPolStates; ++iPol)
		{
			if(!vecHasSFPartner[iPol]) continue;
			if(setPolAlreadySeen.find(iPol) != setPolAlreadySeen.end())
				continue;

			const std::size_t iSF = vecSFIdx[iPol];
			const std::array<t_real, 6>& state = vecPolStates[iPol];

			setPolAlreadySeen.insert(iPol);
			setPolAlreadySeen.insert(iSF);

			t_real dCntsNSF_fg = t_real{0};
			t_real dCntsNSF_bg = t_real{0};
			t_real dCntsSF_fg = t_real{0};
			t_real dCntsSF_bg = t_real{0};

			if(iPtFg*iNumPolStates + iPol < vecCnts.size())
				dCntsNSF_fg = vecCnts[iPtFg*iNumPolStates + iPol];
			if(iPtBg*iNumPolStates + iPol < vecCnts.size())
				dCntsNSF_bg = vecCnts[iPtBg*iNumPolStates + iPol];
			if(iPtFg*iNumPolStates + iSF < vecCnts.size())
				dCntsSF_fg = vecCnts[iPtFg*iNumPolStates + iSF];
			if(iPtBg*iNumPolStates + iSF < vecCnts.size())
				dCntsSF_bg = vecCnts[iPtBg*iNumPolStates + iSF];

			const t_real dCntsNSF = dCntsNSF_fg - dCntsNSF_bg;
			const t_real dCntsSF = dCntsSF_fg - dCntsSF_bg;
			t_real dNSFErr = std::sqrt(dCntsNSF_fg + dCntsNSF_bg);
			t_real dSFErr = std::sqrt(dCntsSF_fg + dCntsSF_bg);

			if(tl::float_equal(dCntsNSF, t_real(0), g_dEps))
				dNSFErr = 1.;
			if(tl::float_equal(dCntsSF, t_real(0), g_dEps))
				dSFErr = 1.;

			bool bInvalid = tl::float_equal(dCntsNSF+dCntsSF, t_real(0), g_dEps);
			t_real dPolElem = 0., dPolErr = 1.;
			if(!bInvalid)
			{
				dPolElem = (dCntsSF-dCntsNSF) / (dCntsSF+dCntsNSF);
				dPolErr = propagate_pol_err(dCntsNSF, dCntsSF, dNSFErr, dSFErr);
			}

			// polarisation matrix elements, e.g. <[100] | P | [010]> = <x|P|y>
			ostrPeakBkgrd << "<tr><td>" << polvec_str(state[0], state[1], state[2]) << "</td>"
				<< "<td>" << polvec_str(state[3], state[4], state[5]) << "</td>"
				<< "<td><b>" << (bInvalid ? "--- ": tl::var_to_str(dPolElem, g_iPrec)) << "</b></td>"
				<< "<td><b>" << (bInvalid ? "--- ": tl::var_to_str(dPolErr, g_iPrec)) << "</b></td>"
				<< "</tr>";
		}

		ostrPeakBkgrd << "</table></p>";
		ostrPeakBkgrd << "</p><br><hr><br>";
	}


	std::string strHtml = "<html><body>" +
		ostrPeakBkgrd.str() + ostrPol.str() + ostrCnts.str() +
		std::string{"</body></html>"};
	editPolMat->setHtml(QString::fromUtf8(strHtml.c_str()));
}
