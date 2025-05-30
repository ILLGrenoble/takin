/**
 * scan viewer
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2025
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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
#include <string>
#include <algorithm>
#include <iterator>

#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/stat.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/file/file.h"
#include "tlibs/log/log.h"
#include "tlibs/helper/misc.h"


using t_real = t_real_glob;
using t_vec = tl::ublas::vector<t_real>;


void ScanViewerDlg::SetupPlotter(unsigned int numCurves)
{
	auto get_colour = [numCurves](unsigned int curve, bool is_curve) -> QColor
	{
		if(numCurves == 2)
		{
			// plotting single data set only
			if(is_curve)  // curve
				return QColor(0, 0, 0x99);
			else          // data points
				return QColor(0xff, 0, 0);
		}

		// plotting multiple data sets
		t_real lerpval = t_real(curve) / t_real(numCurves - 2);
		t_real r = tl::lerp<t_real, t_real>(1., 0., lerpval);
		t_real g = 0.;
		t_real b = tl::lerp<t_real, t_real>(0., 1., lerpval);

		return QColor(int(r * 255.), int(g * 255.), int(b * 255.));
	};

	// already set up with this number of curves?
	if(m_plotwrap && numCurves == m_plotwrap->GetNumCurves())
		return;

	QColor colorBck(240, 240, 240, 255);
	plot->setCanvasBackground(colorBck);

	bool show_legend = (numCurves > 2);  // only when plotting multiple scans
	m_plotwrap.reset(new QwtPlotWrapper(plot, numCurves, true, false, false, show_legend));

	for(unsigned int curve = 0; curve < numCurves; curve += 2)
	{
		// even indices are curves
		QPen penCurve;
		penCurve.setColor(get_colour(curve, true));
		penCurve.setWidth(2);
		m_plotwrap->GetCurve(curve)->setPen(penCurve);
		m_plotwrap->GetCurve(curve)->setStyle(QwtPlotCurve::CurveStyle::Lines);
		m_plotwrap->GetCurve(curve)->setTitle("Curve");
		m_plotwrap->GetCurve(curve)->setItemAttribute(QwtPlotCurve::Legend, false);

		// odd indices are scan points
		QPen penPoints;
		penPoints.setColor(get_colour(curve, false));
		penPoints.setWidth(4);
		m_plotwrap->GetCurve(curve + 1)->setPen(penPoints);
		m_plotwrap->GetCurve(curve + 1)->setStyle(QwtPlotCurve::CurveStyle::Dots);
		m_plotwrap->GetCurve(curve + 1)->setTitle("Data");
		m_plotwrap->GetCurve(curve + 1)->setItemAttribute(QwtPlotCurve::Legend, false);
	}

	m_plotwrap->GetPlot()->updateLegend();
}


void ScanViewerDlg::ClearPlot()
{
	for(tl::FileInstrBase<t_real_glob>* instr : m_instrs)
	{
		if(!instr)
			continue;
		delete instr;
	}

	m_instrs.clear();

	m_vecX.clear();
	m_vecY.clear();
	m_vecYErr.clear();
	m_vecFitX.clear();
	m_vecFitY.clear();

	// set dummy data
	static const std::vector<t_real> dummy;
	for(unsigned int curve = 0; curve < m_plotwrap->GetNumCurves(); ++curve)
		set_qwt_data<t_real>()(*m_plotwrap, dummy, dummy, curve, false);

	m_strX = m_strY = m_strCmd = "";
	plot->setAxisTitle(QwtPlot::xBottom, "");
	plot->setAxisTitle(QwtPlot::yLeft, "");
	plot->setTitle("");

	auto edits = { editA, editB, editC,
		editAlpha, editBeta, editGamma,
		editPlaneX0, editPlaneX1, editPlaneX2,
		editPlaneY0, editPlaneY1, editPlaneY2,
		editTitle, editSample,
		editUser, editContact,
		editKfix, editTimestamp };
	for(auto* pEdit : edits)
		pEdit->setText("");

	comboX->clear();
	comboY->clear();
	comboMon->clear();
	textExportedFile->clear();
	textRawFile->clear();
	spinStart->setValue(0);
	spinStop->setValue(0);
	spinSkip->setValue(0);

	m_plotwrap->GetPlot()->updateLegend();
	m_plotwrap->GetPlot()->replot();
}


void ScanViewerDlg::PlotScan()
{
	if(!m_bDoUpdate)
		return;

	SetupPlotter(m_instrs.size() * 2);  // curves and points

	bool bNormalise = checkNorm->isChecked();

	m_strX = comboX->itemData(comboX->currentIndex(), Qt::UserRole).toString().toStdString();
	m_strY = comboY->itemData(comboY->currentIndex(), Qt::UserRole).toString().toStdString();
	m_strMon = comboMon->itemData(comboMon->currentIndex(), Qt::UserRole).toString().toStdString();
	if(m_strX == "" || m_strY == "" || m_strMon == "")
		return;

	const int iStartIdx = spinStart->value();
	const int iEndSkip = spinStop->value();
	const int iSkipRows = spinSkip->value();
	m_vecX.resize(m_instrs.size());
	m_vecY.resize(m_instrs.size());
	m_vecYErr.resize(m_instrs.size());

	for(std::size_t instr_idx = 0; instr_idx < m_instrs.size(); ++instr_idx)
	{
		const tl::FileInstrBase<t_real_glob> *instr = m_instrs[instr_idx];
		if(!instr)
			continue;

		std::vector<t_real>& vecX = m_vecX[instr_idx];
		std::vector<t_real>& vecY = m_vecY[instr_idx];
		std::vector<t_real>& vecYErr = m_vecYErr[instr_idx];

		// get the data vectors
		vecX = GetCol(instr, m_strX);
		vecY = GetCol(instr, m_strY);
		std::vector<t_real> vecMon = GetCol(instr, m_strMon);

		// TODO: sort the data vectors
		//tl::sort_3(vecX.begin(), vecX.end(), vecY.begin(), vecMon.begin());

		bool bYIsACountVar = (m_strY == instr->GetCountVar() || m_strY == instr->GetMonVar());
		m_plotwrap->GetCurve(instr_idx*2 + 1)->SetShowErrors(bYIsACountVar);

		// see if there's a corresponding error column for the selected counter or monitor
		std::string ctr_err_col, mon_err_col;

		// get counter error if defined
		if(m_strY == instr->GetCountVar())
			ctr_err_col = instr->GetCountErr();
		else if(m_strY == instr->GetMonVar())
			ctr_err_col = instr->GetMonErr();

		if(ctr_err_col != "")
		{
			// use given error column
			vecYErr = GetCol(instr, ctr_err_col);
		}
		else
		{
			vecYErr.clear();
			vecYErr.reserve(vecY.size());

			// calculate error
			for(std::size_t iY = 0; iY < vecY.size(); ++iY)
			{
				t_real err = tl::float_equal(vecY[iY], t_real(0), g_dEps) ? t_real(1) : std::sqrt(std::abs(vecY[iY]));
				vecYErr.push_back(err);
			}
		}

		// get monitor error if defined
		std::vector<t_real> vecMonErr;
		if(m_strMon == instr->GetCountVar())
			mon_err_col = instr->GetCountErr();
		else if(m_strMon == instr->GetMonVar())
			mon_err_col = instr->GetMonErr();

		if(mon_err_col != "")
		{
			// use given error column
			vecMonErr = GetCol(instr, mon_err_col);
		}
		else
		{
			vecMonErr.reserve(vecY.size());

			// calculate error
			for(std::size_t iY = 0; iY < vecMon.size(); ++iY)
			{
				t_real err = tl::float_equal(vecMon[iY], t_real(0), g_dEps) ? t_real(1) : std::sqrt(std::abs(vecMon[iY]));
				vecMonErr.push_back(err);
			}
		}


		// remove points from start
		if(iStartIdx != 0)
		{
			if(std::size_t(iStartIdx) >= vecX.size())
				vecX.clear();
			else
				vecX.erase(vecX.begin(), vecX.begin() + iStartIdx);

			if(std::size_t(iStartIdx) >= vecY.size())
			{
				vecY.clear();
				vecMon.clear();
				vecYErr.clear();
				vecMonErr.clear();
			}
			else
			{
				vecY.erase(vecY.begin(), vecY.begin() + iStartIdx);
				vecMon.erase(vecMon.begin(), vecMon.begin() + iStartIdx);
				vecYErr.erase(vecYErr.begin(), vecYErr.begin() + iStartIdx);
				vecMonErr.erase(vecMonErr.begin(), vecMonErr.begin() + iStartIdx);
			}
		}

		// remove points from end
		if(iEndSkip != 0)
		{
			if(std::size_t(iEndSkip) >= vecX.size())
				vecX.clear();
			else
				vecX.erase(vecX.end() - iEndSkip, vecX.end());

			if(std::size_t(iEndSkip) >= vecY.size())
			{
				vecY.clear();
				vecMon.clear();
				vecYErr.clear();
				vecMonErr.clear();
			}
			else
			{
				vecY.erase(vecY.end() - iEndSkip, vecY.end());
				vecMon.erase(vecMon.end() - iEndSkip, vecMon.end());
				vecYErr.erase(vecYErr.end() - iEndSkip, vecYErr.end());
				vecMonErr.erase(vecMonErr.end() - iEndSkip, vecMonErr.end());
			}
		}

		// interleave rows
		if(iSkipRows != 0)
		{
			std::vector<t_real> vecXNew, vecYNew, vecMonNew, vecYErrNew, vecMonErrNew;

			for(std::size_t iRow = 0; iRow < std::min(vecX.size(), vecY.size()); ++iRow)
			{
				vecXNew.push_back(vecX[iRow]);
				vecYNew.push_back(vecY[iRow]);
				vecMonNew.push_back(vecMon[iRow]);
				vecYErrNew.push_back(vecYErr[iRow]);
				vecMonErrNew.push_back(vecMonErr[iRow]);

				iRow += iSkipRows;
			}

			vecX = std::move(vecXNew);
			vecY = std::move(vecYNew);
			vecMon = std::move(vecMonNew);
			vecYErr = std::move(vecYErrNew);
			vecMonErr = std::move(vecMonErrNew);
		}


		// errors
		if(vecMon.size() != vecY.size() || vecMonErr.size() != vecYErr.size())
		{
			bNormalise = false;
			tl::log_err("Counter and monitor data count do not match, cannot normalise.");
		}

		// normalise to monitor?
		if(bNormalise)
		{
			for(std::size_t iY = 0; iY < vecY.size(); ++iY)
			{
				if(tl::float_equal(vecMon[iY], t_real(0), g_dEps))
				{
					tl::log_warn("Monitor counter is zero for point ", iY + 1, ".");

					vecY[iY] = 0.;
					vecYErr[iY] = 1.;
				}
				else
				{
					std::tie(vecY[iY], vecYErr[iY]) = tl::norm_cnts_to_mon(
						vecY[iY], vecYErr[iY], vecMon[iY], vecMonErr[iY]);
				}
			}
		}


		// show fit (for fist scan file)
		if(m_vecFitX.size() && instr_idx == 0)
			set_qwt_data<t_real>()(*m_plotwrap, m_vecFitX, m_vecFitY, instr_idx*2, false);
		else
			set_qwt_data<t_real>()(*m_plotwrap, vecX, vecY, instr_idx*2, false);
		set_qwt_data<t_real>()(*m_plotwrap, vecX, vecY, instr_idx*2 + 1, false, &vecYErr);

		// legend
		m_plotwrap->GetCurve(instr_idx*2 + 1)->setTitle(instr->GetScanNumber().c_str());
		m_plotwrap->GetCurve(instr_idx*2 + 0)->setItemAttribute(QwtPlotCurve::Legend, false);
		m_plotwrap->GetCurve(instr_idx*2 + 1)->setItemAttribute(QwtPlotCurve::Legend, m_instrs.size() > 1);
	}

	// labels
	QString strY = m_strY.c_str();
	if(bNormalise)
	{
		strY += " / ";
		strY += m_strMon.c_str();
	}
	plot->setAxisTitle(QwtPlot::xBottom, m_strX.c_str());
	plot->setAxisTitle(QwtPlot::yLeft, strY);
	plot->setTitle(m_strCmd.c_str());

	// replot
	set_zoomer_base(m_plotwrap->GetZoomer(), m_vecX, m_vecY, false, m_plotwrap.get(), true, &m_vecYErr);
	m_plotwrap->GetPlot()->updateLegend();
	m_plotwrap->GetPlot()->replot();

	GenerateExternal(comboExport->currentIndex());
}
