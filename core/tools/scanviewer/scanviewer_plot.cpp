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
#include <sstream>
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
		m_plotwrap->GetCurve(curve + 0)->setPen(penCurve);
		m_plotwrap->GetCurve(curve + 0)->setStyle(QwtPlotCurve::CurveStyle::Lines);
		m_plotwrap->GetCurve(curve + 0)->setTitle("Curve");
		m_plotwrap->GetCurve(curve + 0)->setItemAttribute(QwtPlotCurve::Legend, false);

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

	m_strX = m_strY = m_strMon = m_strCmd = "";
	m_strX2 = m_strY2 = m_strMon2 = "";
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

	for(QComboBox* combo : { comboX, comboX2, comboY, comboY2, comboMon, comboMon2 })
		combo->clear();
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

	// sub-curves per scan
	std::size_t num_sub_curves = checkCurve2->isChecked() ? 2 : 1;
	SetupPlotter(m_instrs.size() * num_sub_curves * 2);  // curves and points

	m_strX = comboX->itemData(comboX->currentIndex(), Qt::UserRole).toString().toStdString();
	m_strY = comboY->itemData(comboY->currentIndex(), Qt::UserRole).toString().toStdString();
	m_strMon = comboMon->itemData(comboMon->currentIndex(), Qt::UserRole).toString().toStdString();
	if(m_strX == "" || m_strY == "" || m_strMon == "")
		return;

	if(num_sub_curves > 1)
	{
		m_strX2 = comboX->itemData(comboX2->currentIndex(), Qt::UserRole).toString().toStdString();
		m_strY2 = comboY->itemData(comboY2->currentIndex(), Qt::UserRole).toString().toStdString();
		m_strMon2 = comboMon->itemData(comboMon2->currentIndex(), Qt::UserRole).toString().toStdString();
	}

	bool bNormalise = checkNorm->isChecked();
	const int iStartIdx = spinStart->value();
	const int iEndSkip = spinStop->value();
	const int iSkipRows = spinSkip->value();
	m_vecX.resize(m_instrs.size() * num_sub_curves);
	m_vecY.resize(m_instrs.size() * num_sub_curves);
	m_vecYErr.resize(m_instrs.size() * num_sub_curves);

	// plot the individual scan files
	for(std::size_t scanfile_idx = 0; scanfile_idx < m_instrs.size(); ++scanfile_idx)
	{
		const tl::FileInstrBase<t_real_glob> *instr = m_instrs[scanfile_idx];
		if(!instr)
			continue;

		auto plot_data = [this, instr, bNormalise, iStartIdx, iEndSkip, iSkipRows, num_sub_curves](
			std::size_t plot_idx, std::size_t plot_sub_idx,
			const std::string& strX, const std::string& strY, const std::string& strMon)
		{
			std::vector<t_real>& vecX = m_vecX[plot_idx];
			std::vector<t_real>& vecY = m_vecY[plot_idx];
			std::vector<t_real>& vecYErr = m_vecYErr[plot_idx];

			// get the data vectors
			vecX = GetCol(instr, strX);
			vecY = GetCol(instr, strY);
			std::vector<t_real> vecMon = GetCol(instr, strMon);

			// TODO: sort the data vectors
			//tl::sort_3(vecX.begin(), vecX.end(), vecY.begin(), vecMon.begin());

			bool bYIsACountVar = (strY == instr->GetCountVar() || strY == instr->GetMonVar());
			m_plotwrap->GetCurve(plot_idx*2 + 1)->SetShowErrors(bYIsACountVar);

			// see if there's a corresponding error column for the selected counter or monitor
			std::string ctr_err_col, mon_err_col;

			// get counter error if defined
			if(strY == instr->GetCountVar())
				ctr_err_col = instr->GetCountErr();
			else if(strY == instr->GetMonVar())
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
			if(strMon == instr->GetCountVar())
				mon_err_col = instr->GetCountErr();
			else if(strMon == instr->GetMonVar())
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
				std::size_t new_size = std::min(vecX.size(), vecY.size());

				std::vector<t_real> vecXNew, vecYNew, vecMonNew, vecYErrNew, vecMonErrNew;
				vecXNew.reserve(new_size);
				vecYNew.reserve(new_size);
				vecMonNew.reserve(new_size);
				vecYErrNew.reserve(new_size);
				vecMonErrNew.reserve(new_size);

				for(std::size_t iRow = 0; iRow < new_size; ++iRow)
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
			bool bNorm = bNormalise;
			if(vecMon.size() != vecY.size() || vecMonErr.size() != vecYErr.size())
			{
				bNorm = false;
				tl::log_err("Counter and monitor data count do not match, cannot normalise.");
			}

			// normalise to monitor?
			if(bNorm)
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


			// show fit (for first scan file)
			if(m_vecFitX.size() && plot_idx == 0)
				set_qwt_data<t_real>()(*m_plotwrap, m_vecFitX, m_vecFitY, plot_idx*2, false);
			else
				set_qwt_data<t_real>()(*m_plotwrap, vecX, vecY, plot_idx*2, false);
			set_qwt_data<t_real>()(*m_plotwrap, vecX, vecY, plot_idx*2 + 1, false, &vecYErr);

			// legend
			std::ostringstream ostrLegend;
			ostrLegend << instr->GetScanNumber();
			if(num_sub_curves > 1)
				ostrLegend << " (" << plot_sub_idx + 1 << ")";
			m_plotwrap->GetCurve(plot_idx*2 + 1)->setTitle(ostrLegend.str().c_str());
			m_plotwrap->GetCurve(plot_idx*2 + 0)->setItemAttribute(QwtPlotCurve::Legend, false);
			m_plotwrap->GetCurve(plot_idx*2 + 1)->setItemAttribute(QwtPlotCurve::Legend, num_sub_curves > 1);
		};


		// plot the individual sub-curves
		std::size_t plot_idx = scanfile_idx * num_sub_curves;

		for(std::size_t sub_idx = 0; sub_idx < num_sub_curves; ++sub_idx)
		{
			const std::string& strX = (sub_idx == 1 ? m_strX2 : m_strX);
			const std::string& strY = (sub_idx == 1 ? m_strY2 : m_strY);
			const std::string& strMon = (sub_idx == 1 ? m_strMon2 : m_strMon);

			plot_data(plot_idx + sub_idx, sub_idx, strX, strY, strMon);
		}
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
