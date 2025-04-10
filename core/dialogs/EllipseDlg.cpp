/**
 * Ellipse Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013 - 2024
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "EllipseDlg.h"

#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/time/chrono.h"
#include "tlibs/helper/flags.h"
#include "libs/version.h"

#include <future>
#include <fstream>

#include <QFileDialog>


#define ELLIPSE_DLG_TITLE "Resolution Ellipses"


EllipseDlg::EllipseDlg(QWidget* pParent, QSettings* pSett, Qt::WindowFlags fl)
	: QDialog(pParent, fl), m_pSettings(pSett)
{
	setupUi(this);
	setWindowTitle(ELLIPSE_DLG_TITLE);
	setSizeGripEnabled(true);

	if(m_pSettings)
	{
		// font
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		// window geometry
		if(m_pSettings->contains("reso/ellipse_geo"))
			restoreGeometry(m_pSettings->value("reso/ellipse_geo").toByteArray());

		m_bCenterOn0 = m_pSettings->value("reso/center_around_origin", 1).toInt() != 0;
	}

	m_vecplotwrap.reserve(4);
	m_elliProj.resize(4);
	m_elliSlice.resize(4);
	m_vecXCurvePoints.resize(8);
	m_vecYCurvePoints.resize(8);
	m_vecMCXCurvePoints.resize(4);
	m_vecMCYCurvePoints.resize(4);

	// generate plots
	QwtPlot* pPlots[] = { plot1, plot2, plot3, plot4 };
	for(unsigned int i = 0; i < 4; ++i)
	{
		m_vecplotwrap.push_back(std::unique_ptr<QwtPlotWrapper>(new QwtPlotWrapper(pPlots[i], 3)));
		m_vecplotwrap[i]->GetPlot()->setMinimumSize(200,200);

		m_vecplotwrap[i]->GetCurve(0)->setTitle("Projected Ellipse (HWHM Contour)");
		m_vecplotwrap[i]->GetCurve(1)->setTitle("Sliced Ellipse (HWHM Contour)");

		QPen penProj, penSlice, penPoints;
		penProj.setColor(QColor(0, 0x99, 0));
		penSlice.setColor(QColor(0, 0, 0x99));
		penPoints.setColor(QColor(0xff, 0, 0, 0xff));
		penProj.setWidth(2);
		penSlice.setWidth(2);
		penPoints.setWidth(1);

		m_vecplotwrap[i]->GetCurve(0)->setStyle(QwtPlotCurve::CurveStyle::Dots);
		m_vecplotwrap[i]->GetCurve(1)->setStyle(QwtPlotCurve::CurveStyle::Lines);
		m_vecplotwrap[i]->GetCurve(2)->setStyle(QwtPlotCurve::CurveStyle::Lines);

		m_vecplotwrap[i]->GetCurve(0)->setPen(penPoints);
		m_vecplotwrap[i]->GetCurve(1)->setPen(penProj);
		m_vecplotwrap[i]->GetCurve(2)->setPen(penSlice);

		if(m_vecplotwrap[i]->HasTrackerSignal())
		{
			connect(m_vecplotwrap[i]->GetPicker(), &QwtPlotPicker::moved,
				this, &EllipseDlg::cursorMoved);
		}
	}

	// connections
	QObject::connect(comboCoord, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
		this, &EllipseDlg::Calc);
	QObject::connect(checkCenter, static_cast<void(QCheckBox::*)(bool)>(&QCheckBox::toggled),
		this, &EllipseDlg::SetCenterOn0);
	QObject::connect(btnSave, &QPushButton::clicked, this, &EllipseDlg::SaveEllipses);

	m_bReady = true;
}


EllipseDlg::~EllipseDlg()
{
	m_bReady = false;
	m_vecplotwrap.clear();
}


void EllipseDlg::SetTitle(const char* pcTitle)
{
	QString strTitle = ELLIPSE_DLG_TITLE;
	strTitle += " - ";
	strTitle += pcTitle;
	this->setWindowTitle(strTitle);
}


void EllipseDlg::cursorMoved(const QPointF& pt)
{
	std::string strX = tl::var_to_str(pt.x(), g_iPrecGfx);
	std::string strY = tl::var_to_str(pt.y(), g_iPrecGfx);

	std::ostringstream ostr;
	ostr << "(" << strX << ", " << strY << ")";

	this->labelStatus->setText(ostr.str().c_str());
}


void EllipseDlg::Calc()
{
	if(!m_bReady)
		return;
	const EllipseCoordSys coord = static_cast<EllipseCoordSys>(comboCoord->currentIndex());

	const ublas::matrix<t_real_reso> *pReso = nullptr;
	const ublas::vector<t_real_reso> *pReso_v = nullptr;
	const ublas::vector<t_real_reso> *pQavg = nullptr;
	const std::vector<ublas::vector<t_real_reso>> *pvecMC = nullptr;

	switch(coord)
	{
		case EllipseCoordSys::Q_AVG:	// Q|| Qperp system in 1/A
			pReso = &m_reso;
			pQavg = &m_Q_avg;
			pReso_v = &m_reso_v;
			pvecMC = m_params.vecMC_direct;
			break;
		case EllipseCoordSys::RLU:		// rlu system
			pReso = &m_resoHKL;
			pQavg = &m_Q_avgHKL;
			pReso_v = &m_reso_vHKL;
			pvecMC = m_params.vecMC_HKL;
			break;
		case EllipseCoordSys::RLU_ORIENT:	// rlu system
			pReso = &m_resoOrient;
			pQavg = &m_Q_avgOrient;
			pReso_v = &m_reso_vOrient;
			break;
		default:
			tl::log_err("Unknown coordinate system selected.");
			return;
	}


	const ublas::matrix<t_real_reso>& reso = *pReso;
	const ublas::vector<t_real_reso>& reso_v = *pReso_v;
	const t_real_reso& reso_s = m_reso_s;
	const ublas::vector<t_real_reso>& _Q_avg = *pQavg;

	try
	{
		// parameters are: x, y, project 1, project 2, remove 1, remove 2
		int iParams[2][4][6] =
		{
			{	// projected
				{ 0, 3, 1, -1, 2, -1 },
				{ 1, 3, 0, -1, 2, -1 },
				{ 2, 3, 0, -1, 1, -1 },
				{ 0, 1, 3, -1, 2, -1 }
			},
			/*{	// projected
				{ 0, 3, 2, 1, -1, -1 },
				{ 1, 3, 2, 0, -1, -1 },
				{ 2, 3, 1, 0, -1, -1 },
				{ 0, 1, 3, 2, -1, -1 }
			},*/
			{	// sliced
				{ 0, 3, -1, -1, 2, 1 },
				{ 1, 3, -1, -1, 2, 0 },
				{ 2, 3, -1, -1, 1, 0 },
				{ 0, 1, -1, -1, 2, 3 }
			}
		};

		static const std::string strDeg = tl::get_spec_char_utf8("deg");

		ublas::vector<t_real_reso> Q_avg = _Q_avg;
		if(m_bCenterOn0)
			Q_avg = ublas::zero_vector<t_real_reso>(Q_avg.size());


		std::vector<std::future<Ellipse2d<t_real_reso>>> tasks_ell_proj, tasks_ell_slice;

		for(unsigned int iEll = 0; iEll < 4; ++iEll)
		{
			// load ellipse configuration from settings
			if(m_pSettings)
			{
				for(unsigned int iSubEll = 0; iSubEll < 2; ++iSubEll)
				{
					std::string elli_name = tl::var_to_str(iEll + 1) + char('a' + iSubEll);

					int x = m_pSettings->value(("reso/ellipse_" + elli_name + "_x").c_str(), -2).toInt();
					if(x > -2)
						iParams[iSubEll][iEll][0] = tl::clamp(x, -1, 3);

					int y = m_pSettings->value(("reso/ellipse_" + elli_name + "_y").c_str(), -2).toInt();
					if(y > -2)
						iParams[iSubEll][iEll][1] = tl::clamp(y, -1, 3);

					int proj1 = m_pSettings->value(("reso/ellipse_" + elli_name + "_proj1").c_str(), -2).toInt();
					if(proj1 > -2)
						iParams[iSubEll][iEll][2] = tl::clamp(proj1, -1, 3);

					int proj2 = m_pSettings->value(("reso/ellipse_" + elli_name + "_proj2").c_str(), -2).toInt();
					if(proj2 > -2)
						iParams[iSubEll][iEll][3] = tl::clamp(proj2, -1, 3);

					int rem1 = m_pSettings->value(("reso/ellipse_" + elli_name + "_rem1").c_str(), -2).toInt();
					if(rem1 > -2)
						iParams[iSubEll][iEll][4] = tl::clamp(rem1, -1, 3);

					int rem2 = m_pSettings->value(("reso/ellipse_" + elli_name + "_rem2").c_str(), -2).toInt();
					if(rem2 > -2)
						iParams[iSubEll][iEll][5] = tl::clamp(rem2, -1, 3);
				}
			}

			const int *iP = iParams[0][iEll];
			const int *iS = iParams[1][iEll];

			std::future<Ellipse2d<t_real_reso>> ell_proj =
				std::async(std::launch::deferred|std::launch::async,
				[=, &reso, &Q_avg]()
				{
					return ::calc_res_ellipse<t_real_reso>(
						reso, reso_v, reso_s, Q_avg,
						iP[0], iP[1], iP[2], iP[3], iP[4], iP[5]);
				});
			std::future<Ellipse2d<t_real_reso>> ell_slice =
				std::async(std::launch::deferred|std::launch::async,
				[=, &reso, &Q_avg]()
				{
					return ::calc_res_ellipse<t_real_reso>(
						reso, reso_v, reso_s, Q_avg,
						iS[0], iS[1], iS[2], iS[3], iS[4], iS[5]);
				});

			tasks_ell_proj.push_back(std::move(ell_proj));
			tasks_ell_slice.push_back(std::move(ell_slice));


			// MC neutrons
			if(pvecMC)
			{
				m_vecMCXCurvePoints[iEll].resize(pvecMC->size());
				m_vecMCYCurvePoints[iEll].resize(pvecMC->size());

				for(std::size_t iMC = 0; iMC < pvecMC->size(); ++iMC)
				{
					const ublas::vector<t_real_reso>& vecMC = (*pvecMC)[iMC];
					m_vecMCXCurvePoints[iEll][iMC] = vecMC[iP[0]];
					m_vecMCYCurvePoints[iEll][iMC] = vecMC[iP[1]];

					if(m_bCenterOn0)
					{
						m_vecMCXCurvePoints[iEll][iMC] -= _Q_avg[iP[0]];
						m_vecMCYCurvePoints[iEll][iMC] -= _Q_avg[iP[1]];
					}
				}
			}
			else
			{
				m_vecMCXCurvePoints[iEll].clear();
				m_vecMCYCurvePoints[iEll].clear();
			}
		}

		for(unsigned int iEll = 0; iEll < 4; ++iEll)
		{
			m_elliProj[iEll] = tasks_ell_proj[iEll].get();
			m_elliSlice[iEll] = tasks_ell_slice[iEll].get();

			std::vector<t_real_reso>& vecXProj = m_vecXCurvePoints[iEll*2 + 0];
			std::vector<t_real_reso>& vecYProj = m_vecYCurvePoints[iEll*2 + 0];
			std::vector<t_real_reso>& vecXSlice = m_vecXCurvePoints[iEll*2 + 1];
			std::vector<t_real_reso>& vecYSlice = m_vecYCurvePoints[iEll*2 + 1];
			std::vector<t_real_reso>& vecXMC = m_vecMCXCurvePoints[iEll];
			std::vector<t_real_reso>& vecYMC = m_vecMCYCurvePoints[iEll];

			t_real_reso dBBProj[4], dBBSlice[4];
			m_elliProj[iEll].GetCurvePoints(vecXProj, vecYProj, GFX_NUM_POINTS, dBBProj);
			m_elliSlice[iEll].GetCurvePoints(vecXSlice, vecYSlice, GFX_NUM_POINTS, dBBSlice);

			set_qwt_data<t_real_reso>()(*m_vecplotwrap[iEll], vecXProj, vecYProj, 1, false);
			set_qwt_data<t_real_reso>()(*m_vecplotwrap[iEll], vecXSlice, vecYSlice, 2, false);
			set_qwt_data<t_real_reso>()(*m_vecplotwrap[iEll], vecXMC, vecYMC, 0, false);


			std::ostringstream ostrSlope;
			ostrSlope.precision(g_iPrecGfx);
			ostrSlope << "Projected ellipse (green, HWHM contour):\n";
			ostrSlope << "\tSlope: " << m_elliProj[iEll].slope << "\n";
			ostrSlope << "\tAngle: " << tl::r2d(m_elliProj[iEll].phi) << strDeg << "\n";
			ostrSlope << "\tArea " << m_elliProj[iEll].area << "\n";
			ostrSlope << "Sliced ellipse (blue, HWHM contour):\n";
			ostrSlope << "\tSlope: " << m_elliSlice[iEll].slope << "\n";
			ostrSlope << "\tAngle: " << tl::r2d(m_elliSlice[iEll].phi) << strDeg << "\n";
			ostrSlope << "\tArea " << m_elliSlice[iEll].area;
			m_vecplotwrap[iEll]->GetPlot()->setToolTip(QString::fromUtf8(ostrSlope.str().c_str()));

			const std::string& strLabX = ellipse_labels(iParams[0][iEll][0], coord, m_bCenterOn0);
			const std::string& strLabY = ellipse_labels(iParams[0][iEll][1], coord, m_bCenterOn0);
			m_vecplotwrap[iEll]->GetPlot()->setAxisTitle(QwtPlot::xBottom, strLabX.c_str());
			m_vecplotwrap[iEll]->GetPlot()->setAxisTitle(QwtPlot::yLeft, strLabY.c_str());

			m_vecplotwrap[iEll]->GetPlot()->replot();

			QRectF rect;
			rect.setLeft(std::min(dBBProj[0], dBBSlice[0]));
			rect.setRight(std::max(dBBProj[1], dBBSlice[1]));
			rect.setTop(std::max(dBBProj[2], dBBSlice[2]));
			rect.setBottom(std::min(dBBProj[3], dBBSlice[3]));
			if(m_vecplotwrap[iEll]->GetZoomer())
				m_vecplotwrap[iEll]->GetZoomer()->setZoomBase(rect);

			switch(m_algo)
			{
				case ResoAlgo::CN: SetTitle("Cooper-Nathans Algorithm (Pointlike TAS)"); break;
				case ResoAlgo::POP_CN: SetTitle("Popovici Algorithm (Pointlike TAS)"); break;
				case ResoAlgo::POP: SetTitle("Popovici Algorithm (TAS)"); break;
				case ResoAlgo::ECK: SetTitle("Eckold-Sobolev Algorithm (TAS)"); break;
				case ResoAlgo::ECK_EXT: SetTitle("Extended Eckold-Sobolev Algorithm (TAS)"); break;
				case ResoAlgo::VIO: SetTitle("Violini Algorithm (TOF)"); break;
				case ResoAlgo::SIMPLE: SetTitle("Simple Algorithm"); break;
				case ResoAlgo::MC: SetTitle("MC"); break;
				default: SetTitle("Unknown Resolution Algorithm"); break;
			}
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err("Cannot calculate ellipses.");
		SetTitle("Error");
	}
}


void EllipseDlg::SetCenterOn0(bool bCenter)
{
	m_bCenterOn0 = bCenter;
	Calc();
}


void EllipseDlg::SetParams(const EllipseDlgParams& params)
{
	m_params = params;

	static const ublas::matrix<t_real_reso> mat0 = ublas::zero_matrix<t_real_reso>(4, 4);
	static const ublas::vector<t_real_reso> vec0 = ublas::zero_vector<t_real_reso>(4);

	m_reso = params.reso ? *params.reso : mat0;
	m_reso_v = params.reso_v ? *params.reso_v : vec0;
	m_reso_s = params.reso_s;
	m_Q_avg = params.Q_avg ? *params.Q_avg : vec0;
	m_resoHKL = params.resoHKL ? *params.resoHKL : mat0;
	m_reso_vHKL = params.reso_vHKL ? *params.reso_vHKL : vec0;
	m_Q_avgHKL = params.Q_avgHKL ? *params.Q_avgHKL : vec0;
	m_resoOrient = params.resoOrient ? *params.resoOrient : mat0;
	m_reso_vOrient = params.reso_vOrient ? *params.reso_vOrient : vec0;
	m_Q_avgOrient = params.Q_avgOrient ? *params.Q_avgOrient : vec0;
	m_algo = params.algo;

	Calc();
}


/**
 * export ellipses as gnuplot script
 */
void EllipseDlg::SaveEllipses()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = ".";
	if(m_pSettings)
		strDirLast = m_pSettings->value("reso/last_dir_ellipse", ".").toString();
	QString qstrFile = QFileDialog::getSaveFileName(this,
		"Save resolution ellipses", strDirLast, "Gnuplot files (*.gpl)",
		nullptr, fileopt);

	if(qstrFile == "")
		return;

	std::string strExport = R"RAWSTR(#!gnuplot --persist
#
# Resolution ellipse plot.
# Created with Takin %%TAKIN_VER%% (https://dx.doi.org/10.5281/zenodo.4117437).
# Date: %%DATE%%.
#

# -----------------------------------------------------------------------------
# output to a file
# -----------------------------------------------------------------------------
#set term pdf color enhanced font "NimbusSans-Regular, 54" size 20, 15
#set output "reso.pdf"
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# resolution ellipse
# -----------------------------------------------------------------------------
ellipse_x(t, hwhm_x, hwhm_y, angle, offs_x) = \
	hwhm_x*cos(2*pi*t)*cos(angle) - hwhm_y*sin(2*pi*t)*sin(angle) + offs_x
ellipse_y(t, hwhm_x, hwhm_y, angle, offs_y) = \
	hwhm_x*cos(2*pi*t)*sin(angle) + hwhm_y*sin(2*pi*t)*cos(angle) + offs_y
# -----------------------------------------------------------------------------


set parametric
set multiplot layout 2, 2 #margins 0.15, 0.95, 0.15, 0.95 spacing 0.15, 0.15
set border linewidth 2
set trange [ 0 : 1 ]


# -----------------------------------------------------------------------------
# ellipse parameters
# -----------------------------------------------------------------------------
%%PARAMETERS%%
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plots
# -----------------------------------------------------------------------------
# ellipse colours
col_proj  = "#ff0000"
col_slice = "#0000ff"

linew = 2


%%PLOTS%%
# -----------------------------------------------------------------------------
)RAWSTR";

	std::string strFile = qstrFile.toStdString();
	std::string strDir = tl::get_dir(strFile);
	if(tl::get_fileext(strFile, 1) != "gpl")
		strFile += ".gpl";

	std::string labels_x[] =
	{
		"set xlabel \"Q_{||} - <Q> (Å⁻¹)\"",
		"set xlabel \"Q_{⟂} - <Q> (Å⁻¹)\"",
		"set xlabel \"Q_z - <Q> (Å⁻¹)\"",
		"set xlabel \"Q_{||} - <Q> (Å⁻¹)\"",
	};

	std::string labels_y[] =
	{
		"set ylabel \"E (meV)\"",
		"set ylabel \"E (meV)\"",
		"set ylabel \"E (meV)\"",
		"set ylabel \"Q_{⟂} - <Q> (Å⁻¹)\"",
	};

	std::ostringstream ostrParams, ostrPlots;
	ostrParams.precision(g_iPrec);
	ostrPlots.precision(g_iPrec);
	std::size_t num_ellis = std::min(m_elliProj.size(), m_elliSlice.size());
	for(std::size_t elli_idx = 0; elli_idx < num_ellis; ++elli_idx)
	{
		ostrParams << "hwhm_proj_" << (elli_idx+1) << "_x  = " << m_elliProj[elli_idx].x_hwhm << "\n";
		ostrParams << "hwhm_proj_" << (elli_idx+1) << "_y  = " << m_elliProj[elli_idx].y_hwhm << "\n";
		ostrParams << "angle_proj_" << (elli_idx+1) << "   = " << m_elliProj[elli_idx].phi << "\n";
		ostrParams << "offs_proj_" << (elli_idx+1) << "_x  = " << m_elliProj[elli_idx].x_offs << "\n";
		ostrParams << "offs_proj_" << (elli_idx+1) << "_y  = " << m_elliProj[elli_idx].y_offs << "\n\n";

		ostrParams << "hwhm_slice_" << (elli_idx+1) << "_x = " << m_elliSlice[elli_idx].x_hwhm << "\n";
		ostrParams << "hwhm_slice_" << (elli_idx+1) << "_y = " << m_elliSlice[elli_idx].y_hwhm << "\n";
		ostrParams << "angle_slice_" << (elli_idx+1) << "  = " << m_elliSlice[elli_idx].phi << "\n";
		ostrParams << "offs_slice_" << (elli_idx+1) << "_x = " << m_elliSlice[elli_idx].x_offs << "\n";
		ostrParams << "offs_slice_" << (elli_idx+1) << "_y = " << m_elliSlice[elli_idx].y_offs << "\n";


		ostrPlots << labels_x[elli_idx] << "\n";
		ostrPlots << labels_y[elli_idx] << "\n\n";
		ostrPlots << "plot \\\n";
		ostrPlots << "\tellipse_x(t, hwhm_proj_" << (elli_idx+1)
			<< "_x, hwhm_proj_" << (elli_idx+1)
			<< "_y, angle_proj_" << (elli_idx+1)
			<< ", offs_proj_" << (elli_idx+1)
			<< "_x),\\\n";
		ostrPlots << "\tellipse_y(t, hwhm_proj_" << (elli_idx+1)
			<< "_x, hwhm_proj_" << (elli_idx+1)
			<< "_y, angle_proj_" << (elli_idx+1) << ", offs_proj_" << (elli_idx+1)
			<<  "_y) \\\n";
		ostrPlots << "\t\tlinewidth linew linecolor rgb col_proj notitle, \\\n";
		ostrPlots << "\tellipse_x(t, hwhm_slice_" << (elli_idx+1)
			<< "_x, hwhm_slice_" <<(elli_idx+1)
			<< "_y, angle_slice_" << (elli_idx+1)
			<< ", offs_slice_" << (elli_idx+1)
			<< "_x),\\\n";
		ostrPlots << "\tellipse_y(t, hwhm_slice_" << (elli_idx+1)
			<< "_x, hwhm_slice_" << (elli_idx+1)
			<< "_y, angle_slice_" << (elli_idx+1)
			<< ", offs_slice_" << (elli_idx+1)
			<< "_y) \\\n";
		ostrPlots << "\t\tlinewidth linew linecolor rgb col_slice notitle\n";


		if(elli_idx < num_ellis - 1)
		{
			ostrParams << "\n\n";
			ostrPlots << "\n\n";
		}
	}

	tl::find_all_and_replace<std::string>(strExport, "%%PARAMETERS%%", ostrParams.str());
	tl::find_all_and_replace<std::string>(strExport, "%%PLOTS%%", ostrPlots.str());
	tl::find_all_and_replace<std::string>(strExport, "%%TAKIN_VER%%", TAKIN_VER);
	tl::find_all_and_replace<std::string>(strExport, "%%DATE%%",
		tl::epoch_to_str<t_real_reso>(tl::epoch<t_real_reso>(), "%b %d, %Y; %H:%M:%S (%Z)"));

	std::ofstream ofstr(strFile);
	ofstr << strExport;
	ofstr.flush();

	if(m_pSettings)
		m_pSettings->setValue("reso/last_dir_ellipse", QString(strDir.c_str()));
}


void EllipseDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("reso/ellipse_geo", saveGeometry());

	QDialog::accept();
}


void EllipseDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


void EllipseDlg::closeEvent(QCloseEvent *pEvt)
{
	QDialog::closeEvent(pEvt);
}


#include "moc_EllipseDlg.cpp"
