/**
 * Dynamic Plane Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013, jan-2015
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

#include "DynPlaneDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/phys/neutrons.h"
#include <boost/units/io.hpp>


using t_real = t_real_glob;
static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();


DynPlaneDlg::DynPlaneDlg(QWidget* pParent, QSettings *pSettings)
		: QDialog(pParent), m_pSettings(pSettings)
{
	this->setupUi(this);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}


	m_plotwrap.reset(new QwtPlotWrapper(plot));
	m_plotwrap->GetCurve(0)->setTitle("Kinematic Plane");
	m_plotwrap->GetPlot()->setAxisTitle(QwtPlot::xBottom, "Q (1/A)");
	m_plotwrap->GetPlot()->setAxisTitle(QwtPlot::yLeft, "E (meV)");
	if(m_plotwrap->HasTrackerSignal())
		connect(m_plotwrap->GetPicker(), &QwtPlotPicker::moved, this, &DynPlaneDlg::cursorMoved);


	QObject::connect(comboFixedE, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &DynPlaneDlg::FixedKiKfToggled);

	std::vector<QDoubleSpinBox*> vecSpinBoxes = {spinEiEf, spinMinQ, spinMaxQ, spinAngle};
	for(QDoubleSpinBox* pSpin : vecSpinBoxes)
		QObject::connect(pSpin, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &DynPlaneDlg::Calc);
	QObject::connect(btnSync, &QPushButton::toggled, this, &DynPlaneDlg::Calc);


	if(m_pSettings)
	{
		if(m_pSettings->contains("dyn_plane/geo"))
			restoreGeometry(m_pSettings->value("dyn_plane/geo").toByteArray());

		if(m_pSettings->contains("dyn_plane/kfix"))
			comboFixedE->setCurrentIndex(m_pSettings->value("dyn_plane/kfix").toInt());
	}


	Calc();
}


DynPlaneDlg::~DynPlaneDlg()
{}


/**
 * shows accessible E ranges for cursor Q position
 */
void DynPlaneDlg::cursorMoved(const QPointF& pt)
{
	try
	{
		const bool bFixedKi = (comboFixedE->currentIndex() == 0);
		const tl::t_angle_si<t_real> twotheta = tl::d2r(t_real(spinAngle->value())) * rads;
		const tl::t_energy_si<t_real> EiEf = t_real(spinEiEf->value()) * meV;

		const tl::t_wavenumber_si<t_real> Q = t_real(pt.x()) / angs;
		tl::t_energy_si<t_real> dE0 = tl::kinematic_plane(bFixedKi, false, EiEf, Q, twotheta);
		tl::t_energy_si<t_real> dE1 = tl::kinematic_plane(bFixedKi, true, EiEf, Q, twotheta);

		t_real _dQ = Q * angs;
		t_real _dE0 = dE0 / meV;
		t_real _dE1 = dE1 / meV;

		std::ostringstream ostr;

		// valid accessible energy ranges?
		if(!std::isnan(_dQ) && !std::isnan(_dE0) && !std::isnan(_dE1)
			&& !std::isinf(_dQ) && !std::isinf(_dE0) && !std::isinf(_dE1))
		{
			if(_dE1 < _dE0)
				std::swap(_dE0, _dE1);

			std::string strQ = tl::var_to_str(_dQ, g_iPrecGfx);
			std::string strE0 = tl::var_to_str(_dE0, g_iPrecGfx);
			std::string strE1 = tl::var_to_str(_dE1, g_iPrecGfx);

			ostr << "Q = " << strQ << " /A, E = [ " << strE0 << ", " << strE1 << " ] meV";
		}
		// otherwise just show curser coordinates
		else
		{
			std::string strQ = tl::var_to_str(pt.x(), g_iPrecGfx);
			std::string strE = tl::var_to_str(pt.y(), g_iPrecGfx);

			ostr << "Q = " << strQ << " /A, E = " << strE << " meV [not accessible]";
		}

		this->labelStatus->setText(ostr.str().c_str());
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


/**
 * calculates and plots the accessible ranges
 */
void DynPlaneDlg::Calc()
{
	try
	{
		const bool bFixedKi = (comboFixedE->currentIndex()==0);

		if(btnSync->isChecked())
		{
			spinEiEf->setValue(bFixedKi ? m_dEi : m_dEf);
			spinAngle->setValue(tl::r2d(m_d2Theta));
			spinMinQ->setValue(0.);
			spinMaxQ->setValue(m_dQ*2.);
		}

		const t_real dMinQ = spinMinQ->value();
		const t_real dMaxQ = spinMaxQ->value();
		const t_real dAngle = tl::d2r(t_real(spinAngle->value()));
		const tl::t_energy_si<t_real> EiEf = t_real(spinEiEf->value()) * meV;


		std::vector<t_real> vecQ[2], vecE[2];
		vecQ[0].reserve(GFX_NUM_POINTS); vecE[0].reserve(GFX_NUM_POINTS);
		vecQ[1].reserve(GFX_NUM_POINTS); vecE[1].reserve(GFX_NUM_POINTS);

		tl::t_angle_si<t_real> twotheta = dAngle * rads;

		for(std::size_t iPt = 0; iPt < GFX_NUM_POINTS; ++iPt)
		{
			for(unsigned char iSign = 0; iSign <= 1; ++iSign)
			{
				tl::t_wavenumber_si<t_real> Q = (dMinQ + (dMaxQ - dMinQ)/t_real(GFX_NUM_POINTS)*t_real(iPt)) / angs;
				tl::t_energy_si<t_real> dE = tl::kinematic_plane(bFixedKi, iSign, EiEf, Q, twotheta);

				t_real _dQ = Q * angs;
				t_real _dE = dE / meV;

				if(!std::isnan(_dQ) && !std::isnan(_dE) && !std::isinf(_dQ) && !std::isinf(_dE))
				{
					vecQ[iSign].push_back(Q * angs);
					vecE[iSign].push_back(dE / meV);
				}
			}
		}

		m_vecQ.clear();
		m_vecE.clear();

		m_vecQ.insert(m_vecQ.end(), vecQ[0].rbegin(), vecQ[0].rend());
		m_vecE.insert(m_vecE.end(), vecE[0].rbegin(), vecE[0].rend());

		m_vecQ.insert(m_vecQ.end(), vecQ[1].begin(), vecQ[1].end());
		m_vecE.insert(m_vecE.end(), vecE[1].begin(), vecE[1].end());

		set_qwt_data<t_real>()(*m_plotwrap, m_vecQ, m_vecE);
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


/**
 * received new set of reciprocal coordinates
 */
void DynPlaneDlg::RecipParamsChanged(const RecipParams& params)
{
	m_d2Theta = params.d2Theta;
	m_dEi = tl::k2E(params.dki/angs)/meV;
	m_dEf = tl::k2E(params.dkf/angs)/meV;
	m_dQ = params.dQ;

	Calc();
}


void DynPlaneDlg::FixedKiKfToggled()
{
	if(comboFixedE->currentIndex() == 0)
		labelFixedKiKf->setText("E_i (meV):");
	else
		labelFixedKiKf->setText("E_f (meV):");

	Calc();
}



void DynPlaneDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


void DynPlaneDlg::accept()
{
	if(m_pSettings)
	{
		m_pSettings->setValue("dyn_plane/geo", saveGeometry());
		m_pSettings->setValue("dyn_plane/kfix", comboFixedE->currentIndex());
	}

	QDialog::accept();
}


#include "moc_DynPlaneDlg.cpp"
