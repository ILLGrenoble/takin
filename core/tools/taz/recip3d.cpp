/**
 * Scattering Triangle Tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014 - 2017
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

#include "recip3d.h"
#include "tlibs/math/geo.h"
#include <QGridLayout>


#define DEF_PEAK_SIZE 0.04
#define MIN_PEAK_SIZE 0.015
#define MAX_PEAK_SIZE 0.15


using t_real = t_real_glob;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;


Recip3DDlg::Recip3DDlg(QWidget* pParent, QSettings *pSettings)
	: QDialog(pParent, Qt::Tool), m_pSettings(pSettings),
	m_pStatus(new QStatusBar(this)),
	m_pPlot(new PlotGl(this, pSettings, 0.25))
{
	m_pPlot->SetEnabled(false);

	setWindowTitle("Reciprocal Space");
	m_pStatus->setSizeGripEnabled(1);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("recip3d/geo"))
			restoreGeometry(m_pSettings->value("recip3d/geo").toByteArray());
		else
			resize(800, 600);
	}

	m_pPlot->SetPrec(g_iPrecGfx);
	m_pPlot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

	QGridLayout *gridLayout = new QGridLayout(this);
	gridLayout->setContentsMargins(4, 4, 4, 4);
	gridLayout->addWidget(m_pPlot, 0, 0, 1, 1);
	gridLayout->addWidget(m_pStatus, 1, 0, 1, 1);

	m_pPlot->AddHoverSlot([this](const PlotObjGl* pObj)
	{
		std::string strStatus;
		if(pObj)
			strStatus = pObj->strLabel;
		tl::find_all_and_replace<std::string>(strStatus, "\n", ", ");

		if(strStatus.length())
			m_pStatus->showMessage(strStatus.c_str());
		else
			m_pStatus->clearMessage();
	});

	m_vecQ_rlu = m_vecQ = tl::ublas::zero_vector<t_real>(3);

	m_pPlot->SetLabels("a", "b", "c");
	m_pPlot->SetEnabled(true);
}


Recip3DDlg::~Recip3DDlg()
{
	if(m_pPlot)
	{
		m_pPlot->SetEnabled(false);
		delete m_pPlot;
		m_pPlot = 0;
	}
}


struct Peak3d
{
	std::string strName;
	t_vec vecPeak;
	t_real dF = -1.;
	const std::vector<t_real>* pvecColor = nullptr;
};

void Recip3DDlg::CalcPeaks(const xtl::LatticeCommon<t_real_glob>& recipcommon)
{
	const tl::Lattice<t_real>& recip = recipcommon.recip;
	const xtl::SpaceGroup<t_real>* pSpaceGroup = recipcommon.pSpaceGroup;
	const tl::Plane<t_real>& plane = recipcommon.plane;
	const t_vec &vec0 = plane.GetDir0(),
		&vec1 = plane.GetDir1(),
		&vecNorm = plane.GetNorm();

	static const std::vector<t_real> vecColPeak = { 0., 0., 1., 0.7 };
	static const std::vector<t_real> vecColPeakPlane = { 1., 0., 0., 0.7 };
	static const std::vector<t_real> vecCol000 = { 0., 1., 0., 0.7 };
	static const std::vector<t_real> vecColCurQ = { 0., 0., 0., 1. };
	static const std::vector<t_real> vecColScatPlane = { 0., 0., 0.75, 0.15 };
	static const std::vector<t_real> vecColLine = { 0., 0., 0., 1. };

	bool bShowScatPlane = 1;
	bool bShowCurQ = 1;

	const t_real dLimMax = std::numeric_limits<t_real>::max();

	std::vector<t_real> vecMin = {dLimMax, dLimMax, dLimMax},
		vecMax = {-dLimMax, -dLimMax, -dLimMax};
	t_real dMinF = std::numeric_limits<t_real>::max(), dMaxF = -1.;


	std::vector<Peak3d> vecPeaks;

	for(t_real h = -m_dMaxPeaks; h <= m_dMaxPeaks; h += 1.)
	for(t_real k = -m_dMaxPeaks; k <= m_dMaxPeaks; k += 1.)
	for(t_real l = -m_dMaxPeaks; l <= m_dMaxPeaks; l += 1.)
	{
		int ih = int(h), ik = int(k), il = int(l);
		if(pSpaceGroup)
		{
			// if reflection is not allowed, skip it
			if(!pSpaceGroup->HasReflection(ih, ik, il))
				continue;
		}

		Peak3d peak;
		peak.vecPeak = recip.GetPos(h,k,l);
		for(int i=0; i<3; ++i)
		{
			vecMin[i] = std::min(peak.vecPeak[i], vecMin[i]);
			vecMax[i] = std::max(peak.vecPeak[i], vecMax[i]);
		}


		// -------------------------------------------------------------
		// get structure factor
		std::string strStructfact;
		t_real dFsq = -1.;

		if(pSpaceGroup && recipcommon.CanCalcStructFact())
		{
			std::tie(std::ignore, peak.dF, dFsq) =
				recipcommon.GetStructFact(peak.vecPeak);

			tl::set_eps_0(dFsq, g_dEpsGfx);
			tl::set_eps_0(peak.dF, g_dEpsGfx);

			dMinF = std::min(peak.dF, dMinF);
			dMaxF = std::max(peak.dF, dMaxF);

			std::ostringstream ostrStructfact;
			ostrStructfact.precision(g_iPrecGfx);
			if(g_bShowFsq)
				ostrStructfact << "\nS = " << dFsq;
			else
				ostrStructfact << "\nF = " << peak.dF;
			strStructfact = ostrStructfact.str();
		}
		// -------------------------------------------------------------


		// -------------------------------------------------------------
		// is the reflection in the scattering plane?
		t_real dDist = 0.;
		t_vec vecDropped = plane.GetDroppedPerp(peak.vecPeak, &dDist);

		peak.pvecColor = &vecColPeak;
		if(tl::float_equal<t_real>(dDist, 0., m_dPlaneDistTolerance))
		{
			// (000) peak?
			if(tl::float_equal<t_real>(h, 0., g_dEps) &&
				tl::float_equal<t_real>(k, 0., g_dEps) &&
				tl::float_equal<t_real>(l, 0., g_dEps))
			{
				peak.pvecColor = &vecCol000;
			}
			else
			{
				peak.pvecColor = &vecColPeakPlane;
			}
		}
		// -------------------------------------------------------------


		std::ostringstream ostrLab;
		ostrLab << "(" << ih << " " << ik << " " << il << ")"
			<< strStructfact;
		peak.strName = ostrLab.str();

		vecPeaks.emplace_back(std::move(peak));
	}



	std::size_t iObjCnt = vecPeaks.size();
	if(bShowScatPlane)
		iObjCnt += 2;
	if(bShowCurQ)
		++iObjCnt;

	m_pPlot->SetEnabled(false);
	m_pPlot->clear();
	m_pPlot->SetObjectCount(iObjCnt);

	std::size_t iCurObjIdx = 0;


	// plot all allowed reflections
	for(const Peak3d& peak : vecPeaks)
	{
		t_real dFRad = DEF_PEAK_SIZE;

		// valid structure factors
		if(peak.dF >= 0. && !tl::float_equal(dMinF, dMaxF, g_dEpsGfx))
		{
			t_real dFScale = (peak.dF-dMinF) / (dMaxF-dMinF);
			dFRad = tl::lerp(MIN_PEAK_SIZE, MAX_PEAK_SIZE, dFScale);
		}

		m_pPlot->PlotSphere(peak.vecPeak, dFRad, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, *peak.pvecColor);
		m_pPlot->SetObjectLabel(iCurObjIdx, peak.strName);

		++iCurObjIdx;
	}


	// scattering plane
	if(bShowScatPlane)
	{
		std::vector<t_vec> vecVerts =
		{
			-vec0*m_dMaxPeaks - vec1*m_dMaxPeaks,
			-vec0*m_dMaxPeaks + vec1*m_dMaxPeaks,
			 vec0*m_dMaxPeaks + vec1*m_dMaxPeaks,
			 vec0*m_dMaxPeaks - vec1*m_dMaxPeaks
		};

		m_pPlot->PlotPoly(vecVerts, vecNorm, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColScatPlane);
		m_pPlot->SetObjectCull(iCurObjIdx, 0);
		++iCurObjIdx;

		m_pPlot->PlotLines(vecVerts, 2., iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColLine);
		++iCurObjIdx;
	}


	// current Q position
	if(bShowCurQ)
	{
		m_pPlot->PlotSphere(m_vecQ, 0.04, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColCurQ);
		m_pPlot->SetObjectAnimation(iCurObjIdx, 1);

		// current Q label
		std::ostringstream ostrTip;
		ostrTip.precision(g_iPrecGfx);
		ostrTip << "Q = (" << m_vecQ_rlu[0] << ", " << m_vecQ_rlu[1] << ", " << m_vecQ_rlu[2] << ") rlu";
		ostrTip << "\nQ = (" << m_vecQ[0] << ", " << m_vecQ[1] << ", " << m_vecQ[2] << ") 1/A";
		m_pPlot->SetObjectLabel(iCurObjIdx, ostrTip.str());

		m_iQIdx = iCurObjIdx;
		++iCurObjIdx;
	}

	m_pPlot->SetMinMax(vecMin, vecMax);
	m_pPlot->SetEnabled(true);
}


/**
 * scattering triangle changed
 */
void Recip3DDlg::RecipParamsChanged(const RecipParams& recip)
{
	m_vecQ = tl::make_vec({ -recip.Q[0], -recip.Q[1], -recip.Q[2] });
	m_vecQ_rlu = tl::make_vec({ -recip.Q_rlu[0], -recip.Q_rlu[1], -recip.Q_rlu[2] });

	// if a 3d object is already assigned, update it
	if(m_pPlot && m_iQIdx)
	{
		m_pPlot->SetEnabled(false);
		m_pPlot->PlotSphere(m_vecQ, 0.04, *m_iQIdx);

		// current Q label
		std::ostringstream ostrTip;
		ostrTip.precision(g_iPrecGfx);
		ostrTip << "Q = (" << m_vecQ_rlu[0] << ", " << m_vecQ_rlu[1] << ", " << m_vecQ_rlu[2] << ") rlu";
		ostrTip << "\nQ = (" << m_vecQ[0] << ", " << m_vecQ[1] << ", " << m_vecQ[2] << ") 1/A";
		m_pPlot->SetObjectLabel(*m_iQIdx, ostrTip.str());

		m_pPlot->SetEnabled(true);
	}
}


void Recip3DDlg::keyPressEvent(QKeyEvent* pEvt)
{
	if(m_pPlot)
		m_pPlot->keyPressEvent(pEvt);

	QDialog::keyPressEvent(pEvt);
}


void Recip3DDlg::closeEvent(QCloseEvent *pEvt)
{
	if(m_pSettings)
		m_pSettings->setValue("recip3d/geo", saveGeometry());
	QDialog::closeEvent(pEvt);
}


void Recip3DDlg::hideEvent(QHideEvent *pEvt)
{
	if(m_pPlot) m_pPlot->SetEnabled(false);
	QDialog::hideEvent(pEvt);
}


void Recip3DDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
	if(m_pPlot) m_pPlot->SetEnabled(true);
}


#include "moc_recip3d.cpp"
