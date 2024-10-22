/**
 * 3d Brillouin zone drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2017
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

// hack because of include file ordering
#ifdef __MINGW32__
	#include <boost/asio.hpp>
#endif

#include "bz3d.h"
#include <QGridLayout>

#include "tlibs/math/linalg.h"
#include "tlibs/string/string.h"


using t_real = t_real_glob;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;



BZ3DDlg::BZ3DDlg(QWidget* pParent, QSettings *pSettings)
	: QDialog(pParent, Qt::Tool), m_pSettings(pSettings),
	m_pStatus(new QStatusBar(this)),
	m_pPlot(new PlotGl(this, pSettings, 0.25))
{
	m_pPlot->SetEnabled(false);

	setWindowTitle("Reciprocal Space / Brillouin Zone");
	m_pStatus->setSizeGripEnabled(1);
	m_pStatus->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("bz3d/geo"))
			restoreGeometry(m_pSettings->value("bz3d/geo").toByteArray());
		else
			resize(800, 600);
	}

	m_pPlot->SetPrec(g_iPrecGfx);
	m_pPlot->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

	QGridLayout *gridLayout = new QGridLayout(this);
	gridLayout->setContentsMargins(4, 4, 4, 4);
	gridLayout->addWidget(m_pPlot.get(), 0, 0, 1, 1);
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

	m_vecq_rlu = m_vecq = tl::ublas::zero_vector<t_real>(3);
	m_pPlot->SetLabels("x (1/A)", "y (1/A)", "z (1/A)");
	m_pPlot->SetDrawMinMax(0);
	m_pPlot->SetEnabled(true);
}


BZ3DDlg::~BZ3DDlg()
{
	if(m_pPlot)
		m_pPlot->SetEnabled(false);
}


// ----------------------------------------------------------------------------


/**
 * draw 3d BZ
 * assumes to be centred around (0,0,0): see CALC_BZ_AROUND_ZERO in scattering_triangle.cpp
 */
void BZ3DDlg::RenderBZ(const tl::Brillouin3D<t_real_glob>& bz,
	const xtl::LatticeCommon<t_real_glob>& lattice,
	const std::vector<ublas::vector<t_real_glob>>* pScatPlaneVerts,
	const std::vector<ublas::vector<t_real_glob>>* pvecSymmPts)
{
	if(!bz.IsValid() || !m_pPlot)
		return;

	static const std::vector<t_real> vecColVertices = { 0., 0., 0., 0.85 };
	static const std::vector<t_real> vecColSymmVerts = { 1., 0., 0., 0.85 };
	static const std::vector<t_real> vecColPolys = { 0., 0., 1., 0.85 };
	static const std::vector<t_real> vecColEdges = { 0., 0., 0., 1. };
	static const std::vector<t_real> vecColScatPlane = { 1., 1., 0., 1. };
	static const std::vector<t_real> vecColCurq = { 0., 0., 0., 1. };

	m_pPlot->SetEnabled(false);
	m_pPlot->clear();

	t_mat matBinv = tl::transpose(lattice.matA) / (tl::get_pi<t_real>()*t_real(2));


	const bool bShowVerts = false;
	const bool bShowSymmPts = true;
	const bool bShowCurq = true;

	// all objects: polys + edges
	std::size_t iNumObjs =  2*bz.GetPolys().size();
	if(bShowVerts)
		iNumObjs += bz.GetVertices().size();
	if(pScatPlaneVerts)
		++iNumObjs;
	if(bShowSymmPts && pvecSymmPts)
		iNumObjs += pvecSymmPts->size();
	if(bShowCurq)
		++iNumObjs;
	m_pPlot->SetObjectCount(iNumObjs);


	// minimum and maximum coordinates
	const t_real dLimMax = std::numeric_limits<t_real>::max();
	std::vector<t_real> vecMin = { dLimMax, dLimMax, dLimMax },
		vecMax = { -dLimMax, -dLimMax, -dLimMax };

	std::size_t iCurObjIdx = 0;

	// render vertices
	for(const t_vec& vec : bz.GetVertices())
	{
		t_vec vecRLU = ublas::prod(matBinv, vec);
		tl::set_eps_0(vecRLU);

		if(bShowVerts)
		{
			m_pPlot->PlotSphere(vec, 0.025/**0.1*g_dFontSize*/, iCurObjIdx);
			m_pPlot->SetObjectColor(iCurObjIdx, vecColVertices);

			// label
			std::ostringstream ostrTip;
			ostrTip.precision(g_iPrecGfx);
			ostrTip << "(" << vecRLU[0] << ", " << vecRLU[1] << ", " << vecRLU[2] << ") rlu";
			ostrTip << "\n(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ") 1/A";
			m_pPlot->SetObjectLabel(iCurObjIdx, ostrTip.str());

			++iCurObjIdx;
		}

		// min & max
		for(unsigned int i=0; i<3; ++i)
		{
			vecMin[i] = std::min(vec[i], vecMin[i]);
			vecMax[i] = std::max(vec[i], vecMax[i]);
		}
	}

	// render polygons
	for(std::size_t iPoly=0; iPoly<bz.GetPolys().size(); ++iPoly)
	{
		m_pPlot->PlotPoly(bz.GetPolys()[iPoly], bz.GetPlanes()[iPoly].GetNorm(), iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColPolys);

		++iCurObjIdx;
	}

	// render edges
	for(const std::vector<t_vec>& vecPoly : bz.GetPolys())
	{
		m_pPlot->PlotLines(vecPoly, 2.*0.1*g_dFontSize, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColEdges);

		++iCurObjIdx;
	}

	// render scattering plane if available
	if(pScatPlaneVerts)
	{
		m_pPlot->PlotLines(*pScatPlaneVerts, 4.*0.1*g_dFontSize, iCurObjIdx);
		//m_pPlot->PlotPoly(*pScatPlaneVerts, tl::make_vec({0.,0.,1.}), iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColScatPlane);

		++iCurObjIdx;
	}

	// render points of high symmetry if available
	if(bShowSymmPts && pvecSymmPts)
	{
		for(const t_vec& vec : *pvecSymmPts)
		{
			t_vec vecRLU = ublas::prod(matBinv, vec);
			tl::set_eps_0(vecRLU);

			m_pPlot->PlotSphere(vec, 0.025/**0.1*g_dFontSize*/, iCurObjIdx);
			m_pPlot->SetObjectColor(iCurObjIdx, vecColSymmVerts);

			// label
			std::ostringstream ostrTip;
			ostrTip.precision(g_iPrecGfx);
			ostrTip << "(" << vecRLU[0] << ", " << vecRLU[1] << ", " << vecRLU[2] << ") rlu";
			ostrTip << "\n(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ") 1/A";
			m_pPlot->SetObjectLabel(iCurObjIdx, ostrTip.str());

			++iCurObjIdx;
		}
	}

	// current q position
	if(bShowCurq)
	{
		m_pPlot->PlotSphere(m_vecq, 0.02/**0.1*g_dFontSize*/, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColCurq);
		m_pPlot->SetObjectAnimation(iCurObjIdx, 1);

		// current q label
		std::ostringstream ostrTip;
		ostrTip.precision(g_iPrecGfx);
		ostrTip << "q = (" << m_vecq_rlu[0] << ", " << m_vecq_rlu[1] << ", " << m_vecq_rlu[2] << ") rlu";
		ostrTip << "\nq = (" << m_vecq[0] << ", " << m_vecq[1] << ", " << m_vecq[2] << ") 1/A";
		m_pPlot->SetObjectLabel(iCurObjIdx, ostrTip.str());

		m_iqIdx = iCurObjIdx;
		++iCurObjIdx;
	}

	m_pPlot->SetMinMax(vecMin, vecMax);
	m_pPlot->SetEnabled(true);
}


// ----------------------------------------------------------------------------


/**
 * scattering triangle changed
 */
void BZ3DDlg::RecipParamsChanged(const RecipParams& recip)
{
	m_vecq = tl::make_vec({ -recip.q[0], -recip.q[1], -recip.q[2] });
	m_vecq_rlu = tl::make_vec({ -recip.q_rlu[0], -recip.q_rlu[1], -recip.q_rlu[2] });

	// if a 3d object is already assigned, update it
	if(m_pPlot && m_iqIdx)
	{
		m_pPlot->SetEnabled(false);

		m_pPlot->PlotSphere(m_vecq, 0.02/**0.1*g_dFontSize*/, *m_iqIdx);

		// current q label
		std::ostringstream ostrTip;
		ostrTip.precision(g_iPrecGfx);
		ostrTip << "q = (" << m_vecq_rlu[0] << ", " << m_vecq_rlu[1] << ", " << m_vecq_rlu[2] << ") rlu";
		ostrTip << "\nq = (" << m_vecq[0] << ", " << m_vecq[1] << ", " << m_vecq[2] << ") 1/A";
		m_pPlot->SetObjectLabel(*m_iqIdx, ostrTip.str());

		m_pPlot->SetEnabled(true);
	}
}


// ----------------------------------------------------------------------------


void BZ3DDlg::keyPressEvent(QKeyEvent* pEvt)
{
	if(m_pPlot)
		m_pPlot->keyPressEvent(pEvt);

	QDialog::keyPressEvent(pEvt);
}


void BZ3DDlg::closeEvent(QCloseEvent* pEvt)
{
	if(m_pSettings)
		m_pSettings->setValue("bz3d/geo", saveGeometry());
	QDialog::closeEvent(pEvt);
}


void BZ3DDlg::hideEvent(QHideEvent *pEvt)
{
	if(m_pPlot)
		m_pPlot->SetEnabled(false);
	QDialog::hideEvent(pEvt);
}


void BZ3DDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
	if(m_pPlot)
		m_pPlot->SetEnabled(true);
}


#include "moc_bz3d.cpp"
