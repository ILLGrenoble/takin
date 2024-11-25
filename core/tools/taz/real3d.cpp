/**
 * 3d unit cell drawing
 * @author Tobias Weber <tobias.weber@tum.de>
 * @modif_by Victor Mecoli <mecoli@ill.fr>
 * @date oct-2016
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

#include "real3d.h"
#include <QGridLayout>


using t_real = t_real_glob;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;


Real3DDlg::Real3DDlg(QWidget* pParent, QSettings *pSettings)
	: QDialog(pParent, Qt::Tool), m_pSettings(pSettings),
	m_pStatus(new QStatusBar(this)),
	m_pPerspective(new QPushButton("Perspective projection", pParent)),
	m_pTransparency(new QPushButton("Transparency", pParent)),
	m_pDrawFaces(new QPushButton("Draw Faces", pParent)),
	m_pDrawEdges(new QPushButton("Draw Edges", pParent)),
	m_pDrawSpheres(new QPushButton("Draw Spheres", pParent)),
	m_pPlot(new PlotGl(this, pSettings, 0.25))
{
	m_pPlot->SetEnabled(false);

	setWindowTitle("Real Space / Unit Cell");
	m_pStatus->setSizeGripEnabled(1);
	m_pStatus->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("real3d/geo"))
			restoreGeometry(m_pSettings->value("real3d/geo").toByteArray());
		else
			resize(800, 600);
	}

	m_pPlot->SetPrec(g_iPrecGfx);
	m_pPlot->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

	connect(m_pPerspective, &QPushButton::clicked, this, &Real3DDlg::onPerspectiveClicked);
	connect(m_pTransparency, &QPushButton::clicked, this, &Real3DDlg::onTransparencyClicked);
	connect(m_pDrawFaces, &QPushButton::clicked, this, &Real3DDlg::onDrawFacesClicked);
	connect(m_pDrawEdges, &QPushButton::clicked, this, &Real3DDlg::onDrawEdgesClicked);
	connect(m_pDrawSpheres, &QPushButton::clicked, this, &Real3DDlg::onDrawSpheresClicked);

	QGridLayout *gridLayout = new QGridLayout(this);
	QGridLayout *gridButton = new QGridLayout(this);
	gridLayout->setContentsMargins(4, 4, 4, 4);
	gridLayout->addWidget(m_pPlot.get(), 0, 0, 1, 1);

	gridLayout->addLayout(gridButton, 1, 0, 1, 1);
	gridButton->addWidget(m_pPerspective, 0, 0, 1, 3);
	gridButton->addWidget(m_pTransparency, 0, 3, 1, 3);
	gridButton->addWidget(m_pDrawFaces, 1, 0, 1, 2);
	gridButton->addWidget(m_pDrawEdges, 1, 2, 1, 2);
	gridButton->addWidget(m_pDrawSpheres, 1, 4, 1, 2);

	gridLayout->addWidget(m_pStatus, 2, 0, 1, 1);

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

	m_pPlot->SetLabels("a (A)", "b (A)", "c (A)");
	m_pPlot->SetDrawMinMax(0);
	m_pPlot->SetEnabled(true);
}


Real3DDlg::~Real3DDlg()
{
	if(m_pPlot)
		m_pPlot->SetEnabled(false);
}


void Real3DDlg::CalcPeaks(const tl::Brillouin3D<t_real_glob>& ws,
	const xtl::LatticeCommon<t_real_glob>& lattice)
{
	if(!ws.IsValid() || !m_pPlot)
		return;

	const bool bShowUCVerts = false;
	const bool bShowUCPolys = false;

	m_pPlot->SetEnabled(false);
	m_pPlot->clear();
	m_pPlot->SetObjectCount(lattice.vecAllAtoms.size());

	t_mat matAinv = tl::transpose(lattice.matB) / (tl::get_pi<t_real>()*t_real(2));


	// colours
	t_real dA = 0.7;
	std::vector<std::vector<t_real>> vecColor = {
		{0., 0., 1., dA},
		{1., 0., 0., dA},
		{0., 1., 0., dA},
		{1., 1., 0., dA},
		{0., 1., 1., dA},
		{1., 0., 1., dA},

		{0.5, 0., 0., dA},
		{0., 0.5, 0., dA},
		{0., 0., 0.5, dA},
		{0.5, 0.5, 0., dA},
		{0., 0.5, 0.5, dA},
		{0.5, 0., 0.5, dA},

		{1., 1., 1., dA},
		{0., 0., 0., dA},
		{0.5, 0.5, 0.5, dA},
	};

	static const std::vector<t_real> vecColVertices = { 0., 0., 0., 0.85 };
	static const std::vector<t_real> vecColPolys = { 0., 0., 1., 0.85 };
	static const std::vector<t_real> vecColEdges = { 0., 0., 0., 1. };
	static const std::vector<t_real> vecColScatPlane = { 1., 1., 0., 1. };


	std::size_t iCurObjIdx = 0;
	const t_real dLimMax = std::numeric_limits<t_real>::max();
	std::vector<t_real> vecMin = {dLimMax, dLimMax, dLimMax},
		vecMax = {-dLimMax, -dLimMax, -dLimMax};

	t_real d = t_real(0.5);
	t_vec vecPeaks[] = {
		lattice.lattice.GetPos(-d,0,0), lattice.lattice.GetPos(d,0,0),
		lattice.lattice.GetPos(0,-d,0), lattice.lattice.GetPos(0,d,0),
		lattice.lattice.GetPos(0,0,-d), lattice.lattice.GetPos(0,0,d),
	};


	// minimum and maximum coordinates
	for(const t_vec& vecPeak : vecPeaks)
	{
		for(unsigned int i = 0; i < 3; ++i)
		{
			vecMin[i] = std::min(vecPeak[i], vecMin[i]);
			vecMax[i] = std::max(vecPeak[i], vecMax[i]);
		}
	}


	// iterate atoms
	std::shared_ptr<const xtl::PeriodicSystem<t_real>> pElems =
		xtl::PeriodicSystem<t_real>::GetInstance();

	for(std::size_t iAtom=0; iAtom<lattice.vecAllAtoms.size(); ++iAtom)
	{
		const std::string& strElem = lattice.vecAllNames[iAtom];
		const t_vec& vecThisAtom = lattice.vecAllAtoms[iAtom];
		const t_vec& vecThisAtomFrac = lattice.vecAllAtomsFrac[iAtom];
		std::size_t iCurAtomType = lattice.vecAllAtomTypes[iAtom];

		// use covalent atom radii if available
		t_real dElemRad = 0.1;
		const auto* pElem = pElems->Find(strElem);
		if(pElem)
			dElemRad = pElem->GetRadiusCov()*0.1;

		m_pPlot->PlotSphere(vecThisAtom, dElemRad, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColor[iCurAtomType % vecColor.size()]);

		std::ostringstream ostrTip;
		ostrTip.precision(g_iPrecGfx);
		ostrTip << strElem;
		ostrTip << "\n("
			<< vecThisAtomFrac[0] << ", "
			<< vecThisAtomFrac[1] << ", "
			<< vecThisAtomFrac[2] << ") frac";
		ostrTip << "\n("
			<< vecThisAtom[0] << ", "
			<< vecThisAtom[1] << ", "
			<< vecThisAtom[2] << ") A";

		m_pPlot->SetObjectLabel(iCurObjIdx, ostrTip.str());
		++iCurObjIdx;


		// coordination polyhedra
		if(lattice.vecAllAtomPosAux.size() == lattice.vecAllAtoms.size())
		{
			for(auto vecPoly : lattice.vecAllAtomPosAux[iAtom].vecPolys)
			{
				// add centre to polyhedron
				for(auto& vecVert : vecPoly) vecVert += vecThisAtom;

				t_vec vecPolyNorm = tl::get_face_normal<t_vec, std::vector, t_real>(vecPoly, vecThisAtom);
				m_pPlot->PlotPoly(vecPoly, vecPolyNorm, iCurObjIdx);
				m_pPlot->SetObjectColor(iCurObjIdx, vecColor[iCurAtomType % vecColor.size()]);
				++iCurObjIdx;

				m_pPlot->PlotLines(vecPoly, 2.*0.1*g_dFontSize, iCurObjIdx);
				m_pPlot->SetObjectColor(iCurObjIdx, {0., 0., 0., 1.});
				++iCurObjIdx;
			}
		}
	}


	// render vertices
	for(const t_vec& vec : ws.GetVertices())
	{
		t_vec vecRLU = ublas::prod(matAinv, vec);
		tl::set_eps_0(vecRLU);

		if(bShowUCVerts)
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
	if(bShowUCPolys)
	{
		for(std::size_t iPoly=0; iPoly<ws.GetPolys().size(); ++iPoly)
		{
			m_pPlot->PlotPoly(ws.GetPolys()[iPoly], ws.GetPlanes()[iPoly].GetNorm(), iCurObjIdx);
			m_pPlot->SetObjectColor(iCurObjIdx, vecColPolys);

			++iCurObjIdx;
		}
	}

	// render edges
	for(const std::vector<t_vec>& vecPoly : ws.GetPolys())
	{
		m_pPlot->PlotLines(vecPoly, 2.*0.1*g_dFontSize, iCurObjIdx);
		m_pPlot->SetObjectColor(iCurObjIdx, vecColEdges);

		++iCurObjIdx;
	}


	m_pPlot->SetMinMax(vecMin, vecMax);
	m_pPlot->SetEnabled(true);
}



// ----------------------------------------------------------------------------


void Real3DDlg::keyPressEvent(QKeyEvent* pEvt)
{
	if(m_pPlot)
		m_pPlot->keyPressEvent(pEvt);

	QDialog::keyPressEvent(pEvt);
}


void Real3DDlg::closeEvent(QCloseEvent* pEvt)
{
	if(m_pSettings)
		m_pSettings->setValue("real3d/geo", saveGeometry());
	QDialog::closeEvent(pEvt);
}


void Real3DDlg::hideEvent(QHideEvent *pEvt)
{
	if(m_pPlot)
		m_pPlot->SetEnabled(false);
	QDialog::hideEvent(pEvt);
}


void Real3DDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
	if(m_pPlot)
		m_pPlot->SetEnabled(true);
}


void Real3DDlg::onPerspectiveClicked()
{
	m_pPlot->TogglePerspective();
}

void Real3DDlg::onTransparencyClicked()
{
	m_pPlot->ToggleZTest();
}

void Real3DDlg::onDrawFacesClicked()
{
	m_pPlot->ToggleDrawPolys();
}

void Real3DDlg::onDrawEdgesClicked()
{
	m_pPlot->ToggleDrawLines();
}

void Real3DDlg::onDrawSpheresClicked()
{
	m_pPlot->ToggleDrawSpheres();
}

#include "moc_real3d.cpp"
