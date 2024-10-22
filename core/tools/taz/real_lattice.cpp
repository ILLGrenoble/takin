/**
 * Real crystal lattice
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014 - 2016
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

#include "tlibs/helper/flags.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/phys/atoms.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/log/log.h"
#include "libs/globals.h"
#include "real_lattice.h"

#include <QToolTip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <tuple>


using t_real = t_real_glob;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;


#define REAL_LATTICE_NODE_TYPE_KEY	0

enum LatticeNodeType
{
	NODE_REAL_LATTICE,
	NODE_REAL_LATTICE_ATOM,

	NODE_REAL_LATTICE_OTHER
};



// --------------------------------------------------------------------------------

LatticePoint::LatticePoint()
{
	setFlag(QGraphicsItem::ItemIsMovable, false);
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
}

QRectF LatticePoint::boundingRect() const
{
	return QRectF(-3.5*g_dFontSize, -1.*g_dFontSize,
		7.*g_dFontSize, 5.*g_dFontSize);
}

void LatticePoint::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	pPainter->setFont(g_fontGfx);
	pPainter->setBrush(m_color);
	pPainter->drawEllipse(QRectF(-3.*0.1*g_dFontSize, -3.*0.1*g_dFontSize,
		6.*0.1*g_dFontSize, 6.*0.1*g_dFontSize));

	if(m_strLabel != "")
	{
		pPainter->setPen(m_color);
		QRectF rect = boundingRect();
		rect.setTop(rect.top()+1.65*g_dFontSize);
		pPainter->drawText(rect, Qt::AlignHCenter|Qt::AlignTop, m_strLabel);
	}
}

// --------------------------------------------------------------------------------

LatticeAtom::LatticeAtom()
{
	setFlag(QGraphicsItem::ItemIsMovable, false);
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
}

QRectF LatticeAtom::boundingRect() const
{
	return QRectF(-3.5*g_dFontSize, -1.*g_dFontSize,
		7.*g_dFontSize, 5.*g_dFontSize);
}

void LatticeAtom::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	pPainter->setFont(g_fontGfx);
	pPainter->setBrush(m_color);
	pPainter->setPen(Qt::darkCyan);
	pPainter->drawEllipse(QRectF(-3.*0.1*g_dFontSize, -3.*0.1*g_dFontSize,
		6.*0.1*g_dFontSize, 6.*0.1*g_dFontSize));

	if(m_strElem != "")
	{
		pPainter->setPen(m_color);
		QRectF rect = boundingRect();
		rect.setTop(rect.top()+1.65*g_dFontSize);
		//pPainter->drawRect(rect);
		pPainter->drawText(rect, Qt::AlignHCenter|Qt::AlignTop, m_strElem.c_str());
	}
}


// --------------------------------------------------------------------------------


RealLattice::RealLattice(LatticeScene& scene)
	: m_scene(scene)
{
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
	setAcceptedMouseButtons(Qt::NoButton);
	m_bReady = true;
}

RealLattice::~RealLattice()
{
	m_bReady = false;
	ClearPeaks();
}

QRectF RealLattice::boundingRect() const
{
	return QRectF(-100.*m_dZoom*g_dFontSize, -100.*m_dZoom*g_dFontSize,
		200.*m_dZoom*g_dFontSize, 200.*m_dZoom*g_dFontSize);
}

void RealLattice::SetZoom(t_real dZoom)
{
	m_dZoom = dZoom;
	m_scene.update();
}

void RealLattice::SetWSVisible(bool bVisible)
{
	m_bShowWS = bVisible;
	this->update();
}


/**
 * paint real lattice & unit cell
 */
void RealLattice::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	pPainter->setFont(g_fontGfx);

	if(!m_bShowWS || (!m_ws.IsValid() && !m_ws3.IsValid()))
		return;

	// Brillouin zone
	QPen penOrg = pPainter->pen();
	QPen penGray(Qt::darkGray);
	penGray.setWidthF(g_dFontSize*0.1);
	pPainter->setPen(penGray);

	t_vec vecCentral2d;
	std::vector<QPointF> vecWS3;

	// use 3d BZ code
	if(g_b3dBZ && m_ws3.IsValid())
	{
		// convert vertices to QPointFs
		vecWS3.reserve(m_vecWS3Verts.size());
		for(const auto& vecVert : m_vecWS3Verts)
			vecWS3.push_back(vec_to_qpoint(vecVert * m_dScaleFactor * m_dZoom));
	}
	// use 2d BZ code
	else if(m_ws.IsValid())
	{
		vecCentral2d = m_ws.GetCentralReflex() * m_dScaleFactor*m_dZoom;
	}

	for(const LatticePoint* pPeak : m_vecPeaks)
	{
		QPointF peakPos = pPeak->pos();
		peakPos *= m_dZoom;

		// use 3d BZ code
		if(g_b3dBZ && m_ws3.IsValid())
		{
			std::vector<QPointF> vecWS3_peak = vecWS3;
			for(auto& vecVert : vecWS3_peak)
				vecVert += peakPos;
			pPainter->drawPolygon(vecWS3_peak.data(), vecWS3_peak.size());
		}
		// use 2d BZ code
		else if(m_ws.IsValid())
		{
			const tl::Brillouin2D<t_real>::t_vertices<t_real>& verts = m_ws.GetVertices();
			for(const tl::Brillouin2D<t_real>::t_vecpair<t_real>& vertpair : verts)
			{
				const t_vec& vec1 = vertpair.first * m_dScaleFactor * m_dZoom;
				const t_vec& vec2 = vertpair.second * m_dScaleFactor * m_dZoom;

				QPointF pt1 = vec_to_qpoint(vec1 - vecCentral2d) + peakPos;
				QPointF pt2 = vec_to_qpoint(vec2 - vecCentral2d) + peakPos;

				QLineF lineWS(pt1, pt2);
				pPainter->drawLine(lineWS);
			}
		}
	}

	pPainter->setPen(penOrg);
}


/**
 * calculate real space representation
 */
void RealLattice::CalcPeaks(const xtl::LatticeCommon<t_real>& latticecommon)
{
	ClearPeaks();
	m_kdLattice.Unload();
	m_lattice = latticecommon.lattice;
	m_matPlane = latticecommon.matPlaneReal;
	m_matPlane_inv = latticecommon.matPlaneReal_inv;

	m_ws.SetEpsilon(g_dEps);
	m_ws.SetMaxNN(g_iMaxNN);
	m_ws3.SetEpsilon(g_dEps);
	m_ws3.SetMaxNN(g_iMaxNN);

	// central peak for WS cell calculation
	ublas::vector<int> veciCent = tl::make_vec({0.,0.,0.});

	static const std::string strAA = tl::get_spec_char_utf8("AA");

	// --------------------------------------------------------------------
	// atom positions in unit cell
	std::vector<QColor> colors = {QColor(127,0,0), QColor(0,127,0), QColor(0,0,127),
		QColor(127,127,0), QColor(0,127,127), QColor(127,0,127)};

	for(std::size_t iAtom = 0; iAtom < latticecommon.vecAllAtoms.size(); ++iAtom)
	{
		const std::string& strElem = latticecommon.vecAllNames[iAtom];
		const t_vec& vecThisAtom = latticecommon.vecAllAtoms[iAtom];
		const t_vec& vecThisAtomFrac = latticecommon.vecAllAtomsFrac[iAtom];
		std::size_t iCurAtomType = latticecommon.vecAllAtomTypes[iAtom];

		LatticeAtom *pAtom = new LatticeAtom();
		m_vecAtoms.push_back(pAtom);

		pAtom->m_strElem = strElem;
		pAtom->m_vecPos = std::move(vecThisAtom);
		pAtom->m_vecProj = latticecommon.planeReal.GetDroppedPerp(pAtom->m_vecPos/*, &pAtom->m_dProjDist*/);
		pAtom->m_dProjDist = latticecommon.planeReal.GetDist(pAtom->m_vecPos);

		t_vec vecCoord = ublas::prod(m_matPlane_inv, pAtom->m_vecProj);
		t_real dX = vecCoord[0], dY = -vecCoord[1];

		pAtom->setPos(dX * m_dScaleFactor, dY * m_dScaleFactor);
		pAtom->setData(REAL_LATTICE_NODE_TYPE_KEY, NODE_REAL_LATTICE_ATOM);

		std::ostringstream ostrTip;
		ostrTip.precision(g_iPrecGfx);
		ostrTip << pAtom->m_strElem;
		ostrTip << "\n("
			<< vecThisAtomFrac[0] << ", "
			<< vecThisAtomFrac[1] << ", "
			<< vecThisAtomFrac[2] << ") frac";
		ostrTip << "\n("
			<< vecThisAtom[0] << ", "
			<< vecThisAtom[1] << ", "
			<< vecThisAtom[2] << ") " << strAA;
		ostrTip << "\nDistance to Plane: " << pAtom->m_dProjDist << " " << strAA;
		pAtom->setToolTip(QString::fromUtf8(ostrTip.str().c_str(), ostrTip.str().length()));
		pAtom->SetColor(colors[iCurAtomType % colors.size()]);

		m_scene.addItem(pAtom);
	}
	// --------------------------------------------------------------------



	std::list<std::vector<t_real>> lstPeaksForKd;

	for(int ih = -m_iMaxPeaks; ih <= m_iMaxPeaks; ++ih)
	for(int ik = -m_iMaxPeaks; ik <= m_iMaxPeaks; ++ik)
	for(int il = -m_iMaxPeaks; il <= m_iMaxPeaks; ++il)
	{
		const t_real h = t_real(ih), k = t_real(ik), l = t_real(il);
		const t_vec vecPeakHKL = tl::make_vec<t_vec>({h,k,l});
		t_vec vecPeak = m_lattice.GetPos(h,k,l);

		// 3d unit cell
		if(g_b3dBZ)
		{
			if(ih==veciCent[0] && ik==veciCent[1] && il==veciCent[2])
				m_ws3.SetCentralReflex(vecPeak, &vecPeakHKL);
			else if(std::abs(ih-veciCent[0]) <= 2 && std::abs(ik-veciCent[1]) <= 2 && std::abs(il-veciCent[2]) <= 2)
				m_ws3.AddReflex(vecPeak, &vecPeakHKL);
		}

		// add peak in A and in fractional units
		lstPeaksForKd.push_back(std::vector<t_real>{vecPeak[0],vecPeak[1],vecPeak[2], h,k,l});

		t_real dDist = 0.;
		t_vec vecDropped = latticecommon.planeReal.GetDroppedPerp(vecPeak, &dDist);

		if(tl::float_equal<t_real>(dDist, 0., m_dPlaneDistTolerance))
		{
			t_vec vecCoord = ublas::prod(m_matPlane_inv, vecDropped);
			t_real dX = vecCoord[0], dY = -vecCoord[1];

			LatticePoint *pPeak = new LatticePoint();
			if(ih==0 && ik==0 && il==0)
				pPeak->SetColor(Qt::darkGreen);
			pPeak->setPos(dX * m_dScaleFactor, dY * m_dScaleFactor);
			pPeak->setData(REAL_LATTICE_NODE_TYPE_KEY, NODE_REAL_LATTICE);

			std::ostringstream ostrTip;
			ostrTip.precision(g_iPrecGfx);

			ostrTip << "(" << ih << " " << ik << " " << il << ")";
			pPeak->SetLabel(ostrTip.str().c_str());

			tl::set_eps_0(vecPeak, g_dEps);
			ostrTip << " frac\n";
			ostrTip << "("
				<< vecPeak[0] << ", "
				<< vecPeak[1] << ", "
				<< vecPeak[2] << ") " << strAA;
			pPeak->setToolTip(QString::fromUtf8(ostrTip.str().c_str(), ostrTip.str().length()));

			m_vecPeaks.push_back(pPeak);
			m_scene.addItem(pPeak);


			// 2d unit cell
			if(!g_b3dBZ)
			{
				if(ih==veciCent[0] && ik==veciCent[1] && il==veciCent[2])
				{
					t_vec vecCentral = tl::make_vec({dX, dY});
					m_ws.SetCentralReflex(vecCentral, &vecPeakHKL);
				}
				// TODO: check if 2 next neighbours is sufficient for all space groups
				else if(std::abs(ih-veciCent[0])<=2 && std::abs(ik-veciCent[1])<=2
					&& std::abs(il-veciCent[2])<=2)
				{
					t_vec vecN = tl::make_vec({dX, dY});
					m_ws.AddReflex(vecN, &vecPeakHKL);
				}
			}
		}
	}

	if(g_b3dBZ)
	{
		m_ws3.CalcBZ(get_max_threads());

		// ----------------------------------------------------------------
		// calculate intersection with real plane
		tl::Plane<t_real> planeBZ3 = tl::Plane<t_real>(m_ws3.GetCentralReflex(),
			latticecommon.planeReal.GetNorm());

		std::tie(std::ignore, m_vecWS3VertsUnproj) = m_ws3.GetIntersection(planeBZ3);

		for(const t_vec& _vecWS3Vert : m_vecWS3VertsUnproj)
		{
			t_vec vecWS3Vert = ublas::prod(latticecommon.matPlaneReal_inv, _vecWS3Vert - m_ws3.GetCentralReflex());
			vecWS3Vert.resize(2, true);
			vecWS3Vert[1] = -vecWS3Vert[1];

			m_vecWS3Verts.push_back(vecWS3Vert);
		}
		// ----------------------------------------------------------------

	}
	else
	{
		m_ws.CalcBZ();
	}

	m_kdLattice.Load(lstPeaksForKd, 3);

	this->update();
}

t_vec RealLattice::GetHKLFromPlanePos(t_real x, t_real y) const
{
	if(!HasPeaks())
		return t_vec();

	t_vec vec = x*tl::get_column(m_matPlane, 0)
		+ y*tl::get_column(m_matPlane, 1);
	return m_lattice.GetHKL(vec);
}

void RealLattice::ClearPeaks()
{
	m_ws.Clear();
	m_ws3.Clear();

	m_vecWS3VertsUnproj.clear();
	m_vecWS3Verts.clear();

	for(LatticePoint*& pPeak : m_vecPeaks)
	{
		if(pPeak)
		{
			m_scene.removeItem(pPeak);
			delete pPeak;
			pPeak = nullptr;
		}
	}
	m_vecPeaks.clear();

	for(LatticeAtom*& pAtom : m_vecAtoms)
	{
		if(pAtom)
		{
			m_scene.removeItem(pAtom);
			delete pAtom;
			pAtom = nullptr;
		}
	}
	m_vecAtoms.clear();
}

// --------------------------------------------------------------------------------


LatticeScene::LatticeScene(QObject *pParent)
	: QGraphicsScene(pParent), m_pLatt(new RealLattice(*this))
{
	this->addItem(m_pLatt);
}

LatticeScene::~LatticeScene()
{
	delete m_pLatt;
}

void LatticeScene::scaleChanged(t_real dTotalScale)
{
	if(!m_pLatt)
		return;
	m_pLatt->SetZoom(dTotalScale);
}

void LatticeScene::mousePressEvent(QGraphicsSceneMouseEvent *pEvt)
{
	m_bMousePressed = true;
	QGraphicsScene::mousePressEvent(pEvt);
}


void LatticeScene::drawBackground(QPainter* pPainter, const QRectF& rect)
{
	QGraphicsScene::drawBackground(pPainter, rect);
	// TODO: draw accurate WS cell
}

void LatticeScene::mouseMoveEvent(QGraphicsSceneMouseEvent *pEvt)
{
	bool bHandled = false;
	bool bAllowed = true;


	// tooltip
	if(m_pLatt)
	{
		const t_real dX = pEvt->scenePos().x()/m_pLatt->GetScaleFactor();
		const t_real dY = -pEvt->scenePos().y()/m_pLatt->GetScaleFactor();

		t_vec vecHKL = m_pLatt->GetHKLFromPlanePos(dX, dY);
		tl::set_eps_0(vecHKL, g_dEps);

		if(vecHKL.size() == 3)
		{
			const std::vector<t_real>* pvecNearest = nullptr;

			const tl::Kd<t_real>& kd = m_pLatt->GetKdLattice();
			const tl::Lattice<t_real>& lattice = m_pLatt->GetRealLattice();
			t_vec vecHKLA = lattice.GetPos(vecHKL[0], vecHKL[1], vecHKL[2]);

			if(kd.GetRootNode())
			{
				std::vector<t_real> stdvecHKL{vecHKLA[0], vecHKLA[1], vecHKLA[2]};
				pvecNearest = &kd.GetNearestNode(stdvecHKL);
				if(pvecNearest->size() < 6)
				{
					pvecNearest = nullptr;
					tl::log_warn("Invalid WS node.");
				}
			}

			emit coordsChanged(vecHKL[0], vecHKL[1], vecHKL[2],
				pvecNearest != nullptr,
				pvecNearest?(*pvecNearest)[3]:0.,
				pvecNearest?(*pvecNearest)[4]:0.,
				pvecNearest?(*pvecNearest)[5]:0.);
		}
	}


	// node dragging
	if(m_bMousePressed)
	{
		QGraphicsItem* pCurItem = mouseGrabberItem();
		if(pCurItem)
		{
			const int iNodeType = pCurItem->data(REAL_LATTICE_NODE_TYPE_KEY).toInt();

			// nothing there yet...
		}
	}

	if(!bHandled && bAllowed)
		QGraphicsScene::mouseMoveEvent(pEvt);
}

void LatticeScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt)
{
	m_bMousePressed = false;

	QGraphicsScene::mouseReleaseEvent(pEvt);
}

void LatticeScene::keyPressEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Control)
		m_bSnap = true;

	QGraphicsScene::keyPressEvent(pEvt);
}

void LatticeScene::keyReleaseEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Control)
		m_bSnap = false;

	QGraphicsScene::keyReleaseEvent(pEvt);
}


// --------------------------------------------------------------------------------


LatticeView::LatticeView(QWidget* pParent)
	: QGraphicsView(pParent)
{
	setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
		QPainter::SmoothPixmapTransform | QPainter::HighQualityAntialiasing);
	setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);

	setDragMode(QGraphicsView::ScrollHandDrag);
	setMouseTracking(true);
	setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
}

LatticeView::~LatticeView()
{}


void LatticeView::DoZoom(t_real_glob dDelta)
{
	t_real dScale = std::pow(2., dDelta);
	this->scale(dScale, dScale);

	m_dTotalScale *= dScale;
	emit scaleChanged(m_dTotalScale);
}


void LatticeView::keyPressEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Plus)
		DoZoom(0.02);
	else if(pEvt->key() == Qt::Key_Minus)
		DoZoom(-0.02);

	QGraphicsView::keyPressEvent(pEvt);
}

void LatticeView::keyReleaseEvent(QKeyEvent *pEvt)
{
	QGraphicsView::keyReleaseEvent(pEvt);
}

void LatticeView::wheelEvent(QWheelEvent *pEvt)
{
	const t_real dDelta = pEvt->angleDelta().y()/8. / 150.;
	DoZoom(dDelta);
}


#include "moc_real_lattice.cpp"
