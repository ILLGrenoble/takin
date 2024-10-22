/**
 * Reciprocal Lattice
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014 - 2018
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

#include "tlibs/helper/flags.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/log/log.h"
#include "scattering_triangle.h"

#include <QToolTip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <tuple>


// symbol drawing sizes
#define DEF_PEAK_SIZE 3.
#define MIN_PEAK_SIZE 0.2
#define MAX_PEAK_SIZE 5.5

#define CALC_BZ_AROUND_ZERO


using t_real = t_real_glob;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;
static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();


static inline bool flip_text(t_real _dAngle)
{
	t_real dAngle = std::fmod(_dAngle, 360.);
	return std::abs(dAngle) > 90. && std::abs(dAngle) < 270.;
}


ScatteringTriangleNode::ScatteringTriangleNode(ScatteringTriangle* pSupItem)
	: m_pParentItem(pSupItem)
{
	setFlag(QGraphicsItem::ItemSendsGeometryChanges);
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
	setCursor(Qt::CrossCursor);

	setData(TRIANGLE_NODE_TYPE_KEY, NODE_OTHER);
}


QRectF ScatteringTriangleNode::boundingRect() const
{
	return QRectF(-0.5*g_dFontSize, -0.5*g_dFontSize, g_dFontSize, g_dFontSize);
}


void ScatteringTriangleNode::paint(QPainter *pPainter, const QStyleOptionGraphicsItem* popt, QWidget* pwid)
{
	pPainter->drawEllipse(QRectF(-2., -2., 4., 4.));
}


void ScatteringTriangleNode::mousePressEvent(QGraphicsSceneMouseEvent *pEvt)
{
	//setCursor(Qt::ClosedHandCursor);
	QGraphicsItem::mousePressEvent(pEvt);
}


void ScatteringTriangleNode::mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt)
{
	//setCursor(Qt::OpenHandCursor);
	QGraphicsItem::mouseReleaseEvent(pEvt);
}


QVariant ScatteringTriangleNode::itemChange(GraphicsItemChange change, const QVariant &val)
{
	QVariant var = QGraphicsItem::itemChange(change, val);

	if(change == QGraphicsItem::ItemPositionHasChanged)
		m_pParentItem->nodeMoved(this);

	return var;
}


// --------------------------------------------------------------------------------


RecipPeak::RecipPeak()
{
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
	setFlag(QGraphicsItem::ItemIsMovable, false);
}


QRectF RecipPeak::boundingRect() const
{
	return QRectF(-4*g_dFontSize, -g_dFontSize,	// l, t
		8.*g_dFontSize, 5.*g_dFontSize);	// w, h
}


void RecipPeak::paint(QPainter *pPainter, const QStyleOptionGraphicsItem* pOpt, QWidget* pWid)
{
	QFont font = g_fontGfx;
	font.setStrikeOut(!m_bPeakAllowed);

	pPainter->setFont(font);
	pPainter->setBrush(m_color);
	pPainter->drawEllipse(QRectF(-m_dRadius*0.1*g_dFontSize, -m_dRadius*0.1*g_dFontSize,
		m_dRadius*2.*0.1*g_dFontSize, m_dRadius*2.*0.1*g_dFontSize));

	if(m_strLabel != "")
	{
		pPainter->setPen(m_color);
		QRectF rect = boundingRect();
		rect.setTop(rect.top()+1.65*g_dFontSize);
		rect.setBottom(rect.bottom()+1.65*g_dFontSize);
		pPainter->drawText(rect, Qt::AlignHCenter|Qt::AlignTop, m_strLabel);
	}
}


// --------------------------------------------------------------------------------


ScatteringTriangle::ScatteringTriangle(ScatteringTriangleScene& scene)
	: m_scene(scene),
	m_pNodeKiQ(new ScatteringTriangleNode(this)),
	m_pNodeKiKf(new ScatteringTriangleNode(this)),
	m_pNodeKfQ(new ScatteringTriangleNode(this)),
	m_pNodeGq(new ScatteringTriangleNode(this))
{
	setFlag(QGraphicsItem::ItemIgnoresTransformations);

	m_pNodeKiQ->setData(TRIANGLE_NODE_TYPE_KEY, NODE_KIQ);
	m_pNodeKiKf->setData(TRIANGLE_NODE_TYPE_KEY, NODE_KIKF);
	m_pNodeKfQ->setData(TRIANGLE_NODE_TYPE_KEY, NODE_Q);
	m_pNodeGq->setData(TRIANGLE_NODE_TYPE_KEY, NODE_q);

	AllowMouseMove(true);

	m_pNodeKiQ->setPos(0., 0.);
	m_pNodeKiKf->setPos(80., -150.);
	m_pNodeKfQ->setPos(160., 0.);
	m_pNodeGq->setPos(160., 0.);

	m_scene.addItem(m_pNodeKiQ.get());
	m_scene.addItem(m_pNodeKiKf.get());
	m_scene.addItem(m_pNodeKfQ.get());
	m_scene.addItem(m_pNodeGq.get());

	setAcceptedMouseButtons(Qt::NoButton);
	m_bReady = m_bUpdate = true;
}


ScatteringTriangle::~ScatteringTriangle()
{
	m_bUpdate = m_bReady = false;
	ClearPeaks();
}


void ScatteringTriangle::AllowMouseMove(bool bAllow)
{
	m_pNodeKiKf->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pNodeKfQ->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pNodeGq->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
}


void ScatteringTriangle::nodeMoved(const ScatteringTriangleNode* pNode)
{
	if(!m_bReady) return;

	if(m_scene.getSnapq() && pNode==GetNodeKfQ())
		SnapToNearestPeak(GetNodeGq(), GetNodeKfQ());

	if(m_bUpdate)
	{
		update();
		m_scene.emitUpdate();
		m_scene.emitAllParams();
	}
}


QRectF ScatteringTriangle::boundingRect() const
{
	return QRectF(-100.*m_dZoom*g_dFontSize, -100.*m_dZoom*g_dFontSize,
		200.*m_dZoom*g_dFontSize, 200.*m_dZoom*g_dFontSize);
}


void ScatteringTriangle::SetZoom(t_real dZoom)
{
	m_dZoom = dZoom;
	m_scene.update();
}


void ScatteringTriangle::SetqVisible(bool bVisible)
{
	m_bqVisible = bVisible;
	this->m_pNodeGq->setVisible(bVisible);
	this->update();
}


void ScatteringTriangle::SetCoordAxesVisible(bool bVisible)
{
	m_bCoordAxesVisible = bVisible;
	this->update();
}


void ScatteringTriangle::SetBZVisible(bool bVisible)
{
	m_bShowBZ = bVisible;
	this->update();
}


void ScatteringTriangle::SetAllPeaksVisible(bool bVisible)
{
	m_bShowAllPeaks = bVisible;
	this->update();
}


void ScatteringTriangle::SetEwaldSphereVisible(EwaldSphere iEw)
{
	m_bShowEwaldSphere = (iEw != EWALD_NONE);
	m_bEwaldAroundKi = (iEw == EWALD_KI);
	this->update();
}


QPointF ScatteringTriangle::GetGfxMid() const
{
	QPointF ptKiQ = mapFromItem(m_pNodeKiQ.get(), 0, 0);
	QPointF ptKfQ = mapFromItem(m_pNodeKfQ.get(), 0, 0);
	QPointF ptKiKf = mapFromItem(m_pNodeKiKf.get(), 0, 0);

	return (ptKiQ + ptKfQ + ptKiKf) / 3.;
}


void ScatteringTriangle::paint(QPainter *pPainter, const QStyleOptionGraphicsItem* pOpt, QWidget* pWid)
{
	QPen penOrg = pPainter->pen();
	pPainter->setFont(g_fontGfx);

	QPen penRed(Qt::red);
	penRed.setWidthF(g_dFontSize*0.1);
	QPen penBlack(qApp->palette().color(QPalette::WindowText));
	penBlack.setWidthF(g_dFontSize*0.1);
	QPen penGreen(Qt::darkGreen);
	penGreen.setWidthF(g_dFontSize*0.1);
	QPen penLight(QColor(0xaa, 0xaa, 0xaa));
	penLight.setWidthF(g_dFontSize*0.1);
	QPen penBlue(qApp->palette().color(QPalette::Link));
	penBlue.setWidthF(g_dFontSize*0.1);
	QPen penGray(Qt::darkGray);
	penGray.setWidthF(g_dFontSize*0.1);


	// Brillouin zone
	if(m_bShowBZ && (m_bz.IsValid() || m_bz3.IsValid()))
	{
		pPainter->setPen(penGray);

		t_vec vecCentral2d;
		std::vector<QPointF> vecBZ3;

		// use 3d BZ code
		if(g_b3dBZ && m_bz3.IsValid())
		{
			// convert vertices to QPointFs
			vecBZ3.reserve(m_vecBZ3Verts.size());
			for(const auto& vecVert : m_vecBZ3Verts)
				vecBZ3.push_back(vec_to_qpoint(vecVert * m_dScaleFactor * m_dZoom));
		}
		// use 2d BZ code
		else if(m_bz.IsValid())
		{
			vecCentral2d = m_bz.GetCentralReflex() * m_dScaleFactor*m_dZoom;
		}

		for(const RecipPeak* pPeak : m_vecPeaks)
		{
			QPointF peakPos = pPeak->pos();
			peakPos *= m_dZoom;

			// use 3d BZ code
			if(g_b3dBZ && m_bz3.IsValid())
			{
				std::vector<QPointF> vecBZ3_peak = vecBZ3;
				for(auto& vecVert : vecBZ3_peak)
					vecVert += peakPos;
				pPainter->drawPolygon(vecBZ3_peak.data(), vecBZ3_peak.size());
			}
			// use 2d BZ code
			else if(m_bz.IsValid())
			{
				const tl::Brillouin2D<t_real>::t_vertices<t_real>& verts = m_bz.GetVertices();
				for(const tl::Brillouin2D<t_real>::t_vecpair<t_real>& vertpair : verts)
				{
					const t_vec& vec1 = vertpair.first * m_dScaleFactor * m_dZoom;
					const t_vec& vec2 = vertpair.second * m_dScaleFactor * m_dZoom;

					QPointF pt1 = vec_to_qpoint(vec1 - vecCentral2d) + peakPos;
					QPointF pt2 = vec_to_qpoint(vec2 - vecCentral2d) + peakPos;

					QLineF lineBZ(pt1, pt2);
					pPainter->drawLine(lineBZ);
				}
			}
		}

		pPainter->setPen(penOrg);
	}


	// orientation vectors
	QPointF ptOrient[2];

	for(int iOrient=0; iOrient<2; ++iOrient)
	{
		// vecOrient_rlu is normalised
		t_vec vecOrient_rlu = tl::get_column(m_matPlaneRlu, iOrient);
		t_vec vecOrient_lab = m_recip.GetPos(vecOrient_rlu[0], vecOrient_rlu[1], vecOrient_rlu[2]);

		t_real dDistOrient = 0.;
		t_vec vecDroppedOrient = m_plane.GetDroppedPerp(vecOrient_lab, &dDistOrient);
		//bool bOrient0InPlane = tl::float_equal<t_real>(dDistOrient[iOrient], 0., m_dPlaneDistTolerance);

		t_vec vecCoordOrient = ublas::prod(m_matPlane_inv, vecDroppedOrient);

		ptOrient[iOrient] = QPointF(vecCoordOrient[0], -vecCoordOrient[1]);
		ptOrient[iOrient] *= m_dScaleFactor*m_dZoom * 0.75;
	}



	QPointF ptKiQ = mapFromItem(m_pNodeKiQ.get(), 0, 0) * m_dZoom;
	QPointF ptKfQ = mapFromItem(m_pNodeKfQ.get(), 0, 0) * m_dZoom;
	QPointF ptKiKf = mapFromItem(m_pNodeKiKf.get(), 0, 0) * m_dZoom;
	QPointF ptGq = mapFromItem(m_pNodeGq.get(), 0, 0) * m_dZoom;


	// powder lines
	{
		for(std::size_t iLine=0; iLine<m_vecPowderLines.size(); ++iLine)
		{
			const typename tl::Powder<int,t_real>::t_peak& powderpeak = m_vecPowderLines[iLine];
			t_real dLineWidth = m_vecPowderLineWidths[iLine] * g_dFontSize*0.1;

			const int ih = std::get<0>(powderpeak);
			const int ik = std::get<1>(powderpeak);
			const int il = std::get<2>(powderpeak);
			t_real dF = std::get<4>(powderpeak);
			bool bHasF = true;

			if(ih==0 && ik==0 && il==0) continue;
			if(dF < t_real(0))
			{
				bHasF = false;
				dF = t_real(1);
			}

			std::ostringstream ostrPowderLine;
			ostrPowderLine.precision(g_iPrecGfx);
			ostrPowderLine << "(" << ih << " "<< ik << " " << il << ")";
			if(bHasF)
				ostrPowderLine << ", F=" << dF;

			t_vec vec = m_recip.GetPos(t_real(ih), t_real(ik), t_real(il));
			t_real drad = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
			drad *= m_dScaleFactor*m_dZoom;

			QPen penLine(penRed.color());
			penLine.setWidthF(dLineWidth);
			pPainter->setPen(penLine);

			pPainter->drawEllipse(ptKiQ, drad, drad);
			pPainter->drawText(ptKiQ + QPointF(0., drad), ostrPowderLine.str().c_str());
		}

		pPainter->setPen(penOrg);
	}


	QLineF lineQ(ptKiQ, ptKfQ);
	QLineF lineKi(ptKiQ, ptKiKf);
	QLineF lineKf(ptKiKf, ptKfQ);
	QLineF lineG(ptKiQ, ptGq);
	QLineF lineq(ptKfQ, ptGq);
	QLineF lineOrient0(ptOrient[0], ptKiQ);
	QLineF lineOrient1(ptOrient[1], ptKiQ);


	if(m_bCoordAxesVisible)
	{
		pPainter->setPen(penLight);
		pPainter->drawLine(lineOrient0);
		pPainter->drawLine(lineOrient1);
	}

	pPainter->setPen(penRed);
	pPainter->drawLine(lineQ);

	pPainter->setPen(penBlack);
	pPainter->drawLine(lineKi);
	pPainter->drawLine(lineKf);

	if(m_bqVisible)
	{
		pPainter->setPen(penGreen);

		pPainter->drawLine(lineG);
		pPainter->drawLine(lineq);

		pPainter->setPen(penOrg);
	}

	const t_real dG = lineG.length()/m_dScaleFactor/m_dZoom;

	const std::wstring strAA = tl::get_spec_char_utf16("AA") + tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");
	const std::wstring& strDelta = tl::get_spec_char_utf16("Delta");

	std::wostringstream ostrQ, ostrKi, ostrKf, ostrE, ostrq, ostrG;
	ostrQ.precision(g_iPrecGfx); ostrE.precision(g_iPrecGfx);
	ostrKi.precision(g_iPrecGfx); ostrKf.precision(g_iPrecGfx);
	ostrG.precision(g_iPrecGfx); ostrq.precision(g_iPrecGfx);

	t_real dQ = GetQ();	tl::set_eps_0(dQ, g_dEpsGfx);
	t_real dKi = GetKi();	tl::set_eps_0(dKi, g_dEpsGfx);
	t_real dKf = GetKf();	tl::set_eps_0(dKf, g_dEpsGfx);
	t_real dE = GetE();	tl::set_eps_0(dE, g_dEpsGfx);

	ostrQ << L"Q = " << dQ << " " << strAA;
	ostrKi << L"ki = " << dKi << " " << strAA;
	ostrKf << L"kf = " << dKf << " " << strAA;
	ostrE << strDelta << "E = " << dE << " meV";
	if(m_bqVisible)
	{
		ostrq << L"q = " << Getq() << " " << strAA;
		ostrG << L"G = " << dG << " " << strAA;
	}

	pPainter->save();
	t_real dAngleQ = -lineQ.angle();
	pPainter->rotate(dAngleQ);
	pPainter->setPen(penRed);
	pPainter->translate(QPointF(lineQ.length()/5., 16.));
	if(flip_text(dAngleQ))
	{
		pPainter->translate(QPointF(lineQ.length()/2., -10.));
		pPainter->rotate(180.);
	}
	pPainter->drawText(QPointF(0.,0.5*g_dFontSize),
		QString::fromWCharArray(ostrQ.str().c_str()));
	pPainter->restore();

	pPainter->save();
	t_real dAngleKi = -lineKi.angle();
	pPainter->rotate(dAngleKi);
	pPainter->setPen(penBlack);
	pPainter->translate(QPointF(lineKi.length()/5., -4.));
	if(flip_text(dAngleKi))
	{
		pPainter->translate(QPointF(lineKi.length()/2., -12.));
		pPainter->rotate(180.);
	}
	pPainter->drawText(QPointF(0.,0.),
		QString::fromWCharArray(ostrKi.str().c_str()));
	pPainter->restore();

	pPainter->save();
	t_real dAngleKf = -lineKf.angle();
	pPainter->translate(ptKiKf);
	pPainter->rotate(dAngleKf);
	pPainter->translate(QPointF(lineKf.length()/5., -4.));
	if(flip_text(dAngleKf))
	{
		pPainter->translate(QPointF(lineKf.length()/2., 8.));
		pPainter->rotate(180.);
	}
	pPainter->drawText(QPointF(0.,0.),
		QString::fromWCharArray(ostrKf.str().c_str()));
	pPainter->drawText(QPointF(0.,16.*0.1*g_dFontSize),
		QString::fromWCharArray(ostrE.str().c_str()));
	pPainter->restore();

	if(m_bqVisible)
	{
		pPainter->setPen(penGreen);

		pPainter->save();
		t_real dAngleq = -lineq.angle();
		pPainter->translate(ptKfQ);
		pPainter->rotate(dAngleq);
		pPainter->translate(QPointF(lineq.length()/5.,-4.));
		if(flip_text(dAngleq))
		{
			pPainter->translate(QPointF(lineq.length()/2., 8.));
			pPainter->rotate(180.);
		}
		pPainter->drawText(QPointF(0.,0.),
			QString::fromWCharArray(ostrq.str().c_str()));
		pPainter->restore();

		pPainter->save();
		t_real dAngleG = -lineG.angle();
		pPainter->rotate(dAngleG);
		pPainter->translate(QPointF(lineG.length()/5.,-4.));
		if(flip_text(dAngleG))
		{
			pPainter->translate(QPointF(lineG.length()/2., 8.));
			pPainter->rotate(180.);
		}
		pPainter->drawText(QPointF(0.,0.),
			QString::fromWCharArray(ostrG.str().c_str()));
		pPainter->rotate(lineG.angle());
		pPainter->restore();

		pPainter->setPen(penOrg);
	}


	QLineF lineKi2(ptKiKf, ptKiQ);
	QLineF lineKf2(ptKfQ, ptKiKf);
	QLineF lineQ2(ptKfQ, ptKiQ);

	std::vector<const QLineF*> vecLines1 = {&lineKi2, &lineKi, &lineKf};
	std::vector<const QLineF*> vecLines2 = {&lineQ, &lineKf, &lineQ2};
	std::vector<const QPointF*> vecPoints = {&ptKiQ, &ptKiKf, &ptKfQ};

	std::vector<const QLineF*> vecLinesArrow = {&lineKi2, &lineKf, &lineQ2};
	std::vector<const QPointF*> vecPointsArrow = {&ptKiKf, &ptKiKf, &ptKfQ};

	std::vector<bool> vecDrawAngles = {1,1,1};
	std::vector<QColor> vecColor {penBlack.color(), penBlack.color(), penRed.color()};
	std::vector<QColor> vecColArc {penBlue.color(), penBlue.color(), penBlue.color()};
	std::vector<t_real> vecArcScale {1, 1, 1};

	QLineF lineG2(ptGq, ptKiQ);
	QLineF lineq2(ptGq, ptKfQ);

	if(m_bCoordAxesVisible)
	{
		// angle between Q and orient0
		vecLines1.push_back(&lineOrient0);
		vecLines2.push_back(&lineQ);
		vecPoints.push_back(&ptKiQ);
		vecColArc.push_back(penLight.color());
		vecArcScale.push_back(0.75);

		vecLinesArrow.push_back(&lineOrient0);
		vecLinesArrow.push_back(&lineOrient1);

		vecPointsArrow.push_back(&ptOrient[0]);
		vecPointsArrow.push_back(&ptOrient[1]);

		vecDrawAngles.push_back(1);
		vecDrawAngles.push_back(0);

		vecColor.push_back(penLight.color());
		vecColor.push_back(penLight.color());
	}

	if(m_bqVisible)
	{
		vecLinesArrow.push_back(&lineq);
		vecLinesArrow.push_back(&lineG2);

		vecPointsArrow.push_back(&ptKfQ);
		vecPointsArrow.push_back(&ptGq);

		vecDrawAngles.push_back(0);
		vecDrawAngles.push_back(0);

		vecColor.push_back(penGreen.color());
		vecColor.push_back(penGreen.color());
	}


	const std::wstring& strDEG = tl::get_spec_char_utf16("deg");

	for(std::size_t i=0; i<vecDrawAngles.size(); ++i)
	{
		// arrow heads
		t_real dAng = tl::d2r<t_real>(vecLinesArrow[i]->angle() + 90.);
		t_real dC = std::cos(dAng);
		t_real dS = std::sin(dAng);

		t_real dTriagX = 0.5*g_dFontSize, dTriagY = 1.*g_dFontSize;
		QPointF ptTriag1 = *vecPointsArrow[i] + QPointF(dTriagX*dC + dTriagY*dS, -dTriagX*dS + dTriagY*dC);
		QPointF ptTriag2 = *vecPointsArrow[i] + QPointF(-dTriagX*dC + dTriagY*dS, dTriagX*dS + dTriagY*dC);

		QPainterPath triag;
		triag.moveTo(*vecPointsArrow[i]);
		triag.lineTo(ptTriag1);
		triag.lineTo(ptTriag2);

		QPen penCol(vecColor[i]);
		penCol.setWidthF(g_dFontSize*0.1);

		pPainter->setPen(penCol);
		pPainter->fillPath(triag, vecColor[i]);

		if(vecDrawAngles[i])
		{
			// angle arcs
			t_real dArcSize = vecArcScale[i] *
				(vecLines1[i]->length() + vecLines2[i]->length()) / 2. / 3.;
			t_real dBeginArcAngle = vecLines1[i]->angle() + 180.;
			t_real dArcAngle = vecLines1[i]->angleTo(*vecLines2[i]) - 180.;

			QPen penArc(vecColArc[i]);
			penArc.setWidthF(g_dFontSize*0.1);

			pPainter->setPen(penArc);
			pPainter->drawArc(QRectF(vecPoints[i]->x()-dArcSize/2., vecPoints[i]->y()-dArcSize/2.,
				dArcSize, dArcSize), dBeginArcAngle*16., dArcAngle*16.);

			std::wostringstream ostrAngle;
			ostrAngle.precision(g_iPrecGfx);
			ostrAngle << std::fabs(dArcAngle) << strDEG;


			t_real dTotalAngle = -dBeginArcAngle-dArcAngle*0.5 + 180.;
			t_real dTransScale = 50. * m_dZoom * vecArcScale[i];
			pPainter->save();
			pPainter->translate(*vecPoints[i]);
			pPainter->rotate(dTotalAngle);
			pPainter->translate(-dTransScale, +4.*0.1*0.5*g_dFontSize);
			if(flip_text(dTotalAngle))
			{
				pPainter->translate(dTransScale*0.5, -8.*0.1*0.5*g_dFontSize);
				pPainter->rotate(180.);
			}
			pPainter->drawText(QPointF(0.,0.), QString::fromWCharArray(ostrAngle.str().c_str()));
			pPainter->restore();

		}
	}


	// Ewald sphere
	if(m_bShowEwaldSphere)
	{
		t_real dKLen = m_bEwaldAroundKi ? lineKi.length() : lineKf.length();

		QPen penCyan(Qt::darkCyan);
		penCyan.setWidthF(g_dFontSize*0.1);

		pPainter->setPen(penCyan);
		pPainter->drawEllipse(ptKiKf, dKLen, dKLen);
	}
}


t_real ScatteringTriangle::GetKi() const
{
	QPointF ptKiQ = mapFromItem(m_pNodeKiQ.get(), 0, 0);
	QPointF ptKiKf = mapFromItem(m_pNodeKiKf.get(), 0, 0);

	QLineF lineKi(ptKiQ, ptKiKf);
	const t_real dKi = lineKi.length()/m_dScaleFactor;
	return dKi;
}


t_real ScatteringTriangle::GetKf() const
{
	QPointF ptKfQ = mapFromItem(m_pNodeKfQ.get(), 0, 0);
	QPointF ptKiKf = mapFromItem(m_pNodeKiKf.get(), 0, 0);

	QLineF lineKf(ptKiKf, ptKfQ);
	const t_real dKf = lineKf.length()/m_dScaleFactor;
	return dKf;
}


t_real ScatteringTriangle::GetE() const
{
	const t_real dKi = GetKi();
	const t_real dKf = GetKf();
	const t_real dE = tl::get_energy_transfer(dKi/angs, dKf/angs) / meV;
	return dE;
}


t_real ScatteringTriangle::GetQ() const
{
  	QPointF ptKiQ = mapFromItem(m_pNodeKiQ.get(), 0, 0) * m_dZoom;
	QPointF ptKfQ = mapFromItem(m_pNodeKfQ.get(), 0, 0) * m_dZoom;

	QLineF lineQ(ptKiQ, ptKfQ);
	const t_real dQ = lineQ.length()/m_dScaleFactor/m_dZoom;
	return dQ;
}


t_real ScatteringTriangle::Getq() const
{
	QPointF ptKfQ = mapFromItem(m_pNodeKfQ.get(), 0, 0) * m_dZoom;
	QPointF ptGq = mapFromItem(m_pNodeGq.get(), 0, 0) * m_dZoom;

	QLineF lineq(ptKfQ, ptGq);
	const t_real dq = lineq.length()/m_dScaleFactor/m_dZoom;
	return dq;
}


t_real ScatteringTriangle::GetAngleQVec0() const
{
	t_vec vecQ = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0));
	vecQ /= ublas::norm_2(vecQ);
	vecQ = -vecQ;

	// TODO: Q is in rlu! check angle for non-cubic case!
	return tl::vec_angle(vecQ);
}


t_real ScatteringTriangle::GetAngleKiQ(bool bPosSense) const
{
	try
	{
		t_real kiQ = tl::get_angle_ki_Q(GetKi()/angs, GetKf()/angs, GetQ()/angs,
			bPosSense, false) / rads;

		t_real dTT = GetTwoTheta(bPosSense);
		if(dTT > tl::get_pi<t_real>())
			kiQ = -kiQ;

		return kiQ;
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
		return 0.;
	}
}


t_real ScatteringTriangle::GetAngleKfQ(bool bPosSense) const
{
	try
	{
		t_real kfQ = tl::get_angle_kf_Q(GetKi()/angs, GetKf()/angs, GetQ()/angs,
			bPosSense, false) / rads;

		t_real dTT = GetTwoTheta(bPosSense);
		if(dTT > tl::get_pi<t_real>())
			kfQ = -kfQ;

		return kfQ;
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
		return 0.;
	}
}


t_real ScatteringTriangle::GetTheta(bool bPosSense) const
{
	t_vec vecKi = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	vecKi /= ublas::norm_2(vecKi);

	t_real dTh = tl::vec_angle(vecKi) - tl::get_pi<t_real>()/2.;
	dTh = std::fmod(dTh, 2.*tl::get_pi<t_real>());

	//tl::log_info("theta: ", dTh/M_PI*180.);
	return dTh;
}


t_real ScatteringTriangle::GetTwoTheta(bool bPosSense) const
{
	t_vec vecKi = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	t_vec vecKf = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));

	vecKi /= ublas::norm_2(vecKi);
	vecKf /= ublas::norm_2(vecKf);

	t_real dTT = tl::vec_angle(vecKi) - tl::vec_angle(vecKf);
	if(!bPosSense)
		dTT = -dTT;
	if(dTT < 0.)
		dTT += 2.*tl::get_pi<t_real>();
	dTT = std::fmod(dTT, 2.*tl::get_pi<t_real>());

	//tl::log_info("2theta: ", dTT/M_PI*180.);
	return dTT;
}


t_real ScatteringTriangle::GetMonoTwoTheta(t_real dMonoD, bool bPosSense) const
{
	t_vec vecKi = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(), 0, 0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(), 0, 0));
	t_real dKi = ublas::norm_2(vecKi) / m_dScaleFactor;

	return tl::get_mono_twotheta(dKi/angs, dMonoD*angs, bPosSense) / rads;
}


t_real ScatteringTriangle::GetAnaTwoTheta(t_real dAnaD, bool bPosSense) const
{
	t_vec vecKf = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(), 0, 0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(), 0, 0));
	t_real dKf = ublas::norm_2(vecKf) / m_dScaleFactor;

	return tl::get_mono_twotheta(dKf/angs, dAnaD*angs, bPosSense) / rads;
}


void ScatteringTriangle::SetAnaTwoTheta(t_real dTT, t_real dAnaD)
{
	dTT = std::fabs(dTT);
	t_real dKf  = tl::get_pi<t_real>() / std::sin(dTT/2.) / dAnaD;
	dKf *= m_dScaleFactor;

	const t_vec vecNodeKiKf = qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	const t_vec vecNodeKfQ = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0));
	t_vec vecKf = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0))
		- vecNodeKiKf;

	vecKf /= ublas::norm_2(vecKf);
	t_vec vecKf_new = vecKf * dKf;

	m_bUpdate = m_bReady = false;
	m_pNodeKfQ->setPos(vec_to_qpoint(vecNodeKiKf + vecKf_new));
	m_bUpdate = m_bReady = true;

	nodeMoved(m_pNodeKfQ.get());
}


void ScatteringTriangle::SetMonoTwoTheta(t_real dTT, t_real dMonoD)
{
	dTT = std::fabs(dTT);
	t_real dKi = tl::get_pi<t_real>() / std::sin(dTT/2.) / dMonoD;
	dKi *= m_dScaleFactor;

	const t_vec vecNodeKiKf = qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	const t_vec vecNodeKiQ = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0));
	const t_vec vecNodeKfQ = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0));
	const t_vec vecKi_old = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	const t_vec vecKf = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0))
		- vecNodeKiKf;

	t_vec vecKi = vecKi_old;
	vecKi /= ublas::norm_2(vecKi);
	t_vec vecKi_new = vecKi * dKi;

	t_vec vecNodeKiKf_new = vecNodeKiQ - vecKi_new;

	m_bUpdate = m_bReady = false;
	m_pNodeKiKf->setPos(vec_to_qpoint(vecNodeKiKf_new));
	m_pNodeKfQ->setPos(vec_to_qpoint(vecNodeKiKf_new + vecKf));
	m_bReady = true;

	// don't call update twice!
	nodeMoved(m_pNodeKiKf.get());
	m_bUpdate = true;
	nodeMoved(m_pNodeKfQ.get());
}


void ScatteringTriangle::SetTwoTheta(t_real dTT)
{
	const t_vec vecNodeKiKf = qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	const t_vec vecKi = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	const t_vec vecKf = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0))
		- vecNodeKiKf;

	t_vec vecKf_new = ublas::prod(tl::rotation_matrix_2d(-dTT), vecKi);

	vecKf_new /= ublas::norm_2(vecKf_new);
	vecKf_new *= ublas::norm_2(vecKf);

	m_bUpdate = m_bReady = false;
	m_pNodeKfQ->setPos(vec_to_qpoint(vecNodeKiKf + vecKf_new));
	m_bUpdate = m_bReady = true;

	nodeMoved(m_pNodeKfQ.get());
}


void ScatteringTriangle::RotateKiVec0To(bool bSense, t_real dAngle)
{
	t_real dCurAngle = -GetTheta(bSense) + tl::get_pi<t_real>()/2.;
	//tl::log_info("kivec0 old: ", dCurAngle/M_PI*180., "new: ", dAngle/M_PI*180.);

	t_vec vecNodeKiKf = qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	t_vec vecNodeKfQ = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0));

	t_mat matRot = tl::rotation_matrix_2d(dCurAngle - dAngle);
	vecNodeKiKf = ublas::prod(matRot, vecNodeKiKf);
	vecNodeKfQ = ublas::prod(matRot, vecNodeKfQ);

	m_bUpdate = m_bReady = false;
	m_pNodeKiKf->setPos(vec_to_qpoint(vecNodeKiKf));
	m_pNodeKfQ->setPos(vec_to_qpoint(vecNodeKfQ));
	m_bReady = true;

	// don't call update twice
	nodeMoved(m_pNodeKiKf.get());
	m_bUpdate = true;
	nodeMoved(m_pNodeKfQ.get());
}


template<class T=double>
static inline std::string print_complex(const std::complex<T>& c)
{
	T r = c.real();
	T i = c.imag();

	tl::set_eps_0(r, g_dEps);
	tl::set_eps_0(i, g_dEps);

	std::ostringstream ostr;
	ostr.precision(g_iPrec);

	// imaginary part 0?
	if(tl::float_equal<T>(i, 0))
		ostr << r;
	else
		ostr << "(" << r << (i>=0. ? " +" : " -") << std::abs(i) << " i)";
	return ostr.str();
}


void ScatteringTriangle::CalcPeaks(const xtl::LatticeCommon<t_real>& recipcommon, bool bIsPowder)
{
	ClearPeaks();
	m_vecPowderLines.clear();
	m_vecPowderLineWidths.clear();
	m_kdLattice.Unload();

	m_lattice = recipcommon.lattice;
	m_recip = recipcommon.recip;
	m_plane = recipcommon.plane;
	m_matPlane = recipcommon.matPlane;
	m_matPlaneRlu = recipcommon.matPlaneRLU;
	m_matPlane_inv = recipcommon.matPlane_inv;

	tl::Powder<int, t_real_glob> powder;
	powder.SetRecipLattice(&m_recip);

	m_bz.SetEpsilon(g_dEps);
	m_bz.SetMaxNN(g_iMaxNN);
	m_bz3.SetEpsilon(g_dEps);
	m_bz3.SetMaxNN(g_iMaxNN);

	// -------------------------------------------------------------------------
	// central peak for BZ calculation
#ifdef CALC_BZ_AROUND_ZERO
	ublas::vector<int> veciCent = tl::make_vec({0.,0.,0.});
#else
	std::vector<std::tuple<int, int>> vecPeaksToTry =
	{
		std::make_tuple(1, 2), std::make_tuple(1, 3), std::make_tuple(1, 4), std::make_tuple(1, 5), std::make_tuple(1, 6),
		std::make_tuple(2, 1), std::make_tuple(2, 3), std::make_tuple(2, 4), std::make_tuple(2, 5), std::make_tuple(2, 6),
		std::make_tuple(3, 1), std::make_tuple(3, 2), std::make_tuple(3, 4), std::make_tuple(3, 5), std::make_tuple(3, 6),
	};

	ublas::vector<int> veciCent;
	for(const std::tuple<int, int>& tup : vecPeaksToTry)
	{
		veciCent.clear();

		t_vec vecdCent = std::get<0>(tup) * recipcommon.dir0RLU +
			std::get<1>(tup) * recipcommon.dir1RLU;
		veciCent = tl::convert_vec<t_real, int>(vecdCent);
		if(recipcommon.pSpaceGroup &&
			!recipcommon.pSpaceGroup->HasGenReflection(veciCent[0], veciCent[1], veciCent[2]))
			continue;
		break;
	}

	if(!veciCent.size())
		veciCent = tl::make_vec({0.,0.,0.});
#endif
	// -------------------------------------------------------------------------


	static const std::string strAA = tl::get_spec_char_utf8("AA") +
		tl::get_spec_char_utf8("sup-") +
		tl::get_spec_char_utf8("sup1");
	static const std::string strSup2 = tl::get_spec_char_utf8("sup2");

	static const QColor colPeakAllowed = Qt::red;
	static const QColor colPeakForbidden(0xaa, 0xaa, 0xaa);
	static const QColor colPeakOrigin = Qt::darkGreen;

	std::list<std::vector<t_real>> lstPeaksForKd;
	t_real dMinF = std::numeric_limits<t_real>::max(), dMaxF = -1.;

	const int iMaxNN = g_iMaxNN <= 4 ? 2 : g_iMaxNN-2;	// TODO

	// iterate over all bragg peaks
	const int iMaxPeaks = bIsPowder ? m_iMaxPeaks/2 : m_iMaxPeaks;
	for(int ih = -iMaxPeaks; ih <= iMaxPeaks; ++ih)
	for(int ik = -iMaxPeaks; ik <= iMaxPeaks; ++ik)
	for(int il = -iMaxPeaks; il <= iMaxPeaks; ++il)
	{
		const t_real h = t_real(ih);
		const t_real k = t_real(ik);
		const t_real l = t_real(il);
		const t_vec vecPeakHKL = tl::make_vec<t_vec>({ h, k, l });

		bool bHasRefl = true;
		bool bHasGenRefl = true;

		if(recipcommon.pSpaceGroup)
		{
			bHasRefl = recipcommon.pSpaceGroup->HasReflection(ih, ik, il);
			bHasGenRefl = recipcommon.pSpaceGroup->HasGenReflection(ih, ik, il);
		}

		if(!bHasGenRefl)
			continue;

		t_vec vecPeak = m_recip.GetPos(h, k, l);

		// add peak in 1/A and rlu units (only 1/A vectors are used for kd calculation)
		lstPeaksForKd.push_back(std::vector<t_real>
			{ vecPeak[0],vecPeak[1],vecPeak[2], h,k,l/*, dF*/ });

		// add peaks for 3d calculation of 1st BZ
		if(g_b3dBZ && bHasGenRefl)
		{
			if(ih==veciCent[0] && ik==veciCent[1] && il==veciCent[2])
				m_bz3.SetCentralReflex(vecPeak, &vecPeakHKL);
			else if(std::abs(ih-veciCent[0]) <= iMaxNN &&
				std::abs(ik-veciCent[1]) <= iMaxNN &&
				std::abs(il-veciCent[2]) <= iMaxNN)
				m_bz3.AddReflex(vecPeak, &vecPeakHKL);
		}

		t_real dDist = 0.;
		t_vec vecDropped = m_plane.GetDroppedPerp(vecPeak, &dDist);
		bool bInPlane = tl::float_equal<t_real>(dDist, 0., m_dPlaneDistTolerance);

		// --------------------------------------------------------------------
		// structure factors
		std::complex<t_real> cF(-1., -1.);
		t_real dF = -1., dFsq = -1.;

		if(bHasRefl && recipcommon.CanCalcStructFact() && (bInPlane || bIsPowder))
		{
			std::tie(cF, dF, dFsq) =
				recipcommon.GetStructFact(vecPeak);

			//dFsq *= tl::lorentz_factor(dAngle);
			tl::set_eps_0(dFsq, g_dEpsGfx);

			tl::set_eps_0(dF, g_dEpsGfx);
			dMinF = std::min(dF, dMinF);
			dMaxF = std::max(dF, dMaxF);
		}
		// --------------------------------------------------------------------

		t_vec vecCoord = ublas::prod(m_matPlane_inv, vecDropped);
		t_real dX = vecCoord[0];
		t_real dY = -vecCoord[1];

		if(bIsPowder && (bHasRefl || m_bShowAllPeaks))
			powder.AddPeak(ih, ik, il, dF);

		// ignore peaks that are not in the scattering plane
		if(!bInPlane)
			continue;

		// (000), i.e. direct beam, also needed for powder
		if(bIsPowder && (ih != 0 || ik != 0 || il != 0))
			continue;

		if(bHasRefl || m_bShowAllPeaks)
		{
			RecipPeak *pPeak = new RecipPeak();
			if(ih==0 && ik==0 && il==0)
				pPeak->SetColor(bHasRefl ? colPeakOrigin : colPeakForbidden);
			else
				pPeak->SetColor(bHasRefl ? colPeakAllowed : colPeakForbidden);

			pPeak->setPos(dX * m_dScaleFactor, dY * m_dScaleFactor);
			pPeak->SetRadius(dF >= 0. ? dF : 1.);
			pPeak->setData(TRIANGLE_NODE_TYPE_KEY, NODE_BRAGG);

			std::ostringstream ostrLabel, ostrTip;
			ostrLabel.precision(g_iPrecGfx);
			ostrTip.precision(g_iPrec);

			ostrLabel << "(" << ih << " " << ik << " " << il << ")";
			ostrTip << "G = (" << ih << " " << ik << " " << il << ") rlu";

			tl::set_eps_0(vecPeak, g_dEps);
			ostrTip << "\nG = (" << vecPeak[0] << ", "
				<< vecPeak[1] << ", "
				<< vecPeak[2] << ") " << strAA;

			if(dFsq > -1.)
			{
				if(g_bShowFsq)
					ostrLabel << "\nS = " << dFsq;
				else
					ostrLabel << "\nF = " << dF;

				ostrTip << "\nF = " << print_complex<t_real>(cF) << " fm";
				ostrTip << "\nS = " << dFsq << " fm" << strSup2;
			}
			else if(!bHasRefl)
			{
				pPeak->SetPeakAllowed(0);
				//pPeak->SetRadius(0.);
				ostrTip << "\nStructurally forbidden reflection.";
			}


			pPeak->SetLabel(ostrLabel.str().c_str());
			pPeak->setToolTip(QString::fromUtf8(ostrTip.str().c_str(), ostrTip.str().length()));

			m_vecPeaks.push_back(pPeak);
			m_scene.addItem(pPeak);
		}


		// add peaks for 2d approximation of 1st BZ
		if(!g_b3dBZ && bHasGenRefl)
		{
			t_vec vecN = tl::make_vec({dX, dY});
			if(ih==veciCent[0] && ik==veciCent[1] && il==veciCent[2])
			{
				m_bz.SetCentralReflex(vecN, &vecPeakHKL);
			}
			else if(std::abs(ih-veciCent[0])<=2 && std::abs(ik-veciCent[1])<=2
				&& std::abs(il-veciCent[2])<=2)
			{
				m_bz.AddReflex(vecN, &vecPeakHKL);
			}
		}
	}  // peak iteration

	// single crystal
	if(!bIsPowder)
	{
		if(g_b3dBZ)
		{
			m_bz3.CalcBZ(get_max_threads());

			// ----------------------------------------------------------------
			// calculate points of high symmetry
			std::vector<t_vec> vecSymmDirs = {
				tl::make_vec<t_vec>({1,0,0}), tl::make_vec<t_vec>({0,1,0}), tl::make_vec<t_vec>({0,0,1}),

				tl::make_vec<t_vec>({1,1,0}), tl::make_vec<t_vec>({0,1,1}), tl::make_vec<t_vec>({1,0,1}),
				tl::make_vec<t_vec>({1,-1,0}), tl::make_vec<t_vec>({0,1,-1}), tl::make_vec<t_vec>({1,0,-1}),

				tl::make_vec<t_vec>({1,1,1}), tl::make_vec<t_vec>({1,1,-1}), tl::make_vec<t_vec>({1,-1,-1}),
				tl::make_vec<t_vec>({1,-1,1}),
			};

			for(const t_vec& vecSymmDir : vecSymmDirs)
			{
				const t_vec vecSymmDirInvA = m_recip.GetPos(vecSymmDir[0], vecSymmDir[1], vecSymmDir[2]);
				tl::Line<t_real> lineSymmDir(m_bz3.GetCentralReflex(), vecSymmDirInvA);
				std::vector<t_vec> vecSymmIntersects = m_bz3.GetIntersection(lineSymmDir);
				for(t_vec& vecSymmIntersect : vecSymmIntersects)
					m_vecBZ3SymmPts.emplace_back(std::move(vecSymmIntersect));
			}
			// ----------------------------------------------------------------

			// ----------------------------------------------------------------
			// calculate intersection with scattering plane
			tl::Plane<t_real> planeBZ3 = tl::Plane<t_real>(m_bz3.GetCentralReflex(),
				m_plane.GetNorm());

			std::tie(std::ignore, m_vecBZ3VertsUnproj) = m_bz3.GetIntersection(planeBZ3);

			for(const t_vec& _vecBZ3Vert : m_vecBZ3VertsUnproj)
			{
				t_vec vecBZ3Vert = ublas::prod(m_matPlane_inv, _vecBZ3Vert - m_bz3.GetCentralReflex());
				vecBZ3Vert.resize(2, true);
				vecBZ3Vert[1] = -vecBZ3Vert[1];

				m_vecBZ3Verts.push_back(vecBZ3Vert);
			}
			// ----------------------------------------------------------------
		}
		else
		{
			m_bz.CalcBZ();
		}

		m_kdLattice.Load(lstPeaksForKd, 3);
	}

	// single crystal peaks
	if(dMaxF >= 0.)
	{
		bool bValidStructFacts = !tl::float_equal(dMinF, dMaxF, g_dEpsGfx);
		for(RecipPeak *pPeak : m_vecPeaks)
		{
			if(bValidStructFacts)
			{
				t_real dFScale = (pPeak->GetRadius()-dMinF) / (dMaxF-dMinF);
				pPeak->SetRadius(tl::lerp(MIN_PEAK_SIZE, MAX_PEAK_SIZE, dFScale));
			}
			else
			{
				pPeak->SetRadius(DEF_PEAK_SIZE);
			}
		}
	}

	// powder lines
	if(bIsPowder)
	{
		using t_line = typename decltype(powder)::t_peak;
		m_vecPowderLines = powder.GetUniquePeaksSumF();
		m_vecPowderLineWidths.reserve(m_vecPowderLines.size());

		t_real dMinFLine = 0.;
		t_real dMaxFLine = 0.;

		if(dMaxF >= 0.)
		{
			auto minmaxiters = std::minmax_element(m_vecPowderLines.begin(), m_vecPowderLines.end(),
				[](const t_line& line1, const t_line& line2) -> bool
				{
					return std::get<4>(line1) < std::get<4>(line2);
				});
			dMinFLine = std::get<4>(*minmaxiters.first);
			dMaxFLine = std::get<4>(*minmaxiters.second);
		}

		bool bValidStructFacts = !tl::float_equal(dMinFLine, dMaxFLine, g_dEpsGfx);
		for(t_line& line : m_vecPowderLines)
		{
			if(bValidStructFacts)
			{
				t_real dFScale = (std::get<4>(line)-dMinFLine) / (dMaxFLine-dMinFLine);
				m_vecPowderLineWidths.push_back(tl::lerp(MIN_PEAK_SIZE, MAX_PEAK_SIZE, dFScale));
			}
			else
			{
				m_vecPowderLineWidths.push_back(1.);
			}
		}
	}

	m_scene.emitAllParams();
	this->update();
}


t_vec ScatteringTriangle::GetHKLFromPlanePos(t_real x, t_real y) const
{
	if(!HasPeaks())
		return t_vec();

	t_vec vec = x*tl::get_column(m_matPlane, 0)
		+ y*tl::get_column(m_matPlane, 1);
	return m_recip.GetHKL(vec);
}


t_vec ScatteringTriangle::GetQVecPlane(bool bSmallQ) const
{
	t_vec vecQPlane;

	if(bSmallQ)
		vecQPlane = qpoint_to_vec(mapFromItem(m_pNodeGq.get(),0,0))
			- qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(), 0, 0));
	else
		vecQPlane = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
			- qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0));

	vecQPlane[1] = -vecQPlane[1];
	vecQPlane /= m_dScaleFactor;

	return vecQPlane;
}


t_vec ScatteringTriangle::GetQVec(bool bSmallQ, bool bRLU) const
{
	t_vec vecQPlane = GetQVecPlane(bSmallQ);

	t_vec vecQ;
	if(bRLU)
		vecQ = GetHKLFromPlanePos(vecQPlane[0], vecQPlane[1]);
	else
		vecQ = vecQPlane[0]*tl::get_column(m_matPlane, 0)
				+ vecQPlane[1]*tl::get_column(m_matPlane, 1);

	return vecQ;
}


t_vec ScatteringTriangle::GetKiVecPlane() const
{
	t_vec vecKi = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	vecKi[1] = -vecKi[1];
	vecKi /= m_dScaleFactor;
	return vecKi;
}


t_vec ScatteringTriangle::GetKfVecPlane() const
{
	t_vec vecKf = qpoint_to_vec(mapFromItem(m_pNodeKfQ.get(),0,0))
		- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
	vecKf[1] = -vecKf[1];
	vecKf /= m_dScaleFactor;
	return vecKf;
}


void ScatteringTriangle::ClearPeaks()
{
	m_bz.Clear();
	m_bz3.Clear();
	m_vecBZ3VertsUnproj.clear();
	m_vecBZ3Verts.clear();
	m_vecBZ3SymmPts.clear();

	for(RecipPeak*& pPeak : m_vecPeaks)
	{
		if(pPeak)
		{
			m_scene.removeItem(pPeak);
			delete pPeak;
			pPeak = nullptr;
		}
	}
	m_vecPeaks.clear();
}


std::vector<ScatteringTriangleNode*> ScatteringTriangle::GetNodes()
{
	return std::vector<ScatteringTriangleNode*>
		{ m_pNodeKiQ.get(), m_pNodeKiKf.get(),
		m_pNodeKfQ.get(), m_pNodeGq.get() };
}


std::vector<std::string> ScatteringTriangle::GetNodeNames() const
{
	return std::vector<std::string>{ "kiQ", "kikf", "kfQ", "Gq" };
}


static std::pair<bool, QPointF>
get_nearest_elastic_kikf_pos(const QPointF& ptKiKf, const QPointF& ptKiQ, const QPointF& ptKfQ)
{
	std::pair<bool, QPointF> pairRet;
	bool& bOk = pairRet.first;
	QPointF& ptRet = pairRet.second;

	t_vec vecKiQ = qpoint_to_vec(ptKiQ);
	t_vec vecKfQ = qpoint_to_vec(ptKfQ);
	t_vec vecKiKf = qpoint_to_vec(ptKiKf);
	t_vec vecQ = vecKfQ - vecKiQ;

	t_vec vecQperp = tl::make_vec({vecQ[1], -vecQ[0]});
	t_vec vecQMid = vecKiQ + vecQ*0.5;

	tl::Line<t_real> lineQMidPerp(vecQMid, vecQperp);
	t_vec vecDrop = lineQMidPerp.GetDroppedPerp(vecKiKf);
	ptRet = vec_to_qpoint(vecDrop);

	bOk = true;
	return pairRet;
}


static std::tuple<bool, t_real, QPointF>
get_nearest_node(const QPointF& pt,
	const QGraphicsItem* pCurItem, const QList<QGraphicsItem*>& nodes, t_real dFactor,
	const std::vector<typename ScatteringTriangle::t_powderline>* pPowderLines=nullptr,
	const tl::Lattice<t_real>* pRecip=nullptr)
{
	if(nodes.size()==0)
		return std::tuple<bool, t_real, QPointF>(0, 0., QPointF());

	t_real dMinLen = std::numeric_limits<t_real>::max();
	int iMinIdx = -1;
	bool bHasPowderPeak = false;

	// Bragg peaks
	for(int iNode=0; iNode<nodes.size(); ++iNode)
	{
		const QGraphicsItem *pNode = nodes[iNode];
		if(pNode == pCurItem || pNode->data(TRIANGLE_NODE_TYPE_KEY)!=NODE_BRAGG)
			continue;

		t_real dLen = QLineF(pt, pNode->scenePos()).length();
		if(dLen < dMinLen)
		{
			dMinLen = dLen;
			iMinIdx = iNode;
		}
	}

	const QGraphicsItem* pNodeOrigin = nullptr;
	for(const QGraphicsItem* pNode : nodes)
	{
		if(pNode->data(TRIANGLE_NODE_TYPE_KEY) == NODE_KIQ)
		{
			pNodeOrigin = pNode;
			break;
		}
	}



	// Powder peaks
	QPointF ptPowder;
	if(pNodeOrigin && pPowderLines && pRecip)
	{
		t_vec vecOrigin = qpoint_to_vec(pNodeOrigin->scenePos());
		t_vec vecPt = qpoint_to_vec(pt);
		t_vec vecOriginPt = vecPt-vecOrigin;
		const t_real dDistToOrigin = ublas::norm_2(vecOriginPt);
		vecOriginPt /= dDistToOrigin;

		for(const typename ScatteringTriangle::t_powderline& powderpeak : *pPowderLines)
		{
			const int ih = std::get<0>(powderpeak);
			const int ik = std::get<1>(powderpeak);
			const int il = std::get<2>(powderpeak);
			if(ih==0 && ik==0 && il==0) continue;

			t_vec vec = pRecip->GetPos(t_real(ih), t_real(ik), t_real(il));
			t_real drad = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
			drad *= dFactor;

			if(std::fabs(drad-dDistToOrigin) < dMinLen)
			{
				bHasPowderPeak = true;
				dMinLen = std::fabs(drad-dDistToOrigin);

				t_vec vecPowder = vecOrigin + vecOriginPt*drad;
				ptPowder = vec_to_qpoint(vecPowder);
			}
		}
	}



	if(bHasPowderPeak)
	{
		return std::tuple<bool, t_real, QPointF>(1, dMinLen, ptPowder);
	}
	else
	{
		if(iMinIdx < 0)
			return std::tuple<bool, t_real, QPointF>(0, 0., QPointF());

		return std::tuple<bool, t_real, QPointF>(1, dMinLen, nodes[iMinIdx]->pos());
	}
}


// snap pNode to a peak near pNodeOrg
void ScatteringTriangle::SnapToNearestPeak(ScatteringTriangleNode* pNode,
	const ScatteringTriangleNode* pNodeOrg)
{
	if(!pNode) return;
	if(!pNodeOrg) pNodeOrg = pNode;
	if(!HasPeaks()) return;

	std::tuple<bool, t_real, QPointF> tupNearest =
		get_nearest_node(pNodeOrg->pos(), pNode, m_scene.items(),
			GetScaleFactor(), &GetPowder(), &GetRecipLattice());

	if(std::get<0>(tupNearest))
		pNode->setPos(std::get<2>(tupNearest));
}


bool ScatteringTriangle::KeepAbsKiKf(t_real dQx, t_real dQy)
{
	try
	{
		t_vec vecCurKfQ = tl::make_vec({
			m_pNodeKfQ->scenePos().x(),
			m_pNodeKfQ->scenePos().y() });

		t_vec vecNewKfQ = tl::make_vec({
			m_pNodeKfQ->scenePos().x() + dQx,
			m_pNodeKfQ->scenePos().y() + dQy });

		t_vec vecKiQ = qpoint_to_vec(mapFromItem(m_pNodeKiQ.get(),0,0));
		t_vec vecKi = vecKiQ
			- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));

		t_vec vecKf = vecCurKfQ
			- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));
		t_vec vecNewKf = vecNewKfQ
			- qpoint_to_vec(mapFromItem(m_pNodeKiKf.get(),0,0));

		t_vec vecCurQ = vecKiQ - vecCurKfQ;
		t_vec vecNewQ = vecKiQ - vecNewKfQ;

		t_real dKi = ublas::norm_2(vecKi)/m_dScaleFactor;
		t_real dKf = ublas::norm_2(vecKf)/m_dScaleFactor;
		t_real dNewQ = ublas::norm_2(vecNewQ)/m_dScaleFactor;

		bool bPosSense = m_scene.GetSampleSense();
		t_real dAngKiQ = tl::get_angle_ki_Q(dKi/angs, dKf/angs, dNewQ/angs,
			bPosSense, false) / tl::radians;

		t_mat matRot = tl::rotation_matrix_2d(-dAngKiQ);
		vecKi = ublas::prod(matRot, vecNewQ) / ublas::norm_2(vecNewQ);
		vecKi *= dKi * m_dScaleFactor;
		t_vec vecPtKiKf = vecKiQ - vecKi;

		m_bReady = false;
		m_pNodeKiKf->setPos(vec_to_qpoint(vecPtKiKf));
		m_bReady = true;

		return 1;	// allowed
	}
	catch(const std::exception&)
	{
		return 0;	// not allowed
	}
}


// --------------------------------------------------------------------------------



ScatteringTriangleScene::ScatteringTriangleScene(QObject* pParent)
	: QGraphicsScene(pParent),
	m_pTri(new ScatteringTriangle(*this))
{
	this->addItem(m_pTri.get());
}


ScatteringTriangleScene::~ScatteringTriangleScene()
{}


void ScatteringTriangleScene::SetDs(t_real dMonoD, t_real dAnaD)
{
	m_dMonoD = dMonoD;
	m_dAnaD = dAnaD;

	emitUpdate();
}


void ScatteringTriangleScene::emitUpdate()
{
	if(!m_pTri || !m_pTri->IsReady() || m_bDontEmitChange)
		return;

	TriangleOptions opts;
	opts.bChangedMonoTwoTheta = true;
	opts.bChangedAnaTwoTheta = true;
	opts.bChangedTwoTheta = true;
	opts.bChangedTheta = true;

	try
	{
		opts.dTwoTheta = m_pTri->GetTwoTheta(m_bSamplePosSense);
		opts.dTheta = m_pTri->GetTheta(m_bSamplePosSense);
	}
	catch(const std::exception&)
	{
		opts.dTwoTheta = std::nan("");
	}

	try
	{
		opts.dMonoTwoTheta = m_pTri->GetMonoTwoTheta(m_dMonoD, m_bMonoPosSense);
	}
	catch(const std::exception&)
	{
		opts.dMonoTwoTheta = std::nan("");
	}

	try
	{
		opts.dAnaTwoTheta = m_pTri->GetAnaTwoTheta(m_dAnaD, m_bAnaPosSense);
	}
	catch(const std::exception&)
	{
		opts.dAnaTwoTheta = std::nan("");
	}

	//tl::log_debug("triangle: triangleChanged");
	emit triangleChanged(opts);
}


void ScatteringTriangleScene::emitAllParams()
{
	if(!m_pTri || !m_pTri->IsReady() /*|| m_bDontEmitChange*/) // emit even with m_bDontEmitChange because of Q vector in real space
		return;

	RecipParams parms;
	parms.d2Theta = m_pTri->GetTwoTheta(m_bSamplePosSense);
	parms.dTheta = m_pTri->GetTheta(m_bSamplePosSense);
	parms.dE = m_pTri->GetE();
	parms.dQ = m_pTri->GetQ();
	parms.dq = m_pTri->Getq();
	parms.dki = m_pTri->GetKi();
	parms.dkf = m_pTri->GetKf();
	parms.dKiQ = m_pTri->GetAngleKiQ(m_bSamplePosSense);
	parms.dKfQ = m_pTri->GetAngleKfQ(m_bSamplePosSense);
	parms.dAngleQVec0 = m_pTri->GetAngleQVec0();

	t_vec vecQ = m_pTri->GetQVec(0,0);
	t_vec vecQrlu = m_pTri->GetQVec(0,1);
	t_vec vecq = m_pTri->GetQVec(1,0);
	t_vec vecqrlu = m_pTri->GetQVec(1,1);
	t_vec vecG = vecQ - vecq;
	t_vec vecGrlu = vecQrlu - vecqrlu;

	const t_mat& matPlane = m_pTri->GetPlane(true);
	t_vec vec0 = tl::get_column(matPlane, 0);
	t_vec vec1 = tl::get_column(matPlane, 1);
	t_vec vecUp = tl::get_column(matPlane, 2);

	tl::set_eps_0(vecQ, g_dEps); tl::set_eps_0(vecQrlu, g_dEps);
	tl::set_eps_0(vecq, g_dEps); tl::set_eps_0(vecqrlu, g_dEps);
	tl::set_eps_0(vecG, g_dEps); tl::set_eps_0(vecGrlu, g_dEps);


	for(unsigned i = 0; i < 3; ++i)
	{
		parms.Q[i] = vecQ.size() ? vecQ[i] : 0.;
		parms.Q_rlu[i] = vecQrlu.size() ? vecQrlu[i] : 0.;

		parms.q[i] = vecq.size() ? vecq[i] : 0.;
		parms.q_rlu[i] = vecqrlu.size() ? vecqrlu[i] : 0.;

		parms.G[i] = vecG.size() ? vecG[i] : 0.;
		parms.G_rlu[i] = vecGrlu.size() ? vecGrlu[i] : 0.;

		parms.orient_0[i] = vec0.size() ? vec0[i] : 0.;
		parms.orient_1[i] = vec1.size() ? vec1[i] : 0.;
		parms.orient_up[i] = vecUp.size() ? vecUp[i] : 0.;
	}


	// nearest node (exact G)
	if(vecQrlu.size())
	{
		parms.G_rlu_accurate[0] = parms.G_rlu_accurate[1] = parms.G_rlu_accurate[2] = 0.;
		const tl::Kd<t_real>& kd = m_pTri->GetKdLattice();
		t_vec vecHKLinvA = m_pTri->GetRecipLattice().GetPos(-vecQrlu[0], -vecQrlu[1], -vecQrlu[2]);

		if(kd.GetRootNode())
		{
			std::vector<t_real> stdvecHKL{vecHKLinvA[0], vecHKLinvA[1], vecHKLinvA[2]};
			const std::vector<t_real>* pvecNearest = &kd.GetNearestNode(stdvecHKL);

			if(pvecNearest)
			{
				parms.G_rlu_accurate[0] = (*pvecNearest)[3];
				parms.G_rlu_accurate[1] = (*pvecNearest)[4];
				parms.G_rlu_accurate[2] = (*pvecNearest)[5];
			}
		}
	}

	CheckForSpurions();

	emit paramsChanged(parms);
}


// check for spurions
void ScatteringTriangleScene::CheckForSpurions()
{
	typedef tl::t_energy_si<t_real> energy;

	const t_vec vecq = m_pTri->GetQVecPlane(1);
	const t_vec vecKi = m_pTri->GetKiVecPlane();
	const t_vec vecKf = m_pTri->GetKfVecPlane();
	energy E = m_pTri->GetE() * meV;
	energy Ei = tl::k2E(m_pTri->GetKi()/angs);
	energy Ef = tl::k2E(m_pTri->GetKf()/angs);

	// elastic currat-axe spurions
	tl::ElasticSpurion spuris = tl::check_elastic_spurion(vecKi, vecKf, vecq);

	// inelastic higher-order spurions
	std::vector<tl::InelasticSpurion<t_real>> vecInelCKI = tl::check_inelastic_spurions(1, Ei, Ef, E, 5);
	std::vector<tl::InelasticSpurion<t_real>> vecInelCKF = tl::check_inelastic_spurions(0, Ei, Ef, E, 5);

	emit spurionInfo(spuris, vecInelCKI, vecInelCKF);
}


void ScatteringTriangleScene::tasChanged(const TriangleOptions& opts)
{
	if(!m_pTri || !m_pTri->IsReady())
		return;

	m_bDontEmitChange = true;

	if(opts.bChangedMonoTwoTheta)
	{
		m_pTri->SetMonoTwoTheta(opts.dMonoTwoTheta, m_dMonoD);
	}
	if(opts.bChangedAnaTwoTheta)
	{
		m_pTri->SetAnaTwoTheta(opts.dAnaTwoTheta, m_dAnaD);
	}

	if(opts.bChangedTwoTheta)
	{
		m_pTri->SetTwoTheta(m_bSamplePosSense ? opts.dTwoTheta : -opts.dTwoTheta);
	}

	if(opts.bChangedAngleKiVec0)
	{
		m_pTri->RotateKiVec0To(m_bSamplePosSense, opts.dAngleKiVec0);
	}

	m_bDontEmitChange = false;
}


void ScatteringTriangleScene::SetSampleSense(bool bPos)
{
	t_real tt = m_pTri->GetTwoTheta(m_bSamplePosSense);
	m_bSamplePosSense = bPos;
	m_pTri->SetTwoTheta(m_bSamplePosSense ? tt : -tt);
}


void ScatteringTriangleScene::SetMonoSense(bool bPos)
{
	m_bMonoPosSense = bPos;
	emitUpdate();
	emitAllParams();
}


void ScatteringTriangleScene::SetAnaSense(bool bPos)
{
	m_bAnaPosSense = bPos;
	emitUpdate();
	emitAllParams();
}


void ScatteringTriangleScene::scaleChanged(t_real dTotalScale)
{
	if(!m_pTri) return;
	m_pTri->SetZoom(dTotalScale);
}


void ScatteringTriangleScene::mousePressEvent(QGraphicsSceneMouseEvent *pEvt)
{
	m_bMousePressed = true;
	QGraphicsScene::mousePressEvent(pEvt);
	emit nodeEvent(1);
}


void ScatteringTriangleScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt)
{
	m_bMousePressed = false;
	QGraphicsScene::mouseReleaseEvent(pEvt);
	emit nodeEvent(0);
}


void ScatteringTriangleScene::setSnapq(bool bSnap)
{
	m_bSnapq = bSnap;

	if(m_bSnapq && m_pTri)
		m_pTri->SnapToNearestPeak(m_pTri->GetNodeGq(), m_pTri->GetNodeKfQ());
}


void ScatteringTriangleScene::drawBackground(QPainter* pPainter, const QRectF& rect)
{
	QGraphicsScene::drawBackground(pPainter, rect);
}


void ScatteringTriangleScene::mouseMoveEvent(QGraphicsSceneMouseEvent *pEvt)
{
	bool bHandled = false;
	bool bAllowed = true;

	// tooltip
	if(m_pTri)
	{
		const t_real dX = pEvt->scenePos().x()/m_pTri->GetScaleFactor();
		const t_real dY = -pEvt->scenePos().y()/m_pTri->GetScaleFactor();
		t_vec vecHKL = m_pTri->GetHKLFromPlanePos(dX, dY);
		tl::set_eps_0(vecHKL, g_dEps);

		if(vecHKL.size() == 3)
		{
			const std::vector<t_real>* pvecNearest = nullptr;

			const tl::Kd<t_real>& kd = m_pTri->GetKdLattice();
			const tl::Lattice<t_real>& recip = m_pTri->GetRecipLattice();
			t_vec vecHKLinvA = recip.GetPos(vecHKL[0], vecHKL[1], vecHKL[2]);

			if(kd.GetRootNode())
			{
				std::vector<t_real> stdvecHKL{vecHKLinvA[0], vecHKLinvA[1], vecHKLinvA[2]};
				pvecNearest = &kd.GetNearestNode(stdvecHKL);
				if(pvecNearest->size() < 6)
				{
					pvecNearest = nullptr;
					tl::log_warn("Invalid BZ node.");
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
			const int iNodeType = pCurItem->data(TRIANGLE_NODE_TYPE_KEY).toInt();

			if(m_bSnap || (m_bSnapq && iNodeType == NODE_q))
			{
				QList<QGraphicsItem*> nodes = items();
				std::tuple<bool, t_real, QPointF> tupNearest =
					get_nearest_node(pEvt->scenePos(), pCurItem, nodes,
						m_pTri->GetScaleFactor(), &m_pTri->GetPowder(),
						&m_pTri->GetRecipLattice());

				if(std::get<0>(tupNearest))
				{
					pCurItem->setPos(std::get<2>(tupNearest));
					bHandled = true;
				}
			}
			else if(iNodeType == NODE_KIKF && m_bSnapKiKfToElastic)
			{
				std::pair<bool, QPointF> pairNearest =
					get_nearest_elastic_kikf_pos(pEvt->scenePos() /*m_pTri->GetNodeKiKf()->scenePos()*/,
						m_pTri->GetNodeKiQ()->scenePos(),
						m_pTri->GetNodeKfQ()->scenePos());

				if(pairNearest.first)
				{
					pCurItem->setPos(pairNearest.second);
					bHandled = true;
				}
			}

			// TODO: case for both snapping and abs ki kf fixed
			else if(iNodeType == NODE_Q && m_bKeepAbsKiKf)
			{
				t_real dX = pEvt->scenePos().x()-pEvt->lastScenePos().x();
				t_real dY = pEvt->scenePos().y()-pEvt->lastScenePos().y();

				bAllowed = m_pTri->KeepAbsKiKf(dX, dY);
			}
		}
	}

	if(!bHandled && bAllowed)
		QGraphicsScene::mouseMoveEvent(pEvt);
}


void ScatteringTriangleScene::keyPressEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Control)
		m_bSnap = true;
	if(pEvt->key() == Qt::Key_Shift)
		m_bSnapKiKfToElastic = true;

	QGraphicsScene::keyPressEvent(pEvt);
}


void ScatteringTriangleScene::keyReleaseEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Control)
		m_bSnap = false;
	if(pEvt->key() == Qt::Key_Shift)
		m_bSnapKiKfToElastic = false;

	QGraphicsScene::keyReleaseEvent(pEvt);
}


// --------------------------------------------------------------------------------


ScatteringTriangleView::ScatteringTriangleView(QWidget* pParent)
	: QGraphicsView(pParent)
{
	setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
		QPainter::SmoothPixmapTransform | QPainter::HighQualityAntialiasing);
	setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);
	setDragMode(QGraphicsView::ScrollHandDrag);
	setMouseTracking(1);
	setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
}

ScatteringTriangleView::~ScatteringTriangleView()
{}


void ScatteringTriangleView::DoZoom(t_real_glob dDelta)
{
	const t_real dScale = std::pow(2., dDelta);
	this->scale(dScale, dScale);

	m_dTotalScale *= dScale;
	emit scaleChanged(m_dTotalScale);
}


void ScatteringTriangleView::keyPressEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Plus)
		DoZoom(0.02);
	else if(pEvt->key() == Qt::Key_Minus)
		DoZoom(-0.02);

	QGraphicsView::keyPressEvent(pEvt);
}


void ScatteringTriangleView::keyReleaseEvent(QKeyEvent *pEvt)
{
	QGraphicsView::keyReleaseEvent(pEvt);
}


void ScatteringTriangleView::wheelEvent(QWheelEvent *pEvt)
{
	const t_real dDelta = pEvt->angleDelta().y()/8. / 150.;
	DoZoom(dDelta);
}


bool ScatteringTriangleView::event(QEvent *pEvt)
{
	return QGraphicsView::event(pEvt);
}


#include "moc_scattering_triangle.cpp"
