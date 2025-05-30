/**
 * TOF layout
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2016
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

#include "tof_layout.h"
#include "tlibs/string/spec_char.h"
#include <iostream>

using t_real = t_real_glob;
using t_vec = ublas::vector<t_real>;
using t_mat = ublas::matrix<t_real>;


static inline bool flip_text(t_real _dAngle)
{
	t_real dAngle = std::fmod(_dAngle, 360.);
	return std::abs(dAngle) > 90. && std::abs(dAngle) < 270.;
}


TofLayoutNode::TofLayoutNode(TofLayout* pSupItem) : m_pParentItem(pSupItem)
{
	setFlag(QGraphicsItem::ItemSendsGeometryChanges);
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
	setCursor(Qt::CrossCursor);
}


QRectF TofLayoutNode::boundingRect() const
{
	return QRectF(-0.5*g_dFontSize, -0.5*g_dFontSize,
		1.*g_dFontSize, 1.*g_dFontSize);
}


void TofLayoutNode::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	pPainter->drawEllipse(QRectF(-2.*0.1*g_dFontSize, -2.*0.1*g_dFontSize,
		4.*0.1*g_dFontSize, 4.*0.1*g_dFontSize));
}


QVariant TofLayoutNode::itemChange(GraphicsItemChange change, const QVariant &val)
{
	QVariant var = QGraphicsItem::itemChange(change, val);

	if(change == QGraphicsItem::ItemPositionHasChanged)
		m_pParentItem->nodeMoved(this);

	return var;
}

// --------------------------------------------------------------------------------

TofLayout::TofLayout(TofLayoutScene& scene) : m_scene(scene)
{
	setFlag(QGraphicsItem::ItemIgnoresTransformations);

	m_pSrc = new TofLayoutNode(this);
	m_pSample = new TofLayoutNode(this);
	m_pDet = new TofLayoutNode(this);

	m_pSrc->setToolTip("Source");
	m_pSample->setToolTip("Sample");
	m_pDet->setToolTip("Detector");

	m_pSrc->setPos(0., m_dLenSrcSample*m_dScaleFactor);
	m_pSample->setPos(0., 0.);
	m_pDet->setPos(-m_dLenSampleDet*m_dScaleFactor, 0.);

	AllowMouseMove(true);

	scene.addItem(m_pSrc);
	scene.addItem(m_pSample);
	scene.addItem(m_pDet);

	setAcceptedMouseButtons(Qt::NoButton);
	m_bUpdate = m_bReady = true;
}


TofLayout::~TofLayout()
{
	m_bUpdate = m_bReady = false;

	delete m_pSrc;
	delete m_pSample;
	delete m_pDet;
}


void TofLayout::AllowMouseMove(bool bAllow)
{
	//m_pSrc->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pSample->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pDet->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
}


void TofLayout::nodeMoved(const TofLayoutNode *pNode)
{
	if(!m_bReady)
		return;

	// prevents recursive calling of update
	static bool bAllowUpdate = true;
	if(!bAllowUpdate)
		return;

	const t_vec vecSrc = qpoint_to_vec(mapFromItem(m_pSrc, 0, 0));
	const t_vec vecSample = qpoint_to_vec(mapFromItem(m_pSample, 0, 0));
	const t_vec vecDet = qpoint_to_vec(mapFromItem(m_pDet, 0, 0));

	bAllowUpdate = false;
	if(pNode == m_pSample)
	{
		t_real dTwoTheta = m_dTwoTheta;
		t_vec vecSrcSample = vecSample - vecSrc;

		if(m_bAllowChanges)
			m_dLenSrcSample = ublas::norm_2(vecSrcSample)/m_dScaleFactor;
		vecSrcSample /= ublas::norm_2(vecSrcSample);


		t_vec vecSampleDet =
			ublas::prod(tl::rotation_matrix_2d(-dTwoTheta), vecSrcSample);
		vecSampleDet /= ublas::norm_2(vecSampleDet);
		vecSampleDet *= m_dLenSampleDet*m_dScaleFactor;

		m_pDet->setPos(vec_to_qpoint(vecSample + vecSampleDet));
		vecSampleDet /= ublas::norm_2(vecSampleDet);


		TriangleOptions opts;
		opts.bChangedTwoTheta = true;
		opts.dTwoTheta = dTwoTheta;
		if(!m_scene.GetSampleSense())
			opts.dTwoTheta = -opts.dTwoTheta;
		m_scene.emitUpdate(opts);
	}
	if(pNode == m_pDet)
	{
		t_vec vecSampleDet = vecDet - vecSample;
		if(pNode==m_pDet && m_bAllowChanges)
			m_dLenSampleDet = ublas::norm_2(vecSampleDet)/m_dScaleFactor;
		vecSampleDet /= ublas::norm_2(vecSampleDet);

		t_vec vecSrcSample = vecSample - vecSrc;
		vecSrcSample /= ublas::norm_2(vecSrcSample);

		if(m_bAllowChanges)
		{
			m_dTwoTheta = -(tl::vec_angle(vecSampleDet) - tl::vec_angle(vecSrcSample));
			if(m_dTwoTheta < -tl::get_pi<t_real>()) m_dTwoTheta += 2.*tl::get_pi<t_real>();
			if(m_dTwoTheta > tl::get_pi<t_real>()) m_dTwoTheta -= 2.*tl::get_pi<t_real>();
		}

		TriangleOptions opts;
		opts.bChangedTwoTheta = true;
		opts.dTwoTheta = m_dTwoTheta;
		if(!m_scene.GetSampleSense())
			opts.dTwoTheta = -opts.dTwoTheta;
		m_scene.emitUpdate(opts);
	}

	bAllowUpdate = true;

	if(m_bUpdate)
	{
		this->update();
		m_scene.emitAllParams();
	}
}


QRectF TofLayout::boundingRect() const
{
	return QRectF(-100.*m_dZoom*g_dFontSize, -100.*m_dZoom*g_dFontSize,
		200.*m_dZoom*g_dFontSize, 200.*m_dZoom*g_dFontSize);
}


/**
 * gets the centre of the tof drawing
 */
QPointF TofLayout::GetGfxMid() const
{
	t_real num_nodes = 0.;
	QPointF mid(0., 0.);

	for(const TofLayoutNode* node : GetNodes())
	{
		mid += mapFromItem(node, 0, 0);
		num_nodes += 1.;
	}

	return mid / num_nodes;
}


void TofLayout::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	pPainter->setFont(g_fontGfx);
	const bool bDisplayLengths = false;

	QPointF ptSrc = mapFromItem(m_pSrc, 0, 0) * m_dZoom;
	QPointF ptSample = mapFromItem(m_pSample, 0, 0) * m_dZoom;
	QPointF ptDet = mapFromItem(m_pDet, 0, 0) * m_dZoom;

	QLineF lineKi(ptSrc, ptSample);
	QLineF lineKf(ptSample, ptDet);

	t_vec vecSrc = qpoint_to_vec(ptSrc);
	t_vec vecSample = qpoint_to_vec(ptSample);
	t_vec vecDet = qpoint_to_vec(ptDet);

	t_vec vecSrcSample = vecSample-vecSrc;
	t_vec vecSampleDet = vecDet-vecSample;


	QPen penOrig = pPainter->pen();

	QPen penRed(Qt::red);
	penRed.setWidthF(g_dFontSize*0.1);
	QPen penBlack(qApp->palette().color(QPalette::WindowText));
	penBlack.setWidthF(g_dFontSize*0.1);
	QPen penBlue(qApp->palette().color(QPalette::Link));
	penBlue.setWidthF(g_dFontSize*0.1);
	QPen penGray(Qt::gray);
	penGray.setWidthF(g_dFontSize*0.1);


	pPainter->setPen(penBlack);

	pPainter->drawLine(lineKi);
	pPainter->drawLine(lineKf);


	// draw detector bank
	QPen penArc(penGray.color());
	penArc.setWidthF(1.5*g_dFontSize*0.1);
	pPainter->setPen(penArc);

	t_real dBaseAngle = tl::r2d(tl::vec_angle(vecSrcSample));
	t_real dLenKf = lineKf.length();
	pPainter->drawArc(QRectF(ptSample.x()-dLenKf, ptSample.y()-dLenKf,
		dLenKf*2., dLenKf*2.), (m_dBeginDetArcAngle-dBaseAngle)*16., m_dDetArcAngle*16.);


	// draw choppers
	QPointF ptChopper1 = ptSrc + (ptSample-ptSrc) * 1./3.;
	QPointF ptChopper2 = ptSrc + (ptSample-ptSrc) * 2./3.;
	t_vec vecOrthKi(2);
	vecOrthKi[0] = vecSrcSample[1]; vecOrthKi[1] = -vecSrcSample[0];
	vecOrthKi /= ublas::norm_2(vecOrthKi);
	vecOrthKi *= 0.1*m_dLenSrcSample*m_dZoom*m_dScaleFactor;
	pPainter->drawLine(ptChopper1-vec_to_qpoint(vecOrthKi), ptChopper1+vec_to_qpoint(vecOrthKi));
	pPainter->drawLine(ptChopper2-vec_to_qpoint(vecOrthKi), ptChopper2+vec_to_qpoint(vecOrthKi));


	// write lengths
	QPointF ptMidKi = ptSrc + (ptSample-ptSrc)/2.;
	QPointF ptMidKf = ptSample + (ptDet-ptSample)/2.;

	if(bDisplayLengths)
	{
		std::ostringstream ostrLenKi, ostrLenKf;
		ostrLenKi.precision(g_iPrecGfx);
		ostrLenKf.precision(g_iPrecGfx);

		ostrLenKi << m_dLenSrcSample << " cm";
		ostrLenKf << m_dLenSampleDet << " cm";

		pPainter->drawText(ptMidKi, ostrLenKi.str().c_str());
		pPainter->drawText(ptMidKf, ostrLenKf.str().c_str());
	}


	// sample theta rotation
	t_vec vecRotDir =
		ublas::prod(tl::rotation_matrix_2d(-m_dTheta), vecSrcSample);
	vecRotDir /= ublas::norm_2(vecRotDir);
	vecRotDir *= m_dLenSample*m_dScaleFactor;

	QPointF ptThM = vec_to_qpoint(vecSample - vecRotDir*m_dZoom);
	QPointF ptThP = vec_to_qpoint(vecSample + vecRotDir*m_dZoom);
	QLineF lineRot = QLineF(ptThM, ptThP);

	QPen pen(penRed.color());
	pen.setWidthF(1.5*g_dFontSize*0.1);
	pPainter->setPen(pen);
	pPainter->drawLine(lineRot);


	// component names
	pPainter->setPen(penOrig);
	pPainter->save();
		pPainter->translate(vec_to_qpoint(vecSample));
		t_real dCompAngle = 180. + tl::r2d(tl::vec_angle(vecRotDir));
		pPainter->rotate(dCompAngle);
		pPainter->translate(-4., 16.*0.1*g_dFontSize);
		if(flip_text(dCompAngle))
		{
			pPainter->translate(4., -8.*0.1*g_dFontSize);
			pPainter->rotate(180.);
		}
		pPainter->drawText(QPointF(0., 0.), "S");
	pPainter->restore();

	pPainter->drawText(ptSample - vec_to_qpoint(vecSrcSample*1.1), "Src");
	pPainter->drawText(ptSample + vec_to_qpoint(vecSampleDet*1.1), "D");
	pPainter->drawText(ptChopper1 + vec_to_qpoint(vecOrthKi*2.), "ChP");
	pPainter->drawText(ptChopper2 + vec_to_qpoint(vecOrthKi*2.), "ChM");


	// dashed extended lines
	QPen penDash(Qt::DashLine);
	penDash.setColor(qApp->palette().color(QPalette::WindowText));
	penDash.setWidthF(g_dFontSize*0.1);

	pPainter->setPen(penDash);
	QLineF lineki_ext(ptSample, ptSample + (ptSample-ptSrc)/2.);
	pPainter->drawLine(lineki_ext);


	std::unique_ptr<QLineF> plineQ(nullptr);
	std::unique_ptr<QPointF> pptQ(nullptr);
	// Q vector direction visible?
	if(this->m_bRealQVisible)
	{
		const t_real &dAngleKiQ = m_dAngleKiQ;
		t_mat matRotQ = tl::rotation_matrix_2d(dAngleKiQ);
		t_vec vecKi = vecSample-vecSrc;
		t_vec vecQ = ublas::prod(matRotQ, vecKi);
		vecQ /= ublas::norm_2(vecQ);
		vecQ *= (m_dLenSrcSample + m_dLenSampleDet)/2.;	// some arbitrary length
		vecQ *= m_dScaleFactor * m_dZoom;

		pptQ.reset(new QPointF(vec_to_qpoint(vecSample + vecQ)));
		plineQ.reset(new QLineF(ptSample, *pptQ));

		pPainter->setPen(penRed);
		pPainter->drawLine(*plineQ);
		pPainter->save();
			pPainter->translate(ptSample);
			const t_real dQAngle = -plineQ->angle();
			pPainter->rotate(dQAngle);
			pPainter->translate(QPointF(plineQ->length()/2.,12.));
			if(flip_text(dQAngle))
				pPainter->rotate(180.);
			pPainter->drawText(QPointF(0,0), "Q");
		pPainter->restore();
	}
	pPainter->setPen(penOrig);


	// angle arcs
	const QLineF* pLines1[] = { &lineKi, &lineRot };
	const QLineF* pLines2[] = { &lineKf, &lineKi };
	const QPointF* pPoints[] = { &ptSample, &ptSample };
	const t_real dAngles[] = { m_dTwoTheta, -m_dTheta };
	const t_real dAngleOffs[] = { 0., 180. };

	QPen* arcPens[] = {&penBlue, &penRed};
	const std::wstring strDEG = tl::get_spec_char_utf16("deg");

	for(std::size_t i=0; i<sizeof(pPoints)/sizeof(*pPoints); ++i)
	{
		t_real dArcSize = (pLines1[i]->length() + pLines2[i]->length()) / 2. / 3.;
		t_real dBeginArcAngle = pLines1[i]->angle() + dAngleOffs[i];
		t_real dArcAngle = tl::r2d(dAngles[i]);

		std::wostringstream ostrAngle;
		ostrAngle.precision(g_iPrecGfx);
		if(!tl::is_nan_or_inf<t_real>(dArcAngle))
		{
			ostrAngle << std::fabs(dArcAngle) << strDEG;
		}
		else
		{
			dArcAngle = t_real(180);
			ostrAngle << "invalid";
		}

		pPainter->setPen(*arcPens[i]);
		pPainter->drawArc(QRectF(pPoints[i]->x()-dArcSize/2., pPoints[i]->y()-dArcSize/2.,
			dArcSize, dArcSize), dBeginArcAngle*16., dArcAngle*16.);

		bool bFlip = dAngleOffs[i] > 90.;
		t_real dTotalAngle = -dBeginArcAngle-dArcAngle*0.5 + 180.;
		if(bFlip) dTotalAngle += 180.;
		t_real dTransScale = bFlip ? -40. : 80.;
		dTransScale *= m_dZoom;
		pPainter->save();
			pPainter->translate(*pPoints[i]);
			pPainter->rotate(dTotalAngle);
			pPainter->translate(-dTransScale, +4.*0.1*0.5*g_dFontSize);
			if(flip_text(dTotalAngle))
			{
				if(bFlip)
					pPainter->translate(-dTransScale, -8.*0.1*0.5*g_dFontSize);
				else
					pPainter->translate(dTransScale*0.5, -8.*0.1*0.5*g_dFontSize);
				pPainter->rotate(180.);
			}
			pPainter->drawText(QPointF(0.,0.), QString::fromWCharArray(ostrAngle.str().c_str()));
		pPainter->restore();
	}

	pPainter->setPen(penOrig);


	// arrow heads
	const QLineF* pLines_arrow[] = {&lineKi, &lineKf, plineQ.get()};
	const QPointF* pPoints_arrow[] = {&ptSample, &ptDet, pptQ.get()};
	QColor colArrowHead[] = {penBlack.color(), penBlack.color(), penRed.color()};
	for(std::size_t i=0; i<sizeof(pLines_arrow)/sizeof(*pLines_arrow); ++i)
	{
		if(!pLines_arrow[i] || !pPoints_arrow[i])
			continue;

		t_real dAng = tl::d2r<t_real>(pLines_arrow[i]->angle() - 90.);
		t_real dC = std::cos(dAng);
		t_real dS = std::sin(dAng);

		t_real dTriagX = 5., dTriagY = 10.;
		QPointF ptTriag1 = *pPoints_arrow[i] + QPointF(dTriagX*dC + dTriagY*dS, -dTriagX*dS + dTriagY*dC);
		QPointF ptTriag2 = *pPoints_arrow[i] + QPointF(-dTriagX*dC + dTriagY*dS, dTriagX*dS + dTriagY*dC);

		QPainterPath triag;
		triag.moveTo(*pPoints_arrow[i]);
		triag.lineTo(ptTriag1);
		triag.lineTo(ptTriag2);


		QPen penArrow(colArrowHead[i]);
		penArrow.setWidthF(g_dFontSize*0.1);

		pPainter->setPen(penArrow);
		pPainter->fillPath(triag, colArrowHead[i]);
	}

	pPainter->setPen(penOrig);
}


void TofLayout::SetSampleTwoTheta(t_real dAngle)
{
	m_dTwoTheta = dAngle;
	if(!m_scene.GetSampleSense())
		m_dTwoTheta = -m_dTwoTheta;

	t_vec vecSrc = qpoint_to_vec(mapFromItem(m_pSrc, 0, 0));
	t_vec vecSample = qpoint_to_vec(mapFromItem(m_pSample, 0, 0));
	t_vec vecDet = qpoint_to_vec(mapFromItem(m_pDet, 0, 0));

	t_vec vecKi = vecSample - vecSrc;
	vecKi /= ublas::norm_2(vecKi);
	t_real dLenKf = ublas::norm_2(vecDet-vecSample);

	t_vec vecKf = ublas::prod(tl::rotation_matrix_2d(-dAngle), vecKi);
	vecKf /= ublas::norm_2(vecKf);
	vecKf *= dLenKf;

	m_bUpdate = m_bReady = false;
	m_pDet->setPos(vec_to_qpoint(vecSample + vecKf));
	m_bReady = true;

	// don't call update twice
	nodeMoved(m_pSample);
	m_bUpdate = true;
	nodeMoved(m_pDet);
}


void TofLayout::SetSampleTheta(t_real dAngle)
{
	m_dTheta = dAngle;
	nodeMoved();
}


void TofLayout::SetAngleKiQ(t_real dAngle)
{
	m_dAngleKiQ = dAngle;
	nodeMoved();
}


void TofLayout::SetRealQVisible(bool bVisible)
{
	m_bRealQVisible = bVisible;
	this->update();
}


void TofLayout::SetZoom(t_real dZoom)
{
	m_dZoom = dZoom;
	m_scene.update();
}


std::vector<TofLayoutNode*> TofLayout::GetNodes() const
{
	return std::vector<TofLayoutNode*> { m_pSrc, m_pSample, m_pDet };
}


std::vector<std::string> TofLayout::GetNodeNames() const
{
	return std::vector<std::string> { "source", "sample", "detector" };
}

// --------------------------------------------------------------------------------


TofLayoutScene::TofLayoutScene(QObject *pParent) : QGraphicsScene(pParent)
{
	m_pTof = new TofLayout(*this);
	this->addItem(m_pTof);
}


TofLayoutScene::~TofLayoutScene()
{
	if(m_pTof)
		delete m_pTof;
}


void TofLayoutScene::emitAllParams()
{
	if(!m_pTof || !m_pTof->IsReady() || m_bDontEmitChange)
		return;

	RealParams parms;
	parms.dSampleTT = m_pTof->GetSampleTwoTheta();
	parms.dSampleT = m_pTof->GetSampleTheta();

	if(!GetSampleSense())
		parms.dSampleTT = -parms.dSampleTT;

	emit paramsChanged(parms);
}


void TofLayoutScene::triangleChanged(const TriangleOptions& opts)
{
	if(!m_pTof || !m_pTof->IsReady())
		return;

	m_bDontEmitChange = true;
	m_pTof->AllowChanges(false);

	if(opts.bChangedTheta)
		m_pTof->SetSampleTheta(opts.dTheta);
	if(opts.bChangedTwoTheta)
		m_pTof->SetSampleTwoTheta(opts.dTwoTheta);

	m_pTof->AllowChanges(true);
	m_bDontEmitChange = false;
}


void TofLayoutScene::recipParamsChanged(const RecipParams& params)
{
	m_pTof->SetAngleKiQ(params.dKiQ);
}


void TofLayoutScene::emitUpdate(const TriangleOptions& opts)
{
	if(!m_pTof || !m_pTof->IsReady() || m_bDontEmitChange)
		return;
	emit tasChanged(opts);
}


void TofLayoutScene::scaleChanged(t_real dTotalScale)
{
	if(!m_pTof)
		return;
	m_pTof->SetZoom(dTotalScale);
}


void TofLayoutScene::mousePressEvent(QGraphicsSceneMouseEvent *pEvt)
{
	QGraphicsScene::mousePressEvent(pEvt);
	emit nodeEvent(1);
}


void TofLayoutScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt)
{
	QGraphicsScene::mouseReleaseEvent(pEvt);
	emit nodeEvent(0);
}


// --------------------------------------------------------------------------------


TofLayoutView::TofLayoutView(QWidget* pParent) : QGraphicsView(pParent)
{
	setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
		QPainter::SmoothPixmapTransform | QPainter::HighQualityAntialiasing);
	setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);
	setDragMode(QGraphicsView::ScrollHandDrag);
	setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
}


TofLayoutView::~TofLayoutView()
{}


void TofLayoutView::DoZoom(t_real_glob dDelta)
{
	const t_real dScale = std::pow(2., dDelta);
	this->scale(dScale, dScale);

	m_dTotalScale *= dScale;
	emit scaleChanged(m_dTotalScale);
}


void TofLayoutView::keyPressEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Plus)
		DoZoom(0.02);
	else if(pEvt->key() == Qt::Key_Minus)
		DoZoom(-0.02);

	QGraphicsView::keyPressEvent(pEvt);
}


void TofLayoutView::keyReleaseEvent(QKeyEvent *pEvt)
{
	QGraphicsView::keyReleaseEvent(pEvt);
}


void TofLayoutView::wheelEvent(QWheelEvent *pEvt)
{
	const t_real dDelta = pEvt->angleDelta().y()/8. / 150.;
	DoZoom(dDelta);
}


#include "moc_tof_layout.cpp"
