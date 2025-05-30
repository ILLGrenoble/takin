/**
 * TAS layout
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2014
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

#include "tas_layout.h"
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


TasLayoutNode::TasLayoutNode(TasLayout* pSupItem) : m_pParentItem(pSupItem)
{
	setFlag(QGraphicsItem::ItemSendsGeometryChanges);
	setFlag(QGraphicsItem::ItemIgnoresTransformations);
	setCursor(Qt::CrossCursor);
}


QRectF TasLayoutNode::boundingRect() const
{
	return QRectF(-5.*0.1*g_dFontSize, -5.*0.1*g_dFontSize,
		10.*0.1*g_dFontSize, 10.*0.1*g_dFontSize);
}


void TasLayoutNode::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	pPainter->drawEllipse(QRectF(-2.*0.1*g_dFontSize, -2.*0.1*g_dFontSize,
		4.*0.1*g_dFontSize, 4.*0.1*g_dFontSize));
}


QVariant TasLayoutNode::itemChange(GraphicsItemChange change, const QVariant &val)
{
	QVariant var = QGraphicsItem::itemChange(change, val);

	if(change == QGraphicsItem::ItemPositionHasChanged)
		m_pParentItem->nodeMoved(this);

	return var;
}

// --------------------------------------------------------------------------------

TasLayout::TasLayout(TasLayoutScene& scene) : m_scene(scene),
	m_pSrc(new TasLayoutNode(this)),
	m_pMono(new TasLayoutNode(this)),
	m_pSample(new TasLayoutNode(this)),
	m_pAna(new TasLayoutNode(this)),
	m_pDet(new TasLayoutNode(this))
{
	setFlag(QGraphicsItem::ItemIgnoresTransformations);

	m_pSrc->setToolTip("Source");
	m_pMono->setToolTip("Monochromator");
	m_pSample->setToolTip("Sample");
	m_pAna->setToolTip("Analyser");
	m_pDet->setToolTip("Detector");

	m_pSrc->setPos(200., m_dLenMonoSample*m_dScaleFactor);
	m_pMono->setPos(0., m_dLenMonoSample*m_dScaleFactor);
	m_pSample->setPos(0., 0.);
	m_pAna->setPos(-m_dLenSampleAna*m_dScaleFactor, 0.);
	m_pDet->setPos(-100., -m_dLenAnaDet*m_dScaleFactor);

	AllowMouseMove(true);

	scene.addItem(m_pSrc.get());
	scene.addItem(m_pMono.get());
	scene.addItem(m_pSample.get());
	scene.addItem(m_pAna.get());
	scene.addItem(m_pDet.get());

	setAcceptedMouseButtons(Qt::NoButton);
	m_bUpdate = m_bReady = true;
}


TasLayout::~TasLayout()
{
	m_bUpdate = m_bReady = false;
}


void TasLayout::AllowMouseMove(bool bAllow)
{
	m_pMono->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pSample->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pAna->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
	m_pDet->setFlag(QGraphicsItem::ItemIsMovable, bAllow);
}


void TasLayout::nodeMoved(const TasLayoutNode *pNode)
{
	if(!m_bReady)
		return;

	// prevents recursive calling of update
	static bool bAllowUpdate = true;
	if(!bAllowUpdate)
		return;

	const t_vec vecSrc = qpoint_to_vec(mapFromItem(m_pSrc.get(), 0, 0));
	const t_vec vecMono = qpoint_to_vec(mapFromItem(m_pMono.get(), 0, 0));
	const t_vec vecSample = qpoint_to_vec(mapFromItem(m_pSample.get(), 0, 0));
	const t_vec vecAna = qpoint_to_vec(mapFromItem(m_pAna.get(), 0, 0));
	const t_vec vecDet = qpoint_to_vec(mapFromItem(m_pDet.get(), 0, 0));

	bAllowUpdate = false;
	if(pNode == m_pSample.get())
	{
		t_real dTwoTheta = m_dTwoTheta;
		t_real dAnaTwoTheta = m_dAnaTwoTheta;

		t_vec vecSrcMono = vecMono - vecSrc;
		vecSrcMono /= ublas::norm_2(vecSrcMono);

		t_vec vecMonoSample = vecSample - vecMono;
		if(m_bAllowChanges)
			m_dLenMonoSample = ublas::norm_2(vecMonoSample)/m_dScaleFactor;
		vecMonoSample /= ublas::norm_2(vecMonoSample);

		if(m_bAllowChanges)
		{
			m_dMonoTwoTheta = -(tl::vec_angle(vecMonoSample) - tl::vec_angle(vecSrcMono));
			if(m_dMonoTwoTheta < -tl::get_pi<t_real>()) m_dMonoTwoTheta += 2.*tl::get_pi<t_real>();
			if(m_dMonoTwoTheta > tl::get_pi<t_real>()) m_dMonoTwoTheta -= 2.*tl::get_pi<t_real>();
		}


		t_vec vecSampleAna =
			ublas::prod(tl::rotation_matrix_2d(-dTwoTheta), vecMonoSample);
		vecSampleAna /= ublas::norm_2(vecSampleAna);
		vecSampleAna *= m_dLenSampleAna*m_dScaleFactor;

		t_vec vecAnaNew = vecSample + vecSampleAna;
		m_pAna->setPos(vec_to_qpoint(vecAnaNew));


		vecSampleAna /= ublas::norm_2(vecSampleAna);

		t_vec vecAnaDet =
			ublas::prod(tl::rotation_matrix_2d(-dAnaTwoTheta), vecSampleAna);
		vecAnaDet /= ublas::norm_2(vecAnaDet);
		vecAnaDet *= m_dLenAnaDet*m_dScaleFactor;

		m_pDet->setPos(vec_to_qpoint(vecAnaNew+vecAnaDet));


		TriangleOptions opts;
		opts.bChangedMonoTwoTheta = true;
		opts.dMonoTwoTheta = m_dMonoTwoTheta;
		m_scene.emitUpdate(opts);
	}
	else if(pNode == m_pMono.get())
	{
		t_vec vecSrcMono = vecMono-vecSrc;
		vecSrcMono /= ublas::norm_2(vecSrcMono);

		t_vec vecMonoSample =
			ublas::prod(tl::rotation_matrix_2d(-m_dMonoTwoTheta), vecSrcMono);
		vecMonoSample /= ublas::norm_2(vecMonoSample);
		vecMonoSample *= m_dLenMonoSample*m_dScaleFactor;

		t_vec vecSampleNew = vecMono + vecMonoSample;
		m_pSample->setPos(vec_to_qpoint(vecSampleNew));


		t_vec vecSampleAna =
			ublas::prod(tl::rotation_matrix_2d(-m_dTwoTheta), vecMonoSample);
		vecSampleAna /= ublas::norm_2(vecSampleAna);
		vecSampleAna *= m_dLenSampleAna*m_dScaleFactor;

		t_vec vecAnaNew = vecSampleNew + vecSampleAna;
		m_pAna->setPos(vec_to_qpoint(vecAnaNew));


		t_vec vecAnaDet =
			ublas::prod(tl::rotation_matrix_2d(-m_dAnaTwoTheta), vecSampleAna);
		vecAnaDet /= ublas::norm_2(vecAnaDet);
		vecAnaDet *= m_dLenAnaDet*m_dScaleFactor;

		m_pDet->setPos(vec_to_qpoint(vecAnaNew + vecAnaDet));
	}
	else if(pNode == m_pDet.get())
	{
		t_vec vecSampleAna = vecAna - vecSample;
		vecSampleAna /= ublas::norm_2(vecSampleAna);

		t_vec vecAnaDet = vecDet - vecAna;
		if(m_bAllowChanges)
			m_dLenAnaDet = ublas::norm_2(vecAnaDet)/m_dScaleFactor;
		vecAnaDet /= ublas::norm_2(vecAnaDet);

		if(m_bAllowChanges)
		{
			m_dAnaTwoTheta = -(tl::vec_angle(vecAnaDet) - tl::vec_angle(vecSampleAna));
			if(m_dAnaTwoTheta < -tl::get_pi<t_real>()) m_dAnaTwoTheta += 2.*tl::get_pi<t_real>();
			if(m_dAnaTwoTheta > tl::get_pi<t_real>()) m_dAnaTwoTheta -= 2.*tl::get_pi<t_real>();
		}

		TriangleOptions opts;
		opts.bChangedAnaTwoTheta = true;
		opts.dAnaTwoTheta = m_dAnaTwoTheta;
		m_scene.emitUpdate(opts);
	}

	if(pNode == m_pAna.get())
	{
		t_vec vecSampleAna = vecAna-vecSample;
		if(pNode == m_pAna.get() && m_bAllowChanges)
			m_dLenSampleAna = ublas::norm_2(vecSampleAna)/m_dScaleFactor;
		vecSampleAna /= ublas::norm_2(vecSampleAna);

		t_vec vecAnaDet =
			ublas::prod(tl::rotation_matrix_2d(-m_dAnaTwoTheta), vecSampleAna);
		vecAnaDet /= ublas::norm_2(vecAnaDet);
		vecAnaDet *= m_dLenAnaDet*m_dScaleFactor;

		m_pDet->setPos(vec_to_qpoint(vecAna+vecAnaDet));

		t_vec vecMonoSample = vecSample - vecMono;
		vecMonoSample /= ublas::norm_2(vecMonoSample);

		if(m_bAllowChanges)
		{
			m_dTwoTheta = -(tl::vec_angle(vecSampleAna) - tl::vec_angle(vecMonoSample));
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


QRectF TasLayout::boundingRect() const
{
	return QRectF(-100.*m_dZoom*g_dFontSize, -100.*m_dZoom*g_dFontSize,
		200.*m_dZoom*g_dFontSize, 200.*m_dZoom*g_dFontSize);
}


/**
 * gets the centre of the tas drawing
 */
QPointF TasLayout::GetGfxMid() const
{
	t_real num_nodes = 0.;
	QPointF mid(0., 0.);

	for(const TasLayoutNode* node : GetNodes())
	{
		mid += mapFromItem(node, 0, 0);
		num_nodes += 1.;
	}

	return mid / num_nodes;
}


/**
 * draws the TAS schematic
 */
void TasLayout::paint(QPainter *pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
	const bool bDisplayLengths = false;
	pPainter->setFont(g_fontGfx);

	QPointF ptSrc = mapFromItem(m_pSrc.get(), 0, 0) * m_dZoom;
	QPointF ptMono = mapFromItem(m_pMono.get(), 0, 0) * m_dZoom;
	QPointF ptSample = mapFromItem(m_pSample.get(), 0, 0) * m_dZoom;
	QPointF ptAna = mapFromItem(m_pAna.get(), 0, 0) * m_dZoom;
	QPointF ptDet = mapFromItem(m_pDet.get(), 0, 0) * m_dZoom;

	QLineF lineSrcMono(ptSrc, ptMono);
	QLineF lineKi(ptMono, ptSample);
	QLineF lineKf(ptSample, ptAna);
	QLineF lineAnaDet(ptAna, ptDet);

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

	pPainter->drawLine(lineSrcMono);
	pPainter->drawLine(lineKi);
	pPainter->drawLine(lineKf);
	pPainter->drawLine(lineAnaDet);


	// ------------------------------------------------------------------------
	// write lengths
	QPointF ptMidKi = ptMono + (ptSample-ptMono)/2.;
	QPointF ptMidKf = ptSample + (ptAna-ptSample)/2.;
	QPointF ptMidAnaDet = ptAna + (ptDet-ptAna)/2.;

	if(bDisplayLengths)
	{
		std::ostringstream ostrLenKi, ostrLenKf, ostrLenAnaDet;
		ostrLenKi.precision(g_iPrecGfx);
		ostrLenKf.precision(g_iPrecGfx);
		ostrLenAnaDet.precision(g_iPrecGfx);

		ostrLenKi << m_dLenMonoSample << " cm";
		ostrLenKf << m_dLenSampleAna << " cm";
		ostrLenAnaDet << m_dLenAnaDet << " cm";

		pPainter->drawText(ptMidKi, ostrLenKi.str().c_str());
		pPainter->drawText(ptMidKf, ostrLenKf.str().c_str());
		pPainter->drawText(ptMidAnaDet, ostrLenAnaDet.str().c_str());
	}
	// ------------------------------------------------------------------------



	// ------------------------------------------------------------------------
	t_vec vecSrc = qpoint_to_vec(ptSrc);
	t_vec vecMono = qpoint_to_vec(ptMono);
	t_vec vecSample = qpoint_to_vec(ptSample);
	t_vec vecAna = qpoint_to_vec(ptAna);
	t_vec vecDet = qpoint_to_vec(ptDet);

	t_vec vecSrcMono = vecMono-vecSrc;
	t_vec vecMonoSample = vecSample-vecMono;
	t_vec vecSampleAna = vecAna-vecSample;
	t_vec vecAnaDet = vecDet-vecAna;

	t_real dThetas[] = {-m_dMonoTwoTheta/t_real(2.), -m_dAnaTwoTheta/t_real(2.), -m_dTheta};
	std::vector<const t_vec*> vecPos = {&vecMono, &vecAna, &vecSample};
	std::vector<const t_vec*> vecDirs = {&vecSrcMono, &vecSampleAna, &vecMonoSample};
	QColor colThs[] = {penGray.color(), penGray.color(), penRed.color()};
	const char* pcComp[] = {"M", "A", "S"};
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// mono/ana/sample theta rotation
	QLineF lineRot[3];
	QPointF ptThP[3];

	for(std::size_t iTh=0; iTh<sizeof(dThetas)/sizeof(*dThetas); ++iTh)
	{
		t_vec vecRotDir =
			ublas::prod(tl::rotation_matrix_2d(dThetas[iTh]), *vecDirs[iTh]);
		vecRotDir /= ublas::norm_2(vecRotDir);
		vecRotDir *= m_dLenSample*m_dScaleFactor;

		QPointF ptThM = vec_to_qpoint(*vecPos[iTh] - vecRotDir*m_dZoom);
		ptThP[iTh] = vec_to_qpoint(*vecPos[iTh] + vecRotDir*m_dZoom);
		lineRot[iTh] = QLineF(ptThM, ptThP[iTh]);

		QPen pen(colThs[iTh]);
		pen.setWidthF(1.5*g_dFontSize*0.1);
		pPainter->setPen(pen);
		pPainter->drawLine(lineRot[iTh]);


		// component names
		pPainter->setPen(penOrig);
		pPainter->save();
			pPainter->translate(vec_to_qpoint(*vecPos[iTh]));
			t_real dCompAngle = 180. + tl::r2d(tl::vec_angle(vecRotDir));
			pPainter->rotate(dCompAngle);
			pPainter->translate(-4., 16.*0.1*g_dFontSize);
			if(flip_text(dCompAngle))
			{
				pPainter->translate(4., -8.*0.1*g_dFontSize);
				pPainter->rotate(180.);
			}
			pPainter->drawText(QPointF(0., 0.), pcComp[iTh]);
		pPainter->restore();
	}

	pPainter->drawText(ptMono - vec_to_qpoint(vecSrcMono*1.1), "R");
	pPainter->drawText(ptAna + vec_to_qpoint(vecAnaDet*1.1), "D");
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// dashed extended lines
	QPen penDash(Qt::DashLine);
	penDash.setColor(qApp->palette().color(QPalette::WindowText));
	penDash.setWidthF(g_dFontSize*0.1);

	pPainter->setPen(penDash);
	QLineF lineSrcMono_ext(ptMono, ptMono + (ptMono-ptSrc)/2.);
	QLineF lineki_ext(ptSample, ptSample + (ptSample-ptMono)/2.);
	QLineF linekf_ext(ptAna, ptAna + (ptAna-ptSample)/2.);

	pPainter->drawLine(lineSrcMono_ext);
	pPainter->drawLine(lineki_ext);
	pPainter->drawLine(linekf_ext);
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// dark angles
	bool bBeamObstructed = false;

	struct DarkAngleArc
	{
		t_real dCentre[2];
		t_real dAngleStart, dAngleRange;

		QLineF *pLineIn = nullptr;
		QLineF *pLineOut = nullptr;

		bool intersects() const
		{
			for(const QLineF* pLine : {pLineIn, pLineOut})
			{
				if(!pLine) continue;

				t_real dLineAngle = pLine->angle();
				if(pLine == pLineIn)
				{
					dLineAngle += 180.;
					dLineAngle = std::fmod(dLineAngle, 360.);
				}

				bool bIntersect = tl::is_in_angular_range(
					tl::d2r(dAngleStart), tl::d2r(dAngleRange), tl::d2r(dLineAngle));
				if(bIntersect)
				{
					return true;
				}
			}
			return false;
		}
	};


	if(m_pvecDarkAngles)
	{
		QPen pen(QColor(0x90,0x50,0x20));
		pen.setWidthF(3.*g_dFontSize*0.1);
		pPainter->setPen(pen);

		const t_real dArcSize = (lineKi.length() + lineKf.length()) / 2. / 2.;

		for(const DarkAngle<t_real>& angle : *m_pvecDarkAngles)
		{
			DarkAngleArc darkarc;

			// default case: around sample
			QPointF *pCentre = &ptSample;
			darkarc.pLineIn = &lineKi;
			darkarc.pLineOut = &lineKf;
			t_real dCrystalTheta = tl::r2d(m_dTheta);

			switch(angle.iCentreOn)
			{
				case 0:	// around mono
				{
					pCentre = &ptMono;
					darkarc.pLineIn = &lineSrcMono;
					darkarc.pLineOut = &lineKi;
					dCrystalTheta = tl::r2d(m_dMonoTwoTheta/t_real(2.));
					break;
				}
				case 2:	// around ana
				{
					pCentre = &ptAna;
					darkarc.pLineIn = &lineKf;
					darkarc.pLineOut = &lineAnaDet;
					dCrystalTheta = tl::r2d(m_dAnaTwoTheta/t_real(2.));
					break;
				}
			}

			t_real dAbsOffs = 0.;
			switch(angle.iRelativeTo)
			{
				case 0:	// relative to crystal angle theta
				{
					dAbsOffs = darkarc.pLineIn->angle() + dCrystalTheta;
					break;
				}
				case 1:	// relative to incoming axis
				{
					dAbsOffs = darkarc.pLineIn->angle();
					break;
				}
				case 2:	// relative to outgoing axis
				{
					dAbsOffs = darkarc.pLineOut->angle();
					break;
				}
			}

			darkarc.dAngleStart = angle.dAngleStart + angle.dAngleOffs + dAbsOffs;
			darkarc.dAngleRange = angle.dAngleEnd - angle.dAngleStart;

			pPainter->drawArc(QRectF(pCentre->x()-dArcSize/2., pCentre->y()-dArcSize/2.,
				dArcSize, dArcSize), darkarc.dAngleStart*16., darkarc.dAngleRange*16.);

			if(!bBeamObstructed && darkarc.intersects())
				bBeamObstructed = true;
		}
	}
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	std::unique_ptr<QLineF> plineQ(nullptr);
	std::unique_ptr<QPointF> pptQ(nullptr);
	// Q vector direction visible?
	if(this->m_bRealQVisible)
	{
		t_real dAngleKiQ = m_dAngleKiQ;
		t_mat matRotQ = tl::rotation_matrix_2d(dAngleKiQ);

		t_vec vecKi = vecSample - vecMono;
		t_vec vecQ = ublas::prod(matRotQ, vecKi);
		vecQ /= ublas::norm_2(vecQ);
		vecQ *= (m_dLenMonoSample + m_dLenSampleAna)/2.;	// some arbitrary length
		vecQ *= m_dScaleFactor * m_dZoom;

		pptQ.reset(new QPointF(vec_to_qpoint(vecSample + vecQ)));
		plineQ.reset(new QLineF(ptSample, *pptQ));

		pPainter->setPen(penRed);
		pPainter->drawLine(*plineQ);
		pPainter->save();
			pPainter->translate(ptSample);
			const t_real dQAngle = -plineQ->angle();
			pPainter->rotate(dQAngle);
			pPainter->translate(QPointF(plineQ->length()/2., 12.*0.1*g_dFontSize));
			if(flip_text(dQAngle))
				pPainter->rotate(180.);
			pPainter->drawText(QPointF(0,0), "Q");
		pPainter->restore();
	}
	pPainter->setPen(penOrig);
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// angle arcs
	const QLineF* pLines1[] = {&lineSrcMono, &lineKi, &lineKf, &lineRot[2]};
	const QLineF* pLines2[] = {&lineKi, &lineKf, &lineAnaDet, &lineKi};
	const QPointF* pPoints[] = {&ptMono, &ptSample, &ptAna, &ptSample};
	const t_real dAngles[] = {m_dMonoTwoTheta, m_dTwoTheta, m_dAnaTwoTheta, -m_dTheta};
	const t_real dAngleOffs[] = {0., 0., 0., 180.};

	QPen* arcPens[] = {&penBlue, &penBlue, &penBlue, &penRed};
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
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// arrow heads
	const QLineF* pLines_arrow[] = {&lineKi, &lineKf, plineQ.get(), &lineSrcMono, &lineAnaDet};
	const QPointF* pPoints_arrow[] = {&ptSample, &ptAna, pptQ.get(), &ptMono, &ptDet};
	QColor colArrowHead[] = {penBlack.color(), penBlack.color(), penRed.color(), penGray.color(), penGray.color()};
	for(std::size_t i=0; i<sizeof(pLines_arrow)/sizeof(*pLines_arrow); ++i)
	{
		if(!pLines_arrow[i] || !pPoints_arrow[i])
			continue;

		t_real dAng = tl::d2r(pLines_arrow[i]->angle() - 90.);
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



	// ------------------------------------------------------------------------
	if(bBeamObstructed)
	{
		pPainter->setPen(penRed);
		pPainter->drawText(QPointF(0.,0.), "Beam obstructed!");
	}
	// ------------------------------------------------------------------------


	pPainter->setPen(penOrig);
}


void TasLayout::SetSampleTwoTheta(t_real dAngle)
{
	m_dTwoTheta = dAngle;
	if(!m_scene.GetSampleSense())
		m_dTwoTheta = -m_dTwoTheta;

	t_vec vecMono = qpoint_to_vec(mapFromItem(m_pMono.get(), 0, 0));
	t_vec vecSample = qpoint_to_vec(mapFromItem(m_pSample.get(), 0, 0));
	t_vec vecAna = qpoint_to_vec(mapFromItem(m_pAna.get(), 0, 0));

	t_vec vecKi = vecSample - vecMono;
	vecKi /= ublas::norm_2(vecKi);
	t_real dLenKf = ublas::norm_2(vecAna-vecSample);

	t_vec vecKf = ublas::prod(tl::rotation_matrix_2d(-dAngle), vecKi);
	vecKf /= ublas::norm_2(vecKf);
	vecKf *= dLenKf;


	m_bReady = m_bUpdate = false;
	m_pAna->setPos(vec_to_qpoint(vecSample + vecKf));
	m_bReady = true;

	// don't call update twice
	nodeMoved(m_pSample.get());
	m_bUpdate = true;
	nodeMoved(m_pAna.get());
}


void TasLayout::SetSampleTheta(t_real dAngle)
{
	//tl::log_info("sample theta: ", dAngle/M_PI*180.);
	m_dTheta = dAngle;
	nodeMoved();
}


void TasLayout::SetMonoTwoTheta(t_real dAngle)
{
	m_dMonoTwoTheta = dAngle;
	nodeMoved(m_pMono.get());
}


void TasLayout::SetAnaTwoTheta(t_real dAngle)
{
	m_dAnaTwoTheta = dAngle;
	nodeMoved(m_pAna.get());
}


void TasLayout::SetAngleKiQ(t_real dAngle)
{
	m_dAngleKiQ = dAngle;
	nodeMoved();
}


void TasLayout::SetRealQVisible(bool bVisible)
{
	m_bRealQVisible = bVisible;
	this->update();
}


void TasLayout::SetZoom(t_real dZoom)
{
	m_dZoom = dZoom;
	m_scene.update();
}


std::vector<TasLayoutNode*> TasLayout::GetNodes() const
{
	return std::vector<TasLayoutNode*>
			{ m_pSrc.get(), m_pMono.get(), m_pSample.get(),
			m_pAna.get(), m_pDet.get() };
}


std::vector<std::string> TasLayout::GetNodeNames() const
{
	return std::vector<std::string>
		{ "source", "monochromator", "sample", "analyser", "detector" };
}


void TasLayout::SetDarkAngles(const std::vector<DarkAngle<t_real_glob>> *pvecDarkAngles)
{
	m_pvecDarkAngles = pvecDarkAngles;
	update();
}

// --------------------------------------------------------------------------------


TasLayoutScene::TasLayoutScene(QObject *pParent)
	: QGraphicsScene(pParent), m_pTas(new TasLayout(*this))
{
	this->addItem(m_pTas.get());
}


TasLayoutScene::~TasLayoutScene()
{}


void TasLayoutScene::emitAllParams()
{
	if(!m_pTas || !m_pTas->IsReady() || m_bDontEmitChange)
		return;

	RealParams parms;
	parms.dAnaTT = m_pTas->GetAnaTwoTheta();
	parms.dAnaT = m_pTas->GetAnaTheta();
	parms.dMonoTT = m_pTas->GetMonoTwoTheta();
	parms.dMonoT = m_pTas->GetMonoTheta();
	parms.dSampleTT = m_pTas->GetSampleTwoTheta();
	parms.dSampleT = m_pTas->GetSampleTheta();

	parms.dLenMonoSample = m_pTas->GetLenMonoSample();
	parms.dLenSampleAna = m_pTas->GetLenSampleAna();
	parms.dLenAnaDet = m_pTas->GetLenAnaDet();

	if(!GetSampleSense())
		parms.dSampleTT = -parms.dSampleTT;

	//log_info(parms.dSampleT/M_PI*180.);
	//log_debug("tas: emitAllParams");
	emit paramsChanged(parms);
}


void TasLayoutScene::triangleChanged(const TriangleOptions& opts)
{
	if(!m_pTas || !m_pTas->IsReady())
		return;

	m_bDontEmitChange = true;
	m_pTas->AllowChanges(false);

	if(opts.bChangedMonoTwoTheta)
		m_pTas->SetMonoTwoTheta(opts.dMonoTwoTheta);
	if(opts.bChangedAnaTwoTheta)
		m_pTas->SetAnaTwoTheta(opts.dAnaTwoTheta);
	if(opts.bChangedTheta)
		m_pTas->SetSampleTheta(opts.dTheta);
	if(opts.bChangedTwoTheta)
		m_pTas->SetSampleTwoTheta(opts.dTwoTheta);

	m_pTas->AllowChanges(true);
	m_bDontEmitChange = false;
}


void TasLayoutScene::recipParamsChanged(const RecipParams& params)
{
	m_pTas->SetAngleKiQ(params.dKiQ);
}


void TasLayoutScene::emitUpdate(const TriangleOptions& opts)
{
	if(!m_pTas || !m_pTas->IsReady() || m_bDontEmitChange)
		return;
	emit tasChanged(opts);
}


void TasLayoutScene::scaleChanged(t_real dTotalScale)
{
	if(!m_pTas)
		return;
	m_pTas->SetZoom(dTotalScale);
}


void TasLayoutScene::mousePressEvent(QGraphicsSceneMouseEvent *pEvt)
{
	QGraphicsScene::mousePressEvent(pEvt);
	emit nodeEvent(1);
}


void TasLayoutScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt)
{
	QGraphicsScene::mouseReleaseEvent(pEvt);
	emit nodeEvent(0);
}


// --------------------------------------------------------------------------------


TasLayoutView::TasLayoutView(QWidget* pParent) : QGraphicsView(pParent)
{
	setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
		QPainter::SmoothPixmapTransform | QPainter::HighQualityAntialiasing);
	setViewportUpdateMode(QGraphicsView::BoundingRectViewportUpdate);
	setDragMode(QGraphicsView::ScrollHandDrag);
	setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
}


TasLayoutView::~TasLayoutView()
{}


void TasLayoutView::DoZoom(t_real_glob dDelta)
{
	const t_real dScale = std::pow(2., dDelta);
	this->scale(dScale, dScale);

	m_dTotalScale *= dScale;
	emit scaleChanged(m_dTotalScale);
}


void TasLayoutView::keyPressEvent(QKeyEvent *pEvt)
{
	if(pEvt->key() == Qt::Key_Plus)
		DoZoom(0.02);
	else if(pEvt->key() == Qt::Key_Minus)
		DoZoom(-0.02);

	QGraphicsView::keyPressEvent(pEvt);
}


void TasLayoutView::keyReleaseEvent(QKeyEvent *pEvt)
{
	QGraphicsView::keyReleaseEvent(pEvt);
}


void TasLayoutView::wheelEvent(QWheelEvent *pEvt)
{
	const t_real dDelta = pEvt->angleDelta().y()/8. / 150.;
	DoZoom(dDelta);
}


#include "moc_tas_layout.cpp"
