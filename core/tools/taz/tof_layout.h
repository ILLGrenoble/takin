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

#ifndef __TOF_LAYOUT_H__
#define __TOF_LAYOUT_H__

#include "tlibs/helper/flags.h"
#include "libs/globals.h"
#include "libs/globals_qt.h"

#include <cmath>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsItem>
#include <QGraphicsSceneDragDropEvent>
#include <QGraphicsTextItem>
#include <QWheelEvent>
#include <QtWidgets>

#include "tasoptions.h"
#include "dialogs/RealParamDlg.h"	// for RealParams struct
#include "dialogs/RecipParamDlg.h"	// for RecipParams struct


class TofLayout;
class TofLayoutNode : public QGraphicsItem
{
	protected:
		TofLayout *m_pParentItem = nullptr;

	protected:
		virtual QRectF boundingRect() const override;
		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

		virtual QVariant itemChange(GraphicsItemChange change, const QVariant &val) override;

	public:
		TofLayoutNode(TofLayout* pSupItem);
};


class TofLayoutScene;
class TofLayout : public QGraphicsItem
{
	protected:
		bool m_bReady = false, m_bUpdate = false;
		TofLayoutScene& m_scene;

		TofLayoutNode *m_pSrc = nullptr;
		TofLayoutNode *m_pSample = nullptr;
		TofLayoutNode *m_pDet = nullptr;

		t_real_glob m_dTwoTheta = 3.1415/2.;
		t_real_glob m_dTheta = 3.1415/4.;
		t_real_glob m_dAngleKiQ = 3.1415/4.;

		t_real_glob m_dLenSrcSample = 150.;
		t_real_glob m_dLenSampleDet = 100.;
		t_real_glob m_dLenSample = 25.;

		t_real_glob m_dBeginDetArcAngle = -140.;
		t_real_glob m_dDetArcAngle = -m_dBeginDetArcAngle*2.;

		t_real_glob m_dScaleFactor = 1.4; 	// pixels per cm for zoom == 1
		t_real_glob m_dZoom = 1.;

		bool m_bRealQVisible = true;
		bool m_bAllowChanges = true;

	public:
		t_real_glob GetSampleTwoTheta() const { return m_dTwoTheta; }
		t_real_glob GetSampleTheta() const { return m_dTheta; }
		t_real_glob GetLenSrcSample() const { return m_dLenSrcSample; }
		t_real_glob GetLenSampleDet() const { return m_dLenSampleDet; }

	protected:
		QRectF boundingRect() const override;

	public:
		TofLayout(TofLayoutScene& scene);
		virtual ~TofLayout();

		bool IsReady() const { return m_bReady; }

		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

		void SetSampleTwoTheta(t_real_glob dTT);
		void SetSampleTheta(t_real_glob dT);
		void SetAngleKiQ(t_real_glob dKiQ);

		void SetReady(bool bReady) { m_bReady = bReady; }
		void nodeMoved(const TofLayoutNode* pNode = nullptr);

		void AllowChanges(bool bAllow) { m_bAllowChanges = bAllow; };
		void AllowMouseMove(bool bAllow);

		void SetZoom(t_real_glob dZoom);
		t_real_glob GetZoom() const { return m_dZoom; }

	public:
		std::vector<TofLayoutNode*> GetNodes() const;
		std::vector<std::string> GetNodeNames() const;
		QPointF GetGfxMid() const;

		t_real_glob GetScaleFactor() const { return m_dScaleFactor; }
		void SetScaleFactor(t_real_glob dScale) { m_dScaleFactor = dScale; }

		void SetRealQVisible(bool bVisible);
		bool GetRealQVisible() const { return m_bRealQVisible; }
};


class TofLayoutScene : public QGraphicsScene
{	Q_OBJECT
	protected:
		TofLayout *m_pTof = nullptr;
		bool m_bDontEmitChange = true;
		bool m_bSamplePosSense = true;

	protected:
		virtual void mousePressEvent(QGraphicsSceneMouseEvent *pEvt) override;
		virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt) override;


	public:
		TofLayoutScene(QObject *pParent=nullptr);
		virtual ~TofLayoutScene();

		bool GetSampleSense() const { return m_bSamplePosSense; }
		void SetSampleSense(bool bPos) { m_bSamplePosSense = bPos; }

		void SetEmitChanges(bool bEmit) { m_bDontEmitChange = !bEmit; }
		void emitUpdate(const TriangleOptions& opts);

		TofLayout* GetTofLayout() { return m_pTof; }

		void emitAllParams();

	public slots:
		void triangleChanged(const TriangleOptions& opts);
		void scaleChanged(t_real_glob dTotalScale);

		void recipParamsChanged(const RecipParams& params);

	signals:
		// relevant parameters for triangle view
		void tasChanged(const TriangleOptions& opts);
		// all parameters
		void paramsChanged(const RealParams& parms);

		void nodeEvent(bool bStarted);
};


class TofLayoutView : public QGraphicsView
{
	Q_OBJECT
	protected:
		t_real_glob m_dTotalScale = 1.;

		void DoZoom(t_real_glob delta);
		virtual void wheelEvent(QWheelEvent* pEvt) override;
		virtual void keyPressEvent(QKeyEvent *pEvt) override;
		virtual void keyReleaseEvent(QKeyEvent *pEvt) override;

	public:
		TofLayoutView(QWidget* pParent = nullptr);
		virtual ~TofLayoutView();

	signals:
		void scaleChanged(t_real_glob dTotalScale);
};


#endif
