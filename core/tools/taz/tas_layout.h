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

#ifndef __TAS_LAYOUT_H__
#define __TAS_LAYOUT_H__

#include "tlibs/helper/flags.h"
#include "libs/globals.h"
#include "libs/globals_qt.h"

#include <cmath>
#include <memory>
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
#include "dialogs/DarkAnglesDlg.h"	// for DarkAngle struct


class TasLayout;
class TasLayoutNode : public QGraphicsItem
{
	protected:
		TasLayout *m_pParentItem = nullptr;

	protected:
		virtual QRectF boundingRect() const override;
		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

		virtual QVariant itemChange(GraphicsItemChange change, const QVariant &val) override;

	public:
		TasLayoutNode(TasLayout* pSupItem);
};


class TasLayoutScene;
class TasLayout : public QGraphicsItem
{
	protected:
		bool m_bReady = false, m_bUpdate = false;
		TasLayoutScene& m_scene;

		std::unique_ptr<TasLayoutNode> m_pSrc,
			m_pMono, m_pSample, m_pAna, m_pDet;

		t_real_glob m_dMonoTwoTheta = 3.1415/2.;
		t_real_glob m_dAnaTwoTheta = 3.1415/2.;
		t_real_glob m_dTwoTheta = 3.1415/2.;
		t_real_glob m_dTheta = 3.1415/4.;
		t_real_glob m_dAngleKiQ = 3.1415/4.;

		t_real_glob m_dLenMonoSample = 150.;
		t_real_glob m_dLenSampleAna = 100.;
		t_real_glob m_dLenAnaDet = 50.;
		t_real_glob m_dLenSample = 25.;

		t_real_glob m_dScaleFactor = 1.4; // pixels per cm for zoom == 1
		t_real_glob m_dZoom = 1.;

		bool m_bRealQVisible = true;
		bool m_bAllowChanges = true;

		const std::vector<DarkAngle<t_real_glob>> *m_pvecDarkAngles = nullptr;


	public:
		t_real_glob GetMonoTwoTheta() const { return m_dMonoTwoTheta; }
		t_real_glob GetMonoTheta() const { return m_dMonoTwoTheta/2.; }
		t_real_glob GetAnaTwoTheta() const { return m_dAnaTwoTheta; }
		t_real_glob GetAnaTheta() const { return m_dAnaTwoTheta/2.; }
		t_real_glob GetSampleTwoTheta() const { return m_dTwoTheta; }
		t_real_glob GetSampleTheta() const { return m_dTheta; }
		t_real_glob GetLenMonoSample() const { return m_dLenMonoSample; }
		t_real_glob GetLenSampleAna() const { return m_dLenSampleAna; }
		t_real_glob GetLenAnaDet() const { return m_dLenAnaDet; }

	protected:
		virtual QRectF boundingRect() const override;


	public:
		TasLayout(TasLayoutScene& scene);
		virtual ~TasLayout();

		bool IsReady() const { return m_bReady; }

		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

		void SetSampleTwoTheta(t_real_glob dTT);
		void SetSampleTheta(t_real_glob dT);
		void SetMonoTwoTheta(t_real_glob dTT);
		void SetAnaTwoTheta(t_real_glob dTT);
		void SetAngleKiQ(t_real_glob dKiQ);

		void SetReady(bool bReady) { m_bReady = bReady; }
		void nodeMoved(const TasLayoutNode* pNode=0);

		void AllowChanges(bool bAllow) { m_bAllowChanges = bAllow; };
		void AllowMouseMove(bool bAllow);

		void SetZoom(t_real_glob dZoom);
		t_real_glob GetZoom() const { return m_dZoom; }

		std::vector<TasLayoutNode*> GetNodes() const;
		std::vector<std::string> GetNodeNames() const;
		QPointF GetGfxMid() const;

		t_real_glob GetScaleFactor() const { return m_dScaleFactor; }
		void SetScaleFactor(t_real_glob dScale) { m_dScaleFactor = dScale; }

		void SetRealQVisible(bool bVisible);
		bool GetRealQVisible() const { return m_bRealQVisible; }

		void SetDarkAngles(const std::vector<DarkAngle<t_real_glob>> *pvecDarkAngles);
};


class TasLayoutScene : public QGraphicsScene
{	Q_OBJECT
	protected:
		std::unique_ptr<TasLayout> m_pTas;
		bool m_bDontEmitChange = true;
		bool m_bSamplePosSense = true;

	protected:
		virtual void mousePressEvent(QGraphicsSceneMouseEvent *pEvt) override;
		virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt) override;


	public:
		TasLayoutScene(QObject *pParent = nullptr);
		virtual ~TasLayoutScene();

		bool GetSampleSense() const { return m_bSamplePosSense; }
		void SetSampleSense(bool bPos) { m_bSamplePosSense = bPos; }

		void SetEmitChanges(bool bEmit) { m_bDontEmitChange = !bEmit; }
		void emitUpdate(const TriangleOptions& opts);

		TasLayout* GetTasLayout() { return m_pTas.get(); }

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


class TasLayoutView : public QGraphicsView
{
	Q_OBJECT
	protected:
		t_real_glob m_dTotalScale = 1.;

		void DoZoom(t_real_glob delta);
		virtual void wheelEvent(QWheelEvent* pEvt) override;
		virtual void keyPressEvent(QKeyEvent *pEvt) override;
		virtual void keyReleaseEvent(QKeyEvent *pEvt) override;

	public:
		TasLayoutView(QWidget* pParent = nullptr);
		virtual ~TasLayoutView();

	signals:
		void scaleChanged(t_real_glob dTotalScale);
};


#endif
