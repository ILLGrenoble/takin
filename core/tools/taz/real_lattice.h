/**
 * Real crystal lattice
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014 - 2016
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

#ifndef __TAZ_REAL_LATTICE_H__
#define __TAZ_REAL_LATTICE_H__

#include "tlibs/math/linalg.h"
#include "tlibs/math/kd.h"
#include "tlibs/phys/lattice.h"
#include "tlibs/phys/bz.h"
#include "tlibs/phys/neutrons.h"

#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "libs/spacegroups/spacegroup.h"
#include "libs/spacegroups/latticehelper.h"

#include "tasoptions.h"
#include "dialogs/AtomsDlg.h"

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsItem>
#include <QGraphicsSceneDragDropEvent>
#include <QGraphicsTextItem>
#include <QWheelEvent>
#include <QtWidgets>

namespace ublas = boost::numeric::ublas;


class RealLattice;

class LatticePoint : public QGraphicsItem
{
	protected:
		QColor m_color = Qt::red;
		virtual QRectF boundingRect() const override;
		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

	protected:
		QString m_strLabel;

	public:
		LatticePoint();

		void SetLabel(const QString& str) { m_strLabel = str; }
		void SetColor(const QColor& col) { m_color = col; }
};


class LatticeAtom : public QGraphicsItem
{
	friend class RealLattice;

	protected:
		QColor m_color = Qt::cyan;
		virtual QRectF boundingRect() const override;
		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

	protected:
		std::string m_strElem;
		ublas::vector<t_real_glob> m_vecPos;
		ublas::vector<t_real_glob> m_vecProj;
		t_real_glob m_dProjDist = 0.;

	public:
		LatticeAtom();
		void SetColor(const QColor& col) { m_color = col; }
};


class LatticeScene;
class RealLattice : public QGraphicsItem
{
	protected:
		bool m_bReady = false;
		LatticeScene &m_scene;

		t_real_glob m_dScaleFactor = 48.;	// pixels per A for zoom == 1.
		t_real_glob m_dZoom = 1.;
		t_real_glob m_dPlaneDistTolerance = 0.01;
		int m_iMaxPeaks = 7;

		tl::Lattice<t_real_glob> m_lattice;
		ublas::matrix<t_real_glob> m_matPlane, m_matPlane_inv;
		std::vector<LatticePoint*> m_vecPeaks;
		std::vector<LatticeAtom*> m_vecAtoms;

		tl::Kd<t_real_glob> m_kdLattice;

		bool m_bShowWS = true;
		tl::Brillouin2D<t_real_glob> m_ws;	// "Wigner-Seitz cell"
		tl::Brillouin3D<t_real_glob> m_ws3;	// "Wigner-Seitz cell"
		std::vector<ublas::vector<t_real_glob>> m_vecWS3VertsUnproj, m_vecWS3Verts;

	protected:
		virtual QRectF boundingRect() const override;

	public:
		RealLattice(LatticeScene& scene);
		virtual ~RealLattice();

		virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

		void SetReady(bool bReady) { m_bReady = bReady; }
		bool IsReady() const { return m_bReady; }

		const ublas::matrix<t_real_glob>& GetPlane() const { return m_matPlane; }

	public:
		bool HasPeaks() const { return m_vecPeaks.size()!=0 && m_lattice.IsInited(); }
		void ClearPeaks();
		void CalcPeaks(const xtl::LatticeCommon<t_real_glob>& latticecommon);

		void SetPlaneDistTolerance(t_real_glob dTol) { m_dPlaneDistTolerance = dTol; }
		void SetMaxPeaks(int iMax) { m_iMaxPeaks = iMax; }
		unsigned int GetMaxPeaks() const { return m_iMaxPeaks; }
		void SetZoom(t_real_glob dZoom);
		t_real_glob GetZoom() const { return m_dZoom; }

		void SetWSVisible(bool bVisible);
		const tl::Brillouin3D<t_real_glob>& GetWS3D() const { return m_ws3; }
		const std::vector<ublas::vector<t_real_glob>>& GetWS3DPlaneVerts() const { return m_vecWS3VertsUnproj; }

		const tl::Kd<t_real_glob>& GetKdLattice() const { return m_kdLattice; }

	public:
		t_real_glob GetScaleFactor() const { return m_dScaleFactor; }
		void SetScaleFactor(t_real_glob dScale) { m_dScaleFactor = dScale; }

		ublas::vector<t_real_glob> GetHKLFromPlanePos(t_real_glob x, t_real_glob y) const;
		const tl::Lattice<t_real_glob>& GetRealLattice() const { return m_lattice; }
};


class LatticeScene : public QGraphicsScene
{	Q_OBJECT
	protected:
		RealLattice *m_pLatt = nullptr;
		bool m_bSnap = false;
		bool m_bMousePressed = false;

	public:
		LatticeScene(QObject *pParent = nullptr);
		virtual ~LatticeScene();

		const RealLattice* GetLattice() const { return m_pLatt; }
		RealLattice* GetLattice() { return m_pLatt; }

	public slots:
		void scaleChanged(t_real_glob dTotalScale);

	signals:
		void coordsChanged(t_real_glob dh, t_real_glob dk, t_real_glob dl,
			bool bHasNearest,
			t_real_glob dNearestH, t_real_glob dNearestK, t_real_glob dNearestL);

	protected:
		virtual void mousePressEvent(QGraphicsSceneMouseEvent *pEvt) override;
		virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *pEvt) override;
		virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *pEvt) override;

		virtual void keyPressEvent(QKeyEvent *pEvt) override;
		virtual void keyReleaseEvent(QKeyEvent *pEvt) override;

		virtual void drawBackground(QPainter*, const QRectF&) override;
};


class LatticeView : public QGraphicsView
{
	Q_OBJECT
	protected:
		t_real_glob m_dTotalScale = 1.;

		void DoZoom(t_real_glob delta);
		virtual void wheelEvent(QWheelEvent* pEvt) override;
		virtual void keyPressEvent(QKeyEvent *pEvt) override;
		virtual void keyReleaseEvent(QKeyEvent *pEvt) override;

	public:
		LatticeView(QWidget* pParent = nullptr);
		virtual ~LatticeView();

	signals:
		void scaleChanged(t_real_glob dTotalScale);
};

#endif
