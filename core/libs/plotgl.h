/**
 * gl plotter without threading
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 19-may-2013 -- jan-2019
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TAKIN_PLOT_GL_2__
#define __TAKIN_PLOT_GL_2__

#include "plotgl.h"
#include <QTimer>
#include <QtWidgets>
#include <QOpenGLWidget>

#include <vector>

#include "tlibs/gfx/gl.h"
#include "tlibs/gfx/gl_font.h"

#include "libs/globals.h"
#include "libs/globals_qt.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace ublas = boost::numeric::ublas;

#include <boost/signals2.hpp>
namespace sig = boost::signals2;


/**
 * types of plottable objects
 */
enum PlotTypeGl
{
	PLOT_INVALID,

	PLOT_SPHERE,
	PLOT_ELLIPSOID,

	PLOT_POLY,
	PLOT_LINES,
};


/**
 * plottable object
 */
struct PlotObjGl
{
	PlotTypeGl plttype = PLOT_INVALID;

	ublas::vector<t_real_glob> vecPos;
	ublas::vector<t_real_glob> vecScale;
	std::vector<t_real_glob> vecRotMat;

	std::vector<t_real_glob> vecColor;
	t_real_glob dLineWidth = 2.;

	std::vector<ublas::vector<t_real_glob>> vecVertices;
	ublas::vector<t_real_glob> vecNorm;

	bool bSelected = false;
	bool bUseLOD = true;
	bool bCull = true;

	bool bAnimated = false;
	t_real_glob dScaleMult = 1.;

	std::string strLabel;
};


struct PlotGlSize
{
	int iW = 800, iH = 600;
	t_real_glob dDPIScale = 1.;
	bool bDoResize = true;
};


class PlotGl : public QOpenGLWidget, public QTimer
{
public:
	using t_sigHover = sig::signal<void(const PlotObjGl*)>;


protected:
	QSettings *m_pSettings = nullptr;
	bool m_bEnabled;

	static constexpr t_real_glob m_dFOV = 45./180.*M_PI;
	tl::t_mat4_gen<t_real_glob> m_matProj, m_matView;
	bool m_bPerspective = true;	// perspective or orthogonal projection?
	bool m_bResetPrespective = false;
	ublas::vector<t_real_glob> m_vecCam;
	bool m_bDoZTest = false;
	bool m_bDrawPolys = true;
	bool m_bDrawLines = true;
	bool m_bDrawSpheres = true;

	tl::GlFontMap<t_real_glob> *m_pFont = nullptr;

	std::vector<PlotObjGl> m_vecObjs;
	GLuint m_iLstSphere[4];
	QString m_strLabels[3];

	std::size_t m_iPrec = 6;
	bool m_bDrawMinMax = true;

	t_real_glob m_dTime = 0.;

	t_sigHover m_sigHover;

	t_real_glob m_dXMin=-10., m_dXMax=10.;
	t_real_glob m_dYMin=-10., m_dYMax=10.;
	t_real_glob m_dZMin=-10., m_dZMax=10.;
	t_real_glob m_dXMinMaxOffs, m_dYMinMaxOffs, m_dZMinMaxOffs;

	// mouse stuff
	bool m_bMouseRotateActive = false;
	t_real_glob m_dMouseRot[2];
	t_real_glob m_dMouseBegin[2];

	bool m_bMouseScaleActive = false;
	t_real_glob m_dMouseScale;
	t_real_glob m_dMouseScaleBegin;

	PlotGlSize m_size;


protected:
	virtual void timerEvent(QTimerEvent *pEvt) override;

	void SetColor(t_real_glob r, t_real_glob g, t_real_glob b, t_real_glob a = 1.);
	void SetColor(std::size_t iIdx);

	virtual void initializeGL() override;
	virtual void resizeGL(int w, int h) override;
	virtual void paintGL() override;
	//virtual void paintEvent(QPaintEvent*) override;

	void updateViewMatrix();
	void mouseSelectObj(t_real_glob dX, t_real_glob dY);

	void SetPerspective(int w, int h);
	void freeGL();
	void tick(t_real_glob dTime);

	virtual void mousePressEvent(QMouseEvent*) override;
	virtual void mouseReleaseEvent(QMouseEvent*) override;
	virtual void mouseMoveEvent(QMouseEvent*) override;
	virtual void wheelEvent(QWheelEvent*) override;

	t_real_glob GetCamObjDist(const PlotObjGl& obj) const;
	std::vector<std::size_t> GetObjSortOrder() const;


public:
	PlotGl(QWidget* pParent, QSettings *pSettings=nullptr, t_real_glob dMouseScale = 25.);
	virtual ~PlotGl();

	virtual void AddHoverSlot(const typename t_sigHover::slot_type& conn);

	virtual void clear();
	virtual void TogglePerspective();
	virtual void ToggleZTest() { m_bDoZTest = !m_bDoZTest; }
	virtual void ToggleDrawPolys() { m_bDrawPolys = !m_bDrawPolys; }
	virtual void ToggleDrawLines() { m_bDrawLines = !m_bDrawLines; }
	virtual void ToggleDrawSpheres() { m_bDrawSpheres = !m_bDrawSpheres; }

	virtual void PlotSphere(const ublas::vector<t_real_glob>& vecPos, t_real_glob dRadius, int iObjIdx = -1);
	virtual void PlotEllipsoid(const ublas::vector<t_real_glob>& widths,
		const ublas::vector<t_real_glob>& offsets,
		const ublas::matrix<t_real_glob>& rot,
		int iObjsIdx=-1);
	virtual void PlotPoly(const std::vector<ublas::vector<t_real_glob>>& vecVertices,
		const ublas::vector<t_real_glob>& vecNorm, int iObjIdx=-1);
	virtual void PlotLines(const std::vector<ublas::vector<t_real_glob>>& vecVertices,
		t_real_glob dLW=2., int iObjIdx=-1);

	virtual void SetObjectCount(std::size_t iSize) { m_vecObjs.resize(iSize); }
	virtual void SetObjectColor(std::size_t iObjIdx, const std::vector<t_real_glob>& vecCol);
	virtual void SetObjectLabel(std::size_t iObjIdx, const std::string& strLab);
	virtual void SetObjectUseLOD(std::size_t iObjIdx, bool bLOD);
	virtual void SetObjectCull(std::size_t iObjIdx, bool bCull);
	virtual void SetObjectAnimation(std::size_t iObjIdx, bool bAnimate);

	virtual void SetLabels(const char* pcLabX, const char* pcLabY, const char* pcLabZ);
	virtual void SetDrawMinMax(bool b) { m_bDrawMinMax = b; }

	virtual void SetEnabled(bool b);
	virtual void SetPrec(std::size_t iPrec) { m_iPrec = iPrec; }

	virtual void keyPressEvent(QKeyEvent *pEvt) override;


	template<class t_vec>
	/*virtual*/ void SetMinMax(const t_vec& vecMin, const t_vec& vecMax, const t_vec* pOffs = nullptr)
	{
		m_dXMin = vecMin[0]; m_dXMax = vecMax[0];
		m_dYMin = vecMin[1]; m_dYMax = vecMax[1];
		m_dZMin = vecMin[2]; m_dZMax = vecMax[2];

		m_dXMinMaxOffs =  pOffs ? (*pOffs)[0] : 0.;
		m_dYMinMaxOffs =  pOffs ? (*pOffs)[1] : 0.;
		m_dZMinMaxOffs =  pOffs ? (*pOffs)[2] : 0.;
	}

	template<class t_vec=ublas::vector<t_real_glob>>
	/*virtual*/ void SetMinMax(const t_vec& vec, const t_vec* pOffs = nullptr)
	{
		m_dXMin = -vec[0]; m_dXMax = vec[0];
		m_dYMin = -vec[1]; m_dYMax = vec[1];
		m_dZMin = -vec[2]; m_dZMax = vec[2];

		m_dXMinMaxOffs =  pOffs ? (*pOffs)[0] : 0.;
		m_dYMinMaxOffs =  pOffs ? (*pOffs)[1] : 0.;
		m_dZMinMaxOffs =  pOffs ? (*pOffs)[2] : 0.;
	}
};


#endif
