/**
 * gl plotter, non-threaded
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

#include "plotgl.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/math.h"
#include "tlibs/math/geo_prim.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/flags.h"

#include <time.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>


#define RENDER_FPS 40

using t_real = t_real_glob;
using t_mat4 = tl::t_mat4_gen<t_real>;
using t_mat3 = tl::t_mat3_gen<t_real>;
using t_vec4 = tl::t_vec4_gen<t_real>;
using t_vec3 = tl::t_vec3_gen<t_real>;


#define POS_F localPos


static const int g_iTimerInterval = int(t_real(1e3) / t_real(RENDER_FPS));

// ----------------------------------------------------------------------------


PlotGl::PlotGl(QWidget* pParent, QSettings *pSettings, t_real dMouseScale) :
	QOpenGLWidget(pParent), m_pSettings(pSettings),
	m_bEnabled(true), m_matProj(tl::unit_m<t_mat4>(4)), m_matView(tl::unit_m<t_mat4>(4))
{
	m_dMouseRot[0] = m_dMouseRot[1] = 0.;
	m_dMouseScale = dMouseScale;
	updateViewMatrix();

	setUpdatesEnabled(true);
	QTimer::setSingleShot(false);
}


PlotGl::~PlotGl()
{
	SetEnabled(false);
}
// ----------------------------------------------------------------------------


void PlotGl::SetEnabled(bool b)
{
	m_bEnabled = b;

	// stop updating if widget is disabled
	if(b)
		QTimer::start(g_iTimerInterval);
	else
		QTimer::stop();
}


void PlotGl::SetColor(t_real r, t_real g, t_real b, t_real a)
{
	t_real pfCol[] = {r, g, b, a};
	tl::gl_traits<t_real>::SetMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, pfCol);
}


void PlotGl::SetColor(std::size_t iIdx)
{
	static const t_real cols[4][4] =
	{
		{ 0., 0., 1., 0.7 },
		{ 0., 0.5, 0., 0.7 },
		{ 1., 0., 0., 0.7 },
		{ 0., 0., 0., 0.7 }
	};

	tl::gl_traits<t_real>::SetMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, cols[iIdx % 4]);
}


void PlotGl::initializeGL()
{
	QSurfaceFormat form = format();
	tl::log_debug("Renderer started using GL version ", form.majorVersion(), ".", form.minorVersion(), ".");

	glClearColor(1.,1.,1.,0.);
	glShadeModel(GL_SMOOTH);
	//glShadeModel(GL_FLAT);
	glDisable(GL_NORMALIZE);

	glClearDepth(1.);
	glDepthMask(GL_TRUE);
	glDepthFunc(GL_LEQUAL);
	glDisable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	//glEnable(GL_MULTISAMPLE);
	//glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glFrontFace(GL_CCW);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);

	const t_real vecLightColDiff[] = { 1., 1., 1., 1. };
	const t_real vecLightColAmb[] = { 0.25, 0.25, 0.25, 1. };
	const t_real vecLightColSpec[] = { 0., 0., 0., 1. };

	t_real dLightDir = 1./std::sqrt(t_real(3));
	const t_real vecLight0[] = { dLightDir, dLightDir, dLightDir, t_real(0) };
	tl::gl_traits<t_real>::SetLight(GL_LIGHT0, GL_AMBIENT, vecLightColAmb);
	tl::gl_traits<t_real>::SetLight(GL_LIGHT0, GL_DIFFUSE, vecLightColDiff);
	tl::gl_traits<t_real>::SetLight(GL_LIGHT0, GL_SPECULAR, vecLightColSpec);
	tl::gl_traits<t_real>::SetLight(GL_LIGHT0, GL_POSITION, vecLight0);


	// generate spheres in several detail levels
	for(std::size_t iSphere=0; iSphere<sizeof(m_iLstSphere)/sizeof(*m_iLstSphere); ++iSphere)
	{
		m_iLstSphere[iSphere] = glGenLists(1);
		tl::TesselSphere<t_vec3> prim(t_real(1), iSphere);

		glNewList(m_iLstSphere[iSphere], GL_COMPILE);
			for(std::size_t iPoly=0; iPoly<prim.GetPolyCount(); ++iPoly)
			{
				glBegin(GL_POLYGON);
					//t_vec3 vecFaceNorm = prim.GetPolyNormal(iPoly);
					//tl::gl_traits<t_real>::SetNorm(vecFaceNorm[0], vecFaceNorm[1], vecFaceNorm[2]);
					for(const t_vec3& vec : prim.GetPoly(iPoly))
					{
						t_vec3 vecNorm = vec / ublas::norm_2(vec);
						tl::gl_traits<t_real>::SetNorm(vecNorm[0], vecNorm[1], vecNorm[2]);
						tl::gl_traits<t_real>::SetVertex(vec[0], vec[1], vec[2]);
					}
				glEnd();
			}
		glEndList();
	}

	//glActiveTexture(GL_TEXTURE0);

	if(g_strFontGL == "")
		g_strFontGL = DEF_FONT;
	if(g_iFontGLSize <= 0)
		g_iFontGLSize = DEF_FONT_SIZE;

	std::string strFont = find_resource(g_strFontGL);
	m_pFont = new tl::GlFontMap<t_real>(strFont.c_str(), g_iFontGLSize);

	QTimer::start(g_iTimerInterval);
	QWidget::setMouseTracking(true);
}


void PlotGl::freeGL()
{
	SetEnabled(false);
	QWidget::setMouseTracking(false);

	for(std::size_t iSphere=0; iSphere<sizeof(m_iLstSphere)/sizeof(*m_iLstSphere); ++iSphere)
		glDeleteLists(m_iLstSphere[iSphere], 1);

	if(m_pFont) { delete m_pFont; m_pFont = nullptr; }
}


void PlotGl::SetPerspective(int w, int h)
{
	glMatrixMode(GL_PROJECTION);
	if(m_bPerspective)
		m_matProj = tl::perspective_matrix<t_mat4, t_real>(m_dFOV, t_real(w)/t_real(h), 0.1, 100.);
	else
		m_matProj = tl::ortho_matrix<t_mat4, t_real>(-1.,1.,-1.,1.,0.1,100.);
	t_real glmat[16]; tl::to_gl_array(m_matProj, glmat);
	tl::gl_traits<t_real>::LoadMatrix(glmat);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


void PlotGl::resizeGL(int w, int h)
{
	if(w <= 0)
		w = 1;
	if(h <= 0)
		h = 1;

	m_size.iW = w;
	m_size.iH = h;
	m_size.dDPIScale = devicePixelRatioF();
	//tl::log_debug("gl resize: ", m_size.iW, "x", m_size.iH, ", dpi: ", m_size.dDPIScale);

	glViewport(0, 0, w, h);
	SetPerspective(w, h);
}


/**
 * average distance of object from camera position
 */
t_real PlotGl::GetCamObjDist(const PlotObjGl& obj) const
{
	ublas::vector<t_real> vecPos;
	t_real dShiftPt = 0.;

	if(obj.plttype == PLOT_POLY || obj.plttype == PLOT_LINES)
	{
		// mean position
		vecPos = tl::mean_value(obj.vecVertices);
		//if(obj.plttype == PLOT_LINES)
		//	dShiftPt = obj.dLineWidth*0.25;
	}
	else if(obj.plttype == PLOT_SPHERE || obj.plttype == PLOT_ELLIPSOID)
	{
		vecPos = obj.vecPos;
		// point on sphere closest to camera
		dShiftPt = ublas::norm_2(obj.vecScale);
	}
	else
	{
		return t_real(0);
	}

	// shift position towards camera
	if(!tl::float_equal(dShiftPt, t_real(0)))
	{
		ublas::vector<t_real> vecCamObj = m_vecCam - vecPos;
		vecCamObj /= ublas::norm_2(vecCamObj);
		vecPos += vecCamObj*dShiftPt;
	}

	if(vecPos.size() != m_vecCam.size())
		return t_real(0);
	return ublas::norm_2(vecPos - m_vecCam);
}


/**
 * get sorted object indices based on distance to camera
 */
std::vector<std::size_t> PlotGl::GetObjSortOrder() const
{
	std::vector<std::size_t> vecIdx(m_vecObjs.size());
	std::iota(vecIdx.begin(), vecIdx.end(), 0);

	// sort indices based on obj distance
	std::stable_sort(vecIdx.begin(), vecIdx.end(),
		[this](const std::size_t& iIdx0, const std::size_t& iIdx1) -> bool
		{
			if(iIdx0 < m_vecObjs.size() && iIdx1 < m_vecObjs.size())
			{
				t_real tDist0 = GetCamObjDist(m_vecObjs[iIdx0]);
				t_real tDist1 = GetCamObjDist(m_vecObjs[iIdx1]);
				return m_bDoZTest ? tDist0 < tDist1 : tDist0 >= tDist1;
			}
			else
			{
				tl::log_err("Possible corruption in GL object list.");
				return false;
			}
		});

	return vecIdx;
}


void PlotGl::timerEvent(QTimerEvent *pEvt)
{
	QTimer::timerEvent(pEvt);

	m_dTime += t_real(1)/t_real(RENDER_FPS);
	tick(m_dTime);
}


/**
 * set absolute time to "dTime" seconds
 */
void PlotGl::tick(t_real dTime)
{
	if(!m_bEnabled) return;
	//tl::log_debug("tick: ", dTime);

	// cycle between dNum1 and dNum2
	auto fktCycle = [](t_real dTime, t_real dNum1, t_real dNum2) -> t_real
	{
		dNum2 -= dNum1;
		return (dNum1 + std::abs(std::fmod(dTime,dNum2*2.)-dNum2));
	};


	{	// look for objects with animation
		for(PlotObjGl& obj : m_vecObjs)
		{
			if(!obj.bAnimated)
				continue;
			obj.dScaleMult = fktCycle(2.*dTime, 0.5, 2.);
		}
	}

	if(isVisible())
		update();
}


/**
 * main render function
 */
void PlotGl::paintGL()
{
	if(!m_bEnabled) return;

	//QPainter painter{this};
	//painter.setRenderHint(QPainter::Antialiasing);
	//painter.beginNativePainting();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// does the projection matrix have to be reset?
	if(m_bResetPrespective)
	{
		SetPerspective(m_size.iW, m_size.iH);
		m_bResetPrespective = false;
	}

	glMatrixMode(GL_MODELVIEW);
	t_real glmat[16];
	tl::to_gl_array(m_matView, glmat);
	tl::gl_traits<t_real>::LoadMatrix(glmat);


	// camera position
	tl::Line<t_real> rayMid = tl::screen_ray(t_real(0.5), t_real(0.5), m_matProj, m_matView);
	m_vecCam = rayMid.GetX0();


	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_TEXTURE_2D);

	if(m_bDoZTest)
	{
		glEnable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);
	}
	else
	{
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
	}


	// draw axes
	glPushMatrix();
	glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_LIGHTING_BIT | GL_DEPTH_BUFFER_BIT);
		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glLineWidth(2.);
		tl::gl_traits<t_real>::SetColor(0., 0., 0.);

		const t_real dAxisScale = 1.8;
		glBegin(GL_LINES);
			tl::gl_traits<t_real>::SetVertex(m_dXMin*dAxisScale, 0., 0.);
			tl::gl_traits<t_real>::SetVertex(m_dXMax*dAxisScale, 0., 0.);
			tl::gl_traits<t_real>::SetVertex(0., m_dYMin*dAxisScale, 0.);
			tl::gl_traits<t_real>::SetVertex(0., m_dYMax*dAxisScale, 0.);
			tl::gl_traits<t_real>::SetVertex(0., 0., m_dZMin*dAxisScale);
			tl::gl_traits<t_real>::SetVertex(0., 0., m_dZMax*dAxisScale);
		glEnd();
	glPopAttrib();
	glPopMatrix();


	// draw objects
	for(std::size_t iObjIdx : GetObjSortOrder())
	{
		//tl::log_debug("Plotting object ", iObjIdx);
		if(iObjIdx >= m_vecObjs.size())
			continue;
		const PlotObjGl& obj = m_vecObjs[iObjIdx];
		bool bIsSphereLikeObj = false;

		if(obj.bCull)
			glEnable(GL_CULL_FACE);
		else
			glDisable(GL_CULL_FACE);

		if(obj.bSelected)
		{
			SetColor(0.25, 0.25, 0.25, 0.9);
		}
		else
		{
			if(obj.vecColor.size())
				SetColor(obj.vecColor[0], obj.vecColor[1], obj.vecColor[2], obj.vecColor[3]);
			else
				SetColor(iObjIdx);
		}


		glPushMatrix();
		if(obj.plttype == PLOT_SPHERE)
		{
			if(m_bDrawSpheres)
			{
				tl::gl_traits<t_real>::SetTranslate(obj.vecPos[0], obj.vecPos[1], obj.vecPos[2]);
				tl::gl_traits<t_real>::SetScale(
					obj.vecScale[0] * obj.dScaleMult,
					obj.vecScale[0] * obj.dScaleMult,
					obj.vecScale[0] * obj.dScaleMult);

				bIsSphereLikeObj = true;
			}
		}
		else if(obj.plttype == PLOT_ELLIPSOID)
		{
			if(m_bDrawSpheres)
			{
				tl::gl_traits<t_real>::SetTranslate(obj.vecPos[0], obj.vecPos[1], obj.vecPos[2]);

				t_real dMatRot[] = {obj.vecRotMat[0], obj.vecRotMat[1], obj.vecRotMat[2], 0.,
					obj.vecRotMat[3], obj.vecRotMat[4], obj.vecRotMat[5], 0.,
					obj.vecRotMat[6], obj.vecRotMat[7], obj.vecRotMat[8], 0.,
					0., 0., 0., 1. };
				tl::gl_traits<t_real>::MultMatrix(dMatRot);
				tl::gl_traits<t_real>::SetScale(
					obj.vecScale[0] * obj.dScaleMult,
					obj.vecScale[1] * obj.dScaleMult,
					obj.vecScale[2] * obj.dScaleMult);

				bIsSphereLikeObj = true;
			}
		}
		else if(obj.plttype == PLOT_POLY)
		{
			if(m_bDrawPolys)
			{
				glBegin(GL_POLYGON);
					tl::gl_traits<t_real>::SetNorm(obj.vecNorm[0], obj.vecNorm[1], obj.vecNorm[2]);
					for(const ublas::vector<t_real>& vec : obj.vecVertices)
						tl::gl_traits<t_real>::SetVertex(vec[0], vec[1], vec[2]);
				glEnd();
			}
		}
		else if(obj.plttype == PLOT_LINES)
		{
			if(m_bDrawLines)
			{
				glLineWidth(obj.dLineWidth);
				glBegin(GL_LINE_LOOP);
					for(const ublas::vector<t_real>& vec : obj.vecVertices)
						tl::gl_traits<t_real>::SetVertex(vec[0], vec[1], vec[2]);
 				glEnd();
			}
		}
		else
		{
			tl::log_warn("Unknown plot object at index ", iObjIdx, ".");
		}


		if(bIsSphereLikeObj)
		{
			int iLODMax = sizeof(m_iLstSphere)/sizeof(*m_iLstSphere)-1;
			int iLOD = iLODMax;

			if(obj.bUseLOD)
			{
				t_real dLenDist = tl::gl_proj_sphere_size(/*dRadius*/1.);
				iLOD = dLenDist * 50.;
				iLOD = tl::clamp<int>(iLOD, 0, iLODMax);
			}

			glCallList(m_iLstSphere[iLOD]);
		}
		glPopMatrix();
	}


	// draw object labels (after all objects to make sure they are drawn last)
	for(std::size_t iObjIdx : GetObjSortOrder())
	{
		//tl::log_debug("Plotting object ", iObjIdx);
		if(iObjIdx >= m_vecObjs.size())
			continue;
		const PlotObjGl& obj = m_vecObjs[iObjIdx];

		// draw label of selected object
		if(obj.bSelected && obj.strLabel.length() && m_pFont && m_pFont->IsOk())
		{
			glPushMatrix();
			glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_LIGHTING_BIT | GL_DEPTH_BUFFER_BIT);

			m_pFont->BindTexture();
			tl::gl_traits<t_real>::SetTranslate(obj.vecPos[0], obj.vecPos[1], obj.vecPos[2]);
			tl::gl_traits<t_real>::SetColor(0., 0., 0., 1.);
			m_pFont->DrawText(0., 0., 0., obj.strLabel);

			glPopAttrib();
			glPopMatrix();
		}
	}


	// draw axis labels
	if(m_pFont && m_pFont->IsOk())
	{
		glPushMatrix();
		glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_LIGHTING_BIT | GL_DEPTH_BUFFER_BIT);
		m_pFont->BindTexture();

		tl::gl_traits<t_real>::SetColor(0., 0., 1., 1.);
		m_pFont->DrawText(m_dXMax*dAxisScale, 0., 0., m_strLabels[0].toStdString());
		m_pFont->DrawText(0., m_dYMax*dAxisScale , 0., m_strLabels[1].toStdString());
		m_pFont->DrawText(0., 0., m_dZMax*dAxisScale , m_strLabels[2].toStdString());

		if(m_bDrawMinMax)
		{
			tl::gl_traits<t_real>::SetColor(0., 0., 0., 1.);
			m_pFont->DrawText(m_dXMin, 0., 0., tl::var_to_str(m_dXMin+m_dXMinMaxOffs, m_iPrec));
			m_pFont->DrawText(m_dXMax, 0., 0., tl::var_to_str(m_dXMax+m_dXMinMaxOffs, m_iPrec));
			m_pFont->DrawText(0., m_dYMin, 0., tl::var_to_str(m_dYMin+m_dYMinMaxOffs, m_iPrec));
			m_pFont->DrawText(0., m_dYMax, 0., tl::var_to_str(m_dYMax+m_dYMinMaxOffs, m_iPrec));
			m_pFont->DrawText(0., 0., m_dZMin, tl::var_to_str(m_dZMin+m_dZMinMaxOffs, m_iPrec));
			m_pFont->DrawText(0., 0., m_dZMax, tl::var_to_str(m_dZMax+m_dZMinMaxOffs, m_iPrec));
		}
		glPopAttrib();
		glPopMatrix();
	}

	//painter.endNativePainting();
	//painter.end();
}


// ----------------------------------------------------------------------------


void PlotGl::clear()
{
	m_vecObjs.clear();
}


void PlotGl::SetObjectColor(std::size_t iObjIdx, const std::vector<t_real>& vecCol)
{
	if(m_vecObjs.size() <= iObjIdx)
		return;
	m_vecObjs[iObjIdx].vecColor = vecCol;
}


void PlotGl::SetObjectLabel(std::size_t iObjIdx, const std::string& strLab)
{
	if(m_vecObjs.size() <= iObjIdx)
		return;
	m_vecObjs[iObjIdx].strLabel = strLab;
}


void PlotGl::SetObjectUseLOD(std::size_t iObjIdx, bool bLOD)
{
	if(m_vecObjs.size() <= iObjIdx)
		return;
	m_vecObjs[iObjIdx].bUseLOD = bLOD;
}


void PlotGl::SetObjectCull(std::size_t iObjIdx, bool bCull)
{
	if(m_vecObjs.size() <= iObjIdx)
		return;
	m_vecObjs[iObjIdx].bCull = bCull;
}


void PlotGl::SetObjectAnimation(std::size_t iObjIdx, bool bAnim)
{
	if(m_vecObjs.size() <= iObjIdx)
		return;
	m_vecObjs[iObjIdx].bAnimated = bAnim;
}


// ----------------------------------------------------------------------------


void PlotGl::PlotSphere(const ublas::vector<t_real>& vecPos,
	t_real dRadius, int iObjIdx)
{
	if(iObjIdx < 0)
	{
		clear();
		iObjIdx = 0;
	}

	if(iObjIdx >= int(m_vecObjs.size()))
		m_vecObjs.resize(iObjIdx+1);
	PlotObjGl& obj = m_vecObjs[iObjIdx];

	obj.plttype = PLOT_SPHERE;
	obj.vecPos = vecPos;
	obj.vecScale = tl::make_vec<ublas::vector<t_real>>({ dRadius, dRadius, dRadius });
}


void PlotGl::PlotEllipsoid(const ublas::vector<t_real>& widths,
	const ublas::vector<t_real>& offsets, const ublas::matrix<t_real>& rot,
	int iObjIdx)
{
	if(iObjIdx < 0)
	{
		clear();
		iObjIdx = 0;
	}

	if(iObjIdx >= int(m_vecObjs.size()))
		m_vecObjs.resize(iObjIdx+1);
	PlotObjGl& obj = m_vecObjs[iObjIdx];

	obj.plttype = PLOT_ELLIPSOID;
	obj.vecScale = tl::make_vec<ublas::vector<t_real>>({ widths[0], widths[1], widths[2] });
	obj.vecPos = offsets;

	if(obj.vecRotMat.size() != 9)
		obj.vecRotMat.resize(9);

	std::size_t iNum = 0;
	for(std::size_t i = 0; i < 3; ++i)
		for(std::size_t j = 0; j < 3; ++j)
			obj.vecRotMat[iNum++] = rot(j,i);
}


void PlotGl::PlotPoly(const std::vector<ublas::vector<t_real>>& vecVertices,
	const ublas::vector<t_real>& vecNorm, int iObjIdx)
{
	if(iObjIdx < 0)
	{
		clear();
		iObjIdx = 0;
	}

	if(iObjIdx >= int(m_vecObjs.size()))
		m_vecObjs.resize(iObjIdx+1);
	PlotObjGl& obj = m_vecObjs[iObjIdx];

	obj.plttype = PLOT_POLY;
	obj.vecVertices = vecVertices;
	obj.vecNorm = vecNorm;
}


void PlotGl::PlotLines(const std::vector<ublas::vector<t_real>>& vecVertices,
	t_real_glob dLW, int iObjIdx)
{
	if(iObjIdx < 0)
	{
		clear();
		iObjIdx = 0;
	}

	if(iObjIdx >= int(m_vecObjs.size()))
		m_vecObjs.resize(iObjIdx+1);
	PlotObjGl& obj = m_vecObjs[iObjIdx];

	obj.plttype = PLOT_LINES;
	obj.vecVertices = vecVertices;
	obj.dLineWidth = dLW;
}


// ----------------------------------------------------------------------------


void PlotGl::mousePressEvent(QMouseEvent *event)
{
	if(event->buttons() & Qt::LeftButton)
	{
		m_bMouseRotateActive = true;
		m_dMouseBegin[0] = event->POS_F().x();
		m_dMouseBegin[1] = event->POS_F().y();
	}

	if(event->buttons() & Qt::RightButton)
	{
		m_bMouseScaleActive = true;
		m_dMouseScaleBegin = event->POS_F().y();
	}
}


void PlotGl::mouseReleaseEvent(QMouseEvent *event)
{
	if((event->buttons() & Qt::LeftButton) == 0)
		m_bMouseRotateActive = false;

	if((event->buttons() & Qt::RightButton) == 0)
		m_bMouseScaleActive = false;
}


void PlotGl::mouseMoveEvent(QMouseEvent *pEvt)
{
	t_real dNewX = t_real(pEvt->POS_F().x());
	t_real dNewY = t_real(pEvt->POS_F().y());

	bool bUpdateView = false;
	if(m_bMouseRotateActive)
	{
		m_dMouseRot[0] += dNewX - m_dMouseBegin[0];
		m_dMouseRot[1] += dNewY - m_dMouseBegin[1];

		m_dMouseBegin[0] = dNewX;
		m_dMouseBegin[1] = dNewY;

		bUpdateView = true;
	}

	if(m_bMouseScaleActive)
	{
		m_dMouseScale *= 1.-(dNewY - m_dMouseScaleBegin)/t_real(height()) * 2.;
		m_dMouseScaleBegin = dNewY;

		bUpdateView = true;
	}

	if(bUpdateView)
		updateViewMatrix();


	// scale coordinates to [-1 .. 1] range
	t_real dMouseX = 2.*dNewX /**m_size.dDPIScale*/ /t_real(m_size.iW) - 1.;
	t_real dMouseY = -(2.*dNewY /**m_size.dDPIScale*/ /t_real(m_size.iH) - 1.);
	//tl::log_debug("mouse: ", dNewX, " ", dNewY, ", scaled mouse: ", dMouseX, " ", dMouseY);

	bool bHasSelected = false;
	if(m_bEnabled)
	{
		mouseSelectObj(dMouseX, dMouseY);

		for(PlotObjGl& obj : m_vecObjs)
		{
			if(obj.bSelected)
			{
				m_sigHover(&obj);
				bHasSelected = true;
				break;
			}
		}
	}

	if(!bHasSelected)
		m_sigHover(nullptr);
}


void PlotGl::wheelEvent(QWheelEvent* pEvt)
{
	const t_real dDelta = pEvt->angleDelta().y()/8. / 150.;

	m_dMouseScale *= std::pow(2., dDelta);;
	updateViewMatrix();

	QWidget::wheelEvent(pEvt);
}


void PlotGl::TogglePerspective()
{
	m_bPerspective = !m_bPerspective;
	m_bResetPrespective = true;
}


void PlotGl::ToggleZTest()
{
	m_bDoZTest = !m_bDoZTest;
}


void PlotGl::ToggleDrawPolys()
{
	m_bDrawPolys = !m_bDrawPolys;
}


void PlotGl::ToggleDrawLines()
{
	m_bDrawLines = !m_bDrawLines;
}


void PlotGl::ToggleDrawSpheres()
{
	m_bDrawSpheres = !m_bDrawSpheres;
}



// ----------------------------------------------------------------------------


void PlotGl::updateViewMatrix()
{
	t_mat4 matScale = tl::make_mat<t_mat4>(
		{{m_dMouseScale,              0,             0, 0},
		{            0,   m_dMouseScale,             0, 0},
		{            0,               0, m_dMouseScale, 0},
		{            0,               0,             0, 1}});

	t_mat4 matR0 = tl::rotation_matrix_3d_z(tl::d2r<t_real>(m_dMouseRot[0]));
	t_mat4 matR1 = tl::rotation_matrix_3d_x(tl::d2r<t_real>(-90. + m_dMouseRot[1]));
	matR0.resize(4,4,1); matR1.resize(4,4,1);
	matR0(3,3) = matR1(3,3) = 1.;
	for(short i=0; i<3; ++i)
		matR0(i,3) = matR0(3,i) = matR1(i,3) = matR1(3,i) = 0.;
	t_mat4 matRot0 = matR0, matRot1 = matR1;

	t_mat4 matTrans = tl::make_mat<t_mat4>(
		{{ 1, 0, 0,  0},
		{  0, 1, 0,  0},
		{  0, 0, 1, -2},
		{  0, 0, 0,  1}});

	{
		m_matView = ublas::prod(matTrans, matRot1);
		m_matView = ublas::prod(m_matView, matRot0);
		m_matView = ublas::prod(m_matView, matScale);
	}
}


void PlotGl::AddHoverSlot(const typename PlotGl::t_sigHover::slot_type& conn)
{
	m_sigHover.connect(conn);
}


void PlotGl::mouseSelectObj(t_real dX, t_real dY)
{
	tl::Line<t_real> ray = tl::screen_ray(dX, dY, m_matProj, m_matView);

	for(PlotObjGl& obj : m_vecObjs)
	{
		obj.bSelected = false;

		std::unique_ptr<tl::Quadric<t_real>> pQuad;
		t_vec3 vecOffs = ublas::zero_vector<t_real>(3);

		if(obj.plttype == PLOT_SPHERE)
		{
			pQuad.reset(new tl::QuadSphere<t_real>(
				obj.vecScale[0] * obj.dScaleMult));
			vecOffs = obj.vecPos;
		}
		else if(obj.plttype == PLOT_ELLIPSOID)
		{
			pQuad.reset(new tl::QuadEllipsoid<t_real>(
				obj.vecScale[0] * obj.dScaleMult,
				obj.vecScale[1] * obj.dScaleMult,
				obj.vecScale[2] * obj.dScaleMult));

			vecOffs = obj.vecPos;
			t_mat3 matRot = tl::make_mat<t_mat3>(
				{{obj.vecRotMat[0], obj.vecRotMat[1], obj.vecRotMat[2]},
				 {obj.vecRotMat[3], obj.vecRotMat[4], obj.vecRotMat[5]},
				 {obj.vecRotMat[6], obj.vecRotMat[7], obj.vecRotMat[8]}});
			pQuad->transform(matRot);
		}

		std::vector<t_real> vecT;
		if(pQuad)
		{
			pQuad->SetOffset(vecOffs);
			vecT = pQuad->intersect(ray);
		}

		if(vecT.size() > 0)
		{
			for(t_real t : vecT)
			{
				if(t < 0.) continue; // beyond "near" plane
				if(t > 1.) continue; // beyond "far" plane
				obj.bSelected = true;
			}
		}
	}
}


void PlotGl::SetLabels(const char* pcLabX, const char* pcLabY, const char* pcLabZ)
{
	m_strLabels[0] = pcLabX;
	m_strLabels[1] = pcLabY;
	m_strLabels[2] = pcLabZ;
}


void PlotGl::keyPressEvent(QKeyEvent* pEvt)
{
	if(pEvt->key() == Qt::Key_Space)
		TogglePerspective();
	else if(pEvt->key() == Qt::Key_Z)
		ToggleZTest();
	else if(pEvt->key() == Qt::Key_P)
		ToggleDrawPolys();
	else if(pEvt->key() == Qt::Key_L)
		ToggleDrawLines();
	else if(pEvt->key() == Qt::Key_S)
		ToggleDrawSpheres();

	QOpenGLWidget::keyPressEvent(pEvt);
}
