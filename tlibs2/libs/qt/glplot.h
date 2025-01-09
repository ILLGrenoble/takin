/**
 * tlibs2 -- GL plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2017-2021
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "magtools" (https://github.com/t-weber/magtools).
 *
 * References:
 *   - http://doc.qt.io/qt-5/qopenglwidget.html#details
 *   - http://code.qt.io/cgit/qt/qtbase.git/tree/examples/opengl/threadedqopenglwidget
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * magtools
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __MAG_GL_PLOT_H__
#define __MAG_GL_PLOT_H__

#include <QtCore/QThread>
#include <QtCore/QTimer>
#include <QtWidgets/QDialog>
#include <QtGui/QMouseEvent>

#include <functional>
#include <memory>
#include <chrono>
#include <atomic>

#include "gl.h"
#include "../maths.h"
#include "../cam.h"


namespace tl2 {

// ----------------------------------------------------------------------------
class GlPlot;


/**
 * GL plot renderer
 */
class GlPlotRenderer : public QObject
{ Q_OBJECT
public:
	using t_cam = tl2::Camera<t_mat_gl, t_vec_gl, t_vec3_gl, t_real_gl>;


private:
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
	t_qt_mutex m_mutexObj{};
#else
	t_qt_mutex m_mutexObj{QMutex::Recursive};
#endif


protected:
	GlPlot *m_pPlot = nullptr;
	std::string m_strGlVer{}, m_strGlShaderVer{},
		m_strGlVendor{}, m_strGlRenderer{};

	std::shared_ptr<QOpenGLShaderProgram> m_pShaders{};

	GLint m_attrVertex = -1;
	GLint m_attrVertexNorm = -1;
	GLint m_attrVertexCol = -1;
	GLint m_uniConstCol = -1;
	GLint m_uniLightPos = -1;
	GLint m_uniNumActiveLights = -1;
	GLint m_uniMatrixProj = -1;
	GLint m_uniMatrixCam = -1;
	GLint m_uniMatrixCamInv = -1;
	GLint m_uniMatrixObj = -1;
	GLint m_uniMatrixA = -1;
	GLint m_uniMatrixB = -1;
	GLint m_uniIsRealSpace = -1;
	GLint m_uniCoordSys = -1;

	t_mat_gl m_matA = tl2::unit<t_mat_gl>();
	t_mat_gl m_matB = tl2::unit<t_mat_gl>();
	t_real_gl m_CoordMax = 2.5;       // extent of coordinate axes

	std::atomic<bool> m_platform_supported = true;
	std::atomic<bool> m_initialised = false;
	std::atomic<bool> m_viewport_needs_update = false;
	std::atomic<bool> m_picker_enabled = true;
	std::atomic<bool> m_picker_needs_update = false;
	std::atomic<bool> m_lights_need_update = false;
	std::atomic<bool> m_Btrafo_needs_update = false;
	std::atomic<bool> m_cull = true;
	std::atomic<bool> m_blend = false;
	std::atomic<bool> m_is_real_space = true;
	std::atomic<int> m_coordsys = 0;  // 0: orthogonal, 1: using crystal matrix
	t_real_gl m_pickerSphereRadius = 1;
	bool m_showLabels = true;

	std::vector<t_vec3_gl> m_lights{};
	std::vector<GlPlotObj> m_objs{};
	std::optional<std::size_t> m_coordCrossLab{}, m_coordCrossXtal{};

	QPointF m_posMouse{};
	QPointF m_posMouseRotationStart{}, m_posMouseRotationEnd{};
	bool m_in_rotation = false;
	bool m_restrict_cam_theta = true;

	QTimer m_timer{};
	t_cam m_cam{};


protected:
	qgl_funcs* GetGlFunctions()
	{
		return get_gl_functions((QOpenGLWidget*)m_pPlot);
	}

	void UpdateCam();
	void UpdatePicker();
	void UpdateLights();
	void UpdateBTrafo();
	void RequestPlotUpdate();

	void DoPaintGL(qgl_funcs *pGL);
	void DoPaintNonGL(QPainter &painter);

	void tick(const std::chrono::milliseconds& ms);


public:
	GlPlotRenderer(GlPlot *pPlot = nullptr);
	virtual ~GlPlotRenderer();

	GlPlotRenderer(const GlPlotRenderer&) = delete;
	const GlPlotRenderer& operator=(const GlPlotRenderer&) = delete;

	static constexpr bool m_isthreaded = false;
	static constexpr bool m_usetimer = false;

	QPointF GlToScreenCoords(const t_vec_gl& vec, bool *pVisible = nullptr) const;

	const t_cam& GetCamera() const
	{
		return m_cam;
	}

	t_cam& GetCamera()
	{
		return m_cam;
	}

	std::tuple<std::string, std::string, std::string, std::string>
	GetGlDescr() const
	{
		return std::make_tuple(
			m_strGlVer, m_strGlShaderVer,
			m_strGlVendor, m_strGlRenderer);
	}

	void SetPickerSphereRadius(t_real_gl rad)
	{
		m_pickerSphereRadius = rad;
	}

	GlPlotObj CreateTriangleObject(const std::vector<t_vec3_gl>& verts,
		const std::vector<t_vec3_gl>& triag_verts, const std::vector<t_vec3_gl>& norms,
		const t_vec_gl& colour, bool bUseVertsAsNorm=false);
	GlPlotObj CreateLineObject(const std::vector<t_vec3_gl>& verts, const t_vec_gl& colour);

	std::size_t GetNumObjects() const { return m_objs.size(); }
	void RemoveObject(std::size_t obj);
	std::size_t AddLinkedObject(std::size_t linkTo,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 0,
		t_real_gl r = 1, t_real_gl g = 1, t_real_gl b = 1, t_real_gl a = 1);
	std::size_t AddSphere(t_real_gl rad = 1,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 0,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddCylinder(t_real_gl rad = 1, t_real_gl h = 1,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 0,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddCone(t_real_gl rad = 1, t_real_gl h = 1,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 0,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddArrow(t_real_gl rad = 1, t_real_gl h = 1,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 0,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddPlane(t_real_gl nx = 0, t_real_gl ny = 0, t_real_gl nz = 1,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 1, t_real_gl size = 10,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddPatch(
		std::function<t_real_gl(t_real_gl, t_real_gl)> fkt,
		t_real_gl x = 0, t_real_gl y = 0, t_real_gl z = 0,
		t_real_gl w = 10, t_real_gl h = 10, std::size_t pts_x = 16, std::size_t pts_y = 16,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddTriangleObject(const std::vector<t_vec3_gl>& triag_verts,
		const std::vector<t_vec3_gl>& triag_norms,
		t_real_gl r = 0, t_real_gl g = 0, t_real_gl b = 0, t_real_gl a = 1);
	std::size_t AddCoordinateCross(t_real_gl min, t_real_gl max);

	void SetObjectMatrix(std::size_t idx, const t_mat_gl& mat);
	void SetObjectCol(std::size_t idx, t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a = 1);
	void SetObjectLabel(std::size_t idx, const std::string& label);
	void SetObjectDataString(std::size_t idx, const std::string& data);
	void SetObjectVisible(std::size_t idx, bool visible);
	void SetObjectIntersectable(std::size_t idx, bool intersect);
	void SetObjectHighlight(std::size_t idx, bool highlight);
	void SetObjectPriority(std::size_t idx, int prio);
	void SetObjectsHighlight(bool highlight);

	const t_mat_gl& GetObjectMatrix(std::size_t idx) const;
	const std::string& GetObjectLabel(std::size_t idx) const;
	const std::string& GetObjectDataString(std::size_t idx) const;
	bool GetObjectVisible(std::size_t idx) const;
	bool GetObjectHighlight(std::size_t idx) const;

	void SetScreenDims(int w, int h);

	void SetCoordMax(t_real_gl d)
	{
		m_CoordMax = d;
	}

	void SetLight(std::size_t idx, const t_vec3_gl& pos);

	void SetCull(bool b)
	{
		m_cull = b;
	}

	void SetBlend(bool b)
	{
		m_blend = b;
	}

	void SetRestrictCamTheta(bool b)
	{
		m_restrict_cam_theta = b;
	}

	void SetBTrafo(const t_mat_gl& matB, const t_mat_gl* matA = nullptr,
		bool is_real_space = true);
	void SetCoordSys(int iSys);

	bool IsInitialised() const
	{
		return m_initialised;
	}

	const QPointF& GetMousePosition() const
	{
		return m_posMouse;
	}

	void SetLabelsVisible(bool show)
	{
		m_showLabels = show;
	}

	std::optional<std::size_t> GetCoordCross(bool xtal = false) const
	{
		if(xtal)
			return m_coordCrossXtal;
		return m_coordCrossLab;
	}

	void RequestViewportUpdate();
	void UpdateViewport();


public slots:
	void paintGL();

	void startedThread();
	void stoppedThread();

	void initialiseGL();

	void mouseMoveEvent(const QPointF& pos);
	void zoom(t_real_gl val);
	void ResetZoom();

	void BeginRotation();
	void EndRotation();

	void EnablePicker(bool b);


protected slots:
	void tick();


signals:
	void PickerIntersection(const t_vec3_gl *pos,
		std::size_t objIdx, std::size_t triagIdx,
		const t_vec3_gl *posSphere);
	void CameraHasUpdated();
};



/**
 * GL plotter widget
 */
class GlPlot : public QOpenGLWidget
{ Q_OBJECT
public:
	GlPlot(QWidget *pParent = nullptr);
	virtual ~GlPlot();

	GlPlot(const GlPlot&) = delete;
	const GlPlot& operator=(const GlPlot&) = delete;

	GlPlotRenderer* GetRenderer()
	{
		return m_renderer.get();
	}

	static constexpr bool m_isthreaded = GlPlotRenderer::m_isthreaded;


protected:
	virtual void paintEvent(QPaintEvent*) override;
	virtual void initializeGL() override;
	virtual void paintGL() override;
	virtual void resizeGL(int w, int h) override;

	virtual void mouseMoveEvent(QMouseEvent *pEvt) override;
	virtual void mousePressEvent(QMouseEvent *Evt) override;
	virtual void mouseReleaseEvent(QMouseEvent *Evt) override;
	virtual void wheelEvent(QWheelEvent *pEvt) override;


private:
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
	mutable t_qt_mutex m_mutex{};
#else
	mutable t_qt_mutex m_mutex{QMutex::Recursive};
#endif

	std::unique_ptr<GlPlotRenderer> m_renderer{};
	std::unique_ptr<QThread> m_thread_impl{};
	bool m_mouseMovedBetweenDownAndUp = false;
	bool m_mouseDown[3] = { false, false, false };


public:
	t_qt_mutex* GetMutex() { return &m_mutex; }

	void MoveContextToThread();
	bool IsContextInThread() const;


protected slots:
	void beforeComposing();
	void afterComposing();
	void beforeResizing();
	void afterResizing();


signals:
	void AfterGLInitialisation();
	void GLInitialisationFailed();

	void MouseDown(bool left, bool mid, bool right);
	void MouseUp(bool left, bool mid, bool right);
	void MouseClick(bool left, bool mid, bool right);
};
// ----------------------------------------------------------------------------

}
#endif
