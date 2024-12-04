/**
 * tlibs2 -- common gl functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2017-2021
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - http://doc.qt.io/qt-5/qopenglwidget.html#details
 *   - http://code.qt.io/cgit/qt/qtbase.git/tree/examples/opengl/threadedqopenglwidget
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAG_GL_COMMON_H__
#define __MAG_GL_COMMON_H__


#include <QtCore/QtGlobal>

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
	#include <QtOpenGL/QOpenGLShaderProgram>
	#include <QtOpenGL/QOpenGLVertexArrayObject>
	#include <QtOpenGL/QOpenGLBuffer>
	#include <QtOpenGL/QOpenGLTexture>
	#include <QtOpenGL/QOpenGLFramebufferObject>
	#include <QtOpenGL/QOpenGLFramebufferObjectFormat>
	#include <QtOpenGL/QOpenGLPaintDevice>
	#include <QtOpenGLWidgets/QOpenGLWidget>
#else
	#include <QtGui/QOpenGLShaderProgram>
	#include <QtGui/QOpenGLVertexArrayObject>
	#include <QtGui/QOpenGLBuffer>
	#include <QtGui/QOpenGLTexture>
	#include <QtGui/QOpenGLFramebufferObject>
	#include <QtGui/QOpenGLFramebufferObjectFormat>
	#include <QtGui/QOpenGLPaintDevice>
	#include <QtWidgets/QOpenGLWidget>
#endif

#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
	#include <QtCore/QRecursiveMutex>
	using t_qt_mutex = QRecursiveMutex;
#else
	#include <QtCore/QMutex>
	using t_qt_mutex = QMutex;
#endif

#include <QtGui/QMatrix4x4>
#include <QtGui/QVector4D>
#include <QtGui/QVector3D>
#include <QtGui/QVector2D>

#include <memory>
#include <optional>

#include "../maths.h"



// ----------------------------------------------------------------------------
// GL versions
#if !defined(_GL_MAJ_VER) || !defined(_GL_MIN_VER)
	#define _GL_MAJ_VER 3
	#define _GL_MIN_VER 3
#endif

#if _GL_MAJ_VER<=3 && _GL_MIN_VER<2
	#if !defined(_GL_SUFFIX)
		#define _GL_SUFFIX
	#endif

	#if _GL_MAJ_VER==3 && _GL_MIN_VER==1
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 4
	#elif _GL_MAJ_VER==3 && _GL_MIN_VER==0
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 3
	#elif _GL_MAJ_VER==2 && _GL_MIN_VER==1
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 2
	#elif _GL_MAJ_VER==2 && _GL_MIN_VER==0
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 1
	#endif
#else
	#if !defined(_GL_SUFFIX)
		#define _GL_SUFFIX _Core
	#endif

	#if _GL_MAJ_VER==3 && _GL_MIN_VER==2
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 5
	#else
		#define _GLSL_MAJ_VER _GL_MAJ_VER
		#define _GLSL_MIN_VER _GL_MIN_VER
	#endif
#endif


// GL functions include
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
	#define _GL_INC_IMPL(MAJ, MIN, SUFF) <QtOpenGL/QOpenGLFunctions_ ## MAJ ## _ ## MIN ## SUFF>
#else
	#define _GL_INC_IMPL(MAJ, MIN, SUFF) <QtGui/QOpenGLFunctions_ ## MAJ ## _ ## MIN ## SUFF>
#endif
#define _GL_INC(MAJ, MIN, SUFF) _GL_INC_IMPL(MAJ, MIN, SUFF)
#include _GL_INC(_GL_MAJ_VER, _GL_MIN_VER, _GL_SUFFIX)

// GL functions typedef
#define _GL_FUNC_IMPL(MAJ, MIN, SUFF) QOpenGLFunctions_ ## MAJ ## _ ## MIN ## SUFF
#define _GL_FUNC(MAJ, MIN, SUFF) _GL_FUNC_IMPL(MAJ, MIN, SUFF)
using qgl_funcs = _GL_FUNC(_GL_MAJ_VER, _GL_MIN_VER, _GL_SUFFIX);


// GL error codes: https://www.khronos.org/opengl/wiki/OpenGL_Error
#define LOGGLERR(pGl) { while(true) {	\
		auto err = pGl->glGetError();	\
		if(err == GL_NO_ERROR) break;	\
		std::cerr << "GL error in " << __func__ \
			<< ", file: " << __FILE__ \
			<< ", line " << std::dec <<  __LINE__  \
			<< ": " << std::hex << "0x" << err \
			<< "." << std::endl; \
	}}
// ----------------------------------------------------------------------------


namespace tl2 {

// ----------------------------------------------------------------------------
// types
#ifdef _T_REAL_GL
	using t_real_gl _T_REAL_GL;
#else
	using t_real_gl = GLfloat;
	//using t_real_gl = GLdouble;
#endif

#ifdef _T_VEC2_GL
	using t_vec2_gl _T_VEC2_GL;
#else
	using t_vec2_gl = tl2::qvecN_adapter<int, 2, t_real_gl, QVector2D>;
#endif

#ifdef _T_VEC3_GL
	using t_vec3_gl _T_VEC3_GL;
#else
	using t_vec3_gl = tl2::qvecN_adapter<int, 3, t_real_gl, QVector3D>;
#endif

#ifdef _T_VEC_GL
	using t_vec_gl _T_VEC_GL;
#else
	using t_vec_gl = tl2::qvecN_adapter<int, 4, t_real_gl, QVector4D>;
#endif

#ifdef _T_MAT_GL
	using t_mat_gl _T_MAT_GL;
#else
	using t_mat_gl = tl2::qmatNN_adapter<int, 4, 4, t_real_gl, QMatrix4x4>;
#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// plotter objects
enum class GlRenderObjType
{
	TRIANGLES,
	LINES
};


#ifndef GlPlotObjType
	#define GlPlotObjType GlRenderObjType
#endif


struct GlRenderObj
{
	GlPlotObjType m_type = GlPlotObjType::TRIANGLES;

	std::shared_ptr<QOpenGLVertexArrayObject> m_vertex_array{};
	std::shared_ptr<QOpenGLBuffer> m_vertex_buffer{};
	std::shared_ptr<QOpenGLBuffer> m_normals_buffer{};
	std::shared_ptr<QOpenGLBuffer> m_uv_buffer{};
	std::shared_ptr<QOpenGLBuffer> m_colour_buffer{};

	std::vector<t_vec3_gl> m_vertices{}, m_triangles{}, m_uvs{};

	t_vec_gl m_colour = tl2::create<t_vec_gl>({ 0., 0., 1., 1. });	// rgba
};


struct GlPlotObj : public GlRenderObj
{
	// does not define a geometry itself, but just links to another object
	std::optional<std::size_t> linkedObj{};

	t_mat_gl m_mat = tl2::unit<t_mat_gl>();

	bool m_invariant = false;	// invariant to A, B matrices
	bool m_visible = true;		// object shown?
	bool m_highlighted = false;	// object highlighted?
	bool m_valid = true;		// object deleted?
	bool m_intersect = true;        // object can be intersected by picker
	int m_priority = 1;		// object rendering priority

	t_vec3_gl m_labelPos = tl2::create<t_vec3_gl>({ 0., 0., 0. });
	std::string m_label{};
	std::string m_datastr{};

	t_vec3_gl m_boundingSpherePos = tl2::create<t_vec3_gl>({ 0., 0., 0. });
	t_real_gl m_boundingSphereRad = 0.;
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// functions
// GL surface format
extern void set_gl_format(bool bCore = true,
	int iMajorVer = 3, int iMinorVer = 3,
	int iSamples = 8);

extern QSurfaceFormat gl_format(bool bCore = true,
	int iMajorVer = 3, int iMinorVer = 3, int iSamples = 8,
	QSurfaceFormat surf = QSurfaceFormat::defaultFormat());

// get gl functions
extern qgl_funcs* get_gl_functions(QOpenGLWidget *pGLWidget);

// create a triangle-based object
extern bool create_triangle_object(QOpenGLWidget* pGLWidget, GlRenderObj& obj,
	const std::vector<t_vec3_gl>& verts, const std::vector<t_vec3_gl>& triagverts,
	const std::vector<t_vec3_gl>& norms, const std::vector<t_vec3_gl>& uvs,
	const t_vec_gl& colour, bool bUseVertsAsNorm, GLint attrVertex,
	GLint attrVertexNormal, GLint attrVertexcolour, GLint attrTextureCoords = -1);

// create a line-based object
extern bool create_line_object(QOpenGLWidget* pGLWidget, GlRenderObj& obj,
	const std::vector<t_vec3_gl>& verts, const t_vec_gl& colour,
	GLint attrVertex, GLint attrVertexcolour);


extern void delete_render_object(GlRenderObj& obj);
// ----------------------------------------------------------------------------

}
#endif
