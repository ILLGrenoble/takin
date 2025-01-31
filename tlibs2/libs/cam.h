/**
 * camera for 3d plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date feb-2022
 * @note Code forked from my privately developed "misc" and "magtools" project (https://github.com/t-weber/misc and https://github.com/t-weber/magtools).
 * @note Class moved over to tlibs on 10-jan-2022 from my "TAS-Paths" project (https://code.ill.fr/scientific-software/takin/paths).
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - http://doc.qt.io/qt-5/qopenglwidget.html#details
 *   - http://code.qt.io/cgit/qt/qtbase.git/tree/examples/opengl/threadedqopenglwidget
 *   - http://doc.qt.io/qt-5/qtgui-openglwindow-example.html
 *   - http://doc.qt.io/qt-5/qopengltexture.html
 *   - (Sellers 2014) G. Sellers et al., ISBN: 978-0-321-90294-8 (2014).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * TAS-Paths (part of the Takin software suite)
 * Copyright (C) 2021       Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * magtools
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __GL_RENDERER_CAM_H__
#define __GL_RENDERER_CAM_H__

#include <tuple>
#include <array>
#include <unordered_set>
#include <algorithm>

#include "maths.h"


namespace tl2{

template<class t_mat, class t_vec, class t_vec3, class t_real>
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec3>
class Camera
{
public:
	Camera() = default;
	~Camera() = default;

	Camera(const Camera<t_mat, t_vec, t_vec3, t_real>&) = default;
	Camera<t_mat, t_vec, t_vec3, t_real>& operator=(
		const Camera<t_mat, t_vec, t_vec3, t_real>&) = default;


	/**
	 * centre camera on object matrix
	 */
	void Centre(const t_mat& objmat, bool set_z = true)
	{
		m_matTrans(0,3) = -objmat(0,3);
		m_matTrans(1,3) = -objmat(1,3);
		if(set_z)
			m_matTrans(2,3) = -objmat(2,3);

		m_trafo_needs_update = true;
	}


	/**
	 * set the camera's field of view
	 */
	void SetFOV(t_real angle)
	{
		m_FOV = angle;
		m_persp_needs_update = true;
	}


	/**
	 * get the camera's field of view
	 */
	t_real GetFOV() const
	{
		return m_FOV;
	}


	/**
	 * set the camera zoom
	 */
	void SetZoom(t_real zoom)
	{
		m_zoom = zoom;
		m_trafo_needs_update = true;
	}


	/**
	 * get the camera's zoom factor
	 */
	t_real GetZoom() const
	{
		return m_zoom;
	}


	/**
	 * set the camera's distance from the translation centre
	 */
	void SetDist(t_real dist)
	{
		m_dist = dist;
		m_trafo_needs_update = true;
	}


	/**
	 * get the camera's distance from the translation centre
	 */
	t_real GetDist() const
	{
		return m_dist;
	}


	/**
	 * set the camera frustum's near plane
	 */
	void SetNearPlane(t_real z)
	{
		m_nearPlane = z;
		m_persp_needs_update = true;
	}


	/**
	 * set the camera frustum's near plane
	 */
	void SetFarPlane(t_real z)
	{
		m_farPlane = z;
		m_persp_needs_update = true;
	}


	/**
	 * set the range for parallel projection
	 */
	void SetParalellRange(t_real range)
	{
		m_parallel_range = range;
		if(!GetPerspectiveProjection())
			m_persp_needs_update = true;
	}


	t_real GetParalellRange() const
	{
		return m_parallel_range;
	}


	/**
	 * set the camera position
	 */
	void SetPosition(const t_vec3& pos)
	{
		m_matTrans(0, 3) = pos[0];
		m_matTrans(1, 3) = pos[1];
		m_matTrans(2, 3) = pos[2];

		m_trafo_needs_update = true;
	}


	/**
	 * get the camera position
	 */
	t_vec3 GetPosition() const
	{
		return tl2::create<t_vec3>(
		{
			m_matTrans(0, 3),
			m_matTrans(1, 3),
			m_matTrans(2, 3),
		});
	}


	/**
	 * set the rotation angles
	 */
	void SetRotation(const t_real phi, t_real theta)
	{
		m_phi_saved = m_phi = phi;
		m_theta_saved = m_theta = theta;

		m_trafo_needs_update = true;
	}


	/**
	 * get the camera's rotation matrix
	 */
	std::tuple<t_real, t_real> GetRotation() const
	{
		return std::make_tuple(m_phi, m_theta);
	}


	/**
	 * save the current rotation angles
	 */
	void SaveRotation()
	{
		m_phi_saved = m_phi;
		m_theta_saved = m_theta;
	}


	/**
	 * set transformation matrix to look from "pos" to "target"
	 */
	void SetLookAt(const t_vec3& pos, const t_vec3& target, const t_vec3& up)
	{
		m_mat = tl2::hom_lookat<t_mat, t_vec3>(
			pos, target, up);

		std::tie(m_mat_inv, std::ignore)
			= tl2::inv<t_mat>(m_mat);

		m_trafo_needs_update = false;

		// TODO: extract position and angles corresponding to this matrix
	}


	/**
	 * rotate the camera by the given delta angles
	 */
	void Rotate(t_real dphi, t_real dtheta, bool restrict_theta = true)
	{
		m_phi = dphi + m_phi_saved;
		m_theta = dtheta + m_theta_saved;

		// wrap around phi angle
		m_phi = tl2::mod_pos<t_real>(m_phi, t_real(2)*tl2::pi<t_real>);

		// wrap around theta angle
		m_theta = tl2::mod_posneg<t_real>(m_theta, t_real(2)*tl2::pi<t_real>);

		// restrict theta angle
		if(restrict_theta)
			m_theta = std::clamp<t_real>(m_theta, -t_real(0.5)*tl2::pi<t_real>, 0.);

		m_trafo_needs_update = true;
	}


	/**
	 * translate the camera by the given deltas
	 */
	void Translate(t_real dx, t_real dy, t_real dz)
	{
		t_vec3 xdir = tl2::row<t_mat, t_vec3>(m_matRot, 0);
		t_vec3 ydir = tl2::row<t_mat, t_vec3>(m_matRot, 1);
		t_vec3 zdir = tl2::row<t_mat, t_vec3>(m_matRot, 2);

		t_vec3 xinc = xdir * dx;
		t_vec3 yinc = ydir * dy;
		t_vec3 zinc = zdir * dz;

		m_matTrans(0,3) += xinc[0] + yinc[0] + zinc[0];
		m_matTrans(1,3) += xinc[1] + yinc[1] + zinc[1];
		m_matTrans(2,3) += xinc[2] + yinc[2] + zinc[2];

		m_trafo_needs_update = true;
	}


	/**
	 * zoom disginguising between perspective or parallel projection
	 */
	void Zoom(t_real zoom)
	{
		const t_real factor = std::pow(t_real(2), zoom);

		if(GetPerspectiveProjection())
		{
			SetZoom(GetZoom() * factor);
			m_trafo_needs_update = true;
		}
		else
		{
			SetParalellRange(GetParalellRange() / factor);
			m_persp_needs_update = true;
		}
	}


	/**
	 * get the camera's full transformation matrix
	 */
	const t_mat& GetTransformation() const
	{
		return m_mat;
	}


	/**
	 * get the camera's inverse transformation matrix
	 */
	const t_mat& GetInverseTransformation() const
	{
		return m_mat_inv;
	}


	/**
	 * get the camera's full transformation matrix
	 */
	const t_mat& GetPerspective() const
	{
		return m_matPerspective;
	}


	/**
	 * get the camera's inverse transformation matrix
	 */
	const t_mat& GetInversePerspective() const
	{
		return m_matPerspective_inv;
	}


	/**
	 * get the camera's full transformation matrix
	 */
	const t_mat& GetViewport() const
	{
		return m_matViewport;
	}


	/**
	 * get the camera's inverse transformation matrix
	 */
	const t_mat& GetInverseViewport() const
	{
		return m_matViewport_inv;
	}


	/**
	 * sets perspective or parallel projection
	 */
	void SetPerspectiveProjection(bool proj)
	{
		m_persp_proj = proj;
		m_persp_needs_update = true;
	}


	/**
	 * is the perspective projection active
	 */
	bool GetPerspectiveProjection() const
	{
		return m_persp_proj;
	}


	/**
	 * sets scree aspect ratio, height / width
	 */
	void SetAspectRatio(t_real aspect)
	{
		m_aspect = aspect;
		m_persp_needs_update = true;
	}


	/**
	 * set the screen width and height
	 */
	void SetScreenDimensions(int w, int h)
	{
		m_screenDims[0] = w;
		m_screenDims[1] = h;

		m_viewport_needs_update = true;

		SetAspectRatio(t_real(h) / t_real(w));
	}


	/**
	 * get the screen width and height
	 */
	const std::array<int, 2>& GetScreenDimensions() const
	{
		return m_screenDims;
	}


	/**
	 * get the z buffer depth range
	 */
	std::tuple<t_real, t_real> GetDepthRange() const
	{
		return std::make_tuple(m_z_near, m_z_far);
	}


	/**
	 * is the transformation matrix outdated?
	 */
	bool TransformationNeedsUpdate() const
	{
		return m_trafo_needs_update;
	}


	/**
	 * is the perspective matrix outdated?
	 */
	bool PerspectiveNeedsUpdate() const
	{
		return m_persp_needs_update;
	}


	/**
	 * is the viewport matrix outdated?
	 */
	bool ViewportNeedsUpdate() const
	{
		return m_viewport_needs_update;
	}


	/**
	 * convert a vector into screen coordinates
	 */
	t_vec ToScreenCoords(
		const t_vec& vec4, bool *visible = nullptr) const
	{
		auto [ persp, vec ] =
			tl2::hom_to_screen_coords<t_mat, t_vec>(
				vec4, GetTransformation(), GetPerspective(),
				GetViewport(), true);

		// position not visible -> return a point outside the viewport
		if(persp[2] > 1.)
		{
			if(visible)
				*visible = false;

			return tl2::create<t_vec>(
			{
				t_real(-1.) * GetScreenDimensions()[0],
				t_real(-1.) * GetScreenDimensions()[1]
			});
		}

		if(visible)
			*visible = true;

		return vec;
	}


	/**
	 * get position of object relative to the camera frustum
	 * -1: left of camera frustum
	 * +1: right of camera frustum
	 * -2: below camera frustum
	 * +2: above of camera frustum
	 * -3: in front of near plane
	 * +3: beyond far plane
	 */
	std::unordered_set<int> GetFrustumSides(const t_vec& vec) const
	{
		// projected vector
		t_vec vec_trafo = m_matPerspective * m_mat * vec;
		vec_trafo /= vec_trafo[3];

		std::unordered_set<int> sides;

		if(vec_trafo[0] < t_real(-1.))
			sides.insert(-1);
		else if(vec_trafo[0] > t_real(1.))
			sides.insert(1);

		if(vec_trafo[1] < t_real(-1.))
			sides.insert(-2);
		else if(vec_trafo[1] > t_real(1.))
			sides.insert(2);

		if(vec_trafo[2] < t_real(-1.))
			sides.insert(-3);
		else if(vec_trafo[2] > t_real(1.))
			sides.insert(3);

		return sides;
	}


	/**
	 * test if bounding box is outside frustum
	 */
	bool IsBoundingBoxOutsideFrustum(const t_mat& matObj,
		const std::vector<t_vec>& bbox) const
	{
		bool first_vec = true;
		std::unordered_set<int> set_inters;

		for(const t_vec& _vec : bbox)
		{
			t_vec vec = matObj * _vec;
			std::unordered_set<int> sides = GetFrustumSides(vec);

			// inside the frustum?
			if(sides.size() == 0)
				return false;

			if(first_vec)
			{
				set_inters = sides;
				first_vec = false;
			}
			else
			{
				std::unordered_set<int> new_inters;
				std::set_intersection(
					sides.begin(), sides.end(),
					set_inters.begin(), set_inters.end(),
					std::inserter(new_inters, new_inters.begin()));
				set_inters = std::move(new_inters);

				if(set_inters.size() == 0)
					return false;
			}
		}

		// all outside the same frustum plane?
		return set_inters.size() != 0;
	}


	/**
	 * get a projected bounding rectangle from an object's bounding box
	 * (note that this is larger than the bounding rectangle of the actual object vertices)
	 */
	std::vector<t_vec> GetBoundingRect(const t_mat& matObj,
		const std::vector<t_vec>& bbox) const
	{
		std::vector<t_vec> brect;
		brect.reserve(bbox.size());

		for(const t_vec& _vec : bbox)
		{
			// projected vector
			t_vec vec = m_matPerspective * m_mat * matObj * _vec;
			vec /= vec[3];

			brect.emplace_back(std::move(vec));
		}

		auto [min, max] = tl2::bounding_box(brect);

		return std::vector<t_vec>{{
			tl2::create<t_vec>({min[0], min[1], 0., 1.}),
			tl2::create<t_vec>({min[0], max[1], 0., 1.}),
			tl2::create<t_vec>({max[0], max[1], 0., 1.}),
			tl2::create<t_vec>({max[0], min[1], 0., 1.}),
		}};
	}


	/**
	 * get a ray from screen coordinates
	 */
	std::tuple<t_vec3, t_vec3> GetPickerRay(t_real mouseX, t_real mouseY) const
	{
		// picker ray
		auto [org, dir] = tl2::hom_line_from_screen_coords<t_mat, t_vec>(
			mouseX, mouseY, 0., 1.,
			GetInverseTransformation(),
			GetInversePerspective(),
			GetInverseViewport(),
			&GetViewport(), true);

		t_vec3 org3 = tl2::create<t_vec3>({org[0], org[1], org[2]});
		t_vec3 dir3 = tl2::create<t_vec3>({dir[0], dir[1], dir[2]});

		return std::make_tuple(std::move(org3), std::move(dir3));
	}


//protected:
	/**
	 * update camera matrices
	 */
	void UpdateTransformation()
	{
		const t_vec3 vecCamDir[2] =
		{
			tl2::create<t_vec3>({ 1., 0., 0. }),
			tl2::create<t_vec3>({ 0. ,0., 1. })
		};

		m_matRot = tl2::hom_rotation<t_mat, t_vec3>(
			vecCamDir[0], m_theta, 0);
		m_matRot *= tl2::hom_rotation<t_mat, t_vec3>(
			vecCamDir[1], m_phi, 0);

		t_real zoom = GetPerspectiveProjection() ? m_zoom : 0.5;
		t_mat matDist = hom_translation<t_mat, t_real>(0., 0., -m_dist / zoom);
		m_mat = matDist * m_matRot * m_matTrans;
		std::tie(m_mat_inv, std::ignore) = tl2::inv<t_mat>(m_mat);

		m_trafo_needs_update = false;
	}


	/**
	 * update camera perspective matrices
	 */
	void UpdatePerspective()
	{
		// projection
		if(m_persp_proj)
		{
			m_matPerspective = tl2::hom_perspective<t_mat, t_real>(
				m_nearPlane, m_farPlane, m_FOV, m_aspect);
		}
		else
		{
			m_matPerspective = tl2::hom_ortho_sym<t_mat, t_real>(
				m_nearPlane, m_farPlane,
				m_parallel_range, m_parallel_range * m_aspect);
		}

		std::tie(m_matPerspective_inv, std::ignore) =
			tl2::inv<t_mat>(m_matPerspective);

		m_persp_needs_update = false;
	}


	/**
	 * update camera perspective matrices
	 */
	void UpdateViewport()
	{
		// viewport
		m_matViewport = tl2::hom_viewport<t_mat, t_real>(
			m_screenDims[0], m_screenDims[1],
			m_z_near, m_z_far);
		std::tie(m_matViewport_inv, std::ignore) =
			tl2::inv<t_mat>(m_matViewport);

		m_viewport_needs_update = false;
	}


private:
	// full transformation matrix and its inverse
	t_mat m_mat = tl2::unit<t_mat>();
	t_mat m_mat_inv = tl2::unit<t_mat>();

	// rotation and translation matrices
	t_mat m_matRot = tl2::unit<t_mat>();
	t_mat m_matTrans = tl2::unit<t_mat>();

	// field of view
	t_real m_FOV = tl2::pi<t_real>*t_real(0.5);

	// camera frustum near and far planes
	t_real m_nearPlane = 0.1;
	t_real m_farPlane = 1000.;

	// range for parallel projection
	t_real m_parallel_range = 20.;

	// camera rotation
	t_real m_phi = pi<t_real>*t_real(0.25);
	t_real m_theta = -pi<t_real>*(0.25);
	t_real m_phi_saved = pi<t_real>*t_real(0.25);
	t_real m_theta_saved = -pi<t_real>*(0.25);

	// camera zoom (giving the fraction of m_dist towards the translation centre)
	t_real m_zoom = 1.;

	// distance from camera centre
	t_real m_dist = 15.;

	// perspective matrix and its inverse
	t_mat m_matPerspective = tl2::unit<t_mat>();
	t_mat m_matPerspective_inv = tl2::unit<t_mat>();

	// perspective or parallel projection?
	bool m_persp_proj = true;

	// screen aspect ratio
	t_real m_aspect = 1.;

	// screen viewport
	t_mat m_matViewport = tl2::unit<t_mat>();
	t_mat m_matViewport_inv = tl2::unit<t_mat>();

	// z buffer range
	t_real m_z_near{0}, m_z_far{1};

	// screen dimensions
	std::array<int, 2> m_screenDims = { 800, 600 };

	// does the transformation matrix need an update?
	bool m_trafo_needs_update = true;

	// does the perspective matrix need an update?
	bool m_persp_needs_update = true;

	// does the perspective matrix need an update?
	bool m_viewport_needs_update = true;
};

}

#endif
