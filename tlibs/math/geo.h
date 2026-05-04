/**
 * linalg and geometry helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 30-apr-2013 - 2018
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
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

#ifndef __TLIBS_GEO_H__
#define __TLIBS_GEO_H__

#include "../helper/flags.h"
#include "../helper/exception.h"
#include "linalg.h"
#include "linalg2.h"
#include "stat.h"
#include "../log/log.h"

#include <iostream>
#include <cmath>
#include <utility>


namespace tl {

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;


//------------------------------------------------------------------------------

template<typename T> class Line;


/**
 * analytical geometry of a plane
 * @see (Stoecker 1999), chapter "Analytische Geometrie".
 */
template<typename T> class Plane
{
public:
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;


protected:
	bool m_bValid = false;
	t_vec m_vecX0;
	t_vec m_vecDir0, m_vecDir1;
	t_vec m_vecNorm;
	T m_d;


public:
	/**
	 * plane from a point and a normal
	 */
	Plane(const t_vec& vec0, const t_vec& vecNorm)
		: m_vecX0(vec0), m_vecNorm(vecNorm)
	{
		// normalise normal
		T tLenNorm = veclen(m_vecNorm);
		if(float_equal<T>(tLenNorm, 0.) || tLenNorm!=tLenNorm)
		{
			m_bValid = false;
			return;
		}
		m_vecNorm /= tLenNorm;

		// Hessian form: vecX0*vecNorm - d = 0
		m_d = inner(m_vecX0, m_vecNorm);


		// find direction vectors
		std::vector<t_vec> vecTry =
			{ tl::make_vec({T(1), T(0), T(0)}),
			tl::make_vec({T(0), T(1), T(0)}),
			tl::make_vec({T(0), T(0), T(1)})};

		std::size_t iIdxBest = 0;
		T dDot = T(1);
		for(std::size_t iIdx=0; iIdx<vecTry.size(); ++iIdx)
		{
			const t_vec& vec = vecTry[iIdx];

			T dDotCur = std::abs(inner(vec, m_vecNorm));
			if(dDotCur < dDot)
			{
				iIdxBest = iIdx;
				dDot = dDotCur;
			}
		}
		m_vecDir0 = vecTry[iIdxBest];
		m_vecDir1 = cross_3(m_vecNorm, m_vecDir0);
		m_vecDir0 = cross_3(m_vecDir1, m_vecNorm);

		m_bValid = true;
	}


	/**
	 * plane from a point and two directions on the plane
	 */
	Plane(const t_vec& vec0,
		const t_vec& dir0, const t_vec& dir1)
		: m_vecX0(vec0), m_vecDir0(dir0), m_vecDir1(dir1)
	{
		// calculate normal
		m_vecNorm = cross_3(dir0, dir1);

		// normalise normal
		T tLenNorm = veclen(m_vecNorm);
		if(float_equal<T>(tLenNorm, 0.) || tLenNorm!=tLenNorm)
		{
			m_bValid = false;
			return;
		}
		m_vecNorm /= tLenNorm;

		// Hessian form: vecX0*vecNorm - d = 0
		m_d = inner(m_vecX0, m_vecNorm);
		m_bValid = true;
	}


	Plane() = default;
	~Plane() = default;


	const t_vec& GetX0() const { return m_vecX0; }
	const t_vec& GetDir0() const { return m_vecDir0; }
	const t_vec& GetDir1() const { return m_vecDir1; }
	const t_vec& GetNorm() const { return m_vecNorm; }
	const T& GetD() const { return m_d; }


	T GetDist(const t_vec& vecPt) const
	{
		return inner(vecPt, m_vecNorm) - m_d;
	}


	T GetAngle(const Plane<T>& plane) const
	{
		T dot = inner(GetNorm(), plane.GetNorm());
		return std::acos(dot);
	}


	T GetAngle(const t_vec& _vec) const
	{
		t_vec vec = _vec / veclen(_vec);
		T dot = inner(GetNorm(), vec);
		return std::asin(dot);
	}


	void FlipNormal()
	{
		m_vecNorm = -m_vecNorm;
		m_d = -m_d;

		std::swap(m_vecDir0, m_vecDir1);
	}


	/**
	 * "Lotfußpunkt"
	 * @see e.g. https://de.wikipedia.org/wiki/Lot_(Mathematik)
	 */
	t_vec GetDroppedPerp(const t_vec& vecP, T *pdDist=0) const
	{
		T dist = GetDist(vecP);
		t_vec vecdropped = vecP - dist*m_vecNorm;

		if(pdDist)
		{
			t_vec vecD = vecP - vecdropped;
			*pdDist = std::sqrt(inner(vecD, vecD));
		}

		return vecdropped;
	}


	/**
	 * determine on which side of the plane a point is located
	 */
	bool GetSide(const t_vec& vecP, T *pdDist=0) const
	{
		T dDist = GetDist(vecP);
		if(pdDist) *pdDist = dDist;
		return dDist < T(0);
	}


	/**
	 * determine if a point is on the plane
	 */
	bool IsOnPlane(const t_vec& vecPt, T eps = tl::get_epsilon<T>()) const
	{
		T dDist = GetDist(vecPt);
		return float_equal(dDist, T(0), eps);
	}


	bool IsParallel(const Plane<T>& plane, T eps = tl::get_epsilon<T>()) const
	{
		return vec_is_collinear<t_vec>(GetNorm(), plane.GetNorm(), eps);
	}


	/**
	 * plane-plane intersection
	 * @see http://mathworld.wolfram.com/Plane-PlaneIntersection.html
	 */
	bool intersect(const Plane<T>& plane2, Line<T>& lineRet,
		T eps = tl::get_epsilon<T>()) const
	{
		if(IsParallel(plane2, eps))
			return false;

		const Plane<T>& plane1 = *this;

		// direction vector
		t_vec vecDir = cross_3(plane1.GetNorm(), plane2.GetNorm());

		// find common point in the two planes
		t_mat M = row_matrix( { plane1.GetNorm(), plane2.GetNorm() } );

		t_vec vecD(2);
		vecD[0] = plane1.GetD();
		vecD[1] = plane2.GetD();

		t_vec vec0(3);
		if(!tl::solve_linear(M, vecD, vec0))
			return 0;

		lineRet = Line<T>(vec0, vecDir);
		return true;
	}


	/**
	 * intersection point of three planes
	 * @see http://mathworld.wolfram.com/Plane-PlaneIntersection.html
	 */
	bool intersect(const Plane<T>& plane2, const Plane<T>& plane3, t_vec& ptRet,
		T eps = tl::get_epsilon<T>()) const
	{
		const Plane<T>& plane1 = *this;

		// det
		const t_mat matNorms = row_matrix( { plane1.GetNorm(), plane2.GetNorm(), plane3.GetNorm() } );
		const T detM = determinant(matNorms);
		if(float_equal(detM, T(0), eps))
			return false;

		// direction vectors
		const t_vec vecDir12 = cross_3(plane1.GetNorm(), plane2.GetNorm());
		const t_vec vecDir23 = cross_3(plane2.GetNorm(), plane3.GetNorm());
		const t_vec vecDir31 = cross_3(plane3.GetNorm(), plane1.GetNorm());

		ptRet = (plane1.GetD()*vecDir23 + plane2.GetD()*vecDir31 + plane3.GetD()*vecDir12) / detM;
		if(is_nan_or_inf(ptRet))
			return false;
		return true;
	}


	bool IsValid() const { return m_bValid; }
};



//------------------------------------------------------------------------------



/**
 * analytical geometry of a line
 * @see (Stoecker 1999), chapter "Analytische Geometrie".
 */
template<typename T> class Line
{
public:
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;


protected:
	t_vec m_vecX0;
	t_vec m_vecDir;


public:
	Line() {}


	Line(const t_vec& vec0, const t_vec& dir)
		: m_vecX0(vec0), m_vecDir(dir)
	{}


	~Line() = default;


	t_vec operator()(T t) const
	{
		return m_vecX0 + t*m_vecDir;
	}


	const t_vec& GetX0() const { return m_vecX0; }
	const t_vec& GetDir() const { return m_vecDir; }


	/**
	 * distance to a point
	 * @see e.g.: (Arens 2015), p. 711
	 */
	T GetDist(const t_vec& vecPt) const
	{
		const t_vec& vecX0 = GetX0();
		t_vec vecDir = GetDir() / veclen(GetDir());

		// shift everything so that line goes through the origin
		t_vec vecPtShift = vecPt - vecX0;

		// project point on direction vector
		t_vec vecClosestPt = inner(vecPtShift, vecDir) * vecDir;

		// distance between point and projected point
		return veclen(vecClosestPt - vecPtShift);
	}


	/**
	 * distance to line l1
	 * @see e.g.: (Arens 2015), p. 711
	 */
	T GetDist(const Line<T>& l1) const
	{
		const Line<T>& l0 = *this;

		// vector normal to both directions defining the distance line
		t_vec vecNorm = cross_3<t_vec>(l0.GetDir(), l1.GetDir());
		T tlenNorm = veclen(vecNorm);

		t_vec vec01 = l1.GetX0() - l0.GetX0();

		// if the lines are parallel, any point (e.g. the x0s) can be used
		if(float_equal(tlenNorm, T(0)))
			return GetDist(l1.GetX0());

		// project x0_1 - x0_0 onto vecNorm
		T tdot = std::abs(inner(vec01, vecNorm));
		return tdot / tlenNorm;
	}


	bool IsParallel(const Line<T>& line, T eps = tl::get_epsilon<T>()) const
	{
		return vec_is_collinear<t_vec>(GetDir(), line.GetDir(), eps);
	}


	T GetAngle(const Line<T>& line) const
	{
		t_vec dir1 = GetDir();
		t_vec dir2 = line.GetDir();

		dir1 /= veclen(dir1);
		dir2 /= veclen(dir2);

		T dot = inner(dir1, dir2);
		return std::acos(dot);
	}


	/**
	 * "Lotfußpunkt"
	 * @see e.g. https://de.wikipedia.org/wiki/Lot_(Mathematik)
	 */
	t_vec GetDroppedPerp(const t_vec& vecP, T *pdDist=0) const
	{
		const t_vec& vecDir = GetDir();

		// projection of vecP-x0 onto vecDir
		T t = inner<t_vec>(vecP-GetX0(), vecDir) / inner(vecDir, vecDir);

		t_vec vecdropped = operator()(t);

		if(pdDist)
		{
			t_vec vecD = vecP - vecdropped;
			*pdDist = std::sqrt(inner(vecD, vecD));
		}

		return vecdropped;
	}


	/**
	 * determine on which side of the line a point is located
	 */
	bool GetSide(const t_vec& vecP, T *pdDist=0) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 2)
		{
			log_err("\"Side of line\" only defined for 2d vectors.");
			return false;
		}

		t_vec vecDropped = GetDroppedPerp(vecP, pdDist);


		t_vec vecNorm(2);
		vecNorm[0] = m_vecDir[1];
		vecNorm[1] = -m_vecDir[0];

		T tDot = inner<t_vec>(vecP-vecDropped, vecNorm);
		return tDot < T(0);
	}


	/**
	 * line-plane intersection
	 * @see http://mathworld.wolfram.com/Line-PlaneIntersection.html
	 */
	bool intersect(const Plane<T>& plane, T& t, T eps = tl::get_epsilon<T>()) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 3)
		{
			log_err("Line-plane intersection only implemented for 3d vectors.");
			return false;
		}

		const t_vec& posl = this->GetX0();
		const t_vec& dirl = this->GetDir();

		const t_vec& xp0 = plane.GetX0();
		const t_vec xp1 = plane.GetX0() + plane.GetDir0();
		const t_vec xp2 = plane.GetX0() + plane.GetDir1();

		t_mat matDenom(N+1,N+1);
		matDenom(0,0) = 1;	matDenom(0,1) = 1;	matDenom(0,2) = 1;	matDenom(0,3) = 0;
		matDenom(1,0) = xp0[0];	matDenom(1,1) = xp1[0];	matDenom(1,2) = xp2[0];	matDenom(1,3) = dirl[0];
		matDenom(2,0) = xp0[1];	matDenom(2,1) = xp1[1];	matDenom(2,2) = xp2[1];	matDenom(2,3) = dirl[1];
		matDenom(3,0) = xp0[2];	matDenom(3,1) = xp1[2];	matDenom(3,2) = xp2[2];	matDenom(3,3) = dirl[2];

		T denom = determinant(matDenom);
		if(tl::float_equal<T>(denom, 0., eps))
			return false;

		t_mat matNum(N+1,N+1);
		matNum(0,0) = 1;	matNum(0,1) = 1;	matNum(0,2) = 1;	matNum(0,3) = 1;
		matNum(1,0) = xp0[0];	matNum(1,1) = xp1[0];	matNum(1,2) = xp2[0];	matNum(1,3) = posl[0];
		matNum(2,0) = xp0[1];	matNum(2,1) = xp1[1];	matNum(2,2) = xp2[1];	matNum(2,3) = posl[1];
		matNum(3,0) = xp0[2];	matNum(3,1) = xp1[2];	matNum(3,2) = xp2[2];	matNum(3,3) = posl[2];

		T num = determinant(matNum);

		t = -num / denom;
		return true;
	}


	/**
	 * line-line intersection
	 * @see e.g.: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
	 *
	 * pos0 + t0*dir0 = pos1 + t1*dir1
	 * pos0 - pos1 = t1*dir1 - t0*dir0
	 * pos0 - pos1 = (-dir0 dir1) * (t0 t1)^t
	 * exact solution for the t params: (-dir0 dir1)^(-1) * (pos0 - pos1) = (t0 t1)^t
	 *
	 * generally:
	 * exact: b = M t  ->  M^(-1)*b = t
	 * approx.: M^t b = M^t M t  ->  (M^t M)^(-1) * M^t b = t
	 */
	bool intersect(const Line<T>& line1, T& t0, T eps = tl::get_epsilon<T>(), T *pt1=nullptr) const
	{
		if(IsParallel(line1, eps))
		{
			t0 = 0;
			if(pt1) *pt1 = 0;
			return false;
		}

		const t_vec& pos0 =  this->GetX0();
		const t_vec& pos1 =  line1.GetX0();

		const t_vec& dir0 =  this->GetDir();
		const t_vec& dir1 =  line1.GetDir();

		const t_vec pos = pos0-pos1;
		t_mat M = column_matrix({-dir0, dir1});
		t_mat Mt = transpose(M);
		t_mat MtM = prod_mm(Mt, M);

		t_mat MtMinv;
		if(!tl::inverse(MtM, MtMinv))
			return false;

		t_vec Mtb = prod_mv(Mt, pos);
		t_vec params = prod_mv(MtMinv, Mtb);

		// get parameters
		t0 = params[0];
		T t1 = params[1];
		if(pt1) *pt1 = t1;

		// get closest points between the two lines
		t_vec vecInters0 = (*this)(t0);
		t_vec vecInters1 = line1(t1);

		return veclen(vecInters1-vecInters0) <= eps;
	}


	/**
	 * middle perpendicular line (in 2d)
	 */
	bool GetMiddlePerp(Line<T>& linePerp) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 2)
		{
			log_err("Perpendicular line only implemented for 2d vectors.");
			return false;
		}

		t_vec vecDir(2);
		vecDir[0] = -m_vecDir[1];
		vecDir[1] = m_vecDir[0];

		t_vec vecPos = this->operator()(0.5);

		linePerp = Line<T>(vecPos, vecDir);
		return true;
	}


	/**
	 * middle perpendicular plane (in 3d)
	 */
	bool GetMiddlePerp(Plane<T>& planePerp) const
	{
		const std::size_t N = m_vecDir.size();
		if(N != 3)
		{
			log_err("Perpendicular plane only implemented for 3d vectors.");
			return false;
		}

		t_vec vecPos = this->operator()(0.5);
		planePerp = Plane<T>(vecPos, m_vecDir);
		return true;
	}
};


template<typename T>
std::ostream& operator<<(std::ostream& ostr, const Plane<T>& plane)
{
	ostr << plane.GetX0() << " + s*" << plane.GetDir0()
		<< " + t*" << plane.GetDir1();
	return ostr;
}


template<typename T>
std::ostream& operator<<(std::ostream& ostr, const Line<T>& line)
{
	ostr << line.GetX0() << " + t*" << line.GetDir();
	return ostr;
}



/**
 * intersection of "plane" and polygon (defined by "planePoly" and vertPoly")
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
bool intersect_plane_poly(const Plane<T>& plane,
	const Plane<T>& planePoly, const t_cont<t_vec>& vertPoly,
	Line<T>& lineRes, T eps = tl::get_epsilon<T>())
{
	if(!vertPoly.size())
		return false;

	bool bFirstSide = plane.GetSide(vertPoly[0]);

	// are all vertices on the same side?
	if(std::all_of(vertPoly.begin(), vertPoly.end(),
		[&plane, bFirstSide](const t_vec& vert) -> bool
		{ return plane.GetSide(vert) == bFirstSide; }))
		return false;


	if(!plane.intersect(planePoly, lineRes, eps))
		return false;

	return true;
}


/**
 * intersection of "line" and polygon (defined by "planePoly" and vertPoly")
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
bool intersect_line_poly(const Line<T>& line,
	const Plane<T>& planePoly, const t_cont<t_vec>& vertPoly,
	t_vec& vecIntersect, T eps = tl::get_epsilon<T>())
{
	// point of intersection with plane
	T t;
	if(!line.intersect(planePoly, t, eps))
		return false;
	vecIntersect = line(t);

	// is intersection point within polygon?
	const t_vec vecFaceCentre = mean_value(vertPoly);

	for(std::size_t iVert = 0; iVert < vertPoly.size(); ++iVert)
	{
		std::size_t iNextVert = iVert < (vertPoly.size()-1) ? iVert+1 : 0;

		const t_vec vecEdge = vertPoly[iNextVert] - vertPoly[iVert];
		const t_vec vecEdgeCentre = vertPoly[iVert] + T(0.5)*vecEdge;
		const t_vec vecOut = vecFaceCentre - vecEdgeCentre;
		t_vec vecNorm = cross_3(vecEdge, planePoly.GetNorm());
		if(inner(vecNorm, vecOut) < T(0))
			vecNorm = -vecNorm;

		const Plane<T> planeEdge(vecEdgeCentre, vecNorm);

		if(planeEdge.GetDist(vecIntersect) < -eps)
			return false;
	}

	return true;
}



/**
 * sort vertices in a convex polygon
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
void sort_poly_verts_norm(t_cont<t_vec>& vecPoly, const t_vec& _vecNorm)
{
	if(vecPoly.size() <= 1)
		return;

	// line from centre to vertex
	const t_vec vecCentre = mean_value(vecPoly);
	const t_vec vecNorm = _vecNorm / veclen(_vecNorm);

	t_vec vec0 = vecPoly[0] - vecCentre;

	std::stable_sort(vecPoly.begin(), vecPoly.end(),
		[&vecCentre, &vec0, &vecNorm](const t_vec& vertex1, const t_vec& vertex2) -> bool
		{
			t_vec vec1 = vertex1 - vecCentre;
			t_vec vec2 = vertex2 - vecCentre;

			return vec_angle(vec0, vec1, &vecNorm) < vec_angle(vec0, vec2, &vecNorm);
		});
}


/**
 * sort vertices in a convex polygon using an absolute centre for determining the normal
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
void sort_poly_verts(t_cont<t_vec>& vecPoly, const t_vec& vecAbsCentre)
{
	if(vecPoly.size() <= 1)
		return;

	// line from centre to vertex
	const t_vec vecCentre = mean_value(vecPoly);
	// face normal
	t_vec vecNorm = vecCentre - vecAbsCentre;

	sort_poly_verts_norm<t_vec, t_cont, T>(vecPoly, vecNorm);
}


/**
 * sort vertices in a convex polygon determining normal
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
void sort_poly_verts(t_cont<t_vec>& vecPoly)
{
	if(vecPoly.size() <= 1)
		return;

	// line from centre to vertex
	const t_vec vecCentre = mean_value(vecPoly);
	// face normal
	t_vec vecNormBest;
	T tBestCross = T(0);

	// find non-collinear vectors
	for(std::size_t iVecPoly=1; iVecPoly<vecPoly.size(); ++iVecPoly)
	{
		t_vec vecNorm = cross_3<t_vec>(vecPoly[0]-vecCentre, vecPoly[1]-vecCentre);
		T tCross = veclen(vecNorm);
		if(tCross > tBestCross)
		{
			tBestCross = tCross;
			vecNormBest = vecNorm;
		}
	}

	// nothing found
	if(vecNormBest.size() < vecCentre.size())
		return;

	sort_poly_verts_norm<t_vec, t_cont, T>(vecPoly, vecNormBest);
}


/**
 * get the polygon's face normal vector
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
t_vec get_face_normal(const t_cont<t_vec>& vecVerts, t_vec vecCentre,
	T eps = tl::get_epsilon<T>())
{
	t_vec vecNorm;

	if(vecVerts.size() < 3)
		return vecNorm;

	const t_vec& vec1 = vecVerts[1] - vecVerts[0];
	for(std::size_t iVec2 = 2; iVec2 < vecVerts.size(); ++iVec2)
	{
		const t_vec& vec2 = vecVerts[iVec2] - vecVerts[0];
		t_vec vecCross = tl::cross_3(vec1, vec2);
		if(veclen(vecCross) > eps)
		{
			vecNorm = vecCross;
			break;
		}
	}

	// nothing found
	if(vecNorm.size() == 0)
		return vecNorm;


	// vector pointing "outwards"
	const t_vec vecPolyCentre = mean_value(vecVerts);
	t_vec vecCentreToFace = vecPolyCentre - vecCentre;
	vecCentreToFace /= veclen(vecCentreToFace);

	vecNorm /= veclen(vecNorm);
	if(inner(vecNorm, vecCentreToFace) < T(0))
		vecNorm = -vecNorm;

	return vecNorm;
}


/**
 * creates polyhedron faces out of vertices
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_cont = std::vector,
	class T = typename t_vec::value_type>
t_cont<t_cont<t_vec>> verts_to_polyhedron(
	const t_cont<t_vec>& vecVerts, T eps = tl::get_epsilon<T>())
{
	t_cont<t_cont<t_vec>> vecPolys;
	t_cont<Plane<T>> vecPlanes;

	// too few vertices
	if(vecVerts.size() < 3)
		return vecPolys;


	auto fkt_try_add_new_pt = [eps](t_cont<t_vec>& vecPoly, const t_vec& vecPt) -> void
	{
		for(const t_vec& vecOld : vecPoly)
		{
			// vecPt already in set?
			if(vec_equal(vecOld, vecPt, T(1)*eps))
				return;
		}

		// add new point
		vecPoly.push_back(vecPt);
	};


	const t_vec vecCentre = mean_value(vecVerts);


	for(std::size_t i=0; i<vecVerts.size(); ++i)
	{
		for(std::size_t j=i+1; j<vecVerts.size(); ++j)
		{
			//if(i==j) continue;
			for(std::size_t k=j+1; k<vecVerts.size(); ++k)
			{
				//if(i==k || j==k) continue;

				const t_vec& vert0 = vecVerts[i];
				const t_vec& vert1 = vecVerts[j];
				const t_vec& vert2 = vecVerts[k];

				bool bPointsUsed = false;
				for(std::size_t iPlane=0; iPlane<vecPlanes.size(); ++iPlane)
				{
					const Plane<T>& plane = vecPlanes[iPlane];
					t_cont<t_vec>& vecPoly = vecPolys[iPlane];

					// plane already in set?
					if(plane.IsOnPlane(vert0, T(10)*eps) &&
						plane.IsOnPlane(vert1, T(10)*eps) &&
						plane.IsOnPlane(vert2, T(10)*eps))
					{
						// add new points to poly if not already present
						fkt_try_add_new_pt(vecPoly, vert0);
						fkt_try_add_new_pt(vecPoly, vert1);
						fkt_try_add_new_pt(vecPoly, vert2);

						bPointsUsed = true;
						break;
					}
				}


				// new plane / poly
				if(!bPointsUsed)
				{
					t_cont<t_vec> vecPoly = { vert0, vert1, vert2 };
					sort_poly_verts<t_vec, t_cont, T>(vecPoly, vecCentre);
					Plane<T> plane(vecPoly[0], vecPoly[1]-vecPoly[0], vecPoly[2]-vecPoly[0]);

					// check if all other vertices are on the negative side of the polygon
					bool bIsHull = true;
					for(const t_vec& vecOtherVert : vecVerts)
					{
						if(plane.GetDist(vecOtherVert) > T(10)*eps)
						{
							bIsHull = false;
							break;
						}
					}

					if(bIsHull)
					{
						vecPolys.emplace_back(std::move(vecPoly));
						vecPlanes.emplace_back(std::move(plane));
					}
				}
			}
		}
	}

	// too few polygons => remove polyhedron
	if(vecPolys.size() < 3)
		vecPolys = decltype(vecPolys)();

	// sort vertices
	for(auto& vecPoly : vecPolys)
		sort_poly_verts<t_vec, t_cont, T>(vecPoly, vecCentre);

	return vecPolys;
}


//------------------------------------------------------------------------------


/**
 * quadric
 * @see e.g.: (Arens 2015), ch. 21
 * @see e.g.: (Merziger 1993), p. 224
 */
template<class T = double>
class Quadric
{
public:
	using t_vec = ublas::vector<T>;
	using t_mat = ublas::matrix<T>;


protected:
	// x^T Q x  +  r x  +  s  =  0
	t_mat m_Q = zero_m<t_mat>(3,3);
	t_vec m_r = zero_v<t_vec>(3);
	T m_s = 0;

	t_vec m_vecOffs = zero_v<t_vec>(3);
	bool m_bQSymm = true;


protected:
	void CheckSymm()
	{
		m_bQSymm = is_symmetric(m_Q, std::cbrt(get_epsilon<T>()));
	}


public:
	Quadric() {}


	Quadric(std::size_t iDim)
		: m_Q(zero_m<t_mat>(iDim,iDim)), m_r(zero_v<t_vec>(iDim))
	{ CheckSymm(); }


	Quadric(const t_mat& Q) : m_Q(Q)
	{ CheckSymm(); }


	Quadric(const t_mat& Q, const t_vec& r, T s)
		: m_Q(Q), m_r(r), m_s(s)
	{ CheckSymm(); }
	~Quadric() {}


	void SetDim(std::size_t iDim) { m_Q.resize(iDim, iDim, 1); }


	/**
	 * get the classification of the quadric
	 * @see https://mathworld.wolfram.com/QuadraticSurface.html
	 * @returns [rank, rank_ext, signature, signature_ext]
	 */
	std::tuple<int,int, int,int,int, int,int,int>
	ClassifyQuadric(T eps = tl::get_epsilon<T>()) const
	{
		// extended matrix
		t_mat Qext = m_Q;
		Qext.resize(m_Q.size1()+1, m_Q.size2()+1, true);
		Qext(Qext.size1()-1, Qext.size2()-1) = m_s;
		for(std::size_t i=0; i<m_r.size(); ++i)
			Qext(i, Qext.size2()-1) = Qext(Qext.size1()-1, i) = m_r[i];


		// ------------------------------------------------------------
		// get eigenvalues
		std::vector<t_vec> evecs;
		std::vector<T> evals;
		bool bEV = false;
		if(m_bQSymm)
			bEV = eigenvec_sym(m_Q, evecs, evals);
		else
			bEV = eigenvec_approxsym(m_Q, evecs, evals);

		if(!bEV)
			return std::make_tuple(-1,-1, -1,-1,-1, -1,-1,-1);

		std::vector<t_vec> evecsext;
		std::vector<T> evalsext;
		if(m_bQSymm)
			bEV = eigenvec_sym(Qext, evecsext, evalsext);
		else
			bEV = eigenvec_approxsym(Qext, evecsext, evalsext);

		if(!bEV)
			return std::make_tuple(-1,-1, -1,-1,-1, -1,-1,-1);
		// ------------------------------------------------------------


		// ------------------------------------------------------------
		// get ranks

		auto get_rank = [eps](const std::vector<T>& evals) -> std::tuple<int, int,int,int>
		{
			int pos_evals = 0;
			int neg_evals = 0;
			int zero_evals = 0;
			int rank = 0;

			for(T eval : evals)
			{
				if(float_equal<T>(eval, T(0), eps))
				{
					++zero_evals;
				}
				else if(eval > T(0))
				{
					++pos_evals;
					++rank;
				}
				else if(eval < T(0))
				{
					++neg_evals;
					++rank;
				}
			}

			return std::make_tuple(rank, pos_evals, neg_evals, zero_evals);
		};

		int pos_evals = 0, pos_evalsext = 0;
		int neg_evals = 0, neg_evalsext = 0;
		int zero_evals = 0, zero_evalsext = 0;
		int rank = 0, rankext = 0;

		std::tie(rank, pos_evals, neg_evals, zero_evals) = get_rank(evals);
		std::tie(rankext, pos_evalsext, neg_evalsext, zero_evalsext) = get_rank(evalsext);
		// ------------------------------------------------------------


		return std::make_tuple(rank, rankext,
			pos_evals, neg_evals, zero_evals,
			pos_evalsext, neg_evalsext, zero_evalsext);
	}


	const Quadric<T>& operator=(const Quadric<T>& quad)
	{
		this->m_Q = quad.m_Q;
		this->m_r = quad.m_r;
		this->m_s = quad.m_s;
		this->m_vecOffs = quad.m_vecOffs;
		this->m_bQSymm = quad.m_bQSymm;

		return *this;
	}


	Quadric<T>& operator=(Quadric<T>&& quad)
	{
		this->m_Q = std::move(quad.m_Q);
		this->m_r = std::move(quad.m_r);
		this->m_s = std::move(quad.m_s);
		this->m_vecOffs = std::move(quad.m_vecOffs);
		this->m_bQSymm = quad.m_bQSymm;

		return *this;
	}


	Quadric(const Quadric<T>& quad) { *this = quad; }
	Quadric(Quadric<T>&& quad) { *this = quad; }


	void SetOffset(const t_vec& vec) { m_vecOffs = vec; }
	const t_vec& GetOffset() const { return m_vecOffs; }


	const t_mat& GetQ() const { return m_Q; }
	const t_vec& GetR() const { return m_r; }
	T GetS() const { return m_s; }


	void SetQ(const t_mat& Q) { m_Q = Q; CheckSymm(); }
	void SetR(const t_vec& r) { m_r = r; }
	void SetS(T s) { m_s = s; }


	T operator()(const t_vec& _x) const
	{
		t_vec x = _x-m_vecOffs;

		t_vec vecQ = prod_mv(m_Q, x);
		T dQ = inner(x, vecQ);
		T dR = inner(m_r, x);

		return dQ + dR + m_s;
	}


	/**
	 * remove column and row iIdx
	 */
	void RemoveElems(std::size_t iIdx)
	{
		m_Q = remove_elems(m_Q, iIdx);
		m_r = remove_elem(m_r, iIdx);
		m_vecOffs = remove_elem(m_vecOffs, iIdx);
	}


	void transform(const t_mat& S)
	{
		m_Q = tl::transform<t_mat>(m_Q, S, 1);
		CheckSymm();
	}


	/**
	 * Q = O D O^T
	 * O: eigenvecs, D: eigenvals
	 */
	bool GetPrincipalAxes(t_mat& matEvecs, std::vector<T>& vecEvals,
		Quadric<T>* pquadPrincipal = nullptr) const
	{
		std::vector<t_vec> evecs;

		bool bEV = false;
		if(m_bQSymm)
			bEV = eigenvec_sym(m_Q, evecs, vecEvals);
		else
			bEV = eigenvec_approxsym(m_Q, evecs, vecEvals);
		if(!bEV)
		{
			log_err("Could not exactly determine quadric's principal axes, Q = ", m_Q, ".");
			//return false;
		}

		sort_eigenvecs<T>(evecs, vecEvals, 1,
			[](T d) -> T { return 1./std::sqrt(d); });

		if(determinant(matEvecs) < T(0) && evecs.size() >= 2)
		{
			std::swap(evecs[evecs.size()-2], evecs[evecs.size()-1]);
			std::swap(vecEvals[vecEvals.size()-2], vecEvals[vecEvals.size()-1]);
		}

		matEvecs = column_matrix(evecs);

		if(pquadPrincipal)
		{
			t_mat matEvals = diag_matrix(vecEvals);

			pquadPrincipal->SetDim(vecEvals.size());
			pquadPrincipal->SetQ(matEvals);
			pquadPrincipal->SetS(GetS());

			t_mat matEvecsT = transpose(matEvecs);
			pquadPrincipal->SetR(prod_mv(matEvecsT, GetR()));
		}

		return bEV;
	}


	/**
	 * only valid in principal axis system with diagonal Q:
	 * x^T Q x + rx = 0
	 * q11*x1^2 + r1*x1 + ... = 0
	 * q11*(x1^2 + r1/q11*x1) = 0
	 * completing the square: q11*(x1 + r1/(2*q11))^2 - r1^2/(4*q11)
	 */
	t_vec GetPrincipalOffset() const
	{
		/*
		// component version (only valid in principal axis system)
		t_vec vecOffs = GetR();

		for(std::size_t i = 0; i < vecOffs.size(); ++i)
			vecOffs[i] /= -T(2)*GetQ()(i,i);
		*/

		// matrix version
		t_mat Qinv;
		if(!tl::inverse(GetQ(), Qinv))
		{
			tl::log_err("Could not invert quadric matrix.");
			return tl::zero_vector<t_vec>(GetR().size());
		}

		t_vec vecOffs = -0.5*tl::prod_mv(Qinv, GetR());
		return vecOffs;
	}


	/**
	 * here: only for x^T Q x + s  =  0, i.e. for r=0
	 * quad: x^T Q x + s = 0; line: x = x0 + t d
	 * (x0 + t d)^T Q (x0 + t d) + s = 0
	 * (x0 + t d)^T Q x0 + (x0 + t d)^T Q t d + s = 0
	 * (x0^T + t d^T) Q x0 + (x0^T + t d^T) Q t d + s = 0
	 * x0^T Q x0 + s  +  (d^T Q x0 + x0^T Q d) t  +  d^T Q d t^2 = 0
	 */
	std::vector<T> intersect(const Line<T>& line) const
	{
		const t_mat& Q = GetQ();
		const T& s = m_s;
		const t_vec& d = line.GetDir();
		const t_vec x0 = line.GetX0() - m_vecOffs;;

		// solving at^2 + bt + c = 0 for t
		t_vec vecQd = prod_mv(Q, d);
		T a = inner(d, vecQd);

		t_vec vecQx0 = prod_mv(Q, x0);
		T c = inner(x0, vecQx0) + s;

		T b = inner(x0, vecQd);
		b += inner(d, vecQx0);

		return solve_quadratic(a,b,c);
	}
};


template<class T = double>
std::ostream& operator<<(std::ostream& ostr, const Quadric<T>& quad)
{
	ostr << "Q = " << quad.GetQ() << ", ";
	ostr << "s = " << quad.GetS();
	return ostr;
}


template<class T = double>
class QuadSphere : public Quadric<T>
{
protected:

public:
	QuadSphere() : Quadric<T>()
	{
		this->m_s = T(-1);
	}

	QuadSphere(std::size_t iDim) : Quadric<T>(iDim)
	{
		this->m_s = T(-1);
	}

	QuadSphere(T r) : Quadric<T>(3)
	{
		this->m_Q(0,0) =
			this->m_Q(1,1) =
			this->m_Q(2,2) = T(1.)/(r*r);

		this->m_s = T(-1.);
	}

	QuadSphere(std::size_t iDim, T r) : Quadric<T>(iDim)
	{
		for(std::size_t i=0; i<iDim; ++i)
			this->m_Q(i,i) = T(1.)/(r*r);

		this->m_s = T(-1.);
	}

	/**
	 * only valid in principal axis system
	 */
	T GetRadius() const
	{
		return std::abs(this->m_s) /
			std::sqrt(std::abs(this->m_Q(0,0)));
	}

	T GetVolume() const
	{
		return get_ellipsoid_volume(this->m_Q) /
			std::abs(this->m_s);
	}

	~QuadSphere() {}
};


template<class T = double>
class QuadEllipsoid : public Quadric<T>
{
protected:

public:
	QuadEllipsoid() : Quadric<T>()
	{
		this->m_s = T(-1);
	}

	QuadEllipsoid(std::size_t iDim) : Quadric<T>(iDim)
	{
		this->m_s = T(-1);
	}

	QuadEllipsoid(T a, T b) : Quadric<T>(2)
	{
		this->m_Q(0,0) = T(1)/(a*a);
		this->m_Q(1,1) = T(1)/(b*b);

		this->m_s = T(-1);
	}

	QuadEllipsoid(T a, T b, T c) : Quadric<T>(3)
	{
		this->m_Q(0,0) = T(1)/(a*a);
		this->m_Q(1,1) = T(1)/(b*b);
		this->m_Q(2,2) = T(1)/(c*c);

		this->m_s = T(-1);
	}

	QuadEllipsoid(T a, T b, T c, T d) : Quadric<T>(4)
	{
		this->m_Q(0,0) = T(1)/(a*a);
		this->m_Q(1,1) = T(1)/(b*b);
		this->m_Q(2,2) = T(1)/(c*c);
		this->m_Q(3,3) = T(1)/(d*d);

		this->m_s = T(-1);
	}

	~QuadEllipsoid() {}

	/**
	 * only valid in principal axis system
	 */
	T GetRadius(std::size_t i) const
	{
		return std::abs(this->m_s) /
			std::sqrt(std::abs(this->m_Q(i,i)));
	}

	T GetVolume() const
	{
		return get_ellipsoid_volume(this->m_Q) /
			std::abs(this->m_s);
	}
};


//------------------------------------------------------------------------------


template<typename T = double>
std::vector<std::size_t> find_zeroes(std::size_t N, const T* pIn, T eps = tl::get_epsilon<T>())
{
	using t_vec = ublas::vector<T>;
	std::vector<std::size_t> vecIndices;

	for(std::size_t i=0; i<N-1; ++i)
	{
		t_vec zero(2);
		zero[0] = zero[1] = 0.;
		t_vec xdir(2);
		xdir[0] = 1.; xdir[1] = 0.;
		Line<T> xaxis(zero, xdir);

		t_vec pos0(2);
		pos0[0] = 0.; pos0[1] = pIn[i];
		t_vec pos1(2);
		pos1[0] = 1.; pos1[1] = pIn[i+1];
		Line<T> line(pos0, pos1-pos0);

		T param{};
		if(!line.intersect(xaxis, param, eps))
			continue;

		t_vec posInters = line(param);
		if(posInters[0]>=0. && posInters[0]<=1.)
			vecIndices.push_back(i);
	}

	return vecIndices;
}

}

#endif
