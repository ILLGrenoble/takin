/**
 * resolution ellipse calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 14-may-2013
 * @license GPLv2
 *
 * @desc This is a reimplementation in C++ of the files rc_projs.m and rc_int.m of the
 *   - 'rescal5' package by Zinkin, McMorrow, Tennant, Farhi, and Wildes (ca. 1995-2007):
 *       http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/
 *   - and the 'mcresplot.pl' program from McStas (www.mcstas.org):
 *       https://github.com/McStasMcXtrace/McCode/blob/master/tools/Legacy-Perl/mcresplot.pl
 *   - see also: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __RES_ELLIPSE__
#define __RES_ELLIPSE__

#include <string>
#include <ostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <utility>

#include "tlibs/math/quat.h"
#include "tlibs/math/geo.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/linalg2.h"
#include "tlibs/math/math.h"
#include "defs.h"

namespace ublas = boost::numeric::ublas;


template<class t_real = t_real_reso>
struct Ellipse2d
{
	tl::QuadEllipsoid<t_real> quad;
	ublas::matrix<t_real> rot;

	t_real phi, slope;
	t_real x_hwhm, y_hwhm;
	t_real x_hwhm_bound, y_hwhm_bound;
	t_real x_offs, y_offs;
	t_real area;

	std::string x_lab, y_lab;

	ublas::vector<t_real> operator()(t_real t, bool bOffs=true) const;
	void GetCurvePoints(std::vector<t_real>& x, std::vector<t_real>& y,
		std::size_t iPoints=512, t_real *pLRTB=0);
};


template<class t_real = t_real_reso>
struct Ellipsoid3d
{
	tl::QuadEllipsoid<t_real> quad;
	ublas::matrix<t_real> rot;

	t_real x_hwhm, y_hwhm, z_hwhm;
	t_real x_offs, y_offs, z_offs;
	t_real vol;

	std::string x_lab, y_lab, z_lab;
};


template<class t_real = t_real_reso>
struct Ellipsoid4d
{
	tl::QuadEllipsoid<t_real> quad;
	ublas::matrix<t_real> rot;

	t_real x_hwhm, y_hwhm, z_hwhm, w_hwhm;
	t_real x_offs, y_offs, z_offs, w_offs;
	t_real vol;

	std::string x_lab, y_lab, z_lab, w_lab;
};


enum class EllipseCoordSys : int
{
	AUTO       = -1,

	Q_AVG      = 0,  // Q|| Qperp system (1/A)
	RLU        = 1,  // absolute hkl system (rlu)
	RLU_ORIENT = 2   // system using scattering plane (rlu)
};



// --------------------------------------------------------------------------------


/**
 * project along one axis of the quadratic part of the quadric to remove line and column iIdx
 * this is a 1:1 C++ reimplementation of 'rc_int' from 'mcresplot.pl' and 'rescal5'
 * (see also [eck14], equ. 57)
 * project along row/column iIdx
 *
 * quadric M: <x|M|x> = c
 * projector along |v>: |v><v|, with |v> normalised
 * remove projected contribution: M - s * |v><v|
 * choose column/row i as |v> and scale factor s so that they vanish with the projection
 *
 * cf. projector orthogonal to |v>: (1 - |v><v|) * M
 * cf. also householder mirror along |v>: (1 - 2*|v><v|) * M
 *
 * e.g. projection along 2nd column:
 *     ( m00 m01 m02 )        ( m01 )           ( m01 m01  m01 m11  m01 m12 )
 * M = ( m01 m11 m12 ), |v> = ( m11 ), |v><v| = ( m11 m01  m11 m11  m11 m12 )
 *     ( m02 m12 m22 )        ( m12 )           ( m12 m01  m12 m11  m12 m12 )
 *
 *                           ( m01 m01/m11    m01    m01 m12/m11 )
 * s * |v><v| = |v><v|/m11 = (         m01    m11            m12 )
 *                           ( m12 m01/m11    m12    m12 m12/m11 )
 */
template<class T = t_real_reso>
ublas::matrix<T> quadric_proj(const ublas::matrix<T>& quadric, std::size_t iIdx)
{
	using t_mat = ublas::matrix<T>;
	using t_vec = ublas::vector<T>;

	if(tl::float_equal<T>(quadric(iIdx, iIdx), T{0}))
	{
		tl::log_warn("Cannot project quadric, slicing instead.");
		return tl::remove_elems(quadric, iIdx);
	}

	// symmetric matrix -> col and row are equal to one another and to this average b
	t_vec b = T(0.5)*(tl::get_column(quadric, iIdx) + tl::get_row(quadric, iIdx));
	//T blen = tl::veclen(b);
	//b /= blen;

	t_mat m = quadric;
	//T dscale = blen / b[iIdx];	// == blen*blen / mat(iIdx,iIdx)
	T dscale = 1. / quadric(iIdx, iIdx);
	m -= dscale * tl::outer<t_vec,t_mat>(b, b);

	//tl::log_debug(quadric, " -> ", m);
	m = tl::remove_elems(m, iIdx);

	return m;
}


/**
 * project along one axis of the linear part of the quadric
 * (see [eck14], equ. 57)
 */
template<class T = t_real_reso>
ublas::matrix<T> quadric_proj_mat(const ublas::matrix<T>& mat,
	const ublas::matrix<T>& quadric, std::size_t iIdx)
{
	using t_mat = ublas::matrix<T>;
	using t_vec = ublas::vector<T>;

	if(tl::float_equal<T>(quadric(iIdx, iIdx), T{0}))
	{
		tl::log_warn("Cannot project quadric, slicing instead.");
		return tl::remove_elems(quadric, iIdx);
	}

	// symmetric matrix -> col and row are equal to one another and to this average b
	t_vec b = T(0.5)*(tl::get_column(quadric, iIdx) + tl::get_row(quadric, iIdx));
	b /= quadric(iIdx, iIdx);
	t_mat m = mat - tl::outer<t_vec,t_mat>(b, tl::get_column(mat, iIdx));

	//tl::log_debug(quadric, " -> ", m);
	m = tl::remove_elems(m, iIdx);

	return m;
}


/**
 * project along one axis of the linear part of the quadric
 * (see [eck14], equ. 57)
 */
template<class T = t_real_reso>
ublas::vector<T> quadric_proj(const ublas::vector<T>& vec,
	const ublas::matrix<T>& quadric, std::size_t iIdx)
{
	using t_vec = ublas::vector<T>;

	if(tl::float_equal<T>(quadric(iIdx, iIdx), T{0}))
	{
		tl::log_warn("Cannot project vector part of quadric, slicing instead.");
		return tl::remove_elem(vec, iIdx);
	}

	t_vec b = T(0.5)*(tl::get_column(quadric, iIdx) + tl::get_row(quadric, iIdx));

	t_vec vecProj = vec;
	T dscale = vecProj[iIdx] / quadric(iIdx, iIdx);
	vecProj -= dscale * b;

	vecProj = tl::remove_elem(vecProj, iIdx);
	return vecProj;
}


// --------------------------------------------------------------------------------

template<class t_real = t_real_reso>
std::ostream& operator<<(std::ostream& ostr, const Ellipse2d<t_real>& ell)
{
	ostr << "phi = " << tl::r2d(ell.phi) << " deg \n";
	ostr << "slope = " << ell.slope << " deg \n";
	ostr << "x_hwhm = " << ell.x_hwhm << ", ";
	ostr << "y_hwhm = " << ell.y_hwhm << "\n";
	ostr << "x_hwhm_bound = " << ell.x_hwhm_bound << ", ";
	ostr << "y_hwhm_bound = " << ell.y_hwhm_bound << "\n";
	ostr << "x_offs = " << ell.x_offs << ", ";
	ostr << "y_offs = " << ell.y_offs << "\n";
	ostr << "x_lab = " << ell.x_lab << ", ";
	ostr << "y_lab = " << ell.y_lab << "\n";
	ostr << "area = " << ell.area;

	return ostr;
}


template<class t_real = t_real_reso>
std::ostream& operator<<(std::ostream& ostr, const Ellipsoid4d<t_real>& ell)
{
	ostr << "x_hwhm = " << ell.x_hwhm << ", ";
	ostr << "y_hwhm = " << ell.y_hwhm << ", ";
	ostr << "z_hwhm = " << ell.z_hwhm << ", ";
	ostr << "w_hwhm = " << ell.w_hwhm << "\n";
	ostr << "x_offs = " << ell.x_offs << ", ";
	ostr << "y_offs = " << ell.y_offs << ", ";
	ostr << "z_offs = " << ell.z_offs << ", ";
	ostr << "w_offs = " << ell.w_offs << "\n";
	ostr << "x_lab = " << ell.x_lab << ", ";
	ostr << "y_lab = " << ell.y_lab << ", ";
	ostr << "z_lab = " << ell.z_lab << ", ";
	ostr << "w_lab = " << ell.w_lab << "\n";
	ostr << "volume = " << ell.vol;

	return ostr;
}


// --------------------------------------------------------------------------------


/**
 * get points on ellipse
 */
template<class t_real>
ublas::vector<t_real> Ellipse2d<t_real>::operator()(t_real t, bool bAddOffs) const
{
	ublas::vector<t_real> vec = tl::make_vec<ublas::vector<t_real>>
		({ x_hwhm * std::cos(t_real(2)*tl::get_pi<t_real>()*t),
			y_hwhm * std::sin(t_real(2)*tl::get_pi<t_real>()*t) });

	vec = ublas::prod(rot, vec);

	if(bAddOffs)
	{
		vec[0] += x_offs;
		vec[1] += y_offs;
	}

	return vec;
}


template<class t_real>
void Ellipse2d<t_real>::GetCurvePoints(std::vector<t_real>& x, std::vector<t_real>& y,
	std::size_t iPoints, t_real *pLRTB)
{
	x.resize(iPoints);
	y.resize(iPoints);

	for(std::size_t i=0; i<iPoints; ++i)
	{
		t_real dT = t_real(i)/t_real(iPoints-1);
		ublas::vector<t_real> vec = operator()(dT);

		x[i] = vec[0];
		y[i] = vec[1];
	}

	if(pLRTB)  // bounding rect
	{
		auto pairX = std::minmax_element(x.begin(), x.end());
		auto pairY = std::minmax_element(y.begin(), y.end());

		*(pLRTB+0) = *pairX.first;   // left
		*(pLRTB+1) = *pairX.second;  // right
		*(pLRTB+2) = *pairY.second;  // top
		*(pLRTB+3) = *pairY.first;   // bottom
	}
}


// --------------------------------------------------------------------------------


template<class T = t_real_reso>
static void quad_proj(tl::QuadEllipsoid<T>& quad, std::size_t iIdx)
{
	ublas::vector<T> vecRint = quadric_proj(quad.GetR(), quad.GetQ(), iIdx);
	ublas::matrix<T> matQint = quadric_proj(quad.GetQ(), iIdx);
	quad.RemoveElems(iIdx);
	quad.SetQ(matQint);
	quad.SetR(vecRint);
}


static const std::string g_strLabelsCentre[] = {"Q_{para}-<Q> (1/A)", "Q_{ortho}-<Q> (1/A)", "Q_z-<Q> (1/A)", "E (meV)"};
static const std::string g_strLabels[] = {"Q_{para} (1/A)", "Q_{ortho} (1/A)", "Q_z (1/A)", "E (meV)"};
static const std::string g_strLabelsHKLCentre[] = {"h-<h> (rlu)", "k-<k> (rlu)", "l-<l> (rlu)", "E (meV)"};
static const std::string g_strLabelsHKL[] = {"h (rlu)", "k (rlu)", "l (rlu)", "E (meV)"};
static const std::string g_strLabelsHKLOrient[] = {"Reflex 1 (rlu)", "Reflex 2 (rlu)", "Up (rlu)", "E (meV)"};


static inline const std::string& ellipse_labels(int iCoord, EllipseCoordSys sys, bool bCentre=0)
{
	switch(sys)
	{
		case EllipseCoordSys::RLU:
			if(bCentre)
				return g_strLabelsHKLCentre[iCoord];
			return g_strLabelsHKL[iCoord];
		case EllipseCoordSys::RLU_ORIENT:
			return g_strLabelsHKLOrient[iCoord];
		case EllipseCoordSys::AUTO:
		case EllipseCoordSys::Q_AVG:
		default:
			if(bCentre)
				return g_strLabelsCentre[iCoord];
			return g_strLabels[iCoord];
	}
}


/*
 * this is a 1:1 C++ reimplementation of 'proj_elip' from 'mcresplot.pl' and 'rescal5'
 * iX, iY: dimensions to plot
 * iInt, iInt2: dimensions to integrate
 * iRem1, iRem2: dimensions to remove
 */
template<class t_real = t_real_reso>
Ellipse2d<t_real> calc_res_ellipse(
	const ublas::matrix<t_real>& reso,	// quadratic part of quadric
	const ublas::vector<t_real>& reso_vec,	// linear part
	t_real reso_const,			// const part
	const ublas::vector<t_real>& Q_avg,
	int iX, int iY, int iInt1, int iInt2, int iRem1, int iRem2)
{
	Ellipse2d<t_real> ell;
	ell.quad.SetDim(4);
	ell.quad.SetQ(reso);
	ell.quad.SetR(reso_vec);
	//ell.quad.SetS(reso_const);

	ell.x_offs = ell.y_offs = 0.;

	// labels only valid for non-rotated system
	ell.x_lab = g_strLabels[iX];
	ell.y_lab = g_strLabels[iY];


	ublas::vector<t_real> Q_offs = Q_avg;

	if(iRem1 > -1)
	{
		ell.quad.RemoveElems(iRem1);
		Q_offs = tl::remove_elem(Q_offs, iRem1);

		if(iInt1 >= iRem1) --iInt1;
		if(iInt2 >= iRem1) --iInt2;
		if(iRem2 >= iRem1) --iRem2;
		if(iX >= iRem1) --iX;
		if(iY >= iRem1) --iY;
	}

	if(iRem2 > -1)
	{
		ell.quad.RemoveElems(iRem2);
		Q_offs = tl::remove_elem(Q_offs, iRem2);

		if(iInt1 >= iRem2) --iInt1;
		if(iInt2 >= iRem2) --iInt2;
		if(iX >= iRem2) --iX;
		if(iY >= iRem2) --iY;
	}

	if(iInt1 > -1)
	{
		quad_proj(ell.quad, iInt1);
		Q_offs = tl::remove_elem(Q_offs, iInt1);

		if(iInt2 >= iInt1) --iInt2;
		if(iX >= iInt1) --iX;
		if(iY >= iInt1) --iY;
	}

	if(iInt2 > -1)
	{
		quad_proj(ell.quad, iInt2);
		Q_offs = tl::remove_elem(Q_offs, iInt2);

		if(iX >= iInt2) --iX;
		if(iY >= iInt2) --iY;
	}

	std::vector<t_real> evals;

	tl::QuadEllipsoid<t_real> quad(2);
	ell.quad.GetPrincipalAxes(ell.rot, evals, &quad);

	ell.phi = tl::rotation_angle(ell.rot)[0];

	ell.x_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(0);
	ell.y_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(1);

	ell.x_offs = Q_offs[iX];
	ell.y_offs = Q_offs[iY];


	// calculate bounding rect
	std::tie(ell.x_hwhm_bound, ell.y_hwhm_bound)
		= [&ell]() -> std::pair<t_real, t_real>
	{
		t_real t0deg = ell.phi / (tl::get_pi<t_real>()*t_real(2));
		t_real t90deg = (ell.phi + tl::get_pi<t_real>()/t_real(2))
			/ (tl::get_pi<t_real>()*t_real(2));

		ublas::vector<t_real> vec1 = ell(t0deg, false);
		ublas::vector<t_real> vec2 = ell(t90deg, false);

		t_real dX = std::max(std::abs(vec1[0]), std::abs(vec2[0]));
		t_real dY = std::max(std::abs(vec1[1]), std::abs(vec2[1]));

		return std::pair<t_real, t_real>(dX, dY);
	}();


	// linear part of quadric
	const ublas::vector<t_real> vecTrans
		= ublas::prod(ell.rot, quad.GetPrincipalOffset());

	if(vecTrans.size() == 2)
	{
		ell.x_offs += vecTrans[0];
		ell.y_offs += vecTrans[1];
	}
	else
	{
		tl::log_err("Invalid ellipse shift.");
	}

	ell.area = quad.GetVolume();
	ell.slope = std::tan(ell.phi);

	return ell;
}


// --------------------------------------------------------------------------------


/**
 * transforms resolution ellipsoid from <Qpara Qperp Qup>
 * to crystal hkl coordinate system
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::tuple<t_mat, t_vec, t_vec>
conv_lab_to_rlu(t_real dAngleQVec0,
	const t_mat& matUB, const t_mat& matUBinv,
	const t_mat& reso, const t_vec& reso_v, const t_vec& Q_avg)
{
	// hkl crystal system:
	// Qavg system in 1/A -> rotate back to orient system in 1/A ->
	// transform to hkl rlu system
	t_mat matQVec0 = tl::rotation_matrix_2d(-dAngleQVec0);
	tl::resize_unity(matQVec0, 4);
	const t_mat matQVec0inv = ublas::trans(matQVec0);

	const t_mat matUBinvQVec0 = ublas::prod(matUBinv, matQVec0);
	const t_mat matQVec0invUB = ublas::prod(matQVec0inv, matUB);
	t_mat resoHKL = tl::transform(reso, matQVec0invUB, 1);
	t_vec Q_avgHKL = ublas::prod(matUBinvQVec0, Q_avg);

	t_vec reso_vHKL;
	if(reso_v.size() == 4)
		reso_vHKL = ublas::prod(matUBinvQVec0, reso_v);

	return std::make_tuple(resoHKL, reso_vHKL, Q_avgHKL);
}


/**
 * transforms resolution ellipsoid from <Qpara Qperp Qup>
 * to crystal orientation vector hkl system
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
std::tuple<t_mat, t_vec, t_vec>
conv_lab_to_rlu_orient(t_real dAngleQVec0,
	const t_mat& matUB, const t_mat& matUBinv,
	const t_mat& matUrlu, const t_mat& matUinvrlu,
	const t_mat& reso, const t_vec& reso_v, const t_vec& Q_avg)
{
	// hkl crystal system:
	// Qavg system in 1/A -> rotate back to orient system in 1/A ->
	// transform to hkl rlu system
	t_mat matQVec0 = tl::rotation_matrix_2d(-dAngleQVec0);
	tl::resize_unity(matQVec0, 4);
	const t_mat matQVec0inv = ublas::trans(matQVec0);

	const t_mat matUBinvQVec0 = ublas::prod(matUBinv, matQVec0);
	const t_mat matQVec0invUB = ublas::prod(matQVec0inv, matUB);

	// system of scattering plane: (orient1, orient2, up)
	// Qavg system in 1/A -> rotate back to orient system in 1/A ->
	// transform to hkl rlu system -> rotate forward to orient system in rlu
	const t_mat matToOrient = ublas::prod(matUrlu, matUBinvQVec0);
	const t_mat matToOrientinv = ublas::prod(matQVec0invUB, matUinvrlu);

	t_mat resoOrient = tl::transform(reso, matToOrientinv, 1);
	t_vec Q_avgOrient = ublas::prod(matToOrient, Q_avg);

	t_vec reso_vOrient;
	if(reso_v.size() == 4)
		reso_vOrient = ublas::prod(matToOrient, reso_v);

	return std::make_tuple(resoOrient, reso_vOrient, Q_avgOrient);
}


/**
 * coherent (Bragg) widths
 */
template<class t_mat, class t_real = typename t_mat::value_type>
std::vector<t_real> calc_bragg_fwhms(const t_mat& reso)
{
	static const t_real sig2fwhm = tl::get_SIGMA2FWHM<t_real>();

	std::vector<t_real> vecFwhms;
	vecFwhms.reserve(reso.size1());

	for(std::size_t i = 0; i < reso.size1(); ++i)
		vecFwhms.push_back(sig2fwhm/sqrt(reso(i,i)));

	return vecFwhms;
}


/**
 * incoherent (vanadium) widths
 */
template<class t_mat, class t_real = typename t_mat::value_type>
std::vector<t_real> calc_vanadium_fwhms(const t_mat& reso)
{
	std::vector<t_real> widths;

	// Q_parallel
	{
		ublas::matrix<t_real> M = quadric_proj(reso, 3);
		M = quadric_proj(M, 2);
		M = quadric_proj(M, 1);

		t_real w = tl::get_SIGMA2FWHM<t_real>() / std::sqrt(std::abs(M(0,0)));
		widths.push_back(w);
	}
	// Q_perpendicular
	{
		ublas::matrix<t_real> M = quadric_proj(reso, 3);
		M = quadric_proj(M, 2);
		M = quadric_proj(M, 0);

		t_real w = tl::get_SIGMA2FWHM<t_real>() / std::sqrt(std::abs(M(0,0)));
		widths.push_back(w);
	}
	// Q_up
	{
		ublas::matrix<t_real> M = quadric_proj(reso, 3);
		M = quadric_proj(M, 1);
		M = quadric_proj(M, 0);

		t_real w = tl::get_SIGMA2FWHM<t_real>() / std::sqrt(std::abs(M(0,0)));
		widths.push_back(w);
	}
	// energy -> project all Q axes
	{
		ublas::matrix<t_real> M = quadric_proj(reso, 0);
		M = quadric_proj(M, 0);
		M = quadric_proj(M, 0);

		t_real w = tl::get_SIGMA2FWHM<t_real>() / std::sqrt(std::abs(M(0,0)));
		widths.push_back(w);
	}

	return widths;
}


// --------------------------------------------------------------------------------


template<class t_real = t_real_reso>
Ellipsoid3d<t_real> calc_res_ellipsoid(
	const ublas::matrix<t_real>& reso,
	const ublas::vector<t_real>& reso_vec,
	t_real reso_const,
	const ublas::vector<t_real>& Q_avg,
	int iX, int iY, int iZ, int iInt, int iRem)
{
	Ellipsoid3d<t_real> ell;

	ell.quad.SetDim(4);
	ell.quad.SetQ(reso);
	ell.quad.SetR(reso_vec);
	//ell.quad.SetS(reso_const);

	ell.x_offs = ell.y_offs = ell.z_offs = 0.;

	// labels only valid for non-rotated system
	ell.x_lab = g_strLabels[iX];
	ell.y_lab = g_strLabels[iY];
	ell.z_lab = g_strLabels[iZ];

	ublas::vector<t_real> Q_offs = Q_avg;

	if(iRem >- 1)
	{
		ell.quad.RemoveElems(iRem);
		Q_offs = tl::remove_elem(Q_offs, iRem);

		if(iInt >= iRem) --iInt;
		if(iX >= iRem) --iX;
		if(iY >= iRem) --iY;
		if(iZ >= iRem) --iZ;
	}

	if(iInt >- 1)
	{
		quad_proj(ell.quad, iInt);
		Q_offs = tl::remove_elem(Q_offs, iInt);

		if(iX >= iInt) --iX;
		if(iY >= iInt) --iY;
		if(iZ >= iInt) --iZ;
	}

	std::vector<t_real> evals;
	tl::QuadEllipsoid<t_real> quad(3);
	ell.quad.GetPrincipalAxes(ell.rot, evals, &quad);

	//tl::log_info("Principal axes: ", quad.GetQ());
	ell.x_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(0);
	ell.y_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(1);
	ell.z_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(2);

	ell.x_offs = Q_offs[iX];
	ell.y_offs = Q_offs[iY];
	ell.z_offs = Q_offs[iZ];

	// linear part of quadric
	const ublas::vector<t_real> vecTrans = ublas::prod(ell.rot, quad.GetPrincipalOffset());

	if(vecTrans.size() == 3)
	{
		ell.x_offs += vecTrans[0];
		ell.y_offs += vecTrans[1];
		ell.z_offs += vecTrans[2];
	}
	else
	{
		tl::log_err("Invalid ellipsoid shift.");
	}

	ell.vol = quad.GetVolume();
	return ell;
}


// --------------------------------------------------------------------------------


template<class t_real = t_real_reso>
Ellipsoid4d<t_real> calc_res_ellipsoid4d(
	const ublas::matrix<t_real>& reso,
	const ublas::vector<t_real>& reso_vec,
	t_real reso_const,
	const ublas::vector<t_real>& Q_avg)
{
	Ellipsoid4d<t_real> ell;
	ell.quad.SetDim(4);
	ell.quad.SetQ(reso);
	ell.quad.SetR(reso_vec);
	//ell.quad.SetS(reso_const);

	std::vector<t_real> evals;
	tl::QuadEllipsoid<t_real> quad(4);
	ell.quad.GetPrincipalAxes(ell.rot, evals, &quad);

	ell.x_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(0);
	ell.y_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(1);
	ell.z_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(2);
	ell.w_hwhm = tl::get_SIGMA2HWHM<t_real>() * quad.GetRadius(3);

	ell.x_offs = Q_avg[0];
	ell.y_offs = Q_avg[1];
	ell.z_offs = Q_avg[2];
	ell.w_offs = Q_avg[3];

	// linear part of quadric
	const ublas::vector<t_real> vecTrans = ublas::prod(ell.rot, quad.GetPrincipalOffset());

	if(vecTrans.size() == 4)
	{
		ell.x_offs += vecTrans[0];
		ell.y_offs += vecTrans[1];
		ell.z_offs += vecTrans[2];
		ell.w_offs += vecTrans[3];
	}
	else
	{
		tl::log_err("Invalid ellipsoid shift.");
	}

	// labels only valid for non-rotated system
	ell.x_lab = g_strLabels[0];
	ell.y_lab = g_strLabels[1];
	ell.z_lab = g_strLabels[2];
	ell.w_lab = g_strLabels[3];

	ell.vol = quad.GetVolume();

	return ell;
}

#endif
