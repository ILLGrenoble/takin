/**
 * advanced linalg helpers
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 2013-2018
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

#ifndef __TLIBS_LINALG2_IMPL_H__
#define __TLIBS_LINALG2_IMPL_H__

#include "linalg2.h"
#include <memory>

#ifndef NO_LAPACK
extern "C"
{
	#define lapack_complex_float std::complex<float>
	#define lapack_complex_float_real(c) (c.real())
	#define lapack_complex_float_imag(c) (c.imag())
	#define lapack_complex_double std::complex<double>
	#define lapack_complex_double_real(c) (c.real())
	#define lapack_complex_double_imag(c) (c.imag())

	#include <lapacke.h>
}


namespace tl {

// selects the float or double version of a lapack function
template<class T1, class T2, class F1, class F2>
struct select_func
{
	F1* m_f1 = nullptr;
	F2* m_f2 = nullptr;

	select_func(F1* f1, F2* f2) : m_f1(f1), m_f2(f2) {}

	template<class T>
	typename std::enable_if<std::is_same<T, T1>::value, F1*>::type
		get_func() { return m_f1; }
	template<class T>
	typename std::enable_if<std::is_same<T, T2>::value, F2*>::type
		get_func() { return m_f2; }
};

// ----------------------------------------------------------------------------

template<class T>
bool eigenvec(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs_real,
	std::vector<ublas::vector<T>>& evecs_imag,
	std::vector<T>& evals_real,
	std::vector<T>& evals_imag,
	bool bNorm)
{
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_sgeev), decltype(LAPACKE_dgeev)>
		sfunc(LAPACKE_sgeev, LAPACKE_dgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs_real.resize(iOrder); evecs_imag.resize(iOrder);
	evals_real.resize(iOrder); evals_imag.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
	{
		evecs_real[i].resize(iOrder);
		evecs_imag[i].resize(iOrder);
	}

	std::unique_ptr<T, std::default_delete<T[]>>
		uptrMem(new T[iOrder*iOrder + iOrder*iOrder + iOrder*iOrder]);
	T *pMatrix = uptrMem.get();
	T *pEVs = pMatrix + iOrder*iOrder;

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'V', iOrder,
		pMatrix, iOrder, evals_real.data(), evals_imag.data(),
		nullptr, iOrder, pEVs, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve general real eigenproblem", 
			" (lapack error ", iInfo , ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		bool bIsReal = 0;
		if(float_equal<T>(evals_imag[i], 0.))
			bIsReal = 1;

		if(bIsReal)
		{
			for(std::size_t j=0; j<iOrder; ++j)
			{
				evecs_real[i][j] = pEVs[j*iOrder + i];
				evecs_imag[i][j] = 0.;
			}
		}
		else
		{
			for(std::size_t j=0; j<iOrder; ++j)
			{
				evecs_real[i][j] = pEVs[j*iOrder + i];
				evecs_imag[i][j] = pEVs[j*iOrder + i+1];

				evecs_real[i+1][j] = pEVs[j*iOrder + i];
				evecs_imag[i+1][j] = -pEVs[j*iOrder + i+1];
			}
			++i; // check: (next eigenval) == -(currrent eigenval)
		}
	}


	// normalise
	if(bNorm && bOk)
	{
		for(std::size_t i=0; i<evecs_real.size(); ++i)
		{
			T len = T(0);
			for(std::size_t j=0; j<evecs_real[i].size(); ++j)
				len += evecs_real[i][j]*evecs_real[i][j] + evecs_imag[i][j]*evecs_imag[i][j];
			len = std::sqrt(len);

			evecs_real[i] /= len;
			evecs_imag[i] /= len;
		}
	}

	return bOk;
}


template<class T>
bool eigenval(const ublas::matrix<T>& mat, std::vector<T>& evals_real, std::vector<T>& evals_imag)
{
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_sgeev), decltype(LAPACKE_dgeev)>
		sfunc(LAPACKE_sgeev, LAPACKE_dgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evals_real.resize(iOrder); evals_imag.resize(iOrder);

	std::unique_ptr<T, std::default_delete<T[]>> uptrMem(new T[iOrder*iOrder]);
	T *pMatrix = uptrMem.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'N', iOrder,
		pMatrix, iOrder, evals_real.data(), evals_imag.data(),
		nullptr, iOrder, nullptr, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve general real eigenproblem", 
			" (lapack error ", iInfo , ").");
		bOk = false;
	}

	return bOk;
}


template<class T>
bool eigenvec_cplx(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<std::complex<T>>& evals,
	bool bNorm)
{
	using t_cplx = std::complex<T>;
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_cgeev), decltype(LAPACKE_zgeev)>
		sfunc(LAPACKE_cgeev, LAPACKE_zgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMem(new t_cplx[iOrder*iOrder + iOrder*iOrder + iOrder]);
	t_cplx *pMatrix = uptrMem.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	t_cplx *pEVs = pMatrix + iOrder*iOrder;
	t_cplx *pEVals = pEVs + iOrder*iOrder;

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'V', iOrder,
		pMatrix, iOrder, pEVals,
		nullptr, iOrder, pEVs, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve general complex eigenproblem",
			" (lapack error ", iInfo , ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pEVs[j*iOrder + i];
		evals[i] = pEVals[i];

		if(bNorm && bOk)
			evecs[i] /= veclen(evecs[i]);
	}

	return bOk;
}


template<class T>
bool eigenval_cplx(const ublas::matrix<std::complex<T>>& mat, std::vector<std::complex<T>>& evals)
{
	using t_cplx = std::complex<T>;
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_cgeev), decltype(LAPACKE_zgeev)>
		sfunc(LAPACKE_cgeev, LAPACKE_zgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evals.resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMem(new t_cplx[iOrder*iOrder + iOrder]);
	t_cplx *pMatrix = uptrMem.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	t_cplx *pEVals = pMatrix + iOrder*iOrder;

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'N', iOrder,
		pMatrix, iOrder, pEVals,
		nullptr, iOrder, nullptr, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve general complex eigenproblem",
			" (lapack error ", iInfo , ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
		evals[i] = pEVals[i];

	return bOk;
}


// ----------------------------------------------------------------------------


template<class T>
bool eigenvec_sym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs,
	std::vector<T>& evals,
	bool bNorm)
{
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_ssyev), decltype(LAPACKE_dsyev)>
		sfunc(LAPACKE_ssyev, LAPACKE_dsyev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<T, std::default_delete<T[]>>
		uptrMat(new T[iOrder*iOrder]);
	T *pMatrix = uptrMat.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = (j>=i ? mat(i,j) : T(0));

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'V', 'U',
		iOrder, pMatrix, iOrder, evals.data());

	if(iInfo!=0)
	{
		log_err("Could not solve symmetric eigenproblem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pMatrix[j*iOrder + i];

		if(bNorm && bOk)
			evecs[i] /= veclen(evecs[i]);
	}

	return bOk;
}


template<class T>
bool eigenval_sym(const ublas::matrix<T>& mat, std::vector<T>& evals)
{
	bool bOk = true;
	select_func<float, double, decltype(LAPACKE_ssyev), decltype(LAPACKE_dsyev)>
		sfunc(LAPACKE_ssyev, LAPACKE_dsyev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evals.resize(iOrder);

	std::unique_ptr<T, std::default_delete<T[]>>
		uptrMat(new T[iOrder*iOrder]);
	T *pMatrix = uptrMat.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'U',
		iOrder, pMatrix, iOrder, evals.data());

	if(iInfo!=0)
	{
		log_err("Could not solve symmetric eigenproblem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	return bOk;
}


template<class T>
bool eigenvec_herm(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<T>& evals,
	bool bNorm)
{
	using t_cplx = std::complex<T>;
	bool bOk = true;

	select_func<float, double, decltype(LAPACKE_cheev), decltype(LAPACKE_zheev)>
		sfunc(LAPACKE_cheev, LAPACKE_zheev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMat(new t_cplx[iOrder*iOrder]);
	t_cplx *pMatrix = uptrMat.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = (j>=i ? mat(i,j) : T(0));

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'V', 'U',
		iOrder, pMatrix, iOrder, evals.data());

	if(iInfo != 0)
	{
		log_err("Could not solve hermitian eigenproblem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pMatrix[j*iOrder + i];
		if(bNorm)
			evecs[i] /= veclen(evecs[i]);
	}
	return bOk;
}


template<class T>
bool eigenval_herm(const ublas::matrix<std::complex<T>>& mat, std::vector<T>& evals)
{
	using t_cplx = std::complex<T>;
	bool bOk = true;

	select_func<float, double, decltype(LAPACKE_cheev), decltype(LAPACKE_zheev)>
		sfunc(LAPACKE_cheev, LAPACKE_zheev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evals.resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMat(new t_cplx[iOrder*iOrder]);
	t_cplx *pMatrix = uptrMat.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'U',
		iOrder, pMatrix, iOrder, evals.data());

	if(iInfo != 0)
	{
		log_err("Could not solve hermitian eigenproblem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	return bOk;
}


template<class T>
bool eigenvecsel_herm(const ublas::matrix<std::complex<T>>& mat,
	std::vector<ublas::vector<std::complex<T>>>& evecs,
	std::vector<T>& evals,
	bool bNorm, T minval, T maxval, T eps)
{
	// select needed functions
	select_func<float, double, decltype(LAPACKE_cheevr), decltype(LAPACKE_zheevr)>
	sfunc(LAPACKE_cheevr, LAPACKE_zheevr);
	auto pfunc = sfunc.get_func<T>();

	select_func<float, double, decltype(LAPACKE_slamch), decltype(LAPACKE_dlamch)>
	_lamch(LAPACKE_slamch, LAPACKE_dlamch);
	auto lamch = _lamch.get_func<T>();


	using t_cplx = std::complex<T>;
	bool bOk = true;

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>>
		uptrMat(new t_cplx[iOrder*iOrder + iOrder*iOrder]);
	t_cplx *pMatrix = uptrMat.get();
	t_cplx *pEVsOrtho = pMatrix + iOrder*iOrder;

	std::unique_ptr<int, std::default_delete<int[]>>
		uptrIdxArr(new int[2*iOrder]);
	int *pIdxArr = uptrIdxArr.get();

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = (j>=i ? mat(i,j) : T(0));

	// use maximum precision if none given
	if(eps < T(0))
		eps = lamch('S');

	// if an invalid range is given, select all eigenvalues
	bool bSelectAll = (minval > maxval);

	int minidx = 1, maxidx = iOrder;
	int iNumFound = 0;
	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'V', bSelectAll?'A':'V', 'U',
		iOrder, pMatrix, iOrder, minval, maxval, minidx, maxidx,
		eps, &iNumFound, evals.data(), pEVsOrtho, iOrder, pIdxArr);

	if(iInfo != 0)
	{
		log_err("Could not solve hermitian eigenproblem",
				" (lapack error ", iInfo, ").");
		bOk = false;
	}

	evecs.resize(iNumFound);
	evals.resize(iNumFound);
	for(int i = 0; i < iNumFound; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = /*pMatrix*/pEVsOrtho[j*iOrder + i];
		if(bNorm)
			evecs[i] /= veclen(evecs[i]);
	}
	return bOk;
}


// ----------------------------------------------------------------------------


template<typename T>
bool singvec(const ublas::matrix<T>& mat,
	ublas::matrix<T>& matU, ublas::matrix<T>& matV, std::vector<T>& vecsvals)
{
	select_func<float, double, decltype(LAPACKE_sgesvd), decltype(LAPACKE_dgesvd)>
		sfunc(LAPACKE_sgesvd, LAPACKE_dgesvd);
	auto pfunc = sfunc.get_func<T>();

	const std::size_t iM = mat.size1();
	const std::size_t iN = mat.size2();
	const std::size_t iMin = std::min(iM,iN);

	vecsvals.resize(iMin);
	matU.resize(iM, iM);
	matV.resize(iN, iN);

	std::unique_ptr<T, std::default_delete<T[]>> uptrMat(new T[iM*iN]);
	std::unique_ptr<T, std::default_delete<T[]>> uptrWork(new T[iM*iN]);	// TODO: find correct size
	std::unique_ptr<T, std::default_delete<T[]>> uptrU(new T[iM*iM]);
	std::unique_ptr<T, std::default_delete<T[]>> uptrVt(new T[iN*iN]);

	for(std::size_t i=0; i<iM; ++i)
		for(std::size_t j=0; j<iN; ++j)
			uptrMat.get()[i*iN + j] = mat(i,j);


	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'A', 'A', iM, iN, uptrMat.get(), iN,
		vecsvals.data(), uptrU.get(), iM, uptrVt.get(), iN, uptrWork.get());

	bool bOk = true;
	if(iInfo != 0)
	{
		log_err("Could not solve real singular value problem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iM; ++i)
		for(std::size_t j=0; j<iM; ++j)
			matU(i,j) = uptrU.get()[i*iM + j];
	for(std::size_t i=0; i<iN; ++i)
		for(std::size_t j=0; j<iN; ++j)
			matV(j,i) = uptrVt.get()[i*iN + j];	// transposed

	return bOk;
}


template<typename T>
bool singvec_cplx(const ublas::matrix<std::complex<T>>& mat,
	ublas::matrix<std::complex<T>>& matU, ublas::matrix<std::complex<T>>& matV,
	std::vector<T>& vecsvals)
{
	using t_cplx = std::complex<T>;

	select_func<float, double, decltype(LAPACKE_cgesvd), decltype(LAPACKE_zgesvd)>
		sfunc(LAPACKE_cgesvd, LAPACKE_zgesvd);
	auto pfunc = sfunc.get_func<T>();

	const std::size_t iM = mat.size1();
	const std::size_t iN = mat.size2();
	const std::size_t iMin = std::min(iM,iN);

	vecsvals.resize(iMin);
	matU.resize(iM, iM);
	matV.resize(iN, iN);

	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>> uptrMat(new t_cplx[iM*iN]);
	std::unique_ptr<T, std::default_delete<T[]>> uptrWork(new T[iM*iN]);	// TODO: find correct size
	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>> uptrU(new t_cplx[iM*iM]);
	std::unique_ptr<t_cplx, std::default_delete<t_cplx[]>> uptrVt(new t_cplx[iN*iN]);

	for(std::size_t i=0; i<iM; ++i)
		for(std::size_t j=0; j<iN; ++j)
			uptrMat.get()[i*iN + j] = mat(i,j);


	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'A', 'A', iM, iN, uptrMat.get(), iN,
		vecsvals.data(), uptrU.get(), iM, uptrVt.get(), iN, uptrWork.get());

	bool bOk = true;
	if(iInfo != 0)
	{
		log_err("Could not solve complex singular value problem",
			" (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iM; ++i)
		for(std::size_t j=0; j<iM; ++j)
			matU(i,j) = uptrU.get()[i*iM + j];
	for(std::size_t i=0; i<iN; ++i)
		for(std::size_t j=0; j<iN; ++j)
			matV(j,i) = uptrVt.get()[i*iN + j];	// transposed

	return bOk;
}


// ----------------------------------------------------------------------------


template<class T>
bool qr(const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R)
{
	select_func<float, double, decltype(LAPACKE_sgeqrf), decltype(LAPACKE_dgeqrf)>
		sfunc(LAPACKE_sgeqrf, LAPACKE_dgeqrf);
	auto pfunc = sfunc.get_func<T>();

	const typename ublas::matrix<T>::size_type m = M.size1();
	const typename ublas::matrix<T>::size_type n = M.size2();

	const std::size_t iTauSize = m; //std::min<std::size_t>(m,n);

	std::unique_ptr<T, std::default_delete<T[]>>
		uptrMem(new T[n*m + iTauSize]);
	T *pMem = uptrMem.get();

	T *pMat = pMem;
	T *pTau = pMem + n*m;

	for(std::size_t i=0; i<m; ++i)
		for(std::size_t j=0; j<n; ++j)
			pMat[i*n + j] = M(i,j);

	// see: http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html
	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, m, n, pMat, n, pTau);

	R = ublas::matrix<T>(m,n);
	for(std::size_t i = 0; i < m; ++i)
		for(std::size_t j = 0; j < n; ++j)
		{
			if(j >= i)
				R(i,j) = pMat[i*n + j];
			else
				R(i,j) = 0.;
		}

	ublas::vector<T> v(iTauSize);

	const ublas::matrix<T> ident = unit_m<ublas::matrix<T>>(iTauSize);
	Q = ident;

	for(std::size_t k = 1; k <= iTauSize; ++k)
	{
		T dTau = pTau[k-1];

		for(std::size_t i = 1; i <= k-1; ++i)
			v[i - 1] = 0.;
		v[k - 1] = 1.;

		for(std::size_t i = k+1; i <= iTauSize; ++i)
			v[i-1] = pMat[(i-1)*n + (k-1)];

		ublas::matrix<T> VV = outer(v, v);
		ublas::matrix<T> H = ident - dTau*VV;

		Q = prod_mm(Q, H);
	}

	return (iInfo == 0);
}

}

#endif

#endif
