/**
 * S(Q, E) module example for a pre-calculated grid (file format version 2)
 * @author Tobias Weber <tweber@ill.fr>
 * @date 06-jan-2020
 * @license GPLv2
 */

#ifndef __MCONV_SQW_GRID2_H__
#define __MCONV_SQW_GRID2_H__

#include "tools/monteconvo/sqwbase.h"
#include "tlibs/math/linalg.h"


class SqwMod : public SqwBase
{
	public:
		using SqwBase::t_var;
		using t_real = t_real_reso;
		using t_vec = tl::ublas::vector<t_real>;

	protected:
		// temperature for Bose factor
		t_real m_dT = t_real(100);

		// Bose cutoff
		t_real m_dcut = t_real(0.02);

		// peak width for creation and annihilation
		t_real m_dSigma = t_real(0.05);

		// S(q,E) scaling factor
		t_real m_dS0 = t_real(1.);

		// incoherent amplitude and width
		t_real m_dIncAmp = t_real(0.);
		t_real m_dIncSigma = t_real(0.05);

		// grid data file
		std::string m_strDataFile;
		std::size_t m_indexBlockOffset = 0;

		t_real m_hmin=0., m_hmax=0., m_hstep=0.;
		t_real m_kmin=0., m_kmax=0., m_kstep=0.;
		t_real m_lmin=0., m_lmax=0., m_lstep=0.;
		std::size_t m_hsize=0, m_ksize=0, m_lsize=0;



	public:
		SqwMod();
		SqwMod(const std::string& strDatFile);
		virtual ~SqwMod();

		virtual std::tuple<std::vector<t_real>, std::vector<t_real>>
			disp(t_real dh, t_real dk, t_real dl) const override;
		virtual t_real operator()(t_real dh, t_real dk, t_real dl, t_real dE) const override;

		virtual std::vector<t_var> GetVars() const override;
		virtual void SetVars(const std::vector<t_var>&) override;
		virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override;

		virtual SqwBase* shallow_copy() const override;
};

#endif
