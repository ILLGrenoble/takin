/**
 * pre-defined sqw modules
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015 -- 2018
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

#ifndef __MCONV_SQWMOD_PHONON_H__
#define __MCONV_SQWMOD_PHONON_H__

#include <list>
#include <unordered_map>

#include "tlibs/helper/boost_hacks.h"
#include <boost/numeric/ublas/vector.hpp>

#include "tlibs/math/math.h"
#include "tlibs/math/kd.h"
#include "tlibs/file/loaddat.h"
#include "../../res/defs.h"
#include "../sqwbase.h"

#ifdef USE_RTREE
	#include "tlibs/math/rt.h"
	#define RT_ELEMS 64
#endif


namespace ublas = boost::numeric::ublas;


/**
 * simple phonon model
 */
class SqwPhonon : public SqwBase
{
private:
	SqwPhonon() {};

protected:
	static t_real_reso phonon_disp(t_real_reso dq, t_real_reso da, t_real_reso df);

	void create();
	void destroy();

protected:
#ifdef USE_RTREE
	std::shared_ptr<tl::Rt<t_real_reso, 3, RT_ELEMS>> m_rt;
#else
	std::shared_ptr<tl::Kd<t_real_reso>> m_kd;
#endif
	unsigned int m_iNumqs = 250;
	unsigned int m_iNumArc = 50;
	t_real_reso m_dArcMax = 10.;

	ublas::vector<t_real_reso> m_vecBragg;

	ublas::vector<t_real_reso> m_vecLA;
	ublas::vector<t_real_reso> m_vecTA1;
	ublas::vector<t_real_reso> m_vecTA2;

	t_real_reso m_dLA_amp=20., m_dLA_freq=M_PI, m_dLA_E_HWHM=0.1, m_dLA_q_HWHM=0.1, m_dLA_S0=1.;
	t_real_reso m_dTA1_amp=15., m_dTA1_freq=M_PI, m_dTA1_E_HWHM=0.1, m_dTA1_q_HWHM=0.1, m_dTA1_S0=1.;
	t_real_reso m_dTA2_amp=10., m_dTA2_freq=M_PI, m_dTA2_E_HWHM=0.1, m_dTA2_q_HWHM=0.1, m_dTA2_S0=1.;

	t_real_reso m_dIncAmp=0., m_dIncSig=0.1;
	t_real_reso m_dT = 100.;
	bool m_bUseQDistr = true;

public:
	SqwPhonon(const ublas::vector<t_real_reso>& vecBragg,
		const ublas::vector<t_real_reso>& vecTA1,
		const ublas::vector<t_real_reso>& vecTA2,
		t_real_reso dLA_amp, t_real_reso dLA_freq, t_real_reso dLA_E_HWHM, t_real_reso dLA_q_HWHM, t_real_reso dLA_S0,
		t_real_reso dTA1_amp, t_real_reso dTA1_freq, t_real_reso dTA1_E_HWHM, t_real_reso dTA1_q_HWHM, t_real_reso dTA1_S0,
		t_real_reso dTA2_amp, t_real_reso dTA2_freq, t_real_reso dTA2_E_HWHM, t_real_reso dTA2_q_HWHM, t_real_reso dTA2_S0,
		t_real_reso dT);
	SqwPhonon(const char* pcFile);

	virtual ~SqwPhonon() = default;

	//virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
	//	disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override;
	virtual t_real_reso operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;

	const ublas::vector<t_real_reso>& GetBragg() const { return m_vecBragg; }
	const ublas::vector<t_real_reso>& GetLA() const { return m_vecLA; }
	const ublas::vector<t_real_reso>& GetTA1() const { return m_vecTA1; }
	const ublas::vector<t_real_reso>& GetTA2() const { return m_vecTA2; }

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;
};


/**
 * even simpler phonon model
 */
class SqwPhononSingleBranch : public SqwBase
{
private:
	SqwPhononSingleBranch() {};

protected:
	static t_real_reso phonon_disp(t_real_reso dq, t_real_reso da, t_real_reso df);

protected:
	ublas::vector<t_real_reso> m_vecBragg;

	t_real_reso m_damp=20., m_dfreq=M_PI, m_dHWHM=0.1, m_dS0=1.;
	t_real_reso m_dIncAmp=0., m_dIncSig=0.1;
	t_real_reso m_dT = 100.;

public:
	SqwPhononSingleBranch(const char* pcFile);

	virtual ~SqwPhononSingleBranch() = default;

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override;
	virtual t_real_reso
		operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;

	const ublas::vector<t_real_reso>& GetBragg() const { return m_vecBragg; }

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;
};


#endif
