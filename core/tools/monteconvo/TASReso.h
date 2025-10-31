/**
 * loads reso settings
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jul-2015
 * @license GPLv2
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

#ifndef __DO_RESO_H__
#define __DO_RESO_H__

#include "../res/algos/eck_ext.h"
#include "../res/algos/vio_ext.h"
#include "../res/ellipse.h"
#include "../res/mc.h"

#include <vector>


enum class ResoFocus : unsigned
{
	FOC_UNCHANGED = 0,

	FOC_MONO_FLAT = (1<<0),
	FOC_MONO_H = (1<<1),
	FOC_MONO_V = (1<<2),

	FOC_ANA_FLAT = (1<<3),
	FOC_ANA_H = (1<<4),
	FOC_ANA_V = (1<<5),
};


class TASReso
{
protected:
	ResoAlgo m_algo = ResoAlgo::CN;
	ResoFocus m_foc = ResoFocus::FOC_UNCHANGED;

	McNeutronOpts<ublas::matrix<t_real_reso>> m_opts;
	EckParams m_reso;
	VioExtParams m_tofreso;

	// randomly smear out sample position if vector size >= 1
	std::vector<ResoResults> m_res;

	bool m_bKiFix = 0;
	t_real_reso m_dKFix = 1.4;

	// for better compatibility with older versions
	t_real_reso m_R0_scale = 1.;
	t_real_reso m_dPlaneDistTolerance = std::cbrt(tl::get_epsilon<t_real_reso>());

public:
	TASReso();
	TASReso(const TASReso& res);
	const TASReso& operator=(const TASReso& res);

	virtual ~TASReso() = default;

	bool LoadRes(const char* pcXmlFile);
	bool LoadLattice(const char* pcXmlFile, bool flip_coords = false);

	bool SetLattice(t_real_reso a, t_real_reso b, t_real_reso c,
		t_real_reso alpha, t_real_reso beta, t_real_reso gamma,
		const ublas::vector<t_real_reso>& vec1, const ublas::vector<t_real_reso>& vec2);
	bool SetHKLE(t_real_reso h, t_real_reso k, t_real_reso l, t_real_reso E);
	Ellipsoid4d<t_real_reso> GenerateMC(std::size_t iNum, std::vector<ublas::vector<t_real_reso>>&) const;
	Ellipsoid4d<t_real_reso> GenerateMC_deferred(std::size_t iNum, std::vector<ublas::vector<t_real_reso>>&) const;

	void SetKiFix(bool bKiFix) { m_bKiFix = bKiFix; }
	void SetKFix(t_real_reso dKFix) { m_dKFix = dKFix; }

	void SetAlgo(ResoAlgo algo) { m_algo = algo; }
	void SetOptimalFocus(ResoFocus foc) { m_foc = foc; }
	void SetPlaneDistTolerance(t_real_reso eps) { m_dPlaneDistTolerance = eps; }

	const EckParams& GetResoParams() const { return m_reso; }
	const VioParams& GetTofResoParams() const { return m_tofreso; }
	const McNeutronOpts<ublas::matrix<t_real_reso>>& GetMCOpts() const { return m_opts; }
	EckParams& GetResoParams() { return m_reso; }
	VioParams& GetTofResoParams() { return m_tofreso; }
	t_real_reso GetR0Scale() const { return m_R0_scale; }

	const ResoResults& GetResoResults() const { return m_res[0]; }

	void SetRandomSamplePos(std::size_t iNum) { m_res.resize(iNum); }
};

#endif
