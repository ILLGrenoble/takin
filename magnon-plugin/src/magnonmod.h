/**
 * S(Q, E) module for magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv2, see 'LICENSE' file
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

#ifndef __MAGNON_SQW_MOD_H__
#define __MAGNON_SQW_MOD_H__

#include "core/tools/monteconvo/sqwbase.h"
#include "mag-core/libs/magdyn.h"


class MagnonMod : public SqwBase
{
	public:
		using SqwBase::t_var;

		using t_size = std::size_t;
		using t_real = t_real_reso;
		using t_cplx = std::complex<t_real>;

		using t_vec_cplx = tl2::vec<t_cplx>;
		using t_mat_cplx = tl2::mat<t_cplx>;

		using t_magdyn = magdyn::MagDyn<
			t_mat_cplx, t_vec_cplx,
			t_cplx, t_real, t_size>;

		using t_vec3_real = typename t_magdyn::t_vec3_real;
		using t_mat33_real = typename t_magdyn::t_mat33_real;

		using EnergiesAndWeights = typename t_magdyn::EnergiesAndWeights;


	protected:
		t_magdyn m_dyn{};

		// peak width
		t_real m_sigma { t_real(0.025) };
		int m_lineshape { 0 };  // 0: gaussian, 1: lorentzian, 2: rectangular

		// S(q, E) scaling factor
		t_real m_S0 { 1. };

		// incoherent amplitude and width
		t_real m_incoh_amp { 0. };
		t_real m_incoh_sigma { 0.025 };
		int m_incoh_lineshape { 0 };  // 0: gaussian, 1: lorentzian, 2: rectangular

		// polarisation channel, -1: unpolarised
		int m_channel { -1 };
		// mode index, -1: all
		int m_mode_idx { -1 };

		// use model's implementation of the bose factor
		bool m_use_model_bose{ false };
		// temperature
		t_real m_T { 300 };

		// unit cell from 0 to 1 (or from -0.5 to 0.5)?
		bool m_uc_01 { false };

		// rotate spin-spin correlation matrix
		bool m_use_polcoords { false };

		// use powder calculation (otherwise single crystal)?
		bool m_is_powder { false };

		// number of Q points to use for powder
		unsigned int m_powder_Qs { 512 };

		// twinning
		bool m_use_twinning { false };
		t_vec3_real m_twinning_axis { 0., 0., 1 };
		t_real m_twinning_angle { 5. };
		t_real m_twinning_fraction { 0.5 };  // population factor for main grain


#ifdef MAGNONMOD_ALLOW_QSIGNS
		// for quickly flipping coordinates
		std::vector<t_real> m_Qsigns { 1., 1., 1. };
#endif


	public:
		MagnonMod();
		MagnonMod(const std::string& strCfgFile);
		virtual ~MagnonMod();

		virtual std::tuple<std::vector<t_real>, std::vector<t_real>>
			disp(t_real dh, t_real dk, t_real dl) const override;
		virtual t_real operator()(t_real dh, t_real dk, t_real dl, t_real dE) const override;

		virtual std::vector<t_var> GetVars() const override;
		virtual void SetVars(const std::vector<t_var>&) override;
		virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override;

		virtual SqwBase* shallow_copy() const override;
};

#endif
