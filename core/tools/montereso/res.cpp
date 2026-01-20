/**
 * Monte-Carlo resolution calculation
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date July 2012, Sep. 2014
 * @license GPLv2
 *
 * References:
 *   - based on Mcstas' 'mcresplot.pl' perl program:
 *        - www.mcstas.org, https://github.com/McStasMcXtrace/McCode/blob/master/tools/Legacy-Perl/mcresplot.pl
 *        - 'mcresplot.pl' is first mentioned here: https://doi.org/10.1016/S0921-4526(99)01409-X
 *   - also based on the rescal5 matlab program: http://www.ill.eu/en/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab/)
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

#include "res.h"

#include "../res/ellipse.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/math/stat.h"
#include "tlibs/log/log.h"

#include <algorithm>
#include <boost/algorithm/minmax_element.hpp>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>

using namespace ublas;
using t_real = t_real_reso;


void normalise_P(std::vector<t_real>* pp_vec)
{
	if(!pp_vec)
		return;

	if(pp_vec->size())
	{
		const t_real dpMax = *std::max_element(pp_vec->begin(), pp_vec->end());
		for(t_real& dP : *pp_vec)
			dP /= dpMax;
	}
}


/*
 * this function tries to be a 1:1 C++ reimplementation of the Perl function
 * 'read_mcstas_res' of the McStas 'mcresplot.pl' program
 */
Resolution calc_res(const std::vector<vector<t_real>>& Q_vec, const std::vector<t_real>* pp_vec,
	const ublas::vector<t_real> *qPara, const ublas::vector<t_real> *qPerp)
{
	vector<t_real> Q_avg = pp_vec ? tl::mean_value(*pp_vec, Q_vec) : tl::mean_value(Q_vec);

	vector<t_real> Q_dir, Q_perp;

	// use user-given or average Q vector
	if(qPara)
		Q_dir = *qPara;
	else
		Q_dir = tl::make_vec({Q_avg[0], Q_avg[1], Q_avg[2]});
	Q_dir = Q_dir / norm_2(Q_dir);

	// vector perpendicular to Q, in the scattering plane
	if(qPerp)
		Q_perp = *qPerp;
	else
		Q_perp = tl::make_vec({-Q_dir[1], Q_dir[0], Q_dir[2]});
	Q_perp = Q_perp / norm_2(Q_perp);

	// normal to scattering plane
	vector<t_real> vecUp = tl::cross_3(Q_dir, Q_perp);


	tl::log_info("Q_para = ", Q_dir, ", Q_perp = ", Q_perp, ", Q_up = ", vecUp, ".");


	/*
	 * transformation from the <Q_x, Q_y, Q_z, E> system
	 * into the <Qperp (approx. Q_avg), Q_perp, Q_z, E> system,
	 * i.e. rotate by the (ki,Q) angle
	 */
	matrix<t_real> trafo = identity_matrix<t_real>(4);
	tl::set_column(trafo, 0, Q_dir);
	tl::set_column(trafo, 1, Q_perp);
	tl::set_column(trafo, 2, vecUp);
	//trafo(3,3) = -1;
	tl::log_info("Transformation: ", trafo);


	Resolution reso;
	reso.Q_avg_notrafo = Q_avg;
	reso.Q_avg = prod(trans(trafo), Q_avg);
	tl::log_info("Transformed average Q vector: ", reso.Q_avg);

	reso.res.resize(4, 4, 0);
	reso.cov.resize(4, 4, 0);

	std::tie(reso.cov, std::ignore) = tl::covariance(Q_vec, pp_vec);
	tl::log_info("Covariance matrix (untransformed): ", reso.cov);

	reso.cov = tl::transform<matrix<t_real>>(reso.cov, trafo, true);
	tl::log_info("Covariance matrix: ", reso.cov);

	if(!(reso.bHasRes = tl::inverse(reso.cov, reso.res)))
		tl::log_err("Covariance matrix could not be inverted!");

	if(reso.bHasRes)
	{
		reso.dQ = calc_bragg_fwhms(reso.res);
		reso.dEinc = calc_vanadium_fwhms(reso.res)[3];

		tl::log_info("Resolution matrix: ", reso.res);

		std::ostringstream ostrVals;
		ostrVals << "Coherent / Bragg FWHM values (Qx, Qy, Qz, E): ";
		std::copy(reso.dQ.begin(), reso.dQ.end(),
			std::ostream_iterator<t_real>(ostrVals, ", "));

		std::ostringstream ostrIncVals;
		ostrIncVals << "Incoherent / Vanadium FWHM value (E in meV): " << reso.dEinc;

		std::ostringstream ostrElli;
		ostrElli << "Ellipsoid offsets: ";
		std::copy(reso.Q_avg.begin(), reso.Q_avg.end(),
			std::ostream_iterator<t_real>(ostrElli, ", "));

		reso.vecQ = Q_vec;

		// transform monte carlo events into Q||, Qperp, Qup
		for(vector<t_real>& Q : reso.vecQ)
			Q = prod(tl::transpose(trafo), Q - reso.Q_avg_notrafo) /* + reso.Q_avg*/;


		tl::log_info(ostrVals.str());
		tl::log_info(ostrIncVals.str());
		tl::log_info(ostrElli.str());
	}

	return reso;
}



/*
 * this function tries to be a 1:1 C++ reimplementation of the Perl function
 * 'read_mcstas_res' of the McStas 'mcresplot.pl' program
 */
Resolution calc_res(
	const std::vector<ublas::vector<t_real>>& vecKi, const std::vector<ublas::vector<t_real>>& vecKf,
	const std::vector<t_real>* _p_i, const std::vector<t_real>* _p_f,
	const ublas::vector<t_real> *qPara, const ublas::vector<t_real> *qPerp)
{
	tl::log_info("Calculating resolution...");

	std::size_t uiLen = vecKi.size();
	std::vector<vector<t_real>> Q_vec;
	std::vector<t_real> p_vec;
	Q_vec.reserve(uiLen);
	p_vec.reserve(uiLen);

	const t_real pi_max = _p_i ? *std::max_element(_p_i->begin(), _p_i->end()) : 1.;
	const t_real pf_max = _p_f ? *std::max_element(_p_f->begin(), _p_f->end()) : 1.;
	const t_real p_max = std::abs(pi_max*pf_max);

	vector<t_real> Q_avg(4);
	Q_avg[0] = Q_avg[1] = Q_avg[2] = Q_avg[3] = 0.;

	t_real p_sum = 0.;
	for(std::size_t uiRow=0; uiRow<uiLen; ++uiRow)
	{
		const vector<t_real>& ki = vecKi[uiRow];
		const vector<t_real>& kf = vecKf[uiRow];

		t_real p = (_p_i && _p_f) ? fabs((*_p_i)[uiRow] * (*_p_f)[uiRow]) : 1.;
		p /= p_max;		// normalize p to 0..1
		p_sum += p;

		vector<t_real> Q = ki - kf;
		t_real Ei = tl::get_KSQ2E<t_real>() * inner_prod(ki, ki);
		t_real Ef = tl::get_KSQ2E<t_real>() * inner_prod(kf, kf);
		t_real dE = Ei - Ef;

		// insert the energy into the Q vector
		Q.resize(4, true);
		Q[3] = dE;

		Q_avg += Q*p;

		Q_vec.push_back(Q);
		p_vec.push_back(p);
	}
	Q_avg /= p_sum;
	tl::log_info("Average Q vector: ", Q_avg);

	return calc_res(Q_vec, &p_vec, qPara, qPerp);
}
