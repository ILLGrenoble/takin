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

#ifndef __MONTERES_H__
#define __MONTERES_H__

#include "../res/defs.h"
#include "tlibs/math/linalg.h"
#include <utility>
namespace ublas = boost::numeric::ublas;


struct Resolution
{
	// covariance matrix
	ublas::matrix<t_real_reso> cov;

	// resolution matrix
	bool bHasRes = 0;
	ublas::matrix<t_real_reso> res;

	// full-widths (coh)
	std::vector<t_real_reso> dQ;	// in 1/A and meV

	// full-width (inc)
	t_real_reso dEinc;

	// ellipse origin
	ublas::vector<t_real_reso> Q_avg, Q_avg_notrafo;

	// all MC events
	std::vector<ublas::vector<t_real_reso>> vecQ;
};


void normalise_P(std::vector<t_real_reso>*);


Resolution calc_res(const std::vector<ublas::vector<t_real_reso>>& Q_vec,
	const std::vector<t_real_reso> *pp_vec = nullptr,
	const ublas::vector<t_real_reso> *qPara = nullptr,
	const ublas::vector<t_real_reso> *qPerp = nullptr);


Resolution calc_res(const std::vector<ublas::vector<t_real_reso>>& vecKi,
	const std::vector<ublas::vector<t_real_reso>>& vecKf,
	const std::vector<t_real_reso>* p_i = nullptr,
	const std::vector<t_real_reso>* p_f = nullptr,
	const ublas::vector<t_real_reso> *qPara = nullptr,
	const ublas::vector<t_real_reso> *qPerp = nullptr);

#endif
