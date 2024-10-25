/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc The magdyn library implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#ifndef __TLIBS2_MAGDYN_H__
#define __TLIBS2_MAGDYN_H__

#ifndef USE_LAPACK
	#define USE_LAPACK 1
#endif

// enables debug output
//#define __TLIBS2_MAGDYN_DEBUG_OUTPUT__

// enables ground state minimisation
//#define __TLIBS2_MAGDYN_USE_MINUIT__


#include "magdyn/magdyn.h"
#include "magdyn/structs.h"
#include "magdyn/helpers.h"
#include "magdyn/getters.h"
#include "magdyn/generators.h"
#include "magdyn/file.h"
#include "magdyn/configuration.h"
#include "magdyn/groundstate.h"
#include "magdyn/precalc.h"
#include "magdyn/hamilton.h"
#include "magdyn/correlation.h"
#include "magdyn/polarisation.h"
#include "magdyn/dispersion.h"


#endif
