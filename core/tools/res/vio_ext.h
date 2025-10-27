/**
 * re-implementation of V. Mecoli's extension of Violini's TOF reso algorithm [mec25]
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date oct-2025
 * @license GPLv2
 *
 * @desc for extended method, see: [mec25] V. Mecoli, PhD thesis in preparation
 * @desc for original method, see: [vio14] N. Violini et al., NIM A 736 (2014) pp. 31-39, doi: 10.1016/j.nima.2013.10.042
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TOFRESO_EXT_H__
#define __TOFRESO_EXT_H__

#include "vio.h"


/**
 * TOF parameters in sigma
 */
struct VioExtParams : public VioParams
{
};


extern ResoResults calc_vio_ext(const VioExtParams& params);


#endif
