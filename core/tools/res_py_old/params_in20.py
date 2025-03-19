#
# in20 parameters (without and with flatcone)
#
# @author Tobias Weber <tweber@ill.fr>
# @date jul-2024
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

import numpy as np
import helpers


#
# mono vertical curvature
#
def mono_curv_v_formula(params):
    foclen = helpers.focal_len(params["dist_vsrc_mono"] / helpers.cm2A,
        params["dist_mono_sample"] / helpers.cm2A)

    curv_cm = helpers.foc_curv_2(foclen,
        params["ki"], params["mono_xtal_d"], True)

    # clamp at minimally possible curvature
    if curv_cm < 50.:
        curv_cm = 50.

    return curv_cm


#
# ana vertical curvature
#
def ana_curv_v_formula(params):
    foclen = helpers.focal_len(params["dist_sample_ana"] / helpers.cm2A,
        params["dist_ana_det"] / helpers.cm2A)

    curv_cm = helpers.foc_curv_2(foclen,
        params["kf"], params["ana_xtal_d"], True)

    return curv_cm


#
# mono horizontal curvature
#
def mono_curv_h_formula(params):
    foclen = helpers.focal_len(params["dist_hsrc_mono"] / helpers.cm2A,
        params["dist_mono_sample"] / helpers.cm2A)

    curv_cm = helpers.foc_curv_2(foclen,
        params["ki"], params["mono_xtal_d"], False)

    return curv_cm


#
# ana horizontal curvature
#
def ana_curv_h_formula(params):
    len1 = params["dist_sample_ana"] / helpers.cm2A
    #len2 = params["dist_ana_det"] / helpers.cm2A
    len2 = params["dist_sample_ana"] / helpers.cm2A
    foclen = helpers.focal_len(len1, len2)

    curv_cm = helpers.foc_curv_2(foclen,
        params["kf"], params["ana_xtal_d"], False)

    return curv_cm




#
# pre-defined parameters for IN20
#
params = {
    # options
    "verbose" : True,

    # resolution method, "eck", "pop", or "cn"
    "reso_method" : "eck",

    # scattering triangle
    "ki" : 2.662,
    "kf" : 2.662,
    "E"  : helpers.get_E(2.662, 2.662),
    "Q"  : 2.,

    # d spacings
    "mono_xtal_d" : 3.355,   # PG
    "ana_xtal_d"  : 3.355,   # PG

     # scattering senses
    "mono_sense"   : -1.,
    "sample_sense" :  1.,
    "ana_sense"    : -1.,
    "mirror_Qperp" : False,

    # distances
    "dist_vsrc_mono"   : 702.5 * helpers.cm2A,
    "dist_hsrc_mono"   : 227.5 * helpers.cm2A,
    "dist_mono_sample" : 227.5 * helpers.cm2A,
    "dist_sample_ana"  : 105.7 * helpers.cm2A,
    "dist_ana_det"     : 23.5  * helpers.cm2A,

    # shapes
    "src_shape"    : "rectangular",  # "rectangular" or "circular"
    "sample_shape" : "cylindrical",  # "cuboid" or "cylindrical"
    "det_shape"    : "rectangular",  # "rectangular" or "circular"

    # component sizes
    "src_w"    : 2.   * helpers.cm2A,
    "src_h"    : 17.  * helpers.cm2A,
    "mono_d"   : 0.2  * helpers.cm2A,
    "mono_w"   : 22.5 * helpers.cm2A,
    "mono_h"   : 20.  * helpers.cm2A,
    "sample_d" : 1.   * helpers.cm2A,
    "sample_w" : 1.   * helpers.cm2A,
    "sample_h" : 1.   * helpers.cm2A,
    "ana_d"    : 0.2  * helpers.cm2A,
    "ana_w"    : 23.  * helpers.cm2A,
    "ana_h"    : 9.6  * helpers.cm2A,
    "det_w"    : 2.5  * helpers.cm2A,
    "det_h"    : 2.5  * helpers.cm2A,

    # horizontal collimation
    "coll_h_pre_mono"    : 9999. * helpers.min2rad,
    "coll_h_pre_sample"  : 9999. * helpers.min2rad,
    "coll_h_post_sample" : 9999. * helpers.min2rad,
    "coll_h_post_ana"    : 9999. * helpers.min2rad,

    # vertical collimation
    "coll_v_pre_mono"    : 9999. * helpers.min2rad,
    "coll_v_pre_sample"  : 9999. * helpers.min2rad,
    "coll_v_post_sample" : 9999. * helpers.min2rad,
    "coll_v_post_ana"    : 9999. * helpers.min2rad,

    # horizontal focusing
    "mono_curv_h" : 0.,
    "ana_curv_h"  : 0.,
    "mono_is_curved_h" : True,
    "ana_is_curved_h"  : False,
    "mono_is_optimally_curved_h" : False,
    "ana_is_optimally_curved_h"  : False,
    "mono_curv_h_formula" : mono_curv_h_formula,
    "ana_curv_h_formula" : ana_curv_h_formula,

    # vertical focusing
    "mono_curv_v" : 0.,
    "ana_curv_v"  : 37. * helpers.cm2A,
    "mono_is_curved_v" : True,
    "ana_is_curved_v"  : True,
    "mono_is_optimally_curved_v" : False,
    "ana_is_optimally_curved_v"  : False,
    "mono_curv_v_formula" : mono_curv_v_formula,
    "ana_curv_v_formula" : ana_curv_v_formula,

    # guide before monochromator
    "use_guide"   : False,
    "guide_div_h" : 15. * helpers.min2rad,
    "guide_div_v" : 15. * helpers.min2rad,

    # horizontal mosaics
    "mono_mosaic"   : 27. * helpers.min2rad,
    "sample_mosaic" : 30. * helpers.min2rad,
    "ana_mosaic"    : 27. * helpers.min2rad,

    # vertical mosaics
    "mono_mosaic_v"   : 27. * helpers.min2rad,
    "sample_mosaic_v" : 30. * helpers.min2rad,
    "ana_mosaic_v"    : 27. * helpers.min2rad,

    # calculate R0 factor (not needed if only the ellipses are to be plotted)
    "calc_R0" : True,

    # crystal reflectivities; TODO, so far always 1
    "dmono_refl" : 1.,
    "dana_effic" : 1.,

    # off-center scattering
    # WARNING: while this is calculated, it is not yet considered in the ellipse plots
    "pos_x" : 0. * helpers.cm2A,
    "pos_y" : 0. * helpers.cm2A,
    "pos_z" : 0. * helpers.cm2A,

    # vertical scattering in kf, keep "False" for normal TAS
    "kf_vert" : False,
}


#
# pre-defined parameters for IN20 with FlatCone
#
params_fc = {
    # options
    "verbose" : True,

    # resolution method, "eck", "pop", or "cn"
    "reso_method" : "eck",

    # scattering triangle
    "ki" : 2.981,
    "kf" : 2.981,
    "E"  : helpers.get_E(2.981, 2.981),
    "Q"  : 2.,

    # d spacings
    "mono_xtal_d" : 3.355,   # PG
    "ana_xtal_d"  : 3.1355,  # Si

     # scattering senses
    "mono_sense"   : -1.,
    "sample_sense" :  1.,
    "ana_sense"    : -1.,
    "mirror_Qperp" : False,

    # distances
    "dist_vsrc_mono"   : 702.5 * helpers.cm2A,
    "dist_hsrc_mono"   : 227.5 * helpers.cm2A,
    "dist_mono_sample" : 227.5 * helpers.cm2A,
    "dist_sample_ana"  : 76.3  * helpers.cm2A,
    "dist_ana_det"     : 23.5  * helpers.cm2A,

    # shapes
    "src_shape"    : "rectangular",  # "rectangular" or "circular"
    "sample_shape" : "cylindrical",  # "cuboid" or "cylindrical"
    "det_shape"    : "rectangular",  # "rectangular" or "circular"

    # component sizes
    "src_w"    : 2.   * helpers.cm2A,
    "src_h"    : 17.  * helpers.cm2A,
    "mono_d"   : 0.2  * helpers.cm2A,
    "mono_w"   : 22.5 * helpers.cm2A,
    "mono_h"   : 20.  * helpers.cm2A,
    "sample_d" : 1.   * helpers.cm2A,
    "sample_w" : 1.   * helpers.cm2A,
    "sample_h" : 1.   * helpers.cm2A,
    "ana_d"    : 1.   * helpers.cm2A,
    "ana_w"    : 1.5  * helpers.cm2A,
    "ana_h"    : 13.3 * helpers.cm2A,
    "det_w"    : 2.5  * helpers.cm2A,
    "det_h"    : 6.58 * helpers.cm2A,

    # horizontal collimation
    "coll_h_pre_mono"    : 9999. * helpers.min2rad,
    "coll_h_pre_sample"  : 9999. * helpers.min2rad,
    "coll_h_post_sample" : 9999. * helpers.min2rad,
    "coll_h_post_ana"    : 9999. * helpers.min2rad,

    # vertical collimation
    "coll_v_pre_mono"    : 9999. * helpers.min2rad,
    "coll_v_pre_sample"  : 9999. * helpers.min2rad,
    "coll_v_post_sample" : 9999. * helpers.min2rad,
    "coll_v_post_ana"    : 9999. * helpers.min2rad,

    # horizontal focusing
    "mono_curv_h" : 0.,
    "ana_curv_h"  : 0.,
    "mono_is_curved_h" : True,
    "ana_is_curved_h"  : False,
    "mono_is_optimally_curved_h" : False,
    "ana_is_optimally_curved_h"  : False,
    "mono_curv_h_formula" : mono_curv_h_formula,
    "ana_curv_h_formula" : None,

    # vertical focusing
    "mono_curv_v" : 0.,
    "ana_curv_v"  : 200. * helpers.cm2A,
    "mono_is_curved_v" : True,
    "ana_is_curved_v"  : True,
    "mono_is_optimally_curved_v" : False,
    "ana_is_optimally_curved_v"  : False,
    "mono_curv_v_formula" : mono_curv_v_formula,
    "ana_curv_v_formula" : None,

    # guide before monochromator
    "use_guide"   : False,
    "guide_div_h" : 15. * helpers.min2rad,
    "guide_div_v" : 15. * helpers.min2rad,

    # horizontal mosaics
    "mono_mosaic"   : 27. * helpers.min2rad,
    "sample_mosaic" : 30. * helpers.min2rad,
    "ana_mosaic"    :  1. * helpers.min2rad,

    # vertical mosaics
    "mono_mosaic_v"   : 27. * helpers.min2rad,
    "sample_mosaic_v" : 30. * helpers.min2rad,
    "ana_mosaic_v"    :  1. * helpers.min2rad,

    # calculate R0 factor (not needed if only the ellipses are to be plotted)
    "calc_R0" : True,

    # crystal reflectivities; TODO, so far always 1
    "dmono_refl" : 1.,
    "dana_effic" : 1.,

    # off-center scattering
    # WARNING: while this is calculated, it is not yet considered in the ellipse plots
    "pos_x" : 0. * helpers.cm2A,
    "pos_y" : 0. * helpers.cm2A,
    "pos_z" : 0. * helpers.cm2A,

    # vertical scattering in kf, keep "False" for normal TAS
    "kf_vert" : True,
}
