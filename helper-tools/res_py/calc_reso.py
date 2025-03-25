#!/usr/bin/env python3
#
# calculates the instrumental resolution
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019, jul-2024
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

import libs.tas as tas
import libs.helpers as helpers
import instruments.params_in20 as params_in20

# requires numpy version >= 1.10
import numpy as np
np.set_printoptions(floatmode = "fixed",  precision = 4)


# -----------------------------------------------------------------------------
#
# default resolution parameters
# NOTE: not all parameters are used by all calculation backends
#
params = {
    # options
    "verbose" : True,

    # resolution method, "eck", "eck_ext", "pop", "cn", or "vio"
    "reso_method" : "pop",

    # scattering triangle
    "ki" : 1.4,
    "kf" : 1.4,
    "E"  : tas.get_E(1.4, 1.4),
    #"Q"  : 1.777,  # length in 1/A
    "Q"  : np.array([1., 0., 0.]),  # (hkl) in rlu

    # if Q is defined as (hkl) vector, the following crystal infos are also needed:
    "sample_lattice" : np.array([ 5., 5., 5. ]),
    "sample_angles"  : np.array([ 90., 90., 90. ]) * helpers.deg2rad,
    "sample_plane_1" : np.array([ 1., 0., 0. ]),
    "sample_plane_2" : np.array([ 0., 1., 0. ]),

    # d spacings
    "mono_xtal_d" : 3.355,
    "ana_xtal_d"  : 3.355,

     # scattering senses
    "mono_sense"   : -1.,
    "sample_sense" :  1.,
    "ana_sense"    : -1.,
    "mirror_Qperp" : False,

    # distances
    "dist_vsrc_mono"   :  10. * helpers.cm2A,
    "dist_hsrc_mono"   :  10. * helpers.cm2A,
    "dist_mono_sample" : 200. * helpers.cm2A,
    "dist_sample_ana"  : 115. * helpers.cm2A,
    "dist_ana_det"     :  85. * helpers.cm2A,

    # shapes
    "src_shape"    : "rectangular",  # "rectangular" or "circular"
    "sample_shape" : "cylindrical",  # "cuboid" or "cylindrical"
    "det_shape"    : "rectangular",  # "rectangular" or "circular"

    # component sizes
    "src_w"    : 6.   * helpers.cm2A,
    "src_h"    : 12.  * helpers.cm2A,
    "mono_d"   : 0.15 * helpers.cm2A,
    "mono_w"   : 12.  * helpers.cm2A,
    "mono_h"   : 8.   * helpers.cm2A,
    "sample_d" : 1.   * helpers.cm2A,
    "sample_w" : 1.   * helpers.cm2A,
    "sample_h" : 1.   * helpers.cm2A,
    "ana_d"    : 0.3  * helpers.cm2A,
    "ana_w"    : 12.  * helpers.cm2A,
    "ana_h"    : 8.   * helpers.cm2A,
    "det_w"    : 1.5  * helpers.cm2A,
    "det_h"    : 5.   * helpers.cm2A,

    # horizontal collimation
    "coll_h_pre_mono"    : 30. * helpers.min2rad,
    "coll_h_pre_sample"  : 30. * helpers.min2rad,
    "coll_h_post_sample" : 30. * helpers.min2rad,
    "coll_h_post_ana"    : 30. * helpers.min2rad,

    # vertical collimation
    "coll_v_pre_mono"    : 30. * helpers.min2rad,
    "coll_v_pre_sample"  : 30. * helpers.min2rad,
    "coll_v_post_sample" : 30. * helpers.min2rad,
    "coll_v_post_ana"    : 30. * helpers.min2rad,

    # horizontal focusing
    "mono_curv_h" : 0.,
    "ana_curv_h"  : 0.,
    "mono_is_curved_h" : False,
    "ana_is_curved_h"  : False,
    "mono_is_optimally_curved_h" : False,
    "ana_is_optimally_curved_h"  : False,
    "mono_curv_h_formula" : None,
    "ana_curv_h_formula" : None,

    # vertical focusing
    "mono_curv_v" : 0.,
    "ana_curv_v"  : 0.,
    "mono_is_curved_v" : False,
    "ana_is_curved_v"  : False,
    "mono_is_optimally_curved_v" : False,
    "ana_is_optimally_curved_v"  : False,
    "mono_curv_v_formula" : None,
    "ana_curv_v_formula" : None,

    # guide before monochromator
    "use_guide"   : True,
    "guide_div_h" : 15. * helpers.min2rad,
    "guide_div_v" : 15. * helpers.min2rad,

    # horizontal mosaics
    "mono_mosaic"   : 45. * helpers.min2rad,
    "sample_mosaic" : 30. * helpers.min2rad,
    "ana_mosaic"    : 45. * helpers.min2rad,

    # vertical mosaics
    "mono_mosaic_v"   : 45. * helpers.min2rad,
    "sample_mosaic_v" : 30. * helpers.min2rad,
    "ana_mosaic_v"    : 45. * helpers.min2rad,

    # crystal reflectivities; TODO, so far always 1
    "dmono_refl" : 1.,
    "dana_effic" : 1.,

    # off-center scattering (for eck and eck_ext)
    # WARNING: while this is calculated, it is not yet considered in the ellipse plots
    "pos_x" : 0. * helpers.cm2A,
    "pos_y" : 0. * helpers.cm2A,
    "pos_z" : 0. * helpers.cm2A,

    # vertical scattering in kf (for eck and eck_ext), keep "False" for normal TAS
    "kf_vert" : False,

    # sample integration method (for eck_ext): "gaussian" or "analytical"
    "sample_int" : "gaussian",
}
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# get command-line arguments
# -----------------------------------------------------------------------------
import argparse
argparser = argparse.ArgumentParser(
    description = "Calculates the resolution ellipsoid of a TAS instrument.")

argparser.add_argument("-i", "--instr", default = None, type = str,
    help = "set the parameters to a pre-defined instrument (in20/in20fc)")
argparser.add_argument("--silent", action = "store_true",
    help = "disable output")
argparser.add_argument("-o", "--out_file", default = "", type = str,
    help = "output file for results")
argparser.add_argument("--plot_file", default = "", type = str,
    help = "output file for plots")
argparser.add_argument("-p", "--plot", action = "store_true",
    help = "plot results on screen")
argparser.add_argument("-m", "--reso_method", default = None, type = str,
    help = "resolution method to use (cn/pop/eck/eck_ext/vio)")
argparser.add_argument("--kf_vert", default = None, action = "store_true",
    help = "scatter vertically in the kf axis (only for Eckold-Sobolev method)")
argparser.add_argument("--ki", default = None, type = float,
    help = "incoming wavenumber in 1/A")
argparser.add_argument("--kf", default = None, type = float,
    help = "outgoing wavenumber in 1/A")
argparser.add_argument("-E", "--E", default = None, type = float,
    help = "energy transfer in meV")
argparser.add_argument("-Q", "--Q", default = None, type = float,
    help = "momentum transfer in 1/A")
argparser.add_argument("-Qh", "--Qh", default = None, type = float,
    help = "h Miller index in rlu")
argparser.add_argument("-Qk", "--Qk", default = None, type = float,
    help = "k Miller index in rlu")
argparser.add_argument("-Ql", "--Ql", default = None, type = float,
    help = "l Miller index in rlu")
argparser.add_argument("--twotheta", default = None, type = float,
    help = "sample scattering angle")
argparser.add_argument("--pos_x", default = None, type = float,
    help = "sample x position")
argparser.add_argument("-a", "--a", default = None, type = float, nargs="?",
    help = "lattice constant a")
argparser.add_argument("-b", "--b", default = None, type = float, nargs="?",
    help = "lattice constant b")
argparser.add_argument("-c", "--c", default = None, type = float, nargs="?",
    help = "lattice constant c")
argparser.add_argument("--alpha", default = 90., type = float, nargs="?",
    help = "lattice angle alpha")
argparser.add_argument("--beta", default = 90., type = float, nargs="?",
    help = "lattice angle beta")
argparser.add_argument("--gamma", default = 90., type = float, nargs="?",
    help = "lattice angle gamma")
argparser.add_argument("--ax", default = None, type = float, nargs="?",
    help = "x component of first orientation vector")
argparser.add_argument("--ay", default = None, type = float, nargs="?",
    help = "y component of first orientation vector")
argparser.add_argument("--az", default = None, type = float, nargs="?",
    help = "z component of first orientation vector")
argparser.add_argument("--bx", default = None, type = float, nargs="?",
    help = "x component of second orientation vector")
argparser.add_argument("--by", default = None, type = float, nargs="?",
    help = "y component of second orientation vector")
argparser.add_argument("--bz", default = None, type = float, nargs="?",
    help = "z component of second orientation vector")
argparser.add_argument("--pos_y", default = None, type = float,
    help = "sample y position")
argparser.add_argument("--pos_z", default = None, type = float,
    help = "sample z position")
argparser.add_argument("--mono_sense", default = None, type = float,
    help = "monochromator scattering sense")
argparser.add_argument("--sample_sense", default = None, type = float,
    help = "sample scattering sense")
argparser.add_argument("--ana_sense", default = None, type = float,
    help = "analyser scattering sense")
argparser.add_argument("--mono_curv_v", default = None, type = float,
    help = "monochromator vertical curvature")
argparser.add_argument("--mono_curv_h", default = None, type = float,
    help = "monochromator horizontal curvature")
argparser.add_argument("--ana_curv_v", default = None, type = float,
    help = "analyser vertical curvature")
argparser.add_argument("--ana_curv_h", default = None, type = float,
    help = "analyser horizontal curvature")
argparser.add_argument("--mono_curv_v_opt", action = "store_true",
    help = "set optimal monochromator vertical curvature")
argparser.add_argument("--mono_curv_h_opt", action = "store_true",
    help = "set optimal monochromator horizontal curvature")
argparser.add_argument("--ana_curv_v_opt", action = "store_true",
    help = "set optimal analyser vertical curvature")
argparser.add_argument("--ana_curv_h_opt", action = "store_true",
    help = "set optimal analyser horizontal curvature")

parsedargs = argparser.parse_args()

if parsedargs.instr != None:
    if parsedargs.instr == "in20":
        print("Loaded IN20 default parameters.\n")
        params = params_in20.params
    elif parsedargs.instr == "in20fc":
        print("Loaded IN20/Flatcone default parameters.\n")
        params = params_in20.params_fc

# get parsed command-line arguments
out_file = parsedargs.out_file
show_plots = parsedargs.plot
plot_file = parsedargs.plot_file
params["verbose"] = not parsedargs.silent

if parsedargs.reso_method != None:
    params["reso_method"] = parsedargs.reso_method

if parsedargs.kf_vert != None:
    params["kf_vert"] = parsedargs.kf_vert
    # TODO: flip y and z parameters if this is enabled
if parsedargs.ki != None:
    params["ki"] = parsedargs.ki
if parsedargs.kf != None:
    params["kf"] = parsedargs.kf

if parsedargs.Qh != None and parsedargs.Qk != None and parsedargs.Ql != None:
    # Q is given as (hkl) vector
    params["Q"] = np.array([ parsedargs.Qh, parsedargs.Qk, parsedargs.Ql ])
elif parsedargs.Q != None:
    # Q is given in 1/A
    params["Q"] = parsedargs.Q
elif parsedargs.twotheta != None:
    # Q is defined by the scattering angle
    params["twotheta"] = parsedargs.twotheta * helpers.deg2rad
    params["Q"] = tas.get_Q(params["ki"], params["kf"], params["twotheta"])
#else:
#    params["twotheta"] = tas.get_scattering_angle(params["ki"], params["kf"], params["Q"])

if parsedargs.ax != None and parsedargs.ay != None and parsedargs.az != None:
    params["sample_plane_1"] = np.array([ parsedargs.ax, parsedargs.ay, parsedargs.az ])
if parsedargs.bx != None and parsedargs.by != None and parsedargs.bz != None:
    params["sample_plane_2"] = np.array([ parsedargs.bx, parsedargs.by, parsedargs.bz ])
if parsedargs.a != None and parsedargs.b != None and parsedargs.c != None:
    params["sample_lattice"] = np.array([ parsedargs.a, parsedargs.b, parsedargs.c ])
if parsedargs.alpha != None and parsedargs.beta != None and parsedargs.gamma != None:
    params["sample_angles"] = np.array([ parsedargs.alpha, parsedargs.beta, parsedargs.gamma ]) * helpers.deg2rad

if parsedargs.pos_x != None:
    params["pos_x"] = parsedargs.pos_x * helpers.cm2A
if parsedargs.pos_y != None:
    params["pos_y"] = parsedargs.pos_y * helpers.cm2A
if parsedargs.pos_z != None:
    params["pos_z"] = parsedargs.pos_z * helpers.cm2A

if parsedargs.mono_sense != None:
    params["mono_sense"] = parsedargs.mono_sense
if parsedargs.sample_sense != None:
    params["sample_sense"] = parsedargs.sample_sense
if parsedargs.ana_sense != None:
    params["ana_sense"] = parsedargs.ana_sense

# ensure that the senses are either +1 or -1
if params["mono_sense"] > 0.:
    params["mono_sense"] = 1.
else:
    params["mono_sense"] = -1.
if params["sample_sense"] > 0.:
    params["sample_sense"] = 1.
else:
    params["sample_sense"] = -1.
if params["ana_sense"] > 0.:
    params["ana_sense"] = 1.
else:
    params["ana_sense"] = -1.

# set instrument position
if parsedargs.E != None:
    # calculate ki from E and kf
    params["E"] = parsedargs.E
    params["ki"] = tas.get_ki(params["kf"], params["E"])
else:
    # calculate E from ki and kf
    params["E"] = tas.get_E(params["ki"], params["kf"])

# set fixed curvatures
if parsedargs.mono_curv_v != None:
    params["mono_curv_v"] = parsedargs.mono_curv_v
    params["mono_is_curved_v"] = True
if parsedargs.mono_curv_h != None:
    params["mono_curv_h"] = parsedargs.mono_curv_h
    params["mono_is_curved_h"] = True
if parsedargs.ana_curv_v != None:
    params["ana_curv_v"] = parsedargs.ana_curv_v
    params["ana_is_curved_v"] = True
if parsedargs.ana_curv_h != None:
    params["ana_curv_h"] = parsedargs.ana_curv_h
    params["ana_is_curved_h"] = True

# set optimal curvatures
if parsedargs.mono_curv_v_opt:
    params["mono_is_optimally_curved_v"] = True
    params["mono_is_curved_v"] = True
if parsedargs.mono_curv_h_opt:
    params["mono_is_optimally_curved_h"] = True
    params["mono_is_curved_h"] = True
if parsedargs.ana_curv_v_opt:
    params["ana_is_optimally_curved_v"] = True
    params["ana_is_curved_v"] = True
if parsedargs.ana_curv_h_opt:
    params["ana_is_optimally_curved_h"] = True
    params["ana_is_curved_h"] = True


def log(msg):
    if params["verbose"]:
        print(msg)


# is Q given as (hkl) vector?
if isinstance(params["Q"], type([])) or isinstance(params["Q"], type(np.array([]))):
    log("Sample lattice: a = %g A, b = %g A, c = %g A, alpha = %g deg, beta = %g deg, gamma = %g deg." %
        (params["sample_lattice"][0], params["sample_lattice"][1], params["sample_lattice"][2],
         params["sample_angles"][0] * helpers.rad2deg,
         params["sample_angles"][1] * helpers.rad2deg,
         params["sample_angles"][2] * helpers.rad2deg))
    log("Scattering plane: orient1: %s rlu, orient2 %s rlu." % (params["sample_plane_1"], params["sample_plane_2"]))
    log("ki = %g / A, kf = %g / A, E = %g meV, Q = %s rlu." %
        (params["ki"], params["kf"], params["E"], params["Q"]))
else:
    log("ki = %g / A, kf = %g / A, E = %g meV, Q = %g / A." %
        (params["ki"], params["kf"], params["E"], params["Q"]))
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# calculate resolution ellipsoid using the given backend
# -----------------------------------------------------------------------------
import libs.reso as reso
import algos.pop as pop
import algos.eck as eck
import algos.eck_ext as eck_ext
import algos.vio as vio

if params["reso_method"] == "eck":
    log("\nCalculating using Eckold-Sobolev method. Scattering %s in kf." %
        ("vertically" if params["kf_vert"] else "horizontally"))
    res = eck.calc(params)
elif params["reso_method"] == "eck_ext":
    log("\nCalculating using extended Eckold-Sobolev method. Scattering %s in kf." %
        ("vertically" if params["kf_vert"] else "horizontally"))
    res = eck_ext.calc(params)
elif params["reso_method"] == "pop":
    log("\nCalculating using Popovici method.")
    res = pop.calc(params, False)
elif params["reso_method"] == "cn":
    log("\nCalculating using Cooper-Nathans method.")
    res = pop.calc(params, True)
elif params["reso_method"] == "vio":
    log("\nCalculating using Violini method.")
    res = vio.calc(params)
else:
    raise ValueError("ResPy: Invalid resolution calculation method selected.")

if not res["ok"]:
    print("RESOLUTION CALCULATION FAILED!")
    exit(-1)


log("R0 = %g, Vol = %g" % (res["r0"], res["res_vol"]))
log("Resolution matrix (Q_para [1/A], Q_perp [1/A], Q_up [1/A], E [meV]):\n%s" % res["reso"])
log("Resolution vector (Q_para [1/A], Q_perp [1/A], Q_up [1/A], E [meV]):\n%s" % res["reso_v"])
log("Resolution scalar: %g" % res["reso_s"])
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# describe and plot ellipses
# -----------------------------------------------------------------------------
ellipses = reso.calc_ellipses(res["reso"], params["verbose"])


if show_plots or plot_file != "":
    # plot ellipses
    reso.plot_ellipses(ellipses, verbose = params["verbose"], plot_results = show_plots, plot_file = plot_file)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# write the results to a file
# -----------------------------------------------------------------------------
if out_file != "":
    with open(out_file, "w") as file:
        log("\nWriting results to \"%s\"." % out_file)

        print("[scattering_triangle]", file = file)
        print("\tQ = %g" % res["Q_avg"][0], file = file)
        print("\tE = %g" % res["Q_avg"][3], file = file)
        print("\tki = %g" % res["ki"], file = file)
        print("\tkf = %g" % res["kf"], file = file)
        print("\tQ_ki = %g" % (res["Q_ki"] * helpers.rad2deg), file = file)
        print("\tQ_kf = %g" % (res["Q_kf"] * helpers.rad2deg), file = file)
        print("\t2theta = %g" % (res["twotheta"] * helpers.rad2deg), file = file)
        print("\ttheta_m = %g" % (res["theta_m"] * helpers.rad2deg), file = file)
        print("\ttheta_a = %g" % (res["theta_a"] * helpers.rad2deg), file = file)
        print("", file = file)


        print("\n[resolution]", file = file)

        print("\tR0 = %g" % res["r0"], file = file)
        print("\tV = %g" % res["res_vol"], file = file)
        print("", file = file)

        for i in range(0, 4):
            for j in range(0, 4):
                print("\treso_%d%d = %g" % (i, j, res["reso"][i][j]), file = file)
        print("", file = file)

        for i in range(0, 4):
            print("\treso_v%d = %g" % (i, res["reso_v"][i]), file = file)
        print("", file = file)

        print("\treso_s = %g" % res["reso_s"], file = file)
        print("", file = file)


        print("\n[principal_axes]", file = file)

        for i in range(0, 4):
            for j in range(0, 4):
                print("\teigensys_%d%d = %g" % (i, j, ellipses["rot"][i][j]), file = file)
        print("", file = file)

        for i in range(0, 4):
            print("\teigenval_%d = %g" % (i, ellipses["evals"][i]), file = file)
        print("", file = file)

        for i in range(0, 4):
            print("\tfwhm_%d = %g" % (i, ellipses["fwhms"][i]), file = file)
        print("", file = file)

        for i in range(0, 4):
            print("\tfwhm_coh_%d = %g" % (i, ellipses["fwhms_coh"][i]), file = file)
        print("", file = file)

        for i in range(0, 4):
            print("\tfwhm_inc_%d = %g" % (i, ellipses["fwhms_inc"][i]), file = file)
        print("", file = file)


        print("\n[parameters]", file = file)

        for [ key, value ] in params.items():
            if isinstance(value, float):
                print("\t%s = %g" % (key, value), file = file)
            elif isinstance(value, int):
                print("\t%s = %d" % (key, value), file = file)
            elif callable(value):
                print("\t%s = <func>" % key, file = file)
            else:
                print("\t%s = \"%s\"" % (key, value), file = file)
# -----------------------------------------------------------------------------
