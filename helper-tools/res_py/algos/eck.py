#
# implementation of the eckold-sobolev algo
#
# @author Tobias Weber <tweber@ill.fr>
# @date feb-2015, oct-2019
# @license see 'LICENSE' file
#
# @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
# @desc for alternate R0 normalisation: [mit84] P. W. Mitchell, R. A. Cowley and S. A. Higgins, Acta Cryst. Sec A, 40(2), 152-160 (1984), doi: 10.1107/S0108767384000325
# @desc for vertical scattering modification: [eck20] G. Eckold, personal communication, 2020.
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

import libs.reso as reso
import libs.tas as tas
import libs.helpers as helpers

import numpy as np
import numpy.linalg as la


#
# mono (and ana) resolution calculation
#
def get_mono_vals(src_w, src_h, mono_w, mono_h,
    dist_vsrc_mono, dist_hsrc_mono,
    dist_mono_sample,
    ki, thetam,
    coll_h_pre_mono, coll_h_pre_sample,
    coll_v_pre_mono, coll_v_pre_sample,
    mono_mosaic, mono_mosaic_v,
    inv_mono_curv_h, inv_mono_curv_v,
    sample_pos, refl):

    # A matrix: equ. 26 in [eck14]
    A = np.identity(3)

    A_t0 = 1. / mono_mosaic
    A_tx = inv_mono_curv_h*dist_mono_sample / np.abs(np.sin(thetam))
    A_t1 = A_t0 * A_tx

    # typo in paper?
    A[0, 0] = 0.5*helpers.sig2fwhm**2. / ki**2. * np.tan(thetam)**2. * \
        ( (2./coll_h_pre_mono)**2. + (2*dist_hsrc_mono/src_w)**2. + A_t0**2. )

    A[0, 1] = A[1, 0] = 0.5*helpers.sig2fwhm**2. / ki**2. * np.tan(thetam) * \
        ( + 2.*(1./coll_h_pre_mono)**2. + \
            2.*dist_hsrc_mono*(dist_hsrc_mono-dist_mono_sample)/src_w**2. + \
            A_t0**2. - A_t0*A_t1 )

    A[1, 1] = 0.5*helpers.sig2fwhm**2. / ki**2. * \
        ( 1./coll_h_pre_mono**2. + 1./coll_h_pre_sample**2. \
            + ((dist_hsrc_mono-dist_mono_sample)/src_w)**2. \
            + (dist_mono_sample/(mono_w*np.abs(np.sin(thetam))))**2. \
            + (A_t0 - A_t1)**2. )


    # Av matrix: equ. 38 in [eck14]
    # some typos in paper leading to the (false) result of a better Qz resolution when focusing
    # => trying to match terms in Av with corresponding terms in A
    # corresponding pre-mono terms commented out in Av, as they are not considered there
    Av = np.identity(2)

    Av_t0 = 0.5 / (mono_mosaic_v * np.abs(np.sin(thetam)))
    Av_t1 = inv_mono_curv_v * dist_mono_sample / mono_mosaic_v

    Av[0, 0] = 0.5*helpers.sig2fwhm**2. / ki**2. * \
        ( 1./coll_v_pre_sample**2. + (dist_mono_sample/src_h)**2. + (dist_mono_sample/mono_h)**2. + \
        (Av_t0 - Av_t1)**2. )     # typo/missing in paper?

    Av[0, 1] = Av[1, 0] = 0.5*helpers.sig2fwhm**2. / ki**2. * \
        ( dist_vsrc_mono*dist_mono_sample/src_h**2. - Av_t0**2. + Av_t0*Av_t1 )

    Av[1, 1] = 0.5*helpers.sig2fwhm**2. / ki**2. * \
        ( (1./coll_v_pre_mono)**2. + (dist_vsrc_mono/src_h)**2. + Av_t0**2. )


    # B vector: equ. 27 in [eck14]
    B = np.array([0., 0., 0.])
    B_t0 = inv_mono_curv_h / (mono_mosaic**2. * np.abs(np.sin(thetam)))

    B[0] = helpers.sig2fwhm**2. * sample_pos[1] / ki * np.tan(thetam) * \
        ( 2.*dist_hsrc_mono / src_w**2. + B_t0 )

    B[1] = helpers.sig2fwhm**2. * sample_pos[1] / ki * \
        ( - dist_mono_sample / (mono_w*np.abs(np.sin(thetam)))**2. + \
        B_t0 - B_t0 * A_tx + \
        (dist_hsrc_mono-dist_mono_sample) / src_w**2. )


    # Bv vector: equ. 39 in [eck14]
    Bv = np.array([0., 0.])

    Bv_t0 = inv_mono_curv_v / mono_mosaic_v**2

    # typo in paper?
    Bv[0] = (-1.) * helpers.sig2fwhm**2. * sample_pos[2] / ki * \
        ( dist_mono_sample / mono_h**2. + \
            dist_mono_sample / src_h**2. + \
            Bv_t0 * inv_mono_curv_v*dist_mono_sample - \
            0.5*Bv_t0 / np.abs(np.sin(thetam)) )

    # typo in paper?
    Bv[1] = (-1.) * helpers.sig2fwhm**2. * sample_pos[2] / ki * \
        ( dist_vsrc_mono / (src_h*src_h) + 0.5*Bv_t0/np.abs(np.sin(thetam)) )


    # C scalar: equ. 28 in [eck14]
    C = 0.5*helpers.sig2fwhm**2. * sample_pos[1]**2. * \
        ( 1./src_w**2. + (1./(mono_w*np.abs(np.sin(thetam))))**2. + \
            (inv_mono_curv_h/(mono_mosaic * np.abs(np.sin(thetam))))**2. )

    # Cv scalar: equ. 40 in [eck14]
    Cv = 0.5*helpers.sig2fwhm**2. * sample_pos[2]**2. * \
        ( 1./src_h**2. + 1./mono_h**2. + (inv_mono_curv_v/mono_mosaic_v)**2. )


    # z components, [eck14], equ. 42
    A[2, 2] = Av[0, 0] - Av[0, 1]**2./Av[1, 1]
    B[2] = Bv[0] - Bv[1]*Av[0, 1]/Av[1, 1]
    # typo in paper? (thanks to F. Bourdarot for pointing this out)
    D = Cv - (0.5*Bv[1])**2./Av[1, 1]


    # [eck14], equ. 54
    therefl = refl * np.sqrt(np.pi / Av[1, 1])  # typo in paper?

    return [ A, B, C, D, therefl ]



#
# Eckold algorithm combining the mono and ana resolutions
#
def calc(param):
    helpers.calc_triangle(param)

    ki = param["ki"]
    kf = param["kf"]
    E = param["E"]
    Q = param["Q"]
    twotheta = param["twotheta"]
    thetam = param["thetam"]
    thetaa = param["thetaa"]
    Q_ki = param["Q_ki"]
    Q_kf = param["Q_kf"]
    sample_pos = np.array([ param["pos_x"],  param["pos_y"], param["pos_z"] ])

    # --------------------------------------------------------------------
    # mono/ana focus
    # use fixed values
    mono_curv_h = param["mono_curv_h"]
    mono_curv_v = param["mono_curv_v"]
    ana_curv_h = param["ana_curv_h"]
    ana_curv_v = param["ana_curv_v"]

    # use a user-defined curvature formula, if given
    if param["mono_curv_h_formula"] != None:
        mono_curv_h = param["mono_curv_h_formula"](param) * helpers.cm2A
    if param["mono_curv_v_formula"] != None:
        mono_curv_v = param["mono_curv_v_formula"](param) * helpers.cm2A
    if param["ana_curv_h_formula"] != None:
        ana_curv_h = param["ana_curv_h_formula"](param) * helpers.cm2A
    if param["ana_curv_v_formula"] != None:
        ana_curv_v = param["ana_curv_v_formula"](param) * helpers.cm2A

    if param["mono_is_optimally_curved_h"]:
        mono_curv_h = helpers.foc_curv(param["dist_hsrc_mono"], \
            param["dist_mono_sample"], np.abs(2.*thetam), False)
    if param["mono_is_optimally_curved_v"]:
        mono_curv_v = helpers.foc_curv(param["dist_vsrc_mono"], \
            param["dist_mono_sample"], np.abs(2.*thetam), True)
    if param["ana_is_optimally_curved_h"]:
        ana_curv_h = helpers.foc_curv(param["dist_sample_ana"], \
            param["dist_ana_det"], np.abs(2.*thetaa), False)
    if param["ana_is_optimally_curved_v"]:
        ana_curv_v = helpers.foc_curv(param["dist_sample_ana"], \
            param["dist_ana_det"], np.abs(2.*thetaa), True)

    if param["verbose"]:
        print("Mono curvature radius: vertical: %g cm, horizontal: %g cm." %
                (mono_curv_v/helpers.cm2A, mono_curv_h/helpers.cm2A))
        print("Ana curvature radius: vertical: %g cm, horizontal: %g cm.\n" %
                (ana_curv_v/helpers.cm2A, ana_curv_h/helpers.cm2A))

    inv_mono_curv_h = 0.
    inv_mono_curv_v = 0.
    inv_ana_curv_h = 0.
    inv_ana_curv_v = 0.

    if param["mono_is_curved_h"]:
        inv_mono_curv_h = 1./mono_curv_h
    if param["mono_is_curved_v"]:
        inv_mono_curv_v = 1./mono_curv_v
    if param["ana_is_curved_h"]:
        inv_ana_curv_h = 1./ana_curv_h
    if param["ana_is_curved_v"]:
        inv_ana_curv_v = 1./ana_curv_v
    # --------------------------------------------------------------------


    lam = tas.k_to_lam(ki)

    coll_h_pre_mono = param["coll_h_pre_mono"]
    coll_v_pre_mono = param["coll_v_pre_mono"]

    if param["use_guide"]:
        coll_h_pre_mono = lam*param["guide_div_h"]
        coll_v_pre_mono = lam*param["guide_div_v"]


    # dict with results
    res = {}

    res["Q_avg"] = np.array([ Q, 0., 0., E ])
    res["ki"] = ki
    res["kf"] = kf
    res["Q_ki"] = Q_ki
    res["Q_kf"] = Q_kf
    res["twotheta"] = twotheta
    res["theta_m"] = thetam
    res["theta_a"] = thetaa

    # -------------------------------------------------------------------------

    # - if the instruments works in kf=const mode and the scans are counted for
    #   or normalised to monitor counts no ki^3 or kf^3 factor is needed.
    # - if the instrument works in ki=const mode the kf^3 factor is needed.

    # TODO
    #tupScFact = get_scatter_factors(param.flags, param.thetam, param.ki, param.thetaa, param.kf)
    tupScFact = [1., 1., 1.]

    dmono_refl = param["dmono_refl"] * tupScFact[0]
    dana_effic = param["dana_effic"] * tupScFact[1]
    dxsec = tupScFact[2]
    #if param.mono_refl_curve:
    #    dmono_refl *= (*param.mono_refl_curve)(param.ki)
    #if param.ana_effic_curve:
    #    dana_effic *= (*param.ana_effic_curve)(param.kf)


    #--------------------------------------------------------------------------
    # mono part
    [A, B, C, D, dReflM] = get_mono_vals(
        param["src_w"], param["src_h"],
        param["mono_w"], param["mono_h"],
        param["dist_vsrc_mono"], param["dist_hsrc_mono"],
        param["dist_mono_sample"],
        ki, thetam,
        coll_h_pre_mono, param["coll_h_pre_sample"],
        coll_v_pre_mono, param["coll_v_pre_sample"],
        param["mono_mosaic"], param["mono_mosaic_v"],
        inv_mono_curv_h, inv_mono_curv_v,
        sample_pos, dmono_refl)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # ana part, equ. 43 in [eck14]
    sample_pos_kf = np.dot(helpers.rotation_matrix_nd(-twotheta, 3), sample_pos)

    # vertical scattering in kf axis, formula from [eck20]
    if param["kf_vert"]:
        T_vert = helpers.rotation_matrix_3d_x(-np.pi / 2.)
        sample_pos_kf = np.dot(T_vert, sample_pos_kf)

    [E, F, G, H, dReflA] = get_mono_vals(
        param["det_w"], param["det_h"],
        param["ana_w"], param["ana_h"],
        param["dist_ana_det"], param["dist_ana_det"],
        param["dist_sample_ana"],
        kf, -thetaa,
        param["coll_h_post_ana"], param["coll_h_post_sample"],
        param["coll_v_post_ana"], param["coll_v_post_sample"],
        param["ana_mosaic"], param["ana_mosaic_v"],
        inv_ana_curv_h, inv_ana_curv_v,
        sample_pos_kf, dana_effic)

    # vertical scattering in kf axis, formula from [eck20]
    if param["kf_vert"]:
        E = np.dot(np.dot(np.transpose(T_vert), E), T_vert)
        F = np.dot(np.transpose(T_vert), F)  # possible typo in [eck20]?
    #--------------------------------------------------------------------------


    # equ. 4 & equ. 53 in [eck14]
    dE = (ki**2. - kf**2.) / (2.*Q**2.)
    dEi = 0.5 + dE
    dEf = 0.5 - dE
    kperp = np.sqrt(ki**2. - (Q*dEi)**2.)
    kperp *= param["sample_sense"]


    # trafo, equ. 52 in [eck14]
    T = np.identity(6)
    T[0,3] = T[1,4] = T[2,5] = -1.
    T[3,0] = 2.*tas.k2_to_E * Q*dEi
    T[3,3] = 2.*tas.k2_to_E * Q*dEf
    T[3,1] = 2.*tas.k2_to_E * kperp
    T[3,4] = -2.*tas.k2_to_E * kperp
    T[4,1] = T[5,2] = dEf
    T[4,4] = T[5,5] = dEi

    Tinv = la.inv(T)


    # equ. 54 in [eck14]
    Dalph_i = helpers.rotation_matrix_nd(-Q_ki, 3)
    Dalph_f = helpers.rotation_matrix_nd(-Q_kf, 3)

    matAE = np.zeros((6,6))
    matAE[0:3, 0:3] = np.dot(np.dot(np.transpose(Dalph_i), A), Dalph_i)
    matAE[3:6, 3:6] = np.dot(np.dot(np.transpose(Dalph_f), E), Dalph_f)

    # U1 matrix
    # typo in paper in quadric trafo in equ. 54 (top)?
    U1 = np.dot(np.dot(np.transpose(Tinv), matAE), Tinv)

    # V1 vector
    vecBF = np.zeros(6)
    vecBF[0:3] = np.dot(np.transpose(Dalph_i), B)
    vecBF[3:6] = np.dot(np.transpose(Dalph_f), F)
    V1 = np.dot(vecBF, Tinv)


    # --------------------------------------------------------------------------
    # integrate last 2 vars -> equs. 57 & 58 in [eck14]

    U2 = reso.quadric_proj(U1, 5)
    # careful: factor -0.5*... missing in U matrix compared to normal gaussian!
    U = 2. * reso.quadric_proj(U2, 4)

    V2 = reso.quadric_proj_vec(V1, U1, 5)
    V = reso.quadric_proj_vec(V2, U2, 4)

    W = C + D + G + H

    # squares in Vs missing in paper? (thanks to F. Bourdarot for pointing this out)
    W -= (0.5*V1[5])**2. / U1[5, 5] + (0.5*V2[4])**2. / U2[4, 4]

    R0 = dReflM*dReflA * np.pi * np.sqrt(1. / np.abs(U1[5, 5] * U2[4, 4]))
    # --------------------------------------------------------------------------


    # --------------------------------------------------------------------------
    # include sample mosaic, see cn.cpp
    # add horizontal sample mosaic
    mos_Q_sq = (param["sample_mosaic"] * Q)**2.
    vec1 = U[:, 1] / helpers.sig2fwhm**2.
    U -= helpers.sig2fwhm**2. * np.outer(vec1, vec1) / \
        (1./mos_Q_sq + U[1,1]/helpers.sig2fwhm**2.)

    # add vertical sample mosaic
    mos_v_Q_sq = (param["sample_mosaic_v"] * Q)**2.
    vec2 = U[:, 2] / helpers.sig2fwhm**2.
    U -= helpers.sig2fwhm**2. * np.outer(vec2, vec2) / \
        (1./mos_v_Q_sq + U[2,2]/helpers.sig2fwhm**2.)
    # --------------------------------------------------------------------------


    # quadratic part of quadric (matrix U)
    R = U
    # linear (vector V) and constant (scalar W) part of quadric
    res["reso_v"] = V
    res["reso_s"] = W


    if param["mirror_Qperp"] and param["sample_sense"] < 0.:
        # mirror Q_perp
        matMirror = helpers.mirror_matrix(len(R), 1)
        R = np.dot(np.dot(np.transpose(matMirror), R), matMirror)
        res["reso_v"][1] = -res["reso_v"][1]


    # prefactor and volume
    res_vol = reso.ellipsoid_volume(R)

    # missing volume prefactor to normalise gaussian,
    # cf. equ. 56 in [eck14] to  equ. 1 in [pop75] and equ. A.57 in [mit84]
    R0 *= res_vol * np.pi * 3.
    R0 *= np.exp(-W)
    R0 *= dxsec

    res["reso"] = R
    res["r0"] = R0
    res["res_vol"] = res_vol

    # Bragg widths
    res["coherent_fwhms"] = reso.calc_coh_fwhms(res["reso"])
    res["ok"] = True

    if np.isnan(res["r0"]) or np.isinf(res["r0"]) or np.isnan(res["reso"].any()) or np.isinf(res["reso"].any()):
        res["ok"] = False

    return res
