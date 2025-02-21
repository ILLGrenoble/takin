#!/usr/bin/env python
#
# resolution ellipsoid calculations
#
# @author Tobias Weber <tweber@ill.fr>
# @date mar-2019
# @license GPLv2
#
# @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
# @desc see covariance calculations: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/cov.py
# @desc see also: https://github.com/McStasMcXtrace/McCode/blob/master/tools/Legacy-Perl/mcresplot.pl
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
import numpy.linalg as la
import helpers

g_eps = 1e-8


#
# volume of the ellipsoid
#
def ellipsoid_volume(mat):
    det = np.abs(la.det(mat))
    return 4./3. * np.pi * np.sqrt(1./det)


#
# projects along one axis of the quadric
# (see [eck14], equ. 57)
#
def quadric_proj(quadric, idx):
    if np.abs(quadric[idx, idx]) < g_eps:
        return np.delete(np.delete(quadric, idx, axis = 0), idx, axis = 1)

    # row/column along which to perform the orthogonal projection
    vec = 0.5 * (quadric[idx, :] + quadric[:, idx])   # symmetrise if not symmetric
    vec /= np.sqrt(quadric[idx, idx])                 # normalise to indexed component
    proj_op = np.outer(vec, vec)                      # projection operator
    ortho_proj = quadric - proj_op                    # projected quadric

    # comparing with simple projection
    #rank = len(quadric)
    #vec /= np.sqrt(np.dot(vec, vec))
    #proj_op = np.outer(vec, vec)
    #ortho_proj = np.dot((np.identity(rank) - proj_op), quadric)

    # remove row/column that was projected out
    #print("\nProjected row/column %d:\n%s\n->\n%s.\n" % (idx, str(quadric), str(ortho_proj)))
    return np.delete(np.delete(ortho_proj, idx, axis = 0), idx, axis = 1)


#
# projects along one axis of the quadric
# (see [eck14], equ. 57)
#
def quadric_proj_mat(mat, quadric, idx):
    if np.abs(quadric[idx, idx]) < g_eps:
        return np.delete(np.delete(mat, idx, axis = 0), idx, axis = 1)

    # row/column along which to perform the orthogonal projection
    vec = 0.5 * (quadric[idx, :] + quadric[:, idx])   # symmetrise if not symmetric
    vec /= quadric[idx, idx]                          # normalise to indexed component
    proj_op = np.outer(vec, mat[idx, :])              # projection operator
    ortho_proj = mat - proj_op                        # projected matrix

    # remove row/column that was projected out
    #return np.delete(np.delete(ortho_proj, idx, axis = 0), idx, axis = 1)
    return np.delete(ortho_proj, idx, axis = 0)


#
# projects linear part of the quadric
# (see [eck14], equ. 57)
#
def quadric_proj_vec(vec, quadric, idx):
    _col = quadric[:,idx]
    col = np.delete(_col, idx, axis = 0)
    if np.abs(_col[idx]) < g_eps:
        return col

    v = np.delete(vec, idx, axis = 0)
    v = v - col*vec[idx]/_col[idx]

    return v


#
# coherent fwhm widths
#
def calc_coh_fwhms(reso):
    vecFwhms = []
    for i in range(len(reso)):
        vecFwhms.append(helpers.sig2fwhm / np.sqrt(reso[i,i]))

    return np.array(vecFwhms)


#
# incoherent fwhm width
#
def calc_incoh_fwhms(reso):
    Qres_proj_Qpara = quadric_proj(reso, 3)
    Qres_proj_Qpara = quadric_proj(Qres_proj_Qpara, 2)
    Qres_proj_Qpara = quadric_proj(Qres_proj_Qpara, 1)

    Qres_proj_Qperp = quadric_proj(reso, 3)
    Qres_proj_Qperp = quadric_proj(Qres_proj_Qperp, 2)
    Qres_proj_Qperp = quadric_proj(Qres_proj_Qperp, 0)

    Qres_proj_Qup = quadric_proj(reso, 3)
    Qres_proj_Qup = quadric_proj(Qres_proj_Qup, 1)
    Qres_proj_Qup = quadric_proj(Qres_proj_Qup, 0)

    Qres_proj_E = quadric_proj(reso, 2)
    Qres_proj_E = quadric_proj(Qres_proj_E, 1)
    Qres_proj_E = quadric_proj(Qres_proj_E, 0)

    return np.array([
        1./np.sqrt(np.abs(Qres_proj_Qpara[0,0])) * helpers.sig2fwhm,
        1./np.sqrt(np.abs(Qres_proj_Qperp[0,0])) * helpers.sig2fwhm,
        1./np.sqrt(np.abs(Qres_proj_Qup[0,0])) * helpers.sig2fwhm,
        1./np.sqrt(np.abs(Qres_proj_E[0,0])) * helpers.sig2fwhm ])


#
# calculates the characteristics of a given ellipse by principal axis trafo
#
def descr_ellipse(quadric):
    [ evals, evecs ] = la.eig(quadric)

    fwhms = 1./np.sqrt(np.abs(evals)) * helpers.sig2fwhm

    angles = np.array([])
    if len(quadric) == 2:
        angles = np.array([ np.arctan2(evecs[1][0], evecs[0][0]) ])

    return [fwhms, angles*helpers.rad2deg, evecs, evals]


#
# describes the ellipsoid by a principal axis trafo and by 2d cuts
#
def calc_ellipses(Qres_Q, verbose = True):
    # 4d ellipsoid
    [fwhms, angles, rot, evals] = descr_ellipse(Qres_Q)
    fwhms_coh = calc_coh_fwhms(Qres_Q)
    fwhms_inc = calc_incoh_fwhms(Qres_Q)

    if verbose:
        print()
        print("Eigenvalues: %s" % evals)
        print("Eigensystem (Q_para [1/A], Q_perp [1/A], Q_up [1/A], E [meV]):\n%s" % rot)
        print()
        print("Principal axes fwhms: %s" % fwhms)
        print("Coherent-elastic fwhms: %s" % fwhms_coh)
        print("Incoherent-elastic fwhms: %s" % fwhms_inc)


    # 2d sliced ellipses
    if verbose:
        print()
    Qres_QxE = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
    Qres_QxE = np.delete(np.delete(Qres_QxE, 1, axis=0), 1, axis=1)
    [fwhms_QxE, angles_QxE, rot_QxE, evals_QxE] = descr_ellipse(Qres_QxE)
    if verbose:
        print("2d Qx,E slice fwhms and slope angle: %s, %.4f" % (fwhms_QxE, angles_QxE[0]))

    Qres_QyE = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
    Qres_QyE = np.delete(np.delete(Qres_QyE, 0, axis=0), 0, axis=1)
    [fwhms_QyE, angles_QyE, rot_QyE, evals_QyE] = descr_ellipse(Qres_QyE)
    if verbose:
        print("2d Qy,E slice fwhms and slope angle: %s, %.4f" % (fwhms_QyE, angles_QyE[0]))

    Qres_QzE = np.delete(np.delete(Qres_Q, 1, axis=0), 1, axis=1)
    Qres_QzE = np.delete(np.delete(Qres_QzE, 0, axis=0), 0, axis=1)
    [fwhms_QzE, angles_QzE, rot_QzE, evals_QzE] = descr_ellipse(Qres_QzE)
    if verbose:
        print("2d Qz,E slice fwhms and slope angle: %s, %.4f" % (fwhms_QzE, angles_QzE[0]))

    Qres_QxQy = np.delete(np.delete(Qres_Q, 3, axis=0), 3, axis=1)
    Qres_QxQy = np.delete(np.delete(Qres_QxQy, 2, axis=0), 2, axis=1)
    [fwhms_QxQy, angles_QxQy, rot_QxQy, evals_QxQy] = descr_ellipse(Qres_QxQy)
    if verbose:
        print("2d Qx,Qy slice fwhms and slope angle: %s, %.4f" % (fwhms_QxQy, angles_QxQy[0]))


    # 2d projected ellipses
    if verbose:
        print()
    Qres_QxE_proj = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
    Qres_QxE_proj = quadric_proj(Qres_QxE_proj, 1)
    [fwhms_QxE_proj, angles_QxE_proj, rot_QxE_proj, evals_QxE_proj] = descr_ellipse(Qres_QxE_proj)
    if verbose:
        print("2d Qx,E projection fwhms and slope angle: %s, %.4f" % (fwhms_QxE_proj, angles_QxE_proj[0]))

    Qres_QyE_proj = np.delete(np.delete(Qres_Q, 2, axis=0), 2, axis=1)
    Qres_QyE_proj = quadric_proj(Qres_QyE_proj, 0)
    [fwhms_QyE_proj, angles_QyE_proj, rot_QyE_proj, evals_QyE_proj] = descr_ellipse(Qres_QyE_proj)
    if verbose:
        print("2d Qy,E projection fwhms and slope angle: %s, %.4f" % (fwhms_QyE_proj, angles_QyE_proj[0]))

    Qres_QzE_proj = np.delete(np.delete(Qres_Q, 1, axis=0), 1, axis=1)
    Qres_QzE_proj = quadric_proj(Qres_QzE_proj, 0)
    [fwhms_QzE_proj, angles_QzE_proj, rot_QzE_proj, evals_QzE_proj] = descr_ellipse(Qres_QzE_proj)
    if verbose:
        print("2d Qz,E projection fwhms and slope angle: %s, %.4f" % (fwhms_QzE_proj, angles_QzE_proj[0]))

    Qres_QxQy_proj = quadric_proj(Qres_Q, 3)
    Qres_QxQy_proj = np.delete(np.delete(Qres_QxQy_proj, 2, axis=0), 2, axis=1)
    [fwhms_QxQy_proj, angles_QxQy_proj, rot_QxQy_proj, evals_QxQy_proj] = descr_ellipse(Qres_QxQy_proj)
    if verbose:
        print("2d Qx,Qy projection fwhms and slope angle: %s, %.4f" % (fwhms_QxQy_proj, angles_QxQy_proj[0]))


    results = {
        # 4d ellipsoid
        "fwhms" : fwhms, "rot" : rot, "evals" : evals,
        "fwhms_coh" : fwhms_coh, "fwhms_inc" : fwhms_inc,

        # projected and sliced ellipses
        "fwhms_QxE" : fwhms_QxE, "rot_QxE" : rot_QxE,
        "fwhms_QyE" : fwhms_QyE, "rot_QyE" : rot_QyE,
        "fwhms_QzE" : fwhms_QzE, "rot_QzE" : rot_QzE,
        "fwhms_QxQy" : fwhms_QxQy,  "rot_QxQy" : rot_QxQy,
        "fwhms_QxE_proj" : fwhms_QxE_proj,  "rot_QxE_proj" : rot_QxE_proj,
        "fwhms_QyE_proj" : fwhms_QyE_proj, "rot_QyE_proj" : rot_QyE_proj,
        "fwhms_QzE_proj" : fwhms_QzE_proj, "rot_QzE_proj" : rot_QzE_proj,
        "fwhms_QxQy_proj" : fwhms_QxQy_proj, "rot_QxQy_proj" : rot_QxQy_proj,
    }

    return results



#
# shows the 2d ellipses
#
def plot_ellipses(ellis, verbose = True, plot_results = True, file = "", dpi = 600, ellipse_points = 128, use_tex = False):
    import mpl_toolkits.mplot3d as mplot3d
    import matplotlib
    import matplotlib.pyplot as plot

    matplotlib.rc("text", usetex = use_tex)


    ellfkt = lambda rad, rot, phi : \
        np.dot(rot, np.array([ rad[0]*np.cos(phi), rad[1]*np.sin(phi) ]))


    phi = np.linspace(0, 2.*np.pi, ellipse_points)

    ell_QxE = ellfkt(ellis["fwhms_QxE"]*0.5, ellis["rot_QxE"], phi)
    ell_QyE = ellfkt(ellis["fwhms_QyE"]*0.5, ellis["rot_QyE"], phi)
    ell_QzE = ellfkt(ellis["fwhms_QzE"]*0.5, ellis["rot_QzE"], phi)
    ell_QxQy = ellfkt(ellis["fwhms_QxQy"]*0.5, ellis["rot_QxQy"], phi)

    ell_QxE_proj = ellfkt(ellis["fwhms_QxE_proj"]*0.5, ellis["rot_QxE_proj"], phi)
    ell_QyE_proj = ellfkt(ellis["fwhms_QyE_proj"]*0.5, ellis["rot_QyE_proj"], phi)
    ell_QzE_proj = ellfkt(ellis["fwhms_QzE_proj"]*0.5, ellis["rot_QzE_proj"], phi)
    ell_QxQy_proj = ellfkt(ellis["fwhms_QxQy_proj"]*0.5, ellis["rot_QxQy_proj"], phi)

    labelQpara = "Qpara (1/A)"
    labelQperp = "Qperp (1/A)"
    labelQup = "Qup (1/A)"

    if use_tex:
        labelQpara = "$Q_{\\parallel}$ (\\AA$^{-1}$)"
        labelQperp = "$Q_{\\perp}$ (\\AA$^{-1}$)"
        labelQup = "$Q_{up}$ (\\AA$^{-1}$)"


    # Qpara, E axis
    fig = plot.figure()
    subplot_QxE = fig.add_subplot(221)
    subplot_QxE.set_xlabel(labelQpara)
    subplot_QxE.set_ylabel("E (meV)")
    subplot_QxE.plot(ell_QxE[0], ell_QxE[1], c="black", linestyle="dashed")
    subplot_QxE.plot(ell_QxE_proj[0], ell_QxE_proj[1], c="black", linestyle="solid")

    # Qperp, E axis
    subplot_QyE = fig.add_subplot(222)
    subplot_QyE.set_xlabel(labelQperp)
    subplot_QyE.set_ylabel("E (meV)")
    subplot_QyE.plot(ell_QyE[0], ell_QyE[1], c="black", linestyle="dashed")
    subplot_QyE.plot(ell_QyE_proj[0], ell_QyE_proj[1], c="black", linestyle="solid")

    # Qup, E axis
    subplot_QzE = fig.add_subplot(223)
    subplot_QzE.set_xlabel(labelQup)
    subplot_QzE.set_ylabel("E (meV)")
    subplot_QzE.plot(ell_QzE[0], ell_QzE[1], c="black", linestyle="dashed")
    subplot_QzE.plot(ell_QzE_proj[0], ell_QzE_proj[1], c="black", linestyle="solid")

    # Qpara, Qperp axis
    subplot_QxQy = fig.add_subplot(224)
    subplot_QxQy.set_xlabel(labelQpara)
    subplot_QxQy.set_ylabel(labelQperp)
    subplot_QxQy.plot(ell_QxQy[0], ell_QxQy[1], c="black", linestyle="dashed")
    subplot_QxQy.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], c="black", linestyle="solid")
    plot.tight_layout()


    # 3d plot
    fig3d = plot.figure()
    subplot3d = fig3d.add_subplot(111, projection="3d")

    subplot3d.set_xlabel(labelQpara)
    subplot3d.set_ylabel(labelQperp)
    #subplot3d.set_ylabel(labelQup)
    subplot3d.set_zlabel("E (meV)")

    # xE
    subplot3d.plot(ell_QxE[0], ell_QxE[1], zs=0., zdir="y", c="black", linestyle="dashed")
    subplot3d.plot(ell_QxE_proj[0], ell_QxE_proj[1], zs=0., zdir="y", c="black", linestyle="solid")
    # yE
    subplot3d.plot(ell_QyE[0], ell_QyE[1], zs=0., zdir="x", c="black", linestyle="dashed")
    subplot3d.plot(ell_QyE_proj[0], ell_QyE_proj[1], zs=0., zdir="x", c="black", linestyle="solid")
    # zE
    #subplot3d.plot(ell_QzE[0], ell_QzE[1], zs=0., zdir="x", c="black", linestyle="dashed")
    #subplot3d.plot(ell_QzE_proj[0], ell_QzE_proj[1], zs=0., zdir="x", c="black", linestyle="solid")
    # xy
    subplot3d.plot(ell_QxQy[0], ell_QxQy[1], zs=0., zdir="z", c="black", linestyle="dashed")
    subplot3d.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], zs=0., zdir="z", c="black", linestyle="solid")


    # save plots to files
    if file != "":
        import os
        splitext = os.path.splitext(file)
        file3d = splitext[0] + "_3d" + splitext[1]

        if verbose:
            print("Saving 2d plot to \"%s\"." % file)
            print("Saving 3d plot to \"%s\"." % file3d)
        fig.savefig(file, dpi=dpi)
        fig3d.savefig(file3d, dpi=dpi)


    # show plots
    if plot_results:
        plot.show()
