#!/usr/bin/env python
#
# resolution ellipsoid calculations
#
# @author Tobias Weber <tweber@ill.fr>
# @date mar-2019
# @license see 'LICENSE' file
#
# @desc for algorithm: [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
# @desc see covariance calculations: https://code.ill.fr/scientific-software/takin/mag-core/blob/master/tools/tascalc/cov.py
# @desc see also: https://github.com/McStasMcXtrace/McCode/blob/master/tools/Legacy-Perl/mcresplot.pl
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

import libs.helpers as helpers

import numpy as np
import numpy.linalg as la

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
    [fwhms_4d, angles_4d, rot_4d, evals_4d] = descr_ellipse(Qres_Q)

    axes_to_delete = [ [2, 1], [2, 0], [1, 0], [3, 2] ]
    slice_first = [ True, True, True, False ]
    results = []

    for ellidx in range(len(axes_to_delete)):
        # sliced 2d ellipse
        Qres = np.delete(np.delete(Qres_Q, axes_to_delete[ellidx][0], axis=0), axes_to_delete[ellidx][0], axis=1)
        Qres = np.delete(np.delete(Qres, axes_to_delete[ellidx][1], axis=0), axes_to_delete[ellidx][1], axis=1)
        [fwhms, angles, rot, evals] = descr_ellipse(Qres)

        # projected 2d ellipse
        if slice_first[ellidx]:
            Qres_proj = np.delete(np.delete(Qres_Q, axes_to_delete[ellidx][0], axis=0), axes_to_delete[ellidx][0], axis=1)
            Qres_proj = quadric_proj(Qres_proj, axes_to_delete[ellidx][1])
        else:
            Qres_proj = quadric_proj(Qres_Q, axes_to_delete[ellidx][0])
            Qres_proj = np.delete(np.delete(Qres_proj, axes_to_delete[ellidx][1], axis=0), axes_to_delete[ellidx][1], axis=1)
        [fwhms_proj, angles_proj, rot_proj, evals_proj] = descr_ellipse(Qres_proj)

        results.append({
            "fwhms" : fwhms, "angles" : angles, "rot" : rot, "evals" : evals,
            "fwhms_proj" : fwhms_proj, "angles_proj" : angles_proj, "rot_proj" : rot_proj, "evals_proj" : evals_proj,
            "fwhms_4d" : fwhms_4d, "angles_4d" : angles_4d, "rot_4d" : rot_4d, "evals_4d" : evals_4d, })

    if verbose:
        print()
        print("4d resolution ellipsoid principal axes fwhm lengths:\n%s" % fwhms_4d)
        print("4d resolution ellipsoid diagonal elements fwhm (coherent-elastic scattering) lengths:\n%s" \
            % (1./np.sqrt(np.abs(np.diag(Qres_Q))) * helpers.sig2fwhm))
        print()

        fwhms_coh = calc_coh_fwhms(Qres_Q)
        fwhms_inc = calc_incoh_fwhms(Qres_Q)

        print("Eigenvalues: %s" % evals_4d)
        print("Eigensystem (Q_para [1/A], Q_perp [1/A], Q_up [1/A], E [meV]):\n%s" % rot_4d)
        print()
        print("Principal axes fwhms: %s" % fwhms)
        print("Coherent-elastic fwhms: %s" % fwhms_coh)
        print("Incoherent-elastic fwhms: %s" % fwhms_inc)
        print()

        Qres_proj_Qpara = quadric_proj(quadric_proj(quadric_proj(Qres_Q, 3), 2), 1)
        Qres_proj_Qperp = quadric_proj(quadric_proj(quadric_proj(Qres_Q, 3), 2), 0)
        Qres_proj_Qup = quadric_proj(quadric_proj(quadric_proj(Qres_Q, 3), 1), 0)
        Qres_proj_E = quadric_proj(quadric_proj(quadric_proj(Qres_Q, 2), 1), 0)
        print("Projected fwhms (incoherent-elastic scattering) lengths:\n" \
            "%.4f / A, %.4f / A, %.4f / A, %.4f meV \n" % ( \
            1./np.sqrt(np.abs(Qres_proj_Qpara[0,0])) * helpers.sig2fwhm, \
            1./np.sqrt(np.abs(Qres_proj_Qperp[0,0])) * helpers.sig2fwhm, \
            1./np.sqrt(np.abs(Qres_proj_Qup[0,0])) * helpers.sig2fwhm, \
            1./np.sqrt(np.abs(Qres_proj_E[0,0])) * helpers.sig2fwhm ))

        print("Qx,E sliced ellipse fwhm and slope angle: %s, %.4f" % (results[0]["fwhms"], results[0]["angles"][0]))
        print("Qy,E sliced ellipse fwhm and slope angle: %s, %.4f" % (results[1]["fwhms"], results[1]["angles"][0]))
        print("Qz,E sliced ellipse fwhm and slope angle: %s, %.4f" % (results[2]["fwhms"], results[2]["angles"][0]))
        print("Qx,Qy sliced ellipse fwhm and slope angle: %s, %.4f" % (results[3]["fwhms"], results[3]["angles"][0]))
        print()
        print("Qx,E projected ellipse fwhm and slope angle: %s, %.4f" % (results[0]["fwhms_proj"], results[0]["angles_proj"][0]))
        print("Qy,E projected ellipse fwhm and slope angle: %s, %.4f" % (results[1]["fwhms_proj"], results[1]["angles_proj"][0]))
        print("Qz,E projected ellipse fwhm and slope angle: %s, %.4f" % (results[2]["fwhms_proj"], results[2]["angles_proj"][0]))
        print("Qx,Qy projected ellipse fwhm and slope angle: %s, %.4f" % (results[3]["fwhms_proj"], results[3]["angles_proj"][0]))

    return results


#
# shows the 2d ellipses
#
def plot_ellipses(ellis, Qs = np.array([]), Qmean = None, centre_on_Q = False,
    verbose = True, plot_results = True, plot_file = "",
    symsize = 1., dpi = 600, ellipse_points = 128, use_tex = False):
    try:
        import mpl_toolkits.mplot3d as mplot3d
        import matplotlib
        import matplotlib.pyplot as plot
    except ImportError:
        print("Matplotlib could not be imported!")
        exit(-1)

    matplotlib.rc("text", usetex = use_tex)

    themarker = "."


    ellfkt = lambda rad, rot, phi, Qmean2d : \
        np.dot(rot, np.array([ rad[0]*np.cos(phi), rad[1]*np.sin(phi) ])) + Qmean2d


    # 2d plots
    fig = plot.figure()

    num_ellis = len(ellis)
    coord_axes = [[0, 3], [1, 3], [2, 3], [0, 1]]
    coord_names = ["Qpara (1/A)", "Qperp (1/A)", "Qup (1/A)", "E (meV)"]

    if use_tex:
        coord_names[0] = "$Q_{\\parallel}$ (\\AA$^{-1}$)"
        coord_names[1] = "$Q_{\\perp}$ (\\AA$^{-1}$)"
        coord_names[2] = "$Q_{up}$ (\\AA$^{-1}$)"


    ellplots = []
    for ellidx in range(num_ellis):
        # centre plots on zero or mean Q vector ?
        QxE = np.array([[0], [0]])

        if centre_on_Q and Qmean != None:
            QxE = np.array([[Qmean[coord_axes[ellidx][0]]], [Qmean[coord_axes[ellidx][0]]]])


        phi = np.linspace(0, 2.*np.pi, ellipse_points)

        ell_QxE = ellfkt(ellis[ellidx]["fwhms"]*0.5, ellis[ellidx]["rot"], phi, QxE)
        ell_QxE_proj = ellfkt(ellis[ellidx]["fwhms_proj"]*0.5, ellis[ellidx]["rot_proj"], phi, QxE)
        ellplots.append({"sliced":ell_QxE, "proj":ell_QxE_proj})


        subplot_QxE = fig.add_subplot(221 + ellidx)
        subplot_QxE.set_xlabel(coord_names[coord_axes[ellidx][0]])
        subplot_QxE.set_ylabel(coord_names[coord_axes[ellidx][1]])

        if len(Qs.shape) == 2 and len(Qs) > 0 and len(Qs[0]) == 4:
            subplot_QxE.scatter(Qs[:, coord_axes[ellidx][0]], Qs[:, coord_axes[ellidx][1]],
                marker = themarker, s = symsize)

        subplot_QxE.plot(ell_QxE[0], ell_QxE[1], c="black", linestyle="dashed")
        subplot_QxE.plot(ell_QxE_proj[0], ell_QxE_proj[1], c="black", linestyle="solid")

    plot.tight_layout()


    # 3d plot
    fig3d = plot.figure()
    subplot3d = fig3d.add_subplot(111, projection="3d")

    subplot3d.set_xlabel(coord_names[0])
    subplot3d.set_ylabel(coord_names[1])
    subplot3d.set_zlabel(coord_names[3])

    if len(Qs.shape) == 2 and len(Qs) > 0 and len(Qs[0]) == 4:
        subplot3d.scatter(Qs[:,0], Qs[:,1], Qs[:,3], marker = themarker, s = symsize)

    # xE
    subplot3d.plot(ellplots[0]["sliced"][0], ellplots[0]["sliced"][1], zs=0., zdir="y", c="black", linestyle="dashed")
    subplot3d.plot(ellplots[0]["proj"][0], ellplots[0]["proj"][1], zs=0., zdir="y", c="black", linestyle="solid")
    # yE
    subplot3d.plot(ellplots[1]["sliced"][0], ellplots[1]["sliced"][1], zs=0., zdir="x", c="black", linestyle="dashed")
    subplot3d.plot(ellplots[1]["proj"][0], ellplots[1]["proj"][1], zs=0., zdir="x", c="black", linestyle="solid")
    # xy
    subplot3d.plot(ellplots[3]["sliced"][0], ellplots[3]["sliced"][1], zs=0., zdir="z", c="black", linestyle="dashed")
    subplot3d.plot(ellplots[3]["proj"][0], ellplots[3]["proj"][1], zs=0., zdir="z", c="black", linestyle="solid")


    if plot_file != "":
        splitext = os.path.splitext(plot_file)
        file3d = splitext[0] + "_3d" + splitext[1]

        if verbose:
            print("Saving 2d plot to \"%s\"." % plot_file)
            print("Saving 3d plot to \"%s\"." % file3d)
        fig.savefig(plot_file, dpi = dpi)
        fig3d.savefig(file3d, dpi = dpi)

    if plot_results:
        plot.show()
