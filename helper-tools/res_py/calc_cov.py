#!/usr/bin/env python3
#
# calculate covariance from neutron events
# @author Tobias Weber <tweber@ill.fr>
# @date 30-mar-2019
# @license see 'LICENSE' file
#
# @desc For a good explanation of the covariance matrix method, see (Arens 2015), pp. 795 and 1372.
# @desc reimplements the functionality of https://github.com/McStasMcXtrace/McCode/blob/master/tools/Legacy-Perl/mcresplot.pl
# @desc see also [eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014), doi: 10.1016/j.nima.2014.03.019
#

import libs.tas as tas
import libs.reso as reso
import libs.helpers as helpers

import os

try:
    import numpy as np
    import numpy.linalg as la
except ImportError:
    print("Numpy could not be imported!")
    exit(-1)


try:
    np.set_printoptions(
        precision = 4,
        floatmode = "fixed")
except TypeError:
    print("Warning: Numpy print options could not be set.")


options = {
    "verbose" : True,        # console outputs
    "plot_results" : True,   # show plot window
    "plot_neutrons" : True,  # also plot neutron events
    "centre_on_Q" : False,   # centre plots on Q or zero?
    "ellipse_points" : 128,  # number of points to draw ellipses
    "symsize" : 32.,
    "dpi" : 600,
    "use_tex" : False,

    # column indices in ki,kf files
    "ki_start_idx" : 0,      # start index of ki 3-vector
    "kf_start_idx" : 3,      # start index of kf 3-vector
    "wi_idx" : 9,            # start index of ki weight factor
    "wf_idx" : 10,           # start index of kf weight factor

    # column indices in Q,E files
    "Q_start_idx" : 0,
    "E_idx" : 3,
    "w_idx" : 4,

    "filter_eps" : 1e-4,
}



#
# normalises events and filters out too low probabilities
#
def filter_events(Q, E, w):
    if w.size == 0:
        raise ValueError("No neutron events available.")

    # normalise intensity/probability
    maxw = np.max(w)
    if np.abs(maxw) < np.finfo(w[0].__class__).eps:
        raise ValueError("Neutron probability factors are zero.")
    w /= maxw

    # filter out too low probabilities
    beloweps = lambda d: np.abs(d) <= options["filter_eps"]
    nonzero_idx = [i for i in range(len(w)) if not beloweps(w[i])]

    Q = Q[nonzero_idx]
    E = E[nonzero_idx]
    w = w[nonzero_idx]

    if w.size == 0:
        raise ValueError("No neutron events left after filtering.")

    return [Q, E, w]



#
# loads a list of neutron events in the [ ki_vec, kf_vec, pos_vec, wi, wf ] format
#
def load_events_kikf(filename):
    dat = np.loadtxt(filename)
    ki = dat[:, options["ki_start_idx"]:options["ki_start_idx"] + 3]
    kf = dat[:, options["kf_start_idx"]:options["kf_start_idx"] + 3]
    wi = dat[:, options["wi_idx"]]
    wf = dat[:, options["wf_idx"]]

    w = wi * wf
    Q = ki - kf
    E = tas.k2_to_E * (np.multiply(ki, ki).sum(1) - np.multiply(kf, kf).sum(1))

    return filter_events(Q, E, w)



#
# loads a list of neutron events in the [ h, k, l, E, w ] format
#
def load_events_QE(filename):
    dat = np.loadtxt(filename)

    Q = dat[:, options["Q_start_idx"]:options["Q_start_idx"] + 3]
    E = dat[:, options["E_idx"]]
    w = dat[:, options["w_idx"]]

    return filter_events(Q, E, w)



#
# calculates the covariance matrix of the (Q, E) 4-vectors
#
def calc_covar(Q, E, w, Qpara, Qperp):
    # make a [Q, E] 4-vector
    Q4 = np.insert(Q, 3, E, axis = 1)

    # calculate the mean Q 4-vector
    Qmean = np.array([ np.average(Q4[:, i], weights = w) for i in range(4) ])
    if options["verbose"]:
        print("Mean (Q, E) vector in lab system:\n%s\n" % Qmean)

    # get the weighted covariance matrix
    Qcov = np.cov(Q4, rowvar = False, aweights = np.abs(w), ddof = 0)
    if options["verbose"]:
        print("Covariance matrix in lab system:\n%s\n" % Qcov)

    # the resolution is the inverse of the covariance
    Qres = la.inv(Qcov)
    if options["verbose"]:
        print("Resolution matrix in lab system:\n%s\n" % Qres)


    # create a matrix to transform into the coordinate system with Q along x
    # choose given coordinate system
    if len(Qpara) == 3 and len(Qperp) == 3:
        Qnorm = Qpara / la.norm(Qpara)
        Qside = Qperp / la.norm(Qperp)
        Qup = np.cross(Qnorm, Qside)
    else:
        Qnorm = Qmean[0:3] / la.norm(Qmean[0:3])
        Qup = np.array([ 0., 1., 0. ])
        Qside = np.cross(Qup, Qnorm)

    if options["verbose"]:
        print("Qpara = %s\nQperp = %s\nQup = %s\n" % (Qnorm, Qside, Qup))

    # trafo matrix
    T = np.transpose(np.array([
        np.insert(Qnorm, 3, 0),
        np.insert(Qside, 3, 0),
        np.insert(Qup, 3, 0),
        [0., 0., 0., 1.] ]))

    if options["verbose"]:
        print("Transformation into (Qpara, Qperp, Qup, E) system:\n%s\n" % T)

    # transform mean Q vector
    Qmean_Q = np.dot(np.transpose(T), Qmean)
    if options["verbose"]:
        print("Mean (Q, E) vector in (Qpara, Qperp, Qup, E) system:\n%s\n" % Qmean_Q)

    # transform the covariance matrix
    Qcov_Q = np.dot(np.transpose(T), np.dot(Qcov, T))
    if options["verbose"]:
        print("Covariance matrix in (Qpara, Qperp, Qup, E) system:\n%s\n" % Qcov_Q)

    # the resolution is the inverse of the covariance
    Qres_Q = la.inv(Qcov_Q)
    if options["verbose"]:
        print("Resolution matrix in (Qpara, Qperp, Qup, E) system:\n%s\n" % Qres_Q)

    #[ evals, evecs ] = la.eig(Qcov_Q)
    #print("Ellipsoid fwhm radii:\n%s\n" % (np.sqrt(np.abs(evals)) * helpers.sig2fwhm))

    # transform all neutron events
    Q4_Q = np.array([ ])
    if options["plot_neutrons"]:
        Q4_Q = np.dot(Q4, T)
        if not options["centre_on_Q"]:
            Q4_Q -= Qmean_Q


    return [Qres_Q, Q4_Q, Qmean_Q]



#
# checks versions of needed packages
#
def check_versions():
    npver = np.version.version.split(".")
    if int(npver[0]) >= 2:
        return
    if int(npver[0]) < 1 or int(npver[1]) < 10:
        print("Numpy version >= 1.10 is required, but installed version is %s." % np.version.version)
        exit(-1)



#
# entry point
#
def run_cov():
    print("This is a covariance matrix calculator using neutron events,\n\twritten by T. Weber <tweber@ill.fr>, 30 March 2019.\n")
    check_versions()

    try:
        import argparse as arg
    except ImportError:
        print("Argparse could not be imported!")
        exit(-1)

    args = arg.ArgumentParser(description="Calculates the covariance matrix of neutron scattering events.")
    args.add_argument("file", type=str, help="input file")
    args.add_argument("-s", "--save", default="", type=str, nargs="?", help="save plot to file")
    args.add_argument("--ellipse", default=options["ellipse_points"], type=int, nargs="?", help="number of points to draw ellipses")
    args.add_argument("--ki", default=options["ki_start_idx"], type=int, nargs="?", help="index of ki vector's first column in kikf file")
    args.add_argument("--kf", default=options["kf_start_idx"], type=int, nargs="?", help="index of kf vector's first column in kikf file")
    args.add_argument("--wi", default=options["wi_idx"], type=int, nargs="?", help="index of ki weight factor column in kikf file")
    args.add_argument("--wf", default=options["wf_idx"], type=int, nargs="?", help="index of kf weight factor column in kikf file")
    args.add_argument("--w", default=options["w_idx"], type=int, nargs="?", help="index of neutron weight factor column in QE file")
    args.add_argument("--Q", default=options["Q_start_idx"], type=int, nargs="?", help="index of Q vector's first column in QE file")
    args.add_argument("--E", default=options["E_idx"], type=int, nargs="?", help="index of E column in QE file")
    args.add_argument("--QEfile", action="store_true", help="use the QE file type")
    args.add_argument("--centreonQ", action="store_true", help="centre plots on mean Q")
    args.add_argument("--noverbose", action="store_true", help="don't show console logs")
    args.add_argument("--noplot", action="store_true", help="don't show any plot windows")
    args.add_argument("--noneutrons", action="store_true", help="don't show neutron events in plots")
    args.add_argument("--symsize", default=options["symsize"], type=float, nargs="?", help="size of the symbols in plots")
    args.add_argument("--ax", default=None, type=float, nargs="?", help="x component of first orientation vector")
    args.add_argument("--ay", default=None, type=float, nargs="?", help="y component of first orientation vector")
    args.add_argument("--az", default=None, type=float, nargs="?", help="z component of first orientation vector")
    args.add_argument("--bx", default=None, type=float, nargs="?", help="x component of second orientation vector")
    args.add_argument("--by", default=None, type=float, nargs="?", help="y component of second orientation vector")
    args.add_argument("--bz", default=None, type=float, nargs="?", help="z component of second orientation vector")
    args.add_argument("--a", "--as", "-a", default=None, type=float, nargs="?", help="lattice constant a (only needed in case data is in rlu)")
    args.add_argument("--b", "--bs", "-b", default=None, type=float, nargs="?", help="lattice constant b (only needed in case data is in rlu)")
    args.add_argument("--c", "--cs", "-c", default=None, type=float, nargs="?", help="lattice constant c (only needed in case data is in rlu)")
    args.add_argument("--aa", "--alpha", default=90., type=float, nargs="?", help="lattice angle alpha (only needed in case data is in rlu)")
    args.add_argument("--bb", "--beta", default=90., type=float, nargs="?", help="lattice angle beta (only needed in case data is in rlu)")
    args.add_argument("--cc", "--gamma", default=90., type=float, nargs="?", help="lattice angle gamma (only needed in case data is in rlu)")
    args.add_argument("--filtereps", default=options["filter_eps"], type=float, nargs="?", help="epsilon probability below which neutron events are filtered out")
    args.add_argument("--dpi", default=options["dpi"], type=int, nargs="?", help="DPI of output plot file")
    args.add_argument("--tex", action="store_true", help="use tex in plots")
    argv = args.parse_args()

    options["verbose"] = (argv.noverbose == False)
    options["plot_results"] = (argv.noplot == False)
    options["plot_neutrons"] = (argv.noneutrons == False)
    options["centre_on_Q"] = argv.centreonQ
    options["dpi"] = argv.dpi
    options["use_tex"] = argv.tex

    B = []
    if argv.a != None and argv.b != None and argv.c != None and argv.aa != None and argv.bb != None and argv.cc != None:
        lattice = np.array([ argv.a, argv.b, argv.c,  ])
        angles = np.array([ argv.aa, argv.bb, argv.cc ]) / 180.*np.pi

        B = tas.get_B(lattice, angles)

        if options["verbose"]:
            print("Crystal B matrix:\n%s\n" % B)


    infile = argv.file
    outfile = argv.save
    options["ellipse_points"] = argv.ellipse
    options["ki_start_idx"] = argv.ki
    options["kf_start_idx"] = argv.kf
    options["wi_idx"] = argv.wi
    options["wf_idx"] = argv.wf
    options["filter_eps"] = argv.filtereps
    options["symsize"] = argv.symsize
    avec = [ argv.az, argv.ay, argv.az ]
    bvec = [ argv.bx, argv.by, argv.bz ]


    try:
        # input file is in h k l E w format?
        if argv.QEfile:
            [Q, E, w] = load_events_QE(infile)
            # convert rlu to 1/A
            if len(B) != 0:
                Q = np.dot(Q, np.transpose(B))
        # input file is in the kix kiy kiz kfx kfy kfz wi wf format?
        else:
            [Q, E, w] = load_events_kikf(infile)
    except OSError:
        print("Could not load input file %s." % infile)
        exit(-1)
    except NameError:
        print("Error processing input file %s." % infile)
        exit(-1)

    print("Loaded %d neutron events." % len(Q))


    if avec[0] != None and avec[1] != None and avec[2] != None and bvec[0] != None and bvec[1] != None and bvec[2] != None:
        Qpara = np.array(avec)
        Qperp = np.array(bvec)

        # convert rlu to 1/A
        if len(B) != 0:
            Qpara = np.dot(B, Qpara)
            Qperp = np.dot(B, Qperp)
    else:
        Qpara = np.array([])
        Qperp = np.array([])

    [Qres, Q4, Qmean] = calc_covar(Q, E, w, Qpara, Qperp)
    calcedellis = reso.calc_ellipses(Qres, verbose = options["verbose"])

    if options["plot_results"] or outfile!="":
        reso.plot_ellipses(calcedellis, plot_file = outfile, Qs = Q4, Qmean = Qmean,
            centre_on_Q = options["centre_on_Q"], symsize = options["symsize"] * w,
            ellipse_points = options["ellipse_points"], dpi = options["dpi"],
            use_tex = options["use_tex"], verbose = options["verbose"])


#
# main
#
if __name__ == "__main__":
    run_cov()
