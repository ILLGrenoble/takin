/**
 * Scan viewer
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2020
 * @license GPLv2
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

#ifndef __SCAN_EXPORTERS_H__
#define __SCAN_EXPORTERS_H__

#include "tlibs/string/string.h"
#include "tlibs/time/chrono.h"
#include "libs/version.h"


/**
 * convert scan file to gnuplot
 */
template<class t_vec>
std::string export_scan_to_gnuplot(
	const t_vec& vecX, const t_vec& vecY, const t_vec& vecYErr,
	const std::string& strLabelX = "",
	const std::string& strLabelY = "",
	const std::string& strTitle = "",
	const std::string& strFile = "")
{
	using t_real = typename t_vec::value_type;

	std::string strSrc =
R"RAWSTR(#!gnuplot --persist

#
# Created %%DATE%% with Takin version %%TAKIN_VER%%.
#

# --------------------------------------------------------------------------------
# choose an output terminal
#set term wxt
set term qt
#set term pdf color enhanced font "NimbusSans-Regular, 14" size 4,3.5
#set output "plot.pdf"
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# parameters
maxx = %%MAXX%%
minx = %%MINX%%
maxy = %%MAXY%%
miny = %%MINY%%
rangex = maxx - minx
rangey = maxy - miny

rangex_tics = rangex / 5.
rangey_tics = rangey / 5.
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# data points
$datapoints << ENDPTS
%%POINTS%%
ENDPTS
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# functions for fitting
gauss(x, a, s, x0, y0) = a*exp(-0.5 * ((x-x0)/s)**2.) + y0
gauss2(x, a1,s1,x01, a2,s2,x02, y0) = \
	gauss(x, a1,s1,x01, 0) + \
	gauss(x, a2,s2,x02, 0) + \
	y0
gauss3(x, a1,s1,x01, a2,s2,x02, a3,s3,x03, y0) = \
	gauss(x, a1,s1,x01, 0) + \
	gauss(x, a2,s2,x02, 0) + \
	gauss(x, a3,s3,x03, 0) + \
	y0

lorentz(x, a, h, x0, y0) = a*h**2. / ((x-x0)**2. + h**2.) + y0
lorentz2(x, a1,h1,x01, a2,h2,x02, y0) = \
	lorentz(x, a1,h1,x01, 0) + \
	lorentz(x, a2,h2,x02, 0) + \
	y0
lorentz3(x, a1,h1,x01, a2,h2,x02, a3,h3,x03, y0) = \
	lorentz(x, a1,h1,x01, 0) + \
	lorentz(x, a2,h2,x02, 0) + \
	lorentz(x, a3,h3,x03, 0) + \
	y0

parabola(x, a, x0, y0) = a * (x-x0)**2. + y0
line(x, m, y0) = m*x + y0
sine(x, a, f, p, y0) = a*sin(f*x + p) + y0

coth(x) = x==0 ? 0 : 1/tanh(x)
brillouin(x, J, x0, y0, scale) = scale*((1+1/(2*J)) * coth((1+1/(2*J))*(-x+x0)) - (1/(2*J)) * coth((1/(2*J))*(-x+x0))) + y0
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# fitting

# initial guesses for gauss and others ...
a1 = rangey
s1 = rangex*0.5		# sigma
x01 = minx + rangex*0.5
y0 = miny

# ... and for line
m1 = rangey / rangex

# ... and for sine
f1 = 1.		# frequency
p1 = 0.		# phase

# ... and for lorentz
h1 = rangex*0.5		# HWHM

# actual fitting
#fit line(x, m1,y0) "-" using ($1):($2):($3) yerrors via m1,y0
#fit parabola(x, a1,x01, y0) "-" using ($1):($2):($3) yerrors via a1,x01,y0
#fit sine(x, a1,f1,p1, y0) "-" using ($1):($2):($3) yerrors via a1,h1,x01,y0
#fit lorentz(x, a1,h1,x01, y0) "-" using ($1):($2):($3) yerrors via a1,h1,x01,y0
fit gauss(x, a1,s1,x01, y0) "$datapoints" using ($1):($2):($3) yerrors via a1,s1,x01,y0
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# plotting
line1_col = "#000000"
points1_col = "#000000"

set xlabel "%%LABELX%%"
set ylabel "%%LABELY%%"
set title "%%TITLE%%"
set grid
#set key top right Left width 0 samplen 3 spacing 1.2

set xrange [minx : maxx]
set yrange [miny : maxy]

#set xtics rangex_tics
#set ytics rangey_tics
#set mxtics 2
#set mytics 2

# show fit results
sig2fwhm = 2.*sqrt(2.*log(2.))
fitres = sprintf("x0 = %.4g\nFWHM = %.4g", x01, abs(sig2fwhm*s1))
set label at graph 0.05, graph 0.95 fitres

# use these functions with the plot command below
#	line(x,m1,y0) with lines linewidth 2 linecolor rgb line1_col notitle, \
#	parabola(x,a1,x01,y0) with lines linewidth 2 linecolor rgb line1_col notitle, \
#	sine(x,a1,f1,p1,y0) with lines linewidth 2 linecolor rgb line1_col notitle, \
#	lorentz(x,a1,h1,x01,y0) with lines linewidth 2 linecolor rgb line1_col notitle, \

plot \
	gauss(x,a1,s1,x01,y0) with lines linewidth 2 linecolor rgb line1_col notitle, \
	"$datapoints" using ($1):($2):($3) pointtype 7 pointsize 1 linecolor rgb points1_col with yerrorbars title "Data"
# --------------------------------------------------------------------------------
)RAWSTR";


	auto minmaxX = std::minmax_element(vecX.begin(), vecX.end());
	auto minmaxY = std::minmax_element(vecY.begin(), vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrPoints;
	ostrPoints.precision(g_iPrec);

	for(std::size_t i = 0; i < std::min(vecX.size(), vecY.size()); ++i)
	{
		ostrPoints
			<< std::left << std::setw(g_iPrec*2) << vecX[i] << " "
			<< std::left << std::setw(g_iPrec*2) << vecY[i] << " "
			<< std::left << std::setw(g_iPrec*2) << vecYErr[i] << "\n";
	}

	tl::find_all_and_replace<std::string>(strSrc, "%%MINX%%", tl::var_to_str(*minmaxX.first, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%TITLE%%", strTitle);
	tl::find_all_and_replace<std::string>(strSrc, "%%LABELX%%", strLabelX);
	tl::find_all_and_replace<std::string>(strSrc, "%%LABELY%%", strLabelY);
	tl::find_all_and_replace<std::string>(strSrc, "%%POINTS%%", ostrPoints.str());
	tl::find_all_and_replace<std::string>(strSrc, "%%TAKIN_VER%%", TAKIN_VER);
	tl::find_all_and_replace<std::string>(strSrc, "%%DATE%%", tl::epoch_to_str<t_real>(tl::epoch<t_real>(),
		"on %b %d, %Y at %H:%M:%S (%Z)"));

	return strSrc;
}



/**
 * convert scan file to python
 */
template<class t_vec>
std::string export_scan_to_python(
	const t_vec& vecX, const t_vec& vecY, const t_vec& vecYErr,
	const std::string& strLabelX = "",
	const std::string& strLabelY = "",
	const std::string& strTitle = "",
	const std::string& strFile = "")
{
	using t_real = typename t_vec::value_type;

	std::string strPySrc =
R"RAWSTR(#
# Created %%DATE%% with Takin version %%TAKIN_VER%%.
#

import numpy as np

x = np.array([ %%VECX%% ])
y = np.array([ %%VECY%% ])
yerr = np.array([ %%VECYERR%% ])

min = np.array([ %%MINX%%, %%MINY%% ])
max = np.array([ %%MAXX%%, %%MAXY%% ])
range = max-min
mid = min + range*0.5

yerr = [a if a!=0. else 0.001*range[1] for a in yerr]



import scipy.optimize as opt

def gauss_model(x, x0, sigma, amp, offs):
        return amp * np.exp(-0.5 * ((x-x0) / sigma)**2.) + offs

hints = [mid[0], range[0]*0.5, range[1]*0.5, min[1]]
popt, pcov = opt.curve_fit(gauss_model, x, y, sigma=yerr, absolute_sigma=True, p0=hints)

x_fine = np.linspace(min[0], max[0], 128)
y_fit = gauss_model(x_fine, *popt)



import matplotlib.pyplot as plt

plt.figure()

plt.xlim(min[0], max[0])
plt.ylim(min[1], max[1])

plt.title("%%TITLE%%")
plt.xlabel("%%LABELX%%")
plt.ylabel("%%LABELY%%")

plt.grid(True)
plt.errorbar(x,y,yerr, fmt="o")
plt.plot(x_fine, y_fit)
plt.show())RAWSTR";


	auto minmaxX = std::minmax_element(vecX.begin(), vecX.end());
	auto minmaxY = std::minmax_element(vecY.begin(), vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrX, ostrY, ostrYErr;
	ostrX.precision(g_iPrec);
	ostrY.precision(g_iPrec);
	ostrYErr.precision(g_iPrec);

	for(std::size_t i = 0; i < std::min(vecX.size(), vecY.size()); ++i)
	{
		ostrX << vecX[i] << ", ";
		ostrY << vecY[i] << ", ";
		ostrYErr << vecYErr[i] << ", ";
	}

	tl::find_all_and_replace<std::string>(strPySrc, "%%MINX%%", tl::var_to_str(*minmaxX.first, g_iPrec));
	tl::find_all_and_replace<std::string>(strPySrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second, g_iPrec));
	tl::find_all_and_replace<std::string>(strPySrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strPySrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strPySrc, "%%TITLE%%", strTitle);
	tl::find_all_and_replace<std::string>(strPySrc, "%%LABELX%%", strLabelX);
	tl::find_all_and_replace<std::string>(strPySrc, "%%LABELY%%", strLabelY);
	tl::find_all_and_replace<std::string>(strPySrc, "%%VECX%%", ostrX.str());
	tl::find_all_and_replace<std::string>(strPySrc, "%%VECY%%", ostrY.str());
	tl::find_all_and_replace<std::string>(strPySrc, "%%VECYERR%%", ostrYErr.str());
	tl::find_all_and_replace<std::string>(strPySrc, "%%TAKIN_VER%%", TAKIN_VER);
	tl::find_all_and_replace<std::string>(strPySrc, "%%DATE%%", tl::epoch_to_str<t_real>(tl::epoch<t_real>(),
		"on %b %d, %Y at %H:%M:%S (%Z)"));

	return strPySrc;
}



/**
 * convert scan file to root
 */
template<class t_vec>
std::string export_scan_to_root(
	const t_vec& vecX, const t_vec& vecY, const t_vec& vecYErr,
	const std::string& strLabelX = "",
	const std::string& strLabelY = "",
	const std::string& strTitle = "",
	const std::string& strFile = "")
{
	using t_real = typename t_vec::value_type;

	std::string strRootSrc =
R"RAWSTR(/*
 * Created %%DATE%% with Takin version %%TAKIN_VER%%.
 */

void scan_plot()
{
	const Double_t vecX[] = { %%VECX%% };
	const Double_t vecY[] = { %%VECY%% };
	const Double_t vecYErr[] = { %%VECYERR%% };

	const Double_t dMin[] = { %%MINX%%, %%MINY%% };
	const Double_t dMax[] = { %%MAXX%%, %%MAXY%% };
	const Int_t iSize = sizeof(vecX)/sizeof(*vecX);

	gStyle->SetOptFit(1);
	TCanvas *pCanvas = new TCanvas("canvas0", "Root Canvas", 800, 600);
	pCanvas->SetGrid(1,1);
	pCanvas->SetTicks(1,1);
	//pCanvas->SetLogy();

	TH1F *pFrame = pCanvas->DrawFrame(dMin[0], dMin[1], dMax[0], dMax[1], "");
	pFrame->SetTitle("%%TITLE%%");
	pFrame->SetXTitle("%%LABELX%%");
	pFrame->SetYTitle("%%LABELY%%");

	TGraphErrors *pGraph = new TGraphErrors(iSize, vecX, vecY, 0, vecYErr);
	pGraph->SetMarkerStyle(20);
	pGraph->Draw("P");
})RAWSTR";


	auto minmaxX = std::minmax_element(vecX.begin(), vecX.end());
	auto minmaxY = std::minmax_element(vecY.begin(), vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrX, ostrY, ostrYErr;
	ostrX.precision(g_iPrec);
	ostrY.precision(g_iPrec);
	ostrYErr.precision(g_iPrec);

	for(std::size_t i = 0; i < std::min(vecX.size(), vecY.size()); ++i)
	{
		ostrX << vecX[i] << ", ";
		ostrY << vecY[i] << ", ";
		ostrYErr << vecYErr[i] << ", ";
	}

	tl::find_all_and_replace<std::string>(strRootSrc, "%%MINX%%", tl::var_to_str(*minmaxX.first, g_iPrec));
	tl::find_all_and_replace<std::string>(strRootSrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second, g_iPrec));
	tl::find_all_and_replace<std::string>(strRootSrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strRootSrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strRootSrc, "%%TITLE%%", strTitle);
	tl::find_all_and_replace<std::string>(strRootSrc, "%%LABELX%%", strLabelX);
	tl::find_all_and_replace<std::string>(strRootSrc, "%%LABELY%%", strLabelY);
	tl::find_all_and_replace<std::string>(strRootSrc, "%%VECX%%", ostrX.str());
	tl::find_all_and_replace<std::string>(strRootSrc, "%%VECY%%", ostrY.str());
	tl::find_all_and_replace<std::string>(strRootSrc, "%%VECYERR%%", ostrYErr.str());
	tl::find_all_and_replace<std::string>(strRootSrc, "%%TAKIN_VER%%", TAKIN_VER);
	tl::find_all_and_replace<std::string>(strRootSrc, "%%DATE%%", tl::epoch_to_str<t_real>(tl::epoch<t_real>(),
		"on %b %d, %Y at %H:%M:%S (%Z)"));

	return strRootSrc;
}



/**
 * convert scan file to hermelin
 */
template<class t_vec>
std::string export_scan_to_hermelin(
	const t_vec& vecX, const t_vec& vecY, const t_vec& vecYErr,
	const std::string& strLabelX = "",
	const std::string& strLabelY = "",
	const std::string& strTitle = "",
	const std::string& strFile = "")
{
	using t_real = typename t_vec::value_type;

	std::string strStoatSrc =
R"RAWSTR(#!./hermelin -t

#
# Created %%DATE%% with Takin version %%TAKIN_VER%%.
#

module_init()
{
	import("apps/instr.scr");

	global fit_dbg = 1;
	global theterm = "wxt";
	global norm_to_mon = 1;
}

scan_plot()
{
	scanfile = "%%FILE%%";
	[instr, datx, daty, datyerr, xlab, ylab] = load_instr(scanfile, norm_to_mon);
	title = "\"" + scanfile + "\"";

	maxx = max(datx); minx = min(datx);
	maxy = max(daty); miny = min(daty);
	rangex = maxx-minx; rangey = maxy-miny;
	midx = minx + rangex*0.5;


	gauss_pos = [midx];
	gauss_amp = [rangey*0.5];
	gauss_sig = [rangex*0.5];
	bckgrd = miny;
	fitsteps = ["xxx x"];

	outfile = "";
	#outfile = "scan_plot.pdf";

	thefit = fit_gauss_manual_singlestep(datx, daty, datyerr, gauss_pos, gauss_amp, gauss_sig, bckgrd, fitsteps);


	plotmap = map();
	plotmap += ["xlimits" : minx + " " + maxx];
	#plotmap += ["ylimits" : "0 0.5"];

	plot_gausses(1, thefit, [datx, daty, datyerr], title, xlab, ylab, outfile, plotmap);
}

main(args)
{
	scan_plot();
})RAWSTR";

	tl::find_all_and_replace<std::string>(strStoatSrc, "%%FILE%%", strFile);
	tl::find_all_and_replace<std::string>(strStoatSrc, "%%TAKIN_VER%%", TAKIN_VER);
	tl::find_all_and_replace<std::string>(strStoatSrc, "%%DATE%%", tl::epoch_to_str<t_real>(tl::epoch<t_real>(),
		"on %b %d, %Y at %H:%M:%S (%Z)"));

	return strStoatSrc;
}



/**
 * TODO
 * convert scan file to julia
 */
template<class t_vec>
std::string export_scan_to_julia(
	const t_vec& vecX, const t_vec& vecY, const t_vec& vecYErr,
	const std::string& strLabelX = "",
	const std::string& strLabelY = "",
	const std::string& strTitle = "",
	const std::string& strFile = "")
{
	using t_real = typename t_vec::value_type;

	std::string strSrc =
R"RAWSTR(#
# Created %%DATE%% with Takin version %%TAKIN_VER%%.
#

x = (Array{Float64})([ %%VECX%% ])
y = (Array{Float64})([ %%VECY%% ])
yerr = (Array{Float64})([ %%VECYERR%% ])

min = (Array{Float64})([ %%MINX%%, %%MINY%% ])
max = (Array{Float64})([ %%MAXX%%, %%MAXY%% ])
range = max-min
mid = min + range*0.5

)RAWSTR";


	auto minmaxX = std::minmax_element(vecX.begin(), vecX.end());
	auto minmaxY = std::minmax_element(vecY.begin(), vecY.end());
	t_real dMaxErrY = *std::max_element(vecYErr.begin(), vecYErr.end());

	std::ostringstream ostrX, ostrY, ostrYErr;
	ostrX.precision(g_iPrec);
	ostrY.precision(g_iPrec);
	ostrYErr.precision(g_iPrec);

	for(std::size_t i = 0; i < std::min(vecX.size(), vecY.size()); ++i)
	{
		ostrX << vecX[i] << ", ";
		ostrY << vecY[i] << ", ";
		ostrYErr << vecYErr[i] << ", ";
	}

	tl::find_all_and_replace<std::string>(strSrc, "%%MINX%%", tl::var_to_str(*minmaxX.first, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%MAXX%%", tl::var_to_str(*minmaxX.second, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%MINY%%", tl::var_to_str(*minmaxY.first-dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%MAXY%%", tl::var_to_str(*minmaxY.second+dMaxErrY, g_iPrec));
	tl::find_all_and_replace<std::string>(strSrc, "%%TITLE%%", strTitle);
	tl::find_all_and_replace<std::string>(strSrc, "%%LABELX%%", strLabelX);
	tl::find_all_and_replace<std::string>(strSrc, "%%LABELY%%", strLabelY);
	tl::find_all_and_replace<std::string>(strSrc, "%%VECX%%", ostrX.str());
	tl::find_all_and_replace<std::string>(strSrc, "%%VECY%%", ostrY.str());
	tl::find_all_and_replace<std::string>(strSrc, "%%VECYERR%%", ostrYErr.str());
	tl::find_all_and_replace<std::string>(strSrc, "%%TAKIN_VER%%", TAKIN_VER);
	tl::find_all_and_replace<std::string>(strSrc, "%%DATE%%", tl::epoch_to_str<t_real>(tl::epoch<t_real>(),
		"on %b %d, %Y at %H:%M:%S (%Z)"));

	return strSrc;
}


#endif
