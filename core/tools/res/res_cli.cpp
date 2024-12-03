/**
 * minimalistic command line client
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date apr-2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "res_cli.h"
#include "tlibs/log/log.h"
#include "tlibs/string/string.h"
#include "libs/version.h"
#include "tools/monteconvo/TASReso.h"

#include <map>
#include <future>

namespace ublas = boost::numeric::ublas;

std::istream& istr = std::cin;
std::ostream& ostr = std::cout;

using t_real = t_real_reso;
template<class KEY, class VAL> using t_map = std::/*unordered_*/map<KEY, VAL>;
using t_mat = ublas::matrix<t_real>;
using t_vec = ublas::vector<t_real>;


// ----------------------------------------------------------------------------
// client function declarations
void show_help(const std::vector<std::string>& vecArgs);
void load_sample(const std::vector<std::string>& vecArgs);
void load_instr(const std::vector<std::string>& vecArgs);
void fix(const std::vector<std::string>& vecArgs);
void calc(const std::vector<std::string>& vecArgs);
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// globals
TASReso g_tas;

using t_func = void(*)(const std::vector<std::string>&);
using t_funcmap = t_map<std::string, t_func>;

t_funcmap g_funcmap =
{
	{"help", &show_help},
	{"load_sample", &load_sample},
	{"load_instr", &load_instr},
	{"fix", &fix},
	{"calc", &calc},
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// client functions
void show_help(const std::vector<std::string>& vecArgs)
{
	std::string strHelp = "Available client functions: ";

	for(t_funcmap::const_iterator iter=g_funcmap.begin(); iter!=g_funcmap.end(); ++iter)
	{
		const t_funcmap::value_type& pair = *iter;
		strHelp += pair.first;
		if(std::next(iter) != g_funcmap.end())
			strHelp += ", ";
	}

	ostr << strHelp << ".\n";
}


void load_sample(const std::vector<std::string>& vecArgs)
{
	if(vecArgs.size() < 2)
	{
		ostr << "Error: No filename given.\n";
		return;
	}

	if(g_tas.LoadLattice(vecArgs[1].c_str()))
		ostr << "OK.\n";
	else
		ostr << "Error: Unable to load " << vecArgs[1] << ".\n";
}


void load_instr(const std::vector<std::string>& vecArgs)
{
	if(vecArgs.size() < 2)
	{
		ostr << "Error: No filename given.\n";
		return;
	}

	if(g_tas.LoadRes(vecArgs[1].c_str()))
		ostr << "OK.\n";
	else
		ostr << "Error: Unable to load " << vecArgs[1] << ".\n";
}


void fix(const std::vector<std::string>& vecArgs)
{
	if(vecArgs.size() < 3)
	{
		ostr << "Error: No variable or value given.\n";
		return;
	}

	if(vecArgs[1] == "ki")
		g_tas.SetKiFix(true);
	else if(vecArgs[1] == "kf")
		g_tas.SetKiFix(false);
	else
	{
		ostr << "Error: Unknown variable " << vecArgs[1] << ".\n";
		return;
	}

	const t_real dVal = tl::str_to_var<t_real>(vecArgs[2]);
	g_tas.SetKFix(dVal);
	ostr << "OK.\n";
}


std::ostream& operator<<(std::ostream& ostr, const t_mat& m)
{
	for(std::size_t i=0; i<m.size1(); ++i)
	{
		for(std::size_t j=0; j<m.size2(); ++j)
			ostr << m(i,j) << " ";

		ostr << " ";
	}

	return ostr;
}


std::ostream& operator<<(std::ostream& ostr, const t_vec& v)
{
	for(std::size_t i=0; i<v.size(); ++i)
		ostr << v(i) << " ";

	return ostr;
}


void calc(const std::vector<std::string>& vecArgs)
{
	if(vecArgs.size() < 5)
	{
		ostr << "Error: No hkl and E position given.\n";
		return;
	}

	const t_real dH = tl::str_to_var<t_real>(vecArgs[1]);
	const t_real dK = tl::str_to_var<t_real>(vecArgs[2]);
	const t_real dL = tl::str_to_var<t_real>(vecArgs[3]);
	const t_real dE = tl::str_to_var<t_real>(vecArgs[4]);

	const ResoResults& res = g_tas.GetResoResults();

	//g_tas.GetResoParams().flags |= CALC_R0;

	if(!g_tas.SetHKLE(dH, dK, dL, dE))
	{
		ostr << "Error: At postion Q=("
			<< dH << "," << dK << "," << dL
			<< "), E=" << dE << ": " << res.strErr
			<< ".\n";
		return;
	}

	//Ellipsoid4d<t_real> ell4d = calc_res_ellipsoid4d(res.reso, res.Q_avg);

	// parameters are: x, y, project 1, project 2, remove 1, remove 2
	int iParams[2][4][6] =
	{
		{	// projected
			{ 0, 3, 1, -1, 2, -1 },
			{ 1, 3, 0, -1, 2, -1 },
			{ 2, 3, 0, -1, 1, -1 },
			{ 0, 1, 3, -1, 2, -1 }
		},
		{	// sliced
			{ 0, 3, -1, -1, 2, 1 },
			{ 1, 3, -1, -1, 2, 0 },
			{ 2, 3, -1, -1, 1, 0 },
			{ 0, 1, -1, -1, 2, 3 }
		}
	};


	ostr << "OK.\n";

	ostr << "Reso: " << res.reso << "\n";
	ostr << "R0: " << res.dR0 << "\n";
	ostr << "Vol: " << res.dResVol << "\n";
	ostr << "Q_avg: " << res.Q_avg << "\n";
	ostr << "Bragg_FWHMs: "
		<< res.dBraggFWHMs[0] << " "
		<< res.dBraggFWHMs[1] << " "
		<< res.dBraggFWHMs[2] << " "
		<< res.dBraggFWHMs[3] << "\n";


	std::vector<std::future<Ellipse2d<t_real>>> tasks_ell_proj, tasks_ell_slice;

	for(unsigned int iEll=0; iEll<4; ++iEll)
	{
		const int *iP = iParams[0][iEll];
		const int *iS = iParams[1][iEll];

		const t_vec& Q_avg = res.Q_avg;
		const t_mat& reso = res.reso;
		const t_vec& reso_v = res.reso_v;
		const t_real& reso_s = res.reso_s;

		std::future<Ellipse2d<t_real>> ell_proj =
			std::async(std::launch::deferred|std::launch::async,
			[=, &reso, &Q_avg]()
			{
				return ::calc_res_ellipse<t_real>(
					reso, reso_v, reso_s,
					Q_avg, iP[0], iP[1], iP[2], iP[3], iP[4], iP[5]);
			});
		std::future<Ellipse2d<t_real>> ell_slice =
			std::async(std::launch::deferred|std::launch::async,
			[=, &reso, &Q_avg]()
			{
				return ::calc_res_ellipse<t_real>(
					reso, reso_v, reso_s,
					Q_avg, iS[0], iS[1], iS[2], iS[3], iS[4], iS[5]);
			});

		tasks_ell_proj.push_back(std::move(ell_proj));
		tasks_ell_slice.push_back(std::move(ell_slice));
	}
	for(unsigned int iEll=0; iEll<4; ++iEll)
	{
		Ellipse2d<t_real> elliProj = tasks_ell_proj[iEll].get();
		Ellipse2d<t_real> elliSlice = tasks_ell_slice[iEll].get();
		const std::string& strLabX = ::ellipse_labels(iParams[0][iEll][0], EllipseCoordSys::Q_AVG);
		const std::string& strLabY = ::ellipse_labels(iParams[0][iEll][1], EllipseCoordSys::Q_AVG);

		ostr << "Ellipse_" << iEll << "_labels: " << strLabX << ", " << strLabY << "\n";

		ostr << "Ellipse_" << iEll << "_proj_angle: " << elliProj.phi << "\n";
		ostr << "Ellipse_" << iEll << "_proj_HWHMs: " << elliProj.x_hwhm << " " << elliProj.y_hwhm << "\n";
		ostr << "Ellipse_" << iEll << "_proj_offs: " << elliProj.x_offs << " " << elliProj.y_offs << "\n";
		//ostr << "Ellipse_" << iEll << "_proj_area: " << elliProj.area << "\n";

		ostr << "Ellipse_" << iEll << "_slice_angle: " << elliSlice.phi << "\n";
		ostr << "Ellipse_" << iEll << "_slice_HWHMs: " << elliSlice.x_hwhm << " " << elliSlice.y_hwhm << "\n";
		ostr << "Ellipse_" << iEll << "_slice_offs: " << elliSlice.x_offs << " " << elliSlice.y_offs << "\n";
		//ostr << "Ellipse_" << iEll << "_slice_area: " << elliSlice.area << "\n";
	}

	ostr.flush();
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
int res_main(int argc, char** argv)
{
	try
	{
		std::size_t cmd_idx = 0;
		auto write_prompt = [&cmd_idx]()
		{
			++cmd_idx;
			ostr << "\n" << cmd_idx << "> ";
		};

		tl::log_info("This is Takin/Reso, version " TAKIN_VER " (built on " __DATE__ ").");
		ostr << TAKIN_LICENSE("Takin/Reso") << std::endl;

		ostr << "\n";
		show_help({});
		write_prompt();

		std::string strLine;
		while(std::getline(istr, strLine))
		{
			std::vector<std::string> vecToks;
			tl::get_tokens<std::string, std::string, decltype(vecToks)>
				(strLine, " \t", vecToks);

			for(std::string& strTok : vecToks)
				tl::trim(strTok);

			if(!vecToks.size()) continue;

			if(vecToks[0] == "exit")
				break;

			t_funcmap::const_iterator iter = g_funcmap.find(vecToks[0]);
			if(iter == g_funcmap.end())
			{
				ostr << "Error: No such function: " << vecToks[0] << ".\n" << std::endl;
				write_prompt();
				continue;
			}

			(*iter->second)(vecToks);
			write_prompt();
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

	return 0;
}
// ----------------------------------------------------------------------------
