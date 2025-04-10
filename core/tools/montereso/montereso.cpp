/**
 * Montereso
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2012, 22-sep-2014
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

#include "montereso.h"
#include "res.h"
#include "tlibs/log/log.h"
#include "tlibs/string/string.h"
#include "tlibs/math/rand.h"
#include "libs/qt/qthelper.h"
#include "libs/version.h"
#include "dialogs/EllipseDlg.h"
#include "../res/ellipse.h"

#include <clocale>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>

#include <boost/program_options.hpp>
namespace opts = boost::program_options;

using namespace ublas;
using t_real = t_real_reso;


static void add_param(std::unordered_map<std::string, std::string>& map, const std::string& strLine)
{
	std::string str = strLine.substr(1);

	std::pair<std::string, std::string> pair = tl::split_first<std::string>(str, ":", true);
	if(pair.first == "Param")
	{
		std::pair<std::string, std::string> pairParam = tl::split_first<std::string>(pair.second, "=", true);
		pair.first = pair.first + "_" + pairParam.first;
		pair.second = pairParam.second;
	}

	map.insert(pair);
}


template<class t_map>
static void print_map(std::ostream& ostr, const t_map& map)
{
	for(typename t_map::const_iterator iter=map.begin(); iter!=map.end(); ++iter)
		ostr << (*iter).first << ": " << (*iter).second << "\n";
}

enum class FileType
{
	NEUTRON_Q_LIST,
	NEUTRON_KIKF_LIST,

	RESOLUTION_MATRIX,
	COVARIANCE_MATRIX,

	UNKNOWN
};


static bool load_mat(const char* pcFile, Resolution& reso, FileType ft)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
	{
		tl::log_err("Cannot open \"", pcFile, "\".");
		return false;
	}

	matrix<t_real>& res = reso.res;
	matrix<t_real>& cov = reso.cov;
	res.resize(4,4,0);
	cov.resize(4,4,0);

	if(ft == FileType::COVARIANCE_MATRIX)
	{
		for(unsigned int i=0; i<4; ++i)
			for(unsigned int j=0; j<4; ++j)
				ifstr >> cov(i,j);

		reso.bHasRes = tl::inverse(cov, res);
	}
	else if(ft == FileType::RESOLUTION_MATRIX)
	{
		for(unsigned int i=0; i<4; ++i)
			for(unsigned int j=0; j<4; ++j)
				ifstr >> res(i,j);

		reso.bHasRes = tl::inverse(res, cov);
	}

	tl::log_info("Covariance matrix: ", cov);
	tl::log_info("Resolution matrix: ", res);
	tl::log_info("Matrix valid: ", reso.bHasRes);

	if(reso.bHasRes)
	{
		reso.dQ = calc_bragg_fwhms(reso.res);
		reso.dEinc = calc_vanadium_fwhms(reso.res)[3];

		std::ostringstream ostrVals;
		ostrVals << "Coherent / Bragg FWHM values (Qx, Qy, Qz, E): ";
		std::copy(reso.dQ.begin(), reso.dQ.end(), std::ostream_iterator<t_real>(ostrVals, ", "));

		std::ostringstream ostrIncVals;
		ostrIncVals << "Incoherent / Vanadium FWHM value (E in meV): " << reso.dEinc;

		tl::log_info(ostrVals.str());
		tl::log_info(ostrIncVals.str());
	}

	return reso.bHasRes;
}


static bool load_mc_list(const char* pcFile, Resolution& res,
	const ublas::vector<t_real> *qPara = nullptr,
	const ublas::vector<t_real> *qPerp = nullptr,
	bool bSwapYZ = false)
{
	FileType ft = FileType::NEUTRON_Q_LIST;

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
	{
		tl::log_err("Cannot open \"", pcFile, "\".");
		return false;
	}

	// neutron Q,E list
	std::vector<vector<t_real>> vecQ;
	std::vector<t_real> vecP;

	// neutron ki, kf list
	std::vector<vector<t_real>> vecKi, vecKf;
	std::vector<t_real> vecPi, vecPf;
	std::string strLine;

	std::unordered_map<std::string, std::string> mapParams;
	bool bEndOfHeader = false;

	std::size_t uiNumNeutr = 0;
	while(std::getline(ifstr, strLine))
	{
		tl::trim(strLine);
		if(strLine.length() == 0)
			continue;
		else if(strLine[0] == '#')
		{
			add_param(mapParams, strLine);
			continue;
		}

		if(!bEndOfHeader)
		{
			bEndOfHeader = true;

			try
			{
				if(mapParams.at("variables") == "ki_x ki_y ki_z kf_x kf_y kf_z x y z p_i p_f" ||
					mapParams.at("variables") == "ki_x ki_y ki_z kf_x kf_y kf_z x y z p_i p_f n")
				{
					tl::log_info("File is a ki, kf list.");
					ft = FileType::NEUTRON_KIKF_LIST;
				}
			}
			catch(const std::out_of_range& ex)
			{
				tl::log_info("File is a Q list.");
				ft = FileType::NEUTRON_Q_LIST;
			}
		}

		std::istringstream istr(strLine);

		if(ft == FileType::NEUTRON_Q_LIST)
		{
			t_real dQh=0., dQk=0., dQl=0., dE=0., dP=1.;
			istr >> dQh >> dQk >> dQl >> dE >> dP;

			if(bSwapYZ)
			{
				std::swap(dQk, dQl);
				dQh = -dQh;
			}

			if(!tl::float_equal(dP, t_real{0.}))
			{
				vector<t_real> _vec = tl::make_vec<vector<t_real>>({ dQh, dQk, dQl, dE });

				vecQ.push_back(std::move(_vec));
				vecP.push_back(dP);
			}
		}
		else if(ft == FileType::NEUTRON_KIKF_LIST)
		{
			t_real dKi[3], dKf[3], dPos[3], dPi=0., dPf=0.;

			istr >> dKi[0] >> dKi[1] >> dKi[2];
			istr >> dKf[0] >> dKf[1] >> dKf[2];
			istr >> dPos[0] >> dPos[1] >> dPos[2];
			istr >> dPi >> dPf;

			if(bSwapYZ)
			{
				std::swap(dKi[1], dKi[2]);
				std::swap(dKf[1], dKf[2]);
				dKi[0] = -dKi[0];
				dKf[0] = -dKf[0];
			}

			vecKi.emplace_back(tl::make_vec<vector<t_real>>({ dKi[0], dKi[1], dKi[2] }));
			vecKf.emplace_back(tl::make_vec<vector<t_real>>({ dKf[0], dKf[1], dKf[2] }));
			vecPi.push_back(dPi);
			vecPf.push_back(dPf);
		}

		++uiNumNeutr;
	}

	tl::log_info("Number of neutrons in file: ", uiNumNeutr);
	//print_map(std::cout, mapParams);


	if(ft == FileType::NEUTRON_Q_LIST)
	{
		normalise_P(&vecP);
		res = calc_res(std::forward<decltype(vecQ)&&>(vecQ), &vecP, qPara, qPerp);
	}
	else if(ft == FileType::NEUTRON_KIKF_LIST)
	{
		res = calc_res(vecKi, vecKf, &vecPi, &vecPf, qPara, qPerp);
	}


	if(!res.bHasRes)
	{
		tl::log_err("Cannot calculate resolution matrix.");
		return false;
	}

	return true;
}


static std::unique_ptr<EllipseDlg> show_ellipses(const Resolution& res)
{
	std::unique_ptr<EllipseDlg> ellidlg =
		std::make_unique<EllipseDlg>(nullptr, nullptr, Qt::Window);

	EllipseDlgParams params{};
	params.algo = ResoAlgo::MC;
	params.reso = &res.res;
	params.vecMC_direct = &res.vecQ;
	ellidlg->SetParams(params);

	focus_dlg(ellidlg.get());
	return ellidlg;
}


int montereso_main(int argc, char **argv)
{
	tl::init_rand();

	std::ios_base::sync_with_stdio(0);
	std::setlocale(LC_ALL, "C");


	ublas::vector<t_real>
		vecQPara = tl::zero_v<ublas::vector<t_real>>(3),
		vecQPerp= tl::zero_v<ublas::vector<t_real>>(3);
	std::string strOrient1, strOrient2;
	std::string strFile;
	bool bReso = false;
	bool bCovar = false;
	bool bSwapYZ = false;


	// skip dummy argument when starting from takin
	int start_arg = 0;
	if(argc > 1 && std::string(argv[1]) == "--montereso")
		++start_arg;

	opts::options_description args("program options");
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("in-file",
		opts::value<decltype(strFile)>(&strFile), "input file")));
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("reso",
		opts::bool_switch(&bReso),
		"file contains the resolution matrix")));
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("covar",
		opts::bool_switch(&bCovar),
		"file contains the covariance matrix")));
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("swapYZ",
		opts::bool_switch(&bSwapYZ),
		"swap event y and z coordinate components")));
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("orient1",
		opts::value<decltype(strOrient1)>(&strOrient1),
		"first orientation vector")));
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("orient2",
		opts::value<decltype(strOrient2)>(&strOrient2),
		"second orientation vector")));

	// dummy arg if launched from takin executable
	bool bStartedFromTakin = false;
#ifndef MONTERESO_STANDALONE
	args.add(boost::shared_ptr<opts::option_description>(
		new opts::option_description("montereso",
		opts::bool_switch(&bStartedFromTakin),
		"launch montereso from takin")));
#endif

	opts::positional_options_description args_pos;
	args_pos.add("in-file", -1);

	opts::basic_command_line_parser<char> clparser(argc, argv);
	clparser.options(args);
	clparser.positional(args_pos);
	opts::basic_parsed_options<char> parsedopts = clparser.run();

	opts::variables_map opts_map;
	opts::store(parsedopts, opts_map);
	opts::notify(opts_map);


	std::vector<t_real> _vecOrient1, _vecOrient2;
	if(strOrient1 != "")
		tl::get_tokens<t_real>(strOrient1, std::string(" ,;"), _vecOrient1);
	if(strOrient2 != "")
		tl::get_tokens<t_real>(strOrient2, std::string(" ,;"), _vecOrient2);

	for(std::size_t i=0; i<std::min(vecQPara.size(), _vecOrient1.size()); ++i)
		vecQPara[i] = _vecOrient1[i];
	for(std::size_t i=0; i<std::min(vecQPerp.size(), _vecOrient2.size()); ++i)
		vecQPerp[i] = _vecOrient2[i];


	if(argc < 2 + start_arg)
	{
		std::cerr << args << std::endl;
		return -1;
	}


	FileType ft = FileType::UNKNOWN;
	if(bReso)
		ft = FileType::RESOLUTION_MATRIX;
	else if(bCovar)
		ft = FileType::COVARIANCE_MATRIX;


	Resolution res;

	if(ft == FileType::RESOLUTION_MATRIX || ft == FileType::COVARIANCE_MATRIX)
	{
		tl::log_info("Loading covariance/resolution matrix from \"", strFile, "\".");
		if(!load_mat(strFile.c_str(), res, ft))
			return -1;
	}
	else
	{
		tl::log_info("Loading neutron list from \"", strFile, "\".");
		if(!load_mc_list(strFile.c_str(), res,
			_vecOrient1.size() ? &vecQPara : nullptr,
			_vecOrient2.size() ? &vecQPerp : nullptr,
			bSwapYZ))
			return -1;
	}


	QApplication app(argc, argv);
	app.setApplicationName("Takin/Montereso");
	app.setApplicationVersion(TAKIN_VER);
	app.setQuitOnLastWindowClosed(true);

	QLocale::setDefault(QLocale::English);

	std::unique_ptr<EllipseDlg> ellidlg{show_ellipses(res)};
	return app.exec();
}
