/**
 * Convolution fitting
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
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

#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>

#include "convofit.h"
#include "convofit_cli.h"
#include "../monteconvo/sqwfactory.h"

#include "libs/version.h"
#include "libs/globals.h"

#include "tlibs/log/log.h"
#include "tlibs/log/debug.h"
#include "tlibs/time/stopwatch.h"
#include "tlibs/helper/thread.h"

namespace opts = boost::program_options;

using t_real = t_real_reso;


// ----------------------------------------------------------------------------
// main program

int convofit_main(int argc, char** argv)
{
	try
	{
#ifdef MONTECONVO_STANDALONE	// only show copyright banner if not already displayed from Takin main program
		tl::log_info("--------------------------------------------------------------------------------");
		tl::log_info("This is the Takin command-line convolution fitter (convofit), version " TAKIN_VER ".");
		tl::log_info("Written by Tobias Weber <tweber@ill.fr>, 2014 - 2024.");
		tl::log_info(TAKIN_LICENSE("Takin/Convofit"));
		tl::log_debug("Resolution calculation uses ", sizeof(t_real_reso)*8, " bit ", tl::get_typename<t_real_reso>(), "s.");
		tl::log_debug("Fitting uses ", sizeof(tl::t_real_min)*8, " bit ", tl::get_typename<tl::t_real_min>(), "s.");
		tl::log_info("--------------------------------------------------------------------------------");
#endif

		load_sqw_plugins();


		// --------------------------------------------------------------------
		// get job files and program options
		std::vector<std::string> vecJobs;

		// normal args
		opts::options_description args("convofit options (overriding job file settings)");
		args.add(boost::make_shared<opts::option_description>(
			"job-file", opts::value<decltype(vecJobs)>(&vecJobs),
			"convolution fitting job file"));
		args.add(boost::make_shared<opts::option_description>(
			"verbose", opts::bool_switch(&g_bVerbose),
			"verbose logging"));
		args.add(boost::make_shared<opts::option_description>(
			"neutrons", opts::value<decltype(g_iNumNeutrons)>(&g_iNumNeutrons),
			"neutron count"));
		args.add(boost::make_shared<opts::option_description>(
			"skip-fit", opts::bool_switch(&g_bSkipFit),
			"skip the fitting step"));
		args.add(boost::make_shared<opts::option_description>(
			"keep-model", opts::bool_switch(&g_bUseValuesFromModel),
			"keep the initial values from the model file"));
		args.add(boost::make_shared<opts::option_description>(
			"model-params", opts::value<decltype(g_strSetParams)>(&g_strSetParams),
			"set S(Q, E) model parameters"));
		args.add(boost::make_shared<opts::option_description>(
			"outfile-suffix", opts::value<decltype(g_strOutFileSuffix)>(&g_strOutFileSuffix),
			"suffix to append to output files"));
		args.add(boost::make_shared<opts::option_description>(
			"plot-points", opts::value<decltype(g_iPlotPoints)>(&g_iPlotSkipBegin),
			"number of plot points"));
		args.add(boost::make_shared<opts::option_description>(
			"plot-skip-begin", opts::value<decltype(g_iPlotSkipBegin)>(&g_iPlotSkipBegin),
			"skip plot points in the beginning of the range"));
		args.add(boost::make_shared<opts::option_description>(
			"plot-skip-end", opts::value<decltype(g_iPlotSkipEnd)>(&g_iPlotSkipEnd),
			"skip plot points in the end of the range"));
		args.add(boost::make_shared<opts::option_description>(
			"max-threads", opts::value<decltype(g_iMaxThreads)>(&g_iMaxThreads),
			"maximum number of threads"));
		args.add(boost::make_shared<opts::option_description>(
			"max-processes", opts::value<decltype(g_iMaxProcesses)>(&g_iMaxProcesses),
			"maximum number of processes"));

		// dummy arg if launched from takin executable
		bool bStartedFromTakin = false;
#ifndef CONVOFIT_STANDALONE
		args.add(boost::make_shared<opts::option_description>(
			"convofit", opts::bool_switch(&bStartedFromTakin),
			"launch convofit from takin"));
#endif

		// positional args
		opts::positional_options_description args_pos;
		args_pos.add("job-file", -1);

		opts::basic_command_line_parser<char> clparser(argc, argv);
		clparser.options(args);
		clparser.positional(args_pos);
		opts::basic_parsed_options<char> parsedopts = clparser.run();

		opts::variables_map opts_map;
		opts::store(parsedopts, opts_map);
		opts::notify(opts_map);


		if(vecJobs.size() >= 2)
		{
			for(tl::Log* log : { &tl::log_info, &tl::log_warn, &tl::log_err, &tl::log_crit, &tl::log_debug })
				log->SetShowThread(1);
		}


		int args_to_ignore = 1;		// started with "convofit"
		if(bStartedFromTakin)
			++args_to_ignore;	// started with "takin --convofit"
		if(argc <= args_to_ignore)
		{
			std::ostringstream ostrHelp;
			ostrHelp << "Usage: ";
			for(int argidx=0; argidx<args_to_ignore; ++argidx)
				ostrHelp << argv[argidx] << " ";
			ostrHelp << "[options] <job-file 1> <job-file 2> ...\n";
			ostrHelp << args;
			tl::log_info(ostrHelp.str());
			return -1;
		}

		if(vecJobs.size() == 0)
		{
			tl::log_err("No job files given.");
			return -1;
		}
		// --------------------------------------------------------------------

		tl::Stopwatch<t_real> watch;
		watch.start();

		unsigned int iNumThreads = get_max_threads();
		tl::ThreadPool<bool()> tp(iNumThreads);

		for(std::size_t iJob=0; iJob<vecJobs.size(); ++iJob)
		{
			const std::string& strJob = vecJobs[iJob];
			tp.AddTask([iJob, strJob]() -> bool
			{
				tl::log_info("Executing job file ", iJob+1, ": \"", strJob, "\".");

				Convofit convo;
				return convo.run_job(strJob);
			});
		}

		tp.Start();
		auto& lstFut = tp.GetResults();
		std::size_t iTask = 0;
		for(auto& fut : lstFut)
		{
			bool bOk = fut.get();
			if(!bOk)
				tl::log_err("Job ", iTask+1, " (", vecJobs[iTask], ") failed or fit invalid!");
			++iTask;
		}

		watch.stop();
		tl::log_info("================================================================================");
		tl::log_info("Start time:     ", watch.GetStartTimeStr());
		tl::log_info("Stop time:      ", watch.GetStopTimeStr());
		tl::log_info("Execution time: ", tl::get_duration_str_secs<t_real>(watch.GetDur()));
		tl::log_info("================================================================================");
	}
	catch(const std::exception& ex)
	{
		tl::log_crit(ex.what());
	}

	return 0;
}
// ----------------------------------------------------------------------------
