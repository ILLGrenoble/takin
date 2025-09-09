/**
 * Simple logger
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 12-sep-2014
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <vector>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <utility>

#include "../helper/array.h"


namespace tl {

enum class LogColor
{
	NONE,
	RED, BLUE, GREEN,
	YELLOW, PURPLE, CYAN,
	WHITE, BLACK
};

class Log
{
private:
	int m_iDepth = 0;

protected:
	static std::recursive_mutex s_mtx;

	// pair of ostream and colour flag
	using t_pairOstr = std::pair<std::ostream*, bool>;
	std::vector<t_pairOstr> m_vecOstrs{};
	using t_mapthreadOstrs = std::unordered_map<std::thread::id, std::vector<t_pairOstr>>;
	t_mapthreadOstrs m_mapOstrsTh{};

	std::string m_strInfo = "";
	LogColor m_col = LogColor::NONE;

	bool m_bEnabled = 1;
	bool m_bShowDate = 1;

	bool m_bShowThread = 0;
	unsigned int m_iNumThreads = 0;

	static bool s_bTermCmds;

protected:
	static std::string get_timestamp();
	static std::string get_thread_id();
	static std::string get_color(LogColor col, bool bBold=0);

	std::vector<t_pairOstr>& GetThreadOstrs();

	void begin_log();
	void end_log();

	void inc_depth();
	void dec_depth();

public:
	Log();
	Log(const std::string& strInfo, LogColor col, std::ostream* = nullptr);
	~Log();

	void AddOstr(std::ostream* pOstr, bool bCol=1, bool bThreadLocal=0);
	void RemoveOstr(std::ostream* pOstr);

	template<typename ...t_args>
	void operator()(t_args&&... args)
	{
		if(!m_bEnabled) return;
		begin_log();

		const std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
		std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

		for(t_pairOstr& pair : vecOstrs)
		{
			if(pair.first)
				((*pair.first) << ... << std::forward<t_args>(args));
		}
		end_log();
	}

	void SetEnabled(bool bEnab) { m_bEnabled = bEnab; }
	void SetShowDate(bool bDate) { m_bShowDate = bDate; }
	void SetShowThread(bool bThread) { m_bShowThread = bThread; }

	static void SetUseTermCmds(bool bCmds) { s_bTermCmds = bCmds; }
};


extern Log log_info, log_warn, log_err, log_crit, log_debug;

}
#endif
