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

#include "log.h"
#include <iomanip>
#include <cstring>
#include <chrono>
#include <boost/date_time/c_time.hpp>


namespace tl {
std::recursive_mutex Log::s_mtx;
bool Log::s_bTermCmds = 1;


std::string Log::get_timestamp()
{
	namespace ch = std::chrono;
	using ch::system_clock;
	using boost::date_time::c_time;

	system_clock::time_point now = system_clock::now();

	// milliseconds
	ch::milliseconds msecs = ch::duration_cast<ch::milliseconds>(now.time_since_epoch());
	auto secs = ch::duration_cast<ch::seconds>(msecs);
	msecs -= ch::duration_cast<ch::milliseconds>(secs);
	std::ostringstream ostrmsecs;
	ostrmsecs << std::setw(3) << std::right << std::setfill('0') << msecs.count();

	// time and date
	std::time_t tm = system_clock::to_time_t(now);
	std::tm tmNow;
	c_time::localtime(&tm, &tmNow);

	std::string::value_type cTime[64];
	std::strftime(cTime, sizeof cTime, "%Y-%b-%d %H:%M:%S.", &tmNow);
	if(std::strlen(cTime))
	{	// if a time string is available, return it
		std::strcat(cTime, ostrmsecs.str().c_str());
		return std::string(cTime);
	}
	else
	{	// else return the raw seconds
		std::ostringstream ostrsecs;
		ostrsecs << secs.count();
		ostrsecs << "." << ostrmsecs.str();
		return ostrsecs.str();
	}
}


std::string Log::get_thread_id()
{
	std::ostringstream ostr;
	ostr << std::hex << std::this_thread::get_id();
	return ostr.str();
}


std::string Log::get_color(LogColor col, bool bBold)
{
	if(!s_bTermCmds) return "";

	switch(col)
	{
		case LogColor::RED: return bBold ? "\033[1;31m" : "\033[0;31m";
		case LogColor::YELLOW: return bBold ? "\033[1;33m" : "\033[0;33m";
		case LogColor::BLUE: return bBold ? "\033[1;34m" : "\033[0;34m";
		case LogColor::GREEN: return bBold ? "\033[1;32m" : "\033[0;32m";
		case LogColor::PURPLE: return bBold ? "\033[1;35m" : "\033[0;35m";
		case LogColor::CYAN: return bBold ? "\033[1;36m" : "\033[0;36m";
		case LogColor::WHITE: return bBold ? "\033[1;37m" : "\033[0;37m";
		case LogColor::BLACK: return bBold ? "\033[1;30m" : "\033[0;30m";
		case LogColor::NONE:
		default: return "\033[0m";
	}
}


void Log::begin_log()
{
	s_mtx.lock();

	const std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
	std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

	for(t_pairOstr &pairOstr : vecOstrs)
	{
		std::ostream *pOstr = pairOstr.first;
		bool bCol = pairOstr.second;

		if(!pOstr)
			continue;
		if(bCol)
			(*pOstr) << get_color(m_col, 1);
		if(m_bShowDate)
		{
			std::string strTimeStamp = get_timestamp();
			if(strTimeStamp != "")
				(*pOstr) << strTimeStamp << ", ";
		}
		if(m_bShowThread)
		{
			using t_threadmap = std::unordered_map<std::thread::id, std::string>;
			static t_threadmap s_threadmap;

			using t_mapKey = typename t_threadmap::value_type::first_type;

			t_mapKey idThread = std::this_thread::get_id();
			typename t_threadmap::const_iterator iterMap = s_threadmap.find(idThread);
			if(iterMap == s_threadmap.end())
			{
				++m_iNumThreads;
				std::ostringstream ostrThread;
				ostrThread << "Thread " << m_iNumThreads;

				iterMap = s_threadmap.insert({idThread, ostrThread.str()}).first;
			}

			(*pOstr) << iterMap->second << ", ";
		}
		(*pOstr) << m_strInfo << ": ";
		if(bCol)
			(*pOstr) << get_color(m_col, 0);
	}
}


void Log::end_log()
{
	const std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
	std::vector<t_pairOstr> vecOstrs = arrayunion({m_vecOstrs, vecOstrsTh});

	for(t_pairOstr& pairOstr : vecOstrs)
	{
		std::ostream *pOstr = pairOstr.first;
		bool bCol = pairOstr.second;

		if(!pOstr)
			continue;
		if(bCol)
			(*pOstr) << get_color(LogColor::NONE);
		(*pOstr) << std::endl;
	}

	s_mtx.unlock();
}


void Log::inc_depth()
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);

	if(m_iDepth++ == 0)
		begin_log();
}


void Log::dec_depth()
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	if(--m_iDepth <= 0)
	{
		m_iDepth = 0;
		end_log();
	}
}


Log::Log() : m_vecOstrs{{&std::cerr, 1}}
{}


Log::Log(const std::string& strInfo, LogColor col, std::ostream* pOstr)
	: m_vecOstrs{{pOstr ? pOstr : &std::cerr, 1}},
	  m_strInfo(strInfo), m_col(col)
{}


Log::~Log()
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	//std::cerr << "Removing " << m_strInfo << " logger for thread " << std::hex << std::this_thread::get_id() << std::endl;
	//std::cout << boost::stacktrace::stacktrace{} << std::endl;

	m_mapOstrsTh.clear();
	m_vecOstrs.clear();
}


std::vector<Log::t_pairOstr>& Log::GetThreadOstrs()
{
	t_mapthreadOstrs::iterator iter = m_mapOstrsTh.find(std::this_thread::get_id());
	if(iter == m_mapOstrsTh.end())
		iter = m_mapOstrsTh.insert({std::this_thread::get_id(), std::vector<t_pairOstr>()}).first;

	return iter->second;
}


void Log::AddOstr(std::ostream* pOstr, bool bCol, bool bThreadLocal)
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	if(bThreadLocal)
	{
		std::vector<t_pairOstr>& vecOstrsTh = GetThreadOstrs();
		vecOstrsTh.push_back({pOstr, bCol});
	}
	else
	{
		m_vecOstrs.push_back({pOstr, bCol});
	}
}


void Log::RemoveOstr(std::ostream* pOstr)
{
	std::lock_guard<decltype(s_mtx)> _lck(s_mtx);
	using t_iter = std::vector<t_pairOstr>::iterator;

	auto fktDel = [pOstr](const t_iter::value_type& pairOstr) -> bool
	{
		return pairOstr.first == pOstr;
	};

	t_iter iterNewEnd = std::remove_if(m_vecOstrs.begin(), m_vecOstrs.end(), fktDel);
	if(iterNewEnd != m_vecOstrs.end())
		m_vecOstrs.resize(iterNewEnd-m_vecOstrs.begin());

	for(t_mapthreadOstrs::value_type& pairTh : m_mapOstrsTh)
	{
		t_iter iterNewEndTh = std::remove_if(pairTh.second.begin(), pairTh.second.end(), fktDel);
		if(iterNewEndTh != pairTh.second.end())
			pairTh.second.resize(iterNewEndTh - pairTh.second.begin());
	}
}


// use -fvisibility=hidden to avoid multiple calls to the destructors in loaded external libraries
Log log_info("INFO", LogColor::WHITE, &std::cerr),
	log_warn("WARNING", LogColor::YELLOW, &std::cerr),
	log_err("ERROR", LogColor::RED, &std::cerr),
	log_crit("CRITICAL", LogColor::PURPLE, &std::cerr),
	log_debug("DEBUG", LogColor::CYAN, &std::cerr);
}
