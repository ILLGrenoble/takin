/**
 * S(Q,w) processes
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date sep-2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __SQW_PROC_IMPL_H__
#define __SQW_PROC_IMPL_H__

#include "sqw_proc.h"
#include "tlibs/string/string.h"
#include "tlibs/file/file.h"
#include "tlibs/log/log.h"
#include "tlibs/math/rand.h"

#include <signal.h>
#include <unistd.h>
#include <errno.h>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/string.hpp>

#define MSG_QUEUE_SIZE     512
#define PARAM_MEM          1024*1024
#define WAIT_END_PROCESSES 250


namespace ipr = boost::interprocess;
using t_real = t_real_reso;


// ----------------------------------------------------------------------------
// types and conversions

template<class t_ch=char>
using t_sh_str_alloc_gen = ipr::allocator<t_ch, ipr::managed_shared_memory::segment_manager>;
template<class t_ch=char>
using t_sh_str_gen = ipr::basic_string<t_ch, std::char_traits<t_ch>, t_sh_str_alloc_gen<t_ch>>;

using t_sh_str_alloc = t_sh_str_alloc_gen<char>;
using t_sh_str = t_sh_str_gen<char>;


/**
 * converts the energy and weight vectors of the dispersion to a string
 */
static void disp_to_str(
	t_sh_str& str, const std::tuple<std::vector<t_real>, std::vector<t_real>>& tup)
{
	try
	{
		str.clear();
		std::string strVecE = tl::cont_to_str<std::vector<t_real>>(std::get<0>(tup), ",");
		std::string strVecW = tl::cont_to_str<std::vector<t_real>>(std::get<1>(tup), ",");

		std::size_t iLenTotal = strVecE.length() + strVecW.length() + 3 + 1;
		if(iLenTotal >= std::size_t(PARAM_MEM))
		{
			tl::log_err("Process buffer limit reached. Cannot proceed.");
			return;
		}

		str.append(strVecE.c_str()); str.append("#,#");
		str.append(strVecW.c_str());
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


/**
 * converts a string to the energy and weight vectors of the dispersion
 */
static std::tuple<std::vector<t_real>, std::vector<t_real>>
str_to_disp(const t_sh_str& _str)
{
	std::vector<t_real> vecE, vecW;

	try
	{
		std::string str;
		for(const t_sh_str::value_type& ch : _str) str.push_back(ch);

		std::string strVecE, strVecW;
		std::tie(strVecE, strVecW) = tl::split_first<std::string>(str, "#,#", 1, 1);

		tl::get_tokens<t_real, std::string, decltype(vecE)>(strVecE, ",", vecE);
		tl::get_tokens<t_real, std::string, decltype(vecW)>(strVecW, ",", vecW);
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	return std::make_tuple(vecE, vecW);
}


/**
 * converts the model parameters to a string
 */
static void pars_to_str(t_sh_str& str, const std::vector<SqwBase::t_var>& vec)
{
	try
	{
		str.clear();
		//str.reserve(PARAM_MEM*0.9);

		std::size_t iLenTotal = 0;
		for(const SqwBase::t_var& var : vec)
		{
			iLenTotal += std::get<0>(var).length() +
				std::get<1>(var).length() +
				std::get<2>(var).length() + 3*3 + 1;
			if(iLenTotal >= std::size_t(PARAM_MEM*0.9))
			{
				tl::log_err("Process buffer limit imminent. Truncating parameter list.");
				break;
			}

			str.append(std::get<0>(var).c_str()); str.append("#,#");
			str.append(std::get<1>(var).c_str()); str.append("#,#");
			str.append(std::get<2>(var).c_str()); str.append("#;#");
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


/**
 * converts a string to the model parameters
 */
static std::vector<SqwBase::t_var> str_to_pars(const t_sh_str& _str)
{
	std::vector<SqwBase::t_var> vec;
	try
	{
		std::string str;
		for(const t_sh_str::value_type& ch : _str) str.push_back(ch);

		std::vector<std::string> vecLines;
		tl::get_tokens_seq<std::string, std::string, std::vector>
			(str, "#;#", vecLines, 1);


		vec.reserve(vecLines.size());

		for(const std::string& strLine : vecLines)
		{
			if(strLine.length() == 0)
				continue;

			std::vector<std::string> vecVals;
			tl::get_tokens_seq<std::string, std::string, std::vector>
				(strLine, "#,#", vecVals, 1);

			if(vecVals.size() != 3)
			{
				tl::log_err("Wrong size of parameters: \"", strLine, "\".");
				continue;
			}

			SqwBase::t_var var;
			std::get<0>(var) = vecVals[0];
			std::get<1>(var) = vecVals[1];
			std::get<2>(var) = vecVals[2];
			vec.push_back(var);
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
	return vec;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// messages

enum class ProcMsgTypes
{
	QUIT,
	NOP,

	DISP,
	SQW,
	BCK,
	GET_VARS,
	SET_VARS,

	IS_OK,
	READY,
};


struct ProcMsg
{
	ProcMsgTypes ty = ProcMsgTypes::NOP;

	t_real dParam1, dParam2, dParam3, dParam4;
	t_real dRet;
	bool bRet;
};


static void msg_send(ipr::message_queue& msgqueue, const ProcMsg& msg)
{
	try
	{
		msgqueue.send(&msg, sizeof(msg), 0);
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


static ProcMsg msg_recv(ipr::message_queue& msgqueue)
{
	ProcMsg msg;

	try
	{
		std::size_t iSize = 0;
		unsigned int iPrio = 0;
		msgqueue.receive(&msg, sizeof(msg), iSize, iPrio);

		if(iSize != sizeof(msg))
			tl::log_err("Message size mismatch.");
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}

	return msg;
}

// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// child (worker) process
// ----------------------------------------------------------------------------

template<class t_sqw>
static void child_proc(ipr::message_queue& msgToParent, ipr::message_queue& msgFromParent,
	const char* pcCfg, void* pSharedPars)
{
	std::unique_ptr<t_sqw> pSqw(new t_sqw(pcCfg));

	// tell parent that pSqw is inited
	ProcMsg msgReady;
	msgReady.ty = ProcMsgTypes::READY;
	pid_t pid = getpid();
	msgReady.dRet = *((double*)&pid);
	msgReady.bRet = pSqw->IsOk();
	msg_send(msgToParent, msgReady);

	t_sh_str *pPars = static_cast<t_sh_str*>(pSharedPars);

	// message loop
	while(true)
	{
		ProcMsg msg = msg_recv(msgFromParent);
		ProcMsg msgRet;
		bool bFailure = false;

		switch(msg.ty)
		{
			case ProcMsgTypes::DISP:	// dispersion
			{
				msgRet.ty = msg.ty;
				disp_to_str(*pPars, pSqw->disp(msg.dParam1, msg.dParam2, msg.dParam3));
				msg_send(msgToParent, msgRet);
				break;
			}
			case ProcMsgTypes::SQW:		// structure factor
			{
				msgRet.ty = msg.ty;
				msgRet.dRet = pSqw->operator()(msg.dParam1, msg.dParam2, msg.dParam3, msg.dParam4);
				msg_send(msgToParent, msgRet);
				break;
			}
			case ProcMsgTypes::BCK:		// background
			{
				msgRet.ty = msg.ty;
				msgRet.dRet = pSqw->GetBackground(msg.dParam1, msg.dParam2, msg.dParam3, msg.dParam4);
				msg_send(msgToParent, msgRet);
				break;
			}
			case ProcMsgTypes::GET_VARS:	// get variables
			{
				msgRet.ty = msg.ty;
				pars_to_str(*pPars, pSqw->GetVars());
				msg_send(msgToParent, msgRet);
				break;
			}
			case ProcMsgTypes::SET_VARS:	// set variables
			{
				pSqw->SetVars(str_to_pars(*pPars));
				msgRet.ty = ProcMsgTypes::READY;
				msgRet.bRet = true;
				msg_send(msgToParent, msgRet);
				break;
			}
			case ProcMsgTypes::IS_OK:
			{
				msgRet.ty = msg.ty;
				msgRet.bRet = pSqw->IsOk();
				msg_send(msgToParent, msgRet);

				if(!msgRet.bRet)
					bFailure = true;
				break;
			}
			case ProcMsgTypes::QUIT:
			{
				tl::log_debug("Child process ", getpid(), " received quit request.");
				return;
			}
			default:
			{
				break;
			}
		}

		if(bFailure)
			break;
	}

	tl::log_debug("Child process ", getpid(), " message loop has ended.");
}



// ----------------------------------------------------------------------------
// parent (control) process
// ----------------------------------------------------------------------------

template<class t_sqw>
SqwProc<t_sqw>::SqwProc()
{
	++m_iRefCnt;
}


/**
 * create sub-process
 */
template<class t_sqw>
SqwProc<t_sqw>::SqwProc(const char* pcCfg, SqwProcStartMode mode,
	const char* pcProcMemName, const char* pcProcExecName, unsigned int iNumChildProcesses)
	: m_iNumChildProcesses(iNumChildProcesses), m_strProcBaseName(tl::rand_name<std::string>(8))
{
	++m_iRefCnt;

	// if a process name is given (e.g. for the child process), use it
	if(pcProcMemName)
		m_strProcBaseName = pcProcMemName;

	try
	{
		if(mode == SqwProcStartMode::START_PARENT_CREATE_CHILD || mode == SqwProcStartMode::START_PARENT_FORK_CHILD)
		{
			tl::log_debug("Starting ", m_iNumChildProcesses, " child process(es).");

			// start all child processes
			for(unsigned int iChild = 0; iChild < m_iNumChildProcesses; ++iChild)
			{
				std::string strProcName = m_strProcBaseName + "_" + tl::var_to_str(iChild);
				tl::log_debug("Creating process memory \"", "takin_sqw_proc_*_", strProcName, "\".");

				m_pmtx.push_back(std::make_shared<std::mutex>());

				m_pMem.push_back(std::make_shared<ipr::managed_shared_memory>(ipr::create_only,
					("takin_sqw_proc_mem_" + strProcName).c_str(), PARAM_MEM));
				m_pSharedPars.push_back(static_cast<void*>(m_pMem[iChild]->template construct<t_sh_str>
					(("takin_sqw_proc_params_" + strProcName).c_str())
					(t_sh_str_alloc(m_pMem[iChild]->get_segment_manager()))));

				m_pmsgIn.push_back(std::make_shared<ipr::message_queue>(ipr::create_only,
					("takin_sqw_proc_in_" + strProcName).c_str(), MSG_QUEUE_SIZE, sizeof(ProcMsg)));
				m_pmsgOut.push_back(std::make_shared<ipr::message_queue>(ipr::create_only,
					("takin_sqw_proc_out_" + strProcName).c_str(), MSG_QUEUE_SIZE, sizeof(ProcMsg)));

				// create a child process
				if(mode == SqwProcStartMode::START_PARENT_CREATE_CHILD)
				{
					if(!tl::file_exists(pcProcExecName))
					{
						tl::log_err("Child process file \"", pcProcExecName, "\" does not exist.");
						return;
					}

					// start child process
					if(std::system((std::string(pcProcExecName)
						+ std::string(" \"") + pcCfg + std::string("\" ")
						+ strProcName + " &").c_str()) < 0)
					{
						const int errnum = errno;
						tl::log_err("Could not create child process \"", pcProcExecName, "\".",
							" Error code: ", errnum, ".");
					}
				}

#ifndef __MINGW32__
				// fork a child process from the parent process
				else if(mode == SqwProcStartMode::START_PARENT_FORK_CHILD)
				{
					pid_t pidChild = fork();

					if(pidChild < 0)
					{
						tl::log_err("Could not fork child process.");
						return;
					}
					else if(pidChild == 0)
					{
						// start of child process
						m_iNumChildProcesses = 1;
						m_pidChild.resize(1);
						m_pidChild[0] = 0;

						child_proc<t_sqw>(*m_pmsgIn[iChild], *m_pmsgOut[iChild], pcCfg, m_pSharedPars[iChild]);

						exit(0);
						// end of child process
						return;
					}

					// in control process
					m_pidChild.push_back(pidChild);
				}
#endif

				tl::log_debug("Waiting for child process ", iChild, " to become ready...");

				ProcMsg msgReady = msg_recv(*m_pmsgIn[iChild]);
				if(mode == SqwProcStartMode::START_PARENT_CREATE_CHILD)
					m_pidChild.push_back(*((pid_t*)&msgReady.dRet));
				m_bOk = msgReady.bRet;

				if(!m_bOk)
				{
					tl::log_err("Child process ", m_pidChild[iChild], " reports failure.");
					m_iNumChildProcesses = iChild+1;
					break;
				}
				else
				{
					tl::log_debug("Child process ", m_pidChild[iChild], " is ready.");
				}
			}
		}

		// this process has already been started as a child process
		else if(mode == SqwProcStartMode::START_CHILD)
		{
			// for the child process, the vectors have only one element
			m_iNumChildProcesses = 1;
			m_pidChild.resize(1);
			m_pidChild[0] = 0;

			const std::string& strProcName = m_strProcBaseName;

			m_pMem.push_back(std::make_shared<ipr::managed_shared_memory>(ipr::open_only,
				("takin_sqw_proc_mem_" + strProcName).c_str()));
			m_pSharedPars.push_back(static_cast<void*>(m_pMem[0]->template find<t_sh_str>
				(("takin_sqw_proc_params_" + strProcName).c_str()).first));

			m_pmsgIn.push_back(std::make_shared<ipr::message_queue>(ipr::open_only,
				("takin_sqw_proc_in_" + strProcName).c_str()));
			m_pmsgOut.push_back(std::make_shared<ipr::message_queue>(ipr::open_only,
				("takin_sqw_proc_out_" + strProcName).c_str()));

			child_proc<t_sqw>(*m_pmsgIn[0], *m_pmsgOut[0], pcCfg, m_pSharedPars[0]);
		}
	}
	catch(const std::exception& ex)
	{
		m_bOk = false;
		tl::log_err(ex.what());
	}
}


template<class t_sqw>
SqwProc<t_sqw>::SqwProc(const std::string& strCfg) : SqwProc<t_sqw>::SqwProc(strCfg.c_str())
{
	++m_iRefCnt;
}


/**
 * clean up sub-process
 */
template<class t_sqw>
SqwProc<t_sqw>::~SqwProc()
{
	if(m_pidChild.size() == 0)
	{
		tl::log_err("No process id registered.");
		return;
	}

	// we're in a child process
	if(m_pidChild[0] == 0)
	{
		tl::log_debug("Child process ", getpid(), " ending.");
		return;
	}

	for(unsigned int iChild=0; iChild<m_iNumChildProcesses; ++iChild)
	{
		// make sure that this instance is the last
		if(m_pMem[iChild].use_count() > 1)
			return;
	}

	// shut down the parent process
	try
	{
		for(unsigned int iChild=0; iChild<m_iNumChildProcesses; ++iChild)
		{
			if(m_pmsgOut[iChild])
			{
				ProcMsg msg;
				msg.ty = ProcMsgTypes::QUIT;
				msg_send(*m_pmsgOut[iChild], msg);

				//kill(m_pidChild[iChild], SIGABRT);
			}
		}

		if(--m_iRefCnt == 0)
		{
			// give clients time to end before removing the shared memory
			std::this_thread::sleep_for(std::chrono::milliseconds{WAIT_END_PROCESSES});

			for(unsigned int iChild=0; iChild<m_iNumChildProcesses; ++iChild)
			{
				std::string strProcName = m_strProcBaseName + "_" + tl::var_to_str(iChild);

				ipr::message_queue::remove(("takin_sqw_proc_in_" + strProcName).c_str());
				ipr::message_queue::remove(("takin_sqw_proc_out_" + strProcName).c_str());

				m_pMem[iChild]->template destroy<t_sh_str>(("takin_sqw_proc_params_" + strProcName).c_str());
				ipr::shared_memory_object::remove(("takin_sqw_proc_mem_" + strProcName).c_str());

				tl::log_debug("Removed process memory \"", "takin_sqw_proc_*_",
					strProcName, "\" for child process ", m_pidChild[iChild], ".");
			}
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_debug("Child process unloading exception: ", ex.what());
	}
}


/**
 * query dispersion
 */
template<class t_sqw>
std::tuple<std::vector<t_real>, std::vector<t_real>>
SqwProc<t_sqw>::disp(t_real dh, t_real dk, t_real dl) const
{
	auto query_disp = [this](std::size_t iChild, t_real dh, t_real dk, t_real dl)
		-> std::tuple<std::vector<t_real>, std::vector<t_real>>
	{
		ProcMsg msg;
		msg.ty = ProcMsgTypes::DISP;
		msg.dParam1 = dh;
		msg.dParam2 = dk;
		msg.dParam3 = dl;
		msg_send(*m_pmsgOut[iChild], msg);

		t_sh_str *pPars = static_cast<t_sh_str*>(m_pSharedPars[iChild]);

		ProcMsg msgDisp = msg_recv(*m_pmsgIn[iChild]);
		return str_to_disp(*pPars);
	};


	if(!m_bOk)
		return std::make_tuple(std::vector<t_real>{}, std::vector<t_real>{});

	// find a free child process
	for(unsigned int iChild = 0; iChild < m_iNumChildProcesses; ++iChild)
	{
		std::unique_lock<std::mutex> lock(*m_pmtx[iChild], std::defer_lock);
		if(lock.try_lock())
			return query_disp(iChild, dh, dk, dl);
	}

	// if all processes are occupied, queue at the first one
	std::lock_guard<std::mutex> lock(*m_pmtx[0]);
	return query_disp(0, dh, dk, dl);
}


/**
 * query dynamical structure factor
 */
template<class t_sqw>
t_real SqwProc<t_sqw>::operator()(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	auto query_sqw = [this](std::size_t iChild, t_real dh, t_real dk, t_real dl, t_real dE)
		-> t_real
	{
		ProcMsg msg;
		msg.ty = ProcMsgTypes::SQW;
		msg.dParam1 = dh;
		msg.dParam2 = dk;
		msg.dParam3 = dl;
		msg.dParam4 = dE;
		msg_send(*m_pmsgOut[iChild], msg);

		ProcMsg msgS = msg_recv(*m_pmsgIn[iChild]);
		return msgS.dRet;
	};


	if(!m_bOk)
		return 0.;

	// find a free child process
	for(unsigned int iChild = 0; iChild < m_iNumChildProcesses; ++iChild)
	{
		std::unique_lock<std::mutex> lock(*m_pmtx[iChild], std::defer_lock);
		if(lock.try_lock())
			return query_sqw(iChild, dh, dk, dl, dE);
	}

	// if all processes are occupied, queue at the first one
	std::lock_guard<std::mutex> lock(*m_pmtx[0]);
	return query_sqw(0, dh, dk, dl, dE);
}


/**
 * get background
 */
template<class t_sqw>
t_real SqwProc<t_sqw>::GetBackground(t_real dh, t_real dk, t_real dl, t_real dE) const
{
	auto query_bck = [this](std::size_t iChild, t_real dh, t_real dk, t_real dl, t_real dE)
		-> t_real
	{
		ProcMsg msg;
		msg.ty = ProcMsgTypes::BCK;
		msg.dParam1 = dh;
		msg.dParam2 = dk;
		msg.dParam3 = dl;
		msg.dParam4 = dE;
		msg_send(*m_pmsgOut[iChild], msg);

		ProcMsg msgBck = msg_recv(*m_pmsgIn[iChild]);
		return msgBck.dRet;
	};


	if(!m_bOk)
		return 0.;

	// find a free child process
	for(unsigned int iChild = 0; iChild < m_iNumChildProcesses; ++iChild)
	{
		std::unique_lock<std::mutex> lock(*m_pmtx[iChild], std::defer_lock);
		if(lock.try_lock())
			return query_bck(iChild, dh, dk, dl, dE);
	}

	// if all processes are occupied, queue at the first one
	std::lock_guard<std::mutex> lock(*m_pmtx[0]);
	return query_bck(0, dh, dk, dl, dE);
}


template<class t_sqw>
bool SqwProc<t_sqw>::IsOk() const
{
	if(!m_bOk)
		return false;

	// check all sub-processes
	for(unsigned int iChild=0; iChild<m_iNumChildProcesses; ++iChild)
	{
		std::lock_guard<std::mutex> lock(*m_pmtx[iChild]);

		ProcMsg msg;
		msg.ty = ProcMsgTypes::IS_OK;
		msg_send(*m_pmsgOut[iChild], msg);

		ProcMsg msgRet = msg_recv(*m_pmsgIn[iChild]);
		if(!msgRet.bRet)
			return false;
	}

	return true;
}


/**
 * query variables
 */
template<class t_sqw>
std::vector<SqwBase::t_var> SqwProc<t_sqw>::GetVars() const
{
	if(!m_bOk)
		return std::vector<SqwBase::t_var>{};

	// the variables should be the same for all sub-processes => only query the first
	std::lock_guard<std::mutex> lock(*m_pmtx[0]);

	ProcMsg msg;
	msg.ty = ProcMsgTypes::GET_VARS;
	msg_send(*m_pmsgOut[0], msg);

	t_sh_str *pPars = static_cast<t_sh_str*>(m_pSharedPars[0]);

	ProcMsg msgRet = msg_recv(*m_pmsgIn[0]);
	return str_to_pars(*pPars);
}


/**
 * set variables
 */
template<class t_sqw>
void SqwProc<t_sqw>::SetVars(const std::vector<SqwBase::t_var>& vecVars)
{
	if(!m_bOk)
		return;

	// set the same variables for all child-processes
	for(unsigned int iChild=0; iChild<m_iNumChildProcesses; ++iChild)
	{
		std::lock_guard<std::mutex> lock(*m_pmtx[iChild]);

		t_sh_str *pPars = static_cast<t_sh_str*>(m_pSharedPars[iChild]);

		ProcMsg msg;
		msg.ty = ProcMsgTypes::SET_VARS;
		pars_to_str(*pPars, vecVars);
		msg_send(*m_pmsgOut[iChild], msg);

		ProcMsg msgRet = msg_recv(*m_pmsgIn[iChild]);
		if(!msgRet.bRet)
			tl::log_err("Could not set variables for child process ", iChild, ".");
	}
}


template<class t_sqw>
SqwBase* SqwProc<t_sqw>::shallow_copy() const
{
	SqwProc* pSqw = new SqwProc();
	*static_cast<SqwBase*>(pSqw) = *static_cast<const SqwBase*>(this);

	pSqw->m_pmtx = this->m_pmtx;
	pSqw->m_iNumChildProcesses = this->m_iNumChildProcesses;
	pSqw->m_pMem = this->m_pMem;
	pSqw->m_pmsgIn = this->m_pmsgIn;
	pSqw->m_pmsgOut = this->m_pmsgOut;
	pSqw->m_strProcBaseName = this->m_strProcBaseName;
	pSqw->m_pidChild = this->m_pidChild;
	pSqw->m_pSharedPars = this->m_pSharedPars;
	pSqw->m_iRefCnt = this->m_iRefCnt;

	return pSqw;
}

// ----------------------------------------------------------------------------

#endif
