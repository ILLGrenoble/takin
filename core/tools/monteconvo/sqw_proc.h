/**
 * delegates S(Q,w) functions through processes
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

#ifndef __SQW_PROC_H__
#define __SQW_PROC_H__

// use filesystem based shared memory emulation
#ifdef __USE_SQW_INTERPROC_EMUL__
	#include <boost/interprocess/detail/workaround.hpp>

	// try to avoid certain shared memory syscalls
	#undef BOOST_INTERPROCESS_XSI_SHARED_MEMORY_OBJECTS
	#undef BOOST_INTERPROCESS_POSIX_SHARED_MEMORY_OBJECTS
	#undef BOOST_INTERPROCESS_FILESYSTEM_BASED_POSIX_SHARED_MEMORY
	#undef BOOST_INTERPROCESS_RUNTIME_FILESYSTEM_BASED_POSIX_SHARED_MEMORY

	// force using emulation mode instead:
	//   https://www.boost.org/doc/libs/1_75_0/doc/html/interprocess/sharedmemorybetweenprocesses.html#interprocess.sharedmemorybetweenprocesses.sharedmemory.emulation
	#define BOOST_INTERPROCESS_FORCE_GENERIC_EMULATION
	//#define BOOST_INTERPROCESS_SHARED_DIR_PATH "/tmp/shared_mem"
#endif

// disable xsi syscalls (shmat(), shmctl(), shmdt())
#ifdef __DISABLE_SQW_INTERPROC_XSI__
	#include <boost/interprocess/detail/workaround.hpp>

	// try to avoid certain shared memory syscalls
	#undef BOOST_INTERPROCESS_XSI_SHARED_MEMORY_OBJECTS
#endif

#include "sqwbase.h"
#include <mutex>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>
#include <boost/interprocess/ipc/message_queue.hpp>


enum class SqwProcStartMode
{
	// run the parent process and create child processes
	START_PARENT_CREATE_CHILD,

	// run the parent process and fork child processes
	START_PARENT_FORK_CHILD,

	// run the child process
	START_CHILD
};


template<class t_sqw>
class SqwProc : public SqwBase
{
private:
	std::size_t m_iRefCnt = 0;


protected:
	mutable std::vector<std::shared_ptr<std::mutex>> m_pmtx;

	std::size_t m_iNumChildProcesses = 1;
	std::string m_strProcBaseName;
	std::vector<pid_t> m_pidChild;

	std::vector<std::shared_ptr<boost::interprocess::managed_shared_memory>> m_pMem;
	std::vector<std::shared_ptr<boost::interprocess::message_queue>> m_pmsgIn, m_pmsgOut;
	std::vector<void*> m_pSharedPars;


public:
	SqwProc();
	SqwProc(const char* pcCfg, SqwProcStartMode mode=SqwProcStartMode::START_PARENT_FORK_CHILD,
		const char* pcProcMemName = nullptr, const char* pcProcExecName = nullptr,
		unsigned int iNumChildProcesses = 1);
	explicit SqwProc(const std::string& strCfg);
	virtual ~SqwProc();

	virtual std::tuple<std::vector<t_real_reso>, std::vector<t_real_reso>>
		disp(t_real_reso dh, t_real_reso dk, t_real_reso dl) const override;
	virtual t_real_reso
		operator()(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;
	virtual t_real_reso
		GetBackground(t_real_reso dh, t_real_reso dk, t_real_reso dl, t_real_reso dE) const override;
	virtual bool IsOk() const override;

	virtual std::vector<SqwBase::t_var> GetVars() const override;
	virtual void SetVars(const std::vector<SqwBase::t_var>&) override;

	virtual SqwBase* shallow_copy() const override;
};

#endif
