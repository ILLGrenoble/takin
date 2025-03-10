/**
 * TcpClient
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2014
 * @license GPLv2 or GPLv3
 *
 * @desc based on Boost's example chat client/server (c) 2003-2015 by C. M. Kohlhoff:
 * @desc   - http://www.boost.org/doc/libs/1_61_0/doc/html/boost_asio/example/cpp11/chat/chat_client.cpp
 * @desc   - https://www.boost.org/doc/libs/1_61_0/doc/html/boost_asio/example/cpp11/chat/chat_server.cpp
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

#ifndef __TCP_CLIENT__
#define __TCP_CLIENT__

#include <boost/asio.hpp>
#include <boost/signals2.hpp>
#include <boost/lockfree/queue.hpp>

#include <string>
#include <list>
#include <thread>


namespace tl {

namespace sys = boost::system;
namespace asio = boost::asio;
namespace ip = boost::asio::ip;
namespace sig = boost::signals2;
namespace lf = boost::lockfree;


template<class t_ch=char, class t_str=std::basic_string<t_ch>>
class TcpTxtClient
{
public:
	static bool get_cmd_tokens(const t_str& str, const t_str& strDelim,
		std::vector<t_str>& vecStr, t_str& strRemainder);

protected:
	t_str m_strHost, m_strService;

#if BOOST_VERSION >= 108700
	asio::io_context *m_pservice = nullptr;
#else
	asio::io_service *m_pservice = nullptr;
#endif
	ip::tcp::socket *m_psock = nullptr;
	std::thread* m_pthread = nullptr;

	t_str m_strCmdDelim = "\n";
	lf::queue<const t_str*, lf::fixed_sized<false>> m_listWriteBuffer;

	static constexpr const std::size_t m_iReadBufLen = 512;
	t_ch m_pcReadBuffer[m_iReadBufLen];
	t_str m_strReadBuffer;

	using t_sigRecv = sig::signal<void(const t_str&)>;
	using t_sigDisconn = sig::signal<void(const t_str&, const t_str&)>;
	using t_sigConn = sig::signal<void(const t_str&, const t_str&)>;

	t_sigRecv m_sigRecv;
	t_sigDisconn m_sigDisconn;
	t_sigConn m_sigConn;

public:
	TcpTxtClient();
	virtual ~TcpTxtClient();
	void set_delim(const t_str& strDelim) { m_strCmdDelim = strDelim; }

	void add_receiver(const typename t_sigRecv::slot_type& conn);
	void add_disconnect(const typename t_sigDisconn::slot_type& disconn);
	void add_connect(const typename t_sigConn::slot_type& conn);

	bool connect(const t_str& strHost, const t_str& strService);
	virtual void disconnect(bool bAlwaysSendSignal = false);
	bool is_connected();

	void write(const t_str& str);
	void wait();

protected:
	void flush_write();
	void read_loop();
};



template<class t_ch=char, class t_str=std::basic_string<t_ch>>
class TcpTxtServer : public TcpTxtClient<t_ch, t_str>
{
protected:
	ip::tcp::endpoint* m_pendpoint = nullptr;
	ip::tcp::acceptor* m_pacceptor = nullptr;

protected:
	using t_sigServerStart = sig::signal<void(unsigned short iPort)>;

	t_sigServerStart m_sigServerStart;

public:
	TcpTxtServer();
	virtual ~TcpTxtServer();

	virtual void disconnect(bool bAlwaysSendSignal = false) override;
	bool start_server(unsigned short iPort);

	void add_server_start(const typename t_sigServerStart::slot_type& start);
};

}


#ifdef TLIBS_INC_HDR_IMPLS
	#include "tcp_impl.h"
#endif
#endif
