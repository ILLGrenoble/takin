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

#ifndef __TL_TCP_IMPL_H__
#define __TL_TCP_IMPL_H__

#include "tcp.h"
#include "../log/log.h"
#include "../string/string.h"
#include <boost/tokenizer.hpp>

namespace tl {


template<class t_ch, class t_str>
bool TcpTxtClient<t_ch, t_str>::get_cmd_tokens(const t_str& str, const t_str& strDelim,
	std::vector<t_str>& vecStr, t_str& strRemainder)
{
	boost::char_separator<t_ch> delim(strDelim.c_str(), "", boost::keep_empty_tokens);
	boost::tokenizer<boost::char_separator<t_ch>> tok(str, delim);

	for(const t_str& strTok : tok)
	{
		vecStr.push_back(strTok);
	}

	if(vecStr.size()<=1)
		return false;

	// keep_empty_tokens leads to an empty last element if no remaining string is left
	if(*vecStr.rbegin() == "")
	{
		strRemainder = "";
		vecStr.pop_back();
	}
	else
	{
		strRemainder = *vecStr.rbegin();
		vecStr.pop_back();
	}

	return true;
}


template<class t_ch, class t_str>
TcpTxtClient<t_ch, t_str>::TcpTxtClient() : m_listWriteBuffer(1024)
{}


template<class t_ch, class t_str>
TcpTxtClient<t_ch, t_str>::~TcpTxtClient()
{
	disconnect();

	m_sigRecv.disconnect_all_slots();
	m_sigDisconn.disconnect_all_slots();
	m_sigConn.disconnect_all_slots();
}


template<class t_ch, class t_str>
bool TcpTxtClient<t_ch, t_str>::connect(const t_str& strHost, const t_str& strService)
{
	m_strHost = strHost;
	m_strService = strService;

	try
	{
		disconnect();

#if BOOST_VERSION >= 108700
		m_pservice = new asio::io_context;
#else
		m_pservice = new asio::io_service;
#endif
		m_psock = new ip::tcp::socket(*m_pservice);
		ip::tcp::resolver res(*m_pservice);

#if BOOST_VERSION >= 108700
		using t_iter = typename ip::basic_resolver_results<ip::tcp>::iterator;
		using t_results = typename ip::tcp::resolver::results_type;
		t_results results = res.resolve(strHost, strService);

		asio::async_connect(*m_psock, results.begin(), results.end(),
#else
		using t_iter = typename ip::tcp::resolver::iterator;
		t_iter iter = res.resolve({strHost, strService});

		asio::async_connect(*m_psock, iter,
#endif
		[this](const sys::error_code& err, t_iter iter)
		{
			if(!err)
			{
				read_loop();
				m_sigConn(m_strHost, m_strService);
			}
			else
			{
				log_err("TCP connection error.",
					" Category: ", err.category().name(),
					", message: ", err.message(), ".");
			}
		});

		//tl::log_debug("Starting service thread.");
		m_pthread = new std::thread([this]()
		{
			try
			{
				m_pservice->run();
			}
			catch(const std::exception& ex)
			{
				log_err("TCP client thread exited with error: ", ex.what(), ".");
				m_pthread = nullptr;
				disconnect(1);
			}
		});
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
		return 0;
	}

	return 1;
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::disconnect(bool bAlwaysSendSignal)
{
	const bool bConnected = is_connected();
	if(bConnected)
	{
		m_psock->shutdown(ip::tcp::socket::shutdown_send);
		m_pservice->stop();
		m_psock->close();
	}

	if(m_psock) { delete m_psock; m_psock = 0; }
	if(m_pthread)
	{
		m_pthread->join();
		delete m_pthread;
		m_pthread = 0;
	}
	if(m_pservice) { delete m_pservice; m_pservice = 0; }

	if(bConnected || bAlwaysSendSignal)
	{
		m_sigDisconn(m_strHost, m_strService);
		m_strHost = m_strService = "";
	}

	// clean up write buffer
	const t_str* pstr = nullptr;
	while(m_listWriteBuffer.pop(pstr))
	{
		if(pstr) { delete pstr; pstr = nullptr; }
	}
}


template<class t_ch, class t_str>
bool TcpTxtClient<t_ch, t_str>::is_connected()
{
	if(!m_psock) return 0;
	return m_psock->is_open();
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::wait()
{
	if(m_pthread)
		m_pthread->join();
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::write(const t_str& str)
{
	try
	{
		if(!is_connected())
		{
			log_err("Not connected, cannot write to socket.");
			disconnect();
			return;
		}

		m_listWriteBuffer.push(new t_str(str));
#if BOOST_VERSION >= 108700
		boost::asio::post(*m_pservice, [&](){ flush_write(); });
#else
		m_pservice->post([&](){ flush_write(); });
#endif
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
	}
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::flush_write()
{
	const t_str* pstr = nullptr;
	if(m_listWriteBuffer.empty()) return;
	if(!m_listWriteBuffer.pop(pstr)) return;
	if(!pstr) return;

	asio::async_write(*m_psock, asio::buffer(pstr->data(), pstr->length()),
	[this, pstr](const sys::error_code& err, std::size_t len)
	{
		if(pstr) delete pstr;

		if(err)
		{
			disconnect();
			return;
		}

		// call this function again if buffer is not yet empty
		if(!m_listWriteBuffer.empty())
			flush_write();
	});
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::read_loop()
{
	//tl::log_debug("In read loop.");
	asio::async_read(*m_psock, asio::buffer(m_pcReadBuffer, m_iReadBufLen), asio::transfer_at_least(1),
	[this](const sys::error_code& err, std::size_t len)
	{
		if(err)
		{
			log_err("TCP read error. Category: ", err.category().name(),
				", message: ", err.message(), ".");
			disconnect();
			return;
		}

		t_str strCurMsg(m_pcReadBuffer, len);
		m_strReadBuffer.append(strCurMsg);

		std::vector<t_str> vecCmds;
		if(get_cmd_tokens(m_strReadBuffer, m_strCmdDelim, vecCmds, m_strReadBuffer))
		{
			//tl::log_debug("remainder: ", m_strReadBuffer);
			for(const t_str& strCmd : vecCmds)
			{
				m_sigRecv(strCmd);
			}
		}

		// call this function again
		read_loop();
	});
  }


// --------------------------------------------------------------------------------
// Signals
template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::add_receiver(const typename t_sigRecv::slot_type& conn)
{
	m_sigRecv.connect(conn);
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::add_disconnect(const typename t_sigDisconn::slot_type& disconn)
{
	m_sigDisconn.connect(disconn);
}


template<class t_ch, class t_str>
void TcpTxtClient<t_ch, t_str>::add_connect(const typename t_sigConn::slot_type& conn)
{
	m_sigConn.connect(conn);
}
// --------------------------------------------------------------------------------



// ================================================================================



template<class t_ch, class t_str>
TcpTxtServer<t_ch, t_str>::TcpTxtServer()
	: TcpTxtClient<t_ch, t_str>()
{}


template<class t_ch, class t_str>
TcpTxtServer<t_ch, t_str>::~TcpTxtServer()
{}


template<class t_ch, class t_str>
void TcpTxtServer<t_ch, t_str>::disconnect(bool bAlwaysSendSignal)
{
	//if(!this->is_connected()) { return; }
	if(m_pacceptor) { delete m_pacceptor; m_pacceptor = nullptr; }
	if(m_pendpoint) { delete m_pendpoint; m_pendpoint = nullptr; }

	TcpTxtClient<t_ch, t_str>::disconnect(bAlwaysSendSignal);
}


template<class t_ch, class t_str>
bool TcpTxtServer<t_ch, t_str>::start_server(unsigned short iPort)
{
	this->m_strHost = "localhost";
	this->m_strService = tl::var_to_str(iPort);

	try
	{
		disconnect();

#if BOOST_VERSION >= 108700
		this->m_pservice = new asio::io_context;
#else
		this->m_pservice = new asio::io_service;
#endif
		this->m_psock = new ip::tcp::socket(*this->m_pservice);
		m_pendpoint = new ip::tcp::endpoint(ip::tcp::v4(), iPort);
		m_pacceptor = new ip::tcp::acceptor(*this->m_pservice, *m_pendpoint);

		m_pacceptor->listen();
		m_pacceptor->async_accept(*this->m_psock,
		[this, iPort](const sys::error_code& err)
		{
			if(!err)
			{
				this->read_loop();
				m_sigServerStart(iPort);
			}
			else
			{
				log_err("TCP server error.",
					" Category: ", err.category().name(),
					", message: ", err.message(), ".");
			}
		});

		this->m_pthread = new std::thread([this]()
		{
			try
			{
				this->m_pservice->run();
			}
			catch(const std::exception& ex)
			{
				log_err("TCP server thread exited with error: ", ex.what(), ".");
				this->m_pthread = nullptr;
				disconnect(1);
			}
		});
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
		return 0;
	}
	return 1;
}


// --------------------------------------------------------------------------------
// Signals
template<class t_ch, class t_str>
void TcpTxtServer<t_ch, t_str>::add_server_start(const typename t_sigServerStart::slot_type& conn)
{
	m_sigServerStart.connect(conn);
}
// --------------------------------------------------------------------------------


}
#endif
