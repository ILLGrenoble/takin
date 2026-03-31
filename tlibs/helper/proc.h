/**
 * process helpers
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: feb-2017
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

#ifndef __TLIBS_PROC_H__
#define __TLIBS_PROC_H__

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <memory>


namespace tl {


template<class t_ch = char>
class PipeProc
{
protected:
	bool m_bOk = false;
	FILE *m_pipe = nullptr;
	bool m_bWrite = true;

	std::unique_ptr<boost::iostreams::file_descriptor_sink> m_pfdsOut;
	std::unique_ptr<boost::iostreams::file_descriptor_source> m_pfdsIn;

	std::unique_ptr<boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_sink>> m_psbufOut;
	std::unique_ptr<boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_source>> m_psbufIn;

	std::unique_ptr<std::basic_ostream<t_ch>> m_postr;
	std::unique_ptr<std::basic_istream<t_ch>> m_pistr;


public:
	PipeProc(const char* pcProc, bool bWrite = true) : m_bOk{false}, m_bWrite{bWrite}
	{
		namespace ios = boost::iostreams;

		m_pipe = (FILE*)::/*my_*/popen(pcProc, bWrite ? "w" : "r");
		int fileno = ::fileno(m_pipe);
		if(!m_pipe || fileno == -1)
			return;

		if(bWrite)
		{
			m_pfdsOut.reset(new ios::file_descriptor_sink(fileno, ios::close_handle));
			if(m_pfdsOut && m_pfdsOut->is_open())
				m_psbufOut.reset(new ios::stream_buffer<ios::file_descriptor_sink>(*m_pfdsOut));
			if(m_psbufOut && m_psbufOut->is_open())
				m_postr.reset(new std::basic_ostream<t_ch>(m_psbufOut.get()));
			if(m_postr && m_postr->good())
				m_bOk = true;
		}
		else
		{
			m_pfdsIn.reset(new ios::file_descriptor_source(fileno, ios::close_handle));
			if(m_pfdsIn && m_pfdsIn->is_open())
				m_psbufIn.reset(new ios::stream_buffer<ios::file_descriptor_source>(*m_pfdsIn));
			if(m_psbufIn && m_psbufIn->is_open())
				m_pistr.reset(new std::basic_istream<t_ch>(m_psbufIn.get()));
			if(m_pistr && m_pistr->good())
				m_bOk = true;
		}
	}


	virtual ~PipeProc()
	{
		m_postr.reset();
		m_pistr.reset();

		m_psbufOut.reset();
		m_psbufIn.reset();

		m_pfdsOut.reset();
		m_pfdsIn.reset();

		if(m_pipe) ::/*my_*/pclose(m_pipe);
		m_pipe = nullptr;
	}


	bool IsReady() const { return m_bOk; }
	bool IsWritePipe() const { return m_bWrite; }


	std::basic_ostream<t_ch>& GetOstr() { return *m_postr; }
	std::basic_istream<t_ch>& GetIstr() { return *m_pistr; }


	void flush()
	{
		if(m_postr) m_postr->flush();
	}
};

}


template<class t_ch, class T>
std::basic_ostream<t_ch>& operator<<(tl::PipeProc<t_ch>& proc, const T& t)
{
	return (proc.GetOstr() << t);
}


template<class t_ch, class T>
std::basic_istream<t_ch>& operator>>(T& t, tl::PipeProc<t_ch>& proc)
{
	return (t >> proc.GetIstr());
}


#endif
