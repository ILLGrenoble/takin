/**
 * Py helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 27-aug-2014, 17-feb-2015
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

#ifndef __PY_HELPERS_H__
#define __PY_HELPERS_H__

#include "../string/string.h"

namespace tl {

template<class t_str=std::string, class t_cont=std::vector<double>>
t_cont get_py_array(const t_str& str, bool need_brackets = true, bool allow_nocommas = false)
{
	typedef typename t_cont::value_type t_elems;
	t_cont vecArr;

	std::size_t iStart = str.find('[');
	std::size_t iEnd = str.find(']');

	// search for list instead
	if(iStart == t_str::npos || iEnd == t_str::npos)
	{
		iStart = str.find('(');
		iEnd = str.find(')');
	}

	// invalid array
	if(need_brackets && (iStart == t_str::npos || iEnd == t_str::npos || iEnd < iStart))
		return vecArr;

	t_str strArr = str.substr(iStart + 1, iEnd-iStart - 1);
	tl::get_tokens<t_elems, t_str>(strArr, allow_nocommas ? ", \t" : ",", vecArr);

	return vecArr;
}


template<class t_str=std::string>
t_str get_py_string(const t_str& str)
{
	std::size_t iStart = str.find_first_of("\'\"");
	std::size_t iEnd = str.find_last_of("\'\"");

	// invalid string
	if(iStart == t_str::npos || iEnd == t_str::npos || iEnd < iStart)
		return "";

	return str.substr(iStart + 1, iEnd-iStart - 1);
}

}

#endif
