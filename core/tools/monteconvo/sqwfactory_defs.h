/**
 * definitions for the factory and plugin interface for S(Q, E) models
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2016 -- 2026
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

 #ifndef __SQW_FACTORY_DEFS__
 #define __SQW_FACTORY_DEFS__
 

#include "sqw_proc.h"
#include "sqw_proc_impl.h"
#include "sqwrawdelegate.h"
#include "sqwnull.h"

#include "libs/globals.h"

#include <functional>
#include <unordered_map>


// sqw info function: "takin_sqw_info"
// returns: [takin ver, ident, long name]
using t_pfkt_info = std::tuple<std::string, std::string, std::string, std::string>(*)();
using t_fkt_info = typename std::remove_pointer<t_pfkt_info>::type;

// sqw module creation function: "takin_sqw"
// old interface returning a shared pointer, which might be dangerous for so files,
// see: https://www.boost.org/doc/libs/1_72_0/doc/html/boost_dll/missuses.html
using t_pfkt = std::shared_ptr<SqwBase>(*)(const std::string&);
using t_fkt = typename std::remove_pointer<t_pfkt>::type;

// new raw pointer constructor interface
using t_pfkt_raw_new = SqwBase*(*)(const std::string&);
using t_fkt_raw_new = typename std::remove_pointer<t_pfkt_raw_new>::type;
using t_pfkt_raw_del = void(*)(SqwBase*);
using t_fkt_raw_del = typename std::remove_pointer<t_pfkt_raw_del>::type;


// key: identifier, value: [func, long name, help text]
using t_mapSqw = std::unordered_map<std::string, std::tuple<t_pfkt, std::string, std::string>>;
using t_mapSqwRaw = std::unordered_map<std::string, std::tuple<t_pfkt_raw_new, t_pfkt_raw_del, std::string, std::string>>;

// key: identifier, value: [long name, binary file name, help text]
using t_mapSqwExt = std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>>;


#endif
