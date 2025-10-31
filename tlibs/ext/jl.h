/**
 * julia interface helpers
 *
 * @author Tobias Weber <tweber@ill.fr>
 * @date 23-apr-2017, 17-feb-2020
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

#ifndef __TL_JL_H__
#define __TL_JL_H__

#include <julia.h>
#include <tuple>
#include <limits>
#include <string>
#include <functional>

#include <boost/system/error_code.hpp>
#include <boost/dll/shared_library.hpp>


namespace tl
{

// ----------------------------------------------------------------------------
/**
 * Julia data type traits
 */
template<typename T> struct jl_traits {};

template<> struct jl_traits<int64_t>
{
	using value_type = int64_t;

	static jl_datatype_t* get_type() { return jl_int64_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int64(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int64(val); }
};

template<> struct jl_traits<uint64_t>
{
	using value_type = uint64_t;

	static jl_datatype_t* get_type() { return jl_uint64_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint64(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint64(val); }
};

template<> struct jl_traits<int32_t>
{
	using value_type = int32_t;

	static jl_datatype_t* get_type() { return jl_int32_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int32(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int32(val); }
};

template<> struct jl_traits<uint32_t>
{
	using value_type = uint32_t;

	static jl_datatype_t* get_type() { return jl_uint32_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint32(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint32(val); }
};

template<> struct jl_traits<int16_t>
{
	using value_type = int16_t;

	static jl_datatype_t* get_type() { return jl_int16_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int16(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int16(val); }
};

template<> struct jl_traits<uint16_t>
{
	using value_type = uint16_t;

	static jl_datatype_t* get_type() { return jl_uint16_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint16(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint16(val); }
};

template<> struct jl_traits<int8_t>
{
	using value_type = int8_t;

	static jl_datatype_t* get_type() { return jl_int8_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int8(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int8(val); }
};

template<> struct jl_traits<uint8_t>
{
	using value_type = uint8_t;

	static jl_datatype_t* get_type() { return jl_uint8_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint8(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint8(val); }
};

template<> struct jl_traits<float>
{
	using value_type = float;

	static jl_datatype_t* get_type() { return jl_float32_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_float32(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_float32(val); }
};

template<> struct jl_traits<double>
{
	using value_type = double;

	static jl_datatype_t* get_type() { return jl_float64_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_float64(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_float64(val); }
};

template<> struct jl_traits<std::string>
{
	using value_type = std::string;

	static jl_datatype_t* get_type() { return jl_string_type; }
	static value_type unbox(jl_value_t *pVal)
	{
		const char* pc = jl_string_ptr(pVal);
		return std::string(pc);
	}
	static jl_value_t* box(value_type val) { return jl_cstr_to_string(val.c_str()); }
};
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
/**
 * converts an stl container of containers into a julia 2d array
 */
template<template<class...> class t_cont/*=std::vector*/, class T/*=double*/>
jl_array_t* make_jl_2darr(const t_cont<t_cont<T>>& vecvec)
{
	// number of columns and rows
	std::size_t iCols = vecvec.size();
	std::size_t iRows = std::numeric_limits<std::size_t>::max();
	for(const auto& vec : vecvec)
		iRows = std::min(iRows, vec.size());
	if(!iCols) iRows = 0;

	jl_array_t *pArr = jl_alloc_array_2d(
		jl_apply_array_type(reinterpret_cast<jl_value_t*>(jl_traits<T>::get_type()), 2), 
		iRows, iCols);
	T* pDat = jl_array_data(pArr, T);

	std::size_t iCurCol = 0;
	for(const auto& vec : vecvec)
	{
		std::size_t iCurRow = 0;
		for(const auto& d : vec)
		{
			pDat[iCurCol*iRows + iCurRow] = d;
			++iCurRow;
		}

		++iCurCol;
	}

	return pArr;
}


/**
 * converts an stl container of strings into a julia array of strings
 */
template<template<class...> class t_cont/*=std::vector*/, class t_str/*=std::string*/>
jl_array_t* make_jl_str_arr(const t_cont<t_str>& vecStr)
{
	jl_array_t *pArr = jl_alloc_array_1d(jl_apply_array_type(
		reinterpret_cast<jl_value_t*>(jl_string_type), 1), vecStr.size());
	jl_value_t** pDat = jl_array_data(pArr, jl_value_t*);

	std::size_t iIdx = 0;
	for(const t_str& str : vecStr)
	{
		pDat[iIdx] = jl_cstr_to_string(str.c_str());
		++iIdx;
	}

	return pArr;
}


/**
 * converts a julia array into an stl container
 */
template<template<class...> class t_cont/*=std::vector*/, class t_val>
t_cont<t_val> from_jl_arr(jl_array_t *pArr, std::size_t iSkipFront = 0)
{
	const std::size_t iSize = jl_array_len(pArr);

	t_cont<t_val> vecRet;
	vecRet.reserve(iSize);

	const t_val *arr = jl_array_data(pArr, t_val);

	for(std::size_t iElem = iSkipFront; iElem < iSize; ++iElem)
		vecRet.push_back(arr[iElem]);

	return vecRet;
}


/**
 * converts a map of strings into two julia arrays of strings (key & value)
 */
template<template<class...> class t_cont/*=std::map*/, class t_str/*=std::string*/>
std::tuple<jl_array_t*, jl_array_t*> make_jl_strmap_arr(const t_cont<t_str, t_str>& map)
{
	jl_array_t *pArrKey = jl_alloc_array_1d(
		jl_apply_array_type(reinterpret_cast<jl_value_t*>(jl_string_type), 1), map.size());
	jl_array_t *pArrVal = jl_alloc_array_1d(
		jl_apply_array_type(reinterpret_cast<jl_value_t*>(jl_string_type), 1), map.size());

	jl_value_t** pDatKey = jl_array_data(pArrKey, jl_value_t*);
	jl_value_t** pDatVal = jl_array_data(pArrVal, jl_value_t*);

	std::size_t iIdx = 0;
	for(const auto& pair : map)
	{
		pDatKey[iIdx] = jl_cstr_to_string(pair.first.c_str());
		pDatVal[iIdx] = jl_cstr_to_string(pair.second.c_str());
		++iIdx;
	}

	return std::make_tuple(pArrKey, pArrVal);
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// dynamic loading of libjulia, wraps (a subset of) the Julia native API
// ----------------------------------------------------------------------------
class LibJulia
{
	public:
		LibJulia(const std::string& name = "libjulia") :
			m_lib{name,
			boost::dll::load_mode::search_system_folders | boost::dll::load_mode::append_decorations |
			boost::dll::load_mode::rtld_lazy | boost::dll::load_mode::rtld_global}
		{
			try
			{
				if(!IsLoaded())
					return;

				this->jl_init = m_lib.get<void()>("jl_init__threading");
				if(!this->jl_init)
					this->jl_init = m_lib.get<void()>("jl_init");

				this->jl_ver_string = m_lib.get<decltype(::jl_ver_string)>("jl_ver_string");

				this->jl_get_global = m_lib.get<decltype(::jl_get_global)>("jl_get_global");
				this->jl_symbol = m_lib.get<decltype(::jl_symbol)>("jl_symbol");
				this->jl_typeof_str = m_lib.get<decltype(::jl_typeof_str)>("jl_typeof_str");

				this->jl_call0 = m_lib.get<decltype(::jl_call0)>("jl_call0");
				this->jl_call1 = m_lib.get<decltype(::jl_call1)>("jl_call1");
				this->jl_call2 = m_lib.get<decltype(::jl_call2)>("jl_call2");
				this->jl_call = m_lib.get<decltype(::jl_call)>("jl_call");

				this->jl_eval_string = m_lib.get<decltype(::jl_eval_string)>("jl_eval_string");
				this->jl_string_ptr = m_lib.get<decltype(::jl_string_ptr)>("jl_string_ptr");
				this->jl_cstr_to_string = m_lib.get<decltype(::jl_cstr_to_string)>("jl_cstr_to_string");

				this->jl_get_nth_field = m_lib.get<decltype(::jl_get_nth_field)>("jl_get_nth_field");

				this->jl_box_int8 = m_lib.get<decltype(::jl_box_int8)>("jl_box_int8");
				this->jl_box_int16 = m_lib.get<decltype(::jl_box_int16)>("jl_box_int16");
				this->jl_box_int32 = m_lib.get<decltype(::jl_box_int32)>("jl_box_int32");
				this->jl_box_int64 = m_lib.get<decltype(::jl_box_int64)>("jl_box_int64");
				this->jl_box_uint8 = m_lib.get<decltype(::jl_box_uint8)>("jl_box_uint8");
				this->jl_box_uint16 = m_lib.get<decltype(::jl_box_uint16)>("jl_box_uint16");
				this->jl_box_uint32 = m_lib.get<decltype(::jl_box_uint32)>("jl_box_uint32");
				this->jl_box_uint64 = m_lib.get<decltype(::jl_box_uint64)>("jl_box_uint64");
				this->jl_box_float32 = m_lib.get<decltype(::jl_box_float32)>("jl_box_float32");
				this->jl_box_float64 = m_lib.get<decltype(::jl_box_float64)>("jl_box_float64");
				this->jl_box_voidpointer = m_lib.get<decltype(::jl_box_voidpointer)>("jl_box_voidpointer");
				this->jl_box_bool = m_lib.get<decltype(::jl_box_bool)>("jl_box_bool");
				this->jl_box_char = m_lib.get<decltype(::jl_box_char)>("jl_box_char");

				this->jl_unbox_int8 = m_lib.get<decltype(::jl_unbox_int8)>("jl_unbox_int8");
				this->jl_unbox_int16 = m_lib.get<decltype(::jl_unbox_int16)>("jl_unbox_int16");
				this->jl_unbox_int32 = m_lib.get<decltype(::jl_unbox_int32)>("jl_unbox_int32");
				this->jl_unbox_int64 = m_lib.get<decltype(::jl_unbox_int64)>("jl_unbox_int64");
				this->jl_unbox_uint8 = m_lib.get<decltype(::jl_unbox_uint8)>("jl_unbox_uint8");
				this->jl_unbox_uint16 = m_lib.get<decltype(::jl_unbox_uint16)>("jl_unbox_uint16");
				this->jl_unbox_uint32 = m_lib.get<decltype(::jl_unbox_uint32)>("jl_unbox_uint32");
				this->jl_unbox_uint64 = m_lib.get<decltype(::jl_unbox_uint64)>("jl_unbox_uint64");
				this->jl_unbox_float32 = m_lib.get<decltype(::jl_unbox_float32)>("jl_unbox_float32");
				this->jl_unbox_float64 = m_lib.get<decltype(::jl_unbox_float64)>("jl_unbox_float64");
				this->jl_unbox_voidpointer = m_lib.get<decltype(::jl_unbox_voidpointer)>("jl_unbox_voidpointer");
				this->jl_unbox_bool = m_lib.get<decltype(::jl_unbox_bool)>("jl_unbox_bool");

				m_ok = true;
			}
			catch(const std::exception& ex)
			{
				m_ok = false;
				throw;
			}
		}


		bool IsLoaded() const { return m_lib.is_loaded(); }
		bool IsOk() const { return IsLoaded() && m_ok; }

		std::string GetFileName() const { return m_lib.location().string(); }


	protected:
		bool m_ok = false;
		boost::dll::shared_library m_lib;


	public:
		// wrapped jl functions
		std::function<void()> jl_init = nullptr;
		std::function<decltype(::jl_ver_string)> jl_ver_string = nullptr;

		std::function<decltype(::jl_get_global)> jl_get_global = nullptr;
		std::function<decltype(::jl_symbol)> jl_symbol = nullptr;
		std::function<decltype(::jl_typeof_str)> jl_typeof_str = nullptr;

		std::function<decltype(::jl_call0)> jl_call0 = nullptr;
		std::function<decltype(::jl_call1)> jl_call1 = nullptr;
		std::function<decltype(::jl_call2)> jl_call2 = nullptr;
		std::function<decltype(::jl_call)> jl_call = nullptr;

		std::function<decltype(::jl_eval_string)> jl_eval_string = nullptr;
		std::function<decltype(::jl_string_ptr)> jl_string_ptr = nullptr;
		std::function<decltype(::jl_cstr_to_string)> jl_cstr_to_string = nullptr;

		std::function<decltype(::jl_get_nth_field)> jl_get_nth_field = nullptr;

		std::function<decltype(::jl_box_int8)> jl_box_int8 = nullptr;
		std::function<decltype(::jl_box_int16)> jl_box_int16 = nullptr;
		std::function<decltype(::jl_box_int32)> jl_box_int32 = nullptr;
		std::function<decltype(::jl_box_int64)> jl_box_int64 = nullptr;
		std::function<decltype(::jl_box_uint8)> jl_box_uint8 = nullptr;
		std::function<decltype(::jl_box_uint16)> jl_box_uint16 = nullptr;
		std::function<decltype(::jl_box_uint32)> jl_box_uint32 = nullptr;
		std::function<decltype(::jl_box_uint64)> jl_box_uint64 = nullptr;
		std::function<decltype(::jl_box_float32)> jl_box_float32 = nullptr;
		std::function<decltype(::jl_box_float64)> jl_box_float64 = nullptr;
		std::function<decltype(::jl_box_voidpointer)> jl_box_voidpointer = nullptr;
		std::function<decltype(::jl_box_bool)> jl_box_bool = nullptr;
		std::function<decltype(::jl_box_char)> jl_box_char = nullptr;

		std::function<decltype(::jl_unbox_int8)> jl_unbox_int8 = nullptr;
		std::function<decltype(::jl_unbox_int16)> jl_unbox_int16 = nullptr;
		std::function<decltype(::jl_unbox_int32)> jl_unbox_int32 = nullptr;
		std::function<decltype(::jl_unbox_int64)> jl_unbox_int64 = nullptr;
		std::function<decltype(::jl_unbox_uint8)> jl_unbox_uint8 = nullptr;
		std::function<decltype(::jl_unbox_uint16)> jl_unbox_uint16 = nullptr;
		std::function<decltype(::jl_unbox_uint32)> jl_unbox_uint32 = nullptr;
		std::function<decltype(::jl_unbox_uint64)> jl_unbox_uint64 = nullptr;
		std::function<decltype(::jl_unbox_float32)> jl_unbox_float32 = nullptr;
		std::function<decltype(::jl_unbox_float64)> jl_unbox_float64 = nullptr;
		std::function<decltype(::jl_unbox_voidpointer)> jl_unbox_voidpointer = nullptr;
		std::function<decltype(::jl_unbox_bool)> jl_unbox_bool = nullptr;
};
// ----------------------------------------------------------------------------



}
#endif
