#
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

cmake_minimum_required(VERSION 3.5)
project(bz_py)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra)

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/../../../setup/cmake"
)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

# system specific settings
message("Building for ${CMAKE_SYSTEM_NAME} systems.")

find_package(Boost)
add_definitions(-DDONT_USE_BOOTS_PROGOPTS)
add_compile_options(${Boost_CXX_FLAGS})

find_package(Gemmi REQUIRED)

include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.." "../.."
	gemmi::headers
)

find_package(Qhull REQUIRED)
if(Qhull_FOUND)
	add_definitions(-DUSE_QHULL)
	include_directories("${Qhull_INCLUDE_DIRS}")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

find_package(Python3 COMPONENTS Interpreter Development)
find_package(SWIG COMPONENTS python)
#find_package(SWIG COMPONENTS javascript)

if(SWIG_FOUND AND SWIG_python_FOUND)
	message("Scripting using python version ${Python3_VERSION} enabled; packages: ${Python3_SITEARCH}.")

	cmake_policy(SET CMP0078 NEW)
	cmake_policy(SET CMP0086 NEW)

	set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
	include(${SWIG_USE_FILE})

	set_source_files_properties(bzlib.i PROPERTIES CPLUSPLUS TRUE)
	set_source_files_properties(bzlib.i PROPERTIES SWIG_FLAGS "-I ${PROJECT_SOURCE_DIR}")
	#set_source_files_properties(bzlib.i PROPERTIES SWIG_FLAGS -node)
	set_source_files_properties(
		${PROJECT_BINARY_DIR}/CMakeFiles/bz_py.dir/bzlibPYTHON_wrap.cxx
		#${PROJECT_BINARY_DIR}/CMakeFiles/bz_py.dir/bzlibJAVASCRIPT_wrap.cxx
		PROPERTIES SKIP_AUTOMOC TRUE SKIP_AUTOUIC TRUE)

	swig_add_library(bz_py LANGUAGE python SOURCES bzlib.i)
	#swig_add_library(bz_py LANGUAGE javascript SOURCES bzlib.i)

	target_link_libraries(bz_py
		Python3::Python
		${Boost_LIBRARIES}
		${Qhull_LIBRARIES}
		gemmi::gemmi_cpp
	)

	# installation
	get_property(BZ_AUX_FILES TARGET bz_py PROPERTY SWIG_SUPPORT_FILES)
	install(FILES ${BZ_AUX_FILES} DESTINATION "${Python3_SITEARCH}")
	install(TARGETS bz_py LIBRARY DESTINATION "${Python3_SITEARCH}")
endif()
