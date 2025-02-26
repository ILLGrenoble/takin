#
# magdyn py bindings
# @author Tobias Weber <tweber@ill.fr>
# @date 2023-2024
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
project(magdyn_py)

cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/cmake"
	"${PROJECT_SOURCE_DIR}/../../../setup/cmake")


message("Project: ${PROJECT_NAME}, build type: ${CMAKE_BUILD_TYPE}.")

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif()


set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra)



# -----------------------------------------------------------------------------
# Boost
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
find_package(Boost REQUIRED COMPONENTS system filesystem iostreams REQUIRED)
add_definitions(${Boost_CXX_FLAGS})

message("Using Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}.")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Lapack(e)
find_package(GFortran)
find_package(Lapacke REQUIRED)
add_definitions(-DUSE_LAPACK=1)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Gemmi
find_package(Gemmi REQUIRED)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Python
find_package(Python3 REQUIRED COMPONENTS Interpreter Development REQUIRED)

message("Using Python version ${Python3_VERSION}. Dist packages: ${Python3_SITEARCH}")
# -----------------------------------------------------------------------------



include_directories(
	"${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/libs"
	"${PROJECT_SOURCE_DIR}/../.." "${PROJECT_SOURCE_DIR}/../../tlibs2"
	"${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"${Lapacke_INCLUDE_DIRS}"
	"${Python3_INCLUDE_DIRS}"
	gemmi::headers
)



# -----------------------------------------------------------------------------
# Python module
find_package(SWIG REQUIRED)
set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
include(${SWIG_USE_FILE})

set_source_files_properties(magdyn.i PROPERTIES CPLUSPLUS TRUE)
set_source_files_properties(magdyn.i PROPERTIES SWIG_FLAGS "-I ${PROJECT_SOURCE_DIR}/libs")

swig_add_library(magdyn_py LANGUAGE python SOURCES magdyn.i)
target_link_libraries(magdyn_py
	Python3::Python
	${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES} ${GFortran_LIBRARIES}
	gemmi::gemmi_cpp
)
#target_link_libraries(magdyn_py Python3::Python /usr/local/lib/liblapacke.a /usr/local/lib/liblapack.a /usr/local/lib/libblas.a libgfortran.so)

# installation
get_property(MAGDYN_AUX_FILES TARGET magdyn_py PROPERTY SWIG_SUPPORT_FILES)
install(FILES ${MAGDYN_AUX_FILES} DESTINATION "${Python3_SITEARCH}")
install(TARGETS magdyn_py LIBRARY DESTINATION "${Python3_SITEARCH}")
# -----------------------------------------------------------------------------
