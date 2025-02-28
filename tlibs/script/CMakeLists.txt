#
# Hermelin script
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2013-2017
# @license GPLv2 or GPLv3
#
# ----------------------------------------------------------------------------
# tlibs -- a physical-mathematical C++ template library
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) version 3.
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
project(hermelin)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_BUILD_TYPE Release)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/.." "." ".."
	"${PROJECT_SOURCE_DIR}/../../setup/cmake")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")


add_custom_command(
	OUTPUT ${PROJECT_SOURCE_DIR}/script.cc ${PROJECT_SOURCE_DIR}/tokens.h
	COMMAND bison --defines=${PROJECT_SOURCE_DIR}/tokens.h --output-file=${PROJECT_SOURCE_DIR}/script.cc ${PROJECT_SOURCE_DIR}/lang/script.yy
)


find_package(Boost REQUIRED COMPONENTS iostreams system filesystem REQUIRED)
find_package(Threads REQUIRED)
find_package(Minuit2 REQUIRED)


add_definitions(-std=c++14)
add_definitions(-march=native -s)

add_definitions(-DNO_LAPACK -DNO_QT)
add_definitions(${Boost_CXX_FLAGS})
add_definitions(-DTLIBS_INC_HDR_IMPLS)

add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")



include_directories("${PROJECT_SOURCE_DIR}" "." "..")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${Minuit2_INCLUDE_DIRS}")

#include_directories("/usr/include/root")
#link_directories("/usr/lib64/root")


file(GLOB lang "lang/*.cpp" "lang/*.h")
file(GLOB runtime "runtime/*.cpp" "runtime/*.h")
file(GLOB tlibs
	"../math/linalg2.cpp" "../math/rand.cpp"
	"../string/spec_char.cpp"
	"../gfx/gnuplot.cpp"
	"../file/file.cpp"
	"../log/log.cpp" "../helper/flags.cpp"
)


add_executable(hermelin ${lang} ${runtime} ${tlibs}
	${PROJECT_SOURCE_DIR}/script.cc ${PROJECT_SOURCE_DIR}/tokens.h
)

target_link_libraries(hermelin
	${Boost_LIBRARIES} Threads::Threads ${Minuit2_LIBRARIES}
)


add_custom_command(TARGET hermelin POST_BUILD
	COMMAND strip hermelin
	MAIN_DEPENDENCY hermelin
)



install(TARGETS hermelin DESTINATION bin)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/apps DESTINATION share/hermelin)
