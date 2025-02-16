#
# random distributions
# @author Tobias Weber <tobias.weber@tum.de>
# @date sep-2016
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

project(distr)
cmake_minimum_required(VERSION 3.0)


set(CMAKE_VERBOSE_MAKEFILE TRUE)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "."
	"${PROJECT_SOURCE_DIR}/../../../setup/cmake")


find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qwt REQUIRED)

add_definitions(-I${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(${Boost_CXX_FLAGS})
add_definitions(-std=c++14)


include_directories("${PROJECT_SOURCE_DIR}" "."
	"${PROJECT_SOURCE_DIR}/.." "${PROJECT_SOURCE_DIR}/../..")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${QWT_INCLUDE_DIRS}")


add_executable(distr
	../../log/log.cpp
	distr_main.cpp distr.cpp distr.h)

target_link_libraries(distr
	${Boost_LIBRARIES} Threads::Threads
	Qt5::Core Qt5::Widgets Qt5::Gui
	${QWT_LIBRARIES})
