#
# @author Tobias Weber <tweber@ill.fr>
# @date 17-jun-2019
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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
project(mag-core)

option(BUILD_PY_MODULES "build py modules" TRUE)
option(ONLY_BUILD_FINISHED "only build finished tools" TRUE)


# mingw settings
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	set(Boost_USE_STATIC_LIBS FALSE)
endif()


# cif converter
message("\n================================================================================")
message("Configuring cif2xml...")
message("================================================================================")
add_subdirectory(tools/cif2xml)

# nuclear structure factor tool
message("\n================================================================================")
message("Configuring structfact...")
message("================================================================================")
add_subdirectory(tools/structfact)

# magnetic structure factor tool
message("\n================================================================================")
message("Configuring magstructfact...")
message("================================================================================")
add_subdirectory(tools/magstructfact)

# Brillouin zone tool
message("\n================================================================================")
message("Configuring bz...")
message("================================================================================")
add_subdirectory(tools/bz)

# magnetic dynamics tool
message("\n================================================================================")
message("Configuring magdyn...")
message("================================================================================")
add_subdirectory(tools/magdyn)

# molecular dynamics tool
message("\n================================================================================")
message("Configuring moldyn...")
message("================================================================================")
add_subdirectory(tools/moldyn)


if(NOT ONLY_BUILD_FINISHED)
	# scan browser tool
	message("\n================================================================================")
	message("Configuring scanbrowser...")
	message("================================================================================")
	add_subdirectory(tools/scanbrowser)

	# setup tools
	message("\n================================================================================")
	message("Configuring setup...")
	message("================================================================================")
	add_subdirectory(tools/setup)

	# magnetic structure factor browser
	message("\n================================================================================")
	message("Configuring magsgbrowser...")
	message("================================================================================")
	add_subdirectory(tools/magsgbrowser)
endif()


if(BUILD_PY_MODULES)
	# magnetic dynamics
	message("\n================================================================================")
	message("Configuring py-magdyn...")
	message("================================================================================")
	add_subdirectory(tools_py/magdyn)

	# brillouin zones
	message("\n================================================================================")
	message("Configuring py-bz...")
	message("================================================================================")
	add_subdirectory(tools_py/bz)

	# instrument data file loader
	message("\n================================================================================")
	message("Configuring py-instr...")
	message("================================================================================")
	add_subdirectory(tools_py/instr)
endif()
