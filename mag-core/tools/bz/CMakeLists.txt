#
# @author Tobias Weber <tweber@ill.fr>
# @date may-2022
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# mag-core (part of the Takin software suite)
# Copyright (C) 2018-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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
project(bz)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-std=c++20)
add_compile_options(-Wall -Wextra)

list(APPEND CMAKE_MODULE_PATH
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/../../../setup/cmake"
)

set(CMAKE_VERBOSE_MAKEFILE TRUE)


# options
option(USE_QT "use Qt" TRUE)
option(MAKE_LIB "create a library" FALSE)
option(USE_BOOST_PROGOPTS "use Boost.program_options" TRUE)


# system specific settings
message("Building for ${CMAKE_SYSTEM_NAME} systems.")


if(USE_BOOST_PROGOPTS)
	find_package(Boost REQUIRED COMPONENTS program_options)
else()
	find_package(Boost)

	message("Compiling without Boost.program_options.")
	add_definitions(-DDONT_USE_BOOTS_PROGOPTS)
endif()
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

if(USE_QT)
	set(QT_VER 5)

	message("Trying to select Qt version ${QT_VER}.")
	if(QT_VER EQUAL 6)
		find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGLWidgets Svg)
		include_directories("${Qt6Core_INCLUDE_DIRS}/..")
		list(APPEND QtLibraries Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGLWidgets Qt6::Svg)
	elseif(QT_VER EQUAL 5)
		find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets Svg)
		include_directories("${Qt5Core_INCLUDE_DIRS}/..")
		list(APPEND QtLibraries Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg)
	else()
		message(FATAL_ERROR "Tried to select non-supported Qt version ${QT_VER}.")
	endif()

	set(CMAKE_AUTOUIC TRUE)
	set(CMAKE_AUTOMOC TRUE)


	add_executable(takin_bz
		bz.cpp bz_file.cpp bz_calc.cpp
		bz_main.cpp
		bz_ops.cpp bz.h
		bz_conf.cpp bz_conf.h
		globals.cpp globals.h
		bz_formulas.cpp
		plot.cpp plot.h
		plot_cut.cpp plot_cut.h
		bz_lib.h
		../../tlibs2/libs/qt/recent.cpp ../../tlibs2/libs/qt/recent.h
		../../tlibs2/libs/qt/gl.cpp ../../tlibs2/libs/qt/gl.h
		../../tlibs2/libs/qt/glplot.cpp ../../tlibs2/libs/qt/glplot.h
		../../libs/symops.h
	)

	target_link_libraries(takin_bz
		${Boost_LIBRARIES}
		${Qhull_LIBRARIES}
		${QtLibraries}
		gemmi::gemmi_cpp
	)
else()
	message("Compiling without Qt.")
	add_definitions(-DDONT_USE_QT)

	add_executable(takin_bz
		bz_main.cpp bz_lib.h
		bz_conf.cpp bz_conf.h
		globals.cpp globals.h

		../../libs/symops.h
	)

	target_link_libraries(takin_bz
		${Boost_LIBRARIES}
		${Qhull_LIBRARIES}
		gemmi::gemmi_cpp
	)
endif()


if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_link_options(takin_bz
		# create an __info_plist section in the binary
		PRIVATE LINKER:-sectcreate,__TEXT,__info_plist,${PROJECT_SOURCE_DIR}/../../../setup/build_mac/plists/InfoBZ.plist
	)

	add_custom_command(TARGET takin_bz POST_BUILD
		COMMAND install_name_tool -add_rpath /usr/local/lib $<TARGET_FILE:takin_bz>
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)
endif()


if(MAKE_LIB)
	message("Compiling as library.")

	add_library(takinbz
		bz_lib.cpp bz_lib.h
		bz_conf.cpp bz_conf.h
		globals.cpp globals.h

		../../libs/symops.h
	)

	target_compile_definitions(takinbz PUBLIC -DDONT_USE_QT)

	target_link_libraries(takinbz
		${Boost_LIBRARIES}
		${Qhull_LIBRARIES}
		gemmi::gemmi_cpp
	)
endif()


if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	# pass linker --subsystem option
	target_compile_options(takin_bz PUBLIC -Wl,--subsystem,windows)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_compile_options(takin_bz PUBLIC -mmacosx-version-min=10.15)
endif()
