#
# Takin
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2013-2020
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------
#

cmake_minimum_required(VERSION 3.5)
project(takin)

message("System name: ${CMAKE_SYSTEM_NAME}.")


option(DEBUG "debug build" FALSE)
option(USE_PLUGINS "use the plugin system" TRUE)
option(USE_STRIP "strip debug info from release binaries" TRUE)
option(USE_GENTAB "build setup utilities" TRUE)
option(USE_RtMp "use auxilliary libraries" TRUE)
option(USE_Dl "use dl libraries" TRUE)
option(USE_Py "search for python" TRUE)
option(USE_LAPACKE "search for Lapacke" FALSE)
option(USE_CUSTOM_THREADPOOL "use custom threadpool implementation" FALSE)
option(USE_INTERPROC_EMUL "force shared memory emulation mode" FALSE)
option(DISABLE_INTERPROC_XSI "disable xsi shared memory syscalls" FALSE)

set(USE_TR1_FUNCS FALSE)
set(BUILD_EXT_TOOLS FALSE)


if(DEBUG)
	message("Doing a DEBUG build.")

	set(CMAKE_VERBOSE_MAKEFILE TRUE)
	set(CMAKE_BUILD_TYPE Debug)

	set(Boost_DEBUG TRUE)
	#set(Boost_USE_DEBUG_LIBS TRUE)
	#set(Boost_USE_RELEASE_LIBS FALSE)

	add_definitions(-DDEBUG)
	add_compile_options(-ggdb -rdynamic)
	add_compile_options(-Wall -Wpedantic -Wextra -Wno-unused-local-typedefs -Wno-unused-parameter)
else()
	message("Doing a RELEASE build.")

	set(CMAKE_VERBOSE_MAKEFILE FALSE)
	set(CMAKE_BUILD_TYPE Release)

	set(Boost_DEBUG FALSE)
	set(Boost_USE_DEBUG_LIBS FALSE)
	set(Boost_USE_RELEASE_LIBS TRUE)
	#set(Boost_USE_STATIC_LIBS TRUE)

	add_definitions(-DNDEBUG -DQT_NO_DEBUG)
	#add_compile_options(-Wno-narrowing)
	#add_compile_options(-march=native -s)
endif()


set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_CXX_STANDARD 14)
add_compile_options(-std=c++14)
#add_compile_options(-fPIC)
add_compile_options(-fvisibility=hidden)

# for application plugins
add_definitions(-DBUILD_APPLI)


set(SOCK2 "")

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	message("Setting MinGW options.")

	# stuff needed by mingw
	set(SOCK2 "ws2_32")
	set(USE_GENTAB FALSE)
	set(USE_RtMp FALSE)
	set(USE_Dl FALSE)
	set(USE_Py FALSE)

	set(Boost_USE_STATIC_LIBS FALSE)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	message("Setting Darwin options.")

	add_compile_options(-mmacosx-version-min=10.10)
	set(USE_GENTAB FALSE)
	set(USE_RtMp FALSE)
	set(USE_Dl FALSE)
	set(USE_Py TRUE)
endif()


if(USE_CUSTOM_THREADPOOL)
	message("Using custom threadpool implementation.")
	add_definitions(-DUSE_OWN_THREADPOOL)
endif()


if(USE_INTERPROC_EMUL)
	message("Using interprocess shared memory emulation mode.")
	add_definitions(-D__USE_SQW_INTERPROC_EMUL__)
endif()


if(DISABLE_INTERPROC_XSI)
	message("Disabling xsi shared memory syscalls.")
	add_definitions(-D__DISABLE_SQW_INTERPROC_XSI__)
endif()


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "."
	"${PROJECT_SOURCE_DIR}/cmake" "${PROJECT_SOURCE_DIR}/../setup/cmake")
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/plugins")
set(CMAKE_AUTOUIC_SEARCH_PATHS "${PROJECT_SOURCE_DIR}/ui")

set(Boost_FIND_QUIETLY FALSE)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	#set(BOOST_SUFFIX -x64)
	set(BOOST_SUFFIX)
else()
	set(BOOST_SUFFIX)
endif()

set(BOOST_COMPS program_options${BOOST_SUFFIX} iostreams${BOOST_SUFFIX} filesystem${BOOST_SUFFIX} system${BOOST_SUFFIX})
set(BOOST_OPT_COMPS)

if(USE_TR1_FUNCS)
	list(APPEND BOOST_COMPS regex)
endif()

if(USE_Py)
	find_package(Python3 COMPONENTS Interpreter Development)
endif()

if(Python3_Interpreter_FOUND AND Python3_Development_FOUND)
	list(APPEND BOOST_OPT_COMPS python3 python37 python38 python39 python310 python311 python312)
	set(Py_FOUND TRUE)
else()
	set(Py_FOUND FALSE)
endif()


find_package(Boost REQUIRED COMPONENTS ${BOOST_COMPS} OPTIONAL_COMPONENTS ${BOOST_OPT_COMPS})
find_package(Threads REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui Svg) #OpenGL)
find_package(Qwt REQUIRED)
find_package(Minuit2 REQUIRED)
find_package(PNG REQUIRED)
find_package(Julia)
find_package(ZLIB)
find_package(BZip2)
find_package(HDF5 COMPONENTS CXX)

if(Py_FOUND)
	if(TARGET Boost::python3)
		set(_BOOST_PY_TARGET Boost::python3)
	elseif(TARGET Boost::python312)
		set(_BOOST_PY_TARGET Boost::python312)
	elseif(TARGET Boost::python311)
		set(_BOOST_PY_TARGET Boost::python311)
	elseif(TARGET Boost::python310)
		set(_BOOST_PY_TARGET Boost::python310)
	elseif(TARGET Boost::python39)
		set(_BOOST_PY_TARGET Boost::python39)
	elseif(TARGET Boost::python38)
		set(_BOOST_PY_TARGET Boost::python38)
	elseif(TARGET Boost::python37)
		set(_BOOST_PY_TARGET Boost::python37)
	else()
		message("Boost.Python not found.")
		set(Py_FOUND FALSE)
	endif()
endif()

if(USE_RtMp)
	find_package(RtMp REQUIRED)
else()
	set(Rt_LIBRARIES "")
	set(Mp_LIBRARIES "")
endif()

if(USE_PLUGINS AND USE_Dl AND NOT (CMAKE_SYSTEM_NAME STREQUAL "Windows"))
	find_package(Dl REQUIRED)
else()
	set(Dl_LIBRARIES "")
endif()
#list(APPEND Dl_LIBRARIES -ldl_test)

#if(USE_GENTAB)
#	find_package(Clipper REQUIRED)
#endif()

set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL)
find_package(Freetype)


if("${QT_MOC_EXECUTABLE}" STREQUAL "")
	get_target_property(QT_MOC_EXECUTABLE Qt5::moc IMPORTED_LOCATION)
endif()

if("${QT_UIC_EXECUTABLE}" STREQUAL "")
	get_target_property(QT_UIC_EXECUTABLE Qt5::uic IMPORTED_LOCATION)
endif()

#if("${QT_QCOLLECTIONGENERATOR_EXECUTABLE}" STREQUAL "")
#	get_target_property(QT_QCOLLECTIONGENERATOR_EXECUTABLE Qt5::qcollectiongenerator IMPORTED_LOCATION)
#endif()


message("Project source dir: ${PROJECT_SOURCE_DIR}")
message("CMake source dir: ${CMAKE_SOURCE_DIR}")
message("CMake binary dir: ${CMAKE_BINARY_DIR}")
message("Install path prefix: ${CMAKE_INSTALL_PREFIX}")
message("Using Qt version ${Qt5Core_VERSION}.")
if(NOT "${QWT_MAJOR_VERSION}" STREQUAL "")
	message("Using Qwt version ${QWT_MAJOR_VERSION}.")
endif()
message("Using Boost version ${Boost_VERSION}.")
message("Using Zlib version ${ZLIB_VERSION_STRING}.")
message("Using Bzip2 version ${BZIP2_VERSION_STRING}.")
if(NOT "${Python3_VERSION}" STREQUAL "") #AND Py_FOUND)
	message("Using Py version ${Python3_VERSION}.")
endif()
if(NOT "${Boost_PYTHON_VERSION}" STREQUAL "")
	message("Using Boost.Py version ${Boost_PYTHON_VERSION}.")
endif()
message("Qt moc: ${QT_MOC_EXECUTABLE}\nQt uic: ${QT_UIC_EXECUTABLE}")


add_definitions(-DNO_QHULL)
add_definitions(-DNO_JPEG -DNO_TIFF)
add_definitions(-DINSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
#add_definitions(${Boost_CXX_FLAGS} -DBOOST_FILESYSTEM_NO_DEPRECATED -DBOOST_ERROR_CODE_HEADER_ONLY)
add_definitions(-DTLIBS_INC_HDR_IMPLS)
add_definitions(-DHAS_COMPLEX_ERF)

if(USE_PLUGINS)
	add_definitions(-DUSE_PLUGINS)
endif()

if(USE_TR1_FUNCS)
	add_definitions(-DUSE_BOOST_REX)
endif()

if(USE_LAPACKE)
	find_package(GFortran)
	find_package(Lapacke)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
	add_definitions(-DNO_HELP_ASSISTANT)
	add_definitions(-DNO_TERM_CMDS)
	add_definitions(-DUSE_WIDE_STR)
	#add_definitions(-DNO_IOSTR)

	# pass linker --subsystem option
	add_compile_options(-Wl,--subsystem,windows)
endif()

# hack to prevent math_fwd.h inclusion
add_definitions(-DBOOST_MATH_FWD_HPP)


#add_definitions("${Qt5Core_DEFINITIONS}")
#add_definitions("${Qt5Widgets_DEFINITIONS}")
#add_definitions("${QtGui_DEFINITIONS}")
#add_definitions("${Qt5OpenGL_DEFINITIONS}")

include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/3rdparty")
include_directories("." "./3rdparty")
include_directories("${Boost_INCLUDE_DIRS}/..")
include_directories("${PNG_INCLUDE_DIRS}")
#include_directories("${JPEG_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS}")
include_directories("${Qt5Core_INCLUDE_DIRS}" "${Qt5Gui_INCLUDE_DIRS}" "${Qt5Widgets_INCLUDE_DIRS}")
include_directories("${QWT_INCLUDE_DIRS}")
include_directories("${ZLIB_INCLUDE_DIRS}")
include_directories("${BZIP2_INCLUDE_DIR}")
include_directories("${Minuit2_INCLUDE_DIRS}")
if(USE_GENTAB)
	include_directories("./3rdparty/clipper-master/clipper")
endif()

set(LIBS_FIT ${Minuit2_LIBRARIES})


if(HDF5_FOUND)
	message("Using HDF5 version ${HDF5_VERSION}.")
	include_directories("${HDF5_CXX_INCLUDE_DIRS}")
	add_definitions(-DUSE_HDF5)
else()
	message("Disabling HDF5 file support.")
	set(HDF5_CXX_LIBRARIES "")
endif()



# py scripting support
if(Py_FOUND)
	message("Enabling py scripting support.")
	add_definitions(-DUSE_PY)
	include_directories("${Python3_INCLUDE_DIRS}")
else()
	message("Disabling py scripting support.")
	add_definitions(-DNO_PY)
endif()



# jl scripting support
if(Julia_FOUND)
	message("Enabling jl scripting support.")
	add_definitions(-DUSE_JL)
	include_directories("${Julia_INCLUDE_DIRS}")
else()
	message("Disabling jl scripting support.")
	add_definitions(-DNO_JL)
endif()



# lapacke support
if(Lapacke_FOUND)
	message("Using Lapacke.")
	add_definitions(-DUSE_LAPACK=1)
	set(LIBS_LINALG ${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES} ${GFortran_LIBRARIES})
	include_directories("${Lapacke_INCLUDE_DIRS}")
else()
	message("Disabling Lapacke support.")
	set(LIBS_LINALG "")
	add_definitions(-DNO_LAPACK)
endif()



# 3d support
if(OPENGL_FOUND AND FREETYPE_FOUND)
	message("Enabling 3D support. Inc: \"${OPENGL_INCLUDE_DIR}\".")

	include_directories("${FREETYPE_INCLUDE_DIRS}")
	if(NOT "${OPENGL_INCLUDE_DIR}" STREQUAL "")
		include_directories("${OPENGL_INCLUDE_DIR}")
		include_directories("${OPENGL_INCLUDE_DIR}/GL" "${OPENGL_INCLUDE_DIR}/OpenGL")
	else()
		message("Empty OGL dirs, guessing some defaults.")
		include_directories("/usr/include/GL" "/usr/include/OpenGL" "/usr/x86_64-w64-mingw32/sys-root/mingw/include/GL")
	endif()

	set(SRCS_3D
		tlibs/gfx/gl.cpp
		libs/plotgl.cpp tools/taz/recip3d.cpp
		tools/taz/real3d.cpp tools/taz/bz3d.cpp dialogs/EllipseDlg3D.cpp
	)

	set(LIBS_3D
		${FREETYPE_LIBRARY}
		${OPENGL_LIBRARIES}  # Qt5::OpenGL
	)
else()
	message("Disabling 3D support: OGL: ${OPENGL_FOUND}, FT: ${FREETYPE_FOUND}")
	add_definitions(-DNO_3D)
endif()



set_property(SOURCE 3rdparty/Faddeeva.hh PROPERTY SKIP_AUTOGEN TRUE)
set(CMAKE_AUTOMOC TRUE)
set(CMAKE_AUTOUIC TRUE)
set(CMAKE_AUTORCC FALSE)



# -----------------------------------------------------------------------------
# gentab
# -----------------------------------------------------------------------------
if(USE_GENTAB)
	add_executable(gentab
		tools/gentab/gentab.cpp
		#libs/spacegroups/spacegroup_clp.cpp libs/spacegroups/crystalsys.cpp
		tlibs/log/log.cpp

		# clipper library
		3rdparty/clipper-master/clipper/core/atomsf.cpp
		3rdparty/clipper-master/clipper/core/coords.cpp
		3rdparty/clipper-master/clipper/core/cell.cpp
		3rdparty/clipper-master/clipper/core/rotation.cpp
		3rdparty/clipper-master/clipper/core/symop.cpp
		3rdparty/clipper-master/clipper/core/clipper_types.cpp
		3rdparty/clipper-master/clipper/core/clipper_util.cpp
		3rdparty/clipper-master/clipper/core/clipper_message.cpp
	)

	target_link_libraries(gentab
		Qt5::Core Qt5::Gui Qt5::Widgets
		Boost::iostreams${BOOST_SUFFIX} Boost::system${BOOST_SUFFIX}
		Boost::filesystem${BOOST_SUFFIX} Boost::program_options${BOOST_SUFFIX}
		${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
		#${Clipper_LIBRARIES}
	)

	add_custom_command(TARGET gentab POST_BUILD
		COMMAND cd ${CMAKE_BINARY_DIR} && ./bin/gentab
		MAIN_DEPENDENCY gentab)
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# takin
# -----------------------------------------------------------------------------
add_executable(takin
	# main
	tools/taz/taz_main.cpp
	tools/taz/taz.cpp tools/taz/tas_layout.cpp tools/taz/tof_layout.cpp
	tools/taz/taz_crys.cpp
	tools/taz/taz_file.cpp tools/taz/taz_export.cpp
	tools/taz/tas_layout.cpp tools/taz/scattering_triangle.cpp
	tools/taz/real_lattice.cpp tools/taz/proj_lattice.cpp
	tools/taz/nicos.cpp tools/taz/sics.cpp tools/taz/taz_net.cpp

	dialogs/SettingsDlg.cpp dialogs/FilePreviewDlg.cpp
	dialogs/GotoDlg.cpp dialogs/ElasticDlg.cpp dialogs/GenPosDlg.cpp
	dialogs/ScatteringFactorsDlg.cpp dialogs/DynPlaneDlg.cpp
	dialogs/NeutronDlg.cpp dialogs/TOFDlg.cpp
	dialogs/SpurionDlg.cpp dialogs/PowderDlg.cpp
	dialogs/FavDlg.cpp
	dialogs/RecipParamDlg.cpp dialogs/RealParamDlg.cpp
	dialogs/SrvDlg.cpp dialogs/NetCacheDlg.cpp dialogs/ScanMonDlg.cpp
	dialogs/EllipseDlg.cpp tools/ffact/FormfactorDlg.cpp
	dialogs/AtomsDlg.cpp dialogs/AboutDlg.cpp
	dialogs/DarkAnglesDlg.cpp dialogs/LogDlg.cpp

	# reso
	tools/res/cn.cpp tools/res/pop.cpp tools/res/pop_cn.cpp
	tools/res/eck.cpp tools/res/vio.cpp tools/res/simple.cpp
	tools/res/ResoDlg.cpp tools/res/ResoDlg_file.cpp
	tools/res/res_cli.cpp

	# monteconvo
	tools/monteconvo/ConvoDlg.cpp tools/monteconvo/ConvoDlg_file.cpp
	tools/monteconvo/ConvoDlg_sim.cpp tools/monteconvo/ConvoDlg_fit.cpp
	tools/monteconvo/SqwParamDlg.cpp tools/monteconvo/TASReso.cpp
	tools/monteconvo/modules/elast.cpp
	tools/monteconvo/modules/kdtree.cpp
	tools/monteconvo/modules/simple_magnon.cpp
	tools/monteconvo/modules/simple_phonon.cpp
	tools/monteconvo/modules/table1d.cpp
	tools/monteconvo/modules/uniform_grid.cpp
	tools/monteconvo/sqwbase.cpp tools/monteconvo/sqwfactory.cpp
	tools/monteconvo/monteconvo_cli.cpp tools/monteconvo/monteconvo_common.cpp

	# convofit
	tools/convofit/convofit.cpp tools/convofit/convofit_import.cpp
	tools/convofit/model.cpp tools/convofit/scan.cpp
	tools/convofit/convofit_cli.cpp

	# scanviewer
	tools/scanviewer/scanviewer.cpp tools/scanviewer/scanviewer_fit.cpp
	tools/scanviewer/scanviewer_pol.cpp
	tools/scanviewer/FitParamDlg.cpp
	tools/scanpos/ScanPosDlg.cpp
	tools/powderfit/PowderFitDlg.cpp

	# libs
	libs/spacegroups/spacegroup.cpp tools/sglist/SgListDlg.cpp
	libs/globals.cpp libs/globals_qt.cpp libs/spacegroups/crystalsys.cpp
	libs/formfactors/formfact.cpp libs/qt/qthelper.cpp libs/qt/qwthelper.cpp

	${SRCS_3D}

	tlibs/file/x3d.cpp
	tlibs/log/log.cpp tlibs/log/debug.cpp
	tlibs/math/rand.cpp
	tlibs/string/spec_char.cpp
	tlibs/helper/flags.cpp  tlibs/version.cpp

	# external libs
	./3rdparty/Faddeeva.cc
)

target_link_libraries(takin
	${LIBS_FIT} ${LIBS_LINALG}
	${SOCK2}
	Qt5::Core Qt5::Gui Qt5::Svg
	${QWT_LIBRARIES}
	${LIBS_3D}
	${PNG_LIBRARIES}
	${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
	Boost::iostreams${BOOST_SUFFIX} Boost::system${BOOST_SUFFIX}
	Boost::filesystem${BOOST_SUFFIX} Boost::program_options${BOOST_SUFFIX}
	Threads::Threads ${Mp_LIBRARIES} ${Rt_LIBRARIES} ${Dl_LIBRARIES}
	${HDF5_CXX_LIBRARIES}
)


if(CMAKE_BUILD_TYPE STREQUAL "Release" AND USE_STRIP)
	add_custom_command(TARGET takin POST_BUILD
		COMMAND strip -v $<TARGET_FILE:takin>
		MAIN_DEPENDENCY takin
	)
endif()


if(USE_GENTAB)
	add_dependencies(takin gentab)
endif()
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# build plugins
# -----------------------------------------------------------------------------

if(Py_FOUND)
#add_library(takinmod_py SHARED
add_executable(takinmod_py
	tools/monteconvo/script_modules/sqw_py.cpp
	tools/monteconvo/sqwbase.cpp
	tlibs/log/log.cpp
	tlibs/math/rand.cpp
)


if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_link_options(takinmod_py
		# create an __info_plist section in the binary
		PRIVATE LINKER:-sectcreate,__TEXT,__info_plist,${PROJECT_SOURCE_DIR}/../setup/build_mac/plists/InfoModPy.plist
	)
endif()

target_link_libraries(takinmod_py
	Boost::iostreams${BOOST_SUFFIX} Boost::system Boost::filesystem${BOOST_SUFFIX}
	${_BOOST_PY_TARGET} ${Python3_LIBRARIES}
	Threads::Threads ${Rt_LIBRARIES} ${Dl_LIBRARIES}
)


if(CMAKE_BUILD_TYPE STREQUAL "Release" AND USE_STRIP)
	add_custom_command(TARGET takinmod_py POST_BUILD
		COMMAND strip -v $<TARGET_FILE:takinmod_py>
		MAIN_DEPENDENCY takinmod_py
	)
endif()
endif()



if(Julia_FOUND)
#add_library(takinmod_jl SHARED
add_executable(takinmod_jl
	tools/monteconvo/script_modules/sqw_jl.cpp
	tools/monteconvo/sqwbase.cpp
	tlibs/log/log.cpp
	tlibs/math/rand.cpp
)

target_link_libraries(takinmod_jl
	Boost::iostreams${BOOST_SUFFIX} Boost::system${BOOST_SUFFIX}
	Boost::filesystem${BOOST_SUFFIX}
	${Julia_LIBRARIES}
	Threads::Threads ${Rt_LIBRARIES}
)


if(CMAKE_BUILD_TYPE STREQUAL "Release" AND USE_STRIP)
	add_custom_command(TARGET takinmod_jl POST_BUILD
		COMMAND strip -v $<TARGET_FILE:takinmod_jl>
		MAIN_DEPENDENCY takinmod_jl
	)
endif()
endif()
# -----------------------------------------------------------------------------



if(BUILD_EXT_TOOLS)
	add_definitions(-DCONVOFIT_STANDALONE)

	# -----------------------------------------------------------------------------
	# convofit, convoseries
	# -----------------------------------------------------------------------------
	add_executable(takin_convofit
		tools/res/cn.cpp tools/res/pop.cpp tools/res/pop_cn.cpp
		tools/res/eck.cpp tools/res/vio.cpp
		#tools/res/simple.cpp

		tools/monteconvo/TASReso.cpp
		tools/monteconvo/modules/kdtree.cpp
		tools/monteconvo/modules/simple_magnon.cpp
		tools/monteconvo/modules/simple_phonon.cpp
		tools/monteconvo/modules/table1d.cpp
		tools/monteconvo/modules/uniform_grid.cpp
		tools/monteconvo/sqwbase.cpp tools/monteconvo/sqwfactory.cpp

		tools/convofit/convofit.cpp tools/convofit/convofit_import.cpp
		tools/convofit/model.cpp tools/convofit/scan.cpp
		tools/convofit/convofit_cli.cpp tools/convofit/convofit_cli_main.cpp

		# statically link tlibs externals
		tlibs/log/log.cpp
		tlibs/math/rand.cpp
		tlibs/file/tmp.cpp
		#tlibs/file/loadinstr.cpp
		libs/globals.cpp
	)

	set_target_properties(takin_convofit PROPERTIES COMPILE_FLAGS "-DNO_QT")

	target_link_libraries(takin_convofit
		${LIBS_FIT} ${LIBS_LINALG}
		${SOCK2}
		Threads::Threads ${Mp_LIBRARIES} ${Rt_LIBRARIES} ${Dl_LIBRARIES}
		Boost::iostreams${BOOST_SUFFIX} Boost::system${BOOST_SUFFIX}
		Boost::filesystem${BOOST_SUFFIX} Boost::program_options${BOOST_SUFFIX}
		${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
	)

	if(CMAKE_BUILD_TYPE STREQUAL "Release" AND USE_STRIP)
		add_custom_command(TARGET takin_convofit POST_BUILD
			COMMAND strip -v $<TARGET_FILE:takin_convofit>
			MAIN_DEPENDENCY takin_convofit
		)
	endif()

	# -----------------------------------------------------------------------------

	add_executable(takin_convoseries
		tools/convofit/scanseries.cpp

		# statically link tlibs externals
		tlibs/log/log.cpp
	)

	set_target_properties(takin_convoseries PROPERTIES COMPILE_FLAGS "-DNO_QT")

	target_link_libraries(takin_convoseries
		Boost::iostreams${BOOST_SUFFIX} Boost::system${BOOST_SUFFIX}
		Boost::filesystem${BOOST_SUFFIX} Boost::program_options${BOOST_SUFFIX}
		${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES}
	)

	if(CMAKE_BUILD_TYPE STREQUAL "Release" AND USE_STRIP)
		add_custom_command(TARGET takin_convoseries POST_BUILD
			COMMAND strip -v $<TARGET_FILE:takin_convoseries>
			MAIN_DEPENDENCY takin_convoseries
	)
	endif()
	# -----------------------------------------------------------------------------
endif()
