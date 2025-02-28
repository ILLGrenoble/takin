#
# @author Tobias Weber <tweber@ill.fr>
# @date 11-apr-2020
# @license: see 'LICENSE' file
# @desc Forked on 18/July/2020 from my privately developed "matrix_calc" project (https://github.com/t-weber/matrix_calc).
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
# matrix_calc
# Copyright (C) 2020       Tobias WEBER (privately developed).
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
project(mcalc)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/../../setup/cmake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(FLEX 2 REQUIRED)
find_package(BISON 3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(GFortran)
find_package(Lapacke REQUIRED)

add_compile_options(-Wall -Wextra)
add_compile_options(-std=c++20)

add_definitions(${Boost_CXX_FLAGS})
add_definitions(-DUSE_LAPACK)


include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/libs"
	"${PROJECT_SOURCE_DIR}/.."
	"${FLEX_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}/.."
	"${CMAKE_CURRENT_BINARY_DIR}"
	"${Lapacke_INCLUDE_DIRS}"
)


BISON_TARGET(parser_impl
	src/parser.y ${CMAKE_CURRENT_BINARY_DIR}/parser_impl.cpp
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser_defs.h
	VERBOSE REPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser_report.txt
)


FLEX_TARGET(lexer_impl
	src/lexer.l ${CMAKE_CURRENT_BINARY_DIR}/lexer_impl.cpp
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/lexer_impl.h
)

ADD_FLEX_BISON_DEPENDENCY(lexer_impl parser_impl)



message("Flex outputs: ${FLEX_lexer_impl_OUTPUTS}")
message("Bison output source: ${BISON_parser_impl_OUTPUT_SOURCE}")
message("Bison output header: ${BISON_parser_impl_OUTPUT_HEADER}")



add_executable(mcalc
	src/parser.cpp src/parser.h src/ast.h src/sym.h
	src/llasm.cpp src/llasm_ops.cpp src/llasm_var.cpp
	src/llasm_arr.cpp src/llasm_func.cpp src/llasm.h
	src/printast.cpp src/printast.h
	${FLEX_lexer_impl_OUTPUTS}
	${BISON_parser_impl_OUTPUT_SOURCE} ${BISON_parser_impl_OUTPUT_HEADER}
)

#add_dependencies(parser parser_impl lexer_impl)
target_link_libraries(mcalc ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})




add_library(mcalc_rt SHARED
	src/runtime.cpp
)

set(TARGET mcalc_rt PROPTERTY CMAKE_POSITION_INDEPENDENT_CODE ON)

target_link_libraries(mcalc_rt ${Lapacke_LIBRARIES} ${Lapack_LIBRARIES} ${BLAS_LIBRARIES} ${GFortran_LIBRARIES})

