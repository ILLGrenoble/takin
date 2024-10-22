#
# finds the lapacke libs
# @author Tobias Weber <tweber@ill.fr>
# @date 8-aug-2020
# @license GPLv3
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
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

find_path(Lapacke_INCLUDE_DIRS
	NAMES lapacke.h
	PATH_SUFFIXES lapacke lapacke/include
	HINTS /usr/local/include/ /usr/include/ /opt/local/include /usr/local/opt/lapack/include /usr/local/Cellar/lapack/*/include
	DOC "Lapacke include directories"
)


find_library(Lapacke_LIBRARIES
	NAMES lapacke
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/local/opt/lapack/lib /usr/local/Cellar/lapack/*/lib
	DOC "Lapacke libraries"
)


find_library(Lapack_LIBRARIES
	NAMES lapack
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/local/opt/lapack/lib /usr/local/Cellar/lapack/*/lib
	DOC "Lapack libraries"
)


find_library(BLAS_LIBRARIES
	NAMES blas
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/local/opt/lapack/lib /usr/local/Cellar/lapack/*/lib
	DOC "BLAS libraries"
)


if(Lapacke_INCLUDE_DIRS AND Lapacke_LIBRARIES AND Lapack_LIBRARIES)
	set(Lapacke_FOUND TRUE)

	message("Lapacke include directories: ${Lapacke_INCLUDE_DIRS}")
	message("Lapack library: ${Lapack_LIBRARIES}")
	message("Lapacke library: ${Lapacke_LIBRARIES}")
else()
	set(Lapacke_FOUND FALSE)

	message("Error: Lapacke could not be found!")
endif()


if(BLAS_LIBRARIES)
	set(BLAS_FOUND TRUE)

	message("BLAS library: ${BLAS_LIBRARIES}")
else()
	set(BLAS_FOUND FALSE)

	message("Error: BLAS could not be found!")
endif()
