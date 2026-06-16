#
# finds the julia libs
# @author Tobias Weber <tweber@ill.fr>
# @date 2016 / 2017
# @license GPLv2
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

find_path(Julia_INCLUDE_DIRS
	NAMES julia.h
	PATH_SUFFIXES julia
	HINTS /usr/local/include /usr/include /opt/local/includ ~/build/julia/usr/include
	DOC "Julia include directories"
)

list(APPEND Julia_INCLUDE_DIRS "${Julia_INCLUDE_DIRS}/..")


find_library(Julia_LIBRARIES
	NAMES julia
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32 ~/build/julia/usr/lib
	DOC "Julia libraries"
)

if(Julia_INCLUDE_DIRS AND Julia_LIBRARIES)
	set(Julia_FOUND TRUE)

	message("Julia include directories: ${Julia_INCLUDE_DIRS}")
	message("Julia library: ${Julia_LIBRARIES}")
else()
	set(Julia_FOUND FALSE)

	message("Error: Julia could not be found!")
endif()
