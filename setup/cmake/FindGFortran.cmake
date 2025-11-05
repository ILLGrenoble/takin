#
# finds the gfortran library
# @author Tobias Weber <tweber@ill.fr>
# @date 9-oct-2024
# @license GPLv3
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

find_library(GFortran_LIBRARIES
	NAMES gfortran libgfortran.so.5
	HINTS /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib
		/usr/local/opt/gcc/lib/gcc/current /usr/local/Cellar/gcc/*/lib/gcc/current
		/opt/homebrew/Cellar/gcc/*/lib/gcc/current/
	DOC "GFortran libraries"
)


if(GFortran_LIBRARIES)
	set(GFortran_FOUND TRUE)

	message("GFortran library: ${GFortran_LIBRARIES}")
else()
	set(GFortran_FOUND FALSE)
	set(GFortran_LIBRARIES "")

	message("GFortran could not be found.")
endif()
