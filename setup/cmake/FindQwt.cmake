#
# finds the qwt libs
# @author Tobias Weber <tweber@ill.fr>
# @date 2016
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

find_path(QWT_INCLUDE_DIRS
	NAMES qwt.h
	PATH_SUFFIXES qwt6-qt6 qwt-qt6 qt6/qwt qt6/qwt/qwt qt6/qwt6 qt6 qwt6 qwt
	HINTS /usr/local/qwt-6.3.0/include /usr/local/include /usr/include /opt/local/include /usr/local/Cellar/qwt/*/lib/qwt.framework/Versions/6/Headers/
	DOC "Qwt include directories"
)

list(APPEND QWT_INCLUDE_DIRS "${QWT_INCLUDE_DIRS}/..")


find_library(QWT_LIBRARIES
	NAMES qwt6-qt6 qwt-qt6 qwt6 qwt
	HINTS /usr/local/qwt-6.3.0/lib /usr/local/lib64 /usr/local/lib /usr/lib64 /usr/lib /opt/local/lib /usr/lib32 /usr/local/lib32
	DOC "Qwt libraries"
)

if(QWT_INCLUDE_DIRS AND QWT_LIBRARIES)
	set(QWT_FOUND TRUE)

	message("Qwt include directories: ${QWT_INCLUDE_DIRS}")
	message("Qwt library: ${QWT_LIBRARIES}")
else()
	set(QWT_FOUND FALSE)

	message("Error: Qwt could not be found!")
endif()
