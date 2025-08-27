#!/bin/bash
#
# builds qwt
# @author Tobias Weber <tweber@ill.fr>
# @date aug-2025
# @license see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

NUM_CORES=$(nproc)


BUILD_FOR_MINGW=0
if [ "$1" == "--mingw" ]; then
	BUILD_FOR_MINGW=1
fi


QWT_REMOTE=http://downloads.sourceforge.net/project/qwt/qwt/6.3.0/qwt-6.3.0.tar.bz2
QWT_DIR=qwt-6.3.0


QWT_LOCAL=${QWT_REMOTE##*[/\\]}
rm -f "${QWT_LOCAL}"


if ! wget ${QWT_REMOTE}; then
	echo -e "Could not download ${QWT_REMOTE}."
	exit -1
fi


rm -rf ${QWT_DIR}
tar -xjvf "${QWT_LOCAL}"
cd ${QWT_DIR}


if [ $BUILD_FOR_MINGW -ne 0 ]; then
	echo -e "TODO"
else
	qmake6 qwt.pro
	make -j${NUM_CORES} && sudo make install
fi
