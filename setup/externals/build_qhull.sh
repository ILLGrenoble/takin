#!/bin/bash
#
# builds qhull
# @author Tobias Weber <tweber@ill.fr>
# @date may-2024
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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



#QHULL_REMOTE=https://codeload.github.com/qhull/qhull/zip/refs/heads/master
QHULL_REMOTE=https://github.com/qhull/qhull/archive/refs/tags/v8.0.2.zip
QHULL_LOCAL_ZIP=${QHULL_REMOTE##*[/\\]}
QHULL_LOCAL=qhull-8.0.2


rm -f "${QHULL_LOCAL}"


if ! wget ${QHULL_REMOTE}; then
	echo -e "Could not download ${QHULL_REMOTE}."
	exit -1
fi


rm -rf "${QHULL_LOCAL}"
unzip "${QHULL_LOCAL_ZIP}"
cd "${QHULL_LOCAL}"


if [ $BUILD_FOR_MINGW -ne 0 ]; then
	mkdir build_lib && cd build_lib
	mingw64-cmake -DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
		-DCMAKE_POSITION_INDEPENDENT_CODE=True \
		-DQHULL_ENABLE_TESTING=False -DBUILD_APPLICATIONS=False \
		-DBUILD_STATIC_LIBS=True -DBUILD_SHARED_LIBS=False ..
	mingw64-make -j${NUM_CORES} && sudo mingw64-make install/strip
else
	cmake -DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_POSITION_INDEPENDENT_CODE=True \
		-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
		-DQHULL_ENABLE_TESTING=False -DBUILD_APPLICATIONS=False \
		-DBUILD_STATIC_LIBS=True -DBUILD_SHARED_LIBS=False \
		-B build_lib .
	cmake --build build_lib --parallel ${NUM_CORES} && sudo cmake --install build_lib --strip
fi
