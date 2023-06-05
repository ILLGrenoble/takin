#!/bin/bash
# prebuilding steps
# @author Tobias Weber <tobias.weber@tum.de>
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

BUILD_DEV_DOC=0


echo -e "Takin root dir: ${TAKINROOT}\n"


if [ $BUILD_DEV_DOC -ne 0 ]; then
	echo -e "--------------------------------------------------------------------------------"
	echo -e "building devel docs..."

	cd ${TAKINROOT}
	doxygen takin-doc

	echo -e "--------------------------------------------------------------------------------\n"
fi
