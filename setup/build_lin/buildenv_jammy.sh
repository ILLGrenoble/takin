#!/bin/bash
#
# install all packages needed for building on jammy (and noble)
# @author Tobias Weber <tobias.weber@tum.de>
# @date 28-jul-20
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

#if [[ $(id -u) -gt 0 ]]; then
#	echo -e "Please run this script as root."
#	exit -1
#fi


# -----------------------------------------------------------------------------
# install packages
# -----------------------------------------------------------------------------
if ! sudo apt-get install cmake clang build-essential \
	libboost-all-dev libclipper-dev libhdf5-dev \
	qtchooser qt6-base-dev qt6-base-dev-tools qt6-svg-dev \
	libpython3-dev swig \
	libqcustomplot-dev \
	libgfortran5 \
	libfreetype6-dev libbz2-dev wget coreutils \
	flex bison
then
	echo -e "Error: Could not install packages necessary for building."
	exit -1
fi

# once it exists:
# libqwt-qt6-dev
# -----------------------------------------------------------------------------

#	g++-10 llvm llvm-dev \
#	liblapacke-dev libqhull-dev gemmi-dev
