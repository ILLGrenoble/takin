#!/bin/bash
#
# gets takin extension modules -- to be run in the takin root directory
# @author Tobias Weber <tweber@ill.fr>
# @date 28-feb-2025
# @license GPLv2
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


#
# magpie project
#
function clone_magpie
{
	# test for existing symlink
	if [ -L magpie ]; then
		echo -e "Magpie project already symlinked..."
		return
	fi

	# test for existing directory
	if [ -d magpie ]; then
		echo -e "Updating Magpie project..."

		pushd magpie
		if ! (git pull); then
			echo -e "Error: Cannot pull Magpie project.";
			return -1;
		fi
		popd

		return
	fi

	# clone a new copy
	echo -e "Cloning Magpie project..."

	if ! (git clone --branch ver-0.9.2 https://github.com/ILLGrenoble/magpie); then
		echo -e "Error: Cannot clone Magpie project.";
		return -1;
	fi
}



#
# taspaths project
#
function clone_taspaths
{
	# test for existing symlink
	if [ -L taspaths ]; then
		echo -e "TAS-Paths project already symlinked..."
		return
	fi

	# test for existing directory
	if [ -d taspaths ]; then
		echo -e "Updating TAS-Paths project..."

		pushd taspaths
		if ! (git pull); then
			echo -e "Error: Cannot pull TAS-Paths project.";
			return -1;
		fi
		popd

		return
	fi

	# clone a new copy
	echo -e "Cloning TAS-Paths project..."

	if ! (git clone https://github.com/ILLGrenoble/taspaths); then
		echo -e "Error: Cannot clone TAS-Paths project.";
		return -1;
	fi
}



echo -e "--------------------------------------------------------------------------------"
clone_magpie
echo -e "--------------------------------------------------------------------------------"

#echo -e "--------------------------------------------------------------------------------"
#clone_taspaths
#echo -e "--------------------------------------------------------------------------------"
