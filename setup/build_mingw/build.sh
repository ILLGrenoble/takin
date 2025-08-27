#!/bin/bash
#
# takin mingw build script
# @author Tobias Weber <tweber@ill.fr>
# @date sep-2020
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

# individual building steps
setup_buildenv=1
setup_externals=1
build_externals=1
build_takin=1
build_magpie=1
build_plugins=1
build_package=1


NUM_CORES=$(nproc)
#NUM_CORES=1


# get root dir of takin repos
TAKIN_ROOT=$(dirname $0)/../..
cd "${TAKIN_ROOT}"
TAKIN_ROOT=$(pwd)
echo -e "Takin root dir: ${TAKIN_ROOT}"


if [ $setup_buildenv -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Setting up build environment..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/setup"
		if ! ./build_mingw/buildenv.sh; then
			exit -1
		fi
	popd
fi


if [ $setup_externals -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Getting external dependencies..."
	echo -e "================================================================================\n"

	if ! ./setup/externals/setup_modules.sh; then
		exit -1
	fi

	pushd "${TAKIN_ROOT}/core"
		rm -rf tmp
		if ! ../setup/externals/setup_externals.sh; then
			exit -1
		fi
		if ! ../setup/externals/get_3rdparty_licenses.sh; then
			exit -1
		fi
	popd
fi


if [ $build_externals -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building external libraries..."
	echo -e "================================================================================\n"

	mkdir -p "${TAKIN_ROOT}/tmp"
	pushd "${TAKIN_ROOT}/tmp"
		if ! "${TAKIN_ROOT}"/setup/externals/build_minuit.sh --mingw; then
			exit -1
		fi
		if ! "${TAKIN_ROOT}"/setup/externals/build_lapacke.sh --mingw; then
			exit -1
		fi
		if ! "${TAKIN_ROOT}"/setup/externals/build_qwt.sh --mingw; then
			exit -1
		fi
	popd
fi


if [ $build_takin -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building main Takin binary..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/core"
		../setup/build_general/clean.sh

		mkdir -p build
		cd build

		if ! mingw64-cmake -DCMAKE_BUILD_TYPE=Release -DDEBUG=False -DUSE_CUSTOM_THREADPOOL=True ..; then
			echo -e "Failed configuring core package."
			exit -1
		fi

		if ! mingw64-make -j${NUM_CORES}; then
			echo -e "Failed building core package."
			exit -1
		fi
	popd
fi


if [ $build_magpie -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Magpie tools..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/mag-core"
		# build external libraries
		if ! ./setup/externals/build_qhull.sh --mingw; then
			exit -1
		fi
		cp -v ./setup/externals/CMakeLists_qcp.txt .
		if ! ./setup/externals/build_qcp.sh --mingw; then
			exit -1
		fi
		if ! ./setup/externals/build_gemmi.sh --mingw; then
			exit -1
		fi

		# build tools
		rm -rf build
		mkdir -p build
		cd build

		if ! mingw64-cmake -DCMAKE_BUILD_TYPE=Release \
			-DONLY_BUILD_FINISHED=True -DBUILD_PY_MODULES=False ..; then
			echo -e "Failed configuring mag-core package."
			exit -1
		fi

		if ! mingw64-make -j${NUM_CORES}; then
			echo -e "Failed building mag-core package."
			exit -1
		fi

		# copy tools to Takin main dir
		cp -v build/magpie/magpie.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/cif2xml/takin_cif2xml.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/cif2xml/takin_findsg.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/bz/takin_bz.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/structfact/takin_structfact.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/structfact/takin_magstructfact.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/scanbrowser/takin_scanbrowser.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/magsgbrowser/takin_magsgbrowser.exe "${TAKIN_ROOT}"/core/bin/
		cp -v build/tools/moldyn/takin_moldyn.exe "${TAKIN_ROOT}"/core/bin/
	popd
fi


if [ $build_plugins -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Takin plugins..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}/magnon-plugin"
		rm -rf build
		mkdir -p build
		cd build

		if ! mingw64-cmake -DCMAKE_BUILD_TYPE=Release ..; then
			echo -e "Failed configuring magnon plugin."
			exit -1
		fi

		if ! mingw64-make -j${NUM_CORES}; then
			echo -e "Failed building magnon plugin."
			exit -1
		fi

		# copy plugin to Takin main dir
		cp -v libmagnonmod.dll "${TAKIN_ROOT}"/core/plugins/
	popd
fi


if [ $build_package -ne 0 ]; then
	echo -e "\n================================================================================"
	echo -e "Building Takin package..."
	echo -e "================================================================================\n"

	pushd "${TAKIN_ROOT}"
		rm -rf tmp-mingw
		cd core
		if ! ../setup/build_mingw/cp_mingw_takin.sh "${TAKIN_ROOT}/tmp-mingw/takin"; then
			exit -1
		fi

		cd ../tmp-mingw
		zip -9 -r takin.zip takin
	popd


	if [ -e  "${TAKIN_ROOT}/tmp-mingw/takin.zip" ]; then
		echo -e "\n================================================================================"
		echo -e "The built Takin package can be found here:\n\t${TAKIN_ROOT}/tmp-mingw/takin.zip"
		echo -e "================================================================================\n"
	else
		echo -e "\n================================================================================"
		echo -e "Error: Takin package could not be built!"
		echo -e "================================================================================\n"
	fi
fi
