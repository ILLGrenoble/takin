#!/bin/bash
#
# creates a DEB distro
# @author Tobias Weber <tobias.weber@tum.de>
# @date 14 may 2024
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


# installation directory
INSTDIR="$1"

if [ "${INSTDIR}" = "" ]; then
	INSTDIR=~/tmp/takin
fi


# directories
mkdir -p ${INSTDIR}/usr/local/bin
mkdir -p ${INSTDIR}/usr/local/lib
mkdir -p ${INSTDIR}/usr/local/lib/takin_plugins
mkdir -p ${INSTDIR}/usr/local/lib/python3.12/dist-packages/takin
mkdir -p ${INSTDIR}/usr/local/share/takin/res
mkdir -p ${INSTDIR}/usr/local/share/takin/3rdparty_licenses
mkdir -p ${INSTDIR}/usr/share/applications
mkdir -p ${INSTDIR}/DEBIAN


# control file
echo -e "Package: takin\nVersion: 2.9" > ${INSTDIR}/DEBIAN/control
echo -e "Description: inelastic neutron scattering software" >> ${INSTDIR}/DEBIAN/control
echo -e "Maintainer: n/a" >> ${INSTDIR}/DEBIAN/control
echo -e "Architecture: $(dpkg --print-architecture)" >> ${INSTDIR}/DEBIAN/control
echo -e "Section: science\nPriority: optional" >> ${INSTDIR}/DEBIAN/control
echo -e "Depends:" \
	"libstdc++6," \
	"libboost-system1.83.0 (>=1.83.0)," \
	"libboost-filesystem1.83.0 (>=1.83.0)," \
	"libboost-iostreams1.83.0 (>=1.83.0)," \
	"libboost-regex1.83.0 (>=1.83.0)," \
	"libboost-program-options1.83.0 (>=1.83.0)," \
	"libboost-python1.83.0 (>=1.83.0)," \
	"libqt6core6t64 (>=6.4.0)," \
	"libqt6gui6t64 (>=6.4.0)," \
	"libqt6opengl6t64 (>=6.4.0)," \
	"libqt6xml6t64 (>=6.4.0)," \
	"libqt6svg6 (>=6.4.0)," \
	"libpython3.12 (>=3.12.0)," \
	"python3.12 (>=3.12.0)," \
	"python3-numpy," \
	"python3-scipy," \
	"libfreetype6," \
	"gnuplot," \
	"gnuplot-qt," \
	"libopengl0," \
	"libgfortran5,"
	"liblapacke," \
	"libqhull-r8.0," \
	"libqhullcpp8.0," \
	"libhdf5-103-1t64," \
	"libhdf5-cpp-103-1t64," \
	"libqcustomplot2.1\n" \
		>> ${INSTDIR}/DEBIAN/control

# libqcustomplot2.0 is only needed by the external programs

# once it exists:
# "libqwt-qt6-6 (>=6.5.0)," \


# copy program files
cp -v bin/takin                            ${INSTDIR}/usr/local/bin
cp -v bin/takinmod_py                      ${INSTDIR}/usr/local/bin
cp -v bin/takinmod_jl                      ${INSTDIR}/usr/local/bin

cp -rv res/*                               ${INSTDIR}/usr/local/share/takin/res/
cp -rv demos                               ${INSTDIR}/usr/local/share/takin
cp -rv data/samples                        ${INSTDIR}/usr/local/share/takin
cp -rv data/instruments                    ${INSTDIR}/usr/local/share/takin
cp -v *.txt                                ${INSTDIR}/usr/local/share/takin
cp -v ../LICENSE                           ${INSTDIR}/usr/local/share/takin/LICENSE.txt
cp -v ../LICENSES                          ${INSTDIR}/usr/local/share/takin/LICENSES.txt
cp -v ../LITERATURE                        ${INSTDIR}/usr/local/share/takin/LITERATURE.txt
cp -v ../AUTHORS                           ${INSTDIR}/usr/local/share/takin/AUTHORS.txt
cp -rv 3rdparty_licenses/*                 ${INSTDIR}/usr/local/share/takin/3rdparty_licenses/
cp -v ../setup/build_lin/takin.desktop     ${INSTDIR}/usr/share/applications
cp -v /usr/local/lib/libMinuit2.so         ${INSTDIR}/usr/local/lib
cp -v /usr/local/lib/libMinuit2Math.so     ${INSTDIR}/usr/local/lib
cp -v /usr/local/qwt-6.3.0/lib/libqwt.so   ${INSTDIR}/usr/local/lib

# create some symbolic links for dynamically linked libraries
pushd ${INSTDIR}/usr/local/lib
if [ -e libMinuit2.so ]; then
	ln -sf libMinuit2.so libMinuit2.so.0
	ln -sf libMinuit2.so libMinuit2.so.0.0
	ln -sf libMinuit2.so libMinuit2.so.0.0.0
fi
if [ -e libMinuit2Math.so ]; then
	ln -sf libMinuit2Math.so libMinuit2Math.so.0
	ln -sf libMinuit2Math.so libMinuit2Math.so.0.0
	ln -sf libMinuit2Math.so libMinuit2Math.so.0.0.0
fi
if [ -e libqwt.so ]; then
	ln -sf libqwt.so libqwt.so.6
	ln -sf libqwt.so libqwt.so.6.3
	ln -sf libqwt.so libqwt.so.6.3.0
fi
popd


# copy external programs
cp -v bin/magpie                ${INSTDIR}/usr/local/bin
cp -v bin/takin_cif2xml         ${INSTDIR}/usr/local/bin
cp -v bin/takin_findsg          ${INSTDIR}/usr/local/bin
cp -v bin/takin_bz              ${INSTDIR}/usr/local/bin
cp -v bin/takin_structfact      ${INSTDIR}/usr/local/bin
cp -v bin/takin_magstructfact   ${INSTDIR}/usr/local/bin
cp -v bin/takin_scanbrowser     ${INSTDIR}/usr/local/bin
cp -v bin/takin_magsgbrowser    ${INSTDIR}/usr/local/bin
cp -v bin/takin_moldyn          ${INSTDIR}/usr/local/bin


# copy plugins
cp -v plugins/*.so              ${INSTDIR}/usr/local/lib/takin_plugins

# copy py modules
cp -v pymods/*                  ${INSTDIR}/usr/local/lib/python3.12/dist-packages/takin


# permissions
chmod a+x ${INSTDIR}/usr/local/bin/*

# stripping
strip -v ${INSTDIR}/usr/local/bin/*
strip -v ${INSTDIR}/usr/local/lib/*
strip -v ${INSTDIR}/usr/local/lib/takin_plugins/*
strip -v ${INSTDIR}/usr/local/lib/python3.12/dist-packages/takin/*.so


# startup script
cp -v takin.sh                  ${INSTDIR}/usr/local/bin


# remove unnecessary files
find ${INSTDIR} -type f -name ".dir" -exec rm -fv {} \; -print
find ${INSTDIR} -type f -name "*.pyc" -exec rm -fv {} \; -print
rm -v ${INSTDIR}/usr/local/share/takin/CMakeLists.txt


# build package
cd ${INSTDIR}/..
chmod -R 775 ${INSTDIR}
dpkg --build takin
