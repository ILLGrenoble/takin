#!/bin/bash
#
# creates a DEB distro
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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
mkdir -p ${INSTDIR}/usr/local/share/takin/res
mkdir -p ${INSTDIR}/usr/local/share/takin/3rdparty_licenses
mkdir -p ${INSTDIR}/usr/share/applications
mkdir -p ${INSTDIR}/DEBIAN


# control file
echo -e "Package: takin\nVersion: 2.8.1" > ${INSTDIR}/DEBIAN/control
echo -e "Description: inelastic neutron scattering software" >> ${INSTDIR}/DEBIAN/control
echo -e "Maintainer: n/a" >> ${INSTDIR}/DEBIAN/control
echo -e "Architecture: $(dpkg --print-architecture)" >> ${INSTDIR}/DEBIAN/control
echo -e "Section: science\nPriority: optional" >> ${INSTDIR}/DEBIAN/control
echo -e "Depends: libstdc++6, libboost-system1.65.1 (>=1.65.1), libboost-filesystem1.65.1 (>=1.65.1), libboost-iostreams1.65.1 (>=1.65.1), libboost-regex1.65.1 (>=1.65.1), libboost-program-options1.65.1 (>=1.65.1), libboost-python1.65.1 (>=1.65.1), libqt5core5a (>=5.9.5), libqt5gui5 (>=5.9.5), libqt5opengl5 (>=5.9.5), libqt5svg5 (>=5.9.5), libqt5xml5 (>=5.9.5), libqwt-qt5-6 (>=6.1.3), libpython3.6 (>=3.6.0), python3.6 (>=3.6.0), python3-numpy, python3-scipy, libfreetype6, gnuplot, gnuplot-qt, libopengl0\n" >> ${INSTDIR}/DEBIAN/control


# copy program files
cp -v bin/takin			${INSTDIR}/usr/local/bin
#cp -v bin/takin_convofit	${INSTDIR}/usr/local/bin
#cp -v bin/takin_convoseries	${INSTDIR}/usr/local/bin
cp -v bin/takinmod_py		${INSTDIR}/usr/local/bin
cp -v bin/takinmod_jl		${INSTDIR}/usr/local/bin

cp -rv res/*			${INSTDIR}/usr/local/share/takin/res/
cp -rv 3rdparty_licenses/*	${INSTDIR}/usr/local/share/takin/3rdparty_licenses/
cp -v *.txt			${INSTDIR}/usr/local/share/takin
cp -v ../setup/build_lin/takin.desktop	${INSTDIR}/usr/share/applications
cp -v /usr/local/lib/libMinuit2.so	${INSTDIR}/usr/local/lib

pushd ${INSTDIR}/usr/local/lib
ln -sf libMinuit2.so libMinuit2.so.0
ln -sf libMinuit2.so libMinuit2.so.0.0
ln -sf libMinuit2.so libMinuit2.so.0.0.0
popd


# copy external programs
cp -v bin/magpie          	${INSTDIR}/usr/local/bin
cp -v bin/takin_cif2xml		${INSTDIR}/usr/local/bin
cp -v bin/takin_findsg		${INSTDIR}/usr/local/bin
cp -v bin/takin_bz		${INSTDIR}/usr/local/bin
cp -v bin/takin_structfact      ${INSTDIR}/usr/local/bin
cp -v bin/takin_magstructfact   ${INSTDIR}/usr/local/bin
cp -v bin/takin_scanbrowser     ${INSTDIR}/usr/local/bin
cp -v bin/takin_magsgbrowser    ${INSTDIR}/usr/local/bin
cp -v bin/takin_moldyn          ${INSTDIR}/usr/local/bin


# copy plugins
cp -v plugins/*.so		${INSTDIR}/usr/local/lib/takin_plugins


# permissions
chmod a+x ${INSTDIR}/usr/local/bin/*

# stripping
strip -v ${INSTDIR}/usr/local/bin/*
strip -v ${INSTDIR}/usr/local/lib/*
strip -v ${INSTDIR}/usr/local/lib/takin_plugins/*


# startup script
cp -v takin.sh			${INSTDIR}/usr/local/bin



# build package
cd ${INSTDIR}/..
chmod -R 775 ${INSTDIR}
dpkg --build takin
