#!/bin/bash
#
# creates a distro for mingw
# @author Tobias Weber <tobias.weber@tum.de>
# @date 2016-2020
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


MINGW_BIN_DIR=/usr/x86_64-w64-mingw32/sys-root/mingw/bin
MINGW_LIB_DIR=/usr/x86_64-w64-mingw32/sys-root/mingw/lib


# installation directory
INSTDIR="$1"

if [ "${INSTDIR}" = "" ]; then
	INSTDIR=~/.wine/drive_c/takin
fi

mkdir -p ${INSTDIR}


# main programs
cp -v bin/*.exe			${INSTDIR}/


# info files
cp -v *.txt			${INSTDIR}/
cp -v ../LICENSE		${INSTDIR}/LICENSE.txt
cp -v ../LICENSES		${INSTDIR}/LICENSES.txt
cp -v ../LITERATURE		${INSTDIR}/LITERATURE.txt
cp -v ../AUTHORS		${INSTDIR}/AUTHORS.txt

cp -rv 3rdparty_licenses/	${INSTDIR}/
rm -v ${INSTDIR}/CMakeLists.txt


# examples
cp -rv examples			${INSTDIR}/
cp -rv data/samples		${INSTDIR}/
cp -rv data/instruments		${INSTDIR}/
cp -rv demos			${INSTDIR}/


# resources
mkdir ${INSTDIR}/res
cp -rv res/* ${INSTDIR}/res/
gunzip -v ${INSTDIR}/res/data/*



# libraries
cp -v $MINGW_BIN_DIR/libstdc++-6.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libwinpthread-1.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/libpng16-16.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libgcc_s_sjlj-1.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/libgcc_s_seh-1.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libbz2-1.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/zlib1.dll			${INSTDIR}/
cp -v $MINGW_BIN_DIR/iconv.dll			${INSTDIR}/
cp -v $MINGW_BIN_DIR/libpcre2-16-0.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libpcre2-8-0.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libpcre-1.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libharfbuzz-0.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libglib-2.0-0.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libintl-8.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libssp-0.dll		${INSTDIR}/
cp -v $MINGW_LIB_DIR/libMinuit2.dll		${INSTDIR}/
cp -v $MINGW_LIB_DIR/libMinuit2Math.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libMinuit2.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libMinuit2Math.dll		${INSTDIR}/

cp -v $MINGW_BIN_DIR/libboost_regex-x64.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_system-x64.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_iostreams-x64.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_filesystem-x64.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_program_options-x64.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_regex-mt-x64.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_system-mt-x64.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_iostreams-mt-x64.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_filesystem-mt-x64.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/libboost_program_options-mt-x64.dll	${INSTDIR}/

#cp -v $MINGW_BIN_DIR/libboost_python39.dll	${INSTDIR}/
#cp -v $MINGW_BIN_DIR/libpython3.9.dll		${INSTDIR}/
#cp -v $MINGW_BIN_DIR/libcrypto-1_1-x64.dll	${INSTDIR}/
#cp -v $MINGW_BIN_DIR/libffi-6.dll		${INSTDIR}/

cp -v $MINGW_BIN_DIR/Qt5Core.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/Qt5Gui.dll			${INSTDIR}/
cp -v $MINGW_BIN_DIR/Qt5Widgets.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/Qt5OpenGL.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/Qt5Svg.dll			${INSTDIR}/
cp -v $MINGW_BIN_DIR/Qt5Xml.dll			${INSTDIR}/
cp -v $MINGW_BIN_DIR/Qt5PrintSupport.dll	${INSTDIR}/
cp -v $MINGW_BIN_DIR/qwt-qt5.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libfreetype-6.dll		${INSTDIR}/

cp -v $MINGW_BIN_DIR/liblapack.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/liblapacke.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libblas.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libquadmath-0.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libgfortran-5.dll		${INSTDIR}/

cp -v $MINGW_BIN_DIR/libqcustomplot.dll		${INSTDIR}/
cp -v $MINGW_BIN_DIR/libqhull_r.dll		${INSTDIR}/


# qt plugins
mkdir -p ${INSTDIR}/lib/plugins/platforms/
mkdir -p ${INSTDIR}/lib/plugins/iconengines/
mkdir -p ${INSTDIR}/lib/plugins/imageformats

cp -v $MINGW_LIB_DIR/qt5/plugins/platforms/*.dll		${INSTDIR}/lib/plugins/platforms/
cp -v $MINGW_LIB_DIR/qt5/plugins/iconengines/qsvgicon.dll	${INSTDIR}/lib/plugins/iconengines/
cp -v $MINGW_LIB_DIR/qt5/plugins/imageformats/qsvg.dll		${INSTDIR}/lib/plugins/imageformats/


# convo plugins
cp -v plugins/*.dll ${INSTDIR}/lib/plugins/


# stripping
find ${INSTDIR}/ \( -name "*.exe" -o -name "*.dll" \) -exec strip -v {} \; -print


# access rights
find ${INSTDIR}/ -name "*.dll" -exec chmod a-x {} \; -print
find ${INSTDIR}/ -name "*.exe" -exec chmod a+x {} \; -print


echo -e "[Paths]\nPlugins = lib/plugins\n" > ${INSTDIR}/qt.conf
