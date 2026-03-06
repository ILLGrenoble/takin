#!/bin/bash
#
# @date 12-apr-17
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# copy framework libs
#
# ----------------------------------------------------------------------------
# Takin (inelastic neutron scattering software package)
# Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

PRG="takin.app"

BIN_DIR="${PRG}/Contents/MacOS/"
DST_FRAMEWORK_DIR="${PRG}/Contents/Frameworks/"
DST_LIB_DIR="${PRG}/Contents/Libraries/"
DST_PLUGIN_DIR="${PRG}/Contents/PlugIns/"
DST_SITEPACKAGES_DIR="${PRG}/Contents/site-packages"
BREW_DIR=/opt/homebrew
QT_DIR="${BREW_DIR}"
#QT_DIR=/usr/local/opt/qt
PLUGIN_DIR="${QT_DIR}/share/qt/plugins/"
GCC_DIR="/opt/homebrew/Cellar/gcc/15.2.0_1"


# -----------------------------------------------------------------------------
# frameworks
declare -a SRC_FRAMEWORKS=(
	"${QT_DIR}/lib/QtCore.framework"
	"${QT_DIR}/lib/QtGui.framework"
	"${QT_DIR}/lib/QtWidgets.framework"
	"${QT_DIR}/lib/QtOpenGL.framework"
	"${QT_DIR}/lib/QtOpenGLWidgets.framework"
	"${QT_DIR}/lib/QtConcurrent.framework"
	"${QT_DIR}/lib/QtXml.framework"
	"${QT_DIR}/lib/QtSvg.framework"
	"${QT_DIR}/lib/QtPrintSupport.framework"
	"${QT_DIR}/lib/QtDBus.framework"
	"/usr/local/qwt-6.3.0/lib/qwt.framework"
	"/Library/Frameworks/Python.framework"
)
#	"/usr/local/opt/python/Frameworks/Python.framework"

# libs which need to have their symbolic link resolved
declare -a SRC_LIBS=(
	"${BREW_DIR}/lib/libMinuit2.0.dylib"
	"${BREW_DIR}/lib/libboost_system.dylib"
	"${BREW_DIR}/lib/libboost_filesystem.dylib"
	"${BREW_DIR}/lib/libboost_container.dylib"
	"${BREW_DIR}/lib/libboost_graph.dylib"
	"${BREW_DIR}/lib/libboost_random.dylib"
	"${BREW_DIR}/lib/libboost_regex.dylib"
	"${BREW_DIR}/lib/libboost_atomic.dylib"
	"${BREW_DIR}/lib/libboost_iostreams.dylib"
	"${BREW_DIR}/lib/libboost_program_options.dylib"
	"${BREW_DIR}/lib/libboost_python314.dylib"
	"${BREW_DIR}/lib/libfreetype.6.dylib"
	"${BREW_DIR}/lib/libpng16.16.dylib"
	"${BREW_DIR}/lib/libjpeg.8.dylib"
	"${BREW_DIR}/lib/libtiff.6.dylib"
	"${BREW_DIR}/lib/liblapacke.3.dylib"
	"${BREW_DIR}/lib/liblapack.3.dylib"
	"${BREW_DIR}/lib/libblas.3.dylib"
	"${GCC_DIR}/lib/gcc/current/libgfortran.5.dylib"
	"${GCC_DIR}/lib/gcc/current/libquadmath.0.dylib"
	"${GCC_DIR}/lib/gcc/current/libgomp.1.dylib"
	"${GCC_DIR}/lib/gcc/current/libgcc_s.1.1.dylib"
	"${BREW_DIR}/lib/libgthread-2.0.0.dylib"
	"${BREW_DIR}/lib/libglib-2.0.0.dylib"
	"${BREW_DIR}/lib/libpcre.1.dylib"
	"${BREW_DIR}/lib/libpcre2-16.0.dylib"
	"${BREW_DIR}/lib/libpcre2-8.0.dylib"
	"${BREW_DIR}/lib/libintl.8.dylib"
	"${BREW_DIR}/lib/libzstd.1.dylib"
	"${BREW_DIR}/lib/liblzma.5.dylib"
	"${BREW_DIR}/lib/libqhull_r.8.0.dylib"
	"${BREW_DIR}/lib/libhdf5_cpp.320.dylib"
	"${BREW_DIR}/lib/libhdf5.320.dylib"
	"${BREW_DIR}/lib/libsz.2.dylib"
	"${BREW_DIR}/lib/libaec.0.dylib"
	"${BREW_DIR}/lib/libicudata.78.dylib"
	"${BREW_DIR}/lib/libicui18n.78.dylib"
	"${BREW_DIR}/lib/libicuuc.78.dylib"
	"${BREW_DIR}/lib/libqcustomplot.dylib"
	"${BREW_DIR}/lib/libmd4c.0.dylib"
	"${BREW_DIR}/lib/libgemmi_cpp.dylib"
	"${BREW_DIR}/lib/libharfbuzz.0.dylib"
	"${BREW_DIR}/lib/libgraphite2.3.dylib"
	"${BREW_DIR}/lib/libdbus-1.3.dylib"
	"${BREW_DIR}/lib/libdouble-conversion.3.dylib"
	"${BREW_DIR}/lib/libb2.1.dylib"
)
#	"/usr/local/opt/openblas/lib/libopenblas.0.dylib"


# qt plugins
declare -a SRC_PLUGINS=(
	"printsupport/libcocoaprintersupport.dylib"
	"imageformats/libqsvg.dylib"
	"imageformats/libqicns.dylib"
	"imageformats/libqjpeg.dylib"
	"iconengines/libqsvgicon.dylib"
	"styles/libqmacstyle.dylib"
	"platforms/libqcocoa.dylib"
	"platforms/libqminimal.dylib"
)
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# clean old files
rm -rfv ${BIN_DIR}
rm -rfv ${DST_FRAMEWORK_DIR}
rm -rfv ${DST_LIB_DIR}
rm -rfv ${DST_PLUGIN_DIR}
rm -rfv ${DST_SITEPACKAGES_DIR}
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# create dirs
mkdir -pv "${DST_FRAMEWORK_DIR}"
mkdir -pv "${DST_LIB_DIR}"
mkdir -pv "${DST_PLUGIN_DIR}/printsupport"
mkdir -pv "${DST_PLUGIN_DIR}/imageformats"
mkdir -pv "${DST_PLUGIN_DIR}/iconengines"
mkdir -pv "${DST_PLUGIN_DIR}/platforms"
mkdir -pv "${DST_PLUGIN_DIR}/platformthemes"
mkdir -pv "${DST_PLUGIN_DIR}/styles"
mkdir -pv "${BIN_DIR}"
mkdir -pv "${PRG}/Contents/Resources"
mkdir -pv "${DST_SITEPACKAGES_DIR}"
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# copy frameworks
for srclib in ${SRC_FRAMEWORKS[@]}; do
	echo -e "Copying framework \"${srclib}\"..."
	cp -rv ${srclib} "${DST_FRAMEWORK_DIR}"
done


# copy libs, following links
for srclib in ${SRC_LIBS[@]}; do
	echo -e "Copying library \"${srclib}\"..."
	cp -Lrv ${srclib} "${DST_LIB_DIR}"
done


# copy plugins
for srclib in ${SRC_PLUGINS[@]}; do
	echo -e "Copying plugin \"${srclib}\"..."
	cp -v "${PLUGIN_DIR}${srclib}" "${DST_PLUGIN_DIR}${srclib}"
	chmod a+rx "${DST_PLUGIN_DIR}${srclib}"
done


# copy convo plugins
cp -v plugins/*.dylib ${DST_PLUGIN_DIR}
chmod a+rx ${DST_PLUGIN_DIR}/*.dylib


# copy py modules
cp -rv pymods ${DST_PLUGIN_DIR}
cp -v ../setup/build_general/install_pymods.py  ${DST_PLUGIN_DIR}
chmod a+rx ${DST_PLUGIN_DIR}/pymods/*.so
chmod a+rx ${DST_PLUGIN_DIR}/pymods/*.py
chmod a+rx ${DST_PLUGIN_DIR}/install_pymods.py


# copy binaries
cp -v bin/takin "${BIN_DIR}"
cp -v bin/magpie "${BIN_DIR}"
cp -v bin/takin_cif2xml "${BIN_DIR}"
cp -v bin/takin_findsg "${BIN_DIR}"
cp -v bin/takin_bz "${BIN_DIR}"
cp -v bin/takinmod_py "${BIN_DIR}"
cp -v bin/takinmod_jl "${BIN_DIR}"
cp -v bin/takin_structfact "${BIN_DIR}"
cp -v bin/takin_magstructfact "${BIN_DIR}"
cp -v bin/takin_scanbrowser "${BIN_DIR}"
cp -v bin/takin_magsgbrowser "${BIN_DIR}"
cp -v bin/takin_moldyn "${BIN_DIR}"

cp -v bin/takin_convofit "${BIN_DIR}"
cp -v bin/takin_convoseries "${BIN_DIR}"
cp -v bin/takin_polextract "${BIN_DIR}"

# copy plugin modules
cp -v plugins/*.dylib "${DST_PLUGIN_DIR}"


# data files
cp -rv data/res "${PRG}/Contents/"
cp -rv doc/* "${PRG}/Contents/res/doc/"
cp -rv 3rdparty_licenses "${PRG}/Contents/Resources/"

cp -v data/res/icons/takin.icns "${PRG}/Contents/Resources/"
cp -v ../setup/build_mac/plists/Info.plist "${PRG}/Contents/"
cp -v *.txt "${PRG}/Contents/Resources/"
cp -v ../LICENSE "${PRG}/Contents/Resources/LICENSE.txt"
cp -v ../LICENSES "${PRG}/Contents/Resources/LICENSES.txt"
cp -v ../LITERATURE "${PRG}/Contents/Resources/LITERATURE.txt"
cp -v ../AUTHORS.md "${PRG}/Contents/Resources/AUTHORS.txt"

cp -rv data/instruments "${PRG}/Contents/Resources/"
cp -rv data/samples "${PRG}/Contents/Resources/"
cp -rv ../demos "${PRG}/Contents/Resources/"

rm -v "${PRG}/Contents/Resources/CMakeLists.txt"


# python site-packages
#cp -rv /usr/local/opt/numpy/lib/python3.9/site-packages/numpy "${DST_SITEPACKAGES_DIR}"
#cp -rv /usr/local/opt/scipy/lib/python3.9/site-packages/scipy "${DST_SITEPACKAGES_DIR}"
##ln -sf "../site-packages" "${BIN_DIR}/site-packages"
##ln -sf ./Frameworks/Python.framework/Versions/Current/lib/python3.9/site-packages ${DST_SITEPACKAGES_DIR}
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# attributes
chmod -R u+rw,a+r "${DST_FRAMEWORK_DIR}"
chmod -R u+rw,a+r "${DST_LIB_DIR}"
chmod -R u+rw,a+r "${DST_PLUGIN_DIR}"

# only the binaries in BIN_DIR should have the executable flag
find "${PRG}" -type f -print0 | xargs -0 chmod a-x
chmod a+x ${BIN_DIR}*

find "${DST_FRAMEWORK_DIR}" -type d -print0 | xargs -0 chmod a+x
find "${DST_LIB_DIR}" -type d -print0 | xargs -0 chmod a+x
find "${DST_PLUGIN_DIR}" -type d -print0 | xargs -0 chmod a+x
# -----------------------------------------------------------------------------
