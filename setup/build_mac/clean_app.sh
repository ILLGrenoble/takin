#!/bin/bash
#
# removes unnecessary files from the app dir
# @author Tobias Weber <tobias.weber@tum.de>
# @date 16-dec-2020
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

aggressive_py_cleaning=1


PRG="takin.app"

DST_FRAMEWORK_DIR="${PRG}/Contents/Frameworks/"
DST_SITEPACKAGES_DIR="${PRG}/Contents/site-packages/"
DST_PY_DIR="${DST_FRAMEWORK_DIR}/Python.framework/"


# ----------------------------------------------------------------------------
# general cleaning
# ----------------------------------------------------------------------------
# clean temporary files
find ${PRG} -type d -name "__pycache__" -exec rm -rfv {} \;
find ${PRG} -type f -name ".dir" -exec rm -fv {} \; -print
find ${PRG} -type f -name ".DS_Store" -exec rm -fv {} \; -print
find ${PRG} -type f -name "*.o" -print0 | xargs -0 rm -rfv
find ${PRG} -type f -name "*.ico" -print0 | xargs -0 rm -rfv
find ${PRG} -type f -name ".*" -exec rm -fv {} \;


# clean non-needed files from frameworks
find "${DST_FRAMEWORK_DIR}" -type d -name "Headers" -print0 | xargs -0 rm -rfv
find "${DST_FRAMEWORK_DIR}" -type d -name "Current" -print0 | xargs -0 rm -rfv
find "${DST_FRAMEWORK_DIR}" -name "*.pyc" -print0 | xargs -0 rm -rfv

declare -a QT_FW_LIBS=("QtCore" "QtGui" "QtWidgets" "QtConcurrent" "QtOpenGL"
	"QtSvg" "QtXml" "QtDBus" "QtPrintSupport" "qwt")

for qlib in ${QT_FW_LIBS[@]}; do
	rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/Resources
	rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/Headers
	rm -rfv ${DST_FRAMEWORK_DIR}/${qlib}.framework/${qlib}*
done
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# clean py dir
# ----------------------------------------------------------------------------
# remove old signatures
codesign --remove-signature ${DST_PY_DIR}/Versions/Current/Python


# clean site-packages
rm -rfv ${DST_SITEPACKAGES_DIR}/PyQt5*
rm -rfv ${DST_SITEPACKAGES_DIR}/wheel*
rm -rfv ${DST_SITEPACKAGES_DIR}/setuptools*
rm -rfv ${DST_SITEPACKAGES_DIR}/pip*
rm -rfv ${DST_SITEPACKAGES_DIR}/distutils*
rm -rfv ${DST_SITEPACKAGES_DIR}/_distutils*
rm -rfv ${DST_SITEPACKAGES_DIR}/pkg_*
rm -rfv ${DST_SITEPACKAGES_DIR}/sip*
rm -rfv ${DST_SITEPACKAGES_DIR}/site*
rm -rfv ${DST_SITEPACKAGES_DIR}/easy_*


# remove tcl/tk
rm -fv ${DST_PY_DIR}/Versions/Current/lib/libtk*
rm -fv ${DST_PY_DIR}/Versions/Current/lib/libtcl*
rm -rfv ${DST_PY_DIR}/Versions/Current/lib/tcl8*
rm -rfv ${DST_PY_DIR}/Versions/Current/lib/tk8*
rm -rfv ${DST_PY_DIR}/Versions/Current/lib/itcl4*
rm -fv ${DST_PY_DIR}/Versions/Current/lib/pkgconfig/tcl.pc
rm -fv ${DST_PY_DIR}/Versions/Current/lib/tcl*
rm -fv ${DST_PY_DIR}/Versions/Current/lib/tcl*
rm -fv ${DST_PY_DIR}/Versions/Current/lib/tk*
rm -fv ${DST_PY_DIR}/Versions/Current/lib/pkgconfig/tk.pc
rm -rfv ${DST_PY_DIR}/Versions/Current/lib/python3.13/tkinter
rm -fv ${DST_PY_DIR}/Versions/Current/lib/python3.13/lib-dynload/_tkinter*.so
rm -fv ${DST_PY_DIR}/Versions/Current/lib/Tk*


# remove non-needed site packages
pushd ${PRG}/Contents
	ln -sf Frameworks/Python.framework/Versions/Current/lib/python3.13/site-packages
popd

rm -rfv ${DST_PY_DIR}/Versions/Current/lib/python3.13/test
rm -rfv ${DST_PY_DIR}/Versions/Current/lib/python3.13/idlelib

rm -rfv ${PRG}/Contents/site-packages/setuptools*
rm -rfv ${PRG}/Contents/site-packages/pip*
rm -rfv ${PRG}/Contents/site-packages/easy*
rm -rfv ${PRG}/Contents/site-packages/pkg_resources


# remove app
rm -rfv ${DST_PY_DIR}/Versions/Current/Resources/Python.app


if [ $aggressive_py_cleaning -ne 0 ]; then
	# remove everything except Python and lib
	rm -rfv ${DST_PY_DIR}/Versions/Current/bin
	rm -rfv ${DST_PY_DIR}/Versions/Current/etc
	rm -rfv ${DST_PY_DIR}/Versions/Current/share
	rm -rfv ${DST_PY_DIR}/Versions/Current/include
	rm -rfv ${DST_PY_DIR}/Versions/Current/Headers
	rm -rfv ${DST_PY_DIR}/Versions/Current/Resources
	rm -rfv ${DST_PY_DIR}/Versions/Current/_CodeSignature

	rm -fv ${DST_PY_DIR}/Resources
	rm -fv ${DST_PY_DIR}/Headers

	for file in ${DST_PY_DIR}/Versions/Current/lib/*; do
		# remove everything except file libpython3.xx.dylib and dir python3.xx
		if [[ $file =~ [a-zA-Z0-9\./]*python[a-zA-Z0-9\./]* ]]; then
			echo -e "Keeping \"${file}\"."
			continue
		fi

		rm -rfv $file
	done

	find ${DST_PY_DIR} -type d -name "tests" -exec rm -rfv {} \;
fi
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# strip binaries
# ----------------------------------------------------------------------------
find ${PRG} -name "*.so" -exec strip -v {} \; -print
find ${PRG} -name "*.dylib" -exec strip -v {} \; -print

strip -v ${PRG}/Contents/MacOS/*
# ----------------------------------------------------------------------------
