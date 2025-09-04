#!/bin/bash
#
# gets license files for 3rd party libraries
# @author Tobias Weber <tweber@ill.fr>
# @date jan-2021
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

LICDIR=../core/3rdparty_licenses


echo -e "Downloading license texts...\n"

# boost
if ! wget http://www.boost.org/LICENSE_1_0.txt -O ${LICDIR}/boost_license.txt; then
	echo -e "Error: Cannot download Boost license.";
fi

# python
if ! wget https://raw.githubusercontent.com/python/cpython/master/Doc/license.rst -O ${LICDIR}/python_license.txt; then
	echo -e "Error: Cannot download Python license.";
fi

# numpy
if ! wget https://raw.githubusercontent.com/numpy/numpy/master/LICENSE.txt -O ${LICDIR}/numpy_license.txt; then
	echo -e "Error: Cannot download Numpy license.";
fi

# scipy
if ! wget https://raw.githubusercontent.com/scipy/scipy/master/LICENSE.txt -O ${LICDIR}/scipy_license.txt; then
	echo -e "Error: Cannot download Scipy license.";
fi

# swig
if ! wget https://raw.githubusercontent.com/swig/swig/master/LICENSE -O ${LICDIR}/swig_license.txt; then
	echo -e "Error: Cannot download swig license.";
fi

# ROOT
if ! wget https://raw.githubusercontent.com/root-project/root/master/LICENSE -O ${LICDIR}/ROOT_license.txt; then
	echo -e "Error: Cannot download ROOT license.";
fi

# qwt
if ! wget https://sourceforge.net/p/qwt/git/ci/develop/tree/COPYING\?format=raw -O ${LICDIR}/qwt_license.txt; then
	echo -e "Error: Cannot download Qwt license.";
fi

# freetype
if ! wget https://git.savannah.gnu.org/cgit/freetype/freetype2.git/plain/docs/FTL.TXT -O ${LICDIR}/freetype_license.txt; then
	echo -e "Error: Cannot download Freetype license.";
fi

# lapack(e)
if ! wget http://www.netlib.org/lapack/LICENSE.txt -O ${LICDIR}/lapack_license.txt; then
	echo -e "Error: Cannot download Lapack(e) license.";
fi

# qhull
if ! wget https://raw.githubusercontent.com/qhull/qhull/master/COPYING.txt -O ${LICDIR}/qhull_license.txt; then
	echo -e "Error: Cannot download Qhull license.";
fi

# neutronpy
if ! wget https://raw.githubusercontent.com/neutronpy/neutronpy/master/LICENSE.txt -O ${LICDIR}/neutronpy_license.txt; then
	echo -e "Error: Cannot download NeutronPy license.";
fi

# bodr
if ! wget https://raw.githubusercontent.com/egonw/bodr/master/bodr/COPYING -O ${LICDIR}/blue_obelisk_data_repository_license.txt; then
	echo -e "Error: Cannot download BODR license.";
fi

# dejavu
#if ! wget https://dejavu-fonts.github.io/License.html -O ${LICDIR}/dejavu_license.htm; then
if ! wget https://raw.githubusercontent.com/dejavu-fonts/dejavu-fonts/master/LICENSE -O ${LICDIR}/dejavu_license.txt; then
	echo -e "Error: Cannot download DejaVu license.";
fi

# libjpeg
if ! wget https://raw.githubusercontent.com/LuaDist/libjpeg/master/README -O ${LICDIR}/libjpeg_license.txt; then
#if ! wget https://raw.githubusercontent.com/freedesktop/libjpeg/master/README -O ${LICDIR}/libjpeg_license.txt; then
	echo -e "Error: Cannot download libjpg license.";
fi

# libtiff
#if ! wget http://www.libtiff.org/misc.html -O ${LICDIR}/libtiff_license.htm; then
if ! wget https://raw.githubusercontent.com/vadz/libtiff/master/COPYRIGHT -O ${LICDIR}/libtiff_license.txt; then
	echo -e "Error: Cannot download libtiff license.";
fi

# libpng
if ! wget http://www.libpng.org/pub/png/src/libpng-LICENSE.txt -O ${LICDIR}/libpng_license.txt; then
	echo -e "Error: Cannot download libpng license.";
fi

# gnuplot
if ! wget https://raw.githubusercontent.com/gnuplot/gnuplot/master/Copyright -O ${LICDIR}/gnuplot_license.txt; then
	echo -e "Error: Cannot download Gnuplot license.";
fi

# gemmi
if ! wget https://raw.githubusercontent.com/project-gemmi/gemmi/master/LICENSE.txt -O ${LICDIR}/gemmi_license.txt; then
	echo -e "Error: Cannot download Gemmi license.";
fi

# qcustomplot
if ! wget https://gitlab.com/ecme2/QCustomPlot/-/raw/master/GPL.txt -O ${LICDIR}/qcustomplot_license.txt; then
	echo -e "Error: Cannot download QCustomPlot license.";
fi

# xz
if ! wget https://git.tukaani.org/\?p=xz.git\;a=blob_plain\;f=COPYING -O ${LICDIR}/xz_license.txt; then
	echo -e "Error: Cannot download XZ license.";
fi

# zstd
if ! wget https://raw.githubusercontent.com/facebook/zstd/dev/LICENSE -O ${LICDIR}/zstd_license.txt; then
	echo -e "Error: Cannot download Zstandard license.";
fi

# libintl
if ! wget https://git.savannah.gnu.org/cgit/gettext.git/plain/gettext-runtime/intl/COPYING.LIB -O ${LICDIR}/libintl_license.txt; then
	echo -e "Error: Cannot download libintl license.";
fi

# pcre
if ! wget https://raw.githubusercontent.com/PCRE2Project/pcre2/master/LICENCE -O ${LICDIR}/pcre_license.txt; then
	echo -e "Error: Cannot download pcre license.";
fi

# glib
#if ! wget https://gitlab.gnome.org/GNOME/glib/-/raw/main/COPYING -O ${LICDIR}/glib_license.txt; then
if ! wget https://gitlab.gnome.org/GNOME/glib/-/raw/main/LICENSES/LGPL-2.1-or-later.txt -O ${LICDIR}/glib_license.txt; then
	echo -e "Error: Cannot download glib license.";
fi

# gcc runtime libraries
if ! wget https://raw.githubusercontent.com/gcc-mirror/gcc/master/COPYING.RUNTIME -O ${LICDIR}/gcc_runtime_license.txt; then
	echo -e "Error: Cannot download gcc_runtime license.";
fi

# libquadmath
if ! wget https://raw.githubusercontent.com/gcc-mirror/gcc/master/libquadmath/COPYING.LIB -O ${LICDIR}/libquadmath_license.txt; then
	echo -e "Error: Cannot download libquadmath license.";
fi

# hdf5
if ! wget https://raw.githubusercontent.com/HDFGroup/hdf5/develop/COPYING -O ${LICDIR}/hdf5_license.txt; then
	echo -e "Error: Cannot download hdf5 license.";
fi

# libaec / libsz2
if ! wget https://gitlab.dkrz.de/k202009/libaec/-/raw/master/LICENSE.txt -O ${LICDIR}/libaec_license.txt; then
	echo -e "Error: Cannot download libaec license.";
fi

# libharfbuzz
if ! wget https://raw.githubusercontent.com/harfbuzz/harfbuzz/main/COPYING -O ${LICDIR}/harfbuzz_license.txt; then
	echo -e "Error: Cannot download harfbuzz license.";
fi

# libdouble-conversion
if ! wget https://raw.githubusercontent.com/google/double-conversion/refs/heads/master/LICENSE -O ${LICDIR}/libdouble-conversion_license.txt; then
	echo -e "Error: Cannot download libdouble-conversion license.";
fi

# graphite
if ! wget https://raw.githubusercontent.com/silnrsi/graphite/refs/heads/master/LICENSE -O ${LICDIR}/graphite_license.txt; then
	echo -e "Error: Cannot download graphite license.";
fi

# dbus
if ! wget https://gitlab.freedesktop.org/dbus/dbus/-/raw/main/COPYING -O ${LICDIR}/dbus_license.txt; then
	echo -e "Error: Cannot download dbus license.";
fi

# blake2
if ! wget https://raw.githubusercontent.com/BLAKE2/BLAKE2/refs/heads/master/COPYING -O ${LICDIR}/blake2_license.txt; then
	echo -e "Error: Cannot download blake2 license.";
fi

# libbzip2
if ! wget https://raw.githubusercontent.com/libarchive/bzip2/master/COPYING -O ${LICDIR}/bzip2_license.txt; then
	echo -e "Error: Cannot download bzip2 license.";
fi

# zlib
if ! wget https://raw.githubusercontent.com/libarchive/zlib/master/LICENSE -O ${LICDIR}/zlib_license.txt; then
	echo -e "Error: Cannot download zlib license.";
fi

# libiconv
if ! wget https://git.savannah.gnu.org/cgit/libiconv.git/plain/COPYING.LIB -O ${LICDIR}/libiconv_license.txt; then
	echo -e "Error: Cannot download libiconv license.";
fi

# mingw64
if ! wget https://raw.githubusercontent.com/mingw-w64/mingw-w64/master/COPYING -O ${LICDIR}/mingw64_license.txt; then
	echo -e "Error: Cannot download mingw64 license.";
fi

# winpthreads
if ! wget https://raw.githubusercontent.com/mingw-w64/mingw-w64/master/mingw-w64-libraries/winpthreads/COPYING -O ${LICDIR}/winpthreads_license.txt; then
	echo -e "Error: Cannot download winpthreads license.";
fi
