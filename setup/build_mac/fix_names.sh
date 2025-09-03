#!/bin/bash
#
# @date 22-mar-17
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# makes the framework libraries locally usable
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

PRG="takin.app"

fix_libs=1

NAME_TOOL=install_name_tool
STRIP=strip

PY_VER=3.13
echo -e "Py version: ${PY_VER}"

# files whose linkage is to be changed
declare -a filestochange=(
	"${PRG}/Contents/Frameworks/QtCore.framework/Versions/Current/QtCore"
	"${PRG}/Contents/Frameworks/QtWidgets.framework/Versions/Current/QtWidgets"
	"${PRG}/Contents/Frameworks/QtGui.framework/Versions/Current/QtGui"
	"${PRG}/Contents/Frameworks/QtConcurrent.framework/Versions/Current/QtConcurrent"
	"${PRG}/Contents/Frameworks/QtOpenGL.framework/Versions/Current/QtOpenGL"
	"${PRG}/Contents/Frameworks/QtOpenGLWidgets.framework/Versions/Current/QtOpenGLWidgets"
	"${PRG}/Contents/Frameworks/QtPrintSupport.framework/Versions/Current/QtPrintSupport"
	"${PRG}/Contents/Frameworks/QtDBus.framework/Versions/Current/QtDBus"
	"${PRG}/Contents/Frameworks/QtSvg.framework/Versions/Current/QtSvg"
	"${PRG}/Contents/Frameworks/QtXml.framework/Versions/Current/QtXml"
	"${PRG}/Contents/Frameworks/QtXmlPatterns.framework/Versions/5/QtXmlPatterns"
	"${PRG}/Contents/Frameworks/qwt.framework/Versions/6/qwt"
	"${PRG}/Contents/MacOS/takin"
	"${PRG}/Contents/MacOS/magpie"
	"${PRG}/Contents/MacOS/takin_cif2xml"
	"${PRG}/Contents/MacOS/takin_findsg"
	"${PRG}/Contents/MacOS/takin_bz"
	"${PRG}/Contents/MacOS/takin_structfact"
	"${PRG}/Contents/MacOS/takin_magstructfact"
	"${PRG}/Contents/MacOS/takin_scanbrowser"
	"${PRG}/Contents/MacOS/takin_magsgbrowser"
	"${PRG}/Contents/MacOS/takin_moldyn"
	"${PRG}/Contents/MacOS/takinmod_py"
	"${PRG}/Contents/MacOS/takinmod_jl"
	"${PRG}/Contents/MacOS/takin_convofit"
	"${PRG}/Contents/MacOS/takin_convoseries"
	"${PRG}/Contents/MacOS/takin_polextract"
	"${PRG}/Contents/PlugIns/libmagnonmod.dylib"
)


# original symbols, pattern-matched
declare -a changefrom=(
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtCore\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtGui\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtWidgets\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtOpenGL \""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtOpenGLWidgets\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtConcurrent\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtXml\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtXmlPatterns\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtSvg\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtPrintSupport\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/QtDBus\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"[-_/@.a-zA-Z0-9]*/qwt\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/Python\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_system.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_filesystem.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_container.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_graph.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_random.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_regex.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_atomic.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_iostreams.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_program_options.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libboost_python313.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libMinuit2.0.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libjpeg.9.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libpng16.16.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libtiff.6.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libfreetype.6.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libgcc_s.1.1.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libgomp.1.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libquadmath.0.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libgfortran.5.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/liblapacke.3.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/liblapack.3.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libblas.3.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libopenblas.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libgthread-2.0.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libglib-2.0.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libpcre.1.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libpcre2-16.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libpcre2-8.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libintl.8.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libzstd.1.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/liblzma.5.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libqhull_r.8.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libhdf5_cpp.310.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libhdf5.310.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libsz.2.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libaec.0.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libicudata.77.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libicui18n.77.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libicuuc.77.dylib\""
	"otool -L __BIN_FILE__ | grep -E -o -m1 \"(/|@rpath)[-_/@.a-zA-Z0-9]*/libqcustomplot[_a-zA-Z0-9]*.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libmd4c.0.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libgemmi_cpp.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libharfbuzz.0.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libgraphite2.3.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libdbus-1.3.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libdouble-conversion.3.dylib\""
	"otool -L __BIN_FILE__ | grep -o -m1 \"/[-_/@.a-zA-Z0-9]*/libb2.1.dylib\""
)

# symbols to change into
declare -a changeto=(
	"@executable_path/../Frameworks/QtCore.framework/Versions/Current/QtCore"
	"@executable_path/../Frameworks/QtGui.framework/Versions/Current/QtGui"
	"@executable_path/../Frameworks/QtWidgets.framework/Versions/Current/QtWidgets"
	"@executable_path/../Frameworks/QtOpenGL.framework/Versions/Current/QtOpenGL"
	"@executable_path/../Frameworks/QtOpenGLWidgets.framework/Versions/Current/QtOpenGLWidgets"
	"@executable_path/../Frameworks/QtConcurrent.framework/Versions/Current/QtConcurrent"
	"@executable_path/../Frameworks/QtXml.framework/Versions/Current/QtXml"
	"@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/Current/QtXmlPatterns"
	"@executable_path/../Frameworks/QtSvg.framework/Versions/Current/QtSvg"
	"@executable_path/../Frameworks/QtPrintSupport.framework/Versions/Current/QtPrintSupport"
	"@executable_path/../Frameworks/QtDBus.framework/Versions/Current/QtDBus"
	"@executable_path/../Frameworks/qwt.framework/Versions/6/qwt"
	"@executable_path/../Frameworks/Python.framework/Versions/${PY_VER}/Python"
	"@executable_path/../Libraries/libboost_system.dylib"
	"@executable_path/../Libraries/libboost_filesystem.dylib"
	"@executable_path/../Libraries/libboost_container.dylib"
	"@executable_path/../Libraries/libboost_graph.dylib"
	"@executable_path/../Libraries/libboost_random.dylib"
	"@executable_path/../Libraries/libboost_regex.dylib"
	"@executable_path/../Libraries/libboost_atomic.dylib"
	"@executable_path/../Libraries/libboost_iostreams.dylib"
	"@executable_path/../Libraries/libboost_program_options.dylib"
	"@executable_path/../Libraries/libboost_python313.dylib"
	"@executable_path/../Libraries/libMinuit2.0.dylib"
	"@executable_path/../Libraries/libjpeg.9.dylib"
	"@executable_path/../Libraries/libpng16.16.dylib"
	"@executable_path/../Libraries/libtiff.6.dylib"
	"@executable_path/../Libraries/libfreetype.6.dylib"
	"@executable_path/../Libraries/libgcc_s.1.1.dylib"
	"@executable_path/../Libraries/libgomp.1.dylib"
	"@executable_path/../Libraries/libquadmath.0.dylib"
	"@executable_path/../Libraries/libgfortran.5.dylib"
	"@executable_path/../Libraries/liblapacke.3.dylib"
	"@executable_path/../Libraries/liblapack.3.dylib"
	"@executable_path/../Libraries/libblas.3.dylib"
	"@executable_path/../Libraries/libopenblas.0.dylib"
	"@executable_path/../Libraries/libgthread-2.0.0.dylib"
	"@executable_path/../Libraries/libglib-2.0.0.dylib"
	"@executable_path/../Libraries//libpcre.1.dylib"
	"@executable_path/../Libraries/libpcre2-16.0.dylib"
	"@executable_path/../Libraries/libpcre2-8.0.dylib"
	"@executable_path/../Libraries/libintl.8.dylib"
	"@executable_path/../Libraries/libzstd.1.dylib"
	"@executable_path/../Libraries/liblzma.5.dylib"
	"@executable_path/../Libraries/libqhull_r.8.0.dylib"
	"@executable_path/../Libraries/libhdf5_cpp.310.dylib"
	"@executable_path/../Libraries/libhdf5.310.dylib"
	"@executable_path/../Libraries/libsz.2.dylib"
	"@executable_path/../Libraries/libaec.0.dylib"
	"@executable_path/../Libraries/libicudata.77.dylib"
	"@executable_path/../Libraries/libicui18n.77.dylib"
	"@executable_path/../Libraries/libicuuc.77.dylib"
	"@executable_path/../Libraries/libqcustomplot.dylib"
	"@executable_path/../Libraries/libmd4c.0.dylib"
	"@executable_path/../Libraries/libgemmi_cpp.dylib"
	"@executable_path/../Libraries/libharfbuzz.0.dylib"
	"@executable_path/../Libraries/libgraphite2.3.dylib"
	"@executable_path/../Libraries/libdbus-1.3.dylib"
	"@executable_path/../Libraries/libdouble-conversion.3.dylib"
	"@executable_path/../Libraries/libb2.1.dylib"
)

CNT=$(expr ${#changefrom[*]} - 1)


function fix_name()
{
	cfile=$1

	if [ ! -e ${cfile} ]; then
		echo -e "Warning: ${cfile} does not exist."
		return
	fi

	echo -e "Processing binary \"${cfile}\"..."
	chmod a+rx ${cfile}

	for idx in $(seq 0 ${CNT}); do
		cfrom_interp=${changefrom[$idx]}
		cto=${changeto[$idx]}
		cfrom_cmd="${cfrom_interp/__BIN_FILE__/${cfile}}"
		cfrom=$(eval "${cfrom_cmd}")
		if [ -z "${cfrom}" ]; then
			continue
		fi

		#echo -e "Command: $cfrom_cmd"
		echo -e "\tChanging \"${cfrom}\"\n\t -> \"${cto}\"."
		chmod u+w ${cfile}
		${NAME_TOOL} -change ${cfrom} ${cto} ${cfile}
	done

	${STRIP} ${cfile}
	echo -e ""
}


# fix names for given files
for cfile in ${filestochange[@]}; do
	fix_name $cfile
done


if [ $fix_libs -ne 0 ]; then
	# fix names for libraries
	find ${PRG}/Contents/Libraries \( -name "*.dylib" -o -name "*.so" \) -print0 \
		| while read -d $'\0' dylib; do
		fix_name $dylib
	done
fi
