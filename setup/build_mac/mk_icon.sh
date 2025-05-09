#!/bin/bash
#
# @date oct-17
# @author Tobias Weber <tobias.weber@tum.de>
# @license GPLv2
#
# create icon files
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

FILE=res/icons/takin.svg


echo "svg -> png (512)"
magick convert -resize 512x512 -antialias -channel rgba -background "#ffffff00" -alpha background -transparent "#ffffffff" $FILE ${FILE%\.svg}.png

echo "png -> icns"
makeicns -in ${FILE%\.svg}.png -out ${FILE%\.svg}.icns


echo "svg -> png (128)"
magick convert -resize 128x128 -antialias -channel rgba -background "#ffffff00" -alpha background -transparent "#ffffffff" $FILE ${FILE%\.svg}.png

echo "png -> ico"
magick convert ${FILE%\.svg}.png ${FILE%\.svg}.ico


# remove temporary file
rm ${FILE%\.svg}.png
