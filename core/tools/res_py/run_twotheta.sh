#!/bin/bash
#
# calculates the instrumental resolution for several 2theta values
#
# @author Tobias Weber <tweber@ill.fr>
# @date jul-2024
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


# output directory
mkdir -p results


# 2theta parameter range
twotheta_start=10
twotheta_end=90
twotheta_step=0.5


# run for several 2theta values
for twotheta in $(LC_ALL=C seq ${twotheta_start} ${twotheta_step} ${twotheta_end})
do
    echo -e "--------------------------------------------------------------------------------"
    echo -e "Running calculation for 2theta = ${twotheta}..."
    echo -e "--------------------------------------------------------------------------------"

    # output file name for results
    results_file="results/tt_${twotheta}_deg.dat"

    # run the resolution calculation
    if ! python3 ./calc_reso.py -i in20fc \
        --ki 2.662 --kf 2.662 --twotheta ${twotheta} \
        -o "${results_file}"
    then
        echo -e "Error: Failed to run calculation program."
        exit -1
    fi

    echo -e "--------------------------------------------------------------------------------\n"
done
