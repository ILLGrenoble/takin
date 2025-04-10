/**
 * monte carlo convolution tool -- common functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date sep-2020
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include "monteconvo_common.h"
#include "libs/version.h"
#include "tlibs/time/chrono.h"


ResoFocus get_reso_focus(int iFocMono, int iFocAna)
{
	unsigned ifocMode = unsigned(ResoFocus::FOC_UNCHANGED);
	if(iFocMono == 1)
		ifocMode |= unsigned(ResoFocus::FOC_MONO_FLAT); // flat
	else if(iFocMono == 2)
		ifocMode |= unsigned(ResoFocus::FOC_MONO_H);    // horizontal
	else if(iFocMono == 3)
		ifocMode |= unsigned(ResoFocus::FOC_MONO_V);    // vertical
	else if(iFocMono == 4)
		ifocMode |= unsigned(ResoFocus::FOC_MONO_V)|unsigned(ResoFocus::FOC_MONO_H);  // both
	if(iFocAna == 1)
		ifocMode |= unsigned(ResoFocus::FOC_ANA_FLAT);  // flat
	else if(iFocAna == 2)
		ifocMode |= unsigned(ResoFocus::FOC_ANA_H);     // horizontal
	else if(iFocAna == 3)
		ifocMode |= unsigned(ResoFocus::FOC_ANA_V);     // vertical
	else if(iFocAna == 4)
		ifocMode |= unsigned(ResoFocus::FOC_ANA_V)|unsigned(ResoFocus::FOC_ANA_H);    // both

	return ResoFocus(ifocMode);
}


bool load_scan_file(const std::string& _strFile, Scan& scan,
	bool bFlipAxis, bool bAllowScanMerging,
	const Filter& filter)
{
	std::string strFile = _strFile;
	tl::trim(strFile);
	if(strFile == "")
		return false;

	std::vector<std::string> vecFiles;
	tl::get_tokens<std::string, std::string>(strFile, ";", vecFiles);
	std::for_each(vecFiles.begin(), vecFiles.end(), [](std::string& str){ tl::trim(str); });

	bool bLoaded = ::load_file(vecFiles, scan, true, filter, bFlipAxis, bAllowScanMerging);

	// if file was not found, alternatively look in global paths
	if(!bLoaded)
	{
		const std::vector<std::string>& vecGlobPaths = get_global_paths();
		for(const std::string& strGlobPath : vecGlobPaths)
		{
			std::vector<std::string> _vecFiles;
			for(const std::string& _strFile : vecFiles)
				_vecFiles.push_back(strGlobPath + "/" + _strFile);

			if((bLoaded = ::load_file(_vecFiles, scan, true, filter, bFlipAxis, bAllowScanMerging)))
				break;
		}
	}

    return bLoaded;
}


void write_takin_metadata(std::ostream& ostrOut)
{
	const char* pcUser = std::getenv("USER");
	if(!pcUser)
		pcUser = "";

	ostrOut << "# Takin/Monteconvo version " << TAKIN_VER << "\n";
	ostrOut << "# DOI: https://dx.doi.org/10.5281/zenodo.4117437\n";
	ostrOut << "# URL: https://github.com/ILLGrenoble/takin\n";
	ostrOut << "# Timestamp: " << tl::epoch_to_str(tl::epoch()) << "\n";
	ostrOut << "# User: " << pcUser << "\n";
}
