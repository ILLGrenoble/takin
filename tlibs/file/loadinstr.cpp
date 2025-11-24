/**
 * Loads instrument-specific data files
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
 * @license GPLv2 or GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include "loadinstr.h"
#include "loadinstr_impl.h"


namespace tl
{
	template FileInstrBase<double>* FileInstrBase<double>::LoadInstr(const char* pcFile);

	template class FilePsi<double>;
	template class FileFrm<double>;
	template class FileMacs<double>;
	template class FileTrisp<double>;
	template class FileRaw<double>;
#ifdef USE_HDF5
	template class FileILLH5<double>;
#endif


	template FileInstrBase<float>* FileInstrBase<float>::LoadInstr(const char* pcFile);

	template class FilePsi<float>;
	template class FileFrm<float>;
	template class FileMacs<float>;
	template class FileTrisp<float>;
	template class FileRaw<float>;
#ifdef USE_HDF5
	template class FileILLH5<float>;
#endif
}
