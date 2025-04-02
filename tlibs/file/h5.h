/**
 * hdf5 loader
 * @author Tobias Weber <tweber@ill.fr>
 * @date 18-oct-2022
 * @license GPLv2 or GPLv3
 *
 * References:
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/readdata_8cpp-example.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/classH5_1_1CommonFG.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/classH5_1_1DataSet.html
 *   @see http://davis.lbl.gov/Manuals/HDF5-1.8.7/cpplus_RM/classH5_1_1DataSpace.html
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TLIBS_H5__
#define __TLIBS_H5__


#include <H5Cpp.h>

#include <type_traits>
#include <memory>
#include <string>
#include <vector>


namespace tl {

/**
 * get the entries of a given directory from an hdf5 file
 */
template<class t_str = std::string>
bool get_h5_entries(H5::H5File& file, const std::string& path, std::vector<t_str>& entries)
{
	// H5G_iterate_t callback, see: /usr/include/hdf5/serial/H5Dpublic.h
	H5G_iterate_t collect_elems = [](hid_t group, const char* name, void* _vec) -> herr_t
	{
		((std::vector<t_str>*)(_vec))->push_back(name);
		return 0;
	};

	if(!file.nameExists(path))
		return false;

	file.iterateElems(path, 0, collect_elems, &entries);
	return true;
}


/**
 * get an entry's attribute
 */
template<class t_str = std::string>
t_str get_h5_attr(H5::H5File& file, const t_str& path, const t_str& attr_name, bool only_group = false)
{
	if(!file.nameExists(path))
		return "";

	H5::H5Object *obj = nullptr;
	H5::DataSet set;
	H5::Group grp;

	// get entry's type
	H5G_stat_t stat;
	file.getObjinfo(path, true, stat);
	if(stat.type == H5G_GROUP)
	{
		grp = file.openGroup(path);
		obj = &grp;
	}
	else if(!only_group && stat.type == H5G_DATASET)
	{
		set = file.openDataSet(path);
		obj = &set;
	}

	if(!obj || !obj->attrExists(attr_name))
		return "";

	// get entry's attribute
	t_str attr_val;
	try
	{
		H5::Attribute attr = obj->openAttribute(attr_name);

		//hsize_t len = attr.getStorageSize();
		attr.read(attr.getStrType(), attr_val);
		attr.close();
	}
	catch(const H5::AttributeIException&)
	{}

	obj->close();
	return attr_val;
}


/**
 * get a scalar value from an hdf5 file
 */
template<class T>
bool get_h5_scalar(H5::H5File& file, const std::string& path, T& val)
{
	if(!file.nameExists(path))
		return false;

	H5::DataSet dset = file.openDataSet(path);
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank = dspace.getSimpleExtentNdims();
	if(rank != 1)
		return false;

	hsize_t dim = 0;
	dspace.getSimpleExtentDims(&dim, 0);
	if(dim != 1)
		return false;

	std::unique_ptr<H5::PredType> ty;
	if(std::is_integral<T>::value)
		ty.reset(new H5::PredType(H5::PredType::NATIVE_INT));
	else if(std::is_floating_point<T>::value)
		ty.reset(new H5::PredType(H5::PredType::NATIVE_DOUBLE));
	else
		return false;

	dset.read(&val, *ty);

	dspace.close();
	dset.close();
	return true;
}


/**
 * get a string value from an hdf5 file
 */
template<class t_str = std::string>
bool get_h5_string(H5::H5File& file, const std::string& path, t_str& val)
{
	if(!file.nameExists(path))
		return false;

	H5::DataSet dset = file.openDataSet(path);
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank = dspace.getSimpleExtentNdims();
	if(rank != 1)
		return false;

	hsize_t dim = 0;
	dspace.getSimpleExtentDims(&dim, 0);
	if(dim != 1)
		return false;

	dset.read(val, dset.getStrType());

	dspace.close();
	dset.close();
	return true;
}


/**
 * get a vector of value from an hdf5 file
 */
template<class T>
bool get_h5_vector(H5::H5File& file, const std::string& path, std::vector<T>& vals)
{
	if(!file.nameExists(path))
		return false;

	H5::DataSet dset = file.openDataSet(path);
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank = dspace.getSimpleExtentNdims();
	if(rank != 1)
		return false;

	hsize_t dim = 0;
	dspace.getSimpleExtentDims(&dim, 0);

	std::unique_ptr<H5::PredType> ty;
	if(std::is_integral<T>::value)
		ty.reset(new H5::PredType(H5::PredType::NATIVE_INT));
	else if(std::is_floating_point<T>::value)
		ty.reset(new H5::PredType(H5::PredType::NATIVE_DOUBLE));
	else
		return false;

	vals.resize(dim);
	dset.read(vals.data(), *ty);

	dspace.close();
	dset.close();
	return true;
}


/**
 * get a vector of string value from an hdf5 file
 */
template<class t_str = std::string>
bool get_h5_string_vector(H5::H5File& file, const std::string& path, std::vector<t_str>& vals)
{
	if(!file.nameExists(path))
		return false;

	H5::DataSet dset = file.openDataSet(path);
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank = dspace.getSimpleExtentNdims();
	if(rank != 1)
		return false;

	hsize_t dim = 0;
	dspace.getSimpleExtentDims(&dim, 0);
	vals.resize(dim);

	const char **chs = new const char*[dim];
	dset.read(chs, H5::StrType(H5::PredType::C_S1, H5T_VARIABLE));
	for(hsize_t i=0; i<dim; ++i)
		vals[i] = chs[i];
	for(hsize_t i=0; i<dim; ++i)
		delete[] chs[i];
	delete[] chs;

	dspace.close();
	dset.close();
	return true;
}


/**
 * get a matrix of value from an hdf5 file
 */
template<class T>
bool get_h5_matrix(H5::H5File& file, const std::string& path, std::vector<std::vector<T>>& vals)
{
	if(!file.nameExists(path))
		return false;

	H5::DataSet dset = file.openDataSet(path);
	H5::DataSpace dspace = dset.getSpace();

	hsize_t rank = dspace.getSimpleExtentNdims();
	if(rank != 2)
		return false;

	hsize_t dims[2] = { 0, 0 };
	dspace.getSimpleExtentDims(dims, 0);

	std::unique_ptr<H5::PredType> ty;
	if(std::is_integral<T>::value)
		ty.reset(new H5::PredType(H5::PredType::NATIVE_INT));
	else if(std::is_floating_point<T>::value)
		ty.reset(new H5::PredType(H5::PredType::NATIVE_DOUBLE));
	else
		return false;

	T* buf = new T[dims[0]*dims[1]];
	dset.read(buf, *ty);

	vals.resize(dims[0]);
	for(hsize_t row=0; row<dims[0]; ++row)
	{
		vals[row].resize(dims[1]);
		for(hsize_t col=0; col<dims[1]; ++col)
			vals[row][col] = buf[row*dims[1] + col];
	}

	delete[] buf;

	dspace.close();
	dset.close();
	return true;
}

}

#endif
