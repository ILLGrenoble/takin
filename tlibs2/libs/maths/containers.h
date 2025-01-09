/**
 * tlibs2 maths library -- containers and adapters
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
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

#ifndef __TLIBS2_MATHS_CONTS_H__
#define __TLIBS2_MATHS_CONTS_H__

#include <cmath>
#include <vector>

#include "decls.h"
#include "operators.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// adapters
// ----------------------------------------------------------------------------
/**
 * vector-like access adapter to a matrix
 */
template<class t_mat> requires is_basic_mat<t_mat>
class matvec_adapter
{
public:
	using value_type = typename t_mat::value_type;
	using size_type = decltype(t_mat{}.size1());

public:
	matvec_adapter(const t_mat &mat) : m_mat{mat} {}
	~matvec_adapter() = default;

	size_type size() const { return m_mat.size1() * m_mat.size2(); }

	const value_type& operator[](size_type i) const
	{
		size_type row = i/m_mat.size2();
		size_type col = i%m_mat.size2();

		return m_mat(row, col);
	}

private:
	const t_mat& m_mat;
};


/**
 * adapter for a qvector
 */
template<typename size_t, size_t N, typename T,
	template<size_t, size_t, class...> class t_mat_base>
class qvec_adapter : public t_mat_base<1, N, T>
{
public:
	// types
	using base_type = t_mat_base<1, N, T>;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qvec_adapter(const base_type& vec) : base_type{vec} {}

	static constexpr size_t size() { return N; }

	T& operator[](size_t i)
	{
		return base_type::operator()(i,0);
	}

	const T operator[](size_t i) const
	{
		return base_type::operator()(i,0);
	}
};


/**
 * adapter for a qmatrix
 */
template<typename size_t, size_t ROWS, size_t COLS, typename T,
	template<size_t, size_t, class...> class t_mat_base>
class qmat_adapter : public t_mat_base<COLS, ROWS, T>
{
public:
	// types
	using base_type = t_mat_base<COLS, ROWS, T>;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qmat_adapter(const base_type& mat) : base_type{mat} {}

	static constexpr size_t size1() { return ROWS; }
	static constexpr size_t size2() { return COLS; }
};


/**
 * adapter for a qvector
 */
template<typename size_t, size_t N, typename T, class t_vec_base>
class qvecN_adapter : public t_vec_base
{
public:
	// types
	using base_type = t_vec_base;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qvecN_adapter(const base_type& vec) : base_type{vec} {}

	static constexpr size_t size() { return N; }

	T& operator[](size_t i)
	{
		return static_cast<base_type&>(*this)[i];
	}

	const T operator[](size_t i) const
	{
		return static_cast<const base_type&>(*this)[i];
	}
};


/**
 * adapter for a qmatrix
 */
template<typename size_t, size_t ROWS, size_t COLS, typename T, class t_mat_base>
class qmatNN_adapter : public t_mat_base
{
public:
	// types
	using base_type = t_mat_base;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qmatNN_adapter(const base_type& mat) : base_type{mat} {}

	// convert from a different matrix type
	template<class t_matOther> qmatNN_adapter(const t_matOther& matOther)
		requires is_basic_mat<t_matOther>
	{
		const std::size_t minRows = std::min(
			static_cast<std::size_t>(size1()),
			static_cast<std::size_t>(matOther.size1()));
		const std::size_t minCols = std::min(
			static_cast<std::size_t>(size2()),
			static_cast<std::size_t>(matOther.size2()));

		for(std::size_t i = 0; i < minRows; ++i)
			for(std::size_t j = 0; j < minCols; ++j)
				(*this)(i, j) = static_cast<value_type>(matOther(i, j));
	}

	static constexpr size_t size1() { return ROWS; }
	static constexpr size_t size2() { return COLS; }
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// vector and matrix containers
// ----------------------------------------------------------------------------

/**
 * generic vector container
 */
template<class T = double, template<class...> class t_cont = std::vector>
requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
class vec : public t_cont<T>
{
public:
	using value_type = T;
	using container_type = t_cont<T>;

	using container_type::container_type;
	using container_type::size;
	using container_type::operator[];

	using typename container_type::iterator;
	using typename container_type::const_iterator;
	using typename container_type::size_type;
	using typename container_type::difference_type;
	using typename container_type::allocator_type;

	~vec() = default;

	vec() : container_type{}
	{}

	vec(const vec<T, t_cont>& other)
		: container_type{other}
	{}

	vec(vec<T, t_cont>&& other)
		: container_type{std::forward<vec<T, t_cont>&&>(other)}
	{}

	template<class T_other, template<class...> class t_cont_other>
	vec(const vec<T_other, t_cont_other>& other)
	{
		this->operator=<T_other, t_cont_other>(other);
	}

	vec<T, t_cont>& operator=(const vec<T, t_cont>& other)
	{
		*static_cast<container_type*>(this) = other;
		return *this;
	}

	const vec<T, t_cont>& operator=(const vec<T, t_cont>& other) const
	{
		*static_cast<container_type*>(this) = other;
		return *this;
	}

	template<class T_other, template<class...> class t_cont_other>
	vec<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other)
	{
		*this = convert<vec<T, t_cont>, vec<T_other, t_cont_other>>(other);
		return *this;
	}

	template<class T_other, template<class...> class t_cont_other>
	const vec<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other) const
	{
		*this = convert<vec<T, t_cont>, vec<T_other, t_cont_other>>(other);
		return *this;
	}

	vec(std::size_t SIZE, const T* arr = nullptr) : container_type(SIZE)
	{
		if(arr)
			from_array(arr);
	}


	const value_type& operator()(std::size_t i) const
	{
		return this->operator[](i);
	}

	value_type& operator()(std::size_t i)
	{
		return this->operator[](i);
	}


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(std::size_t i = 0; i < size(); ++i)
			this->operator[](i) = arr[i];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(std::size_t i = 0; i < size(); ++i)
			arr[i] = this->operator[](i);
	}

	friend vec operator+(const vec& vec1, const vec& vec2) { return tl2_ops::operator+(vec1, vec2); }
	friend vec operator-(const vec& vec1, const vec& vec2) { return tl2_ops::operator-(vec1, vec2); }
	friend const vec& operator+(const vec& vec1) { return tl2_ops::operator+(vec1); }
	friend vec operator-(const vec& vec1) { return tl2_ops::operator-(vec1); }

	friend vec operator*(value_type d, const vec& vec1) { return tl2_ops::operator*(d, vec1); }
	friend vec operator*(const vec& vec1, value_type d) { return tl2_ops::operator*(vec1, d); }
	friend vec operator/(const vec& vec1, value_type d) { return tl2_ops::operator/(vec1, d); }

	vec& operator*=(const vec& vec2) { return tl2_ops::operator*=(*this, vec2); }
	vec& operator+=(const vec& vec2) { return tl2_ops::operator+=(*this, vec2); }
	vec& operator-=(const vec& vec2) { return tl2_ops::operator-=(*this, vec2); }
	vec& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	vec& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }
};



/**
 * generic matrix container
 */
template<class T=double, template<class...> class t_cont = std::vector>
requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
class mat
{
public:
	using value_type = T;
	using container_type = t_cont<T>;

	mat() = default;
	~mat() = default;

	mat(std::size_t ROWS, std::size_t COLS, const T* arr = nullptr)
		: m_data(ROWS*COLS), m_rowsize{ROWS}, m_colsize{COLS}
	{
		if(arr)
			from_array(arr);
	}


	template<class T_other, template<class...> class t_cont_other>
	mat(const mat<T_other, t_cont_other>& other)
	{
		this->operator=<T_other, t_cont_other>(other);
	}

	template<class T_other, template<class...> class t_cont_other>
	mat<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other)
	{
		*this = convert<mat<T, t_cont>, mat<T_other, t_cont_other>>(other);

		this->m_rowsize = other.m_rowsize;
		this->m_colsize = other.m_colsize;

		return *this;
	}

	template<class T_other, template<class...> class t_cont_other>
	const mat<T, t_cont>& operator=(const vec<T_other, t_cont_other>& other) const
	{
		*this = convert<mat<T, t_cont>, mat<T_other, t_cont_other>>(other);

		this->m_rowsize = other.m_rowsize;
		this->m_colsize = other.m_colsize;

		return *this;
	}


	std::size_t size1() const { return m_rowsize; }
	std::size_t size2() const { return m_colsize; }


	// element access
	const T& operator()(std::size_t row, std::size_t col) const
	{
		return m_data[row*m_colsize + col];
	}

	T& operator()(std::size_t row, std::size_t col)
	{
		return m_data[row*m_colsize + col];
	}


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(std::size_t i = 0; i < m_rowsize; ++i)
			for(std::size_t j = 0; j < m_colsize; ++j)
				this->operator()(i, j) = arr[i*m_colsize + j];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(std::size_t i = 0; i < m_rowsize; ++i)
			for(std::size_t j = 0; j < m_colsize; ++j)
				arr[i*m_colsize + j] = this->operator()(i,j);
	}

	friend mat operator+(const mat& mat1, const mat& mat2) { return tl2_ops::operator+(mat1, mat2); }
	friend mat operator-(const mat& mat1, const mat& mat2) { return tl2_ops::operator-(mat1, mat2); }
	friend const mat& operator+(const mat& mat1) { return tl2_ops::operator+(mat1); }
	friend mat operator-(const mat& mat1) { return tl2_ops::operator-(mat1); }

	friend mat operator*(const mat& mat1, const mat& mat2) { return tl2_ops::operator*(mat1, mat2); }
	friend mat operator*(const mat& mat1, value_type d) { return tl2_ops::operator*(mat1, d); }
	friend mat operator*(value_type d, const mat& mat1) { return tl2_ops::operator*(d, mat1); }
	friend mat operator/(const mat& mat1, value_type d) { return tl2_ops::operator/(mat1, d); }

	template<class t_vec> requires is_basic_vec<t_cont<T>> && is_dyn_vec<t_cont<T>>
	friend t_vec operator*(const mat& mat1, const t_vec& vec2) { return tl2_ops::operator*(mat1, vec2); }

	mat& operator*=(const mat& mat2) { return tl2_ops::operator*=(*this, mat2); }
	mat& operator+=(const mat& mat2) { return tl2_ops::operator+=(*this, mat2); }
	mat& operator-=(const mat& mat2) { return tl2_ops::operator-=(*this, mat2); }
	mat& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	mat& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }

private:
	container_type m_data{};
	std::size_t m_rowsize{0};
	std::size_t m_colsize{0};
};

// ----------------------------------------------------------------------------

}

#endif
