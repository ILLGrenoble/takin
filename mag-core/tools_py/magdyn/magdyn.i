/**
 * swig interface for magdyn
 * @author Tobias Weber <tweber@ill.fr>
 * @date 12-oct-2023
 * @license see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

%module magdyn
%{
	#include "../../tlibs2/libs/magdyn.h"
%}


%include "std_vector.i"
%include "std_list.i"
%include "std_array.i"
%include "std_string.i"
%include "std_complex.i"
%include "std_shared_ptr.i"


%template(CplxD) std::complex<double>;

%template(VecD) std::vector<double>;
%template(VecCplx) std::vector<std::complex<double>>;

%template(LstStr) std::list<std::string>;

%template(ArrD3) std::array<double, 3>;
%template(ArrD4) std::array<double, 4>;
%template(ArrCplx3) std::array<std::complex<double>, 3>;
%template(ArrCplx4) std::array<std::complex<double>, 4>;
%template(ArrStr3) std::array<std::string, 3>;
%template(ArrStr33) std::array<std::array<std::string, 3>, 3>;


// ----------------------------------------------------------------------------
// matrix and vector containers
// ----------------------------------------------------------------------------
//%include "maths.h"

// TODO: add these as soon as swig supports concepts

//%template(VectorD) tl2::vec<double>;
//%template(VectorCplx) tl2::vec<std::complex<double>>;
//%template(MatrixD) tl2::mat<double>;
//%template(MatrixCplx) tl2::mat<std::complex<double>>;
// ----------------------------------------------------------------------------


%include "../../tlibs2/libs/magdyn.h"
%include "../../tlibs2/libs/magdyn/magdyn.h"
%include "../../tlibs2/libs/magdyn/structs.h"
%include "../../tlibs2/libs/magdyn/helpers.h"
%include "../../tlibs2/libs/magdyn/getters.h"
%include "../../tlibs2/libs/magdyn/generators.h"
%include "../../tlibs2/libs/magdyn/file.h"
%include "../../tlibs2/libs/magdyn/configuration.h"
%include "../../tlibs2/libs/magdyn/groundstate.h"
%include "../../tlibs2/libs/magdyn/precalc.h"
%include "../../tlibs2/libs/magdyn/hamilton.h"
%include "../../tlibs2/libs/magdyn/correlation.h"
%include "../../tlibs2/libs/magdyn/polarisation.h"
%include "../../tlibs2/libs/magdyn/dispersion.h"


// ----------------------------------------------------------------------------
// input- and output structs (and vectors of them)
// ----------------------------------------------------------------------------
%template(MagneticSite) tl2_mag::t_MagneticSite<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::vec<double>,
	std::size_t,
	double>;
%template(VecMagneticSite) std::vector<
	tl2_mag::t_MagneticSite<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::vec<double>,
		std::size_t,
		double>
	>;
%template(VecMagneticSitePtr) std::vector<
	const tl2_mag::t_MagneticSite<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::vec<double>,
		std::size_t,
		double>*
	>;

%template(ExchangeTerm) tl2_mag::t_ExchangeTerm<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::vec<double>,
	std::size_t,
	std::complex<double>,
	double>;
%template(VecExchangeTerm) std::vector<
	tl2_mag::t_ExchangeTerm<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		tl2::vec<double>,
		std::size_t,
		std::complex<double>,
		double>
	>;

%template(ExternalField) tl2_mag::t_ExternalField<
	tl2::vec<double>,
	double>;

%template(EnergyAndWeight) tl2_mag::t_EnergyAndWeight<
	tl2::mat<std::complex<double>>,
	double, std::complex<double>>;
%template(VecEnergyAndWeight) std::vector<
	tl2_mag::t_EnergyAndWeight<
		tl2::mat<std::complex<double>>,
		double, std::complex<double>>
	>;

%template(SofQE) tl2_mag::t_SofQE<
	tl2::mat<std::complex<double>>,
	double, std::complex<double>>;
%template(VecSofQE) std::vector<
	tl2_mag::t_SofQE<
		tl2::mat<std::complex<double>>,
		double, std::complex<double>>
	>;

%template(Variable) tl2_mag::t_Variable<
	std::complex<double>>;
%template(VecVariable) std::vector<
	tl2_mag::t_Variable<
		std::complex<double>>
	>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
/**
 * main magdyn class
 */
%template(MagDyn) tl2_mag::MagDyn<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::mat<double>,
	tl2::vec<double>,
	std::complex<double>,
	double,
	std::size_t>;
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// additional functions
// ----------------------------------------------------------------------------
%inline
%{
	#include "../../libs/loadcif.h"


	/**
	 * types
	 */
	using t_str = std::string;
	using t_real = double;
	using t_mat = tl2::mat<std::complex<double>>;
	using t_vec = tl2::vec<std::complex<double>>;
	using t_mat_real = tl2::mat<double>;
	using t_vec_real = tl2::vec<double>;
	using t_cplx = std::complex<double>;

	using t_MagDyn = tl2_mag::MagDyn<
		t_mat, t_vec, t_mat_real, t_vec_real,
		t_cplx, t_real, std::size_t>;


	/**
	 * adds a variable
	 */
	void add_variable(t_MagDyn& magdyn,
		const t_str& name,
		const t_cplx& val)
	{
		typename t_MagDyn::Variable var{};

		var.name = name;
		var.value = val;

		// set existing or add new variable
		magdyn.SetVariable(std::move(var));
	}


	/**
	 * sets an external magnetic field
	 */
	void set_field(t_MagDyn& magdyn,
		t_real Bx, t_real By, t_real Bz, t_real Bmag = 1.,
		bool align_spins = true)
	{
		typename t_MagDyn::ExternalField field{};

		field.dir = tl2::create<t_vec_real>({ Bx, By, Bz });
		field.mag = Bmag;
		field.align_spins = align_spins;

		magdyn.SetExternalField(field);
		magdyn.CalcExternalField();
	}


	/**
	 * sets the temperature
	 */
	void set_temperature(t_MagDyn& magdyn, t_real T)
	{
		magdyn.SetTemperature(T);
	}


	/**
	 * sets up an ordering wave vector
	 */
	void set_ordering(t_MagDyn& magdyn, t_real x, t_real y, t_real z)
	{
		t_vec_real vec = tl2::create<t_vec_real>({ x, y, z });
		magdyn.SetOrderingWavevector(vec);
	}



	/**
	 * sets up the rotation axis for the ordering wave vector
	 */
	void set_rotation(t_MagDyn& magdyn, t_real x, t_real y, t_real z)
	{
		t_vec_real vec = tl2::create<t_vec_real>({ x, y, z });
		magdyn.SetRotationAxis(vec);
	}



	/**
	 * adds a magnetic site
	 * (using strings with expressions)
	 */
	void add_site(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& x, const t_str& y, const t_str& z,
		const t_str& sx = "0", const t_str& sy = "0", const t_str& sz = "1",
		const t_str& S = "1")
	{
		typename t_MagDyn::MagneticSite site{};

		site.name = name;

		site.pos[0] = x;
		site.pos[1] = y;
		site.pos[2] = z;

		site.spin_dir[0] = sx;
		site.spin_dir[1] = sy;
		site.spin_dir[2] = sz;
		site.spin_mag = S;

		magdyn.CalcMagneticSite(site);
		magdyn.AddMagneticSite(std::move(site));
	}


	/**
	 * adds a magnetic site
	 * (using scalar types)
	 * TODO: remove this function and use std::variant<std::string, t_real> as soon as swig supports this
	 */
	void add_site(t_MagDyn& magdyn,
		const t_str& name,
		t_real x, t_real y, t_real z,
		t_real sx = 0., t_real sy = 0., t_real sz = 1.,
		t_real S = 1.)
	{
		typename t_MagDyn::MagneticSite site{};

		site.name = name;

		site.pos[0] = tl2::var_to_str(x);
		site.pos[1] = tl2::var_to_str(y);
		site.pos[2] = tl2::var_to_str(z);

		site.spin_dir[0] = tl2::var_to_str(sx);
		site.spin_dir[1] = tl2::var_to_str(sy);
		site.spin_dir[2] = tl2::var_to_str(sz);
		site.spin_mag = tl2::var_to_str(S);

		magdyn.CalcMagneticSite(site);
		magdyn.AddMagneticSite(std::move(site));
	}


	/**
	 * adds a coupling term between two magnetic sites
	 * (using strings with expressions)
	 */
	void add_coupling(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& site1, const t_str& site2,
		const t_str& dx, const t_str& dy, const t_str& dz,
		const t_str& J,
		const t_str& dmix = "", const t_str& dmiy = "", const t_str& dmiz = "",
		const t_str& Jxx = "", const t_str& Jxy = "", const t_str& Jxz = "",
		const t_str& Jyx = "", const t_str& Jyy = "", const t_str& Jyz = "",
		const t_str& Jzx = "", const t_str& Jzy = "", const t_str& Jzz = "")
	{
		typename t_MagDyn::ExchangeTerm coupling{};

		coupling.name = name;

		coupling.site1 = site1;
		coupling.site2 = site2;

		coupling.dist[0] = dx;
		coupling.dist[1] = dy;
		coupling.dist[2] = dz;

		coupling.J = J;

		coupling.dmi[0] = dmix;
		coupling.dmi[1] = dmiy;
		coupling.dmi[2] = dmiz;

		coupling.Jgen[0][0] = Jxx;
		coupling.Jgen[0][1] = Jxy;
		coupling.Jgen[0][2] = Jxz;
		coupling.Jgen[1][0] = Jyx;
		coupling.Jgen[1][1] = Jyy;
		coupling.Jgen[1][2] = Jyz;
		coupling.Jgen[2][0] = Jzx;
		coupling.Jgen[2][1] = Jzy;
		coupling.Jgen[2][2] = Jzz;

		magdyn.CalcExchangeTerm(coupling);
		magdyn.AddExchangeTerm(std::move(coupling));
	}


	/**
	 * adds a coupling term between two magnetic sites
	 * (using scalar types)
	 * TODO: remove this function and use std::variant<std::string, t_real> as soon as swig supports this
	 */
	void add_coupling(t_MagDyn& magdyn,
		const t_str& name,
		const t_str& site1, const t_str& site2,
		t_real dx, t_real dy, t_real dz,
		t_real J,
		t_real dmix = 0., t_real dmiy = 0., t_real dmiz = 0.,
		t_real Jxx = 0., t_real Jxy = 0., t_real Jxz = 0.,
		t_real Jyx = 0., t_real Jyy = 0., t_real Jyz = 0.,
		t_real Jzx = 0., t_real Jzy = 0., t_real Jzz = 0.)
	{
		typename t_MagDyn::ExchangeTerm coupling{};

		coupling.name = name;

		coupling.site1 = site1;
		coupling.site2 = site2;

		coupling.dist[0] = tl2::var_to_str(dx);
		coupling.dist[1] = tl2::var_to_str(dy);
		coupling.dist[2] = tl2::var_to_str(dz);

		coupling.J = tl2::var_to_str(J);

		coupling.dmi[0] = tl2::var_to_str(dmix);
		coupling.dmi[1] = tl2::var_to_str(dmiy);
		coupling.dmi[2] = tl2::var_to_str(dmiz);

		coupling.Jgen[0][0] = tl2::var_to_str(Jxx);
		coupling.Jgen[0][1] = tl2::var_to_str(Jxy);
		coupling.Jgen[0][2] = tl2::var_to_str(Jxz);
		coupling.Jgen[1][0] = tl2::var_to_str(Jyx);
		coupling.Jgen[1][1] = tl2::var_to_str(Jyy);
		coupling.Jgen[1][2] = tl2::var_to_str(Jyz);
		coupling.Jgen[2][0] = tl2::var_to_str(Jzx);
		coupling.Jgen[2][1] = tl2::var_to_str(Jzy);
		coupling.Jgen[2][2] = tl2::var_to_str(Jzz);

		magdyn.CalcExchangeTerm(coupling);
		magdyn.AddExchangeTerm(std::move(coupling));
	}


	/**
	 * calculates sites and couplings
	 */
	void calc(t_MagDyn& magdyn)
	{
		magdyn.CalcExternalField();
		magdyn.CalcMagneticSites();
		magdyn.CalcExchangeTerms();
	}


	/**
	 * get a list of all space groups
	 */
	std::list<std::string> get_spacegroups()
	{
		std::list<std::string> sg_names;

		for(const auto& sg : get_sgs<t_mat_real>(true, false))
			sg_names.push_back(std::get<1>(sg));

		return sg_names;
	}


	/**
	 * create symmetry equivalent site positions using a space group
	 */
	 void symmetrise_sites(t_MagDyn& magdyn, const std::string& spacegroup)
	 {
		std::vector<t_mat_real> ops = get_sg_ops<t_mat_real>(spacegroup);
		if(ops.size() == 0)
		{
			std::cerr << "Error: Space group \"" << spacegroup << "\" was not found. "
				<< "Please use get_spacegroups() to get a list."
				<< std::endl;
			return;
		}

		magdyn.SymmetriseMagneticSites(ops);
	 }


	/**
	 * create symmetry equivalent exchange couplings using a space group
	 */
	 void symmetrise_couplings(t_MagDyn& magdyn, const std::string& spacegroup)
	 {
		std::vector<t_mat_real> ops = get_sg_ops<t_mat_real>(spacegroup);
		if(ops.size() == 0)
		{
			std::cerr << "Error: Space group \"" << spacegroup << "\" was not found. "
				<< "Please use get_spacegroups() to get a list."
				<< std::endl;
			return;
		}

		magdyn.SymmetriseExchangeTerms(ops);
	 }
%}
// ----------------------------------------------------------------------------
