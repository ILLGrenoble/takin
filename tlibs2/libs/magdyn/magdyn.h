/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc The magdyn library implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#ifndef __TLIBS2_MAGDYN_DECL_H__
#define __TLIBS2_MAGDYN_DECL_H__

#include <vector>
#include <array>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <string>
#include <string_view>

#include <boost/container_hash/hash.hpp>
#include <boost/property_tree/ptree.hpp>

#include "../maths.h"
#include "../expr.h"

#include "structs.h"
#include "helpers.h"



#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
#define MAGDYN_TEMPL                                                \
	template<                                                   \
		class t_mat, class t_vec,                           \
		class t_mat_real, class t_vec_real,                 \
		class t_cplx, class t_real, class t_size>           \
	requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> &&        \
		tl2::is_mat<t_mat_real> && tl2::is_vec<t_vec_real>
#else  // SWIG
#define MAGDYN_TEMPL                                          \
	template<                                             \
		class t_mat, class t_vec,                     \
		class t_mat_real, class t_vec_real,           \
		class t_cplx, class t_real, class t_size>
#endif  // SWIG

#define MAGDYN_INST                                           \
	tl2_mag::MagDyn<t_mat, t_vec, t_mat_real, t_vec_real, \
		t_cplx, t_real, t_size>

#define MAGDYN_TYPE typename MAGDYN_INST


// only print if it's not set to silent mode
#define CERR_OPT if(!m_silent) std::cerr


namespace tl2_mag {


/**
 * calculates magnon dynamics,
 * implementing the formalism given in (Toth 2015)
 */
template<
	class t_mat, class t_vec,
	class t_mat_real, class t_vec_real,
	class t_cplx = typename t_mat::value_type,
	class t_real = typename t_mat_real::value_type,
	class t_size = std::size_t>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> &&
	tl2::is_mat<t_mat_real> && tl2::is_vec<t_vec_real>
#endif
class MagDyn
{
public:
	// --------------------------------------------------------------------
	// structs and types
	// --------------------------------------------------------------------
	using MagneticSite = t_MagneticSite<t_mat, t_vec, t_vec_real, t_size, t_real>;
	using MagneticSites = std::vector<MagneticSite>;

	using ExchangeTerm = t_ExchangeTerm<t_mat, t_vec, t_vec_real, t_size, t_cplx, t_real>;
	using ExchangeTerms = std::vector<ExchangeTerm>;

	using Variable = t_Variable<t_cplx>;
	using Variables = std::vector<Variable>;

	using ExternalField = t_ExternalField<t_vec_real, t_real>;

	using EnergyAndWeight = t_EnergyAndWeight<t_mat, t_real, t_cplx>;
	using EnergiesAndWeights = std::vector<EnergyAndWeight>;

	using SofQE = t_SofQE<t_mat, t_real, t_cplx>;
	using SofQEs = std::vector<SofQE>;

	using t_indices = std::pair<t_size, t_size>;
	using t_Jmap = std::unordered_map<t_indices, t_mat, boost::hash<t_indices>>;
	// --------------------------------------------------------------------



public:
	MagDyn() = default;
	~MagDyn() = default;

	MagDyn(const MagDyn&) = default;
	MagDyn& operator=(const MagDyn&) = default;


	// --------------------------------------------------------------------
	// cleanup functions
	// --------------------------------------------------------------------
	/**
	 * clear all
	 */
	void Clear();

	/**
	 * clear all parser variables
	 */
	void ClearVariables();
	void ClearMagneticSites();
	void ClearExchangeTerms();
	void ClearExternalField();
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// getter
	// --------------------------------------------------------------------
	const Variables& GetVariables() const;
	const MagneticSites& GetMagneticSites() const;
	MagneticSites& GetMagneticSites();
	t_size GetMagneticSitesCount() const;
	const ExchangeTerms& GetExchangeTerms() const;
	ExchangeTerms& GetExchangeTerms();
	t_size GetExchangeTermsCount() const;

	const ExternalField& GetExternalField() const;
	const t_vec_real& GetRotationAxis() const;
	const t_vec_real& GetOrderingWavevector() const;

	t_real GetTemperature() const;
	t_real GetBoseCutoffEnergy() const;

	const std::string& GetMagneticFormFactor() const;

	const t_mat_real& GetCrystalATrafo() const;
	const t_mat_real& GetCrystalBTrafo() const;
	const t_mat_real& GetCrystalUBTrafo() const;

	const MagneticSite& GetMagneticSite(t_size idx) const;
	const ExchangeTerm& GetExchangeTerm(t_size idx) const;

	bool IsIncommensurate() const;
	bool GetSilent() const;

	/**
	 * get number of magnetic sites with the given name (to check if the name is unique)
	 */
	std::vector<const MagneticSite*> FindMagneticSites(const std::string& name) const;

	/**
	 * get magnetic site with the given name
	 */
	const MagneticSite* FindMagneticSite(const std::string& name) const;

	/**
	 * get the index of a magnetic site from its name
	 */
	t_size GetMagneticSiteIndex(const std::string& name) const;

	/**
	 * get the index of an exchange term from its name
	 */
	t_size GetExchangeTermIndex(const std::string& name) const;

	std::vector<t_vec_real> GetMagneticSitePositions(bool homogeneous = false) const;

	t_vec_real GetCrystalLattice() const;
	const t_vec_real* GetScatteringPlane() const;

	/**
	 * get the needed supercell ranges from the exchange terms
	 */
	std::tuple<t_vec_real, t_vec_real> GetSupercellMinMax() const;
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// setter
	// --------------------------------------------------------------------
	void SetEpsilon(t_real eps);
	void SetPrecision(int prec);

	void SetTemperature(t_real T);
	void SetBoseCutoffEnergy(t_real E);
	void SetUniteDegenerateEnergies(bool b);
	void SetForceIncommensurate(bool b);
	void SetPerformChecks(bool b);
	void SetSilent(bool b);

	void SetPhaseSign(t_real sign);
	void SetCholeskyMaxTries(t_size max_tries);
	void SetCholeskyInc(t_real delta);

	void SetMagneticFormFactor(const std::string& ffact);

	void SetExternalField(const ExternalField& field);
	void RotateExternalField(const t_vec_real& axis, t_real angle);
	void RotateExternalField(t_real x, t_real y, t_real z, t_real angle);

	/**
	 * set the ordering wave vector (e.g., the helix pitch) for incommensurate structures
	 */
	void SetOrderingWavevector(const t_vec_real& ordering);

	/**
	 * set the rotation axis for the ordering wave vector
	 */
	void SetRotationAxis(const t_vec_real& axis);

	void SetCalcHamiltonian(bool H, bool Hp, bool Hm);

	void AddVariable(Variable&& var);
	void SetVariable(Variable&& var);

	void AddMagneticSite(MagneticSite&& site);
	void AddExchangeTerm(ExchangeTerm&& term);

	/**
	 * calculate the B matrix from the crystal lattice
	 */
	void SetCrystalLattice(t_real a, t_real b, t_real c,
		t_real alpha, t_real beta, t_real gamma);

	/**
	 * calculate the UB matrix from the scattering plane and the crystal lattice
	 * note: SetCrystalLattice() has to be called before this function
	 */
	void SetScatteringPlane(t_real ah, t_real ak, t_real al,
		t_real bh, t_real bk, t_real bl);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	/**
	 * get an expression parser object with registered variables
	 */
	tl2::ExprParser<t_cplx> GetExprParser() const;
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// sanity checks
	// --------------------------------------------------------------------
	/**
	 * check if the site index is valid
	 */
	bool CheckMagneticSite(t_size idx, bool print_error = true) const;

	/**
	 * check if the term index is valid
	 */
	bool CheckExchangeTerm(t_size idx, bool print_error = true) const;

	/**
	 * check if imaginary weights remain
	 */
	bool CheckImagWeights(const t_vec_real& Q_rlu, const EnergiesAndWeights& Es_and_S) const;
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// symmetrisation and generation functions
	// --------------------------------------------------------------------
	/**
	 * generate symmetric positions based on the given symops
	 */
	void SymmetriseMagneticSites(const std::vector<t_mat_real>& symops);

	/**
	 * generate symmetric exchange terms based on the given symops
	 */
	void SymmetriseExchangeTerms(const std::vector<t_mat_real>& symops);

	/**
	 * generate possible couplings up to a certain distance
	 */
	void GeneratePossibleExchangeTerms(
		t_real dist_max, t_size _sc_max,
		std::optional<t_size> couplings_max);

	/**
	 * extend the magnetic structure
	 */
	void ExtendStructure(t_size x_size, t_size y_size, t_size z_size,
		bool remove_duplicates = true, bool flip_spin = true);

	/**
	 * modify exchange term whose sites point to sc positions that are also available in the uc
	 */
	void FixExchangeTerms(t_size x_size = 1, t_size y_size = 1, t_size z_size = 1);

	/**
	 * remove literal duplicate sites (not symmetry-equivalent ones)
	 */
	void RemoveDuplicateMagneticSites();

	/**
	 * remove literal duplicate couplings (not symmetry-equivalent ones)
	 */
	void RemoveDuplicateExchangeTerms();

	/**
	 * are two sites equivalent with respect to the given symmetry operators?
	 */
	bool IsSymmetryEquivalent(const MagneticSite& site1, const MagneticSite& site2,
		const std::vector<t_mat_real>& symops) const;

	/**
	 * are two couplings equivalent with respect to the given symmetry operators?
	 */
	bool IsSymmetryEquivalent(const ExchangeTerm& term1, const ExchangeTerm& term2,
		const std::vector<t_mat_real>& symops) const;

	/**
	 * assign symmetry group indices to sites and couplings
	 */
	void CalcSymmetryIndices(const std::vector<t_mat_real>& symops);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// calculation functions
	// --------------------------------------------------------------------
	/**
	 * calculate the rotation matrix for the external field
	 */
	void CalcExternalField();

	/**
	 * calculate the spin rotation trafo for the magnetic sites
	 * and parse any given expressions
	 */
	void CalcMagneticSite(MagneticSite& site);

	/**
	 * calculate the spin rotation trafo for the magnetic sites
	 * and parse any given expressions
	 */
	void CalcMagneticSites();

	/**
	 * parse the exchange term expressions and calculate all properties
	 */
	void CalcExchangeTerm(ExchangeTerm& term);

	/**
	 * parse all exchange term expressions and calculate all properties
	 */
	void CalcExchangeTerms();

	/**
	 * calculate the real-space interaction matrix J of
	 * equations (10) - (13) from (Toth 2015)
	 */
	t_mat CalcRealJ(const ExchangeTerm& term) const;

	/**
	 * calculate the reciprocal interaction matrices J(Q) and J(-Q) of
	 * equations (12) and (14) from (Toth 2015)
	 */
	std::tuple<t_Jmap, t_Jmap> CalcReciprocalJs(const t_vec_real& Qvec) const;

	/**
	 * get the hamiltonian at the given momentum
	 * @note implements the formalism given by (Toth 2015)
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_mat CalcHamiltonian(const t_vec_real& Qvec) const;

	/**
	 * get the energies from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	EnergiesAndWeights CalcEnergiesFromHamiltonian(
		t_mat _H, const t_vec_real& Qvec,
		bool only_energies = false) const;

	/**
	 * get the dynamical structure factor from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	bool CalcCorrelationsFromHamiltonian(EnergiesAndWeights& energies_and_correlations,
		const t_mat& H_mat, const t_mat& chol_mat, const t_mat& g_sign,
		const t_vec_real& Qvec, const std::vector<t_vec>& evecs) const;

	/**
	 * applies projectors, form and weight factors to get neutron intensities
	 * @note implements the formalism given by (Toth 2015)
	 */
	void CalcIntensities(const t_vec_real& Q_rlu, EnergiesAndWeights& E_and_S) const;

	/**
	 * calculates the polarisation matrix
	 */
	void CalcPolarisation(const t_vec_real& Q_rlu, EnergyAndWeight& E_and_S) const;

	/**
	 * unite degenerate energies and their corresponding eigenstates
	 */
	EnergiesAndWeights UniteEnergies(const EnergiesAndWeights&
		energies_and_correlations) const;

	/**
	 * get the energies and the spin-correlation at the given momentum
	 * (also calculates incommensurate contributions and applies weight factors)
	 * @note implements the formalism given by (Toth 2015)
	 */
	EnergiesAndWeights CalcEnergies(const t_vec_real& Q_rlu,
		bool only_energies = false) const;

	EnergiesAndWeights CalcEnergies(t_real h, t_real k, t_real l,
		bool only_energies = false) const;

	/**
	 * generates the dispersion plot along the given q path
	 */
	SofQEs CalcDispersion(t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end, t_size num_Qs = 128,
		t_size num_threads = 4,
		std::function<bool(int, int)> *progress_fkt = nullptr) const;

	/**
	 * get the energy minimum
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_real CalcMinimumEnergy() const;

	/**
	 * get the ground-state energy
	 * @note zero-operator term in expansion of equation (20) in (Toth 2015)
	 */
	t_real CalcGroundStateEnergy() const;

	/**
	 * minimise energy to found ground state
	 */
	bool CalcGroundState(const std::unordered_set<std::string>* fixed_params = nullptr,
		bool verbose = false, const bool *stop_request = nullptr);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// loading and saving
	// --------------------------------------------------------------------
	/**
	 * generates the dispersion plot along the given Q path
	 */
	bool SaveDispersion(const std::string& filename,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_Qs = 128, t_size num_threads = 4,
		bool as_py = false,
		std::function<bool(int, int)> *progress_fkt = nullptr) const;

	/**
	 * generates the dispersion plot along the given Q path
	 */
	bool SaveDispersion(std::ostream& ostr,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_Qs = 128, t_size num_threads = 4,
		bool as_py = false,
		std::function<bool(int, int)> *progress_fkt = nullptr,
		bool write_header = true) const;

	/**
	 * generates the dispersion plot along multiple Q paths
	 */
	bool SaveMultiDispersion(const std::string& filename,
		const std::vector<t_vec_real>& Qs,
		t_size num_Qs = 128, t_size num_threads = 4,
		bool as_py = false,
		std::function<bool(int, int)> *progress_fkt = nullptr,
		const std::vector<std::string>* Q_names = nullptr) const;

	/**
	 * generates the dispersion plot along multiple Q paths
	 */
	bool SaveMultiDispersion(std::ostream& ostr,
		const std::vector<t_vec_real>& Qs,
		t_size num_Qs = 128, t_size num_threads = 4,
		bool as_py = false,
		std::function<bool(int, int)> *progress_fkt = nullptr,
		const std::vector<std::string>* Q_names = nullptr) const;

	/**
	 * load a configuration from a file
	 */
	bool Load(const std::string& filename);

	/**
	 * save the configuration to a file
	 */
	bool Save(const std::string& filename) const;

	/**
	 * load a configuration from a property tree
	 */
	bool Load(const boost::property_tree::ptree& node);

	/**
	 * save the configuration to a property tree
	 */
	bool Save(boost::property_tree::ptree& node) const;
	// --------------------------------------------------------------------


protected:
	/**
	 * converts the rotation matrix rotating the local spins to ferromagnetic
	 * [001] directions into the vectors comprised of the matrix columns
	 * @see equation (9) and (51) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> rot_to_trafo(const t_mat& R);

	/**
	 * rotate local spin to ferromagnetic [001] direction
	 * @see equations (7) and (9) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> spin_to_trafo(const t_vec_real& spin_dir);


private:
	// magnetic sites
	MagneticSites m_sites{};

	// magnetic couplings
	ExchangeTerms m_exchange_terms{};

	// open variables in expressions
	Variables m_variables{};

	// external field
	ExternalField m_field{};
	// matrix to rotate field into the [001] direction
	t_mat m_rot_field{ tl2::unit<t_mat>(3) };

	// ordering wave vector for incommensurate structures
	t_vec_real m_ordering{ tl2::zero<t_vec_real>(3) };

	// helix rotation axis for incommensurate structures
	t_vec_real m_rotaxis{ tl2::create<t_vec_real>({ 1., 0., 0. }) };

	// calculate the hamiltonian for Q, Q+ordering, and Q-ordering
	bool m_calc_H{ true };
	bool m_calc_Hp{ true };
	bool m_calc_Hm{ true };

	// direction to rotation spins into, usually [001]
	t_vec_real m_zdir{ tl2::create<t_vec_real>({ 0., 0., 1. }) };

	// temperature (-1: disable bose factor)
	t_real m_temperature{ -1. };

	// bose cutoff energy to avoid infinities
	t_real m_bose_cutoff{ 0.025 };

	// formula for the magnetic form factor
	std::string m_magffact_formula{};
	tl2::ExprParser<t_cplx> m_magffact{};

	// crystal lattice
	t_real m_xtallattice[3]{ 5., 5., 5. };
	t_real m_xtalangles[3]
	{
		t_real(0.5) * tl2::pi<t_real>,
		t_real(0.5) * tl2::pi<t_real>,
		t_real(0.5) * tl2::pi<t_real>
	};
	t_mat_real m_xtalA{ tl2::unit<t_mat_real>(3) };
	t_mat_real m_xtalB{ tl2::unit<t_mat_real>(3) };
	t_mat_real m_xtalUB{ tl2::unit<t_mat_real>(3) };
	t_mat_real m_xtalUBinv{ tl2::unit<t_mat_real>(3) };

	//scattering plane
	t_vec_real m_scatteringplane[3]
	{
		tl2::create<t_vec_real>({ 1., 0., 0. }),  // in-plane, x
		tl2::create<t_vec_real>({ 0., 1., 0. }),  // in-plane, y
		tl2::create<t_vec_real>({ 0., 0., 1. }),  // out-of-plane, z
	};

	// settings
	bool m_is_incommensurate{ false };
	bool m_force_incommensurate{ false };
	bool m_unite_degenerate_energies{ true };
	bool m_perform_checks{ true };
	bool m_silent { false };

	// settings for cholesky decomposition
	t_size m_tries_chol{ 50 };
	t_real m_delta_chol{ 0.0025 };

	// precisions
	t_real m_eps{ 1e-6 };
	int m_prec{ 6 };

	// conventions
	t_real m_phase_sign{ -1. };

	// constants
	static constexpr const t_cplx s_imag { t_real(0), t_real(1) };
	static constexpr const t_real s_twopi { t_real(2)*tl2::pi<t_real> };
	static constexpr const std::array<std::string_view, 3> g_comp_names{{ "x", "y", "z" }};
};

}
#endif
