/**
 * tlibs2 -- magnetic dynamics -- saving of dispersions
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
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

#ifndef __TLIBS2_MAGDYN_FILE_H__
#define __TLIBS2_MAGDYN_FILE_H__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/algorithm/string/replace.hpp>

#include "../algos.h"
#include "magdyn.h"


namespace tl2_mag {

// ============================================================================
// py plotting script template
// ============================================================================
template<class t_str = std::string>
t_str g_pyscr = R"RAW(import sys
import numpy
import matplotlib.pyplot as pyplot
pyplot.rcParams.update({
	"font.sans-serif" : "DejaVu Sans",
	"font.family" : "sans-serif",
	"font.size" : 12,
})


# -----------------------------------------------------------------------------
# options
# -----------------------------------------------------------------------------
show_dividers  = False  # show vertical bars between dispersion branches
plot_file      = ""     # file to save plot to
only_pos_E     = True   # ignore magnon annihilation?

S_filter_min   = 1e-5   # cutoff minimum spectral weight
S_scale        = 10.    # spectral weight scaling factor
S_clamp_min    = 0.1    # spectral weight minimum clamping factor
S_clamp_max    = 100.   # spectral weight maximum clamping factor

# settings for dispersion branches
branch_labels  = %%LABELS%%
width_ratios   = %%RATIOS%%
branch_colours = None
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# plot the dispersion branches
# -----------------------------------------------------------------------------
def plot_disp(data):
	num_branches = len(data)

	(plt, axes) = pyplot.subplots(nrows = 1, ncols = num_branches,
		width_ratios = width_ratios, sharey = True)

	# in case there's only one sub-plot
	if type(axes) != numpy.ndarray:
		axes = [ axes ]

	for branch_idx in range(len(data)):
		branch_data = numpy.array(data[branch_idx]).transpose()

		data_h = branch_data[0]
		data_k = branch_data[1]
		data_l = branch_data[2]
		data_E = branch_data[3]
		data_S = branch_data[4]

		if only_pos_E:
			# ignore magnon annihilation
			data_h = numpy.array([ h for (h, E) in zip(data_h, data_E) if E >= 0. ])
			data_k = numpy.array([ k for (k, E) in zip(data_k, data_E) if E >= 0. ])
			data_l = numpy.array([ l for (l, E) in zip(data_l, data_E) if E >= 0. ])
			data_S = numpy.array([ S for (S, E) in zip(data_S, data_E) if E >= 0. ])
			data_E = numpy.array([ E for E in data_E if E >= 0. ])

		if S_filter_min >= 0.:
			# filter weights below cutoff
			data_h = numpy.array([ h for (h, S) in zip(data_h, data_S) if S >= S_filter_min ])
			data_k = numpy.array([ k for (k, S) in zip(data_k, data_S) if S >= S_filter_min ])
			data_l = numpy.array([ l for (l, S) in zip(data_l, data_S) if S >= S_filter_min ])
			data_E = numpy.array([ E for (E, S) in zip(data_E, data_S) if S >= S_filter_min ])
			data_S = numpy.array([ S for S in data_S if S >= S_filter_min ])

		# branch start and end point
		start_Q = ( data_h[0], data_k[0], data_l[0] )
		end_Q = ( data_h[-1], data_k[-1], data_l[-1] )

		# find scan axis
		Q_diff = [
			numpy.abs(start_Q[0] - end_Q[0]),
			numpy.abs(start_Q[1] - end_Q[1]),
			numpy.abs(start_Q[2] - end_Q[2]) ]

		plot_idx = 0
		data_x = data_h
		if Q_diff[1] > Q_diff[plot_idx]:
			plot_idx = 1
			data_x = data_k
		elif Q_diff[2] > Q_diff[plot_idx]:
			plot_idx = 2
			data_x = data_l

		# ticks and labels
		axes[branch_idx].set_xlim(data_x[0], data_x[-1])

		if branch_colours != None and len(branch_colours) != 0:
			axes[branch_idx].set_facecolor(branch_colours[branch_idx])

		if branch_labels != None and len(branch_labels) != 0:
			tick_labels = [
				branch_labels[branch_idx],
				branch_labels[branch_idx + 1] ]
		else:
			tick_labels = [
				"(%.4g %.4g %.4g)" % (start_Q[0], start_Q[1], start_Q[2]),
				"(%.4g %.4g %.4g)" % (end_Q[0], end_Q[1], end_Q[2]) ]

		if branch_idx == 0:
			axes[branch_idx].set_ylabel("E (meV)")
		else:
			axes[branch_idx].get_yaxis().set_visible(False)
			if not show_dividers:
				axes[branch_idx].spines["left"].set_visible(False)

			tick_labels[0] = ""

		if not show_dividers and branch_idx != num_branches - 1:
			axes[branch_idx].spines["right"].set_visible(False)

		axes[branch_idx].set_xticks([data_x[0], data_x[-1]], labels = tick_labels)

		if branch_idx == num_branches / 2 - 1:
			axes[branch_idx].set_xlabel("Q (rlu)")

		# scale and clamp S
		data_S *= S_scale
		if S_clamp_min < S_clamp_max:
			data_S = numpy.clip(data_S, a_min = S_clamp_min, a_max = S_clamp_max)

		# plot the dispersion branch
		axes[branch_idx].scatter(data_x, data_E, marker = '.', s = data_S)

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0)

	if plot_file != "":
		pyplot.savefig(plot_file)
	pyplot.show()
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# data
# -----------------------------------------------------------------------------
all_data = %%DATA%%
# -----------------------------------------------------------------------------


if __name__ == "__main__":
	plot_disp(all_data)
)RAW";
// ============================================================================

}  // tl2_mag


// --------------------------------------------------------------------
// saving of the dispersion data
// --------------------------------------------------------------------
/**
 * generates the dispersion plot along the given Q path
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveDispersion(const std::string& filename,
	t_real h_start, t_real k_start, t_real l_start,
	t_real h_end, t_real k_end, t_real l_end,
	t_size num_Qs, t_size num_threads, bool as_py,
	std::function<bool(int, int)> *progress_fkt) const
{
	std::ofstream ofstr{filename};
	if(!ofstr)
		return false;

	return SaveDispersion(ofstr,
		h_start, k_start, l_start,
		h_end, k_end, l_end, num_Qs,
		num_threads, as_py, progress_fkt,
		true);
}



/**
 * generates the dispersion plot along multiple Q paths
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveMultiDispersion(const std::string& filename,
	const std::vector<t_vec_real>& Qs,
	t_size num_Qs, t_size num_threads, bool as_py,
	std::function<bool(int, int)> *progress_fkt,
	const std::vector<std::string> *Q_names) const
{
	std::ofstream ofstr{filename};
	if(!ofstr)
		return false;

	return SaveMultiDispersion(ofstr,
		Qs, num_Qs, num_threads, as_py,
		progress_fkt, Q_names);
}



/**
 * generates the dispersion plot along the given Q path
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveDispersion(std::ostream& ostr,
	t_real h_start, t_real k_start, t_real l_start,
	t_real h_end, t_real k_end, t_real l_end,
	t_size num_Qs, t_size num_threads, bool as_py,
	std::function<bool(int, int)> *progress_fkt,
	bool write_header) const
{
	ostr.precision(m_prec);
	int field_len = m_prec * 2.5;

	// data for script export
	std::ostringstream all_data;
	all_data.precision(m_prec);

	if(write_header)
	{
		ostr << "#\n# Created with Takin/Magdyn.\n";
		ostr << "# DOI: https://doi.org/10.5281/zenodo.4117437\n";
		ostr << "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n";
		ostr << "#\n\n";
	}

	if(!as_py)  // save as text file
	{
		ostr	<< std::setw(field_len) << std::left << "# h" << " "
			<< std::setw(field_len) << std::left << "k" << " "
			<< std::setw(field_len) << std::left << "l" << " "
			<< std::setw(field_len) << std::left << "E" << " "
			<< std::setw(field_len) << std::left << "S(Q,E)" << " "
			<< std::setw(field_len) << std::left << "S_xx" << " "
			<< std::setw(field_len) << std::left << "S_yy" << " "
			<< std::setw(field_len) << std::left << "S_zz"
			<< "\n";
	}

	SofQEs results = CalcDispersion(h_start, k_start, l_start,
		h_end, k_end, l_end, num_Qs, num_threads, progress_fkt);

	// print results
	for(const auto& result : results)
	{
		if(progress_fkt && !(*progress_fkt)(-1, -1))
			return false;

		// get results
		for(const EnergyAndWeight& E_and_S : result.E_and_S)
		{
			if(!as_py)  // save as text file
			{
				ostr	<< std::setw(field_len) << std::left << result.h << " "
					<< std::setw(field_len) << std::left << result.k << " "
					<< std::setw(field_len) << std::left << result.l << " "
					<< std::setw(field_len) << E_and_S.E << " "
					<< std::setw(field_len) << E_and_S.weight << " "
					<< std::setw(field_len) << E_and_S.S_perp(0, 0).real() << " "
					<< std::setw(field_len) << E_and_S.S_perp(1, 1).real() << " "
					<< std::setw(field_len) << E_and_S.S_perp(2, 2).real()
					<< "\n";
			}
			else        // save as py script
			{
				all_data << "\t"
					<< "[ " << result.h
					<< ", " << result.k
					<< ", " << result.l
					<< ", " << E_and_S.E
					<< ", " << E_and_S.weight
					<< " ],\n";
			}
		}
	}

	if(as_py)  // save as py script
	{
		std::string pyscr = g_pyscr<std::string>;

		namespace algo = boost::algorithm;
		algo::replace_all(pyscr, "%%LABELS%%", "None");
		algo::replace_all(pyscr, "%%RATIOS%%", "None");
		algo::replace_all(pyscr, "%%DATA%%", "[[\n" + all_data.str() + "\n]]");

		ostr << pyscr << "\n";
	}

	ostr.flush();
	return true;
}



/**
 * generates the dispersion plot along multiple Q paths
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveMultiDispersion(std::ostream& ostr,
	const std::vector<t_vec_real>& Qs,
	t_size num_Qs, t_size num_threads, bool as_py,
	std::function<bool(int, int)> *progress_fkt,
	const std::vector<std::string> *Q_names) const
{
	bool ok = true;
	ostr.precision(m_prec);

	const t_size N = Qs.size();

	ostr << "#\n# Created with Takin/Magdyn.\n";
	ostr << "# DOI: https://doi.org/10.5281/zenodo.4117437\n";
	ostr << "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n";
	ostr << "#\n\n";

	// data for script export
	bool has_Q_names = false;
	std::ostringstream all_data, Q_ratios, Q_labels;
	for(std::ostringstream* theostr : { &all_data, &Q_ratios, &Q_labels })
		theostr->precision(m_prec);

	for(t_size i = 0; i < N - 1; ++i)
	{
		const t_vec_real& Q1 = Qs[i];
		const t_vec_real& Q2 = Qs[i + 1];

		// length from Q1 to Q2
		t_vec_real Q1_lab = m_xtalUB * Q1;
		t_vec_real Q2_lab = m_xtalUB * Q2;
		Q_ratios << tl2::norm<t_vec_real>(Q2_lab - Q1_lab) << ", ";

		if(!as_py)  // save as text file
		{
			ostr << "# ";
			if(Q_names && (*Q_names)[i] != "" && (*Q_names)[i + 1] != "")
			{
				ostr << (*Q_names)[i] << " -> " << (*Q_names)[i + 1]
					<< ": ";
			}
			ostr << "(" << Q1[0] << ", " << Q1[1] << ", " << Q1[2]
				<< ") -> (" << Q2[0] << ", " << Q2[1] << ", " << Q2[2]
				<< ")\n";

			if(!SaveDispersion(ostr, Q1[0], Q1[1], Q1[2],
				Q2[0], Q2[1], Q2[2], num_Qs,
				num_threads, false, progress_fkt, false))
			{
				ok = false;
				break;
			}

			ostr << "\n";
		}
		else        // save as py script
		{
			// get Q names
			if(i == 0 && Q_names && (*Q_names)[i] != "")
			{
				Q_labels << "\"" << (*Q_names)[i] << "\", ";
				has_Q_names = true;
			}
			if(Q_names && (*Q_names)[i + 1] != "")
			{
				Q_labels << "\"" << (*Q_names)[i + 1] << "\", ";
				has_Q_names = true;
			}

			SofQEs results = CalcDispersion(Q1[0], Q1[1], Q1[2],
				Q2[0], Q2[1], Q2[2], num_Qs,
				num_threads, progress_fkt);

			if(progress_fkt && !(*progress_fkt)(-1, -1))
				break;

			// get results
			all_data << "[";
			for(const auto& result : results)
			{
				for(const EnergyAndWeight& E_and_S : result.E_and_S)
				{
					all_data << "\t"
						<< "[ " << result.h
						<< ", " << result.k
						<< ", " << result.l
						<< ", " << E_and_S.E
						<< ", " << E_and_S.weight
						<< " ],\n";
				}
			}
			all_data << "],\n";
		}
	}

	if(as_py)  // save as py script
	{
		std::string pyscr = g_pyscr<std::string>;

		namespace algo = boost::algorithm;
		algo::replace_all(pyscr, "%%LABELS%%",
			has_Q_names ? "[ " + Q_labels.str() + " ]" : "None");
		algo::replace_all(pyscr, "%%RATIOS%%", "[ " + Q_ratios.str() + " ]");
		algo::replace_all(pyscr, "%%DATA%%", "[\n" + all_data.str() + "\n]");

		ostr << pyscr << "\n";
	}

	ostr.flush();
	return ok;
}
// --------------------------------------------------------------------

#endif
