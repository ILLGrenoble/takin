/**
 * magnetic dynamics -- calculations for dispersion plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
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

// these need to be included before all other things on mingw
#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include "magdyn.h"
#include "helper.h"

#include <QtWidgets/QApplication>

#include <sstream>
#include <thread>
#include <future>
#include <mutex>

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"

using namespace tl2_ops;



/**
 * clears the dispersion graph
 */
void MagDynDlg::ClearDispersion(bool replot)
{
	m_graphs.clear();

	if(m_plot)
	{
		m_plot->clearPlottables();
		if(replot)
			m_plot->replot();
	}

	m_qs_data.clear();
	m_Es_data.clear();
	m_ws_data.clear();

	for(int i = 0; i < 3; ++i)
	{
		m_qs_data_channel[i].clear();
		m_Es_data_channel[i].clear();
		m_ws_data_channel[i].clear();
	}

	m_Q_idx = 0;
}



/**
 * draw the calculated dispersion curve
 */
void MagDynDlg::PlotDispersion()
{
	if(!m_plot)
		return;

	m_plot->clearPlottables();
	m_graphs.clear();

	const bool plot_channels = m_plot_channels->isChecked();

	// create plot curves
	if(plot_channels)
	{
		const QColor colChannel[3]
		{
			QColor(0xff, 0x00, 0x00),
			QColor(0x00, 0xff, 0x00),
			QColor(0x00, 0x00, 0xff),
		};

		for(int i = 0; i < 3; ++i)
		{
			if(!m_plot_channel[i]->isChecked())
				continue;

			GraphWithWeights *graph = new GraphWithWeights(m_plot->xAxis, m_plot->yAxis);
			QPen pen = graph->pen();
			pen.setColor(colChannel[i]);
			pen.setWidthF(1.);
			graph->setPen(pen);
			graph->setBrush(QBrush(pen.color(), Qt::SolidPattern));
			graph->setLineStyle(QCPGraph::lsNone);
			graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, m_weight_scale->value()));
			graph->setAntialiased(true);
			graph->setData(m_qs_data_channel[i], m_Es_data_channel[i], true /*already sorted*/);
			graph->SetWeights(m_ws_data_channel[i]);
			graph->SetWeightScale(m_weight_scale->value(), m_weight_min->value(), m_weight_max->value());
			graph->SetWeightAsPointSize(m_plot_weights_pointsize->isChecked());
			graph->SetWeightAsAlpha(m_plot_weights_alpha->isChecked());
			m_graphs.push_back(graph);
		}
	}
	else
	{
		GraphWithWeights *graph = new GraphWithWeights(m_plot->xAxis, m_plot->yAxis);
		QPen pen = graph->pen();

		// colour
		int col_comp[3] = { 0xff, 0, 0 };
		get_colour<int>(g_colPlot, col_comp);
		const QColor colFull(col_comp[0], col_comp[1], col_comp[2]);
		pen.setColor(colFull);

		pen.setWidthF(1.);
		graph->setPen(pen);
		graph->setBrush(QBrush(pen.color(), Qt::SolidPattern));
		graph->setLineStyle(QCPGraph::lsNone);
		graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, m_weight_scale->value()));
		graph->setAntialiased(true);
		graph->setData(m_qs_data, m_Es_data, true /*already sorted*/);
		graph->SetWeights(m_ws_data);
		graph->SetWeightScale(m_weight_scale->value(), m_weight_min->value(), m_weight_max->value());
		graph->SetWeightAsPointSize(m_plot_weights_pointsize->isChecked());
		graph->SetWeightAsAlpha(m_plot_weights_alpha->isChecked());
		m_graphs.push_back(graph);
	}

	// set labels
	const char* Q_label[]{ "h (rlu)", "k (rlu)", "l (rlu)" };
	m_plot->xAxis->setLabel(Q_label[m_Q_idx]);

	// set ranges
	auto [min_E_iter, max_E_iter] = std::minmax_element(m_Es_data.begin(), m_Es_data.end());
	m_plot->xAxis->setRange(m_Q_min, m_Q_max);
	if(min_E_iter != m_Es_data.end() && max_E_iter != m_Es_data.end())
	{
		t_real E_range = *max_E_iter - *min_E_iter;

		t_real Emin = *min_E_iter - E_range*0.05;
		t_real Emax = *max_E_iter + E_range*0.05;
		//if(ignore_annihilation)
		//	Emin = t_real(0);

		m_plot->yAxis->setRange(Emin, Emax);
	}
	else
	{
		m_plot->yAxis->setRange(0., 1.);
	}

	m_plot->replot();
}



/**
 * calculate the dispersion branches
 */
void MagDynDlg::CalcDispersion()
{
	if(m_ignoreCalc)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableInput();
	} BOOST_SCOPE_EXIT_END
	DisableInput();

	// nothing to calculate?
	if(m_dyn.GetMagneticSitesCount() == 0 || m_dyn.GetExchangeTermsCount() == 0)
	{
		ClearDispersion(true);
		return;
	}
	else
	{
		ClearDispersion(false);
	}

	t_real Q_start[]
	{
		(t_real)m_Q_start[0]->value(),
		(t_real)m_Q_start[1]->value(),
		(t_real)m_Q_start[2]->value(),
	};

	t_real Q_end[]
	{
		(t_real)m_Q_end[0]->value(),
		(t_real)m_Q_end[1]->value(),
		(t_real)m_Q_end[2]->value(),
	};

	const t_real Q_range[]
	{
		std::abs(Q_end[0] - Q_start[0]),
		std::abs(Q_end[1] - Q_start[1]),
		std::abs(Q_end[2] - Q_start[2]),
	};

	// Q component with maximum range
	m_Q_idx = 0;
	if(Q_range[1] > Q_range[m_Q_idx])
		m_Q_idx = 1;
	if(Q_range[2] > Q_range[m_Q_idx])
		m_Q_idx = 2;

	t_size num_pts = m_num_points->value();

	m_qs_data.reserve(num_pts*10);
	m_Es_data.reserve(num_pts*10);
	m_ws_data.reserve(num_pts*10);

	for(int i = 0; i < 3; ++i)
	{
		m_qs_data_channel[i].reserve(num_pts*10);
		m_Es_data_channel[i].reserve(num_pts*10);
		m_ws_data_channel[i].reserve(num_pts*10);
	}

	m_Q_min = Q_start[m_Q_idx];
	m_Q_max = Q_end[m_Q_idx];

	// options
	const bool is_comm = !m_dyn.IsIncommensurate();
	const bool use_min_E = false;
	const bool unite_degeneracies = m_unite_degeneracies->isChecked();
	const bool ignore_annihilation = m_ignore_annihilation->isChecked();
	const bool use_weights = m_use_weights->isChecked();
	const bool use_projector = m_use_projector->isChecked();
	const bool force_incommensurate = m_force_incommensurate->isChecked();

	t_real E0 = use_min_E ? m_dyn.CalcMinimumEnergy() : 0.;
	m_dyn.SetUniteDegenerateEnergies(unite_degeneracies);
	m_dyn.SetForceIncommensurate(force_incommensurate);
	m_dyn.SetCalcHamiltonian(
		m_hamiltonian_comp[0]->isChecked() || is_comm,  // always calculate commensurate case
		m_hamiltonian_comp[1]->isChecked(),
		m_hamiltonian_comp[2]->isChecked());

	// tread pool
	asio::thread_pool pool{g_num_threads};

	// mutex to protect m_qs_data, m_Es_data, and m_ws_data
	std::mutex mtx;

	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(num_pts);

	// keep the scanned Q component in ascending order
	if(Q_start[m_Q_idx] > Q_end[m_Q_idx])
	{
		std::swap(Q_start[0], Q_end[0]);
		std::swap(Q_start[1], Q_end[1]);
		std::swap(Q_start[2], Q_end[2]);
	}

	m_stopRequested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(num_pts);
	m_progress->setValue(0);
	m_status->setText("Starting calculation.");
	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	for(t_size i = 0; i < num_pts; ++i)
	{
		auto task = [this, &mtx, i, num_pts, E0,
			use_projector, use_weights, ignore_annihilation, &Q_start, &Q_end]()
		{
			const t_vec_real Q = num_pts > 1
				? tl2::create<t_vec_real>(
				{
					std::lerp(Q_start[0], Q_end[0], t_real(i) / t_real(num_pts - 1)),
					std::lerp(Q_start[1], Q_end[1], t_real(i) / t_real(num_pts - 1)),
					std::lerp(Q_start[2], Q_end[2], t_real(i) / t_real(num_pts - 1)),
				})
				: tl2::create<t_vec_real>({ Q_start[0], Q_start[1], Q_start[2] });

			auto energies_and_correlations = m_dyn.CalcEnergies(Q, !use_weights);

			for(const auto& E_and_S : energies_and_correlations)
			{
				if(m_stopRequested)
					break;

				t_real E = E_and_S.E - E0;
				if(std::isnan(E) || std::isinf(E))
					continue;
				if(ignore_annihilation && E < t_real(0))
					continue;

				std::lock_guard<std::mutex> _lck{mtx};

				// weights
				if(use_weights)
				{
					t_real weight = use_projector ? E_and_S.weight : E_and_S.weight_full;
					if(std::isnan(weight) || std::isinf(weight))
						continue;

					m_ws_data.push_back(weight);

					for(int channel = 0; channel < 3; ++channel)
					{
						t_real weight_channel = use_projector
							? E_and_S.S_perp(channel, channel).real()
							: E_and_S.S(channel, channel).real();

						weight_channel = std::abs(weight_channel);

						if(!tl2::equals_0<t_real>(weight_channel, g_eps))
						{
							m_Es_data_channel[channel].push_back(E);
							m_qs_data_channel[channel].push_back(Q[m_Q_idx]);
							m_ws_data_channel[channel].push_back(weight_channel);
						}
					}
				}

				m_qs_data.push_back(Q[m_Q_idx]);
				m_Es_data.push_back(E);
			}
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText("Calculating dispersion.");

	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		qApp->processEvents();  // process events to see if the stop button was clicked
		if(m_stopRequested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	auto sort_data = [](QVector<qreal>& qvec, QVector<qreal>& Evec, QVector<qreal>& wvec)
	{
		// sort vectors by q component
		std::vector<std::size_t> perm = tl2::get_perm(qvec.size(),
			[&qvec, &Evec](std::size_t idx1, std::size_t idx2) -> bool
		{
			// if q components are equal, sort by E
			if(tl2::equals(qvec[idx1], qvec[idx2], (qreal)g_eps))
				return Evec[idx1] < Evec[idx2];
			return qvec[idx1] < qvec[idx2];
		});

		qvec = tl2::reorder(qvec, perm);
		Evec = tl2::reorder(Evec, perm);
		wvec = tl2::reorder(wvec, perm);
	};

	sort_data(m_qs_data, m_Es_data, m_ws_data);
	for(int i = 0; i < 3; ++i)
		sort_data(m_qs_data_channel[i], m_Es_data_channel[i], m_ws_data_channel[i]);

	PlotDispersion();
}



/**
 * calculate the hamiltonian for a single Q value
 */
void MagDynDlg::CalcHamiltonian()
{
	if(m_ignoreCalc)
		return;

	// options
	const bool only_energies = !m_use_weights->isChecked();
	const bool use_projector = m_use_projector->isChecked();
	const bool ignore_annihilation = m_ignore_annihilation->isChecked();
	const bool unite_degeneracies = m_unite_degeneracies->isChecked();
	const bool force_incommensurate = m_force_incommensurate->isChecked();

	m_dyn.SetUniteDegenerateEnergies(unite_degeneracies);
	m_dyn.SetForceIncommensurate(force_incommensurate);

	m_hamiltonian->clear();

	const t_vec_real Q = tl2::create<t_vec_real>(
	{
		(t_real)m_Q[0]->value(),
		(t_real)m_Q[1]->value(),
		(t_real)m_Q[2]->value(),
	});

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);

	// print hamiltonian
	auto print_H = [&ostr](const t_mat& H, const t_vec_real& Qvec,
		const std::string& Qstr = "Q", const std::string& mQstr = "-Q")
	{
		ostr << "<p><h3>Hamiltonian at " << Qstr <<  " = ("
			<< Qvec[0] << ", " << Qvec[1] << ", " << Qvec[2] << ")</h3>";
		ostr << "<table style=\"border:0px\">";

		// horizontal header
		ostr << "<tr><th/>";
		for(std::size_t col = 0; col < H.size2()/2; ++col)
		{
			ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (col + 1)
				<< "</sub>(" << Qstr << ")" << "</th>";
		}
		for(std::size_t col = H.size2()/2; col < H.size2(); ++col)
		{
			ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (col - H.size2()/2 + 1)
				<< "</sub><sup>&#x2020;</sup>(" << mQstr << ")" << "</th>";
		}
		ostr << "</tr>";

		for(std::size_t row = 0; row < H.size1(); ++row)
		{
			ostr << "<tr>";

			// vertical header
			if(row < H.size1() / 2)
			{
				ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (row + 1)
					<< "</sub><sup>&#x2020;</sup>(" << Qstr << ")" << "</th>";
			}
			else
			{
				ostr << "<th style=\"padding-right:8px\">" << "b<sub>" << (row - H.size1()/2 + 1)
					<< "</sub>(" << mQstr << ")" << "</th>";
			}

			// components
			for(std::size_t col = 0; col < H.size2(); ++col)
			{
				t_cplx elem = H(row, col);
				tl2::set_eps_0<t_cplx, t_real>(elem, g_eps);
				ostr << "<td style=\"padding-right:8px\">"
					<< elem << "</td>";
			}

			ostr << "</tr>";
		}
		ostr << "</table></p>";
	};

	// get hamiltonian at Q
	t_mat H = m_dyn.CalcHamiltonian(Q);
	const bool is_comm = !m_dyn.IsIncommensurate();
	if(m_hamiltonian_comp[0]->isChecked() || is_comm)  // always calculate commensurate case
		print_H(H, Q, "Q", "-Q");

	// print shifted hamiltonians for incommensurate case
	if(!is_comm)
	{
		const t_vec_real O = tl2::create<t_vec_real>(
		{
			(t_real)m_ordering[0]->value(),
			(t_real)m_ordering[1]->value(),
			(t_real)m_ordering[2]->value(),
		});

		if(!tl2::equals_0<t_vec_real>(O, g_eps))
		{
			if(m_hamiltonian_comp[1]->isChecked())
			{
				// get hamiltonian at Q + ordering vector
				t_mat H_p = m_dyn.CalcHamiltonian(Q + O);
				print_H(H_p, Q + O, "Q + O", "Q - O");
			}

			if(m_hamiltonian_comp[2]->isChecked())
			{
				// get hamiltonian at Q - ordering vector
				t_mat H_m = m_dyn.CalcHamiltonian(Q - O);
				print_H(H_m, Q - O, "Q - O", "Q + O");
			}
		}
	}

	// get energies and correlation functions
	using t_E_and_S = typename decltype(m_dyn)::EnergyAndWeight;
	std::vector<t_E_and_S> energies_and_correlations;

	if(is_comm)
	{
		// commensurate case
		energies_and_correlations = m_dyn.CalcEnergiesFromHamiltonian(H, Q, only_energies);
		if(!only_energies)
			m_dyn.CalcIntensities(Q, energies_and_correlations);
		if(unite_degeneracies)
			energies_and_correlations = m_dyn.UniteEnergies(energies_and_correlations);
	}
	else
	{
		// incommensurate case
		energies_and_correlations = m_dyn.CalcEnergies(Q, only_energies);
	}

	if(only_energies)
	{
		// split into positive and negative energies
		std::vector<t_magdyn::EnergyAndWeight> Es_neg, Es_pos;
		for(const t_E_and_S& E_and_S : energies_and_correlations)
		{
			t_real E = E_and_S.E;

			if(E < t_real(0))
				Es_neg.push_back(E_and_S);
			else
				Es_pos.push_back(E_and_S);
		}

		std::stable_sort(Es_neg.begin(), Es_neg.end(),
			[](const t_E_and_S& E_and_S_1, const t_E_and_S& E_and_S_2) -> bool
		{
			t_real E1 = E_and_S_1.E;
			t_real E2 = E_and_S_2.E;
			return std::abs(E1) < std::abs(E2);
		});

		std::stable_sort(Es_pos.begin(), Es_pos.end(),
			[](const t_E_and_S& E_and_S_1, const t_E_and_S& E_and_S_2) -> bool
		{
			t_real E1 = E_and_S_1.E;
			t_real E2 = E_and_S_2.E;
			return std::abs(E1) < std::abs(E2);
		});

		ostr << "<p><h3>Energies</h3>";
		ostr << "<table style=\"border:0px\">";
		ostr << "<tr>";
		ostr << "<th style=\"padding-right:8px\">Creation</th>";
		for(const t_E_and_S& E_and_S : Es_pos)
		{
			t_real E = E_and_S.E;
			tl2::set_eps_0(E);

			ostr << "<td style=\"padding-right:8px\">"
				<< E << " meV" << "</td>";
		}
		ostr << "</tr>";

		if(!ignore_annihilation)
		{
			ostr << "<tr>";
			ostr << "<th style=\"padding-right:8px\">Annihilation</th>";
			for(const t_E_and_S& E_and_S : Es_neg)
			{
				t_real E = E_and_S.E;
				tl2::set_eps_0(E);

				ostr << "<td style=\"padding-right:8px\">"
					<< E << " meV" << "</td>";
			}
			ostr << "</tr>";
		}

		ostr << "</table></p>";
	}
	else
	{
		std::stable_sort(energies_and_correlations.begin(), energies_and_correlations.end(),
			[](const t_E_and_S& E_and_S_1, const t_E_and_S& E_and_S_2) -> bool
		{
			t_real E1 = E_and_S_1.E;
			t_real E2 = E_and_S_2.E;
			return E1 < E2;
		});

		ostr << "<p><h3>Spectrum</h3>";
		ostr << "<table style=\"border:0px\">";
		ostr << "<tr>";
		ostr << "<th style=\"padding-right:16px\">Energy E</td>";
		ostr << "<th style=\"padding-right:16px\">Correlation S(Q, E)</td>";
		ostr << "<th style=\"padding-right:16px\">Projection S<sub>&#x27C2;</sub>(Q, E)</td>";
		ostr << "<th style=\"padding-right:16px\">Weight</td>";
		ostr << "</tr>";

		for(const t_E_and_S& E_and_S : energies_and_correlations)
		{
			t_real E = E_and_S.E;
			if(ignore_annihilation && E < t_real(0))
				continue;

			const t_mat& S = E_and_S.S;
			const t_mat& S_perp = E_and_S.S_perp;
			t_real weight = E_and_S.weight;
			if(!use_projector)
				weight = tl2::trace<t_mat>(S).real();

			tl2::set_eps_0(E);
			tl2::set_eps_0(weight);

			// E
			ostr << "<tr>";
			ostr << "<td style=\"padding-right:16px\">"
				<< E << " meV" << "</td>";

			// S(Q, E)
			ostr << "<td style=\"padding-right:16px\">";
			ostr << "<table style=\"border:0px\">";
			for(std::size_t i=0; i<S.size1(); ++i)
			{
				ostr << "<tr>";
				for(std::size_t j=0; j<S.size2(); ++j)
				{
					t_cplx elem = S(i, j);
					tl2::set_eps_0<t_cplx, t_real>(elem, g_eps);
					ostr << "<td style=\"padding-right:8px\">"
						<< elem << "</td>";
				}
				ostr << "</tr>";
			}
			ostr << "</table>";
			ostr << "</td>";

			// S_perp(Q, E)
			ostr << "<td style=\"padding-right:16px\">";
			ostr << "<table style=\"border:0px\">";
			for(std::size_t i=0; i<S_perp.size1(); ++i)
			{
				ostr << "<tr>";
				for(std::size_t j=0; j<S_perp.size2(); ++j)
				{
					t_cplx elem = S_perp(i, j);
					tl2::set_eps_0<t_cplx, t_real>(elem, g_eps);
					ostr << "<td style=\"padding-right:8px\">"
						<< elem << "</td>";
				}
				ostr << "</tr>";
			}
			ostr << "</table>";
			ostr << "</td>";

			// tr(S_perp(Q, E))
			ostr << "<td style=\"padding-right:16px\">" << weight << "</td>";

			ostr << "</tr>";
		}
		ostr << "</table></p>";
	}

	m_hamiltonian->setHtml(ostr.str().c_str());
}



/**
 * set the number of Q points on the dispersion to calculate
 */
void MagDynDlg::SetNumQPoints(t_size num_Q_pts)
{
	m_num_points->setValue((int)num_Q_pts);
}



/**
 * set the current dispersion path and the hamiltonian to the given one
 */
void MagDynDlg::SetCoordinates(const t_vec_real& Qi, const t_vec_real& Qf, bool calc_dynamics)
{
	m_ignoreCalc = true;

	BOOST_SCOPE_EXIT(this_, calc_dynamics)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked() && calc_dynamics)
		{
			this_->CalcDispersion();
			this_->CalcHamiltonian();
		}
	} BOOST_SCOPE_EXIT_END

	// calculate the dispersion from Qi to Qf
	m_Q_start[0]->setValue(Qi[0]);
	m_Q_start[1]->setValue(Qi[1]);
	m_Q_start[2]->setValue(Qi[2]);
	m_Q_end[0]->setValue(Qf[0]);
	m_Q_end[1]->setValue(Qf[1]);
	m_Q_end[2]->setValue(Qf[2]);

	// calculate the hamiltonian for Qi
	m_Q[0]->setValue(Qi[0]);
	m_Q[1]->setValue(Qi[1]);
	m_Q[2]->setValue(Qi[2]);
}



/**
 * set the selected coordinate path as the current one
 */
void MagDynDlg::SetCurrentCoordinate(int which)
{
	using t_item = tl2::NumericTableWidgetItem<t_real>;

	int idx_i = m_coordinates_cursor_row;
	if(idx_i < 0 || idx_i >= m_coordinatestab->rowCount())
		return;

	const auto* hi = static_cast<t_item*>(m_coordinatestab->item(idx_i, COL_COORD_H));
	const auto* ki = static_cast<t_item*>(m_coordinatestab->item(idx_i, COL_COORD_K));
	const auto* li = static_cast<t_item*>(m_coordinatestab->item(idx_i, COL_COORD_L));

	// set dispersion start and end coordinates
	if(which == 0)
	{
		int idx_f = idx_i + 1;

		// wrap around
		if(idx_f >= m_coordinatestab->rowCount())
			idx_f = 0;

		if(idx_f == idx_i)
			return;
		if(idx_f < 0 || idx_f >= m_coordinatestab->rowCount())
			return;

		const auto* hf = static_cast<t_item*>(m_coordinatestab->item(idx_f, COL_COORD_H));
		const auto* kf = static_cast<t_item*>(m_coordinatestab->item(idx_f, COL_COORD_K));
		const auto* lf = static_cast<t_item*>(m_coordinatestab->item(idx_f, COL_COORD_L));

		if(!hi || !ki || !li || !hf || !kf || !lf)
			return;

		m_ignoreCalc = true;

		BOOST_SCOPE_EXIT(this_)
		{
			this_->m_ignoreCalc = false;
			if(this_->m_autocalc->isChecked())
				this_->CalcDispersion();
		} BOOST_SCOPE_EXIT_END

		m_Q_start[0]->setValue(hi->GetValue());
		m_Q_start[1]->setValue(ki->GetValue());
		m_Q_start[2]->setValue(li->GetValue());
		m_Q_end[0]->setValue(hf->GetValue());
		m_Q_end[1]->setValue(kf->GetValue());
		m_Q_end[2]->setValue(lf->GetValue());
	}

	// send initial Q coordinates to hamiltonian calculation
	else if(which == 1)
	{
		if(!hi || !ki || !li)
			return;

		m_ignoreCalc = true;

		BOOST_SCOPE_EXIT(this_)
		{
			this_->m_ignoreCalc = false;
			if(this_->m_autocalc->isChecked())
				this_->CalcHamiltonian();
		} BOOST_SCOPE_EXIT_END

		m_Q[0]->setValue(hi->GetValue());
		m_Q[1]->setValue(ki->GetValue());
		m_Q[2]->setValue(li->GetValue());
	}
}



/**
 * mouse move event of the plot
 */
void MagDynDlg::PlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot->xAxis->pixelToCoord(evt->pos().x());
	t_real E = m_plot->yAxis->pixelToCoord(evt->pos().y());

	QString status("Q = %1 rlu, E = %2 meV.");
	status = status.arg(Q, 0, 'g', g_prec_gui).arg(E, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * mouse button has been pressed in the plot
 */
void MagDynDlg::PlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuDisp)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuDisp->popup(pos);
		evt->accept();
	}
}
