/**
 * magnetic dynamics -- export of the calculated S(Q, E)
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

#include <QtCore/QString>

#include <iostream>
#include <fstream>
#include <sstream>
#include <future>
#include <vector>
#include <deque>
#include <cstdlib>

#include "tlibs2/libs/str.h"
#include "tlibs2/libs/algos.h"

#ifdef USE_HDF5
	#include "tlibs2/libs/h5file.h"
#endif


// for debugging: write individual data chunks in hdf5 file
//#define WRITE_HDF5_CHUNKS

// precision
extern int g_prec;



/**
 * show dialog and export S(Q, E) into a grid
 */
void MagDynDlg::ExportSQE()
{
	QString extension;
	switch(m_exportFormat->currentData().toInt())
	{
		case EXPORT_HDF5: extension = "HDF5 Files (*.hdf)"; break;
		case EXPORT_GRID: extension = "Takin Grid Files (*.bin)"; break;
		case EXPORT_TEXT: extension = "Text Files (*.txt)"; break;
	}

	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Export S(Q,E)", dirLast, extension);
	if(filename == "")
		return;

	if(ExportSQE(filename))
		m_sett->setValue("dir", QFileInfo(filename).path());
}



/**
 * export S(Q, E) into a grid
 */
bool MagDynDlg::ExportSQE(const QString& filename)
{
#ifdef USE_HDF5
	std::unique_ptr<H5::H5File> h5file;
#endif
	std::unique_ptr<std::ofstream> ofstr;
	bool file_opened = false;

	const int format = m_exportFormat->currentData().toInt();
	if(format == EXPORT_GRID || format == EXPORT_TEXT)
	{
		ofstr = std::make_unique<std::ofstream>(filename.toStdString());
		ofstr->precision(g_prec);
		file_opened = ofstr->operator bool();
	}

#ifdef USE_HDF5
	else if(format == EXPORT_HDF5)
	{
		h5file = std::make_unique<H5::H5File>(filename.toStdString().c_str(), H5F_ACC_TRUNC);
		h5file->createGroup("meta_infos");
		h5file->createGroup("infos");
		h5file->createGroup("data");
#ifdef WRITE_HDF5_CHUNKS
		h5file->createGroup("chunks");
#endif

		file_opened = true;
	}
#endif

	if(!file_opened)
	{
		QMessageBox::critical(this, "Magnetic Dynamics", "File could not be opened.");
		return false;
	}

	const t_vec_real Qstart = tl2::create<t_vec_real>({
		(t_real)m_exportStartQ[0]->value(),
		(t_real)m_exportStartQ[1]->value(),
		(t_real)m_exportStartQ[2]->value() });
	const t_vec_real Qend = tl2::create<t_vec_real>({
		(t_real)m_exportEndQ[0]->value(),
		(t_real)m_exportEndQ[1]->value(),
		(t_real)m_exportEndQ[2]->value() });

	const t_size num_pts_h = m_exportNumPoints[0]->value();
	const t_size num_pts_k = m_exportNumPoints[1]->value();
	const t_size num_pts_l = m_exportNumPoints[2]->value();

	t_magdyn dyn = m_dyn;
	dyn.SetUniteDegenerateEnergies(m_unite_degeneracies->isChecked());
	bool use_weights = m_use_weights->isChecked();
	bool use_projector = m_use_projector->isChecked();

	const t_vec_real dir = Qend - Qstart;
	const t_real inc_h = dir[0] / t_real(num_pts_h);
	const t_real inc_k = dir[1] / t_real(num_pts_k);
	const t_real inc_l = dir[2] / t_real(num_pts_l);
	const t_vec_real Qstep = tl2::create<t_vec_real>({inc_h, inc_k, inc_l});

	// tread pool
	asio::thread_pool pool{g_num_threads};


	using t_taskret = std::deque<
		std::tuple<t_real, t_real, t_real,          // h, k, l
		std::vector<t_real>, std::vector<t_real>,   // E, S
		std::size_t, std::size_t, std::size_t>>;    // h_idx, k_idx, l_idx

	// calculation task
	auto task = [this, use_weights, use_projector, &dyn, inc_l, num_pts_l]
		(t_real h_pos, t_real k_pos, t_real l_pos, std::size_t h_idx, std::size_t k_idx) -> t_taskret
	{
		t_taskret ret;

		// iterate last Q dimension
		for(std::size_t l_idx=0; l_idx<num_pts_l; ++l_idx)
		{
			t_real l = l_pos + inc_l*t_real(l_idx);
			auto energies_and_correlations = dyn.CalcEnergies(
				h_pos, k_pos, l, !use_weights);

			std::vector<t_real> Es, weights;
			Es.reserve(energies_and_correlations.size());
			weights.reserve(energies_and_correlations.size());

			for(const auto& E_and_S : energies_and_correlations)
			{
				if(m_stopRequested)
					break;

				t_real E = E_and_S.E;
				if(std::isnan(E) || std::isinf(E))
					continue;

				const t_mat& S = E_and_S.S;
				t_real weight = E_and_S.weight;

				if(!use_projector)
					weight = tl2::trace<t_mat>(S).real();

				if(std::isnan(weight) || std::isinf(weight))
					weight = 0.;

				Es.push_back(E);
				weights.push_back(weight);
			}

			ret.emplace_back(std::make_tuple(h_pos, k_pos, l, Es, weights, h_idx, k_idx, l_idx));
		}

		return ret;
	};


	// tasks
	using t_task = std::packaged_task<t_taskret(t_real, t_real, t_real, std::size_t, std::size_t)>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::deque<t_taskptr> tasks;
	std::deque<std::future<t_taskret>> futures;

	m_stopRequested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(num_pts_h * num_pts_k /** num_pts_l*/);
	m_progress->setValue(0);
	m_status->setText("Starting calculation.");
	DisableInput();
	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	std::size_t task_idx = 0;
	// iterate first two Q dimensions
	for(std::size_t h_idx=0; h_idx<num_pts_h; ++h_idx)
	{
		if(m_stopRequested)
			break;

		for(std::size_t k_idx=0; k_idx<num_pts_k; ++k_idx)
		{
			if(m_stopRequested)
				break;

			qApp->processEvents();  // process events to see if the stop button was clicked
			if(m_stopRequested)
				break;

			t_vec_real Q = Qstart;
			Q[0] += inc_h*t_real(h_idx);
			Q[1] += inc_k*t_real(k_idx);
			//Q[2] += inc_l*t_real(l_idx);

			// create tasks
			t_taskptr taskptr = std::make_shared<t_task>(task);
			tasks.push_back(taskptr);
			futures.emplace_back(taskptr->get_future());
			asio::post(pool, [taskptr, Q, h_idx, k_idx]()
			{
				(*taskptr)(Q[0], Q[1], Q[2], h_idx, k_idx);
			});

			++task_idx;
			m_progress->setValue(task_idx);
		}
	}

	m_progress->setValue(0);
	m_status->setText("Calculating grid.");

	std::deque<std::uint64_t> hklindices;
	if(format == EXPORT_GRID)  // Takin grid format
	{
		std::uint64_t dummy = 0;  // to be filled by index block index
		ofstr->write(reinterpret_cast<const char*>(&dummy), sizeof(dummy));

		for(int i = 0; i<3; ++i)
		{
			ofstr->write(reinterpret_cast<const char*>(&Qstart[i]), sizeof(Qstart[i]));
			ofstr->write(reinterpret_cast<const char*>(&Qend[i]), sizeof(Qend[i]));
			ofstr->write(reinterpret_cast<const char*>(&Qstep[i]), sizeof(Qstep[i]));
		}

		(*ofstr) << "Takin/Magdyn Grid File Version 2 (doi: https://doi.org/10.5281/zenodo.4117437).";
	}

#ifdef USE_HDF5
	std::vector<t_real> data_energies, data_weights;
	std::vector<std::size_t> data_indices, data_num_branches;
	data_energies.reserve(num_pts_h * num_pts_k * num_pts_l * 32);
	data_weights.reserve(num_pts_h * num_pts_k * num_pts_l * 32);
	data_indices.reserve(num_pts_h * num_pts_k * num_pts_l);
	data_num_branches.reserve(num_pts_h * num_pts_k * num_pts_l);
#endif

	for(std::size_t future_idx=0; future_idx<futures.size(); ++future_idx)
	{
		qApp->processEvents();  // process events to see if the stop button was clicked
		if(m_stopRequested)
		{
			pool.stop();
			break;
		}

		auto results = futures[future_idx].get();

		for(std::size_t result_idx=0; result_idx<results.size(); ++result_idx)
		{
			const auto& result = results[result_idx];
			std::size_t num_branches = std::get<3>(result).size();

			const t_real qh = std::get<0>(result);
			const t_real qk = std::get<1>(result);
			const t_real ql = std::get<2>(result);
			const auto& energies = std::get<3>(result);
			const auto& weights = std::get<4>(result);

			if(format == EXPORT_GRID)           // Takin grid format
			{
				hklindices.push_back(ofstr->tellp());

				// write number of branches
				std::uint32_t _num_branches = std::uint32_t(num_branches);
				ofstr->write(reinterpret_cast<const char*>(&_num_branches), sizeof(_num_branches));
			}
			else if(format == EXPORT_TEXT)      // text format
			{
				(*ofstr) << "Q = " << qh << " " << qk << " " << ql << ":\n";
			}
#ifdef USE_HDF5
			else if(format == EXPORT_HDF5)  // hdf5 format
			{
#ifdef WRITE_HDF5_CHUNKS
				const std::size_t h_idx = std::get<5>(result);
				const std::size_t k_idx = std::get<6>(result);
				const std::size_t l_idx = std::get<7>(result);

				std::ostringstream chunk_name;
				chunk_name << std::hex << h_idx << "_" << k_idx << "_" << l_idx;
				H5::Group data_group = h5file->openGroup("chunks");
				data_group.createGroup(chunk_name.str());

				tl2::set_h5_scalar(*h5file, "chunks/" + chunk_name.str() + "/h", qh);
				tl2::set_h5_scalar(*h5file, "chunks/" + chunk_name.str() + "/k", qk);
				tl2::set_h5_scalar(*h5file, "chunks/" + chunk_name.str() + "/l", ql);

				tl2::set_h5_vector(*h5file, "chunks/" + chunk_name.str() + "/E", energies);
				tl2::set_h5_vector(*h5file, "chunks/" + chunk_name.str() + "/S", weights);
#endif

				// TODO: write this incrementally to the hdf5 file, without these large temporary vectors
				data_num_branches.push_back(energies.size());
				data_indices.push_back(data_energies.size());
				data_energies.insert(data_energies.end(), energies.begin(), energies.end());
				data_weights.insert(data_weights.end(), weights.begin(), weights.end());
			}
#endif

			// iterate energies and weights
			for(std::size_t j=0; j<num_branches; ++j)
			{
				t_real energy = energies[j];
				t_real weight = weights[j];

				if(format == EXPORT_GRID)       // Takin grid format
				{
					// write energies and weights
					ofstr->write(reinterpret_cast<const char*>(&energy), sizeof(energy));
					ofstr->write(reinterpret_cast<const char*>(&weight), sizeof(weight));
				}
				else if(format == EXPORT_TEXT)  // text format
				{
					(*ofstr)
						<< "\tE = " << energy
						<< ", S = " << weight
						<< std::endl;
				}
			}
		}

		m_progress->setValue(future_idx+1);
	}

	pool.join();
	EnableInput();

	if(format == EXPORT_GRID)  // Takin grid format
	{
		std::uint64_t idxblock = ofstr->tellp();

		// write hkl indices
		for(std::uint64_t idx : hklindices)
			ofstr->write(reinterpret_cast<const char*>(&idx), sizeof(idx));

		// write index into index block
		ofstr->seekp(0, std::ios_base::beg);
		ofstr->write(reinterpret_cast<const char*>(&idxblock), sizeof(idxblock));
		ofstr->flush();
	}
#ifdef USE_HDF5
	else if(format == EXPORT_HDF5)
	{
		const char* user = std::getenv("USER");
		if(!user)
			user = "";
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/type", "takin_grid");
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/description", "Takin/Magdyn grid format");
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/user", user);
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/date", tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()));
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/url", "https://github.com/ILLGrenoble/takin");
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/doi", "https://doi.org/10.5281/zenodo.4117437");
		tl2::set_h5_string<std::string>(*h5file, "meta_infos/doi_tlibs", "https://doi.org/10.5281/zenodo.5717779");

		tl2::set_h5_string<std::string>(*h5file, "infos/shape", "cuboid");
		tl2::set_h5_vector(*h5file, "infos/Q_start", static_cast<const std::vector<t_real>&>(Qstart));
		tl2::set_h5_vector(*h5file, "infos/Q_end", static_cast<const std::vector<t_real>&>(Qend));
		tl2::set_h5_vector(*h5file, "infos/Q_steps", static_cast<const std::vector<t_real>&>(Qstep));
		tl2::set_h5_vector(*h5file, "infos/Q_dimensions", std::vector<std::size_t>{{ num_pts_h, num_pts_k, num_pts_l }});

		std::vector<std::string> labels{{"h", "k", "l", "E", "S_perp"}};
		tl2::set_h5_string_vector(*h5file, "infos/labels", labels);

		std::vector<std::string> units{{"rlu", "rlu", "rlu", "meV", "a.u."}};
		tl2::set_h5_string_vector(*h5file, "infos/units", units);

		hsize_t index_dims[] = { num_pts_h, num_pts_k, num_pts_l };
		tl2::set_h5_multidim(*h5file, "data/indices", 3, index_dims, data_indices.data());
		tl2::set_h5_multidim(*h5file, "data/branches", 3, index_dims, data_num_branches.data());
		//tl2::set_h5_vector(*h5file, "data/indices", data_indices);
		//tl2::set_h5_vector(*h5file, "data/branches", data_num_branches);

		tl2::set_h5_vector(*h5file, "data/energies", data_energies);
		tl2::set_h5_vector(*h5file, "data/weights", data_weights);

		h5file->close();
	}
#endif

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

	return true;
}
