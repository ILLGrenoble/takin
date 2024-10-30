/**
 * magnetic dynamics -- saving of dispersion data
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

#include "magdyn.h"

#include <QtCore/QString>
#include <QtWidgets/QApplication>

#include <boost/scope_exit.hpp>

#include <vector>

#include "tlibs2/libs/algos.h"


// --------------------------------------------------------------------------------
/**
 * save the plot as pdf
 */
void MagDynDlg::SavePlotFigure()
{
	if(!m_plot)
		return;

	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	m_plot->savePdf(filename);
}



/**
 * save the data for a dispersion direction
 */
void MagDynDlg::SaveDispersion(bool as_scr)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableInput();
	} BOOST_SCOPE_EXIT_END
	DisableInput();
	m_stopRequested = false;

	QString dirLast = m_sett->value("dir", "").toString();

	QString filename;
	if(as_scr)
	{
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data As Script",
			dirLast, "Py Files (*.py)");
	}
	else
	{
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data",
			dirLast, "Data Files (*.dat)");
	}

	if(filename == "")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	const t_real Q_start[]
	{
		(t_real)m_Q_start[0]->value(),
		(t_real)m_Q_start[1]->value(),
		(t_real)m_Q_start[2]->value(),
	};

	const t_real Q_end[]
	{
		(t_real)m_Q_end[0]->value(),
		(t_real)m_Q_end[1]->value(),
		(t_real)m_Q_end[2]->value(),
	};

	const t_size num_pts = m_num_points->value();

	// function to track progress and request stopping
	std::function<bool(int, int)> progress_fkt =
		[this](int progress, int total) -> bool
	{
		if(total >= 0)
			m_progress->setMaximum(total);
		if(progress >= 0)
			m_progress->setValue(progress);

		qApp->processEvents();
		return !m_stopRequested;
	};

	// start calculation
	m_status->setText("Calculating dispersion.");
	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	bool ok = m_dyn.SaveDispersion(filename.toStdString(),
		Q_start[0], Q_start[1], Q_start[2],
		Q_end[0], Q_end[1], Q_end[2],
		num_pts, g_num_threads, as_scr,
		&progress_fkt);

	// print timing information
        stopwatch.stop();
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested)
		ostrMsg << " stopped ";
	else if(ok)
		ostrMsg << " finished ";
	else
		ostrMsg << " failed ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());
}



/**
 * save the data for multiple dispersion directions
 * given by the saved coordinates
 */
void MagDynDlg::SaveMultiDispersion(bool as_scr)
{
	using t_item = tl2::NumericTableWidgetItem<t_real>;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableInput();
	} BOOST_SCOPE_EXIT_END
	DisableInput();
	m_stopRequested = false;

	// get all Qs from the coordinates table
	std::vector<t_vec_real> Qs;
	std::vector<std::string> Q_names;
	Qs.reserve(m_coordinatestab->rowCount());
	Q_names.reserve(m_coordinatestab->rowCount());

	for(int coord_row = 0; coord_row < m_coordinatestab->rowCount(); ++coord_row)
	{
		std::string name = m_coordinatestab->item(
			coord_row, COL_COORD_NAME)->text().toStdString();
		t_real h = static_cast<t_item*>(
			m_coordinatestab->item(coord_row, COL_COORD_H))->GetValue();
		t_real k = static_cast<t_item*>(
			m_coordinatestab->item(coord_row, COL_COORD_K))->GetValue();
		t_real l = static_cast<t_item*>(
			m_coordinatestab->item(coord_row, COL_COORD_L))->GetValue();

		Q_names.emplace_back(std::move(name));
		Qs.emplace_back(tl2::create<t_vec_real>({ h, k, l }));
	}

	if(Qs.size() <= 1)
	{
		QMessageBox::critical(this, "Magnetic Dynamics",
			"Not enough Q coordinates available, "
			"please define them in the \"Coordinates\" tab.");
		return;
	}

	// get file name to save dispersions to
	QString dirLast = m_sett->value("dir", "").toString();
	QString filename;
	if(as_scr)
	{
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data As Script",
			dirLast, "Py Files (*.py)");
	}
	else
	{
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data",
			dirLast, "Data Files (*.dat)");
	}
	if(filename == "")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	const t_size num_pts = m_num_points->value();
	const t_size num_disps = Qs.size() - 1;
	t_size cur_disp = 1;
	int prev_progress = -1;
	int total_progress = -1;

	// function to track progress and request stopping
	std::function<bool(int, int)> progress_fkt =
		[this, num_disps, &cur_disp, &prev_progress, &total_progress](
			int progress, int total) -> bool
	{
		if(progress == 0 && prev_progress == total)
		{
			++cur_disp;
			m_status->setText(QString(
				"Calculating dispersion %1/%2.").arg(cur_disp).arg(num_disps));
		}

		if(total >= 0)
		{
			m_progress->setMaximum(total);
			total_progress = total;
		}
		if(progress >= 0)
		{
			m_progress->setValue(progress);
			prev_progress = progress;
		}

		qApp->processEvents();
		return !m_stopRequested;
	};

	// start calculation
	m_status->setText(QString("Calculating dispersion 1/%1.").arg(num_disps));
	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	bool ok = m_dyn.SaveMultiDispersion(filename.toStdString(),
		Qs, num_pts, g_num_threads, as_scr,
		&progress_fkt, &Q_names);

	// print timing information
        stopwatch.stop();
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested)
		ostrMsg << " stopped ";
	else if(ok)
		ostrMsg << " finished ";
	else
		ostrMsg << " failed ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());
}
// --------------------------------------------------------------------------------
