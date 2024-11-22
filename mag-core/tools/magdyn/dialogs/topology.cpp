/**
 * magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
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

#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <limits>
#include <mutex>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QMessageBox>

#include "topology.h"
#include "helper.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



// ============================================================================
// topology dialog
// ============================================================================

/**
 * sets up the topology dialog
 */
TopologyDlg::TopologyDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Topology");
	setSizeGripEnabled(true);

	// tab widget
	m_tabs = new QTabWidget(this);

	// status bar
	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	// close button
	QPushButton *btnOk = new QPushButton("OK", this);

	// main grid
	QGridLayout *maingrid = new QGridLayout(this);
	maingrid->setSpacing(4);
	maingrid->setContentsMargins(8, 8, 8, 8);
	maingrid->addWidget(m_tabs, 0, 0, 1, 4);
	maingrid->addWidget(m_status, 1, 0, 1, 3);
	maingrid->addWidget(btnOk, 1, 3, 1, 1);

	// tab panels
	m_tabs->addTab(CreateBerryCurvaturePanel(), "Berry Curvature");

	// connections
	connect(btnOk, &QAbstractButton::clicked, this, &QDialog::accept);
}



TopologyDlg::~TopologyDlg()
{
}



/**
 * set a pointer to the main magdyn kernel
 */
void TopologyDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;
}



/**
 * set the Q start and end points from the main window's dispersion
 */
void TopologyDlg::SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend)
{
	m_Qstart = Qstart;
	m_Qend = Qend;
}



void TopologyDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * dialog is closing
 */
void TopologyDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("topology/geo", saveGeometry());
		m_sett->setValue("topology/splitter", m_split_plot->saveState());
	}

	QDialog::accept();
}
// ============================================================================



// ============================================================================
// calculate berry curvature
// ============================================================================

/**
 * column indices in magnon band table for the berry curvature
 */
enum : int
{
	COL_BC_BAND = 0,
	NUM_COLS_BC,
};



/**
 * create the panel for the berry curvature tab
 */
QWidget* TopologyDlg::CreateBerryCurvaturePanel()
{
	QWidget *panelBerryCurvature = new QWidget(this);

	// plotter
	m_plot_bc = new QCustomPlot(panelBerryCurvature);
	m_plot_bc->setFont(font());
	m_plot_bc->xAxis->setLabel("Momentum Transfer Q (rlu)");
	m_plot_bc->yAxis->setLabel("Berry Curvature B");
	m_plot_bc->setInteraction(QCP::iRangeDrag, true);
	m_plot_bc->setInteraction(QCP::iRangeZoom, true);
	m_plot_bc->setSelectionRectMode(QCP::srmZoom);
	m_plot_bc->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// magnon band table
	m_tableBands = new QTableWidget(panelBerryCurvature);
	m_tableBands->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_tableBands->setShowGrid(true);
	m_tableBands->setSortingEnabled(false);
	m_tableBands->setSelectionBehavior(QTableWidget::SelectRows);
	m_tableBands->setSelectionMode(QTableWidget::SingleSelection);
	m_tableBands->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_tableBands->verticalHeader()->setVisible(false);
	m_tableBands->setColumnCount(NUM_COLS_BC);

	m_tableBands->setHorizontalHeaderItem(COL_BC_BAND, new QTableWidgetItem{"Band"});
	m_tableBands->setColumnWidth(COL_BC_BAND, 70);

	// splitter for plot and magnon band list
	m_split_plot = new QSplitter(panelBerryCurvature);
	m_split_plot->setOrientation(Qt::Horizontal);
	m_split_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot->addWidget(m_plot_bc);
	m_split_plot->addWidget(m_tableBands);
	m_split_plot->setCollapsible(0, false);
	m_split_plot->setCollapsible(1, true);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(m_plot_bc), 24);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(m_tableBands), 1);

	// context menu for plotter
	m_menuPlot_bc = new QMenu("Plotter", panelBerryCurvature);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot_bc);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuPlot_bc);
	QAction *acSaveData = new QAction("Save Data...", m_menuPlot_bc);

	acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));

	m_E_positive = new QAction("Ignore Magnon Annihilation", m_menuPlot_bc);
	m_E_positive->setCheckable(true);
	m_E_positive->setChecked(false);

	m_imag_bc = new QAction("Show Imaginary B Component", m_menuPlot_bc);
	m_imag_bc->setCheckable(true);
	m_imag_bc->setChecked(false);
	m_imag_bc->setToolTip("Show the imaginary component of the Berry curvature.");

	m_menuPlot_bc->addAction(acRescalePlot);
	m_menuPlot_bc->addSeparator();
	m_menuPlot_bc->addAction(acSaveFigure);
	m_menuPlot_bc->addAction(acSaveData);
	m_menuPlot_bc->addSeparator();
	m_menuPlot_bc->addAction(m_E_positive);
	m_menuPlot_bc->addAction(m_imag_bc);

	// start and stop coordinates
	m_Q_start_bc[0] = new QDoubleSpinBox(panelBerryCurvature);
	m_Q_start_bc[1] = new QDoubleSpinBox(panelBerryCurvature);
	m_Q_start_bc[2] = new QDoubleSpinBox(panelBerryCurvature);
	m_Q_end_bc[0] = new QDoubleSpinBox(panelBerryCurvature);
	m_Q_end_bc[1] = new QDoubleSpinBox(panelBerryCurvature);
	m_Q_end_bc[2] = new QDoubleSpinBox(panelBerryCurvature);

	m_Q_start_bc[0]->setToolTip("Dispersion initial momentum transfer, h_i (rlu).");
	m_Q_start_bc[1]->setToolTip("Dispersion initial momentum transfer, k_i (rlu).");
	m_Q_start_bc[2]->setToolTip("Dispersion initial momentum transfer, l_i (rlu).");
	m_Q_end_bc[0]->setToolTip("Dispersion final momentum transfer, h_f (rlu).");
	m_Q_end_bc[1]->setToolTip("Dispersion final momentum transfer, k_f (rlu).");
	m_Q_end_bc[2]->setToolTip("Dispersion final momentum transfer, l_f (rlu).");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_start_bc[i]->setDecimals(4);
		m_Q_start_bc[i]->setMinimum(-99.9999);
		m_Q_start_bc[i]->setMaximum(+99.9999);
		m_Q_start_bc[i]->setSingleStep(0.01);
		m_Q_start_bc[i]->setValue(i == 0 ? -1. : 0.);
		m_Q_start_bc[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q_start_bc[i]->setPrefix(hklPrefix[i]);

		m_Q_end_bc[i]->setDecimals(4);
		m_Q_end_bc[i]->setMinimum(-99.9999);
		m_Q_end_bc[i]->setMaximum(+99.9999);
		m_Q_end_bc[i]->setSingleStep(0.01);
		m_Q_end_bc[i]->setValue(i == 0 ? 1. : 0.);
		m_Q_end_bc[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q_end_bc[i]->setPrefix(hklPrefix[i]);
	}

	// number of Q points in the plot
	m_num_Q_bc = new QSpinBox(panelBerryCurvature);
	m_num_Q_bc->setMinimum(1);
	m_num_Q_bc->setMaximum(99999);
	m_num_Q_bc->setValue(128);
	m_num_Q_bc->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_num_Q_bc->setToolTip("Number of Q points to calculate.");

	// dispersion Q button
	QPushButton *btnQ = new QPushButton("Set Main Q", panelBerryCurvature);
	btnQ->setToolTip("Set the Q start and end points from the dispersion in the main window.");

	// coordinate components
	m_coords_bc[0] = new QSpinBox(panelBerryCurvature);
	m_coords_bc[0]->setMinimum(0);
	m_coords_bc[0]->setMaximum(2);
	m_coords_bc[0]->setValue(0);
	m_coords_bc[0]->setPrefix("i = ");
	m_coords_bc[0]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_coords_bc[0]->setToolTip("First component index of B_ij matrix.");

	m_coords_bc[1] = new QSpinBox(panelBerryCurvature);
	m_coords_bc[1]->setMinimum(0);
	m_coords_bc[1]->setMaximum(2);
	m_coords_bc[1]->setValue(1);
	m_coords_bc[1]->setPrefix("j = ");
	m_coords_bc[1]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_coords_bc[1]->setToolTip("Second component index of B_ij matrix.");

	// maximum cutoff for filtering numerical artefacts in berry curvature
	m_B_filter_enable_bc = new QCheckBox("Maximum B:", panelBerryCurvature);
	m_B_filter_enable_bc->setChecked(true);
	m_B_filter_enable_bc->setToolTip("Enable maximum cutoff Berry curvature for filtering numerical artefacts.");

	m_B_filter_bc = new QDoubleSpinBox(panelBerryCurvature);
	m_B_filter_bc->setDecimals(2);
	m_B_filter_bc->setMinimum(0.);
	m_B_filter_bc->setMaximum(999999.99);
	m_B_filter_bc->setSingleStep(1.);
	m_B_filter_bc->setValue(m_B_filter_bc->maximum());
	m_B_filter_bc->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_B_filter_bc->setToolTip("Maximum cutoff Berry curvature for filtering numerical artefacts.");

	// minimum cutoff for filtering S(Q, E)
	m_S_filter_enable_bc = new QCheckBox("Minimum S(Q, E):", panelBerryCurvature);
	m_S_filter_enable_bc->setChecked(false);
	m_S_filter_enable_bc->setToolTip("Enable minimum S(Q, E).");

	m_S_filter_bc = new QDoubleSpinBox(panelBerryCurvature);
	m_S_filter_bc->setDecimals(5);
	m_S_filter_bc->setMinimum(0.);
	m_S_filter_bc->setMaximum(9999.99999);
	m_S_filter_bc->setSingleStep(0.01);
	m_S_filter_bc->setValue(0.01);
	m_S_filter_bc->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_S_filter_bc->setToolTip("Minimum S(Q, E) to keep.");

	// progress bar
	m_progress_bc = new QProgressBar(panelBerryCurvature);
	m_progress_bc->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	// start/stop button
	m_btnStartStop_bc = new QPushButton("Calculate", panelBerryCurvature);

	// component grid
	auto grid = new QGridLayout(panelBerryCurvature);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_split_plot, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Start Q (rlu):", panelBerryCurvature), y, 0, 1, 1);
	grid->addWidget(m_Q_start_bc[0], y, 1, 1, 1);
	grid->addWidget(m_Q_start_bc[1], y, 2, 1, 1);
	grid->addWidget(m_Q_start_bc[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("End Q (rlu):", panelBerryCurvature), y, 0, 1, 1);
	grid->addWidget(m_Q_end_bc[0], y, 1, 1, 1);
	grid->addWidget(m_Q_end_bc[1], y, 2, 1, 1);
	grid->addWidget(m_Q_end_bc[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Count:", panelBerryCurvature), y, 0, 1, 1);
	grid->addWidget(m_num_Q_bc, y, 1, 1, 1);
	grid->addWidget(btnQ, y++, 3, 1, 1);
	grid->addWidget(new QLabel("B Component:", panelBerryCurvature), y, 0, 1, 1);
	grid->addWidget(m_coords_bc[0], y, 1, 1, 1);
	grid->addWidget(m_coords_bc[1], y++, 2, 1, 1);
	grid->addWidget(m_B_filter_enable_bc, y, 0, 1, 1);
	grid->addWidget(m_B_filter_bc, y, 1, 1, 1);
	grid->addWidget(m_S_filter_enable_bc, y, 2, 1, 1);
	grid->addWidget(m_S_filter_bc, y++, 3, 1, 1);
	grid->addWidget(m_progress_bc, y, 0, 1, 3);
	grid->addWidget(m_btnStartStop_bc, y++, 3, 1, 1);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("topology/geo"))
			restoreGeometry(m_sett->value("topology/geo").toByteArray());
		else
			resize(640, 640);

		if(m_sett->contains("topology/splitter"))
			m_split_plot->restoreState(m_sett->value("topology/splitter").toByteArray());
	}

	// connections
	connect(m_plot_bc, &QCustomPlot::mouseMove, this, &TopologyDlg::BerryCurvaturePlotMouseMove);
	connect(m_plot_bc, &QCustomPlot::mousePress, this, &TopologyDlg::BerryCurvaturePlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &TopologyDlg::RescaleBerryCurvaturePlot);
	connect(acSaveFigure, &QAction::triggered, this, &TopologyDlg::SaveBerryCurvaturePlotFigure);
	connect(acSaveData, &QAction::triggered, this, &TopologyDlg::SaveBerryCurvatureData);
	connect(btnQ, &QAbstractButton::clicked, this, &TopologyDlg::SetBerryCurvatureQ);
	connect(m_B_filter_enable_bc, &QCheckBox::toggled, m_B_filter_bc, &QDoubleSpinBox::setEnabled);
	connect(m_S_filter_enable_bc, &QCheckBox::toggled, m_S_filter_bc, &QDoubleSpinBox::setEnabled);

	// calculation
	connect(m_btnStartStop_bc, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled_bc)
			CalculateBerryCurvature();
		else
			m_stopRequested_bc = true;
	});

	// replotting
	connect(m_E_positive, &QAction::toggled, this, &TopologyDlg::PlotBerryCurvature);
	connect(m_imag_bc, &QAction::toggled, this, &TopologyDlg::PlotBerryCurvature);
	connect(m_B_filter_enable_bc, &QCheckBox::toggled, this, &TopologyDlg::PlotBerryCurvature);
	connect(m_S_filter_enable_bc, &QCheckBox::toggled, this, &TopologyDlg::PlotBerryCurvature);
	connect(m_B_filter_bc, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		this, &TopologyDlg::PlotBerryCurvature);
	connect(m_S_filter_bc, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		this, &TopologyDlg::PlotBerryCurvature);

	m_B_filter_bc->setEnabled(m_B_filter_enable_bc->isChecked());
	m_S_filter_bc->setEnabled(m_S_filter_enable_bc->isChecked());
	EnableBerryCurvatureCalculation();

	return panelBerryCurvature;
}



/**
 * clears the table of magnon bands
 */
void TopologyDlg::ClearBerryCurvatureBands()
{
	m_tableBands->clearContents();
	m_tableBands->setRowCount(0);
}



/**
 * adds a magnon band to the table
 */
void TopologyDlg::AddBerryCurvatureBand(const std::string& name, const QColor& colour)
{
	int row = m_tableBands->rowCount();
	m_tableBands->insertRow(row);

	QTableWidgetItem *item = new QTableWidgetItem{name.c_str()};
	item->setFlags(item->flags() & ~Qt::ItemIsEditable);

	QBrush bg = item->background();
	bg.setColor(colour);
	bg.setStyle(Qt::SolidPattern);
	item->setBackground(bg);

	m_tableBands->setItem(row, COL_BC_BAND, item);
}



/**
 * calculate the filtered data sets and plot the berry curvature
 */
void TopologyDlg::PlotBerryCurvature()
{
	if(!m_plot_bc)
		return;

	ClearBerryCurvatureBands();
	ClearBerryCurvaturePlot(false);

	if(m_data_bc.size() == 0)
	{
		m_plot_bc->replot();
		return;
	}

	// get settings
	t_real max_B = m_B_filter_bc->value();
	if(!m_B_filter_enable_bc->isChecked())
		max_B = -1.;  // disable B filter

	t_real min_S = m_S_filter_bc->value();
	if(!m_S_filter_enable_bc->isChecked())
		min_S = -1.;  // disable S(Q, E) filter

	bool show_imag_comp = m_imag_bc->isChecked();
	bool only_creation = m_E_positive->isChecked();

	t_size num_Q = m_data_bc.size();
	t_size num_bands = m_data_bc[0].curvatures.size();

	// filtered momentum transfer and berry curvature per band
	std::vector<QVector<qreal>> Qs_data_bc{num_bands};
	std::vector<QVector<qreal>> Bs_data_bc{num_bands};

	for(t_size Q_idx = 0; Q_idx < num_Q; ++Q_idx)
	{
		const t_vec_real& Q = m_data_bc[Q_idx].momentum;

		for(t_size band = 0; band < num_bands; ++band)
		{
			t_real berry_comp = show_imag_comp
				? m_data_bc[Q_idx].curvatures[band].imag()
				: m_data_bc[Q_idx].curvatures[band].real();

			// filter numerical artefacts in B
			if(max_B >= 0. && std::abs(berry_comp) > max_B)
				continue;

			// filter minimum S(Q, E)
			if(min_S >= 0. && std::abs(m_data_bc[Q_idx].weights[band]) <= min_S)
				continue;

			// filter magnon annihilation
			if(only_creation && m_data_bc[Q_idx].energies[band] < 0.)
				continue;

			Qs_data_bc[band].push_back(Q[m_Q_idx_bc]);
			Bs_data_bc[band].push_back(berry_comp);
		}
	}

	// sort filtered data by Q
	auto sort_data = [](QVector<qreal>& Qvec, QVector<qreal>& Bvec)
	{
		// sort vectors by Q component
		std::vector<std::size_t> perm = tl2::get_perm(Qvec.size(),
			[&Qvec](std::size_t idx1, std::size_t idx2) -> bool
		{
			return Qvec[idx1] < Qvec[idx2];
		});

		Qvec = tl2::reorder(Qvec, perm);
		Bvec = tl2::reorder(Bvec, perm);
	};

	for(t_size band = 0; band < Bs_data_bc.size(); ++band)
		sort_data(Qs_data_bc[band], Bs_data_bc[band]);

	// berry curvature range
	t_real B_min_bc = std::numeric_limits<t_real>::max();
	t_real B_max_bc = -B_min_bc;

	// get berry curvature range
	for(const QVector<t_real>& Bs_data : Bs_data_bc)
	{
		auto [min_B_iter, max_B_iter] = std::minmax_element(Bs_data.begin(), Bs_data.end());
		if(min_B_iter != Bs_data.end() && max_B_iter != Bs_data.end())
		{
			t_real B_range = *max_B_iter - *min_B_iter;

			B_max_bc = std::max(B_max_bc, *max_B_iter + B_range*0.05);
			B_min_bc = std::min(B_min_bc, *min_B_iter - B_range*0.05);
		}
	}

	// how many bands do actually have data?
	t_size num_effective_bands = 0;
	for(t_size band = 0; band < num_bands; ++band)
	{
		if(Bs_data_bc[band].size() != 0)
			++num_effective_bands;
	}

	// plot berry curvatures per band
	t_size effective_band = 0;
	for(t_size band = 0; band < num_bands; ++band)
	{
		// ignore bands with no data
		if(Bs_data_bc[band].size() == 0)
			continue;

		QCPCurve *curve = new QCPCurve(m_plot_bc->xAxis, m_plot_bc->yAxis);

		// colour for this magnon band
		QPen pen = curve->pen();
		int col[3] = {
			int(std::lerp(1., 0., t_real(effective_band) / t_real(num_effective_bands - 1)) * 255.),
			0x00,
			int(std::lerp(0., 1., t_real(effective_band) / t_real(num_effective_bands - 1)) * 255.),
		};

		//get_colour<int>(g_colPlot, col);
		const QColor colFull(col[0], col[1], col[2]);
		pen.setColor(colFull);
		pen.setWidthF(2.);

		curve->setPen(pen);
		curve->setLineStyle(QCPCurve::lsLine);
		curve->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 1));
		curve->setAntialiased(true);
		curve->setData(Qs_data_bc[band], Bs_data_bc[band]);

		m_curves_bc.push_back(curve);
		AddBerryCurvatureBand("#" + tl2::var_to_str(effective_band + 1), colFull);
		++effective_band;
	}

	// set labels
	const char* Q_label[]{ "h (rlu)", "k (rlu)", "l (rlu)" };
	m_plot_bc->xAxis->setLabel(QString("Momentum Transfer ") + Q_label[m_Q_idx_bc]);

	// set ranges
	m_plot_bc->xAxis->setRange(m_Q_min_bc, m_Q_max_bc);
	m_plot_bc->yAxis->setRange(B_min_bc, B_max_bc);

	// set font
	m_plot_bc->setFont(font());
	m_plot_bc->xAxis->setLabelFont(font());
	m_plot_bc->yAxis->setLabelFont(font());
	m_plot_bc->xAxis->setTickLabelFont(font());
	m_plot_bc->yAxis->setTickLabelFont(font());

	m_plot_bc->replot();
}



/**
 * calculate the berry curvature
 */
void TopologyDlg::CalculateBerryCurvature()
{
	if(!m_dyn)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableBerryCurvatureCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableBerryCurvatureCalculation(false);

	ClearBerryCurvaturePlot(false);

	// get coordinates
	t_vec_real Q_start = tl2::create<t_vec_real>(
	{
		m_Q_start_bc[0]->value(),
		m_Q_start_bc[1]->value(),
		m_Q_start_bc[2]->value(),
	});

	t_vec_real Q_end = tl2::create<t_vec_real>(
	{
		m_Q_end_bc[0]->value(),
		m_Q_end_bc[1]->value(),
		m_Q_end_bc[2]->value(),
	});

	// get Q component with maximum range
	t_vec_real Q_range = Q_end - Q_start;
	m_Q_idx_bc = 0;
	if(std::abs(Q_range[1]) > std::abs(Q_range[m_Q_idx_bc]))
		m_Q_idx_bc = 1;
	if(std::abs(Q_range[2]) > std::abs(Q_range[m_Q_idx_bc]))
		m_Q_idx_bc = 2;

	// keep the scanned Q component in ascending order
	if(Q_start[m_Q_idx_bc] > Q_end[m_Q_idx_bc])
		std::swap(Q_start, Q_end);

	// Q range
	m_Q_min_bc = Q_start[m_Q_idx_bc];
	m_Q_max_bc = Q_end[m_Q_idx_bc];

	// get settings
	t_size Q_count = m_num_Q_bc->value();
	std::vector<t_size> *perm = nullptr;
	t_size dim1 = m_coords_bc[0]->value();
	t_size dim2 = m_coords_bc[1]->value();

	// calculate berry curvature
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(false);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested_bc = false;
	m_progress_bc->setMinimum(0);
	m_progress_bc->setMaximum(Q_count);
	m_progress_bc->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count);

	m_data_bc.clear();
	m_data_bc.reserve(Q_count);

	for(t_size Q_idx = 0; Q_idx < Q_count; ++Q_idx)
	{
		auto task = [this, &mtx, &dyn, &Q_start, &Q_end, Q_idx, Q_count, perm, dim1, dim2]()
		{
			const t_vec_real Q = Q_count > 1
				? tl2::lerp(Q_start, Q_end, t_real(Q_idx) / t_real(Q_count - 1))
				: Q_start;

			// calculate berry curvatures per band
			BerryCurvatureData data_bc;
			data_bc.momentum = Q;
			typename t_magdyn::SofQE S;
			std::tie(data_bc.curvatures, S) = dyn.CalcBerryCurvatures(
				Q, g_delta_diff, perm, dim1, dim2, g_evecs_ortho != 0);
			t_size num_bands = data_bc.curvatures.size();
			data_bc.energies.reserve(num_bands);
			data_bc.weights.reserve(num_bands);

			// calculate energies per band
			assert(S.E_and_S.size() == num_bands);
			for(t_size band = 0; band < num_bands; ++band)
			{
				data_bc.energies.push_back(S.E_and_S[band].E);
				data_bc.weights.push_back(S.E_and_S[band].weight);
			}

			std::lock_guard<std::mutex> _lck{mtx};
			m_data_bc.emplace_back(std::move(data_bc));
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText(QString("Calculating in %1 threads...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		qApp->processEvents();  // process events to see if the stop button was clicked
		if(m_stopRequested_bc)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress_bc->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested_bc)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	// sort raw unfiltered data by Q
	std::vector<std::size_t> perm_all = tl2::get_perm(m_data_bc.size(),
		[this](std::size_t idx1, std::size_t idx2) -> bool
	{
		return m_data_bc[idx1].momentum[m_Q_idx_bc]
			< m_data_bc[idx2].momentum[m_Q_idx_bc];
	});

	m_data_bc = tl2::reorder(m_data_bc, perm_all);

	PlotBerryCurvature();
}



/**
 * clears the dispersion graph
 */
void TopologyDlg::ClearBerryCurvaturePlot(bool replot)
{
	m_curves_bc.clear();

	if(m_plot_bc)
	{
		m_plot_bc->clearPlottables();
		if(replot)
			m_plot_bc->replot();
	}
}



/**
 * show current cursor coordinates
 */
void TopologyDlg::BerryCurvaturePlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot_bc->xAxis->pixelToCoord(evt->pos().x());
	t_real berry = m_plot_bc->yAxis->pixelToCoord(evt->pos().y());

	QString status("Q = %1 rlu, B = %2.");
	status = status.arg(Q, 0, 'g', g_prec_gui).arg(berry, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



/**
 * show plot context menu
 */
void TopologyDlg::BerryCurvaturePlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot_bc)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot_bc->popup(pos);
		evt->accept();
	}
}



/**
 * rescale plot axes to fit the content
 */
void TopologyDlg::RescaleBerryCurvaturePlot()
{
	if(!m_plot_bc)
		return;

	m_plot_bc->rescaleAxes();
	m_plot_bc->replot();
}



/**
 * save plot as image file
 */
void TopologyDlg::SaveBerryCurvaturePlotFigure()
{
	if(!m_plot_bc)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("topology/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("topology/dir", QFileInfo(filename).path());

	if(!m_plot_bc->savePdf(filename))
		ShowError(QString("Could not save figure to file \"%1\".").arg(filename).toStdString().c_str());
}



/**
 * save plot as data file
 */
void TopologyDlg::SaveBerryCurvatureData()
{
	if(m_data_bc.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("topology/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Data", dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("topology/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".").arg(filename).toStdString().c_str());
		return;
	}

	t_size num_bands = m_data_bc[0].curvatures.size();

	ofstr.precision(g_prec);
	int field_len = g_prec * 2.5;

	// write meta header
	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr << "#\n"
		<< "# Created by Takin/Magdyn\n"
		<< "# URL: https://github.com/ILLGrenoble/takin\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.4117437\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n# Number of energy bands: " << num_bands << "\n"
		<< "#\n\n";

	// write column header
	ofstr << std::setw(field_len) << std::left << "# h" << " ";
	ofstr << std::setw(field_len) << std::left << "k" << " ";
	ofstr << std::setw(field_len) << std::left << "l" << " ";

	for(t_size band = 0; band < num_bands; ++band)
	{
		std::string E = "E_" + tl2::var_to_str(band);
		std::string S = "Sperp_" + tl2::var_to_str(band);
		std::string ReB = "Re{B_" + tl2::var_to_str(band) + "}";
		std::string ImB = "Im{B_" + tl2::var_to_str(band) + "}";

		ofstr << std::setw(field_len) << std::left << E << " ";
		ofstr << std::setw(field_len) << std::left << S << " ";
		ofstr << std::setw(field_len) << std::left << ReB << " ";
		ofstr << std::setw(field_len) << std::left << ImB << " ";
	}
	ofstr << "\n";

	// write data
	for(const BerryCurvatureData& data: m_data_bc)
	{
		ofstr << std::setw(field_len) << std::left << data.momentum[0] << " ";
		ofstr << std::setw(field_len) << std::left << data.momentum[1] << " ";
		ofstr << std::setw(field_len) << std::left << data.momentum[2] << " ";

		assert(num_bands == data.curvatures.size());
		for(t_size band = 0; band < num_bands; ++band)
		{
			ofstr << std::setw(field_len) << std::left << data.energies[band] << " ";
			ofstr << std::setw(field_len) << std::left << data.weights[band] << " ";
			ofstr << std::setw(field_len) << std::left << data.curvatures[band].real() << " ";
			ofstr << std::setw(field_len) << std::left << data.curvatures[band].imag() << " ";
		}

		ofstr << "\n";
	}

	ofstr.flush();
}



/**
 * toggle between "calculate" and "stop" button
 */
void TopologyDlg::EnableBerryCurvatureCalculation(bool enable)
{
	m_calcEnabled_bc = enable;

	if(enable)
	{
		m_btnStartStop_bc->setText("Calculate");
		m_btnStartStop_bc->setToolTip("Start calculation.");
		m_btnStartStop_bc->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop_bc->setText("Stop");
		m_btnStartStop_bc->setToolTip("Stop running calculation.");
		m_btnStartStop_bc->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * set the berry curvature's Q positions to the main window dispersion Qs
 */
void TopologyDlg::SetBerryCurvatureQ()
{
	if(m_Qstart.size() < 3 || m_Qend.size() < 3)
		return;

	for(int i = 0; i < 3; ++i)
	{
		m_Q_start_bc[i]->setValue(m_Qstart[i]);
		m_Q_end_bc[i]->setValue(m_Qend[i]);
	}
}
// ============================================================================
