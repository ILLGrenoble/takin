/**
 * magnetic dynamics -- 3d dispersion plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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
#include <boost/algorithm/string/replace.hpp>
#include <boost/asio.hpp>
namespace algo = boost::algorithm;
namespace asio = boost::asio;

#include <cstdlib>
#include <mutex>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <limits>

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QFileDialog>

#include "dispersion3d.h"
#include "helper.h"

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/str.h"



/**
 * column indices in magnon band table
 */
enum : int
{
	COL_BC_BAND = 0,
	COL_BC_ACTIVE,
	NUM_COLS_BC,
};



/**
 * sets up the topology dialog
 */
Dispersion3DDlg::Dispersion3DDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("3D Dispersion");
	setSizeGripEnabled(true);

	// create gl plotter
	m_dispplot = new tl2::GlPlot(this);
	m_dispplot->GetRenderer()->SetRestrictCamTheta(false);
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 50, 50, 50 }));
	m_dispplot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -50, -50, -50 }));
	m_dispplot->GetRenderer()->SetCoordMax(50.);
	m_dispplot->GetRenderer()->GetCamera().SetParalellRange(100.);
	m_dispplot->GetRenderer()->GetCamera().SetFOV(tl2::d2r<t_real>(g_structplot_fov));
	m_dispplot->GetRenderer()->GetCamera().SetDist(40.);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// magnon band table
	QWidget *bands_panel = new QWidget(this);
	m_table_bands = new QTableWidget(bands_panel);
	m_table_bands->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_table_bands->setShowGrid(true);
	m_table_bands->setSortingEnabled(false);
	m_table_bands->setSelectionBehavior(QTableWidget::SelectRows);
	m_table_bands->setSelectionMode(QTableWidget::SingleSelection);
	m_table_bands->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_table_bands->verticalHeader()->setVisible(false);
	m_table_bands->setColumnCount(NUM_COLS_BC);
	m_table_bands->setHorizontalHeaderItem(COL_BC_BAND, new QTableWidgetItem{"Band"});
	m_table_bands->setHorizontalHeaderItem(COL_BC_ACTIVE, new QTableWidgetItem{"Act."});
	m_table_bands->setColumnWidth(COL_BC_BAND, 40);
	m_table_bands->setColumnWidth(COL_BC_ACTIVE, 25);
	m_table_bands->resizeColumnsToContents();

	m_only_pos_E = new QCheckBox("E ≥ 0", bands_panel);
	//m_only_pos_E->setChecked(true);
	m_only_pos_E->setToolTip("Ignore magnon annihilation.");

	// splitter for plot and magnon band list
	m_split_plot = new QSplitter(this);
	m_split_plot->setOrientation(Qt::Horizontal);
	m_split_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	m_split_plot->addWidget(m_dispplot);
	m_split_plot->addWidget(bands_panel);
	m_split_plot->setCollapsible(0, false);
	m_split_plot->setCollapsible(1, true);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(m_dispplot), 24);
	m_split_plot->setStretchFactor(m_split_plot->indexOf(bands_panel), 1);

	// general plot context menu
	m_context = new QMenu(this);
	QAction *acCentre = new QAction("Centre Camera", m_context);
	QAction *acSaveData = new QAction("Save Data...", m_context);
	QAction *acSaveScript = new QAction("Save Data As Script...", m_context);
	acSaveData->setIcon(QIcon::fromTheme("text-x-generic"));
	acSaveScript->setIcon(QIcon::fromTheme("text-x-script"));
	m_context->addAction(acCentre);
	m_context->addSeparator();
	m_context->addAction(acSaveData);
	m_context->addAction(acSaveScript);

	// context menu for sites
	m_context_band = new QMenu(this);
	QAction *acCentreOnObject = new QAction("Centre Camera on Band", m_context_band);
	m_context_band->addAction(acCentre);
	m_context_band->addAction(acCentreOnObject);
	m_context_band->addSeparator();
	m_context_band->addAction(acSaveData);
	m_context_band->addAction(acSaveScript);

	// Q coordinates
	QGroupBox *groupQ = new QGroupBox("Dispersion", this);
	m_Q_origin[0] = new QDoubleSpinBox(groupQ);
	m_Q_origin[1] = new QDoubleSpinBox(groupQ);
	m_Q_origin[2] = new QDoubleSpinBox(groupQ);
	m_Q_dir1[0] = new QDoubleSpinBox(groupQ);
	m_Q_dir1[1] = new QDoubleSpinBox(groupQ);
	m_Q_dir1[2] = new QDoubleSpinBox(groupQ);
	m_Q_dir2[0] = new QDoubleSpinBox(groupQ);
	m_Q_dir2[1] = new QDoubleSpinBox(groupQ);
	m_Q_dir2[2] = new QDoubleSpinBox(groupQ);
	m_num_Q_points[0] = new QSpinBox(groupQ);
	m_num_Q_points[1] = new QSpinBox(groupQ);
	m_num_Q_points[0]->setToolTip("Number of grid points along the first momentum axis.");
	m_num_Q_points[1]->setToolTip("Number of grid points along the second momentum axis.");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_origin[i]->setDecimals(4);
		m_Q_origin[i]->setMinimum(-99.9999);
		m_Q_origin[i]->setMaximum(+99.9999);
		m_Q_origin[i]->setSingleStep(0.01);
		m_Q_origin[i]->setValue(0.);
		//m_Q_origin[i]->setSuffix(" rlu");
		m_Q_origin[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_origin[i]->setPrefix(hklPrefix[i]);
		m_Q_origin[i]->setToolTip("Starting momentum transfer.");

		m_Q_dir1[i]->setDecimals(4);
		m_Q_dir1[i]->setMinimum(-99.9999);
		m_Q_dir1[i]->setMaximum(+99.9999);
		m_Q_dir1[i]->setSingleStep(0.01);
		m_Q_dir1[i]->setValue(i == 0 ? 1. : 0.);
		//m_Q_dir1[i]->setSuffix(" rlu");
		m_Q_dir1[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir1[i]->setPrefix(hklPrefix[i]);
		m_Q_dir1[i]->setToolTip("Direction of momentum transfer along the first axis.");

		m_Q_dir2[i]->setDecimals(4);
		m_Q_dir2[i]->setMinimum(-99.9999);
		m_Q_dir2[i]->setMaximum(+99.9999);
		m_Q_dir2[i]->setSingleStep(0.01);
		m_Q_dir2[i]->setValue(i == 1 ? 1. : 0.);
		//m_Q_dir2[i]->setSuffix(" rlu");
		m_Q_dir2[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir2[i]->setPrefix(hklPrefix[i]);
		m_Q_dir2[i]->setToolTip("Direction of momentum transfer along the second axis.");
	}

	for(int i = 0; i < 2; ++i)
	{
		m_num_Q_points[i]->setMinimum(1);
		m_num_Q_points[i]->setMaximum(9999);
		m_num_Q_points[i]->setSingleStep(1);
		m_num_Q_points[i]->setValue(64);
	}

	// main dispersion button
	QPushButton *btnMainQ = new QPushButton("From Main Q", groupQ);
	btnMainQ->setToolTip("Set the Q origin and directions from the dispersion in the main window.");

	// minimum cutoff for filtering S(Q, E)
	m_S_filter_enable = new QCheckBox("Minimum S(Q, E):", groupQ);
	m_S_filter_enable->setChecked(false);
	m_S_filter_enable->setToolTip("Enable minimum S(Q, E).");

	m_S_filter = new QDoubleSpinBox(groupQ);
	m_S_filter->setDecimals(5);
	m_S_filter->setMinimum(0.);
	m_S_filter->setMaximum(9999.99999);
	m_S_filter->setSingleStep(0.01);
	m_S_filter->setValue(0.01);
	m_S_filter->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	m_S_filter->setToolTip("Minimum S(Q, E) to keep.");

	// unite degenerate energies
	m_unite_degeneracies = new QCheckBox("Unite degenerate Es.", groupQ);
	m_unite_degeneracies->setChecked(false);
	m_unite_degeneracies->setToolTip("Unite degenerate energies.");

	// Q and E scale for plot
	QGroupBox *groupPlotOptions = new QGroupBox("Plot Options", this);
	m_Q_scale1 = new QDoubleSpinBox(groupPlotOptions);
	m_Q_scale2 = new QDoubleSpinBox(groupPlotOptions);
	m_E_scale = new QDoubleSpinBox(groupPlotOptions);

	m_Q_scale1->setToolTip("Scaling factor along the first momentum axis.");
	m_Q_scale2->setToolTip("Scaling factor along the second momentum axis.");
	m_E_scale->setToolTip("Scaling factor along the energy axis.");

	for(QDoubleSpinBox *box : { m_Q_scale1, m_Q_scale2, m_E_scale })
	{
		box->setDecimals(3);
		box->setMinimum(0.001);
		box->setMaximum(999.99);
		box->setSingleStep(box == m_E_scale ? 0.1 : 0.5);
		box->setValue(box == m_E_scale ? 4. : 32.);
		box->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
	}

	// camera view angle
	m_cam_phi = new QDoubleSpinBox(this);
	m_cam_phi->setRange(0., 360.);
	m_cam_phi->setSingleStep(1.);
	m_cam_phi->setDecimals(std::max(g_prec_gui - 2, 2));
	m_cam_phi->setPrefix("φ = ");
	m_cam_phi->setSuffix("°");
	m_cam_phi->setToolTip("Camera polar rotation angle φ.");

	m_cam_theta = new QDoubleSpinBox(this);
	m_cam_theta->setRange(-180., 180.);
	m_cam_theta->setSingleStep(1.);
	m_cam_theta->setDecimals(std::max(g_prec_gui - 2, 2));
	m_cam_theta->setPrefix("θ = ");
	m_cam_theta->setSuffix("°");
	m_cam_theta->setToolTip("Camera azimuthal rotation angle θ.");

	// camera perspective projection
	m_perspective = new QCheckBox("Perspective Projection", this);
	m_perspective->setToolTip("Switch between perspective and parallel projection.");
	m_perspective->setChecked(true);

	// progress bar
	m_progress = new QProgressBar(this);
	m_progress->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// start/stop button
	m_btn_start_stop = new QPushButton("Calculate", this);

	// status bar
	m_status = new QLabel(this);
	m_status->setFrameShape(QFrame::Panel);
	m_status->setFrameShadow(QFrame::Sunken);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);

	// close button
	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	QPushButton *btnSaveData = btnbox->addButton("Save Data...", QDialogButtonBox::ActionRole);
	QPushButton *btnSaveScript = btnbox->addButton("Save Script...", QDialogButtonBox::ActionRole);
	btnSaveData->setIcon(QIcon::fromTheme("text-x-generic"));
	btnSaveScript->setIcon(QIcon::fromTheme("text-x-script"));
	btnbox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// bands panel grid
	int y = 0;
	QGridLayout *grid_bands = new QGridLayout(bands_panel);
	grid_bands->setSpacing(4);
	grid_bands->setContentsMargins(6, 6, 6, 6);
	grid_bands->addWidget(m_table_bands, y++, 0, 1, 1);
	grid_bands->addWidget(m_only_pos_E, y++, 0, 1, 1);

	// Q coordinates grid
	y = 0;
	QGridLayout *Qgrid = new QGridLayout(groupQ);
	Qgrid->setSpacing(4);
	Qgrid->setContentsMargins(6, 6, 6, 6);
	Qgrid->addWidget(new QLabel("Q Origin:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_origin[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_origin[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_origin[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Q Direction 1:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_dir1[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_dir1[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_dir1[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Q Direction 2:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_dir2[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_dir2[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_dir2[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Q Grid Points:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_num_Q_points[0], y, 1, 1, 1);
	Qgrid->addWidget(m_num_Q_points[1], y, 2, 1, 1);
	Qgrid->addWidget(btnMainQ, y++, 3, 1, 1);
	Qgrid->addWidget(m_S_filter_enable, y, 0, 1, 1);
	Qgrid->addWidget(m_S_filter, y, 1, 1, 1);
	Qgrid->addWidget(m_unite_degeneracies, y++, 2, 1, 1);

	// plot options grid
	y = 0;
	QGridLayout *plot_options_grid = new QGridLayout(groupPlotOptions);
	plot_options_grid->setSpacing(4);
	plot_options_grid->setContentsMargins(6, 6, 6, 6);
	plot_options_grid->addWidget(new QLabel("Q Scale:", this), y, 0, 1, 1);
	plot_options_grid->addWidget(m_Q_scale1, y, 1, 1, 1);
	plot_options_grid->addWidget(m_Q_scale2, y, 2, 1, 1);
	plot_options_grid->addWidget(new QLabel("E Scale:", this), y, 3, 1, 1);
	plot_options_grid->addWidget(m_E_scale, y++, 4, 1, 1);
	plot_options_grid->addWidget(new QLabel("Camera Angles:", this), y, 0, 1, 1);
	plot_options_grid->addWidget(m_cam_phi, y, 1, 1, 1);
	plot_options_grid->addWidget(m_cam_theta, y, 2, 1, 1);
	plot_options_grid->addWidget(m_perspective, y++, 3, 1, 2);

	// status grid
	QWidget *status_panel = new QWidget(this);
	QGridLayout *status_grid = new QGridLayout(status_panel);
	status_grid->setSpacing(0);
	status_grid->setContentsMargins(0, 0, 0, 0);
	status_grid->addWidget(m_status, y, 0, 1, 3);
	status_grid->addWidget(btnbox, y++, 3, 1, 1);

	// main grid
	y = 0;
	QGridLayout *maingrid = new QGridLayout(this);
	maingrid->setSpacing(4);
	maingrid->setContentsMargins(8, 8, 8, 8);
	maingrid->addWidget(m_split_plot, y++, 0, 1, 4);
	maingrid->addWidget(groupQ, y++, 0, 1, 4);
	maingrid->addWidget(groupPlotOptions, y++, 0, 1, 4);
	maingrid->addWidget(m_progress, y, 0, 1, 3);
	maingrid->addWidget(m_btn_start_stop, y++, 3, 1, 1);
	maingrid->addWidget(status_panel, y++, 0, 1, 4);

	// restore settings
	if(m_sett)
	{
		if(m_sett->contains("dispersion3d/geo"))
			restoreGeometry(m_sett->value("dispersion3d/geo").toByteArray());
		else
			resize(640, 640);

		if(m_sett->contains("dispersion3d/splitter"))
			m_split_plot->restoreState(m_sett->value("dispersion3d/splitter").toByteArray());
	}

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &Dispersion3DDlg::accept);
	connect(acCentreOnObject, &QAction::triggered, this, &Dispersion3DDlg::CentrePlotCameraOnObject);
	connect(acCentre, &QAction::triggered, this, &Dispersion3DDlg::CentrePlotCamera);
	connect(acSaveData, &QAction::triggered, this, &Dispersion3DDlg::SaveData);
	connect(acSaveScript, &QAction::triggered, this, &Dispersion3DDlg::SaveScript);
	connect(btnSaveData, &QAbstractButton::clicked, this, &Dispersion3DDlg::SaveData);
	connect(btnSaveScript, &QAbstractButton::clicked, this, &Dispersion3DDlg::SaveScript);

	connect(m_dispplot, &tl2::GlPlot::AfterGLInitialisation,
		this, &Dispersion3DDlg::AfterPlotGLInitialisation);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
		this, &Dispersion3DDlg::PlotPickerIntersection);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::CameraHasUpdated,
		this, &Dispersion3DDlg::PlotCameraHasUpdated);
	connect(m_dispplot, &tl2::GlPlot::MouseClick, this, &Dispersion3DDlg::PlotMouseClick);
	connect(m_dispplot, &tl2::GlPlot::MouseDown, this, &Dispersion3DDlg::PlotMouseDown);
	connect(m_dispplot, &tl2::GlPlot::MouseUp, this, &Dispersion3DDlg::PlotMouseUp);
	connect(m_perspective, &QCheckBox::toggled, this, &Dispersion3DDlg::SetPlotPerspectiveProjection);
	connect(m_only_pos_E, &QCheckBox::toggled, [this]() { Plot(true); });
	connect(btnMainQ, &QAbstractButton::clicked, this, &Dispersion3DDlg::FromMainQ);
	connect(m_S_filter_enable, &QCheckBox::toggled, m_S_filter, &QDoubleSpinBox::setEnabled);

	for(QDoubleSpinBox *box : { m_Q_scale1, m_Q_scale2, m_E_scale })
	{
		connect(box, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]() { Plot(false); });
	}

	connect(m_cam_phi,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl phi) -> void
	{
		this->SetPlotCameraRotation(phi, m_cam_theta->value());
	});

	connect(m_cam_theta,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl theta) -> void
	{
		this->SetPlotCameraRotation(m_cam_phi->value(), theta);
	});

	// calculation
	connect(m_btn_start_stop, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calc_enabled)
			Calculate();
		else
			m_stop_requested = true;
	});

	m_S_filter->setEnabled(m_S_filter_enable->isChecked());
	EnableCalculation(true);
}



Dispersion3DDlg::~Dispersion3DDlg()
{
}



/**
 * set a pointer to the main magdyn kernel
 */
void Dispersion3DDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;
}



/**
 * save the Q start and end points from the main window's dispersion
 */
void Dispersion3DDlg::SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend)
{
	m_Qstart = Qstart;
	m_Qend = Qend;
}



/**
 * set the Q position and directions from the main window's Q start and end points
 */
void Dispersion3DDlg::FromMainQ()
{
	if(m_Qstart.size() < 3 || m_Qend.size() < 3 || !m_dyn)
		return;

	const t_mat_real& xtalB = m_dyn->GetCrystalBTrafo();
	const t_vec_real* plane = m_dyn->GetScatteringPlane();
	if(!plane)
		return;

	// direction 1 is from the start to the end point
	t_vec_real Qdir1 = m_Qend - m_Qstart;
	// direction 2 is perpendicular to direction 1 inside the scattering plane
	t_vec_real Qdir2 = tl2::cross(xtalB, plane[2], Qdir1);

	for(int i = 0; i < 3; ++i)
	{
		m_Q_origin[i]->setValue(m_Qstart[i]);
		m_Q_dir1[i]->setValue(Qdir1[i]);
		m_Q_dir2[i]->setValue(Qdir2[i]);
	}
}



void Dispersion3DDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * get the Q origin and direction vectors
 */
std::tuple<t_vec_real, t_vec_real, t_vec_real> Dispersion3DDlg::GetQVectors() const
{
	t_vec_real Q_origin = tl2::create<t_vec_real>(
	{
		m_Q_origin[0]->value(),
		m_Q_origin[1]->value(),
		m_Q_origin[2]->value(),
	});

	t_vec_real Q_dir_1 = tl2::create<t_vec_real>(
	{
		m_Q_dir1[0]->value(),
		m_Q_dir1[1]->value(),
		m_Q_dir1[2]->value(),
	});

	t_vec_real Q_dir_2 = tl2::create<t_vec_real>(
	{
		m_Q_dir2[0]->value(),
		m_Q_dir2[1]->value(),
		m_Q_dir2[2]->value(),
	});

	return std::make_tuple(std::move(Q_origin), std::move(Q_dir_1), std::move(Q_dir_2));
}



/**
 * calculate the dispersion
 */
void Dispersion3DDlg::Calculate()
{
	if(!m_dyn)
		return;

	m_max_E = -std::numeric_limits<t_real>::max();
	m_data.clear();

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	// get coordinates
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	m_Q_count_1 = m_num_Q_points[0]->value();
	m_Q_count_2 = m_num_Q_points[1]->value();

	t_vec_real Q_step_1 = Q_dir_1 / t_real(m_Q_count_1);
	t_vec_real Q_step_2 = Q_dir_2 / t_real(m_Q_count_2);

	t_real min_S = m_S_filter->value();
	bool use_weights = m_S_filter_enable->isChecked();
	bool use_projector = true;
	bool unite_degen = m_unite_degeneracies->isChecked();

	// calculate the dispersion
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(unite_degen);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stop_requested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(m_Q_count_1 * m_Q_count_2);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(m_Q_count_1 * m_Q_count_2);

	t_size expected_bands = dyn.GetMagneticSitesCount() * 2;
	if(dyn.IsIncommensurate())
		expected_bands *= 3;
	m_data.resize(expected_bands);

	for(t_size Q_idx_1 = 0; Q_idx_1 < m_Q_count_1; ++Q_idx_1)
	for(t_size Q_idx_2 = 0; Q_idx_2 < m_Q_count_2; ++Q_idx_2)
	{
		auto task = [this, &mtx, &dyn, Q_idx_1, Q_idx_2,
			&Q_origin, &Q_step_1, &Q_step_2, expected_bands,
			unite_degen, use_weights, use_projector, min_S]()
		{
			// calculate the dispersion at the given Q point
			t_vec_real Q = Q_origin + Q_step_1*t_real(Q_idx_1) + Q_step_2*t_real(Q_idx_2);
			auto Es_and_S = dyn.CalcEnergies(Q, !use_weights).E_and_S;

			// iterate the energies for this Q point
			t_size data_band_idx = 0;
			for(t_size band_idx = 0; band_idx < Es_and_S.size() && data_band_idx < expected_bands; ++band_idx, ++data_band_idx)
			{
				const auto& E_and_S = Es_and_S[band_idx];

				bool valid = true;
				t_real E = E_and_S.E;
				if(std::isnan(E) || std::isinf(E))
					valid = false;

				t_real weight = -1;
				if(use_weights)
				{
					weight = E_and_S.weight;

					if(!use_projector)
					{
						const t_mat& S = E_and_S.S;
						weight = tl2::trace<t_mat>(S).real();
					}

					// filter invalid S(Q, E)
					if(std::isnan(weight) || std::isinf(weight))
						weight = 0.;

					// filter minimum S(Q, E)
					if(min_S >= 0. && std::abs(weight) <= min_S)
						valid = false;
				}

				if(valid)
					m_max_E = std::max(m_max_E, E);

				// count energy degeneracy
				t_size degeneracy = E_and_S.degeneracy;
				for(t_size band_idx2 = 0; band_idx2 < Es_and_S.size(); ++band_idx2)
				{
					if(band_idx2 == band_idx)
						continue;

					if(tl2::equals(E, Es_and_S[band_idx2].E, g_eps))
						degeneracy += Es_and_S[band_idx2].degeneracy;
				}

				/*if(degeneracy > 1)
				{
					std::cout << "degenerate point: Q indices: " << Q_idx_1 << " " << Q_idx_2
						<< ", band index: " << band_idx << " (" << data_band_idx
						<< "): " << E << " meV (" << degeneracy << "x)" << std::endl;
				}*/

				// generate and add data point
				t_data_Q dat{std::make_tuple(Q, E, weight, Q_idx_1, Q_idx_2, degeneracy, valid)};

				std::lock_guard<std::mutex> _lck{mtx};
				m_data[data_band_idx].emplace_back(std::move(dat));

				if(unite_degen && degeneracy > 1)
				{
					// skip degeneracies to keep matching bands together
					for(t_size band_idx2 = data_band_idx + 1; band_idx2 < data_band_idx + degeneracy; ++band_idx2)
						m_data[band_idx2].emplace_back(std::make_tuple(Q, 0., 0., Q_idx_1, Q_idx_2, 1, false));

					data_band_idx += degeneracy - 1;
				}
			}

			// fill up band data in case some indices were skipped due to invalid hamiltonians
			for(; data_band_idx < expected_bands; ++data_band_idx)
			{
				t_data_Q dat{std::make_tuple(Q, 0., 0., Q_idx_1, Q_idx_2, 1, false)};

				std::lock_guard<std::mutex> _lck{mtx};
				m_data[data_band_idx].emplace_back(std::move(dat));
			}
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText(QString("Calculating dispersion in %1 threads...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		// process events to see if the stop button was clicked
		// only do this for a fraction of the points to avoid gui overhead
		if(task_idx % std::max<t_size>(tasks.size() / std::sqrt(g_stop_check_fraction), 1) == 0)
			qApp->processEvents();

		if(m_stop_requested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress->setValue(task_idx + 1);
	}

	// finish parallel calculations
	pool.join();

	// get sorting of data by Q
	for(t_size band_idx = 0; band_idx < m_data.size(); ++band_idx)
	{
		std::vector<std::size_t> perm = tl2::get_perm(m_data[band_idx].size(),
			[this, band_idx](std::size_t idx1, std::size_t idx2) -> bool
		{
			/*
			// sorting by Q components
			t_real h1 = std::get<0>(m_data[0][idx1])[0];
			t_real k1 = std::get<0>(m_data[0][idx1])[1];
			t_real l1 = std::get<0>(m_data[0][idx1])[2];

			t_real h2 = std::get<0>(m_data[0][idx2])[0];
			t_real k2 = std::get<0>(m_data[0][idx2])[1];
			t_real l2 = std::get<0>(m_data[0][idx2])[2];

			if(!tl2::equals(h1, h2, g_eps))
				return h1 < h2;
			if(!tl2::equals(k1, k2, g_eps))
				return k1 < k2;

			return l1 < l2;*/

			// sorting by Q indices
			t_size Q1_idx_1 = std::get<3>(m_data[band_idx][idx1]);
			t_size Q1_idx_2 = std::get<4>(m_data[band_idx][idx1]);
			t_size Q2_idx_1 = std::get<3>(m_data[band_idx][idx2]);
			t_size Q2_idx_2 = std::get<4>(m_data[band_idx][idx2]);

			if(Q1_idx_1 != Q2_idx_1)
				return Q1_idx_1 < Q2_idx_1;
			return Q1_idx_2 < Q2_idx_2;
		});

		m_data[band_idx] = tl2::reorder(m_data[band_idx], perm);
	}

	// move degenerate points to the bands where most of the other points are
	if(unite_degen && m_data.size())
	{
		for(t_size band_idx = 0; band_idx < m_data.size() - 1; ++band_idx)
		{
			t_data_Qs& band1 = m_data[band_idx];
			t_data_Qs& band2 = m_data[band_idx + 1];
			if(band1.size() != band2.size())
				continue;

			auto [_valid1, invalid1] = NumValid(band1);
			auto [_valid2, invalid2] = NumValid(band2);

			for(t_size Q_idx = 0; Q_idx < band1.size(); ++Q_idx)
			{
				t_size degen1 = std::get<5>(band1[Q_idx]);
				bool valid1 = std::get<6>(band1[Q_idx]);
				bool valid2 = std::get<6>(band2[Q_idx]);

				if(degen1 > 1 && valid1 && !valid2 && invalid1 > invalid2)
					std::swap(band1[Q_idx], band2[Q_idx]);
			}
		}
	}

	// remove fully invalid bands
	for(auto iter = m_data.begin(); iter != m_data.end();)
	{
		if(!IsValid(*iter))
			iter = m_data.erase(iter);
		else
			++iter;
	}

	// sort band energies in descending order
	std::stable_sort(m_data.begin(), m_data.end(), [this](const t_data_Qs& dat1, const t_data_Qs& dat2)
	{
		return GetMeanEnergy(dat1) >= GetMeanEnergy(dat2);
	});

	// show elapsed time
	stopwatch.stop();
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Dispersion calculation";
	if(m_stop_requested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	Plot(true);
}



/**
 * count number of bands with positive energies (should be half of all bands)
 */
t_size Dispersion3DDlg::NumPositive() const
{
	t_size num_pos = 0;

	for(const t_data_Qs& band : m_data)
	{
		if(IsPositive(band))
			++num_pos;
	}

	return num_pos;
}



/**
 * determine if a band only contains positive energies
 */
bool Dispersion3DDlg::IsPositive(const t_data_Qs& data) const
{
	for(const t_data_Q& data : data)
	{
		t_real E = std::get<1>(data);
		if(E < 0.)
			return false;
	}

	return true;
}



/**
 * determine if a band is valid or only contains invalid points
 */
bool Dispersion3DDlg::IsValid(const t_data_Qs& data) const
{
	for(const t_data_Q& data : data)
	{
		// data point valid?
		if(std::get<6>(data))
			return true;
	}

	// all points invalid
	return false;
}



/**
 * count the number of valid and invalid points in a band
 */
std::pair<t_size, t_size> Dispersion3DDlg::NumValid(const t_data_Qs& data) const
{
	t_size valid = 0, invalid = 0;

	for(const t_data_Q& data : data)
	{
		// data point valid?
		if(std::get<6>(data))
			++valid;
		else
			++invalid;
	}

	return std::make_pair(valid, invalid);
}



/**
 * calculate the mean band energy
 */
t_real Dispersion3DDlg::GetMeanEnergy(const Dispersion3DDlg::t_data_Qs& data) const
{
	t_real E_mean = 0.;
	t_size num_pts = 0;

	for(const t_data_Q& data : data)
	{
		// data point invalid?
		if(!std::get<6>(data))
			continue;

		E_mean += std::get<1>(data);
		++num_pts;
	}

	if(num_pts)
		E_mean /= static_cast<t_real>(num_pts);
	return E_mean;
}



/**
 * calculate the mean band energy
 */
t_real Dispersion3DDlg::GetMeanEnergy(t_size band_idx) const
{
	if(band_idx >= m_data.size())
		return 0.;

	return GetMeanEnergy(m_data[band_idx]);
}



/**
 * get a unique colour for a given magnon band
 */
std::array<int, 3> Dispersion3DDlg::GetBranchColour(t_size branch_idx, t_size num_branches) const
{
	if(num_branches <= 1)
		return std::array<int, 3>{{ 0xff, 0x00, 0x00 }};

	return std::array<int, 3>{{
		int(std::lerp(1., 0., t_real(branch_idx) / t_real(num_branches - 1)) * 255.),
		0x00,
		int(std::lerp(0., 1., t_real(branch_idx) / t_real(num_branches - 1)) * 255.),
	}};
}



/**
 * plot the calculated dispersion
 */
void Dispersion3DDlg::Plot(bool clear_settings)
{
	if(!m_dispplot)
		return;

	// keep some settings from previous plot, e.g. the band visibility flags
	std::vector<bool> active_bands;
	if(!clear_settings)
	{
		active_bands.reserve(m_table_bands->rowCount());
		for(int row = 0; row < m_table_bands->rowCount(); ++row)
			active_bands.push_back(IsBandEnabled(t_size(row)));
	}

	t_real E_scale = m_E_scale->value();
	t_real Q_scale1 = m_Q_scale1->value();
	t_real Q_scale2 = m_Q_scale2->value();

	// reset plot
	ClearBands();
	m_cam_centre = tl2::zero<t_vec_gl>(3);
	m_dispplot->GetRenderer()->RemoveObjects();
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>(
		{ static_cast<t_real_gl>(1.5 * Q_scale1), static_cast<t_real_gl>(1.5 * Q_scale2),
			static_cast<t_real_gl>(2.5 * (m_max_E + 4.) * E_scale) }));
	m_dispplot->GetRenderer()->SetLight(1, -tl2::create<t_vec3_gl>(
		{ static_cast<t_real_gl>(1.5 * Q_scale1), static_cast<t_real_gl>(1.5 * Q_scale2),
			static_cast<t_real_gl>(2.5 * (m_max_E + 4.) * E_scale) }));

	// plot the magnon bands
	t_size num_bands = m_data.size();
	if(m_only_pos_E->isChecked())
	{
		num_bands = NumPositive();
		//num_bands /= 2;
	}
	t_size num_active_bands = 0;
	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		bool band_active = band_idx < active_bands.size() ? active_bands[band_idx] : true;

		// colour for this magnon band
		std::array<int, 3> col = GetBranchColour(band_idx, num_bands);
		const QColor colFull(col[0], col[1], col[2]);

		if(band_active)
		{
			t_data_Qs& data = m_data[band_idx];

			auto patch_fkt = [this, &data, E_scale](
				t_real_gl /*x2*/, t_real_gl /*x1*/, t_size idx_1, t_size idx_2)
					-> std::pair<t_real_gl, bool>
			{
				t_size idx = idx_1 * m_Q_count_2 + idx_2;
				if(idx >= data.size())
					return std::make_pair(0., false);

				bool valid = std::get<6>(data[idx]);
				if(std::get<3>(data[idx]) != idx_1 || std::get<4>(data[idx]) != idx_2)
				{
					std::cerr << "Error: Patch index mismatch: "
						<< "Expected " << std::get<3>(data[idx]) << " for x, but got " << idx_1
						<< "; expected " << std::get<4>(data[idx]) << " for y, but got " << idx_2
						<< "." << std::endl;
					valid = false;
				}

				t_real_gl E = std::get<1>(data[idx]);
				return std::make_pair(E * E_scale, valid);
			};

			std::ostringstream objLabel;
			objLabel << "Band #" << (band_idx + 1);

			t_real_gl r = t_real_gl(col[0]) / t_real_gl(255.);
			t_real_gl g = t_real_gl(col[1]) / t_real_gl(255.);
			t_real_gl b = t_real_gl(col[2]) / t_real_gl(255.);

			std::size_t obj = m_dispplot->GetRenderer()->AddPatch(patch_fkt, 0., 0., 0.,
				Q_scale1, Q_scale2, m_Q_count_1, m_Q_count_2, r, g, b);
			m_dispplot->GetRenderer()->SetObjectLabel(obj, objLabel.str());
			m_band_objs.insert(std::make_pair(obj, band_idx));

			m_cam_centre[2] += GetMeanEnergy(band_idx);
			++num_active_bands;
		}

		AddBand("#" + tl2::var_to_str(band_idx + 1), colFull, band_active);
	}

	if(num_active_bands)
		m_cam_centre[2] /= static_cast<t_real>(num_active_bands);

	if(clear_settings)
		CentrePlotCamera();    // centre camera and update
	else
		m_dispplot->update();  // only update renderer
}



/**
 * clears the table of magnon bands
 */
void Dispersion3DDlg::ClearBands()
{
	m_band_objs.clear();
	m_cur_obj = std::nullopt;

	m_table_bands->clearContents();
	m_table_bands->setRowCount(0);
}



/**
 * adds a magnon band to the table
 */
void Dispersion3DDlg::AddBand(const std::string& name, const QColor& colour, bool enabled)
{
	if(!m_table_bands)
		return;

	int row = m_table_bands->rowCount();
	m_table_bands->insertRow(row);

	QTableWidgetItem *item = new QTableWidgetItem{name.c_str()};
	item->setFlags(item->flags() & ~Qt::ItemIsEditable);

	QBrush bg = item->background();
	bg.setColor(colour);
	bg.setStyle(Qt::SolidPattern);
	item->setBackground(bg);

	QBrush fg = item->foreground();
	fg.setColor(QColor{0xff, 0xff, 0xff});
	fg.setStyle(Qt::SolidPattern);
	item->setForeground(fg);

	QCheckBox *checkBand = new QCheckBox(m_table_bands);
	checkBand->setChecked(enabled);
	connect(checkBand, &QCheckBox::toggled, [this]() { Plot(false); });

	m_table_bands->setItem(row, COL_BC_BAND, item);
	m_table_bands->setCellWidget(row, COL_BC_ACTIVE, checkBand);
}



/**
 * verifies if the band's checkbox is checked
 */
bool Dispersion3DDlg::IsBandEnabled(t_size idx) const
{
	if(!m_table_bands || int(idx) >= m_table_bands->rowCount())
		return true;

	QCheckBox* box = reinterpret_cast<QCheckBox*>(m_table_bands->cellWidget(int(idx), COL_BC_ACTIVE));
	if(!box)
		return true;

	return box->isChecked();
}



/**
 * toggle between "calculate" and "stop" button
 */
void Dispersion3DDlg::EnableCalculation(bool enable)
{
	m_calc_enabled = enable;

	if(enable)
	{
		m_btn_start_stop->setText("Calculate");
		m_btn_start_stop->setToolTip("Start dispersion calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btn_start_stop->setText("Stop");
		m_btn_start_stop->setToolTip("Stop running dispersion calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * mouse intersection with dispersion band
 */
void Dispersion3DDlg::PlotPickerIntersection(
	[[maybe_unused]] const t_vec3_gl* pos,
	[[maybe_unused]] std::size_t objIdx,
	[[maybe_unused]] std::size_t triagIdx,
	[[maybe_unused]] const t_vec3_gl* posSphere)
{
	m_status->setText("");
	m_cur_obj = std::nullopt;

	m_dispplot->GetRenderer()->SetObjectsHighlight(false);

	if(!pos)
		return;

	m_cur_obj = objIdx;

	// get Q and E position at cursor intersection
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	// reconstruct Q position
	t_real_gl Q1param = (0.5*m_Q_scale1->value() + (*pos)[0]) / m_Q_scale1->value();
	t_real_gl Q2param = (0.5*m_Q_scale2->value() + (*pos)[1]) / m_Q_scale2->value();
	t_vec_real Q = Q_origin + Q_dir_1*Q1param + Q_dir_2*Q2param;

	t_real_gl E = (*pos)[2] / m_E_scale->value();

	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr
		<< "Q = (" << Q[0] << ", " << Q[1] << ", " << Q[2] << ") rlu, "
		<< "E = " << E << " meV";

	const std::string& label = m_dispplot->GetRenderer()->GetObjectLabel(objIdx);
	if(label != "")
		ostr << ", " << label;
	ostr << ".";

	m_status->setText(ostr.str().c_str());
}



/**
 * the dispersion plot's camera properties have been updated
 */
void Dispersion3DDlg::PlotCameraHasUpdated()
{
	auto [phi, theta] = m_dispplot->GetRenderer()->GetCamera().GetRotation();

	phi = tl2::r2d<t_real>(phi);
	theta = tl2::r2d<t_real>(theta);

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_cam_phi->blockSignals(false);
		this_->m_cam_theta->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_cam_phi->blockSignals(true);
	m_cam_theta->blockSignals(true);

	m_cam_phi->setValue(phi);
	m_cam_theta->setValue(theta);
}



/**
 * dispersion plot mouse button clicked
 */
void Dispersion3DDlg::PlotMouseClick(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(right)
	{
		const QPointF& _pt = m_dispplot->GetRenderer()->GetMousePosition();
		QPoint pt = m_dispplot->mapToGlobal(_pt.toPoint());

		if(m_cur_obj && m_band_objs.find(*m_cur_obj) != m_band_objs.end())
		{
			// band selected
			m_context_band->popup(pt);
		}
		else
		{
			// no band selected
			m_context->popup(pt);
		}
	}
}



/**
 * dispersion plot mouse button pressed
 */
void Dispersion3DDlg::PlotMouseDown(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * dispersion plot mouse button released
 */
void Dispersion3DDlg::PlotMouseUp(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * dispersion plot has initialised
 */
void Dispersion3DDlg::AfterPlotGLInitialisation()
{
	if(!m_dispplot)
		return;

	// GL device infos
	auto [ver, shader_ver, vendor, renderer]
		= m_dispplot->GetRenderer()->GetGlDescr();

	emit GlDeviceInfos(ver, shader_ver, vendor, renderer);

	m_dispplot->GetRenderer()->SetCull(false);

	ShowPlotCoordCross(false);
	ShowPlotLabels(false);
	SetPlotPerspectiveProjection(m_perspective->isChecked());
	//SetPlotCoordinateSystem(m_coordsys->currentIndex());

	PlotCameraHasUpdated();
}



/**
 * show or hide the coordinate system
 */
void Dispersion3DDlg::ShowPlotCoordCross(bool show)
{
	if(auto obj = m_dispplot->GetRenderer()->GetCoordCross(); obj)
	{
		m_dispplot->GetRenderer()->SetObjectVisible(*obj, show);
		m_dispplot->update();
	}
}



/**
 * show or hide the object labels
 */
void Dispersion3DDlg::ShowPlotLabels(bool show)
{
	m_dispplot->GetRenderer()->SetLabelsVisible(show);
	m_dispplot->update();
}



/**
 * choose between perspective or parallel projection
 */
void Dispersion3DDlg::SetPlotPerspectiveProjection(bool proj)
{
	m_dispplot->GetRenderer()->GetCamera().SetPerspectiveProjection(proj);
	m_dispplot->GetRenderer()->RequestViewportUpdate();
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * sets the camera's rotation angles
 */
void Dispersion3DDlg::SetPlotCameraRotation(t_real_gl phi, t_real_gl theta)
{
	phi = tl2::d2r<t_real>(phi);
	theta = tl2::d2r<t_real>(theta);

	m_dispplot->GetRenderer()->GetCamera().SetRotation(phi, theta);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	PlotCameraHasUpdated();
	m_dispplot->update();
}



/**
 * centre camera on currently selected object
 */
void Dispersion3DDlg::CentrePlotCameraOnObject()
{
	if(!m_cur_obj)
		return;

	t_mat_gl mat = m_dispplot->GetRenderer()->GetObjectMatrix(*m_cur_obj);

	// selected a band?
	auto band_iter = m_band_objs.find(*m_cur_obj);
	if(band_iter != m_band_objs.end())
	{
		// translate camera to the band's mean energy position
		t_real E_mean = GetMeanEnergy(band_iter->second) * m_E_scale->value();
		mat(2, 3) = E_mean;
	}

	m_dispplot->GetRenderer()->GetCamera().Centre(mat);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * centre camera on central position
 */
void Dispersion3DDlg::CentrePlotCamera()
{
	t_mat_gl matCentre = tl2::hom_translation<t_mat_gl>(
		m_cam_centre[0] * m_Q_scale2->value(),
		m_cam_centre[1] * m_Q_scale1->value(),
		m_cam_centre[2] * m_E_scale->value());

	m_dispplot->GetRenderer()->GetCamera().Centre(matCentre);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->update();
}



/**
 * switch between crystal and lab coordinates
 */
void Dispersion3DDlg::SetPlotCoordinateSystem(int which)
{
	m_dispplot->GetRenderer()->SetCoordSys(which);
}



/**
 * dialog is closing
 */
void Dispersion3DDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("dispersion3d/geo", saveGeometry());
		m_sett->setValue("dispersion3d/splitter", m_split_plot->saveState());
	}

	QDialog::accept();
}



/**
 * write the meta data header for text file exports
 */
void Dispersion3DDlg::WriteHeader(std::ostream& ostr) const
{
	using namespace tl2_ops;

	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	const t_size num_bands = m_data.size();
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	ostr << "#\n"
		<< "# Created by Takin/Magdyn\n"
		<< "# URL: https://github.com/ILLGrenoble/takin\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.4117437\n"
		<< "# User: " << user << "\n"
		<< "# Date: " << tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()) << "\n"
		<< "#\n# branch_count: " << num_bands << "\n"
		<< "# Q_origin: " << Q_origin << "\n"
		<< "# Q_direction_1: " << Q_dir_1 << "\n"
		<< "# Q_direction_2: " << Q_dir_2 << "\n"
		<< "#\n\n";
}



/**
 * save the dispersion as a text data file
 */
void Dispersion3DDlg::SaveData()
{
	bool skip_invalid_points = true;
	bool use_weights = m_S_filter_enable->isChecked();

	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("dispersion3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Dispersion Data",
		dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("dispersion3d/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".")
			.arg(filename).toStdString().c_str());
		return;
	}

	ofstr.precision(g_prec);
	WriteHeader(ofstr);

	// write column header
	int field_len = g_prec * 2.5;
	ofstr << std::setw(field_len) << std::left << "# h" << " ";
	ofstr << std::setw(field_len) << std::left << "k" << " ";
	ofstr << std::setw(field_len) << std::left << "l" << " ";
	ofstr << std::setw(field_len) << std::left << "E" << " ";
	if(!skip_invalid_points)
		ofstr << std::setw(field_len) << std::left << "valid" << " ";
	if(use_weights)
		ofstr << std::setw(field_len) << std::left << "S" << " ";
	ofstr << std::setw(field_len) << std::left << "band" << " ";
	ofstr << std::setw(field_len) << std::left << "Qidx1" << " ";
	ofstr << std::setw(field_len) << std::left << "Qidx2" << " ";
	ofstr << std::setw(field_len) << std::left << "degen" << "\n";

	const t_size num_bands = m_data.size();
	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		for(t_data_Q& data : m_data[band_idx])
		{
			const t_vec_real& Q = std::get<0>(data);
			t_real E = std::get<1>(data);
			t_real S = std::get<2>(data);
			t_size Qidx1 = std::get<3>(data);
			t_size Qidx2 = std::get<4>(data);
			t_size degen = std::get<5>(data);
			bool valid = std::get<6>(data);

			if(skip_invalid_points && !valid)
				continue;

			ofstr << std::setw(field_len) << std::left << Q[0] << " ";
			ofstr << std::setw(field_len) << std::left << Q[1] << " ";
			ofstr << std::setw(field_len) << std::left << Q[2] << " ";
			ofstr << std::setw(field_len) << std::left << E << " ";
			if(!skip_invalid_points)
				ofstr << std::setw(field_len) << std::left << valid << " ";
			if(use_weights)
				ofstr << std::setw(field_len) << std::left << S << " ";
			ofstr << std::setw(field_len) << std::left << band_idx << " ";
			ofstr << std::setw(field_len) << std::left << Qidx1 << " ";
			ofstr << std::setw(field_len) << std::left << Qidx2 << " ";
			ofstr << std::setw(field_len) << std::left << degen << "\n";
		}
	}

	ofstr.flush();
}



/**
 * save the dispersion as a script file
 */
void Dispersion3DDlg::SaveScript()
{
	using namespace tl2_ops;

	bool skip_invalid_points = true;
	bool use_weights = m_S_filter_enable->isChecked();

	if(m_data.size() == 0)
		return;

	QString dirLast;
	if(m_sett)
		dirLast = m_sett->value("dispersion3d/dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Dispersion Data As Script",
		dirLast, "Py Files (*.py)");
	if(filename == "")
		return;
	if(m_sett)
		m_sett->setValue("dispersion3d/dir", QFileInfo(filename).path());

	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		ShowError(QString("Could not save data to file \"%1\".")
			.arg(filename).toStdString().c_str());
		return;
	}

	ofstr.precision(g_prec);
	WriteHeader(ofstr);

	// ------------------------------------------------------------------------
	// plot script
	std::string pyscr = R"RAW(import sys
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
plot_file      = ""     # file to save plot to
only_pos_E     = %%ONLY_POS_E%%  # ignore magnon annihilation?
S_filter_min   = 1e-5   # cutoff minimum spectral weight

h_column       =  0     # h column index in data files
k_column       =  1     # k column index
l_column       =  2     # l column index
E_column       =  3     # E column index
S_column       = %%S_INDEX%%     # S column index, -1: not available
# -----------------------------------------------------------------------------

def plot_disp(data, branch_data, degen_data, branch_colours, Q_idx1 = 0, Q_idx2 = 1):
	(plt, axis) = pyplot.subplots(nrows = 1, ncols = 1,
		width_ratios = None, sharey = True,
		subplot_kw = { "projection" : "3d", "proj_type" : "%%PROJ_TYPE%%",
			"computed_zorder" : False,
			"azim" : %%AZIMUTH%%, "elev" : %%ELEVATION%% })

	# iterate energy branches
	E_branch_eff_idx = 0             # effective index of actually plotted bands
	E_branch_max = max(branch_data)  # maximum branch index
	for E_branch_idx in range(0, E_branch_max + 1):
		# filter data for given branch
		data_Q = [
			[ row[h_column + Q_idx1] for (row, branch_idx) in zip(data, branch_data) \
				if branch_idx == E_branch_idx ],
			[ row[h_column + Q_idx2] for (row, branch_idx) in zip(data, branch_data) \
				if branch_idx == E_branch_idx ]
		]
		data_E = [ row[E_column] for (row, branch_idx) in zip(data, branch_data) \
			if branch_idx == E_branch_idx ]
		if S_column >= 0:
			data_S = [ row[S_column] for (row, branch_idx) in zip(data, branch_data) \
				if branch_idx == E_branch_idx ]

		colour_idx = E_branch_eff_idx
		if only_pos_E:
			# ignore magnon annihilation
			data_Q[0] = [ Q for (Q, E) in zip(data_Q[0], data_E) if E >= 0. ]
			data_Q[1] = [ Q for (Q, E) in zip(data_Q[1], data_E) if E >= 0. ]
			if S_column >= 0:
				data_S = [ S for (S, E) in zip(data_S, data_E) if E >= 0. ]
			data_E = [ E for E in data_E if E >= 0. ]
			colour_idx *= 2  # skip every other colour if E >= 0

		if S_column >= 0 and S_filter_min >= 0.:
			# filter weights below cutoff
			data_Q[0] = [ Q for (Q, S) in zip(data_Q[0], data_S) if S >= S_filter_min ]
			data_Q[1] = [ Q for (Q, S) in zip(data_Q[1], data_S) if S >= S_filter_min ]
			data_E = [ E for (E, S) in zip(data_E, data_S) if S >= S_filter_min ]
			data_S = [ S for S in data_S if S >= S_filter_min ]

		if(len(data_E) < 1):
			continue

		axis.plot_trisurf(data_Q[0], data_Q[1], data_E,
			color = branch_colours[colour_idx], alpha = 1., shade = True,
			antialiased = False, zorder = E_branch_max - E_branch_idx)
		E_branch_eff_idx += 1

	labels = [ "h (rlu)", "k (rlu)", "l (rlu)" ]
	axis.set_xlabel(labels[Q_idx1])
	axis.set_ylabel(labels[Q_idx2])
	axis.set_zlabel("E (meV)")

	plt.tight_layout()
	plt.subplots_adjust(wspace = 0)

	if plot_file != "":
		pyplot.savefig(plot_file)
	pyplot.show()

if __name__ == "__main__":
	h_data = %%H_DATA%%
	k_data = %%K_DATA%%
	l_data = %%L_DATA%%
	E_data = %%E_DATA%%
	S_data = %%S_DATA%%
	branch_data = %%BRANCH_DATA%%
	degen_data = %%DEGEN_DATA%%
	colours = %%BRANCH_COLOURS%%

	if len(S_data) == len(E_data):
		data = numpy.array([ h_data, k_data, l_data, E_data, S_data ]).T
	else:
		data = numpy.array([ h_data, k_data, l_data, E_data ]).T
	plot_disp(data, branch_data, degen_data, colours, %%Q_IDX_1%%, %%Q_IDX_2%%))RAW";
	// ------------------------------------------------------------------------

	const t_size num_bands = m_data.size();
	auto [Q_origin, Q_dir_1, Q_dir_2] = GetQVectors();

	// create data arrays
	std::ostringstream h_data, k_data, l_data, E_data, S_data, bandidx_data, degen_data, colours;
	for(std::ostringstream* ostr : { &h_data, &k_data, &l_data, &E_data, &S_data,
		&bandidx_data, &degen_data, &colours })
	{
		ostr->precision(g_prec);
	}

	for(t_size band_idx = 0; band_idx < num_bands; ++band_idx)
	{
		// band data
		for(t_data_Q& data : m_data[band_idx])
		{
			const t_vec_real& Q = std::get<0>(data);
			t_real E = std::get<1>(data);
			t_real S = std::get<2>(data);
			//t_size Qidx1 = std::get<3>(data);
			//t_size Qidx2 = std::get<4>(data);
			t_size degen = std::get<5>(data);
			bool valid = std::get<6>(data);

			if(skip_invalid_points && !valid)
				continue;

			h_data << Q[0] << ", ";
			k_data << Q[1] << ", ";
			l_data << Q[2] << ", ";
			E_data << E << ", ";
			if(use_weights)
				S_data << S << ", ";

			bandidx_data << band_idx << ", ";
			degen_data << degen << ", ";
		}

		// band colour
		std::array<int, 3> col = GetBranchColour(band_idx, num_bands);
		colours << "\"#" << std::hex
			<< std::setw(2) << std::setfill('0') << col[0]
			<< std::setw(2) << std::setfill('0') << col[1]
			<< std::setw(2) << std::setfill('0') << col[2]
			<< "\", ";
	}

	// find first Q axis index for plot
	int Q_idx_1 = 0;
	t_real cur_comp1 = Q_dir_1[Q_idx_1];
	for(int i = 0; i < 3; ++i)
	{
		// use largest absolute component as index
		if(std::abs(Q_dir_1[i]) > std::abs(cur_comp1))
		{
			Q_idx_1 = i;
			cur_comp1 = Q_dir_1[i];
		}
	}

	// find second Q axis index for plot
	int Q_idx_2 = (Q_idx_1 + 1) % 3;
	t_real cur_comp2 = Q_dir_2[Q_idx_2];
	for(int i = 0; i < 3; ++i)
	{
		if(i == Q_idx_1)
			continue;

		// use largest absolute component as index
		if(std::abs(Q_dir_2[i]) > std::abs(cur_comp2))
		{
			Q_idx_2 = i;
			cur_comp2 = Q_dir_2[i];
		}
	}

	// TODO: for the moment the azimuth only matches for the [100], [010] Q directions -> add a trafo
	algo::replace_all(pyscr, "%%H_DATA%%", "[ " + h_data.str() + "]");
	algo::replace_all(pyscr, "%%K_DATA%%", "[ " + k_data.str() + "]");
	algo::replace_all(pyscr, "%%L_DATA%%", "[ " + l_data.str() + "]");
	algo::replace_all(pyscr, "%%E_DATA%%", "[ " + E_data.str() + "]");
	algo::replace_all(pyscr, "%%S_DATA%%", "[ " + S_data.str() + "]");
	algo::replace_all(pyscr, "%%S_INDEX%%", use_weights ? " 4" : "-1");
	algo::replace_all(pyscr, "%%BRANCH_DATA%%", "[ " + bandidx_data.str() + "]");
	algo::replace_all(pyscr, "%%DEGEN_DATA%%", "[ " + degen_data.str() + "]");
	algo::replace_all(pyscr, "%%BRANCH_COLOURS%%", "[ " + colours.str() + "]");
	algo::replace_all(pyscr, "%%ONLY_POS_E%%", m_only_pos_E->isChecked() ? "True " : "False");
	algo::replace_all(pyscr, "%%PROJ_TYPE%%", m_perspective->isChecked() ? "persp" : "ortho");
	algo::replace_all(pyscr, "%%AZIMUTH%%", tl2::var_to_str(-90. - m_cam_phi->value(), g_prec));
	algo::replace_all(pyscr, "%%ELEVATION%%", tl2::var_to_str(90. + m_cam_theta->value(), g_prec));
	algo::replace_all(pyscr, "%%Q_IDX_1%%", tl2::var_to_str(Q_idx_1, g_prec));
	algo::replace_all(pyscr, "%%Q_IDX_2%%", tl2::var_to_str(Q_idx_2, g_prec));

	ofstr << pyscr << std::endl;
}
