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
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <cstdlib>
#include <mutex>
#include <sstream>
#include <iomanip>

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QDialogButtonBox>

#include "dispersion3d.h"
#include "helper.h"



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
	m_dispplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 5, 5, 5 }));
	m_dispplot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -5, -5, -5 }));
	m_dispplot->GetRenderer()->SetCoordMax(1.);
	m_dispplot->GetRenderer()->GetCamera().SetParalellRange(4.);
	m_dispplot->GetRenderer()->GetCamera().SetFOV(tl2::d2r<t_real>(g_structplot_fov));
	m_dispplot->GetRenderer()->GetCamera().SetDist(1.5);
	m_dispplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_dispplot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// Q coordinates
	QGroupBox *groupQ = new QGroupBox("Q Coordinates", this);
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

		m_Q_dir1[i]->setDecimals(4);
		m_Q_dir1[i]->setMinimum(-99.9999);
		m_Q_dir1[i]->setMaximum(+99.9999);
		m_Q_dir1[i]->setSingleStep(0.01);
		m_Q_dir1[i]->setValue(i == 0 ? 1. : 0.);
		//m_Q_dir1[i]->setSuffix(" rlu");
		m_Q_dir1[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir1[i]->setPrefix(hklPrefix[i]);

		m_Q_dir2[i]->setDecimals(4);
		m_Q_dir2[i]->setMinimum(-99.9999);
		m_Q_dir2[i]->setMaximum(+99.9999);
		m_Q_dir2[i]->setSingleStep(0.01);
		m_Q_dir2[i]->setValue(i == 1 ? 1. : 0.);
		//m_Q_dir2[i]->setSuffix(" rlu");
		m_Q_dir2[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Preferred});
		m_Q_dir2[i]->setPrefix(hklPrefix[i]);
	}

	for(int i = 0; i < 2; ++i)
	{
		m_num_Q_points[i]->setMinimum(1);
		m_num_Q_points[i]->setMaximum(9999);
		m_num_Q_points[i]->setSingleStep(1);
		m_num_Q_points[i]->setValue(32);
	}

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
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	// close button
	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);

	// Q coordinates grid
	int y = 0;
	QGridLayout *Qgrid = new QGridLayout(groupQ);
	Qgrid->addWidget(new QLabel("Origin:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_origin[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_origin[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_origin[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Direction 1:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_dir1[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_dir1[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_dir1[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Direction 2:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_Q_dir2[0], y, 1, 1, 1);
	Qgrid->addWidget(m_Q_dir2[1], y, 2, 1, 1);
	Qgrid->addWidget(m_Q_dir2[2], y++, 3, 1, 1);
	Qgrid->addWidget(new QLabel("Number of Points:", this), y, 0, 1, 1);
	Qgrid->addWidget(m_num_Q_points[0], y, 1, 1, 1);
	Qgrid->addWidget(m_num_Q_points[1], y++, 2, 1, 1);

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
	maingrid->addWidget(m_dispplot, y++, 0, 1, 4);
	maingrid->addWidget(groupQ, y++, 0, 1, 4);
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
	}

	// connections
	connect(btnbox, &QDialogButtonBox::accepted, this, &Dispersion3DDlg::accept);

	connect(m_dispplot, &tl2::GlPlot::AfterGLInitialisation,
		this, &Dispersion3DDlg::AfterPlotGLInitialisation);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
		this, &Dispersion3DDlg::PlotPickerIntersection);
	connect(m_dispplot->GetRenderer(), &tl2::GlPlotRenderer::CameraHasUpdated,
		this, &Dispersion3DDlg::PlotCameraHasUpdated);
	connect(m_dispplot, &tl2::GlPlot::MouseClick, this, &Dispersion3DDlg::PlotMouseClick);
	connect(m_dispplot, &tl2::GlPlot::MouseDown, this, &Dispersion3DDlg::PlotMouseDown);
	connect(m_dispplot, &tl2::GlPlot::MouseUp, this, &Dispersion3DDlg::PlotMouseUp);

	// calculation
	connect(m_btn_start_stop, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calc_enabled)
			Calculate();
		else
			m_stop_requested = true;
	});

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



void Dispersion3DDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, windowTitle() + " -- Error", msg);
}



/**
 * calculate the dispersion
 */
void Dispersion3DDlg::Calculate()
{
	if(!m_dyn)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	// get coordinates
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

	t_size Q_count_1 = m_num_Q_points[0]->value();
	t_size Q_count_2 = m_num_Q_points[1]->value();

	t_vec_real Q_step_1 = Q_dir_1 / t_real(Q_count_1);
	t_vec_real Q_step_2 = Q_dir_2 / t_real(Q_count_2);

	bool use_weights = false;
	bool use_projector = true;

	// calculate the dispersion
	t_magdyn dyn = *m_dyn;
	//dyn.SetUniteDegenerateEnergies(false);

	// tread pool and mutex to protect the data vectors
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stop_requested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(Q_count_1 * Q_count_2);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 threads.").arg(g_num_threads));

	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count_1 * Q_count_2);

	for(t_size Q_idx_1 = 0; Q_idx_1 < Q_count_1; ++Q_idx_1)
	for(t_size Q_idx_2 = 0; Q_idx_2 < Q_count_2; ++Q_idx_2)
	{
		auto task = [this, &mtx, &dyn,
			Q_idx_1, Q_idx_2, Q_count_1, Q_count_2,
			&Q_origin, &Q_step_1, &Q_step_2,
			use_weights, use_projector]()
		{
			// calculate the dispersion at the given Q point
			t_vec_real Q = Q_origin + Q_step_1*t_real(Q_idx_1) + Q_step_2*t_real(Q_idx_2);
			auto E_and_S = dyn.CalcEnergies(Q[0], Q[1], Q[2], !use_weights);

			std::vector<t_real> Es, weights;
			Es.reserve(E_and_S.size());
			weights.reserve(E_and_S.size());

			// iterate the energies for this Q point
			for(const auto& E_and_S : E_and_S)
			{
				t_real E = E_and_S.E;
				if(std::isnan(E) || std::isinf(E))
					continue;

				t_real weight = -1;
				if(use_weights)
				{
					weight = E_and_S.weight;

					if(!use_projector)
					{
						const t_mat& S = E_and_S.S;
						weight = tl2::trace<t_mat>(S).real();
					}

					if(std::isnan(weight) || std::isinf(weight))
						weight = 0.;
				}

				// TODO
			}

			std::lock_guard<std::mutex> _lck{mtx};
			//m_data.emplace_back(std::move(data));
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
		if(m_stop_requested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stop_requested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());


	Plot();
}



/**
 * plot the calculated dispersion
 */
void Dispersion3DDlg::Plot()
{
	// TODO
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
		m_btn_start_stop->setToolTip("Start calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btn_start_stop->setText("Stop");
		m_btn_start_stop->setToolTip("Stop running calculation.");
		m_btn_start_stop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



/**
 * dispersion plot picker intersection
 */
void Dispersion3DDlg::PlotPickerIntersection(
	[[maybe_unused]] const t_vec3_gl* pos,
	[[maybe_unused]] std::size_t objIdx,
	[[maybe_unused]] std::size_t triagIdx,
	[[maybe_unused]] const t_vec3_gl* posSphere)
{
	m_dispplot->GetRenderer()->SetObjectsHighlight(false);
}



/**
 * the dispersion plot's camera properties have been updated
 */
void Dispersion3DDlg::PlotCameraHasUpdated()
{
	//auto [phi, theta] = m_dispplot->GetRenderer()->GetCamera().GetRotation();
}



/**
 * dispersion plot mouse button clicked
 */
void Dispersion3DDlg::PlotMouseClick(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	/*if(right)
	{
		const QPointF& _pt = m_dispplot->GetRenderer()->GetMousePosition();
		QPoint pt = m_dispplot->mapToGlobal(_pt.toPoint());

		m_context->popup(pt);
	}*/
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

	//ShowCoordCross(m_coordcross->isChecked());
	//ShowLabels(m_labels->isChecked());
	//SetPerspectiveProjection(m_perspective->isChecked());
	//SetCoordinateSystem(m_coordsys->currentIndex());
	PlotCameraHasUpdated();
}



/**
 * dialog is closing
 */
void Dispersion3DDlg::accept()
{
	if(m_sett)
	{
		m_sett->setValue("dispersion3d/geo", saveGeometry());
	}

	QDialog::accept();
}
