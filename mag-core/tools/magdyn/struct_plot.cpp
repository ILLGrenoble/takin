/**
 * magnetic dynamics -- magnetic structure plotting
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
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

#include "struct_plot.h"
#include "helper.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>

#include <unordered_set>
#include <functional>

#include <boost/functional/hash.hpp>
#include <boost/scope_exit.hpp>

using namespace tl2_ops;



/**
 * shows the 3d view of the magnetic structure
 */
StructPlotDlg::StructPlotDlg(QWidget *parent, QSettings *sett, InfoDlg *info)
	: QDialog{parent}, m_sett{sett}, m_info_dlg{info}
{
	setWindowTitle("Magnetic Structure");
	setSizeGripEnabled(true);

	// create gl plotter
	m_structplot = new tl2::GlPlot(this);
	m_structplot->GetRenderer()->SetRestrictCamTheta(false);
	m_structplot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 5, 5, 5 }));
	m_structplot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -5, -5, -5 }));
	m_structplot->GetRenderer()->SetCoordMax(1.);
	 m_structplot->GetRenderer()->GetCamera().SetParalellRange(4.);
	m_structplot->GetRenderer()->GetCamera().SetFOV(tl2::d2r<t_real>(g_structplot_fov));
	m_structplot->GetRenderer()->GetCamera().SetDist(1.5);
	m_structplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_structplot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	m_coordcross = new QCheckBox("Show Coordinates", this);
	m_coordcross->setToolTip("Show the coordinate system cross.");
	m_coordcross->setChecked(false);

	m_labels = new QCheckBox("Show Labels", this);
	m_labels->setToolTip("Show magnetic site and coupling labels.");
	m_labels->setChecked(false);

	m_perspective = new QCheckBox("Perspective Projection", this);
	m_perspective->setToolTip("Switch between perspective and parallel projection.");
	m_perspective->setChecked(true);

	QPushButton *btn_100 = new QPushButton("[100] View", this);
	QPushButton *btn_010 = new QPushButton("[010] View", this);
	QPushButton *btn_001 = new QPushButton("[001] View", this);
	btn_100->setToolTip("View along [100] axis.");
	btn_010->setToolTip("View along [010] axis.");
	btn_001->setToolTip("View along [001] axis.");

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

	m_coordsys = new QComboBox(this);
	m_coordsys->addItem("Fractional Units (rlu)");
	m_coordsys->addItem("Lab Units (\xe2\x84\xab)");
	m_coordsys->setCurrentIndex(0);
	m_coordsys->setEnabled(false);

	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
	m_status->setFrameShape(QFrame::Panel);
	m_status->setFrameShadow(QFrame::Sunken);

	QPushButton *btnOk = new QPushButton("OK", this);

	// general context menu
	m_context = new QMenu(this);
	QAction *acCentre = new QAction("Centre Camera", m_context);
	m_context->addAction(acCentre);

	// context menu for sites
	m_context_site = new QMenu(this);
	QAction *acDelSite = new QAction("Delete Site", m_context_site);
	QAction *acFlipSpin = new QAction("Flip Spin", m_context_site);
	QAction *acCentreOnObject = new QAction("Centre Camera on Object", m_context_site);
	m_context_site->addAction(acDelSite);
	m_context_site->addAction(acFlipSpin);
	m_context_site->addSeparator();
	m_context_site->addAction(acCentre);
	m_context_site->addAction(acCentreOnObject);

	// context menu for terms
	m_context_term = new QMenu(this);
	QAction *acDelTerm = new QAction("Delete Coupling", m_context_term);
	m_context_term->addAction(acDelTerm);
	m_context_term->addSeparator();
	m_context_term->addAction(acCentre);
	m_context_term->addAction(acCentreOnObject);

	int y = 0;
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);
	grid->addWidget(m_structplot, y++, 0, 1, 6);
	grid->addWidget(m_coordcross, y, 0, 1, 2);
	grid->addWidget(m_labels, y, 2, 1, 2);
	grid->addWidget(m_perspective, y++, 4, 1, 2);
	grid->addWidget(btn_100, y, 0, 1, 2);
	grid->addWidget(btn_010, y, 2, 1, 2);
	grid->addWidget(btn_001, y++, 4, 1, 2);
	grid->addWidget(new QLabel("Camera Angles:", this), y, 0, 1, 2);
	grid->addWidget(m_cam_phi, y, 2, 1, 2);
	grid->addWidget(m_cam_theta, y++, 4, 1, 2);
	grid->addWidget(new QLabel("Coordinate System:", this), y, 0, 1, 2);
	grid->addWidget(m_coordsys, y++, 2, 1, 4);
	grid->addWidget(m_status, y, 0, 1, 5);
	grid->addWidget(btnOk, y++, 5, 1, 1);

	connect(btnOk, &QAbstractButton::clicked, this, &StructPlotDlg::accept);
	connect(m_structplot, &tl2::GlPlot::AfterGLInitialisation,
		this, &StructPlotDlg::AfterGLInitialisation);
	connect(m_structplot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
		this, &StructPlotDlg::PickerIntersection);
	connect(m_structplot->GetRenderer(), &tl2::GlPlotRenderer::CameraHasUpdated,
		this, &StructPlotDlg::CameraHasUpdated);
	connect(m_structplot, &tl2::GlPlot::MouseClick, this, &StructPlotDlg::MouseClick);
	connect(m_structplot, &tl2::GlPlot::MouseDown, this, &StructPlotDlg::MouseDown);
	connect(m_structplot, &tl2::GlPlot::MouseUp, this, &StructPlotDlg::MouseUp);
	connect(acFlipSpin, &QAction::triggered, this, &StructPlotDlg::FlipSpin);
	connect(acDelSite, &QAction::triggered, this, &StructPlotDlg::DeleteItem);
	connect(acDelTerm, &QAction::triggered, this, &StructPlotDlg::DeleteItem);
	connect(acCentreOnObject, &QAction::triggered, this, &StructPlotDlg::CentreCameraOnObject);
	connect(acCentre, &QAction::triggered, this, &StructPlotDlg::CentreCamera);
	connect(m_coordcross, &QCheckBox::toggled, this, &StructPlotDlg::ShowCoordCross);
	connect(m_labels, &QCheckBox::toggled, this, &StructPlotDlg::ShowLabels);
	connect(m_perspective, &QCheckBox::toggled, this, &StructPlotDlg::SetPerspectiveProjection);
	connect(m_coordsys, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
		this, &StructPlotDlg::SetCoordinateSystem);
	connect(btn_100, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(t_real_gl(90.), -t_real_gl(90.));
	});
	connect(btn_010, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(0., -t_real_gl(90.));
	});
	connect(btn_001, &QAbstractButton::clicked, [this]
	{
		this->SetCameraRotation(0., t_real_gl(180.));
	});

	connect(m_cam_phi,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl phi) -> void
	{
		this->SetCameraRotation(phi, m_cam_theta->value());
	});

	connect(m_cam_theta,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		[this](t_real_gl theta) -> void
	{
		this->SetCameraRotation(m_cam_phi->value(), theta);
	});

	if(m_sett && m_sett->contains("struct_view/geo"))
		restoreGeometry(m_sett->value("struct_view/geo").toByteArray());
	else
		resize(800, 800);
}



/**
 * dialog is closing
 */
void StructPlotDlg::accept()
{
	if(m_sett)
		m_sett->setValue("struct_view/geo", saveGeometry());

	QDialog::accept();
}



/**
 * structure plot picker intersection
 */
void StructPlotDlg::PickerIntersection(
        const t_vec3_gl* pos, std::size_t objIdx,
        [[maybe_unused]] const t_vec3_gl* posSphere)
{
	m_status->setText("");
	m_cur_obj = std::nullopt;
	m_cur_site = std::nullopt;
	m_cur_term = std::nullopt;

	m_structplot->GetRenderer()->SetObjectsHighlight(false);

	if(!pos)
		return;

	m_cur_obj = objIdx;

	// look for magnetic sites
	if(auto iter_sites = m_sites.find(objIdx);
		iter_sites != m_sites.end())
	{
		m_cur_site = iter_sites->second.name;
		HighlightSite(*m_cur_site);

		std::ostringstream ostr;
		ostr.precision(g_prec_gui);
		ostr << "Site: " << *m_cur_site
			<< " (position: " << iter_sites->second.sc_pos;
		if(!tl2::equals<t_vec_real>(iter_sites->second.uc_pos,
			iter_sites->second.sc_pos, g_eps))
			ostr << ", unit cell: " << iter_sites->second.uc_pos;
		ostr << ")";

		m_status->setText(ostr.str().c_str());
	}

	// look for exchange terms
	else if(auto iter_terms = m_terms.find(objIdx);
		iter_terms != m_terms.end())
	{
		m_cur_term = iter_terms->second.name;
		HighlightTerm(*m_cur_term);

		std::ostringstream ostr;
		ostr.precision(g_prec_gui);
		ostr << "Coupling: " << *m_cur_term
			<< " (length: " << iter_terms->second.length
			<< " \xe2\x84\xab)";

		m_status->setText(ostr.str().c_str());
	}
}



void StructPlotDlg::HighlightSite(const std::string& name)
{
	bool highlighted = false;

	// highlight all gl objects that are part of this site
	for(const auto& pair : m_sites)
	{
		if(pair.second.name == name)
		{
			m_structplot->GetRenderer()->SetObjectHighlight(pair.first, true);
			highlighted = true;
		}
	}

	if(highlighted)
		m_structplot->update();
}



void StructPlotDlg::HighlightTerm(const std::string& name)
{
	bool highlighted = false;

	// highlight all gl objects that are part of this term
	for(const auto& pair : m_terms)
	{
		if(pair.second.name == name)
		{
			m_structplot->GetRenderer()->SetObjectHighlight(pair.first, true);
			highlighted = true;
		}
	}

	if(highlighted)
		m_structplot->update();
}



/**
 * delete currently selected magnetic site or coupling
 */
void StructPlotDlg::DeleteItem()
{
	if(m_cur_site)
	{
		emit DeleteSite(*m_cur_site);
		m_cur_site = std::nullopt;
	}
	else if(m_cur_term)
	{
		emit DeleteTerm(*m_cur_term);
		m_cur_term = std::nullopt;
	}
}



/**
 * invert the currently selected site's spin
 */
void StructPlotDlg::FlipSpin()
{
	if(m_cur_site)
		emit FlipSiteSpin(*m_cur_site);
}



/**
 * show or hide the coordinate system
 */
void StructPlotDlg::ShowCoordCross(bool show)
{
	if(auto obj = m_structplot->GetRenderer()->GetCoordCross(); obj)
	{
		m_structplot->GetRenderer()->SetObjectVisible(*obj, show);
		m_structplot->update();
	}
}



/**
 * show or hide the object labels
 */
void StructPlotDlg::ShowLabels(bool show)
{
	m_structplot->GetRenderer()->SetLabelsVisible(show);
	m_structplot->update();
}



/**
 * choose between perspective or parallel projection
 */
void StructPlotDlg::SetPerspectiveProjection(bool proj)
{
	m_structplot->GetRenderer()->GetCamera().SetPerspectiveProjection(proj);
	m_structplot->GetRenderer()->RequestViewportUpdate();
	m_structplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_structplot->update();
}



/**
 * sets the camera's rotation angles
 */
void StructPlotDlg::SetCameraRotation(t_real_gl phi, t_real_gl theta)
{
	phi = tl2::d2r<t_real>(phi);
	theta = tl2::d2r<t_real>(theta);

	m_structplot->GetRenderer()->GetCamera().SetRotation(phi, theta);
	m_structplot->GetRenderer()->GetCamera().UpdateTransformation();
	CameraHasUpdated();
	m_structplot->update();
}



/**
 * the camera's properties have been updated
 */
void StructPlotDlg::CameraHasUpdated()
{
	auto [phi, theta] = m_structplot->GetRenderer()->GetCamera().GetRotation();

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
 * switch between crystal and lab coordinates
 */
void StructPlotDlg::SetCoordinateSystem(int which)
{
	m_structplot->GetRenderer()->SetCoordSys(which);
}



/**
 * centre camera on currently selected object
 */
void StructPlotDlg::CentreCameraOnObject()
{
	if(!m_cur_obj)
		return;

	const t_mat_gl& mat = m_structplot->GetRenderer()->GetObjectMatrix(*m_cur_obj);
	m_structplot->GetRenderer()->GetCamera().Centre(mat);
	m_structplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_structplot->update();
}



/**
 * centre camera on central position
 */
void StructPlotDlg::CentreCamera()
{
	t_mat_gl matCentre = tl2::hom_translation<t_mat_gl>(m_centre[0], m_centre[1], m_centre[2]);
	m_structplot->GetRenderer()->GetCamera().Centre(matCentre);
	m_structplot->GetRenderer()->GetCamera().UpdateTransformation();
	m_structplot->update();
}



/**
 * structure plot mouse button clicked
 */
void StructPlotDlg::MouseClick(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(right)
	{
		const QPointF& _pt = m_structplot->GetRenderer()->GetMousePosition();
		QPoint pt = m_structplot->mapToGlobal(_pt.toPoint());

		if(m_cur_site)
			m_context_site->popup(pt);
		else if(m_cur_term)
			m_context_term->popup(pt);
		else
			m_context->popup(pt);
	}
}



/**
 * structure plot mouse button pressed
 */
void StructPlotDlg::MouseDown(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(left && m_cur_site)
		emit SelectSite(*m_cur_site);

	if(left && m_cur_term)
		emit SelectTerm(*m_cur_term);
}



/**
 * structure plot mouse button released
 */
void StructPlotDlg::MouseUp(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}



/**
 * after structure plot initialisation
 */
void StructPlotDlg::AfterGLInitialisation()
{
	if(!m_structplot)
		return;

	// reference sphere for linked objects
	m_sphere = m_structplot->GetRenderer()->AddSphere(
		g_structplot_site_rad, 0.,0.,0., 1.,1.,1.,1.);
	m_structplot->GetRenderer()->SetObjectVisible(m_sphere, false);

	// reference arrow for linked objects
	m_arrow = m_structplot->GetRenderer()->AddArrow(
		g_structplot_dmi_rad, g_structplot_dmi_len, 0.,0.,0.5,  1.,1.,1.,1.);
	m_structplot->GetRenderer()->SetObjectVisible(m_arrow, false);

	// reference cylinder for linked objects
	m_cyl = m_structplot->GetRenderer()->AddCylinder(
		g_structplot_term_rad, 1., 0.,0.,0.5,  1.,1.,1.,1.);
	m_structplot->GetRenderer()->SetObjectVisible(m_cyl, false);

	// GL device info
	if(m_info_dlg)
	{
		auto [strGlVer, strGlShaderVer, strGlVendor, strGlRenderer]
			= m_structplot->GetRenderer()->GetGlDescr();

		m_info_dlg->SetGlInfo(0,
			QString("GL Version: %1.").arg(strGlVer.c_str()));
		m_info_dlg->SetGlInfo(1,
			QString("GL Shader Version: %1.").arg(strGlShaderVer.c_str()));
		m_info_dlg->SetGlInfo(2,
			QString("GL Vendor: %1.").arg(strGlVendor.c_str()));
		m_info_dlg->SetGlInfo(3,
			QString("GL Device: %1.").arg(strGlRenderer.c_str()));
	}

	ShowCoordCross(m_coordcross->isChecked());
	ShowLabels(m_labels->isChecked());
	SetPerspectiveProjection(m_perspective->isChecked());
	SetCoordinateSystem(m_coordsys->currentIndex());
	CameraHasUpdated();

	Sync();
}



/**
 * get the sites and exchange terms and
 * transfer them to the structure plotter
 */
void StructPlotDlg::Sync()
{
	if(!m_structplot || !m_dyn)
		return;

	// get sites and terms
	const auto& sites = m_dyn->GetMagneticSites();
	const auto& terms = m_dyn->GetExchangeTerms();
	const auto& field = m_dyn->GetExternalField();
	const auto& ordering = m_dyn->GetOrderingWavevector();
	const auto& rotaxis = m_dyn->GetRotationAxis();
	const bool is_incommensurate = m_dyn->IsIncommensurate();

	// clear old magnetic sites
	for(const auto& [site_idx, site] : m_sites)
		m_structplot->GetRenderer()->RemoveObject(site_idx);
	m_sites.clear();

	// clear old terms
	for(const auto& [term_idx, term] : m_terms)
		m_structplot->GetRenderer()->RemoveObject(term_idx);
	m_terms.clear();

	// clear centre position for camera
	m_centre = tl2::zero<t_vec_gl>(3);
	t_size total_sites = 0;

	// crystal matrix
	t_mat_gl matA = tl2::convert<t_mat_gl>(m_dyn->GetCrystalATrafo());
	t_mat_gl matB = tl2::convert<t_mat_gl>(m_dyn->GetCrystalBTrafo());

	m_structplot->GetRenderer()->SetBTrafo(matB, &matA);


	// hashes of already seen magnetic sites
	std::unordered_set<t_size> site_hashes;

	// calculate the hash of a magnetic site
	auto get_site_hash = [](const t_magdyn::MagneticSite& site,
		const t_vec_real* sc_dist = nullptr) -> t_size
	{
		// super cell index
		int _sc_x = 0, _sc_y = 0, _sc_z = 0;
		if(sc_dist)
		{
			_sc_x = int(std::round((*sc_dist)[0]));
			_sc_y = int(std::round((*sc_dist)[1]));
			_sc_z = int(std::round((*sc_dist)[2]));
		}

		t_size hash = 0;
		boost::hash_combine(hash, std::hash<std::string>{}(site.name));
		boost::hash_combine(hash, std::hash<int>{}(_sc_x));
		boost::hash_combine(hash, std::hash<int>{}(_sc_y));
		boost::hash_combine(hash, std::hash<int>{}(_sc_z));
		return hash;
	};


	// check if the magnetic site has already been seen
	auto site_not_yet_seen = [&site_hashes, &get_site_hash](
		const t_magdyn::MagneticSite& site,
		const t_vec_real& sc_dist)
	{
		t_size hash = get_site_hash(site, &sc_dist);
		return site_hashes.find(hash) == site_hashes.end();
	};


	// add a magnetic site to the plot
	auto add_site = [this, &site_hashes, &get_site_hash,
		is_incommensurate, &ordering, &rotaxis, &total_sites](
			t_size site_idx,
			const t_magdyn::MagneticSite& site,
			const t_magdyn::ExternalField& field,
			const t_vec_real* sc_dist = nullptr)
	{
		// super cell index
		int _sc_x = 0, _sc_y = 0, _sc_z = 0;
		if(sc_dist)
		{
			_sc_x = int(std::round((*sc_dist)[0]));
			_sc_y = int(std::round((*sc_dist)[1]));
			_sc_z = int(std::round((*sc_dist)[2]));
		}

		// default colour for unit cell magnetic sites
		t_real_gl rgb[3] { 0., 0., 1. };

		// get user-defined colour
		bool user_col = false;
		if(m_sitestab && site_idx < t_size(m_sitestab->rowCount()))
		{
			user_col = get_colour<t_real_gl>(
				m_sitestab->item(site_idx, COL_SITE_RGB)->
					text().toStdString(), rgb);
		}

		// no user-defined colour -> use default for super-cell magnetic sites
		if(!user_col)
		{
			if(_sc_x == 0 && _sc_y == 0 && _sc_z == 0)
			{
				rgb[0] = t_real_gl(1.);
				rgb[1] = t_real_gl(0.);
				rgb[2] = t_real_gl(0.);
			}
		}

		t_real_gl scale = 1.;

		std::size_t obj = m_structplot->GetRenderer()->AddLinkedObject(
			m_sphere, 0,0,0, rgb[0], rgb[1], rgb[2], 1);

		std::size_t arrow = m_structplot->GetRenderer()->AddLinkedObject(
			m_arrow, 0,0,0, rgb[0], rgb[1], rgb[2], 1);

		t_vec_real _pos_vec = site.pos_calc;
		if(sc_dist)
			_pos_vec += *sc_dist;
		t_vec_gl pos_vec = tl2::convert<t_vec_gl>(_pos_vec);
		t_vec_gl spin_vec;

		{
			MagneticSiteInfo siteinfo;
			siteinfo.name = site.name;
			siteinfo.uc_pos = site.pos_calc;
			siteinfo.sc_pos = _pos_vec;

			m_sites.emplace(std::make_pair(obj, siteinfo));
			m_sites.emplace(std::make_pair(arrow, std::move(siteinfo)));
		}

		// align spin to external field?
		if(field.align_spins)
		{
			spin_vec = tl2::create<t_vec_gl>({
				t_real_gl(-field.dir[0] * site.spin_mag_calc),
				t_real_gl(-field.dir[1] * site.spin_mag_calc),
				t_real_gl(-field.dir[2] * site.spin_mag_calc),
			});
		}
		else
		{
			spin_vec = tl2::create<t_vec_gl>({
				t_real_gl(site.spin_dir_calc[0] * site.spin_mag_calc),
				t_real_gl(site.spin_dir_calc[1] * site.spin_mag_calc),
				t_real_gl(site.spin_dir_calc[2] * site.spin_mag_calc),
			});

			if(is_incommensurate)
			{
				// rotate spin vector for incommensurate structures
				t_vec_gl sc_vec;
				if(sc_dist)
					sc_vec = tl2::convert<t_vec_gl>(*sc_dist);
				else
					sc_vec = tl2::zero<t_vec_gl>(3);

				tl2_mag::rotate_spin_incommensurate<t_mat_gl, t_vec_gl, t_real_gl>(
					spin_vec, sc_vec,
					tl2::convert<t_vec_gl>(ordering),
					tl2::convert<t_vec_gl>(rotaxis),
					g_eps);
			}
		}

		m_structplot->GetRenderer()->SetObjectMatrix(obj,
			tl2::hom_translation<t_mat_gl>(
				pos_vec[0], pos_vec[1], pos_vec[2]) *
			tl2::hom_scaling<t_mat_gl>(scale, scale, scale));

		m_structplot->GetRenderer()->SetObjectMatrix(arrow,
			tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
				spin_vec,                           // to
				1,                                  // post-scale
				tl2::create<t_vec_gl>({ 0, 0, 0 }), // post-translate
				tl2::create<t_vec_gl>({ 0, 0, 1 }), // from
				scale,                              // pre-scale
				pos_vec));                          // pre-translate

		m_structplot->GetRenderer()->SetObjectLabel(obj, site.name);
		//m_structplot->GetRenderer()->SetObjectLabel(arrow, site.name);

		// mark the magnetic site as already seen
		t_size hash = get_site_hash(site, sc_dist);
		site_hashes.insert(hash);

		m_centre += pos_vec;
		++total_sites;
	};


	// iterate and add unit cell magnetic sites
	for(t_size site_idx = 0; site_idx < sites.size(); ++site_idx)
		add_site(site_idx, sites[site_idx], field);


	// iterate and add exchange terms
	for(t_size term_idx = 0; term_idx < terms.size(); ++term_idx)
	{
		const t_term& term = terms[term_idx];
		if(term.site1_calc >= sites.size() || term.site2_calc >= sites.size())
			continue;

		// TODO: handle self-couplings (e.g. single-ion anisotropy)
		if(term.site1_calc == term.site2_calc &&
			tl2::equals_0<t_vec_real>(term.dist_calc, g_eps))
			continue;

		const t_site& site1 = sites[term.site1_calc];
		const t_site& site2 = sites[term.site2_calc];

		const t_vec_real& sc_dist = term.dist_calc;

		// get colour
		t_real_gl rgb[3] {0., 0.75, 0.};
		if(m_termstab && term_idx < t_size(m_termstab->rowCount()))
		{
			get_colour<t_real_gl>(
				m_termstab->item(term_idx, COL_XCH_RGB)->text().toStdString(), rgb);
		}

		t_real_gl scale = 1.;

		std::size_t obj = m_structplot->GetRenderer()->AddLinkedObject(
			m_cyl, 0, 0, 0, rgb[0], rgb[1], rgb[2], 1);

		{
			ExchangeTermInfo terminfo;
			terminfo.name = term.name;
			terminfo.length = term.length_calc;

			m_terms.emplace(std::make_pair(obj, std::move(terminfo)));
		}

		// connection from unit cell magnetic site...
		const t_vec_gl pos1_vec = tl2::convert<t_vec_gl>(site1.pos_calc);

		// ... to magnetic site in super cell
		const t_vec_gl pos2_vec = tl2::convert<t_vec_gl>(site2.pos_calc + sc_dist);

		// add the supercell site if it hasn't been inserted yet
		if(site_not_yet_seen(site2, sc_dist))
			add_site(term.site2_calc, site2, field, &sc_dist);

		t_vec_gl dir_vec = pos2_vec - pos1_vec;
		t_real_gl dir_len = tl2::norm<t_vec_gl>(dir_vec);
		t_vec_gl zero_vec = tl2::create<t_vec_gl>({ 0, 0, 0 });
		t_vec_gl z_vec = tl2::create<t_vec_gl>({ 0, 0, 1 });

		// coupling bond
		m_structplot->GetRenderer()->SetObjectMatrix(obj,
			tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
				dir_vec, 1, zero_vec,    // to, post-scale and post-translate
				z_vec, scale, pos1_vec)  // from, pre-scale and pre-translate
			* tl2::hom_translation<t_mat_gl>(
				t_real_gl(0), t_real_gl(0), dir_len*t_real_gl(0.5))
			* tl2::hom_scaling<t_mat_gl>(
				t_real_gl(1), t_real_gl(1), dir_len));

		m_structplot->GetRenderer()->SetObjectLabel(obj, term.name);


		// dmi vector
		t_vec_gl dmi_vec = tl2::zero<t_vec_gl>(3);
		if(term.dmi_calc.size() >= 3)
		{
			dmi_vec[0] = t_real_gl(term.dmi_calc[0].real());
			dmi_vec[1] = t_real_gl(term.dmi_calc[1].real());
			dmi_vec[2] = t_real_gl(term.dmi_calc[2].real());
		}

		if(tl2::norm<t_vec_gl>(dmi_vec) > g_eps)
		{
			std::size_t objDmi = m_structplot->GetRenderer()->AddLinkedObject(
				m_arrow, 0,0,0, rgb[0], rgb[1], rgb[2], 1);

			{
				ExchangeTermInfo terminfo;
				terminfo.name = term.name;
				terminfo.length = term.length_calc;

				m_terms.emplace(std::make_pair(
					objDmi, std::move(terminfo)));
			}

			t_real_gl scale_dmi = 0.5;

			m_structplot->GetRenderer()->SetObjectMatrix(objDmi,
				tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
					dmi_vec, 1, zero_vec,                // to, post-scale and post-translate
					z_vec, scale_dmi,                    // from and pre-scale
					(pos1_vec + pos2_vec)/t_real_gl(2))  // pre-translate
				);

			//m_structplot->GetRenderer()->SetObjectLabel(objDmi, term.name);
		}
	} // terms

	m_centre /= t_real_gl(total_sites);

	CentreCamera();
	m_structplot->update();
}
