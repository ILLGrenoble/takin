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

#include <QtWidgets/QGridLayout>
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

	// main grid
	QGridLayout *maingrid = new QGridLayout(this);
	maingrid->setSpacing(4);
	maingrid->setContentsMargins(8, 8, 8, 8);
	maingrid->addWidget(m_dispplot, 0, 0, 1, 4);
	maingrid->addWidget(m_status, 1, 0, 1, 3);
	maingrid->addWidget(btnbox, 1, 3, 1, 1);

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
 * dispersion plot picker intersection
 */
void Dispersion3DDlg::PlotPickerIntersection(
	const t_vec3_gl* pos, std::size_t objIdx,
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
