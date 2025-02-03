/**
 * Calculation of polarisation vector
 * @author Tobias Weber <tweber@ill.fr>
 * @date Oct-2018
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "magtools" project
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
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

#include "pol.h"
#include "../defs.h"

#include <QtCore/QDir>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QTabWidget>

#include <locale>
#include <iostream>
#include <optional>

#include <boost/version.hpp>
#include <boost/config.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;


using namespace tl2_ops;

using t_vec = tl2::vec<t_cplx, std::vector>;
using t_mat = tl2::mat<t_cplx, std::vector>;
using t_matvec = std::vector<t_mat>;

using t_real_gl = tl2::t_real_gl;
using t_vec2_gl = tl2::t_vec2_gl;
using t_vec3_gl = tl2::t_vec3_gl;
using t_vec_gl = tl2::t_vec_gl;
using t_mat_gl = tl2::t_mat_gl;


// ----------------------------------------------------------------------------
void PolDlg::closeEvent(QCloseEvent *)
{
	// save window size and position
	m_sett->setValue("pol/geo", saveGeometry());

	// save values
	m_sett->setValue("pol/n_re", m_editNRe->text().toDouble());
	m_sett->setValue("pol/n_im", m_editNIm->text().toDouble());
	m_sett->setValue("pol/mx_re", m_editMPerpReX->text().toDouble());
	m_sett->setValue("pol/my_re", m_editMPerpReY->text().toDouble());
	m_sett->setValue("pol/mz_re", m_editMPerpReZ->text().toDouble());
	m_sett->setValue("pol/mx_im", m_editMPerpImX->text().toDouble());
	m_sett->setValue("pol/my_im", m_editMPerpImY->text().toDouble());
	m_sett->setValue("pol/mz_im", m_editMPerpImZ->text().toDouble());
	m_sett->setValue("pol/pix", m_editPiX->text().toDouble());
	m_sett->setValue("pol/piy", m_editPiY->text().toDouble());
	m_sett->setValue("pol/piz", m_editPiZ->text().toDouble());
}


/**
 * called after the plotter has initialised
 */
void PolDlg::AfterGLInitialisation()
{
	// GL device info
	auto [ strGlVer, strGlShaderVer, strGlVendor, strGlRenderer ]
		= m_plot->GetRenderer()->GetGlDescr();
	emit GlDeviceInfos(strGlVer, strGlShaderVer, strGlVendor, strGlRenderer);

	// create 3d objects
	if(!m_3dobjsReady)
	{
		m_arrow_pi = m_plot->GetRenderer()->AddArrow(0.05, 1., 0.,0.,0.5,  0.,0.,0.85,1.);
		m_arrow_pf = m_plot->GetRenderer()->AddArrow(0.05, 1., 0.,0.,0.5,  0.,0.5,0.,1.);
		m_arrow_M_Re = m_plot->GetRenderer()->AddArrow(0.05, 1., 0.,0.,0.5,  0.85,0.,0.,1.);
		m_arrow_M_Im = m_plot->GetRenderer()->AddArrow(0.05, 1., 0.,0.,0.5,  0.85,0.25,0.,1.);

		m_plot->GetRenderer()->SetObjectLabel(m_arrow_pi, "P_i");
		m_plot->GetRenderer()->SetObjectLabel(m_arrow_pf, "P_f");
		m_plot->GetRenderer()->SetObjectLabel(m_arrow_M_Re, "Re{M_perp}");
		m_plot->GetRenderer()->SetObjectLabel(m_arrow_M_Im, "Im{M_perp}");

		m_3dobjsReady = true;
		CalcPol();
	}
}


/**
 * get the length of a vector
 */
t_real PolDlg::GetArrowLen(std::size_t objIdx) const
{
	if(objIdx == m_arrow_pi)
	{
		return std::sqrt(std::pow(m_editPiX->text().toDouble(), 2.)
			+ std::pow(m_editPiY->text().toDouble(), 2.)
			+ std::pow(m_editPiZ->text().toDouble(), 2.));
	}
	else if(objIdx == m_arrow_pf)
	{
		return std::sqrt(std::pow(m_editPfX->text().toDouble(), 2.)
			+ std::pow(m_editPfY->text().toDouble(), 2.)
			+ std::pow(m_editPfZ->text().toDouble(), 2.));
	}
	else if(objIdx == m_arrow_M_Re)
	{
		return std::sqrt(std::pow(m_editMPerpReX->text().toDouble(), 2.)
			+ std::pow(m_editMPerpReY->text().toDouble(), 2.)
			+ std::pow(m_editMPerpReZ->text().toDouble(), 2.));
	}
	else if(objIdx == m_arrow_M_Im)
	{
		return std::sqrt(std::pow(m_editMPerpImX->text().toDouble(), 2.)
			+ std::pow(m_editMPerpImY->text().toDouble(), 2.)
			+ std::pow(m_editMPerpImZ->text().toDouble(), 2.));
	}

	return -1.;
}


/**
 * called when the mouse hovers over an object
 */
void PolDlg::PickerIntersection(const t_vec3_gl* pos,
	std::size_t objIdx, std::size_t /*triagIdx*/,
	const t_vec3_gl* posSphere)
{
	m_curPickedObj.reset();

	if(pos)
	{	// object selected?
		m_curPickedObj = objIdx;

		if(objIdx == m_arrow_pi)
			m_labelStatus->setText("P_i");
		else if(objIdx == m_arrow_pf)
			m_labelStatus->setText("P_f");
		else if(objIdx == m_arrow_M_Re)
			m_labelStatus->setText("Re{M_perp}");
		else if(objIdx == m_arrow_M_Im)
			m_labelStatus->setText("Im{M_perp}");
		else
			m_curPickedObj.reset();
	}

	if(m_curPickedObj)
	{
		setCursor(Qt::CrossCursor);
	}
	else
	{
		m_labelStatus->setText("");
		setCursor(Qt::ArrowCursor);
	}


	if(posSphere && m_mouseDown[0] && m_curDraggedObj)
	{	// picker intersecting unit sphere and mouse dragged?

		t_vec3_gl posSph = *posSphere;

		if(*m_curDraggedObj == m_arrow_pi)
		{
			m_editPiX->setText(tl2::var_to_str(posSph[0], g_prec).c_str());
			m_editPiY->setText(tl2::var_to_str(posSph[1], g_prec).c_str());
			m_editPiZ->setText(tl2::var_to_str(posSph[2], g_prec).c_str());
			CalcPol();
		}
		else if(*m_curDraggedObj == m_arrow_M_Re)
		{
			m_editMPerpReX->setText(tl2::var_to_str(posSph[0], g_prec).c_str());
			m_editMPerpReY->setText(tl2::var_to_str(posSph[1], g_prec).c_str());
			m_editMPerpReZ->setText(tl2::var_to_str(posSph[2], g_prec).c_str());
			CalcPol();
		}
		else if(*m_curDraggedObj == m_arrow_M_Im)
		{
			m_editMPerpImX->setText(tl2::var_to_str(posSph[0], g_prec).c_str());
			m_editMPerpImY->setText(tl2::var_to_str(posSph[1], g_prec).c_str());
			m_editMPerpImZ->setText(tl2::var_to_str(posSph[2], g_prec).c_str());
			CalcPol();
		}
	}
}


/**
 * mouse button pressed
 */
void PolDlg::MouseDown(bool left, bool mid, bool right)
{
	if(left)
		m_mouseDown[0] = true;
	if(mid)
		m_mouseDown[1] = true;
	if(right)
		m_mouseDown[2] = true;

	if(m_mouseDown[0])
	{
		if((m_curDraggedObj = m_curPickedObj))
		{
			auto lenVec = GetArrowLen(*m_curDraggedObj);
			if(lenVec > 0.)
				m_plot->GetRenderer()->SetPickerSphereRadius(lenVec);
		}
	}
}


/**
 * mouse button released
 */
void PolDlg::MouseUp(bool left, bool mid, bool right)
{
	if(left)
		m_mouseDown[0] = false;
	if(mid)
		m_mouseDown[1] = false;
	if(right)
		m_mouseDown[2] = false;

	if(!m_mouseDown[0])
		m_curDraggedObj.reset();
}


/**
 * create UI
 */
PolDlg::PolDlg(QWidget* pParent, QSettings *sett)
	: QDialog{pParent, Qt::Window}, m_sett{sett}
{
	setWindowTitle("Polarisation Vectors");
	setSizeGripEnabled(true);


	// plot panel
	auto plotpanel = new QWidget(this);

	auto labelN = new QLabel("Re{N}, Im{N}:", plotpanel);
	auto labelMPerpRe = new QLabel("Re{M_perp}:", plotpanel);
	auto labelMPerpIm = new QLabel("Im{M_perp}:", plotpanel);
	auto labelPi = new QLabel("P_i:", plotpanel);
	auto labelPf = new QLabel("P_f:", plotpanel);


	for(auto* label : {labelMPerpRe, labelMPerpIm, labelPi, labelPf})
		label->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	for(auto* editPf : {m_editPfX, m_editPfY, m_editPfZ})
		editPf->setReadOnly(true);

	m_labelStatus->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
	m_labelStatus->setFrameStyle(int(QFrame::Sunken) | int(QFrame::Panel));
	m_labelStatus->setLineWidth(1);


	// connections
	for(auto* edit : {m_editNRe, m_editNIm,
		m_editMPerpReX, m_editMPerpReY, m_editMPerpReZ,
		m_editMPerpImX, m_editMPerpImY, m_editMPerpImZ,
		m_editPiX, m_editPiY, m_editPiZ,
		m_editPfX, m_editPfY, m_editPfZ})
		connect(edit, &QLineEdit::textEdited, this, &PolDlg::CalcPol);

	connect(m_plot.get(), &tl2::GlPlot::AfterGLInitialisation, this, &PolDlg::AfterGLInitialisation);
	connect(m_plot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection, this, &PolDlg::PickerIntersection);

	connect(m_plot.get(), &tl2::GlPlot::MouseDown, this, &PolDlg::MouseDown);
	connect(m_plot.get(), &tl2::GlPlot::MouseUp, this, &PolDlg::MouseUp);


	auto pGrid = new QGridLayout(plotpanel);
	pGrid->setSpacing(0);
	pGrid->setContentsMargins(0, 0, 0, 0);

	pGrid->addWidget(m_plot.get(), 0,0,1,4);

	pGrid->addWidget(labelN, 1,0,1,1);
	pGrid->addWidget(labelMPerpRe, 2,0,1,1);
	pGrid->addWidget(labelMPerpIm, 3,0,1,1);
	pGrid->addWidget(labelPi, 4,0,1,1);
	pGrid->addWidget(labelPf, 5,0,1,1);

	pGrid->addWidget(m_editNRe, 1,1,1,1);
	pGrid->addWidget(m_editNIm, 1,2,1,1);

	pGrid->addWidget(m_editMPerpReX, 2,1,1,1);
	pGrid->addWidget(m_editMPerpReY, 2,2,1,1);
	pGrid->addWidget(m_editMPerpReZ, 2,3,1,1);
	pGrid->addWidget(m_editMPerpImX, 3,1,1,1);
	pGrid->addWidget(m_editMPerpImY, 3,2,1,1);
	pGrid->addWidget(m_editMPerpImZ, 3,3,1,1);

	pGrid->addWidget(m_editPiX, 4,1,1,1);
	pGrid->addWidget(m_editPiY, 4,2,1,1);
	pGrid->addWidget(m_editPiZ, 4,3,1,1);
	pGrid->addWidget(m_editPfX, 5,1,1,1);
	pGrid->addWidget(m_editPfY, 5,2,1,1);
	pGrid->addWidget(m_editPfZ, 5,3,1,1);

	pGrid->addWidget(m_labelStatus, 6,0,1,4);

	// restore last values
	if(m_sett->contains("pol/n_re"))
		m_editNRe->setText(m_sett->value("pol/n_re").toString());
	if(m_sett->contains("pol/n_im"))
		m_editNIm->setText(m_sett->value("pol/n_im").toString());
	if(m_sett->contains("pol/mx_re"))
		m_editMPerpReX->setText(m_sett->value("pol/mx_re").toString());
	if(m_sett->contains("pol/my_re"))
		m_editMPerpReY->setText(m_sett->value("pol/my_re").toString());
	if(m_sett->contains("pol/mz_re"))
		m_editMPerpReZ->setText(m_sett->value("pol/mz_re").toString());
	if(m_sett->contains("pol/mx_im"))
		m_editMPerpImX->setText(m_sett->value("pol/mx_im").toString());
	if(m_sett->contains("pol/my_im"))
		m_editMPerpImY->setText(m_sett->value("pol/my_im").toString());
	if(m_sett->contains("pol/mz_im"))
		m_editMPerpImZ->setText(m_sett->value("pol/mz_im").toString());
	if(m_sett->contains("pol/pix"))
		m_editPiX->setText(m_sett->value("pol/pix").toString());
	if(m_sett->contains("pol/piy"))
		m_editPiY->setText(m_sett->value("pol/piy").toString());
	if(m_sett->contains("pol/piz"))
		m_editPiZ->setText(m_sett->value("pol/piz").toString());



	auto pmainGrid = new QGridLayout(this);
	pmainGrid->setSpacing(4);
	pmainGrid->setContentsMargins(4, 4, 4, 4);
	pmainGrid->addWidget(plotpanel, 0, 0, 1, 1);


	// restory window size and position
	if(m_sett->contains("pol/geo"))
		restoreGeometry(m_sett->value("pol/geo").toByteArray());
	else
		resize(800, 800);

	// have scattering plane in horizontal plane
	m_plot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 5, 5, 5 }));
	m_plot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -5, -5, -5 }));
	m_plot->GetRenderer()->SetCoordMax(5.);
	m_plot->GetRenderer()->GetCamera().SetDist(2.5);
	m_plot->GetRenderer()->GetCamera().UpdateTransformation();

	CalcPol();
}


/**
 * calculate final polarisation vector
 */
void PolDlg::CalcPol()
{
	// get values from line edits
	t_real NRe = t_real(m_editNRe->text().toDouble());
	t_real NIm = t_real(m_editNIm->text().toDouble());

	t_real MPerpReX = t_real(m_editMPerpReX->text().toDouble());
	t_real MPerpReY = t_real(m_editMPerpReY->text().toDouble());
	t_real MPerpReZ = t_real(m_editMPerpReZ->text().toDouble());
	t_real MPerpImX = t_real(m_editMPerpImX->text().toDouble());
	t_real MPerpImY = t_real(m_editMPerpImY->text().toDouble());
	t_real MPerpImZ = t_real(m_editMPerpImZ->text().toDouble());

	t_real PiX = t_real(m_editPiX->text().toDouble());
	t_real PiY = t_real(m_editPiY->text().toDouble());
	t_real PiZ = t_real(m_editPiZ->text().toDouble());

	const t_cplx N(NRe, NIm);
	const t_vec Mperp = tl2::create<t_vec>({
		t_cplx(MPerpReX,MPerpImX),
		t_cplx(MPerpReY,MPerpImY),
		t_cplx(MPerpReZ,MPerpImZ) });
	const t_vec Pi = tl2::create<t_vec>({PiX, PiY, PiZ});

	// calculate final polarisation vector and intensity
	auto [I, P_f] = tl2::blume_maleev_indir<t_mat, t_vec, t_cplx>(Pi, Mperp, N);
	//auto [I, P_f] = tl2::blume_maleev<t_vec, t_cplx>(Pi, Mperp, N);

	// set final polarisation
	m_editPfX->setText(tl2::var_to_str(P_f[0].real(), g_prec).c_str());
	m_editPfY->setText(tl2::var_to_str(P_f[1].real(), g_prec).c_str());
	m_editPfZ->setText(tl2::var_to_str(P_f[2].real(), g_prec).c_str());


	// update 3d objects
	if(m_3dobjsReady)
	{
		// P_i
		t_mat_gl matPi = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
			tl2::create<t_vec_gl>({t_real_gl(PiX), t_real_gl(PiY), t_real_gl(PiZ)}),  // to
			1.,                                     // scale
			tl2::create<t_vec_gl>({0,0,0.5}),       // translate
			tl2::create<t_vec_gl>({0,0,1}));        // from
		m_plot->GetRenderer()->SetObjectMatrix(m_arrow_pi, matPi);

		// P_f
		t_mat_gl matPf = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
			tl2::create<t_vec_gl>({t_real_gl(P_f[0].real()), t_real_gl(P_f[1].real()), t_real_gl(P_f[2].real())}),  // to
			1.,                                     // scale
			tl2::create<t_vec_gl>({0,0,0.5}),       // translate
			tl2::create<t_vec_gl>({0,0,1}));        // from
		m_plot->GetRenderer()->SetObjectMatrix(m_arrow_pf, matPf);

		// Re(M)
		const t_real_gl lenReM = t_real_gl(std::sqrt(MPerpReX*MPerpReX + MPerpReY*MPerpReY + MPerpReZ*MPerpReZ));
		t_mat_gl matMRe = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
			tl2::create<t_vec_gl>({t_real_gl(MPerpReX), t_real_gl(MPerpReY), t_real_gl(MPerpReZ)}),  // to
			lenReM,                                 // scale
			tl2::create<t_vec_gl>({0,0,0.5}),       // translate
			tl2::create<t_vec_gl>({0,0,1}));        // from
		m_plot->GetRenderer()->SetObjectMatrix(m_arrow_M_Re, matMRe);
		m_plot->GetRenderer()->SetObjectVisible(m_arrow_M_Re, !tl2::equals(lenReM, t_real_gl(0)));

		// Im(M)
		const t_real_gl lenImM = t_real_gl(std::sqrt(MPerpImX*MPerpImX + MPerpImY*MPerpImY + MPerpImZ*MPerpImZ));
		t_mat_gl matMIm = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
			tl2::create<t_vec_gl>({t_real_gl(MPerpImX), t_real_gl(MPerpImY), t_real_gl(MPerpImZ)}),  // to
			lenImM,                                 // scale
			tl2::create<t_vec_gl>({0,0,0.5}),       // translate
			tl2::create<t_vec_gl>({0,0,1}));        // from
		m_plot->GetRenderer()->SetObjectMatrix(m_arrow_M_Im, matMIm);
		m_plot->GetRenderer()->SetObjectVisible(m_arrow_M_Im, !tl2::equals(lenImM, t_real_gl(0)));

		m_plot->update();
	}
}
// ----------------------------------------------------------------------------
