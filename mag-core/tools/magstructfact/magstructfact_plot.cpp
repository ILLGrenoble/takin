/**
 * magnetic structure factor tool
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2019
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#include "magstructfact.h"

#include <iostream>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/qt/helper.h"

using namespace tl2_ops;


/**
 * add 3d object
 */
void MagStructFactDlg::Add3DItem(int row)
{
	if(!m_plot)
		return;

	// add all items
	if(row < 0)
	{
		for(int row=0; row<m_nuclei->rowCount(); ++row)
			Add3DItem(row);
		return;
	}

	auto objSphere = m_plot->GetRenderer()->AddLinkedObject(m_sphere, 0,0,0, 1,1,1,1);
	auto objArrowRe = m_plot->GetRenderer()->AddLinkedObject(m_arrow, 0,0,0, 1,1,1,1);
	auto objArrowIm = m_plot->GetRenderer()->AddLinkedObject(m_arrow, 0,0,0, 1,1,1,1);

	m_nuclei->item(row, COL_NAME)->setData(Qt::UserRole + 0,
		unsigned(objSphere));	// atomic position
	m_nuclei->item(row, COL_NAME)->setData(Qt::UserRole + 1,
		unsigned(objArrowRe));	// real part of Fourier comp
	m_nuclei->item(row, COL_NAME)->setData(Qt::UserRole + 2,
		unsigned(objArrowIm));	// imag part of Fourier comp

	Sync3DItem(row);
}


/**
 * sync the properties of a 3d object
 */
void MagStructFactDlg::Sync3DItem(int row)
{
	if(!m_plot)
		return;

	// sync all items
	if(row < 0)
	{
		for(int row=0; row<m_nuclei->rowCount(); ++row)
			Sync3DItem(row);
		return;
	}

	std::size_t objSphere = m_nuclei->item(row, COL_NAME)->
		data(Qt::UserRole + 0).toUInt();
	std::size_t objArrowRe = m_nuclei->item(row, COL_NAME)->
		data(Qt::UserRole + 1).toUInt();
	std::size_t objArrowIm = m_nuclei->item(row, COL_NAME)->
		data(Qt::UserRole + 2).toUInt();
	if(!objSphere || !objArrowRe || !objArrowIm)
		return;

	auto *itemName = m_nuclei->item(row, COL_NAME);
	auto *itemx = m_nuclei->item(row, COL_X);
	auto *itemy = m_nuclei->item(row, COL_Y);
	auto *itemz = m_nuclei->item(row, COL_Z);
	auto *itemM = m_nuclei->item(row, COL_M_MAG);
	auto *itemReMX = m_nuclei->item(row, COL_ReM_X);
	auto *itemReMY = m_nuclei->item(row, COL_ReM_Y);
	auto *itemReMZ = m_nuclei->item(row, COL_ReM_Z);
	auto *itemImMX = m_nuclei->item(row, COL_ImM_X);
	auto *itemImMY = m_nuclei->item(row, COL_ImM_Y);
	auto *itemImMZ = m_nuclei->item(row, COL_ImM_Z);
	auto *itemsc = m_nuclei->item(row, COL_RAD);
	auto *itemCol = m_nuclei->item(row, COL_COL);

	t_real_gl posx = tl2::stoval<t_real_gl>(itemx->text().toStdString());
	t_real_gl posy = tl2::stoval<t_real_gl>(itemy->text().toStdString());
	t_real_gl posz = tl2::stoval<t_real_gl>(itemz->text().toStdString());
	t_real_gl M = tl2::stoval<t_real_gl>(itemM->text().toStdString());
	t_real_gl ReMX = tl2::stoval<t_real_gl>(itemReMX->text().toStdString());
	t_real_gl ReMY = tl2::stoval<t_real_gl>(itemReMY->text().toStdString());
	t_real_gl ReMZ = tl2::stoval<t_real_gl>(itemReMZ->text().toStdString());
	t_real_gl ImMX = tl2::stoval<t_real_gl>(itemImMX->text().toStdString());
	t_real_gl ImMY = tl2::stoval<t_real_gl>(itemImMY->text().toStdString());
	t_real_gl ImMZ = tl2::stoval<t_real_gl>(itemImMZ->text().toStdString());
	t_real_gl scale = tl2::stoval<t_real_gl>(itemsc->text().toStdString());

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	qreal r=1, g=1, b=1;
#else
	float r=1, g=1, b=1;
#endif
	QColor col{itemCol->text()};
	col.getRgbF(&r, &g, &b);

	t_mat_gl matSphere = tl2::hom_translation<t_mat_gl>(posx, posy, posz) *
		tl2::hom_scaling<t_mat_gl>(M*scale, M*scale, M*scale);

	auto vecReM = tl2::create<t_vec_gl>({t_real_gl(ReMX), t_real_gl(ReMY), t_real_gl(ReMZ)});
	auto vecImM = tl2::create<t_vec_gl>({t_real_gl(ImMX), t_real_gl(ImMY), t_real_gl(ImMZ)});
	auto normReM = tl2::norm<t_vec_gl>(vecReM);
	auto normImM = tl2::norm<t_vec_gl>(vecImM);

	t_mat_gl matArrowRe = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
		vecReM,                                    // to
		1,                                         // post-scale
		tl2::create<t_vec_gl>({0, 0, 0}),          // post-translate
		tl2::create<t_vec_gl>({0, 0, 1}),          // from
		M * scale,                                 // pre-scale
		tl2::create<t_vec_gl>({posx, posy, posz})  // pre-translate
	);

	t_mat_gl matArrowIm = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
		vecImM,                                    // to
		1,                                         // post-scale
		tl2::create<t_vec_gl>({0, 0, 0}),          // post-translate
		tl2::create<t_vec_gl>({0, 0, 1}),          // from
		M * scale,                                 // pre-scale
		tl2::create<t_vec_gl>({posx, posy, posz})  // pre-translate
	);

	m_plot->GetRenderer()->SetObjectMatrix(objSphere, matSphere);
	m_plot->GetRenderer()->SetObjectMatrix(objArrowRe, matArrowRe);
	m_plot->GetRenderer()->SetObjectMatrix(objArrowIm, matArrowIm);
	m_plot->GetRenderer()->SetObjectLabel(objSphere, itemName->text().toStdString());
	m_plot->GetRenderer()->SetObjectCol(objSphere, r, g, b, 1);
	m_plot->GetRenderer()->SetObjectCol(objArrowRe, r, g, b, 1);
	m_plot->GetRenderer()->SetObjectCol(objArrowIm, 1.-r, 1.-g, 1.-b, 1);
	m_plot->GetRenderer()->SetObjectVisible(objArrowRe,
		!tl2::equals<t_real_gl>(normReM, 0, g_eps));
	m_plot->GetRenderer()->SetObjectVisible(objArrowIm,
		!tl2::equals<t_real_gl>(normImM, 0, g_eps));
	m_plot->update();
}


// ----------------------------------------------------------------------------
/**
 * mouse hovers over 3d object in unit cell view
 */
void MagStructFactDlg::PickerIntersection(const t_vec3_gl* pos,
	std::size_t objIdx, std::size_t /*triagIdx*/,
	const t_vec3_gl* /*posSphere*/)
{
	if(pos && m_plot)
		m_curPickedObj = long(objIdx);
	else
		m_curPickedObj = -1;

	if(m_curPickedObj > 0)
	{
		// find corresponding nucleus in table
		for(int row = 0; row < m_nuclei->rowCount(); ++row)
		{
			std::size_t objSphere = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole + 0).toUInt();
			std::size_t objArrowRe = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole + 1).toUInt();
			std::size_t objArrowIm = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole + 2).toUInt();

			if(long(objSphere)==m_curPickedObj || long(objArrowRe)==m_curPickedObj
				|| long(objArrowIm)==m_curPickedObj)
			{
				auto *itemname = m_nuclei->item(row, COL_NAME);
				auto *itemX = m_nuclei->item(row, COL_X);
				auto *itemY = m_nuclei->item(row, COL_Y);
				auto *itemZ = m_nuclei->item(row, COL_Z);

				t_vec r = tl2::create<t_vec>({ 0, 0, 0 });
				r[0] = tl2::stoval<t_real>(itemX->text().toStdString());
				r[1] = tl2::stoval<t_real>(itemY->text().toStdString());
				r[2] = tl2::stoval<t_real>(itemZ->text().toStdString());
				t_vec rlab = m_crystA * r;

				std::ostringstream ostr; ostr.precision(g_prec);
				ostr << itemname->text().toStdString();
				ostr << "; r = (" << r[0] << ", " << r[1] << ", " << r[2] << ") rlu";
				ostr << "; r = (" << rlab[0] << ", " << rlab[1] << ", " << rlab[2] << ") A";

				m_status3D->setText(ostr.str().c_str());
				break;
			}
		}
	}
	else
		m_status3D->setText("");
}


/**
 * mouse hovers over 3d object in super cell view
 */
void MagStructFactDlg::PickerIntersectionSC(const t_vec3_gl* pos,
	std::size_t objIdx, std::size_t /*triagIdx*/,
	const t_vec3_gl* /*posSphere*/)
{
	if(pos && m_plotSC)
	{
		const std::string& str = m_plotSC->GetRenderer()->GetObjectDataString(objIdx);
		m_status3DSC->setText(str.c_str());
	}
	else
		m_status3DSC->setText("");
}


/**
 * mouse button pressed
 */
void MagStructFactDlg::PlotMouseDown(bool left, bool /*mid*/, bool /*right*/)
{
	if(left && m_curPickedObj > 0)
	{
		// find corresponding nucleus in table
		for(int row = 0; row < m_nuclei->rowCount(); ++row)
		{
			std::size_t objSphere = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole + 0).toUInt();
			std::size_t objArrowRe = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole + 1).toUInt();
			std::size_t objArrowIm = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole + 2).toUInt();

			if(long(objSphere)==m_curPickedObj || long(objArrowRe)==m_curPickedObj
				|| long(objArrowIm)==m_curPickedObj)
			{
				m_nuclei->setCurrentCell(row, 0);
				break;
			}
		}
	}
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
/**
 * initialise objects for unit cell view
 */
void MagStructFactDlg::AfterGLInitialisation()
{
	if(!m_plot)
		return;
	SetGLInfos();

	// reference sphere and arrow for linked objects
	m_sphere = m_plot->GetRenderer()->AddSphere(0.05, 0.,0.,0., 1.,1.,1.,1.);
	m_arrow = m_plot->GetRenderer()->AddArrow(0.015, 0.25, 0.,0.,0.5,  1.,1.,1.,1.);
	m_plot->GetRenderer()->SetObjectVisible(m_sphere, false);
	m_plot->GetRenderer()->SetObjectVisible(m_arrow, false);

	// B matrix
	m_plot->GetRenderer()->SetBTrafo(m_crystB);

	// add all 3d objects
	Add3DItem(-1);
}


/**
 * initialise objects for super cell view
 */
void MagStructFactDlg::AfterGLInitialisationSC()
{
	if(!m_plotSC)
		return;
	SetGLInfos();

	// reference sphere and arrow for linked objects
	m_sphereSC = m_plotSC->GetRenderer()->AddSphere(0.05, 0.,0.,0., 1.,1.,1.,1.);
	m_arrowSC = m_plotSC->GetRenderer()->AddArrow(0.015, 0.25, 0.,0.,0.5,  1.,1.,1.,1.);
	m_plotSC->GetRenderer()->SetObjectVisible(m_sphereSC, false);
	m_plotSC->GetRenderer()->SetObjectVisible(m_arrowSC, false);

	// B matrix
	m_plotSC->GetRenderer()->SetBTrafo(m_crystB);

	// add all 3d objects (generated in calc)
	Calc();
}


/**
 * set descriptions of the gl device in the info tab
 */
void MagStructFactDlg::SetGLInfos()
{
	static bool already_set = false;
	if(already_set)
		return;

	// try whichever gl plotter is available first
	for(auto* plot : { m_plot.get(), m_plotSC.get() })
	{
		if(!plot) continue;

		auto [strGlVer, strGlShaderVer, strGlVendor, strGlRenderer] =
			plot->GetRenderer()->GetGlDescr();
		m_labelGlInfos[0]->setText(QString("GL Version: ") +
			strGlVer.c_str() + QString("."));
		m_labelGlInfos[1]->setText(QString("GL Shader Version: ") +
			strGlShaderVer.c_str() + QString("."));
		m_labelGlInfos[2]->setText(QString("GL Vendor: ") +
			strGlVendor.c_str() + QString("."));
		m_labelGlInfos[3]->setText(QString("GL Device: ") +
			strGlRenderer.c_str() + QString("."));

		already_set = true;
		break;
	}
}
// ----------------------------------------------------------------------------
