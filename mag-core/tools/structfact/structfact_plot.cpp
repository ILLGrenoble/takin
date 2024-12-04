/**
 * structure factor tool
 * @author Tobias Weber <tweber@ill.fr>
 * @date Dec-2018
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

#include <QtWidgets/QGridLayout>

#include "structfact.h"

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"

using namespace tl2_ops;


void StructFactDlg::ShowStructPlot()
{
	// plot widget
	if(!m_dlgPlot)
	{
		m_dlgPlot = new QDialog(this);
		m_dlgPlot->setWindowTitle("Unit Cell - 3D View");
		m_dlgPlot->setFont(this->font());

		m_plot = std::make_shared<tl2::GlPlot>(this);
		m_plot->GetRenderer()->SetRestrictCamTheta(false);
		m_plot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 5, 5, 5 }));
		m_plot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -5, -5, -5 }));
		m_plot->GetRenderer()->SetCoordMax(1.);
		m_plot->GetRenderer()->GetCamera().SetDist(1.5);
		m_plot->GetRenderer()->GetCamera().UpdateTransformation();

		auto labCoordSys = new QLabel("Coordinate System:", /*m_dlgPlot*/ this);
		auto comboCoordSys = new QComboBox(/*m_dlgPlot*/ this);
		m_status3D = new QLabel(/*m_dlgPlot*/ this);

		comboCoordSys->addItem("Fractional Units (rlu)");
		comboCoordSys->addItem("Lab Units (\xe2\x84\xab)");

		m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
		labCoordSys->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

		auto grid = new QGridLayout(m_dlgPlot);
		grid->setSpacing(2);
		grid->setContentsMargins(4,4,4,4);
		grid->addWidget(m_plot.get(), 0,0,1,2);
		grid->addWidget(labCoordSys, 1,0,1,1);
		grid->addWidget(comboCoordSys, 1,1,1,1);
		grid->addWidget(m_status3D, 2,0,1,2);

		connect(m_plot.get(), &tl2::GlPlot::AfterGLInitialisation,
			this, &StructFactDlg::AfterGLInitialisation);
		connect(m_plot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection,
			this, &StructFactDlg::PickerIntersection);
		connect(m_plot.get(), &tl2::GlPlot::MouseDown,
			this, &StructFactDlg::PlotMouseDown);
		connect(m_plot.get(), &tl2::GlPlot::MouseUp,
			this, &StructFactDlg::PlotMouseUp);
		connect(comboCoordSys, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
			[this](int val)
		{
			if(this->m_plot)
				this->m_plot->GetRenderer()->SetCoordSys(val);
		});


		if(m_sett && m_sett->contains("geo_3dview"))
			m_dlgPlot->restoreGeometry(m_sett->value("geo_3dview").toByteArray());
		else
			m_dlgPlot->resize(500,500);
	}

	m_dlgPlot->show();
	m_dlgPlot->raise();
	m_dlgPlot->focusWidget();
}


/**
 * add 3d object
 */
void StructFactDlg::Add3DItem(int row)
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

	auto *itemName = m_nuclei->item(row, COL_NAME);
	auto *itemx = m_nuclei->item(row, COL_X);
	auto *itemy = m_nuclei->item(row, COL_Y);
	auto *itemz = m_nuclei->item(row, COL_Z);
	auto *itemsc = m_nuclei->item(row, COL_RAD);
	auto *itemCol = m_nuclei->item(row, COL_COL);

	t_real_gl posx = tl2::stoval<t_real_gl>(itemx->text().toStdString());
	t_real_gl posy = tl2::stoval<t_real_gl>(itemy->text().toStdString());
	t_real_gl posz = tl2::stoval<t_real_gl>(itemz->text().toStdString());
	t_real_gl scale = tl2::stoval<t_real_gl>(itemsc->text().toStdString());

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	qreal r=1, g=1, b=1;
#else
	float r=1, g=1, b=1;
#endif
	QColor col{itemCol->text()};
	col.getRgbF(&r, &g, &b);

	auto obj = m_plot->GetRenderer()->AddLinkedObject(m_sphere, 0,0,0, r,g,b,1);
	m_plot->GetRenderer()->SetObjectMatrix(obj,
		tl2::hom_translation<t_mat_gl>(posx, posy, posz) *
		tl2::hom_scaling<t_mat_gl>(scale, scale, scale));
	m_plot->GetRenderer()->SetObjectLabel(obj, itemName->text().toStdString());
	m_plot->update();

	m_nuclei->item(row, COL_NAME)->setData(Qt::UserRole, unsigned(obj));
}


/**
 * mouse hovers over 3d object
 */
void StructFactDlg::PickerIntersection(
	const t_vec3_gl* pos, std::size_t objIdx,
	[[maybe_unused]] std::size_t triagIdx,
	[[maybe_unused]] const t_vec3_gl* posSphere)
{
	if(pos)
		m_curPickedObj = long(objIdx);
	else
		m_curPickedObj = -1;


	if(m_curPickedObj > 0)
	{
		// find corresponding nucleus in table
		for(int row = 0; row < m_nuclei->rowCount(); ++row)
		{
			if(std::size_t obj = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole).toUInt(); long(obj) == m_curPickedObj)
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

				Set3DStatusMsg(ostr.str().c_str());
				break;
			}
		}
	}
	else
		Set3DStatusMsg("");
}



/**
 * set status label text in 3d dialog
 */
void StructFactDlg::Set3DStatusMsg(const std::string& msg)
{
	m_status3D->setText(msg.c_str());
}



/**
 * mouse button pressed
 */
void StructFactDlg::PlotMouseDown(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
	if(left && m_curPickedObj > 0)
	{
		// find corresponding nucleus in table
		for(int row=0; row<m_nuclei->rowCount(); ++row)
		{
			if(std::size_t obj = m_nuclei->item(row, COL_NAME)->
				data(Qt::UserRole).toUInt(); long(obj)==m_curPickedObj)
			{
				m_nuclei->setCurrentCell(row, 0);
				break;
			}
		}
	}
}


/**
 * mouse button released
 */
void StructFactDlg::PlotMouseUp(
	[[maybe_unused]] bool left,
	[[maybe_unused]] bool mid,
	[[maybe_unused]] bool right)
{
}


void StructFactDlg::AfterGLInitialisation()
{
	if(!m_plot)
		return;

	// reference sphere for linked objects
	m_sphere = m_plot->GetRenderer()->AddSphere(0.05, 0.,0.,0., 1.,1.,1.,1.);
	m_plot->GetRenderer()->SetObjectVisible(m_sphere, false);

	// B matrix
	m_plot->GetRenderer()->SetBTrafo(m_crystB);

	// add all 3d objects
	Add3DItem(-1);

	// GL device info
	auto [strGlVer, strGlShaderVer, strGlVendor, strGlRenderer]
		= m_plot->GetRenderer()->GetGlDescr();
	m_labelGlInfos[0]->setText(QString("GL Version: ") +
		strGlVer.c_str() + QString("."));
	m_labelGlInfos[1]->setText(QString("GL Shader Version: ") +
		strGlShaderVer.c_str() + QString("."));
	m_labelGlInfos[2]->setText(QString("GL Vendor: ") +
		strGlVendor.c_str() + QString("."));
	m_labelGlInfos[3]->setText(QString("GL Device: ") +
		strGlRenderer.c_str() + QString("."));
}
