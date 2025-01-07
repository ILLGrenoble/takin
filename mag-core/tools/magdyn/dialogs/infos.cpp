/**
 * magnon dynamics -- info dialog
 * @author Tobias Weber <tweber@ill.fr>
 * @date june-2024
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

#include "infos.h"
#include "defs.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QLabel>

#include <boost/version.hpp>
#include <boost/config.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;



/**
 * set up the gui
 */
InfoDlg::InfoDlg(QWidget* parent, QSettings *sett)
	: QDialog{parent}, m_sett(sett)
{
	setWindowTitle("About");
	setSizeGripEnabled(true);

	auto infopanel = new QWidget(this);
	auto grid = new QGridLayout(infopanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	std::string strBoost = BOOST_LIB_VERSION;
	algo::replace_all(strBoost, "_", ".");

	auto labelTitle = new QLabel("Takin / Magnetic Dynamics Calculator", infopanel);
	auto fontTitle = labelTitle->font();
	fontTitle.setBold(true);
	labelTitle->setFont(fontTitle);
	labelTitle->setAlignment(Qt::AlignHCenter);

	auto labelVersion = new QLabel("Version " MAGCORE_VER ".", infopanel);
	labelVersion->setAlignment(Qt::AlignHCenter);

	auto labelAuthor = new QLabel("Written by Tobias Weber <tweber@ill.fr>.", infopanel);
	labelAuthor->setAlignment(Qt::AlignHCenter);

	auto labelDate = new QLabel("2022 - 2024.", infopanel);
	labelDate->setAlignment(Qt::AlignHCenter);

	auto labelPaper = new QLabel(
		"Paper DOI: "
		"<a href=\"https://doi.org/10.1016/j.softx.2023.101471\">10.1016/j.softx.2023.101471</a>."
		"<br>This program implements the formalism from "
		"<a href=\"https://doi.org/10.1088/0953-8984/27/16/166002\">this paper</a> "
		"(which is also available <a href=\"https://doi.org/10.48550/arXiv.1402.6069\">here</a>).",
		infopanel);
	labelPaper->setWordWrap(true);
	labelPaper->setOpenExternalLinks(true);

	auto labelLicense = new QLabel(
		"<p>This program is free software: you can redistribute it and/or modify "
		"it under the terms of the <u>GNU General Public License</u> as published by "
		"the Free Software Foundation, <u>version 3</u> of the License.</p>"

		"<p>This program is distributed in the hope that it will be useful, "
		"but WITHOUT ANY WARRANTY; without even the implied warranty of "
		"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "
		"See the GNU General Public License for more details.</p>"

		"<p>You should have received a copy of the GNU General Public License "
		"along with this program. If not, see "
		"<a href=\"http://www.gnu.org/licenses/\">&lt;http://www.gnu.org/licenses/&gt;</a>.</p>",
		infopanel);
		labelLicense->setWordWrap(true);
		labelLicense->setOpenExternalLinks(true);

	// renderer infos
	for(int i = 0; i < 4; ++i)
	{
		m_labelGlInfos[i] = new QLabel("", infopanel);
		m_labelGlInfos[i]->setSizePolicy(
			QSizePolicy::Ignored,
			m_labelGlInfos[i]->sizePolicy().verticalPolicy());
	}

	int y = 0;
	grid->addWidget(labelTitle, y++,0, 1,1);
	grid->addWidget(labelVersion, y++,0, 1,1);
	grid->addWidget(labelAuthor, y++,0, 1,1);
	grid->addWidget(labelDate, y++,0, 1,1);

	auto sep1 = new QFrame(infopanel);
	sep1->setFrameStyle(QFrame::HLine);
	auto sep2 = new QFrame(infopanel);
	sep2->setFrameStyle(QFrame::HLine);
	auto sep3 = new QFrame(infopanel);
	sep3->setFrameStyle(QFrame::HLine);
	auto sep4 = new QFrame(infopanel);
	sep4->setFrameStyle(QFrame::HLine);
	auto sep5 = new QFrame(infopanel);
	sep5->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep1, y++,0, 1,1);

	grid->addWidget(new QLabel(
		QString("Compiler: ") +
		QString(BOOST_COMPILER) + ".",
		infopanel), y++,0, 1,1);
	grid->addWidget(new QLabel(
		QString("C++ Library: ") +
		QString(BOOST_STDLIB) + ".",
		infopanel), y++,0, 1,1);
	grid->addWidget(new QLabel(
		QString("Build Date: ") +
		QString(__DATE__) + ", " +
		QString(__TIME__) + ".",
		infopanel), y++,0, 1,1);
	grid->addWidget(new QLabel(
		QString("Using %1-bit real and %2-bit integer type.").arg(sizeof(t_real)*8).arg(sizeof(t_size)*8),
		infopanel), y++,0, 1,1);

	grid->addWidget(sep2, y++,0, 1,1);

	grid->addWidget(new QLabel(
		QString("Qt Version: ") +
		QString(QT_VERSION_STR) + ".",
		infopanel), y++,0, 1,1);
	grid->addWidget(new QLabel(
		QString("Boost Version: ") +
		strBoost.c_str() + ".",
		infopanel), y++,0, 1,1);

	grid->addWidget(sep3, y++,0, 1,1);
	grid->addWidget(labelPaper, y++,0, 1,1);
	grid->addWidget(sep4, y++,0, 1,1);
	grid->addWidget(labelLicense, y++,0, 1,1);
	grid->addWidget(sep5, y++,0, 1,1);

	for(int i = 0; i < 4; ++i)
		grid->addWidget(m_labelGlInfos[i], y++,0, 1,1);

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++,0, 1,1);

	QDialogButtonBox *btnbox = new QDialogButtonBox(this);
	btnbox->addButton(QDialogButtonBox::Ok);
	btnbox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	connect(btnbox, &QDialogButtonBox::accepted, this, &InfoDlg::accept);

	auto dlgGrid = new QGridLayout(this);
	dlgGrid->setSpacing(4);
	dlgGrid->setContentsMargins(8, 8, 8, 8);
	dlgGrid->addWidget(infopanel, 0,0, 1,4);
	dlgGrid->addWidget(btnbox, 1,3, 1,1);

	// restore settings
	if(m_sett)
	{
		// restore dialog geometry
		if(m_sett->contains("infos/geo"))
			restoreGeometry(m_sett->value("infos/geo").toByteArray());
		else
			resize(700, 700);
        }
}



void InfoDlg::SetGlDeviceInfos(const std::string& ver, const std::string& shader_ver,
	const std::string& vendor, const std::string& renderer)
{
	m_labelGlInfos[0]->setText(QString("GL Version: %1.").arg(ver.c_str()));
	m_labelGlInfos[1]->setText(QString("GL Shader Version: %1.").arg(shader_ver.c_str()));
	m_labelGlInfos[2]->setText(QString("GL Vendor: %1.").arg(vendor.c_str()));
	m_labelGlInfos[3]->setText(QString("GL Device: %1.").arg(renderer.c_str()));
}



/**
 * close the dialog
 */
void InfoDlg::accept()
{
	if(m_sett)
	{
		// save dialog geometry
		m_sett->setValue("infos/geo", saveGeometry());
	}

	QDialog::accept();
}
