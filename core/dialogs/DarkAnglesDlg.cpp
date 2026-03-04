/**
 * dark angles dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2017, 28-jul-2022
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include "DarkAnglesDlg.h"
#include "tlibs/math/linalg.h"
#include "tlibs/string/string.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QComboBox>

using t_real = t_real_glob;


enum class AngleInfo : int
{
	START_ANGLE  = 0,
	STOP_ANGLE   = 1,
	OFFSET_ANGLE = 2,

	CENTRE       = 3,
	RELATIVE     = 4,
};


DarkAnglesDlg::DarkAnglesDlg(QWidget* pParent, QSettings *pSettings)
	: QDialog(pParent), m_pSettings(pSettings)
{
	setupUi(this);
	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 2);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") &&
			font.fromString(m_pSettings->value("main/font_gen", "").toString()))
		{
			setFont(font);
		}
	}

	btnAddAngle->setIcon(load_icon("res/icons/list-add.svg"));
	btnDelAngle->setIcon(load_icon("res/icons/list-remove.svg"));

	btnAdd->setIcon(load_icon("res/icons/list-add.svg"));
	btnDel->setIcon(load_icon("res/icons/list-remove.svg"));
	btnSave->setIcon(load_icon("res/icons/document-save.svg"));
	btnLoad->setIcon(load_icon("res/icons/document-open.svg"));

	QObject::connect(btnAddAngle, &QAbstractButton::clicked,
		this, &DarkAnglesDlg::AddAngle);
	QObject::connect(btnDelAngle, &QAbstractButton::clicked,
		this, &DarkAnglesDlg::RemoveAngle);
	QObject::connect(buttonBox, &QDialogButtonBox::clicked,
		this, &DarkAnglesDlg::ButtonBoxClicked);

	// list
	QObject::connect(btnAdd, &QAbstractButton::clicked,
		this, static_cast<void(DarkAnglesDlg::*)()>(&DarkAnglesDlg::AddAnglesToList));
	QObject::connect(btnDel, &QAbstractButton::clicked,
		this, &DarkAnglesDlg::RemAnglesFromList);
	QObject::connect(btnLoad, &QAbstractButton::clicked,
		this, &DarkAnglesDlg::LoadList);
	QObject::connect(btnSave, &QAbstractButton::clicked,
		this, &DarkAnglesDlg::SaveList);
	QObject::connect(listSeq, &QListWidget::itemSelectionChanged,
		this, &DarkAnglesDlg::ListItemSelected);
	QObject::connect(listSeq, &QListWidget::itemDoubleClicked,
		this, &DarkAnglesDlg::ListItemDoubleClicked);

	if(m_pSettings && m_pSettings->contains("darkangles/geo"))
		restoreGeometry(m_pSettings->value("darkangles/geo").toByteArray());
}


DarkAnglesDlg::~DarkAnglesDlg()
{
	ClearList();
}


/**
 * removes the currently selected items from the dark angles list
 */
void DarkAnglesDlg::RemoveAngle()
{
	const bool bSort = tableAngles->isSortingEnabled();
	tableAngles->setSortingEnabled(false);

	bool bNothingRemoved = true;

	// remove selected rows
	QList<QTableWidgetSelectionRange> lstSel = tableAngles->selectedRanges();
	for(QTableWidgetSelectionRange& range : lstSel)
	{
		for(int iRow=range.bottomRow(); iRow>=range.topRow(); --iRow)
		{
			tableAngles->removeRow(iRow);
			bNothingRemoved = false;
		}
	}

	// remove last row if nothing is selected
	if(bNothingRemoved)
		tableAngles->removeRow(tableAngles->rowCount()-1);

	tableAngles->setSortingEnabled(bSort);
}


/**
 * adds a single new item to the dark angles list
 */
void DarkAnglesDlg::AddAngle()
{
	const bool bSort = tableAngles->isSortingEnabled();
	tableAngles->setSortingEnabled(0);

	int iRow = tableAngles->rowCount();
	tableAngles->insertRow(iRow);

	tableAngles->setItem(iRow, static_cast<int>(AngleInfo::START_ANGLE),
		new QTableWidgetItem("0"));
	tableAngles->setItem(iRow, static_cast<int>(AngleInfo::STOP_ANGLE),
		new QTableWidgetItem("0"));
	tableAngles->setItem(iRow, static_cast<int>(AngleInfo::OFFSET_ANGLE),
		new QTableWidgetItem("0"));

	QComboBox *pComboCentreOn = new QComboBox(tableAngles);
	pComboCentreOn->addItem("Monochromator");
	pComboCentreOn->addItem("Sample");
	pComboCentreOn->addItem("Analyser");
	pComboCentreOn->setCurrentIndex(1);
	tableAngles->setCellWidget(iRow, static_cast<int>(AngleInfo::CENTRE), pComboCentreOn);

	QComboBox *pComboRelativeTo = new QComboBox(tableAngles);
	pComboRelativeTo->addItem("Xtal Angle");
	pComboRelativeTo->addItem("In Axis");
	pComboRelativeTo->addItem("Out Axis");
	tableAngles->setCellWidget(iRow, static_cast<int>(AngleInfo::RELATIVE), pComboRelativeTo);

	tableAngles->setSortingEnabled(bSort);
}


/**
 * loads the items in the dark angles list from 'vecAngles'
 */
void DarkAnglesDlg::SetDarkAngles(const std::vector<DarkAngle<t_real>>& vecAngles)
{
	const bool bSort = tableAngles->isSortingEnabled();
	tableAngles->setSortingEnabled(0);

	// add missing rows
	while(tableAngles->rowCount() < int(vecAngles.size()))
		AddAngle();

	// remove superfluous rows
	while(tableAngles->rowCount() > int(vecAngles.size()))
		RemoveAngle();


	for(std::size_t iRow=0; iRow<vecAngles.size(); ++iRow)
	{
		const DarkAngle<t_real>& angle = vecAngles[iRow];

		tableAngles->item(iRow, static_cast<int>(AngleInfo::START_ANGLE))->
			setText(tl::var_to_str(angle.dAngleStart).c_str());
		tableAngles->item(iRow, static_cast<int>(AngleInfo::STOP_ANGLE))->
			setText(tl::var_to_str(angle.dAngleEnd).c_str());
		tableAngles->item(iRow, static_cast<int>(AngleInfo::OFFSET_ANGLE))->
			setText(tl::var_to_str(angle.dAngleOffs).c_str());

		QComboBox* pComboCentreOn = (QComboBox*)tableAngles->
			cellWidget(iRow, static_cast<int>(AngleInfo::CENTRE));
		QComboBox* pComboRelativeTo = (QComboBox*)tableAngles->
			cellWidget(iRow, static_cast<int>(AngleInfo::RELATIVE));

		pComboCentreOn->setCurrentIndex(angle.iCentreOn);
		pComboRelativeTo->setCurrentIndex(angle.iRelativeTo);
	}

	tableAngles->setSortingEnabled(bSort);
}


/**
 * gets the dark angles from the list
 */
std::vector<DarkAngle<t_real>> DarkAnglesDlg::DarkAnglesDlg::GetDarkAngles() const
{
	std::vector<DarkAngle<t_real>> vecAngles;
	vecAngles.reserve(tableAngles->rowCount());

	for(int iRow=0; iRow<tableAngles->rowCount(); ++iRow)
	{
		DarkAngle<t_real> angle;
		angle.dAngleStart =
			tl::str_to_var_parse<t_real>(tableAngles->item(
				iRow, static_cast<int>(AngleInfo::START_ANGLE))->text().toStdString());
		angle.dAngleEnd =
			tl::str_to_var_parse<t_real>(tableAngles->item(
				iRow, static_cast<int>(AngleInfo::STOP_ANGLE))->text().toStdString());
		angle.dAngleOffs =
			tl::str_to_var_parse<t_real>(tableAngles->item(
				iRow, static_cast<int>(AngleInfo::OFFSET_ANGLE))->text().toStdString());

		QComboBox* pComboCentreOn = (QComboBox*)tableAngles->
			cellWidget(iRow, static_cast<int>(AngleInfo::CENTRE));
		QComboBox* pComboRelativeTo = (QComboBox*)tableAngles->
			cellWidget(iRow, static_cast<int>(AngleInfo::RELATIVE));

		angle.iCentreOn = pComboCentreOn->currentIndex();
		angle.iRelativeTo = pComboRelativeTo->currentIndex();

		vecAngles.emplace_back(std::move(angle));
	}

	return vecAngles;
}


/**
 * emits the current list of dark angles
 */
void DarkAnglesDlg::SendApplyDarkAngles()
{
	std::vector<DarkAngle<t_real>> angles = GetDarkAngles();
	emit ApplyDarkAngles(angles);
}


void DarkAnglesDlg::AddAnglesToList(const std::vector<DarkAngle<t_real>>& _angles)
{
	std::wostringstream ostrCaption;
	ostrCaption.precision(g_iPrecGfx);

	for(std::size_t angle_idx = 0; angle_idx < _angles.size(); ++angle_idx)
	{
		t_real start = _angles[angle_idx].dAngleStart;
		t_real end = _angles[angle_idx].dAngleEnd;
		t_real offs = _angles[angle_idx].dAngleOffs;

		tl::set_eps_0(start, g_dEps);
		tl::set_eps_0(end, g_dEps);
		tl::set_eps_0(offs, g_dEps);

		ostrCaption << "[" << start << ", " << end << "] + " << offs;
		if(angle_idx != _angles.size()-1)
			ostrCaption << "\n";
	}

	QString qstr = QString::fromWCharArray(ostrCaption.str().c_str());
	QListWidgetItem* item = new QListWidgetItem(qstr, listSeq);
	std::vector<DarkAngle<t_real>>* angles = new std::vector<DarkAngle<t_real>>{_angles};
	item->setData(Qt::UserRole, QVariant::fromValue<void*>(angles));
}


void DarkAnglesDlg::AddAnglesToList()
{
	std::vector<DarkAngle<t_real>> angles = GetDarkAngles();

	if(angles.size() == 0)
	{
		QMessageBox::critical(this, "Error", "No dark angles defined. Please define at least one using \"Add Angle\".");
		return;
	}

	AddAnglesToList(angles);
}


void DarkAnglesDlg::RemAnglesFromList()
{
	QListWidgetItem *item = listSeq->currentItem();
	if(item)
	{
		std::vector<DarkAngle<t_real>>* angles =
			(std::vector<DarkAngle<t_real>>*)item->data(
				Qt::UserRole).value<void*>();
		if(angles) delete angles;
		delete item;
	}
}


void DarkAnglesDlg::ClearList()
{
	while(listSeq->count())
	{
		QListWidgetItem *item = listSeq->item(0);
		if(!item) break;

		std::vector<DarkAngle<t_real>>* angles =
			(std::vector<DarkAngle<t_real>>*)item->data(
				Qt::UserRole).value<void*>();
		if(angles) delete angles;
		delete item;
        }
}


void DarkAnglesDlg::LoadList()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = "~";
	if(m_pSettings)
		strDirLast = m_pSettings->value("darkangles/last_dir", "~").toString();
	QString qstrFile = QFileDialog::getOpenFileName(this,
		"Load Positions", strDirLast,
		"TAZ files (*.taz *.TAZ)", nullptr,
		fileopt);
	if(qstrFile == "")
		return;


	std::string strFile = qstrFile.toStdString();
	std::string strDir = tl::get_dir(strFile);

	tl::Prop<std::string> xml;
	if(!xml.Load(strFile.c_str(), tl::PropType::XML))
	{
		QMessageBox::critical(this, "Error", "Could not load positions.");
		return;
	}

	Load(xml, strXmlRoot);
	if(m_pSettings)
		m_pSettings->setValue("darkangles/last_dir", QString(strDir.c_str()));
}


void DarkAnglesDlg::SaveList()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = "~";
	if(m_pSettings)
		m_pSettings->value("darkangles/last_dir", "~").toString();
	QString qstrFile = QFileDialog::getSaveFileName(this,
		"Save Dark Angles", strDirLast,
		"TAZ files (*.taz *.TAZ)", nullptr,
		fileopt);

	if(qstrFile == "")
		return;

	std::string strFile = qstrFile.toStdString();
	std::string strDir = tl::get_dir(strFile);
	if(tl::get_fileext(strFile,1) != "taz")
		strFile += ".taz";

	std::map<std::string, std::string> mapConf;
	Save(mapConf, strXmlRoot);

	tl::Prop<std::string> xml;
	xml.Add(mapConf);
	if(!xml.Save(strFile.c_str(), tl::PropType::XML))
	{
		QMessageBox::critical(this, "Error", "Could not save dark angles.");
		return;
	}

	if(m_pSettings)
		m_pSettings->setValue("darkangles/last_dir", QString(strDir.c_str()));
}


void DarkAnglesDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	// save current configuration
	std::vector<DarkAngle<t_real>> angles = GetDarkAngles();

	mapConf[strXmlRoot + "darkangles/num"] = tl::var_to_str(angles.size());
	for(std::size_t angle_idx=0; angle_idx<angles.size(); ++angle_idx)
	{
		const DarkAngle<t_real>& angle = angles[angle_idx];

		std::string strCfgNr = tl::var_to_str(angle_idx);
		mapConf[strXmlRoot + "darkangles/" + strCfgNr + "/start"] =
			tl::var_to_str(angle.dAngleStart);
		mapConf[strXmlRoot + "darkangles/" + strCfgNr + "/end"] =
			tl::var_to_str(angle.dAngleEnd);
		mapConf[strXmlRoot + "darkangles/" + strCfgNr + "/offs"] =
			tl::var_to_str(angle.dAngleOffs);
		mapConf[strXmlRoot + "darkangles/" + strCfgNr + "/centreon"] =
			tl::var_to_str(angle.iCentreOn);
		mapConf[strXmlRoot + "darkangles/" + strCfgNr + "/relativeto"] =
			tl::var_to_str(angle.iRelativeTo);
        }


	// save stored configurations
	int num_stored = listSeq->count();
	mapConf[strXmlRoot + "darkangles/num_stored"] = tl::var_to_str(num_stored);
	for(int stored=0; stored<num_stored; ++stored)
	{
		std::vector<DarkAngle<t_real>>* stored_angles =
			(std::vector<DarkAngle<t_real>>*)listSeq->item(stored)->data(
				Qt::UserRole).value<void*>();

		std::string strStoredNr = "stored_" + tl::var_to_str(stored) + "/";

		for(std::size_t angle_idx=0; angle_idx<stored_angles->size(); ++angle_idx)
		{
			const DarkAngle<t_real>& angle = (*stored_angles)[angle_idx];

			std::string strCfgNr = tl::var_to_str(angle_idx);
			mapConf[strXmlRoot + "darkangles/" + strStoredNr + strCfgNr + "/start"] =
				tl::var_to_str(angle.dAngleStart);
			mapConf[strXmlRoot + "darkangles/" + strStoredNr + strCfgNr + "/end"] =
				tl::var_to_str(angle.dAngleEnd);
			mapConf[strXmlRoot + "darkangles/" + strStoredNr + strCfgNr + "/offs"] =
				tl::var_to_str(angle.dAngleOffs);
			mapConf[strXmlRoot + "darkangles/" + strStoredNr + strCfgNr + "/centreon"] =
				tl::var_to_str(angle.iCentreOn);
			mapConf[strXmlRoot + "darkangles/" + strStoredNr + strCfgNr + "/relativeto"] =
				tl::var_to_str(angle.iRelativeTo);
		}
	}
}


void DarkAnglesDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	// load current configuration
	bool ok;
	unsigned int num_angles = xml.Query<unsigned int>(strXmlRoot + "darkangles/num", 0, &ok);
	if(ok)
	{
		std::vector<DarkAngle<t_real>> angles;
		angles.reserve(num_angles);

		for(unsigned int angle_idx=0; angle_idx<num_angles; ++angle_idx)
		{
			DarkAngle<t_real> angle;

			std::string strNr = tl::var_to_str(angle_idx);
			angle.dAngleStart = xml.Query<t_real>(strXmlRoot + "darkangles/" + strNr + "/start", 0.);
			angle.dAngleEnd = xml.Query<t_real>(strXmlRoot + "darkangles/" + strNr + "/end", 0.);
			angle.dAngleOffs = xml.Query<t_real>(strXmlRoot + "darkangles/" + strNr + "/offs", 0.);
			angle.iCentreOn = xml.Query<int>(strXmlRoot + "darkangles/" + strNr + "/centreon", 1);
			angle.iRelativeTo = xml.Query<int>(strXmlRoot + "darkangles/" + strNr + "/relativeto", 0);

			angles.emplace_back(std::move(angle));
		}

		SetDarkAngles(angles);
	}


	// load stored configurations
	ClearList();
	unsigned int num_stored = xml.Query<unsigned int>(strXmlRoot + "darkangles/num_stored", 0, &ok);
	if(ok)
	{
		for(unsigned int stored_idx=0; stored_idx<num_stored; ++stored_idx)
		{
			std::string strStoredNr = "stored_" + tl::var_to_str(stored_idx) + "/";

			std::vector<DarkAngle<t_real>> angles;
			angles.reserve(num_angles);

			for(unsigned int angle_idx=0; angle_idx<num_angles; ++angle_idx)
			{
				DarkAngle<t_real> angle;

				std::string strNr = tl::var_to_str(angle_idx);
				angle.dAngleStart = xml.Query<t_real>(strXmlRoot + "darkangles/" +  strStoredNr + strNr + "/start", 0.);
				angle.dAngleEnd = xml.Query<t_real>(strXmlRoot + "darkangles/" + strStoredNr + strNr + "/end", 0.);
				angle.dAngleOffs = xml.Query<t_real>(strXmlRoot + "darkangles/" + strStoredNr + strNr + "/offs", 0.);
				angle.iCentreOn = xml.Query<int>(strXmlRoot + "darkangles/" + strStoredNr + strNr + "/centreon", 1);
				angle.iRelativeTo = xml.Query<int>(strXmlRoot + "darkangles/" + strStoredNr + strNr + "/relativeto", 0);

				angles.emplace_back(std::move(angle));
			}

			AddAnglesToList(angles);
		}
	}
}


void DarkAnglesDlg::SetAnglesFromList(QListWidgetItem* item)
{
	if(!item) return;
	std::vector<DarkAngle<t_real>>* angles =
		(std::vector<DarkAngle<t_real>>*)item->data(
			Qt::UserRole).value<void*>();
	if(angles)
		SetDarkAngles(*angles);
}


void DarkAnglesDlg::ListItemSelected()
{
	QListWidgetItem *item = listSeq->currentItem();
	SetAnglesFromList(item);
}


void DarkAnglesDlg::ListItemDoubleClicked(QListWidgetItem* item)
{
	SetAnglesFromList(item);
	SendApplyDarkAngles();
}


void DarkAnglesDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
		buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		SendApplyDarkAngles();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("darkangles/geo", saveGeometry());

		QDialog::accept();
	}
}

void DarkAnglesDlg::closeEvent(QCloseEvent* pEvt)
{
	QDialog::closeEvent(pEvt);
}


#include "moc_DarkAnglesDlg.cpp"
