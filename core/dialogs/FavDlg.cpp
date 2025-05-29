/**
 * Favourite Positions Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "FavDlg.h"
#include "tlibs/math/linalg.h"
#include "tlibs/string/spec_char.h"
#include <sstream>

using t_real = t_real_glob;


FavDlg::FavDlg(QWidget* pParent, QSettings* pSett) : QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}

	btnAdd->setIcon(load_icon("res/icons/list-add.svg"));
	btnDel->setIcon(load_icon("res/icons/list-remove.svg"));

	connect(buttonBox, &QDialogButtonBox::clicked, this, &FavDlg::ButtonBoxClicked);
	connect(btnAdd, &QToolButton::clicked, this, static_cast<void (FavDlg::*)()>(&FavDlg::AddPosToList));
	connect(btnDel, &QToolButton::clicked, this, static_cast<void (FavDlg::*)()>(&FavDlg::RemPosFromList));
	connect(listSeq, &QListWidget::itemSelectionChanged, this, &FavDlg::ListItemSelected);
	connect(listSeq, &QListWidget::itemDoubleClicked, this, &FavDlg::ListItemDoubleClicked);

	if(m_pSettings)
	{
		// load saved positions
		QString strLst = m_pSettings->value("fav_pos/pos_list").toString();
		tl::Prop<std::string> prop;

		std::istringstream istr(strLst.toStdString());
		prop.Load(istr, tl::PropType::XML);
		Load(prop, "");

		if(m_pSettings->contains("fav_pos/geo"))
			restoreGeometry(m_pSettings->value("fav_pos/geo").toByteArray());
	}
}

FavDlg::~FavDlg()
{
	ClearList();
}


void FavDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
	   buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		const FavHklPos* pPos = nullptr;
		QListWidgetItem* pItem = listSeq->currentItem();
		if(pItem)
			pPos = (FavHklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(pPos)
			ApplyPos(pPos);
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
		{
			// save positions
			std::map<std::string, std::string> mapConf;
			Save(mapConf, "");
			tl::Prop<std::string> prop;
			prop.Add(mapConf);
			std::ostringstream ostr;
			prop.Save(ostr, tl::PropType::XML);

			m_pSettings->setValue("fav_pos/pos_list", QString(ostr.str().c_str()));

			m_pSettings->setValue("fav_pos/geo", saveGeometry());
		}
		QDialog::accept();
	}
}


void FavDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


//------------------------------------------------------------------------------

void FavDlg::ListItemSelected()
{}


void FavDlg::ApplyPos(const FavHklPos* pPos)
{
	if(!pPos)
		return;
	emit ChangePos(*pPos);
}


void FavDlg::ListItemDoubleClicked(QListWidgetItem* pItem)
{
	if(!pItem)
		return;
	const FavHklPos* pPos = (FavHklPos*)pItem->data(Qt::UserRole).value<void*>();
	ApplyPos(pPos);
}


void FavDlg::AddPosToList(const FavHklPos& pos)
{
	FavHklPos *pPos = new FavHklPos(pos);

	for(t_real* d : {&pPos->dhstart, &pPos->dhstop,
		&pPos->dkstart, &pPos->dkstop,
		&pPos->dlstart, &pPos->dlstop,
		&pPos->dEstart, &pPos->dEstop})
		tl::set_eps_0(*d, g_dEps);

	std::wostringstream ostrHKL;
	ostrHKL.precision(g_iPrecGfx);
	ostrHKL << "Q = (" << pPos->dhstart << ", " << pPos->dkstart << ", " << pPos->dlstart << ") rlu";
	ostrHKL << ", E = " << pPos->dEstart << " meV ->\n"; //<< tl::get_spec_char_utf16("rightarrow") << "\n";
	ostrHKL << "Q = (" << pPos->dhstop << ", " << pPos->dkstop << ", " << pPos->dlstop << ") rlu";
	ostrHKL << ", E = " << pPos->dEstop << " meV";

	QString qstr = QString::fromWCharArray(ostrHKL.str().c_str());
	QListWidgetItem* pItem = new QListWidgetItem(qstr, listSeq);
	pItem->setData(Qt::UserRole, QVariant::fromValue<void*>(pPos));
}


void FavDlg::AddPosToList() { AddPosToList(m_curPos); }


void FavDlg::RemPosFromList()
{
	QListWidgetItem *pItem = listSeq->currentItem();
	if(pItem)
	{
		FavHklPos* pPos = (FavHklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(pPos) delete pPos;
		delete pItem;
	}
	else
	{
		ClearList();
	}
}

void FavDlg::ClearList()
{
	while(listSeq->count())
	{
		QListWidgetItem *pItem = listSeq->item(0);
		if(!pItem)
			break;

		FavHklPos* pPos = (FavHklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(pPos)
			delete pPos;
		delete pItem;
	}
}


//------------------------------------------------------------------------------


void FavDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	// favlist
	for(int iItem = 0; iItem < listSeq->count(); ++iItem)
	{
		const QListWidgetItem *pItem = listSeq->item(iItem);
		if(!pItem)
			continue;
		const FavHklPos* pPos = (FavHklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(!pPos)
			continue;

		std::ostringstream ostrItemBase;
		ostrItemBase << "fav_pos/pos_" << iItem << "/";
		std::string strItemBase = ostrItemBase.str();

		mapConf[strXmlRoot + strItemBase + "h_start"] = tl::var_to_str(pPos->dhstart);
		mapConf[strXmlRoot + strItemBase + "k_start"] = tl::var_to_str(pPos->dkstart);
		mapConf[strXmlRoot + strItemBase + "l_start"] = tl::var_to_str(pPos->dlstart);
		mapConf[strXmlRoot + strItemBase + "E_start"] = tl::var_to_str(pPos->dEstart);
		mapConf[strXmlRoot + strItemBase + "h_stop"] = tl::var_to_str(pPos->dhstop);
		mapConf[strXmlRoot + strItemBase + "k_stop"] = tl::var_to_str(pPos->dkstop);
		mapConf[strXmlRoot + strItemBase + "l_stop"] = tl::var_to_str(pPos->dlstop);
		mapConf[strXmlRoot + strItemBase + "E_stop"] = tl::var_to_str(pPos->dEstop);
	}
}


void FavDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	ClearList();

	unsigned int iItem = 0;
	while(true)
	{
		std::ostringstream ostrItemBase;
		ostrItemBase << "fav_pos/pos_" << iItem << "/";
		std::string strItemBase = ostrItemBase.str();

		if(!xml.Exists(strXmlRoot + strItemBase))
			break;

		FavHklPos pos;
		pos.dhstart = xml.Query<t_real>(strXmlRoot + strItemBase + "h_start", 0.);
		pos.dkstart = xml.Query<t_real>(strXmlRoot + strItemBase + "k_start", 0.);
		pos.dlstart = xml.Query<t_real>(strXmlRoot + strItemBase + "l_start", 0.);
		pos.dEstart = xml.Query<t_real>(strXmlRoot + strItemBase + "E_start", 0.);
		pos.dhstop = xml.Query<t_real>(strXmlRoot + strItemBase + "h_stop", 0.);
		pos.dkstop = xml.Query<t_real>(strXmlRoot + strItemBase + "k_stop", 0.);
		pos.dlstop = xml.Query<t_real>(strXmlRoot + strItemBase + "l_stop", 0.);
		pos.dEstop = xml.Query<t_real>(strXmlRoot + strItemBase + "E_stop", 0.);

		AddPosToList(pos);
		++iItem;
	}
}


#include "moc_FavDlg.cpp"
