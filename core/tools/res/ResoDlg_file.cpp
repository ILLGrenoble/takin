/**
 * resolution calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013 - 2025
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "ResoDlg.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "tlibs/string/string.h"
#include "tlibs/helper/flags.h"
#include "tlibs/helper/misc.h"
#include "tlibs/time/chrono.h"

#include "libs/globals.h"
#include "libs/version.h"
#include "ellipse.h"

#include <QFileDialog>
#include <QMessageBox>


static const auto angs = tl::get_one_angstrom<t_real_reso>();
static const auto rads = tl::get_one_radian<t_real_reso>();
static const auto meV = tl::get_one_meV<t_real_reso>();
static const auto cm = tl::get_one_centimeter<t_real_reso>();
static const auto meters = tl::get_one_meter<t_real_reso>();
static const auto sec = tl::get_one_second<t_real_reso>();


void ResoDlg::SaveRes()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = ".";
	if(m_pSettings)
		strDirLast = m_pSettings->value("reso/last_dir", ".").toString();
	QString qstrFile = QFileDialog::getSaveFileName(this,
		"Save resolution configuration",
		strDirLast,
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
	bool bOk = xml.Save(strFile, tl::PropType::XML);
	if(!bOk)
		QMessageBox::critical(this, "Error", "Could not save resolution file.");

	if(bOk && m_pSettings)
		m_pSettings->setValue("reso/last_dir", QString(strDir.c_str()));
}


void ResoDlg::LoadRes()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = ".";
	if(m_pSettings)
		strDirLast = m_pSettings->value("reso/last_dir", ".").toString();
	QString qstrFile = QFileDialog::getOpenFileName(this,
		"Open resolution configuration...",
		strDirLast,
		"TAZ files (*.taz *.TAZ)", nullptr,
		fileopt);
	if(qstrFile == "")
		return;


	std::string strFile = qstrFile.toStdString();
	std::string strDir = tl::get_dir(strFile);

	tl::Prop<std::string> xml;
	if(!xml.Load(strFile, tl::PropType::XML))
	{
		QMessageBox::critical(this, "Error", "Could not load resolution file.");
		return;
	}

	Load(xml, strXmlRoot);
	m_strCurDir = strDir;

	if(m_pSettings)
		m_pSettings->setValue("reso/last_dir", QString(strDir.c_str()));

}


void ResoDlg::LoadMonoRefl()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = ".";
	if(m_pSettings)
		strDirLast = m_pSettings->value("reso/last_dir_refl", ".").toString();
	QString qstrFile = QFileDialog::getOpenFileName(this,
		"Open reflectivity file...",
		strDirLast,
		"Data files (*.dat *.DAT);;All files (*.*)", nullptr,
		fileopt);
	if(qstrFile == "")
		return;

	editMonoRefl->setText(qstrFile);
}

void ResoDlg::LoadAnaEffic()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = ".";
	if(m_pSettings)
		strDirLast = m_pSettings->value("reso/last_dir_refl", ".").toString();
	QString qstrFile = QFileDialog::getOpenFileName(this,
		"Open reflectivity file...",
		strDirLast,
		"Data files (*.dat *.DAT);;All files (*.*)", nullptr,
		fileopt);
	if(qstrFile == "")
		return;

	editAnaEffic->setText(qstrFile);
}


// -----------------------------------------------------------------------------


void ResoDlg::WriteLastConfig()
{
	if(!m_pSettings)
		return;

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecSpinBoxes.size(); ++iSpinBox)
		m_pSettings->setValue(m_vecSpinNames[iSpinBox].c_str(), m_vecSpinBoxes[iSpinBox]->value());
	for(std::size_t iSpinBox = 0; iSpinBox < m_vecIntSpinBoxes.size(); ++iSpinBox)
		m_pSettings->setValue(m_vecIntSpinNames[iSpinBox].c_str(), m_vecIntSpinBoxes[iSpinBox]->value());
	for(std::size_t iEditBox = 0; iEditBox < m_vecEditBoxes.size(); ++iEditBox)
		m_pSettings->setValue(m_vecEditNames[iEditBox].c_str(), m_vecEditBoxes[iEditBox]->text());
	for(std::size_t iEditBox = 0; iEditBox < m_vecPosEditBoxes.size(); ++iEditBox)
		m_pSettings->setValue(m_vecPosEditNames[iEditBox].c_str(), m_vecPosEditBoxes[iEditBox]->text().toDouble());
	for(std::size_t iRadio = 0; iRadio < m_vecRadioPlus.size(); ++iRadio)
		m_pSettings->setValue(m_vecRadioNames[iRadio].c_str(), m_vecRadioPlus[iRadio]->isChecked());
	for(std::size_t iCheck = 0; iCheck < m_vecCheckBoxes.size(); ++iCheck)
		m_pSettings->setValue(m_vecCheckNames[iCheck].c_str(), m_vecCheckBoxes[iCheck]->isChecked());
	for(std::size_t iCombo = 0; iCombo < m_vecComboBoxes.size(); ++iCombo)
		m_pSettings->setValue(m_vecComboNames[iCombo].c_str(), m_vecComboBoxes[iCombo]->currentIndex());

	const int iAlgo = static_cast<int>(ResoDlg::GetSelectedAlgo());
		m_pSettings->setValue("reso/algo_idx", iAlgo);

	m_pSettings->setValue("reso/use_guide", groupGuide->isChecked());
}


void ResoDlg::ReadLastConfig()
{
	if(!m_pSettings)
		return;

	bool bOldDontCalc = m_bDontCalc;
	m_bDontCalc = true;

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecSpinBoxes.size(); ++iSpinBox)
	{
		if(!m_pSettings->contains(m_vecSpinNames[iSpinBox].c_str()))
			continue;
		m_vecSpinBoxes[iSpinBox]->setValue(m_pSettings->value(m_vecSpinNames[iSpinBox].c_str()).value<t_real_reso>());
	}

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecIntSpinBoxes.size(); ++iSpinBox)
	{
		if(!m_pSettings->contains(m_vecIntSpinNames[iSpinBox].c_str()))
			continue;
		m_vecIntSpinBoxes[iSpinBox]->setValue(m_pSettings->value(m_vecIntSpinNames[iSpinBox].c_str()).value<int>());
	}

	for(std::size_t iEditBox = 0; iEditBox < m_vecEditBoxes.size(); ++iEditBox)
	{
		if(!m_pSettings->contains(m_vecEditNames[iEditBox].c_str()))
			continue;
		QString strEdit = m_pSettings->value(m_vecEditNames[iEditBox].c_str()).value<QString>();
		m_vecEditBoxes[iEditBox]->setText(strEdit);
	}

	for(std::size_t iEditBox = 0; iEditBox < m_vecPosEditBoxes.size(); ++iEditBox)
	{
		if(!m_pSettings->contains(m_vecPosEditNames[iEditBox].c_str()))
			continue;
		t_real_reso dEditVal = m_pSettings->value(m_vecPosEditNames[iEditBox].c_str()).value<t_real_reso>();
		m_vecPosEditBoxes[iEditBox]->setText(tl::var_to_str(dEditVal, g_iPrec).c_str());
	}

	for(std::size_t iCheckBox = 0; iCheckBox < m_vecCheckBoxes.size(); ++iCheckBox)
	{
		if(!m_pSettings->contains(m_vecCheckNames[iCheckBox].c_str()))
			continue;
		m_vecCheckBoxes[iCheckBox]->setChecked(m_pSettings->value(m_vecCheckNames[iCheckBox].c_str()).value<bool>());
	}

	for(std::size_t iRadio = 0; iRadio < m_vecRadioPlus.size(); ++iRadio)
	{
		if(!m_pSettings->contains(m_vecRadioNames[iRadio].c_str()))
			continue;

		bool bChecked = m_pSettings->value(m_vecRadioNames[iRadio].c_str()).value<bool>();
		if(bChecked)
		{
			m_vecRadioPlus[iRadio]->setChecked(1);
			//m_vecRadioMinus[iRadio]->setChecked(0);;
		}
		else
		{
			//m_vecRadioPlus[iRadio]->setChecked(0);
			m_vecRadioMinus[iRadio]->setChecked(1);;
		}
	}

	for(std::size_t iCombo = 0; iCombo < m_vecComboBoxes.size(); ++iCombo)
	{
		if(!m_pSettings->contains(m_vecComboNames[iCombo].c_str()))
			continue;
		m_vecComboBoxes[iCombo]->setCurrentIndex(
			m_pSettings->value(m_vecComboNames[iCombo].c_str()).value<int>());
	}

	if(m_pSettings->contains("reso/algo_idx"))
		SetSelectedAlgo(static_cast<ResoAlgo>(m_pSettings->value("reso/algo_idx").value<int>()));

	groupGuide->setChecked(m_pSettings->value("reso/use_guide").value<bool>());

	m_bDontCalc = bOldDontCalc;
	RefreshQEPos();
	//Calc();
}


// -----------------------------------------------------------------------------


void ResoDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	{
		// save horizontal src-mono distance as default value for older version
		std::ostringstream ostrVal;
		ostrVal.precision(g_iPrec);
		ostrVal << std::scientific << spinDistHSrcMono->value();

		mapConf[strXmlRoot + std::string("reso/pop_dist_src_mono")] = ostrVal.str();
	}

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecSpinBoxes.size(); ++iSpinBox)
	{
		std::ostringstream ostrVal;
		ostrVal.precision(g_iPrec);
		ostrVal << std::scientific;
		ostrVal << m_vecSpinBoxes[iSpinBox]->value();

		mapConf[strXmlRoot + m_vecSpinNames[iSpinBox]] = ostrVal.str();
	}

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecIntSpinBoxes.size(); ++iSpinBox)
	{
		std::ostringstream ostrVal;
		ostrVal << std::scientific;
		ostrVal << m_vecIntSpinBoxes[iSpinBox]->value();

		mapConf[strXmlRoot + m_vecIntSpinNames[iSpinBox]] = ostrVal.str();
	}

	for(std::size_t iEditBox = 0; iEditBox < m_vecEditBoxes.size(); ++iEditBox)
	{
		std::string strVal = m_vecEditBoxes[iEditBox]->text().toStdString();
		mapConf[strXmlRoot + m_vecEditNames[iEditBox]] = strVal;
	}

	for(std::size_t iEditBox = 0; iEditBox < m_vecPosEditBoxes.size(); ++iEditBox)
	{
		std::string strVal = m_vecPosEditBoxes[iEditBox]->text().toStdString();
		mapConf[strXmlRoot + m_vecPosEditNames[iEditBox]] = strVal;
	}

	for(std::size_t iCheckBox = 0; iCheckBox < m_vecCheckBoxes.size(); ++iCheckBox)
		mapConf[strXmlRoot + m_vecCheckNames[iCheckBox]] = (m_vecCheckBoxes[iCheckBox]->isChecked() ? "1" : "0");

	for(std::size_t iRadio = 0; iRadio < m_vecRadioPlus.size(); ++iRadio)
		mapConf[strXmlRoot + m_vecRadioNames[iRadio]] = (m_vecRadioPlus[iRadio]->isChecked() ? "1" : "0");

	for(std::size_t iCombo = 0; iCombo < m_vecComboBoxes.size(); ++iCombo)
		mapConf[strXmlRoot + m_vecComboNames[iCombo]] = tl::var_to_str<int>(m_vecComboBoxes[iCombo]->currentIndex());

	ResoAlgo algo = ResoDlg::GetSelectedAlgo();
	const int iAlgo = static_cast<int>(algo);
	mapConf[strXmlRoot + "reso/algo_idx"] = tl::var_to_str<int>(static_cast<int>(algo));

	switch(algo)
	{
		case ResoAlgo::CN:
			mapConf[strXmlRoot + "reso/algo"] = "cn";
			break;
		case ResoAlgo::POP_CN:
			mapConf[strXmlRoot + "reso/algo"] = "pop_cn";
			break;
		case ResoAlgo::POP:
			mapConf[strXmlRoot + "reso/algo"] = "pop";
			break;
		case ResoAlgo::ECK:
			mapConf[strXmlRoot + "reso/algo"] = "eck";
			break;
		case ResoAlgo::ECK_EXT:
			mapConf[strXmlRoot + "reso/algo"] = "eck_ext";
			break;
		case ResoAlgo::VIO:
			mapConf[strXmlRoot + "reso/algo"] = "vio";
			break;
		case ResoAlgo::SIMPLE:
			mapConf[strXmlRoot + "reso/algo"] = "simple";
			break;
		default:
			break;
	}

	mapConf[strXmlRoot + "reso/use_guide"] = groupGuide->isChecked() ? "1" : "0";

	const char* pcUser = std::getenv("USER");
	if(!pcUser)
		pcUser = "";
	mapConf[strXmlRoot + "meta/comment"] = textComment->toPlainText().toStdString();
	mapConf[strXmlRoot + "meta/timestamp"] = tl::var_to_str<t_real_reso>(tl::epoch<t_real_reso>());
	mapConf[strXmlRoot + "meta/version"] = TAKIN_VER;
	mapConf[strXmlRoot + "meta/info"] = "Created with Takin/Reso.";
	mapConf[strXmlRoot + "meta/url"] = "https://github.com/ILLGrenoble/takin";
	mapConf[strXmlRoot + "meta/doi"] = "https://dx.doi.org/10.5281/zenodo.4117437";
	mapConf[strXmlRoot + "meta/module"] = "takin/res";
	mapConf[strXmlRoot + "meta/user"] = pcUser;
}


void ResoDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	if(!xml.PathExists(strXmlRoot + "reso"))
	{
		QMessageBox::critical(this, "Error", "Cannot load the selected file as is does not seem to be a Takin/Resolution file.");
		return;
	}

	bool bOldDontCalc = m_bDontCalc;
	m_bDontCalc = true;

	// load src-mono distance value from older version as default before overriding with the new ones
	boost::optional<t_real_reso> odSpinVal = xml.QueryOpt<t_real_reso>(
		strXmlRoot + std::string("reso/pop_dist_src_mono"));
	if(odSpinVal)
	{
		spinDistVSrcMono->setValue(*odSpinVal);
		spinDistHSrcMono->setValue(*odSpinVal);
	}

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecSpinBoxes.size(); ++iSpinBox)
	{
		boost::optional<t_real_reso> odSpinVal = xml.QueryOpt<t_real_reso>(strXmlRoot+m_vecSpinNames[iSpinBox]);
		if(odSpinVal) m_vecSpinBoxes[iSpinBox]->setValue(*odSpinVal);
	}

	for(std::size_t iSpinBox = 0; iSpinBox < m_vecIntSpinBoxes.size(); ++iSpinBox)
	{
		boost::optional<int> odSpinVal = xml.QueryOpt<int>(strXmlRoot+m_vecIntSpinNames[iSpinBox]);
		if(odSpinVal) m_vecIntSpinBoxes[iSpinBox]->setValue(*odSpinVal);
	}

	for(std::size_t iEditBox = 0; iEditBox < m_vecEditBoxes.size(); ++iEditBox)
	{
		boost::optional<std::string> odEditVal = xml.QueryOpt<std::string>(strXmlRoot+m_vecEditNames[iEditBox]);
		if(odEditVal) m_vecEditBoxes[iEditBox]->setText(odEditVal->c_str());
	}

	for(std::size_t iEditBox = 0; iEditBox < m_vecPosEditBoxes.size(); ++iEditBox)
	{
		boost::optional<t_real_reso> odEditVal = xml.QueryOpt<t_real_reso>(strXmlRoot+m_vecPosEditNames[iEditBox]);
		if(odEditVal) m_vecPosEditBoxes[iEditBox]->setText(tl::var_to_str(*odEditVal, g_iPrec).c_str());
	}

	for(std::size_t iCheck = 0; iCheck < m_vecCheckBoxes.size(); ++iCheck)
	{
		boost::optional<int> obChecked = xml.QueryOpt<int>(strXmlRoot+m_vecCheckNames[iCheck]);
		if(obChecked) m_vecCheckBoxes[iCheck]->setChecked(*obChecked);
	}

	for(std::size_t iRadio = 0; iRadio < m_vecRadioPlus.size(); ++iRadio)
	{
		boost::optional<int> obChecked = xml.QueryOpt<int>(strXmlRoot+m_vecRadioNames[iRadio]);
		if(obChecked)
		{
			if(*obChecked)
				m_vecRadioPlus[iRadio]->setChecked(1);
			else
				m_vecRadioMinus[iRadio]->setChecked(1);
		}
	}

	for(std::size_t iCombo = 0; iCombo < m_vecComboBoxes.size(); ++iCombo)
	{
		boost::optional<int> oiComboIdx = xml.QueryOpt<int>(strXmlRoot+m_vecComboNames[iCombo]);
		if(oiComboIdx) m_vecComboBoxes[iCombo]->setCurrentIndex(*oiComboIdx);
	}

	// get reso algo by name
	boost::optional<std::string> opAlgo = xml.QueryOpt<std::string>(strXmlRoot+"reso/algo");
	if(opAlgo)
	{
		if(*opAlgo == "cn")
			SetSelectedAlgo(ResoAlgo::CN);
		else if(*opAlgo == "pop_cn")
			SetSelectedAlgo(ResoAlgo::POP_CN);
		else if(*opAlgo == "pop")
			SetSelectedAlgo(ResoAlgo::POP);
		else if(*opAlgo == "eck")
			SetSelectedAlgo(ResoAlgo::ECK);
		else if(*opAlgo == "eck_ext")
			SetSelectedAlgo(ResoAlgo::ECK_EXT);
		else if(*opAlgo == "vio" || *opAlgo == "viol")
			SetSelectedAlgo(ResoAlgo::VIO);
		else if(*opAlgo == "simple")
			SetSelectedAlgo(ResoAlgo::SIMPLE);
	}
	else
	{
		// get reso algo by index
		boost::optional<int> oiAlgo = xml.QueryOpt<int>(strXmlRoot + "reso/algo_idx");
		if(oiAlgo)
			SetSelectedAlgo(static_cast<ResoAlgo>(*oiAlgo));
	}

	boost::optional<int> obGroupVal = xml.QueryOpt<int>(strXmlRoot+"reso/use_guide");
	if(obGroupVal)
		groupGuide->setChecked(*obGroupVal);

	textComment->setText(xml.Query<std::string>(strXmlRoot+"meta/comment").c_str());
	t_real_reso dTimestamp = xml.Query<t_real_reso>(strXmlRoot+"meta/timestamp");
	editTimestamp->setText(tl::epoch_to_str(dTimestamp).c_str());

	m_bDontCalc = bOldDontCalc;
	RefreshQEPos();
	//Calc();
}



// -----------------------------------------------------------------------------


void ResoDlg::RefreshSimCmd()
{
	const t_real_reso dMin = t_real_reso(tl::get_pi<t_real_reso>()/180./60.);

	std::ostringstream ostrCmd;
	ostrCmd.precision(g_iPrec);

	ostrCmd << "./templateTAS -n 1e6 verbose=1 ";

	ostrCmd << "KI=" << t_real_reso(m_tasparams.ki * angs) << " ";
	ostrCmd << "KF=" << t_real_reso(m_tasparams.kf * angs) << " ";
	ostrCmd << "QM=" << t_real_reso(m_tasparams.Q * angs) << " ";
	ostrCmd << "EN=" << t_real_reso(m_tasparams.E / meV) << " ";
	//ostrCmt << "FX=" << (m_tasparams.bki_fix ? "1" : "2") << " ";

	ostrCmd << "L1=" << t_real_reso(m_tasparams.dist_hsrc_mono / meters) << " ";
	ostrCmd << "L2=" << t_real_reso(m_tasparams.dist_mono_sample / meters) << " ";
	ostrCmd << "L3=" << t_real_reso(m_tasparams.dist_sample_ana / meters) << " ";
	ostrCmd << "L4=" << t_real_reso(m_tasparams.dist_ana_det / meters) << " ";

	ostrCmd << "SM=" << m_tasparams.dmono_sense << " ";
	ostrCmd << "SS=" << m_tasparams.dsample_sense << " ";
	ostrCmd << "SA=" << m_tasparams.dana_sense << " ";

	ostrCmd << "DM=" << t_real_reso(m_tasparams.mono_d / angs) << " ";
	ostrCmd << "DA=" << t_real_reso(m_tasparams.ana_d / angs) << " ";

	ostrCmd << "RMV=" << t_real_reso(m_tasparams.mono_curvv / meters) << " ";
	ostrCmd << "RMH=" << t_real_reso(m_tasparams.mono_curvh / meters) << " ";
	ostrCmd << "RAV=" << t_real_reso(m_tasparams.ana_curvv / meters) << " ";
	ostrCmd << "RAH=" << t_real_reso(m_tasparams.ana_curvh / meters) << " ";

	ostrCmd << "ETAM=" << t_real_reso(m_tasparams.mono_mosaic/rads/dMin) << " ";
	ostrCmd << "ETAA=" << t_real_reso(m_tasparams.ana_mosaic/rads/dMin) << " ";

	ostrCmd << "ALF1=" << t_real_reso(m_tasparams.coll_h_pre_mono/rads/dMin) << " ";
	ostrCmd << "ALF2=" << t_real_reso(m_tasparams.coll_h_pre_sample/rads/dMin) << " ";
	ostrCmd << "ALF3=" << t_real_reso(m_tasparams.coll_h_post_sample/rads/dMin) << " ";
	ostrCmd << "ALF4=" << t_real_reso(m_tasparams.coll_h_post_ana/rads/dMin) << " ";
	ostrCmd << "BET1=" << t_real_reso(m_tasparams.coll_v_pre_mono/rads/dMin) << " ";
	ostrCmd << "BET2=" << t_real_reso(m_tasparams.coll_v_pre_sample/rads/dMin) << " ";
	ostrCmd << "BET3=" << t_real_reso(m_tasparams.coll_v_post_sample/rads/dMin) << " ";
	ostrCmd << "BET4=" << t_real_reso(m_tasparams.coll_v_post_ana/rads/dMin) << " ";

	ostrCmd << "WM=" << t_real_reso(m_tasparams.mono_w / meters) << " ";
	ostrCmd << "HM=" << t_real_reso(m_tasparams.mono_h / meters) << " ";
	ostrCmd << "WA=" << t_real_reso(m_tasparams.ana_w / meters) << " ";
	ostrCmd << "HA=" << t_real_reso(m_tasparams.ana_h / meters) << " ";

	ostrCmd << "NVM=" << m_tasparams.mono_numtiles_v << " ";
	ostrCmd << "NHM=" << m_tasparams.mono_numtiles_h << " ";
	ostrCmd << "NVA=" << m_tasparams.ana_numtiles_v << " ";
	ostrCmd << "NHA=" << m_tasparams.ana_numtiles_h << " ";

	editSim->setPlainText(ostrCmd.str().c_str());
}
