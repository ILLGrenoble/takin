/**
 * TAS tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "taz.h"
#include "tlibs/string/string.h"
#include "tlibs/time/chrono.h"
#include "libs/globals.h"
#include "libs/version.h"
#include "libs/qt/recent.h"
#include "dialogs/FilePreviewDlg.h"

#include <cstdlib>
#include <QMessageBox>
#include <QFileDialog>
#include <boost/scope_exit.hpp>

using t_real = t_real_glob;


/**
 * gets space group index from combo box
 */
static int find_sg_from_combo(QComboBox* pCombo, const std::string& str)
{
	for(int iIdx=0; iIdx<pCombo->count(); ++iIdx)
	{
		xtl::SpaceGroup<t_real> *pSG = reinterpret_cast<xtl::SpaceGroup<t_real>*>
			(pCombo->itemData(iIdx).value<void*>());
		if(pSG && pSG->GetName() == str)
			return iIdx;
	}

	return -1;
}



//--------------------------------------------------------------------------------
// loading/saving

void TazDlg::New()
{
	CrystalOptions crys;
	TriangleOptions triag;

	crys.dLattice[0] = crys.dLattice[1] = crys.dLattice[2] = 5.;
	crys.dLatticeAngles[0] = crys.dLatticeAngles[1] = crys.dLatticeAngles[2] = 90.;
	crys.bChangedLattice = crys.bChangedLatticeAngles = true;

	crys.dPlane1[0] = 1.; crys.dPlane1[1] = 0.; crys.dPlane1[2] = 0.;
	crys.dPlane2[0] = 0.; crys.dPlane2[1] = 1.; crys.dPlane2[2] = 0.;
	crys.bChangedPlane1 = crys.bChangedPlane2 = true;

	crys.strSampleName = " ";
	crys.strSpacegroup = "";
	crys.bChangedSpacegroup = true;

	triag.dAnaD = triag.dMonoD = 3.355;
	triag.bChangedAnaD = triag.bChangedMonoD = true;
	triag.dAnaTwoTheta = triag.dMonoTwoTheta = tl::get_pi<t_real>()/2.;
	triag.bChangedAnaTwoTheta = triag.bChangedMonoTwoTheta = true;

	triag.dTwoTheta = tl::get_pi<t_real>()/2.;
	triag.dMonoTwoTheta = triag.dAnaTwoTheta = tl::get_pi<t_real>()/2.;
	triag.dAngleKiVec0 = tl::get_pi<t_real>()/4.;
	triag.bChangedTwoTheta = triag.bChangedAngleKiVec0 = true;
	triag.bChangedMonoTwoTheta = triag.bChangedAnaTwoTheta = true;

	m_vecAtoms.clear();
	m_vecDarkAngles.clear();
	if(m_sceneReal.GetTasLayout())
		m_sceneReal.GetTasLayout()->SetDarkAngles(&m_vecDarkAngles);

	m_strCurFile = "";
	setWindowTitle(s_strTitle.c_str());

	clear_global_paths();
	DeleteDialogs();
	Disconnect();

	checkSenseM->setChecked(false);
	checkSenseS->setChecked(true);
	checkSenseA->setChecked(false);

	VarsChanged(crys, triag);
}



bool TazDlg::Load()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = m_settings.value("main/last_dir", "~").toString();
	QString strFile = QFileDialog::getOpenFileName(this,
		"Open TAS Configuration...", strDirLast, "Takin files (*.taz *.TAZ)",
		nullptr, fileopt);
	if(strFile == "")
		return false;

	return Load(strFile.toStdString().c_str());
}



bool TazDlg::LoadFile(const QString& strFile)
{
	return Load(strFile.toStdString().c_str());
}



bool TazDlg::Load(const char* pcFile)
{
	if(!tl::file_exists(pcFile))
	{
		tl::log_err("File \"", pcFile, "\" does not exist.");
		return 0;
	}

	m_bReady = false;
	BOOST_SCOPE_EXIT(&m_bReady, &m_sceneReal, &m_sceneRecip)
	{
		m_bReady = true;
		m_sceneReal.GetTasLayout()->SetReady(true);
		m_sceneReal.SetEmitChanges(true);

		m_sceneRecip.GetTriangle()->SetReady(true);
		m_sceneRecip.SetEmitChanges(true);
	} BOOST_SCOPE_EXIT_END

	clear_global_paths();
	Disconnect();

	const std::string strXmlRoot("taz/");

	std::string strFile1 = pcFile;
	std::string strDir = tl::get_dir(strFile1);


	tl::Prop<std::string> xml;
	if(!xml.Load(strFile1.c_str(), tl::PropType::XML))
	{
		std::string strErr = "Could not load file \""
			+ std::string(pcFile) + "\".";
		QMessageBox::critical(this, "Error", strErr.c_str());
		return false;
	}

	if(!xml.PathExists(strXmlRoot))
	{
		QMessageBox::critical(this, "Error", "The selected file is not in the Takin format.");
		return false;
	}

	if(!xml.PathExists(strXmlRoot + "sample"))
		QMessageBox::warning(this, "Warning", "This file does not have a sample definition, attempting to load anyway.");


	m_settings.setValue("main/last_dir", QString(strDir.c_str()));


	bool bOk = false;
	m_sceneReal.SetEmitChanges(false);
	m_sceneReal.GetTasLayout()->SetReady(false);
	m_sceneRecip.SetEmitChanges(false);
	m_sceneRecip.GetTriangle()->SetReady(false);


	// edit boxes
	std::vector<std::vector<QLineEdit*>*> vecEdits
		= {&m_vecEdits_real, &m_vecEdits_recip,
			&m_vecEdits_plane, &m_vecEdits_monoana};
	std::vector<std::vector<std::string>*> vecEditNames
		= {&m_vecEditNames_real, &m_vecEditNames_recip,
			&m_vecEditNames_plane, &m_vecEditNames_monoana};
	for(std::size_t iIdxEdit=0; iIdxEdit<vecEdits.size(); ++iIdxEdit)
	{
		const std::vector<QLineEdit*>* pVec = vecEdits[iIdxEdit];
		const std::vector<std::string>* pvecName = vecEditNames[iIdxEdit];

		for(std::size_t iEditBox=0; iEditBox<pVec->size(); ++iEditBox)
		{
			std::string str = xml.Query<std::string>(strXmlRoot+(*pvecName)[iEditBox], "0", &bOk);
			tl::trim(str);
			if(bOk)
				(*pVec)[iEditBox]->setText(str.c_str());
		}
	}

	std::string strDescr = xml.Query<std::string>(strXmlRoot+"sample/descr", "", &bOk);
	if(bOk)
		editDescr->setText(strDescr.c_str());


	// check boxes
	for(std::size_t iCheckBox=0; iCheckBox<m_vecCheckBoxesSenses.size(); ++iCheckBox)
	{
		int iVal = xml.Query<int>(strXmlRoot+m_vecCheckBoxNamesSenses[iCheckBox], 0, &bOk);
		if(bOk)
			m_vecCheckBoxesSenses[iCheckBox]->setChecked(iVal != 0);
	}


	// TAS Layout
	if(xml.Exists(strXmlRoot + "real"))
	{
		t_real dRealScale = xml.Query<t_real>(strXmlRoot + "real/pixels_per_cm", 0., &bOk);
		if(bOk)
			m_sceneReal.GetTasLayout()->SetScaleFactor(dRealScale);

		for(std::size_t iNodeReal=0; iNodeReal<m_sceneReal.GetTasLayout()->GetNodes().size(); ++iNodeReal)
		{
			TasLayoutNode *pNode = m_sceneReal.GetTasLayout()->GetNodes()[iNodeReal];
			std::string strNode = m_sceneReal.GetTasLayout()->GetNodeNames()[iNodeReal];

			bool bOkX = false, bOkY = false;
			t_real dValX = xml.Query<t_real>(strXmlRoot + "real/" + strNode + "_x", 0., &bOkX);
			t_real dValY = xml.Query<t_real>(strXmlRoot + "real/" + strNode + "_y", 0., &bOkY);

			pNode->setPos(dValX, dValY);
		}

		int bWSEnabled = xml.Query<int>(strXmlRoot + "real/enable_ws", 0, &bOk);
		if(bOk)
			m_pWS->setChecked(bWSEnabled!=0);

		int bRealQEnabled = xml.Query<int>(strXmlRoot + "real/enable_realQDir", 0, &bOk);
		if(bOk)
			m_pShowRealQDir->setChecked(bRealQEnabled!=0);
	}


	// scattering triangle
	if(xml.Exists(strXmlRoot + "recip"))
	{
		t_real dRecipScale = xml.Query<t_real>(strXmlRoot + "recip/pixels_per_A-1", 0., &bOk);
		if(bOk)
			m_sceneRecip.GetTriangle()->SetScaleFactor(dRecipScale);

		for(std::size_t iNodeRecip=0; iNodeRecip<m_sceneRecip.GetTriangle()->GetNodes().size(); ++iNodeRecip)
		{
			ScatteringTriangleNode *pNode = m_sceneRecip.GetTriangle()->GetNodes()[iNodeRecip];
			std::string strNode = m_sceneRecip.GetTriangle()->GetNodeNames()[iNodeRecip];

			bool bOkX=0, bOkY=0;
			t_real dValX = xml.Query<t_real>(strXmlRoot + "recip/" + strNode + "_x", 0., &bOkX);
			t_real dValY = xml.Query<t_real>(strXmlRoot + "recip/" + strNode + "_y", 0., &bOkY);

			pNode->setPos(dValX, dValY);
		}



		int bSmallqEnabled = xml.Query<int>(strXmlRoot + "recip/enable_q", 0, &bOk);
		if(bOk)
			m_pSmallq->setChecked(bSmallqEnabled!=0);

		int bCoordAxesEnabled = xml.Query<int>(strXmlRoot + "recip/enable_axes", 0, &bOk);
		if(bOk)
			m_pCoordAxes->setChecked(bCoordAxesEnabled!=0);

		int bSmallqSnapped = xml.Query<int>(strXmlRoot + "recip/snap_q", 1, &bOk);
		if(bOk)
			m_pSnapSmallq->setChecked(bSmallqSnapped!=0);

		int bBZEnabled = xml.Query<int>(strXmlRoot + "recip/enable_bz", 0, &bOk);
		if(bOk)
			m_pBZ->setChecked(bBZEnabled!=0);

		int bAllPeaks = xml.Query<int>(strXmlRoot + "recip/show_all_peaks", 0, &bOk);
		if(bOk)
			m_pAllPeaks->setChecked(bAllPeaks!=0);

		int iEwald = xml.Query<int>(strXmlRoot + "recip/ewald_sphere", 0, &bOk);
		if(bOk)
		{
			if(iEwald == EWALD_NONE) m_pEwaldSphereNone->setChecked(1);
			else if(iEwald == EWALD_KI) m_pEwaldSphereKi->setChecked(1);
			else if(iEwald == EWALD_KF) m_pEwaldSphereKf->setChecked(1);
		}
	}


	// sample definitions
	if(xml.Exists(strXmlRoot + "sample"))
	{
		std::string strSpaceGroup = xml.Query<std::string>(strXmlRoot + "sample/spacegroup", "", &bOk);
		tl::trim(strSpaceGroup);
		if(bOk)
		{
			editSpaceGroupsFilter->clear();
			RepopulateSpaceGroups();

			int iSGIdx = find_sg_from_combo(comboSpaceGroups, strSpaceGroup);
			if(iSGIdx >= 0)
				comboSpaceGroups->setCurrentIndex(iSGIdx);
			else
				comboSpaceGroups->setCurrentIndex(0);
		}


		m_vecAtoms.clear();
		std::size_t iNumAtoms = xml.Query<std::size_t>(strXmlRoot + "sample/atoms/num", 0, &bOk);
		if(bOk)
		{
			m_vecAtoms.reserve(iNumAtoms);

			for(std::size_t iAtom=0; iAtom<iNumAtoms; ++iAtom)
			{
				xtl::AtomPos<t_real> theatom;
				theatom.vecPos.resize(3,0);

				std::string strNr = tl::var_to_str(iAtom);
				theatom.strAtomName = xml.Query<std::string>(strXmlRoot + "sample/atoms/" + strNr + "/name", "");
				theatom.vecPos[0] = xml.Query<t_real>(strXmlRoot + "sample/atoms/" + strNr + "/x", 0.);
				theatom.vecPos[1] = xml.Query<t_real>(strXmlRoot + "sample/atoms/" + strNr + "/y", 0.);
				theatom.vecPos[2] = xml.Query<t_real>(strXmlRoot + "sample/atoms/" + strNr + "/z", 0.);

				m_vecAtoms.emplace_back(std::move(theatom));
			}

			ShowAtomsDlg(1);
			if(m_pAtomsDlg)
			{
				m_pAtomsDlg->SetAtoms(m_vecAtoms);
				if(!m_pAtomsDlg->ShowPossibleErrorDlg())
					ShowAtomsDlg(0);
			}
		}
	}


	// dark angles
	m_vecDarkAngles.clear();
	if(xml.Exists(strXmlRoot + "darkangles"))
	{
		InitDarkAngles();
		m_pDarkAnglesDlg->Load(xml, strXmlRoot);
		m_vecDarkAngles = m_pDarkAnglesDlg->GetDarkAngles();
		if(m_sceneReal.GetTasLayout())
			m_sceneReal.GetTasLayout()->SetDarkAngles(&m_vecDarkAngles);
	}


	// goto dialog
	if(m_pGotoDlg)
		m_pGotoDlg->ClearList();

	if(xml.Exists(strXmlRoot + "goto_favlist") ||
		xml.Exists(strXmlRoot + "goto_pos"))
	{
		InitGoto();
		m_pGotoDlg->Load(xml, strXmlRoot);
	}


	// elastic positions dialog
	if(xml.Exists(strXmlRoot + "elastic_pos"))
	{
		InitElasticDlg();
		m_pElasticDlg->Load(xml, strXmlRoot);
	}


	// reso dialog
	if(xml.Exists(strXmlRoot + "reso"))
	{
		InitReso();
		m_pReso->SetUpdateOn(0,0);
		m_pReso->Load(xml, strXmlRoot);
	}


	// convo dialog
	if(xml.Exists(strXmlRoot + "monteconvo"))
	{
		InitResoConv();
		m_pConvoDlg->Load(xml, strXmlRoot);
	}


	m_strCurFile = strFile1;
	setWindowTitle((s_strTitle + " - " + m_strCurFile).c_str());

	RecentFiles recent(&m_settings, "main/recent");
	recent.AddFile(strFile1.c_str());
	recent.SaveList();
	recent.FillMenu(m_pMenuRecent, [this](const std::string& str){ LoadFile(str.c_str()); });


	m_bReady = true;

	UpdateDs();
	UpdateMonoSense();
	UpdateAnaSense();
	UpdateSampleSense();

	m_sceneReal.GetTasLayout()->SetReady(true);
	m_sceneReal.SetEmitChanges(true);
	//m_sceneReal.emitUpdate();

	m_sceneRecip.GetTriangle()->SetReady(true);
	m_sceneRecip.SetEmitChanges(true);

	CalcPeaks();
	m_sceneRecip.emitUpdate();

	if(m_pReso)
	{
		m_pReso->SetUpdateOn(1,1);
		m_sceneRecip.emitAllParams();
	}

	return true;
}



bool TazDlg::Save()
{
	if(m_strCurFile == "")
		return SaveAs();

	const std::string strXmlRoot("taz/");
	typedef std::map<std::string, std::string> tmap;
	tmap mapConf;


	// edit boxes
	std::vector<const std::vector<QLineEdit*>*> vecEdits
		= {&m_vecEdits_real, &m_vecEdits_recip, &m_vecEdits_plane, &m_vecEdits_monoana};
	std::vector<const std::vector<std::string>*> vecEditNames
		= {&m_vecEditNames_real, &m_vecEditNames_recip, &m_vecEditNames_plane, &m_vecEditNames_monoana};
	for(std::size_t iIdxEdit=0; iIdxEdit<vecEdits.size(); ++iIdxEdit)
	{
		const std::vector<QLineEdit*>* pVec = vecEdits[iIdxEdit];
		const std::vector<std::string>* pvecName = vecEditNames[iIdxEdit];

		for(std::size_t iEditBox=0; iEditBox<pVec->size(); ++iEditBox)
			mapConf[strXmlRoot+(*pvecName)[iEditBox]]
			        = (*pVec)[iEditBox]->text().toStdString();
	}

	mapConf[strXmlRoot + "sample/descr"] = editDescr->text().toStdString();


	// check boxes
	for(std::size_t iCheckBox=0; iCheckBox<m_vecCheckBoxesSenses.size(); ++iCheckBox)
		mapConf[strXmlRoot+m_vecCheckBoxNamesSenses[iCheckBox]]
			= (m_vecCheckBoxesSenses[iCheckBox]->isChecked() ? "1" : "0");


	// TAS layout
	for(std::size_t iNodeReal=0; iNodeReal<m_sceneReal.GetTasLayout()->GetNodes().size(); ++iNodeReal)
	{
		const TasLayoutNode *pNode = m_sceneReal.GetTasLayout()->GetNodes()[iNodeReal];
		std::string strNode = m_sceneReal.GetTasLayout()->GetNodeNames()[iNodeReal];
		std::string strValX = tl::var_to_str(pNode->pos().x());
		std::string strValY = tl::var_to_str(pNode->pos().y());

		mapConf[strXmlRoot + "real/" + strNode + "_x"] = strValX;
		mapConf[strXmlRoot + "real/" + strNode + "_y"] = strValY;
	}
	t_real dRealScale = m_sceneReal.GetTasLayout()->GetScaleFactor();
	mapConf[strXmlRoot + "real/pixels_per_cm"] = tl::var_to_str(dRealScale);


	// scattering triangle
	for(std::size_t iNodeRecip=0; iNodeRecip<m_sceneRecip.GetTriangle()->GetNodes().size(); ++iNodeRecip)
	{
		const ScatteringTriangleNode *pNode = m_sceneRecip.GetTriangle()->GetNodes()[iNodeRecip];
		std::string strNode = m_sceneRecip.GetTriangle()->GetNodeNames()[iNodeRecip];
		std::string strValX = tl::var_to_str(pNode->pos().x());
		std::string strValY = tl::var_to_str(pNode->pos().y());

		mapConf[strXmlRoot + "recip/" + strNode + "_x"] = strValX;
		mapConf[strXmlRoot + "recip/" + strNode + "_y"] = strValY;
	}
	t_real dRecipScale = m_sceneRecip.GetTriangle()->GetScaleFactor();
	mapConf[strXmlRoot + "recip/pixels_per_A-1"] = tl::var_to_str(dRecipScale);


	bool bSmallqEnabled = m_pSmallq->isChecked();
	mapConf[strXmlRoot + "recip/enable_q"] = (bSmallqEnabled ? "1" : "0");

	bool bCoordAxesEnabled = m_pCoordAxes->isChecked();
	mapConf[strXmlRoot + "recip/enable_axes"] = (bCoordAxesEnabled ? "1" : "0");

	bool bSmallqSnapped = m_sceneRecip.getSnapq();
	mapConf[strXmlRoot + "recip/snap_q"] = (bSmallqSnapped ? "1" : "0");

	bool bBZEnabled = m_pBZ->isChecked();
	mapConf[strXmlRoot + "recip/enable_bz"] = (bBZEnabled ? "1" : "0");

	bool bAllPeaks = m_pAllPeaks->isChecked();
	mapConf[strXmlRoot + "recip/show_all_peaks"] = (bAllPeaks ? "1" : "0");

	int iEw = EWALD_NONE;
	if(m_pEwaldSphereKi->isChecked()) iEw = EWALD_KI;
	else if(m_pEwaldSphereKf->isChecked()) iEw = EWALD_KF;
	mapConf[strXmlRoot + "recip/ewald_sphere"] = tl::var_to_str(iEw);

	bool bWSEnabled = m_pWS->isChecked();
	mapConf[strXmlRoot + "real/enable_ws"] = (bWSEnabled ? "1" : "0");

	bool bRealQDir = m_pShowRealQDir->isChecked();
	mapConf[strXmlRoot + "real/enable_realQDir"] = (bRealQDir ? "1" : "0");


	// space group
	xtl::SpaceGroup<t_real> *pSpaceGroup = nullptr;
	std::string strSG = "-1";
	int iSpaceGroupIdx = comboSpaceGroups->currentIndex();
	if(iSpaceGroupIdx != 0)
		pSpaceGroup = (xtl::SpaceGroup<t_real>*)comboSpaceGroups->itemData(iSpaceGroupIdx).value<void*>();
	if(pSpaceGroup)
		strSG = pSpaceGroup->GetName();
	mapConf[strXmlRoot + "sample/spacegroup"] = strSG;


	// atom positions
	mapConf[strXmlRoot + "sample/atoms/num"] = tl::var_to_str(m_vecAtoms.size());
	for(std::size_t iAtom=0; iAtom<m_vecAtoms.size(); ++iAtom)
	{
		const xtl::AtomPos<t_real>& atom = m_vecAtoms[iAtom];

		std::string strAtomNr = tl::var_to_str(iAtom);
		mapConf[strXmlRoot + "sample/atoms/" + strAtomNr + "/name"] =
			atom.strAtomName;
		mapConf[strXmlRoot + "sample/atoms/" + strAtomNr + "/x"] =
			tl::var_to_str(atom.vecPos[0]);
		mapConf[strXmlRoot + "sample/atoms/" + strAtomNr + "/y"] =
			tl::var_to_str(atom.vecPos[1]);
		mapConf[strXmlRoot + "sample/atoms/" + strAtomNr + "/z"] =
			tl::var_to_str(atom.vecPos[2]);
	}


	// meta data
	const char* pcUser = std::getenv("USER");
	if(!pcUser)
		pcUser = "";

	mapConf[strXmlRoot + "meta/timestamp"] = tl::var_to_str<t_real>(tl::epoch<t_real>());
	mapConf[strXmlRoot + "meta/version"] = TAKIN_VER;
	mapConf[strXmlRoot + "meta/info"] = "Created with Takin.";
	mapConf[strXmlRoot + "meta/url"] = "https://github.com/ILLGrenoble/takin";
	mapConf[strXmlRoot + "meta/doi"] = "https://dx.doi.org/10.5281/zenodo.4117437";
	mapConf[strXmlRoot + "meta/module"] = "takin";
	mapConf[strXmlRoot + "meta/user"] = pcUser;


	// dialogs
	if(m_pReso) m_pReso->Save(mapConf, strXmlRoot);
	if(m_pConvoDlg) m_pConvoDlg->Save(mapConf, strXmlRoot);
	if(m_pGotoDlg) m_pGotoDlg->Save(mapConf, strXmlRoot);
	if(m_pElasticDlg) m_pElasticDlg->Save(mapConf, strXmlRoot);
	if(m_pDarkAnglesDlg) m_pDarkAnglesDlg->Save(mapConf, strXmlRoot);
	//if(m_pPowderDlg) m_pPowderDlg->Save(mapConf, strXmlRoot);


	tl::Prop<std::string> xml;
	xml.Add(mapConf);
	if(!xml.Save(m_strCurFile.c_str(), tl::PropType::XML))
	{
		QMessageBox::critical(this, "Error", "Could not save configuration file.");
		return false;
	}

	RecentFiles recent(&m_settings, "main/recent");
	recent.AddFile(m_strCurFile.c_str());
	recent.SaveList();
	recent.FillMenu(m_pMenuRecent, [this](const std::string& str){ LoadFile(str.c_str()); });

	return true;
}



bool TazDlg::SaveAs()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = m_settings.value("main/last_dir", "~").toString();
	QString strFile = QFileDialog::getSaveFileName(this,
		"Save TAS Configuration", strDirLast, "Takin files (*.taz *.TAZ)",
		nullptr, fileopt);

	if(strFile != "")
	{
		std::string strFile1 = strFile.toStdString();
		std::string strDir = tl::get_dir(strFile1);
		if(tl::get_fileext(strFile1,1) != "taz")
			strFile1 += ".taz";

		m_strCurFile = strFile1;
		setWindowTitle((s_strTitle + " - " + m_strCurFile).c_str());
		bool bOk = Save();

		if(bOk)
			m_settings.setValue("main/last_dir", QString(strDir.c_str()));

		return bOk;
	}

	return false;
}




//--------------------------------------------------------------------------------
// data file importing

#include "tlibs/file/loadinstr.h"

bool TazDlg::Import()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	const bool bShowPreview = m_settings.value("main/dlg_previews", true).toBool();
	QString strDirLast = m_settings.value("main/last_import_dir", "~").toString();

	std::unique_ptr<QFileDialog> pdlg;
	if(bShowPreview)
		pdlg.reset(new FilePreviewDlg(this, "Import Data File...", &m_settings));
	else
		pdlg.reset(new QFileDialog(this, "Import Data File..."));

	pdlg->setOptions(fileopt);
	pdlg->setDirectory(strDirLast);
	pdlg->setFileMode(QFileDialog::ExistingFile);
	pdlg->setViewMode(QFileDialog::Detail);
#if !defined NO_IOSTR
	QString strFilter = "Data files (*.dat *.scn *.DAT *.SCN *.ng0 *.NG0 *.log *.LOG *.scn.gz *.SCN.GZ *.dat.gz *.DAT.GZ *.ng0.gz *.NG0.GZ *.log.gz *.LOG.GZ *.scn.bz2 *.SCN.BZ2 *.dat.bz2 *.DAT.BZ2 *.ng0.bz2 *.NG0.BZ2 *.log.bz2 *.LOG.BZ2);;All files (*.* *)";
#else
	QString strFilter = "Data files (*.dat *.scn *.DAT *.SCN *.NG0 *.ng0 *.log *.LOG);;All files (*.* *)";
#endif
	pdlg->setNameFilter(strFilter);
	if(!pdlg->exec())
		return false;
	if(!pdlg->selectedFiles().size())
		return false;

	QString strFile = pdlg->selectedFiles()[0];
	if(strFile == "")
		return false;

	return Import(strFile.toStdString().c_str());
}


bool TazDlg::ImportFile(const QString& strFile)
{
	return Import(strFile.toStdString().c_str());
}


bool TazDlg::Import(const char* pcFile)
{
	try
	{
		Disconnect();

		std::string strFile1 = pcFile;
		std::string strDir = tl::get_dir(strFile1);

		std::size_t iScanNum = 0;


		std::unique_ptr<tl::FileInstrBase<t_real>> ptrDat(
			tl::FileInstrBase<t_real>::LoadInstr(pcFile));
		tl::FileInstrBase<t_real>* pdat = ptrDat.get();
		if(!pdat)
			return false;

		std::array<t_real, 3> arrLatt = pdat->GetSampleLattice();
		std::array<t_real, 3> arrAng = pdat->GetSampleAngles();
		std::array<bool, 3> arrSenses = pdat->GetScatterSenses();
		std::array<t_real, 2> arrD = pdat->GetMonoAnaD();
		std::array<t_real, 3> arrPeak0 = pdat->GetScatterPlane0();
		std::array<t_real, 3> arrPeak1 = pdat->GetScatterPlane1();

		editA->setText(tl::var_to_str(arrLatt[0]).c_str());
		editB->setText(tl::var_to_str(arrLatt[1]).c_str());
		editC->setText(tl::var_to_str(arrLatt[2]).c_str());

		editAlpha ->setText(tl::var_to_str(tl::r2d(arrAng[0])).c_str());
		editBeta->setText(tl::var_to_str(tl::r2d(arrAng[1])).c_str());
		editGamma->setText(tl::var_to_str(tl::r2d(arrAng[2])).c_str());

		editMonoD->setText(tl::var_to_str(arrD[0]).c_str());
		editAnaD->setText(tl::var_to_str(arrD[1]).c_str());

		checkSenseM->setChecked(arrSenses[0]);
		checkSenseS->setChecked(arrSenses[1]);
		checkSenseA->setChecked(arrSenses[2]);

		editScatX0->setText(tl::var_to_str(arrPeak0[0]).c_str());
		editScatX1->setText(tl::var_to_str(arrPeak0[1]).c_str());
		editScatX2->setText(tl::var_to_str(arrPeak0[2]).c_str());

		editScatY0->setText(tl::var_to_str(arrPeak1[0]).c_str());
		editScatY1->setText(tl::var_to_str(arrPeak1[1]).c_str());
		editScatY2->setText(tl::var_to_str(arrPeak1[2]).c_str());

		// spacegroup
		editSpaceGroupsFilter->clear();
		RepopulateSpaceGroups();

		std::string strSpaceGroup = pdat->GetSpacegroup();
		tl::trim(strSpaceGroup);

		int iSGIdx = find_sg_from_combo(comboSpaceGroups, strSpaceGroup);
		if(iSGIdx >= 0)
			comboSpaceGroups->setCurrentIndex(iSGIdx);
		else
			comboSpaceGroups->setCurrentIndex(0);

		// descr
		std::string strExp = pdat->GetTitle();
		std::string strSample = pdat->GetSampleName();
		if(strSample != "")
			strExp += std::string(" - ") + strSample;
		editDescr->setText(strExp.c_str());

		iScanNum = pdat->GetScanCount();
		if(iScanNum)
		{
			InitGoto();
			m_pGotoDlg->ClearList();

			for(std::size_t iScan=0; iScan<iScanNum; ++iScan)
			{
				std::array<t_real, 5> arrScan = pdat->GetScanHKLKiKf(iScan);
				m_pGotoDlg->AddPosToList(arrScan[0],arrScan[1],arrScan[2],arrScan[3],arrScan[4]);
			}
		}


		m_settings.setValue("main/last_import_dir", QString(strDir.c_str()));
		m_strCurFile = /*strFile1*/ "";		// prevents overwriting imported file on saving
		setWindowTitle((s_strTitle + " - " + strFile1).c_str());

		RecentFiles recent(&m_settings, "main/recent_import");
		recent.AddFile(strFile1.c_str());
		recent.SaveList();
		recent.FillMenu(m_pMenuRecentImport, [this](const std::string& str){ ImportFile(str.c_str()); });

		CalcPeaks();

		if(iScanNum && m_pGotoDlg)
		{
			if(m_pGotoDlg->GotoPos(0))
				focus_dlg(m_pGotoDlg);
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
		return false;
	}

	return true;
}

//--------------------------------------------------------------------------------



//--------------------------------------------------------------------------------
// CIF importing

#include "libs/spacegroups/xtl_xml.h"
#include "tlibs/helper/proc.h"


bool TazDlg::ImportCIF()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	const bool bShowPreview = m_settings.value("main/dlg_previews", true).toBool();
	QString strDirLast = m_settings.value("main/last_import_cif_dir", "~").toString();

	std::unique_ptr<QFileDialog> pdlg{new QFileDialog(this, "Import CIF...")};

	pdlg->setOptions(fileopt);
	pdlg->setDirectory(strDirLast);
	pdlg->setFileMode(QFileDialog::ExistingFile);
	pdlg->setViewMode(QFileDialog::Detail);
#if !defined NO_IOSTR
	QString strFilter = "CIFs (*.cif *.CIF *.cif.gz *.CIF.GZ *.cif.bz2 *.CIF.BZ2);;All files (*.*)";
#else
	QString strFilter = "CIFs (*.cif *.CIF);;All files (*.*)";
#endif
	pdlg->setNameFilter(strFilter);
	if(!pdlg->exec())
		return false;
	if(!pdlg->selectedFiles().size())
		return false;

	QString strFile = pdlg->selectedFiles()[0];
	if(strFile == "")
		return false;

	return ImportCIF(strFile.toStdString().c_str());
}



bool TazDlg::ImportCIFFile(const QString& strFile)
{
	return ImportCIF(strFile.toStdString().c_str());
}



static inline std::pair<std::string, std::string> get_ciftool_version()
{
	std::string cifbin = find_program_binary("takin_cif2xml");
	if(cifbin == "")
		return std::make_pair("", cifbin);

	tl::PipeProc<char> proc(("\"" + cifbin + "\" 2>/dev/null").c_str(), false);
	if(!proc.IsReady())
		return std::make_pair("", cifbin);

	std::string strVer;
	std::getline(proc.GetIstr(), strVer);
	tl::trim(strVer);

	return std::make_pair(strVer, cifbin);
}



bool TazDlg::ImportCIF(const char* pcFile)
{
	try
	{
		Disconnect();

		std::string strFile1 = pcFile;
		std::string strDir = tl::get_dir(strFile1);


		// check cif tool
		const std::string strToolNotFound =
			"The external Cif2Xml tool is not correctly installed. "
			"Please specify its path in the program settings.";

		std::string strCifVer;
		std::string strCifBin;
		std::tie(strCifVer, strCifBin) = get_ciftool_version();
		if(strCifVer == "")
		{
			QMessageBox::critical(this, "Error", strToolNotFound.c_str());
			return false;
		}


		// open pipe to external cif2xml tool
		tl::log_info("Invoking ", strCifVer);
		tl::PipeProc<char> proc((strCifBin + " 2>/dev/null \"" + strFile1 + "\"").c_str(), false);
		if(!proc.IsReady())
		{
			QMessageBox::critical(this, "Error", strToolNotFound.c_str());
			return false;
		}


		// load xml representation of cif
		using t_vec = ublas::vector<t_real>;
		using t_mat = ublas::matrix<t_real>;

		bool bOk = false;
		t_real a = 5, b = 5, c = 5;
		t_real alpha = M_PI/2, beta = M_PI/2, gamma = M_PI/2;
		std::vector<std::string> vecAtomNames;
		std::vector<t_vec> vecAtomPos;
		std::vector<std::vector<t_vec>> vecAllAtomPos;
		std::string strSpaceGroup = "";

		//std::ifstream istrCIF(strFile1);
		std::basic_istream<char>& istrCIF = proc.GetIstr();
		std::tie(bOk, a,b,c, alpha,beta,gamma, vecAtomNames, vecAtomPos, vecAllAtomPos, strSpaceGroup) =
			xtl::load_xml<t_real, t_vec, t_mat>(istrCIF);
		tl::trim(strSpaceGroup);

		if(!bOk)
		{
			QMessageBox::critical(this, "Error", "CIF import failed");
			return false;
		}


		// TODO: clean up atom names
		//for(const std::string& strAtom : vecAtomNames)
		//{
		//}


		// update atomic positions
		m_vecAtoms.clear();

		if(strSpaceGroup == "")
		{
			QMessageBox::warning(this, "Warning",
				"No suitable space group could be found which matches the CIF, "
				"using P1 instead and generating all atomic positions.");

			strSpaceGroup = "P1";
			for(std::size_t iAtomType=0; iAtomType<vecAllAtomPos.size(); ++iAtomType)
			{
				for(std::size_t iAtom=0; iAtom<vecAllAtomPos[iAtomType].size(); ++iAtom)
				{
					xtl::AtomPos<t_real> theatom;
					theatom.vecPos = vecAllAtomPos[iAtomType][iAtom];
					theatom.strAtomName = vecAtomNames[iAtomType];

					m_vecAtoms.emplace_back(std::move(theatom));
				}
			}
		}
		else
		{
			for(std::size_t iAtom=0; iAtom<vecAtomPos.size(); ++iAtom)
			{
				xtl::AtomPos<t_real> theatom;
				theatom.vecPos = vecAtomPos[iAtom];
				theatom.strAtomName = vecAtomNames[iAtom];

				m_vecAtoms.emplace_back(std::move(theatom));
			}
		}

		ShowAtomsDlg(1);
		if(m_pAtomsDlg)
		{
			m_pAtomsDlg->SetAtoms(m_vecAtoms);
			if(!m_pAtomsDlg->ShowPossibleErrorDlg())
				ShowAtomsDlg(0);
		}


		// lattice
		editA->setText(tl::var_to_str(a).c_str());
		editB->setText(tl::var_to_str(b).c_str());
		editC->setText(tl::var_to_str(c).c_str());

		editAlpha ->setText(tl::var_to_str(tl::r2d(alpha)).c_str());
		editBeta->setText(tl::var_to_str(tl::r2d(beta)).c_str());
		editGamma->setText(tl::var_to_str(tl::r2d(gamma)).c_str());


		// spacegroup
		editSpaceGroupsFilter->clear();
		RepopulateSpaceGroups();

		int iSGIdx = find_sg_from_combo(comboSpaceGroups, strSpaceGroup);
		if(iSGIdx >= 0)
			comboSpaceGroups->setCurrentIndex(iSGIdx);
		else
			comboSpaceGroups->setCurrentIndex(0);


		// descr
		editDescr->setText("");


		m_settings.setValue("main/last_import_cif_dir", QString(strDir.c_str()));
		m_strCurFile = /*strFile1*/ "";		// prevents overwriting imported file on saving
		setWindowTitle((s_strTitle + " - " + strFile1).c_str());

		RecentFiles recent(&m_settings, "main/recent_import_cif");
		recent.AddFile(strFile1.c_str());
		recent.SaveList();
		recent.FillMenu(m_pMenuRecentImportCIF, [this](const std::string& str){ ImportCIFFile(str.c_str()); });

		CalcPeaks();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
		return false;
	}

	return true;
}
