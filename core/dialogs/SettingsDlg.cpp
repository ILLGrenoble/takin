/**
 * settings
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 5-dec-2014
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

#include "SettingsDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/log/log.h"
#include "tlibs/file/file.h"
#include "tlibs/math/math.h"
#ifndef NO_3D
	#include "tlibs/gfx/gl.h"
#endif
#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "libs/qt/qthelper.h"

#include <QFileDialog>
#include <QFontDialog>
#include <QStyleFactory>
#include <QPushButton>

#include <iostream>
#include <limits>
#include <thread>

using t_real = t_real_glob;


// ellipse table rows
#define ELLI_1A       0
#define ELLI_1B       1
#define ELLI_2A       2
#define ELLI_2B       3
#define ELLI_3A       4
#define ELLI_3B       5
#define ELLI_4A       6
#define ELLI_4B       7

// ellipse table columns
#define ELLI_X        0
#define ELLI_Y        1
#define ELLI_PROJ1    2
#define ELLI_REM1     3
#define ELLI_PROJ2    4
#define ELLI_REM2     5

// ellipsoid table rows
#define ELLO_1        0
#define ELLO_2        1

// ellipsoid table columns
#define ELLO_X        0
#define ELLO_Y        1
#define ELLO_Z        2
#define ELLO_PROJREM  3
// -----------------------------------------------------------------------------



SettingsDlg::SettingsDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	setupUi(this);

	g_fontGen.setStyleHint(QFont::SansSerif);
	g_fontGfx.setStyleHint(QFont::SansSerif, QFont::PreferAntialias);
	setFont(g_fontGen);

	// get possible GUI styles
	for(const QString& strStyle : QStyleFactory::keys())
		comboGUI->addItem(strStyle);

	connect(buttonBox, &QDialogButtonBox::clicked, this, &SettingsDlg::ButtonBoxClicked);
	connect(btnGLFont, &QAbstractButton::clicked, this, &SettingsDlg::SelectGLFont);
	connect(btnGfxFont, &QAbstractButton::clicked, this, &SettingsDlg::SelectGfxFont);
	connect(btnGenFont, &QAbstractButton::clicked, this, &SettingsDlg::SelectGenFont);
	connect(btnGpl, &QAbstractButton::clicked, this, &SettingsDlg::SelectGplTool);

	m_vecEdits =
	{
		// nicos devices
		t_tupEdit("net/sample_name", "nicos/sample/samplename", editSampleName),
		t_tupEdit("net/lattice", "nicos/sample/lattice", editSampleLattice),
		t_tupEdit("net/angles", "nicos/sample/angles", editSampleAngles),
		t_tupEdit("net/orient1", "nicos/sample/orient1", editSampleOrient1),
		t_tupEdit("net/orient2", "nicos/sample/orient2", editSampleOrient2),
		t_tupEdit("net/spacegroup", "nicos/sample/spacegroup", editSampleSG),
		t_tupEdit("net/psi0", "nicos/sample/psi0", editSamplePsi0),

		t_tupEdit("net/stheta", "nicos/sth/value", editSampleTheta),
		t_tupEdit("net/s2theta", "nicos/stt/value", editSample2Theta),

		t_tupEdit("net/mtheta", "nicos/mth/value", editMonoTheta),
		t_tupEdit("net/m2theta", "nicos/mtt/value", editMono2Theta),
		t_tupEdit("net/mono_d", "nicos/mono/dvalue", editMonoD),

		t_tupEdit("net/atheta", "nicos/ath/value", editAnaTheta),
		t_tupEdit("net/a2theta", "nicos/att/value", editAna2Theta),
		t_tupEdit("net/ana_d", "nicos/ana/dvalue", editAnaD),

		t_tupEdit("net/timer", "nicos/timer/value", editCurTime),
		t_tupEdit("net/preset", "nicos/timer/preselection", editPreset),
		t_tupEdit("net/counter", "nicos/ctr1/value", editCounter),

		// sics devices
		t_tupEdit("net/lattice_a_6", "AS", edit6_AS),
		t_tupEdit("net/lattice_b_6", "BS", edit6_BS),
		t_tupEdit("net/lattice_c_6", "CS", edit6_CS),

		t_tupEdit("net/angle_a_6", "AA", edit6_AA),
		t_tupEdit("net/angle_b_6", "BB", edit6_BB),
		t_tupEdit("net/angle_c_6", "CC", edit6_CC),

		t_tupEdit("net/orient1_x_6", "AX", edit6_AX),
		t_tupEdit("net/orient1_y_6", "AY", edit6_AY),
		t_tupEdit("net/orient1_z_6", "AZ", edit6_AZ),

		t_tupEdit("net/orient2_x_6", "BX", edit6_BX),
		t_tupEdit("net/orient2_y_6", "BY", edit6_BY),
		t_tupEdit("net/orient2_z_6", "BZ", edit6_BZ),

		t_tupEdit("net/stheta_6", "A3", edit6_A3),
		t_tupEdit("net/s2theta_6", "A4", edit6_A4),

		t_tupEdit("net/mtheta_6", "A1", edit6_A1),
		t_tupEdit("net/m2theta_6", "A2", edit6_A2),
		t_tupEdit("net/mono_d_6", "DM", edit6_DM),

		t_tupEdit("net/atheta_6", "A5", edit6_A5),
		t_tupEdit("net/a2theta_6", "A6", edit6_A6),
		t_tupEdit("net/ana_d_6", "DA", edit6_DA),

		t_tupEdit("net/timer_6", "counter getmonitor 1", edit6_CurTime),
		t_tupEdit("net/preset_6", "counter getpreset", edit6_Preset),
		t_tupEdit("net/counter_6", "counter getcounts", edit6_Counter),
		t_tupEdit("net/xdat_6", "iscan getvardata", edit6_xDat),
		t_tupEdit("net/ydat_6", "iscan getcounts", edit6_yDat),

		// misc
		t_tupEdit("gl/font", g_strFontGL.c_str(), editGLFont),
		t_tupEdit("main/font_gfx", g_fontGfx.toString().toStdString().c_str(), editGfxFont),
		t_tupEdit("main/font_gen", g_fontGen.toString().toStdString().c_str(), editGenFont),

		// external tools
		t_tupEdit("tools/gpl", "gnuplot", editGpl)
	};

	m_vecChecks =
	{
		t_tupCheck("main/dlg_previews", 1, checkPreview),
		t_tupCheck("main/native_dialogs", 1, checkNativeDlg),
		t_tupCheck("net/flip_orient2", 1, checkFlipOrient2),
		t_tupCheck("net/sth_stt_corr", 0, checkSthSttCorr),
		t_tupCheck("main/ignore_xtal_restrictions", 0, checkIgnoreXtalRestrictions),
		t_tupCheck("main/allow_scan_merging", 0, checkAllowScanMerging),
	};

	m_vecSpins =
	{
		t_tupSpin("main/prec", g_iPrec, spinPrecGen),
		t_tupSpin("main/prec_gfx", g_iPrecGfx, spinPrecGfx),
		t_tupSpin("main/points_gfx", GFX_NUM_POINTS, spinPtsGfx),
		t_tupSpin("main/max_neighbours", g_iMaxNN, spinMaxNN),
		t_tupSpin("main/max_peaks", 10, spinBragg),
		t_tupSpin("main/max_threads", g_iMaxThreads, spinThreads),
		t_tupSpin("main/max_processes", g_iMaxProcesses, spinProcesses),
		t_tupSpin("gl/font_size", 24, spinGLFont),
		t_tupSpin("net/poll", 750, spinNetPoll),
	};

	m_vecCombos =
	{
		t_tupCombo("main/sfact_sq", 0, comboSFact),
		t_tupCombo("main/calc_3d_bz", 0, comboBZ),
		t_tupCombo("main/gui_style", 0, comboGUI),
	};

	spinPrecGen->setMaximum(std::numeric_limits<t_real>::max_digits10);
	spinPrecGfx->setMaximum(std::numeric_limits<t_real>::max_digits10);
	spinThreads->setMaximum(std::thread::hardware_concurrency());
	spinProcesses->setMaximum(std::thread::hardware_concurrency());


	// tables
	m_elli_rows = { ELLI_1A, ELLI_1B, ELLI_2A, ELLI_2B, ELLI_3A, ELLI_3B, ELLI_4A, ELLI_4B };
	m_ello_rows = { ELLO_1, ELLO_2 };

	tabEllipses->setColumnWidth(ELLI_X, 32);
	tabEllipses->setColumnWidth(ELLI_Y, 32);
	tabEllipses->setColumnWidth(ELLI_PROJ1, 64);
	tabEllipses->setColumnWidth(ELLI_REM1, 64);
	tabEllipses->setColumnWidth(ELLI_PROJ2, 64);
	tabEllipses->setColumnWidth(ELLI_REM2, 64);

	tabEllipsoids->setColumnWidth(ELLO_X, 48);
	tabEllipsoids->setColumnWidth(ELLO_Y, 48);
	tabEllipsoids->setColumnWidth(ELLO_Z, 48);
	tabEllipsoids->setColumnWidth(ELLO_PROJREM, 64);

	for(int elli_row : m_elli_rows)
	{
		tabEllipses->setRowHeight(elli_row,
			tabEllipses->verticalHeader()->minimumSectionSize() + 4);
	}

	for(int ello_row : m_ello_rows)
	{
		tabEllipsoids->setRowHeight(ello_row,
			tabEllipsoids->verticalHeader()->minimumSectionSize() + 4);
	}


	if(m_pSettings && m_pSettings->contains("settings/geo"))
		restoreGeometry(m_pSettings->value("settings/geo").toByteArray());

	SetDefaults(false);
	LoadSettings();
}



SettingsDlg::~SettingsDlg()
{}



void SettingsDlg::SetDefaults(bool bOverwrite)
{
	if(!m_pSettings)
		return;

	for(const t_tupEdit& tup : m_vecEdits)
	{
		const std::string& strKey = std::get<0>(tup);
		const std::string& strDef = std::get<1>(tup);

		bool bKeyExists = m_pSettings->contains(strKey.c_str());
		if(bKeyExists && !bOverwrite)
			continue;

		m_pSettings->setValue(strKey.c_str(), strDef.c_str());
	}

	for(const t_tupCheck& tup : m_vecChecks)
	{
		const std::string& strKey = std::get<0>(tup);
		const bool bDef = std::get<1>(tup);

		bool bKeyExists = m_pSettings->contains(strKey.c_str());
		if(bKeyExists && !bOverwrite)
			continue;

		m_pSettings->setValue(strKey.c_str(), bDef);
	}

	for(const t_tupSpin& tup : m_vecSpins)
	{
		const std::string& strKey = std::get<0>(tup);
		const int iDef = std::get<1>(tup);

		bool bKeyExists = m_pSettings->contains(strKey.c_str());
		if(bKeyExists && !bOverwrite)
			continue;

		m_pSettings->setValue(strKey.c_str(), iDef);
	}

	for(const t_tupCombo& tup : m_vecCombos)
	{
		const std::string& strKey = std::get<0>(tup);
		const int iDef = std::get<1>(tup);

		bool bKeyExists = m_pSettings->contains(strKey.c_str());
		if(bKeyExists && !bOverwrite)
			continue;

		m_pSettings->setValue(strKey.c_str(), iDef);
	}


	// ellipse table
	if(bOverwrite)
	{
		m_pSettings->setValue("reso/ellipse_1a_x", 0); m_pSettings->setValue("reso/ellipse_1b_x", 0);
		m_pSettings->setValue("reso/ellipse_2a_x", 1); m_pSettings->setValue("reso/ellipse_2b_x", 1);
		m_pSettings->setValue("reso/ellipse_3a_x", 2); m_pSettings->setValue("reso/ellipse_3b_x", 2);
		m_pSettings->setValue("reso/ellipse_4a_x", 0); m_pSettings->setValue("reso/ellipse_4b_x", 0);

		m_pSettings->setValue("reso/ellipse_1a_y", 3); m_pSettings->setValue("reso/ellipse_1b_y", 3);
		m_pSettings->setValue("reso/ellipse_2a_y", 3); m_pSettings->setValue("reso/ellipse_2b_y", 3);
		m_pSettings->setValue("reso/ellipse_3a_y", 3); m_pSettings->setValue("reso/ellipse_3b_y", 3);
		m_pSettings->setValue("reso/ellipse_4a_y", 1); m_pSettings->setValue("reso/ellipse_4b_y", 1);

		m_pSettings->setValue("reso/ellipse_1a_proj1", 1); m_pSettings->setValue("reso/ellipse_1b_proj1", -1);
		m_pSettings->setValue("reso/ellipse_2a_proj1", 0); m_pSettings->setValue("reso/ellipse_2b_proj1", -1);
		m_pSettings->setValue("reso/ellipse_3a_proj1", 0); m_pSettings->setValue("reso/ellipse_3b_proj1", -1);
		m_pSettings->setValue("reso/ellipse_4a_proj1", 3); m_pSettings->setValue("reso/ellipse_4b_proj1", -1);

		m_pSettings->setValue("reso/ellipse_1a_rem1", 2); m_pSettings->setValue("reso/ellipse_1b_rem1", 2);
		m_pSettings->setValue("reso/ellipse_2a_rem1", 2); m_pSettings->setValue("reso/ellipse_2b_rem1", 2);
		m_pSettings->setValue("reso/ellipse_3a_rem1", 1); m_pSettings->setValue("reso/ellipse_3b_rem1", 1);
		m_pSettings->setValue("reso/ellipse_4a_rem1", 2); m_pSettings->setValue("reso/ellipse_4b_rem1", 2);

		m_pSettings->setValue("reso/ellipse_1a_proj2", -1); m_pSettings->setValue("reso/ellipse_1b_proj2", -1);
		m_pSettings->setValue("reso/ellipse_2a_proj2", -1); m_pSettings->setValue("reso/ellipse_2b_proj2", -1);
		m_pSettings->setValue("reso/ellipse_3a_proj2", -1); m_pSettings->setValue("reso/ellipse_3b_proj2", -1);
		m_pSettings->setValue("reso/ellipse_4a_proj2", -1); m_pSettings->setValue("reso/ellipse_4b_proj2", -1);

		m_pSettings->setValue("reso/ellipse_1a_rem2", -1); m_pSettings->setValue("reso/ellipse_1b_rem2", 1);
		m_pSettings->setValue("reso/ellipse_2a_rem2", -1); m_pSettings->setValue("reso/ellipse_2b_rem2", 0);
		m_pSettings->setValue("reso/ellipse_3a_rem2", -1); m_pSettings->setValue("reso/ellipse_3b_rem2", 0);
		m_pSettings->setValue("reso/ellipse_4a_rem2", -1); m_pSettings->setValue("reso/ellipse_4b_rem2", 3);
	}


	// ellipsoid table
	if(bOverwrite)
	{
		m_pSettings->setValue("reso/ellipsoid3d_1_x", 0);
		m_pSettings->setValue("reso/ellipsoid3d_2_x", 0);

		m_pSettings->setValue("reso/ellipsoid3d_1_y", 1);
		m_pSettings->setValue("reso/ellipsoid3d_2_y", 1);

		m_pSettings->setValue("reso/ellipsoid3d_1_z", 3);
		m_pSettings->setValue("reso/ellipsoid3d_2_z", 2);

		m_pSettings->setValue("reso/ellipsoid3d_1_proj_or_rem", 2);
		m_pSettings->setValue("reso/ellipsoid3d_2_proj_or_rem", 3);
	}
}



bool SettingsDlg::LoadSettings(QSettings *sett)
{
	// use global settings if none given
	if(!sett)
		sett = m_pSettings;
	if(!sett)
		return false;


	for(const t_tupEdit& tup : m_vecEdits)
	{
		const std::string& strKey = std::get<0>(tup);
		const std::string& strDef = std::get<1>(tup);
		QLineEdit* pEdit = std::get<2>(tup);

		QString strVal = sett->value(strKey.c_str(), strDef.c_str()).toString();
		pEdit->setText(strVal);
	}

	for(const t_tupCheck& tup : m_vecChecks)
	{
		const std::string& strKey = std::get<0>(tup);
		bool bDef = std::get<1>(tup);
		QCheckBox* pCheck = std::get<2>(tup);

		bool bVal = sett->value(strKey.c_str(), bDef).toBool();
		pCheck->setChecked(bVal);
	}

	for(const t_tupSpin& tup : m_vecSpins)
	{
		const std::string& strKey = std::get<0>(tup);
		int iDef = std::get<1>(tup);
		QSpinBox* pSpin = std::get<2>(tup);

		int iVal = sett->value(strKey.c_str(), iDef).toInt();
		pSpin->setValue(iVal);
	}

	for(const t_tupCombo& tup : m_vecCombos)
	{
		const std::string& strKey = std::get<0>(tup);
		int iDef = std::get<1>(tup);
		QComboBox* pCombo = std::get<2>(tup);

		int iVal = sett->value(strKey.c_str(), iDef).toInt();
		pCombo->setCurrentIndex(iVal);
	}


	// ellipse table
	std::vector<std::string> elli_names = { "1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b" };

	for(std::size_t row = 0; row < elli_names.size(); ++row)
	{
		const std::string& elli_name = elli_names[row];
		int elli_row = m_elli_rows[row];

		int x = sett->value(("reso/ellipse_" + elli_name + "_x").c_str(), -2).toInt();
		if(x > -2)
			tabEllipses->item(elli_row, ELLI_X)->setText(tl::var_to_str(tl::clamp(x, -1, 3)).c_str());

		int y = sett->value(("reso/ellipse_" + elli_name + "_y").c_str(), -2).toInt();
		if(y > -2)
			tabEllipses->item(elli_row, ELLI_Y)->setText(tl::var_to_str(tl::clamp(y, -1, 3)).c_str());

		int proj1 = sett->value(("reso/ellipse_" + elli_name + "_proj1").c_str(), -2).toInt();
		if(proj1 > -2)
			tabEllipses->item(elli_row, ELLI_PROJ1)->setText(tl::var_to_str(tl::clamp(proj1, -1, 3)).c_str());

		int rem1 = sett->value(("reso/ellipse_" + elli_name + "_rem1").c_str(), -2).toInt();
		if(rem1 > -2)
			tabEllipses->item(elli_row, ELLI_REM1)->setText(tl::var_to_str(tl::clamp(rem1, -1, 3)).c_str());

		int proj2 = sett->value(("reso/ellipse_" + elli_name + "_proj2").c_str(), -2).toInt();
		if(proj2 > -2)
			tabEllipses->item(elli_row, ELLI_PROJ2)->setText(tl::var_to_str(tl::clamp(proj2, -1, 3)).c_str());

		int rem2 = sett->value(("reso/ellipse_" + elli_name + "_rem2").c_str(), -2).toInt();
		if(rem2 > -2)
			tabEllipses->item(elli_row, ELLI_REM2)->setText(tl::var_to_str(tl::clamp(rem2, -1, 3)).c_str());
	}


	// ellipsoid table
	std::vector<std::string> ello_names = { "1", "2"  };

	for(std::size_t row = 0; row < ello_names.size(); ++row)
	{
		const std::string& ello_name = ello_names[row];
		int ello_row = m_ello_rows[row];

		int x = sett->value(("reso/ellipsoid3d_" + ello_name + "_x").c_str(), -2).toInt();
		if(x > -2)
			tabEllipsoids->item(ello_row, ELLO_X)->setText(tl::var_to_str(tl::clamp(x, -1, 3)).c_str());

		int y = sett->value(("reso/ellipsoid3d_" + ello_name + "_y").c_str(), -2).toInt();
		if(y > -2)
			tabEllipsoids->item(ello_row, ELLO_Y)->setText(tl::var_to_str(tl::clamp(y, -1, 3)).c_str());

		int z = sett->value(("reso/ellipsoid3d_" + ello_name + "_z").c_str(), -2).toInt();
		if(z > -2)
			tabEllipsoids->item(ello_row, ELLO_Z)->setText(tl::var_to_str(tl::clamp(z, -1, 3)).c_str());

		int proj = sett->value(("reso/ellipsoid3d_" + ello_name + "_proj_or_rem").c_str(), -2).toInt();
		if(proj > -2)
			tabEllipsoids->item(ello_row, ELLO_PROJREM)->setText(tl::var_to_str(tl::clamp(proj, -1, 3)).c_str());
	}


	SetGlobals(sett);
	return true;
}



bool SettingsDlg::SaveSettings(QSettings *sett) const
{
	// use global settings if none given
	if(!sett)
		sett = m_pSettings;
	if(!sett)
		return false;


	for(const t_tupEdit& tup : m_vecEdits)
	{
		const std::string& strKey = std::get<0>(tup);
		QLineEdit* pEdit = std::get<2>(tup);

		sett->setValue(strKey.c_str(), pEdit->text());
	}

	for(const t_tupCheck& tup : m_vecChecks)
	{
		const std::string& strKey = std::get<0>(tup);
		QCheckBox* pCheck = std::get<2>(tup);

		sett->setValue(strKey.c_str(), pCheck->isChecked());
	}

	for(const t_tupSpin& tup : m_vecSpins)
	{
		const std::string& strKey = std::get<0>(tup);
		QSpinBox* pSpin = std::get<2>(tup);

		sett->setValue(strKey.c_str(), pSpin->value());
	}

	for(const t_tupCombo& tup : m_vecCombos)
	{
		const std::string& strKey = std::get<0>(tup);
		QComboBox* pCombo = std::get<2>(tup);

		// save both combo-box index and value
		sett->setValue(strKey.c_str(), pCombo->currentIndex());
		sett->setValue((strKey+"_value").c_str(), pCombo->currentText());
	}


	// ellipse table
	std::vector<std::string> elli_names = { "1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b" };

	for(std::size_t row = 0; row < elli_names.size(); ++row)
	{
		const std::string& elli_name = elli_names[row];
		int elli_row = m_elli_rows[row];

		sett->setValue(("reso/ellipse_" + elli_name + "_x").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipses->item(elli_row, ELLI_X)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipse_" + elli_name + "_y").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipses->item(elli_row, ELLI_Y)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipse_" + elli_name + "_proj1").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipses->item(elli_row, ELLI_PROJ1)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipse_" + elli_name + "_rem1").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipses->item(elli_row, ELLI_REM1)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipse_" + elli_name + "_proj2").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipses->item(elli_row, ELLI_PROJ2)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipse_" + elli_name + "_rem2").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipses->item(elli_row, ELLI_REM2)->text().toStdString()), -1, 3));
	}


	// ellipsoid table
	std::vector<std::string> ello_names = { "1", "2" };

	for(std::size_t row = 0; row < ello_names.size(); ++row)
	{
		const std::string& ello_name = ello_names[row];
		int ello_row = m_ello_rows[row];

		sett->setValue(("reso/ellipsoid3d_" + ello_name + "_x").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipsoids->item(ello_row, ELLO_X)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipsoid3d_" + ello_name + "_y").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipsoids->item(ello_row, ELLO_Y)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipsoid3d_" + ello_name + "_z").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipsoids->item(ello_row, ELLO_Z)->text().toStdString()), -1, 3));
		sett->setValue(("reso/ellipsoid3d_" + ello_name + "_proj_or_rem").c_str(),
			tl::clamp(tl::str_to_var<int>(tabEllipsoids->item(ello_row, ELLO_PROJREM)->text().toStdString()), -1, 3));
	}


	SetGlobals();
	return true;
}



void SettingsDlg::SetGlobals(QSettings *sett) const
{
	// use global settings if none given
	if(!sett)
		sett = m_pSettings;


	// precisions
	g_iPrec = spinPrecGen->value();
	g_iPrecGfx = spinPrecGfx->value();

	g_iMaxThreads = spinThreads->value();
	g_iMaxProcesses = spinProcesses->value();

	g_dEps = std::pow(10., -t_real(g_iPrec));
	g_dEpsGfx = std::pow(10., -t_real(g_iPrecGfx));

	GFX_NUM_POINTS = spinPtsGfx->value();
	g_iMaxNN = spinMaxNN->value();


	g_bShowFsq = (comboSFact->currentIndex() == 1);
	g_b3dBZ = (comboBZ->currentIndex() == 0);


	// fonts
	QString strGfxFont = editGfxFont->text();
	if(strGfxFont.length() != 0)
	{
		QFont font;
		if(font.fromString(strGfxFont))
		{
			g_fontGfx = font;

			g_dFontSize = g_fontGfx.pointSizeF();
			if(g_dFontSize <= 0.) g_dFontSize = 10.;
		}
	}

	if(editGLFont->text().length() != 0)
		g_strFontGL = editGLFont->text().toStdString();

	if(spinGLFont->value() > 0)
		g_iFontGLSize = spinGLFont->value();

	QString strGenFont = editGenFont->text();
	if(strGenFont.length() != 0)
	{
		QFont font;
		if(font.fromString(strGenFont))
			g_fontGen = font;
	}


	// external tools
	if(editGpl->text().length() != 0)
		g_strGplTool = find_program_binary(editGpl->text().toStdString());


	// if no GUI style has been set, use a safe default
	if(sett && !sett->contains("main/gui_style_value"))
	{
		int iGUIDefault = comboGUI->findText("Fusion", Qt::MatchContains);
		if(iGUIDefault >= 0)
			comboGUI->setCurrentIndex(iGUIDefault);
	}

	// GUI style
	QString strStyle = comboGUI->currentText();
	QStyle *pStyle = QStyleFactory::create(strStyle);
	if(strStyle!="" && pStyle)
		QApplication::setStyle(pStyle);
	else
		tl::log_err("Style \"", strStyle.toStdString(), "\" was not found.");

	emit SettingsChanged();
}



void SettingsDlg::SelectGLFont()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	// find a default font directory
	std::string strFontDir;
	std::vector<std::string> vecFontDir = get_qt_std_path(QtStdPath::FONTS);
	if(vecFontDir.size() == 0)
		tl::log_warn("Could not determine font directory.");
	else
		strFontDir = vecFontDir[0];

	// get font dir either from font file or from default
	std::string strPath;
	if(g_strFontGL != "")
		strPath = tl::get_dir(g_strFontGL);
	else
		strPath = strFontDir;

	QString strFile = QFileDialog::getOpenFileName(this,
		"Open Font File...", strPath.c_str(),
		"Font Files (*.ttf *.TTF)", nullptr, fileopt);
	if(strFile == "")
		return;

	editGLFont->setText(strFile);
}



void SettingsDlg::SelectGfxFont()
{
	bool bOk = false;
	QFont fontNew = QFontDialog::getFont(&bOk, g_fontGfx, this);
	if(bOk)
	{
		g_fontGfx = fontNew;
		g_dFontSize = g_fontGfx.pointSizeF();
		if(g_dFontSize <= 0.) g_dFontSize = 10.;

		editGfxFont->setText(fontNew.toString());
	}
}



void SettingsDlg::SelectGenFont()
{
	bool bOk = false;
	QFont fontNew = QFontDialog::getFont(&bOk, g_fontGen, this);
	if(bOk)
	{
		g_fontGen = fontNew;
		editGenFont->setText(fontNew.toString());
	}
}



void SettingsDlg::SelectGplTool()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strFile = QFileDialog::getOpenFileName(this,
		"Select Gnuplot Tool...", "",
		"Executable Files (* *.exe *.EXE)", nullptr, fileopt);
	if(strFile == "")
		return;

	editGpl->setText(strFile);
}



void SettingsDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}



void SettingsDlg::ButtonBoxClicked(QAbstractButton *pBtn)
{
	// load setting from file
	if(pBtn == static_cast<QAbstractButton*>(buttonBox->button(QDialogButtonBox::Open)))
	{
		QFileDialog::Option fileopt = QFileDialog::Option(0);
		if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
			fileopt = QFileDialog::DontUseNativeDialog;

		QString strDirLast = "";
		if(m_pSettings)
			strDirLast = m_pSettings->value("settings/last_dir", "~").toString();
		QString strFile = QFileDialog::getOpenFileName(this,
			"Load Configuration...", strDirLast,
			"Settings Files (*.ini *.INI)", nullptr, fileopt);
		if(strFile == "")
			return;

		QSettings sett(strFile, QSettings::IniFormat, this);
		bool ok = LoadSettings(&sett);

		if(ok && m_pSettings)
		{
			std::string strDir = tl::get_dir(strFile.toStdString());
			m_pSettings->setValue("settings/last_dir", QString(strDir.c_str()));
		}

		return;
	}

	// save setting to file
	else if(pBtn == static_cast<QAbstractButton*>(buttonBox->button(QDialogButtonBox::Save)))
	{
		QFileDialog::Option fileopt = QFileDialog::Option(0);
		if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
			fileopt = QFileDialog::DontUseNativeDialog;

		QString strDirLast = "";
		if(m_pSettings)
			strDirLast = m_pSettings->value("settings/last_dir", "~").toString();
		QString strFile = QFileDialog::getSaveFileName(this,
			"Save Configuration", strDirLast,
			"Settings FIles (*.ini *.INI)", nullptr, fileopt);
		if(strFile == "")
			return;

		QSettings sett(strFile, QSettings::IniFormat, this);
		bool ok = SaveSettings(&sett);

		if(ok && m_pSettings)
		{
			std::string strDir = tl::get_dir(strFile.toStdString());
			m_pSettings->setValue("settings/last_dir", QString(strDir.c_str()));
		}

		return;
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
	   buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		SaveSettings();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ResetRole)
	{
		SetDefaults(true);
		LoadSettings();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		QDialog::reject();
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("settings/geo", saveGeometry());

		QDialog::accept();
	}
}


#include "moc_SettingsDlg.cpp"
