/**
 * TAS tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2014
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

// mingw hack
#if defined(__MINGW32__) || defined(__MINGW64__)
	#if !defined(__kernel_entry)
		#define __kernel_entry
	#endif
#endif

#include "taz.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <QApplication>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>
#include <QFileDialog>
#include <QScrollBar>
#include <QMessageBox>
#include <QUrl>
#include <QDesktopServices>

#include "tlibs/phys/lattice.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/flags.h"
#include "tlibs/helper/proc.h"
#include "tlibs/log/log.h"
#include "tlibs/file/prop.h"

#include "libs/qt/recent.h"

namespace algo = boost::algorithm;
namespace fs = boost::filesystem;

using t_real = t_real_glob;
const std::string TazDlg::s_strTitle = "Takin";
const t_real_glob TazDlg::s_dPlaneDistTolerance = std::cbrt(tl::get_epsilon<t_real_glob>());


TazDlg::TazDlg(QWidget* pParent, const std::string& strLogFile)
	: QMainWindow(pParent), m_settings("takin", "core"),
		m_pSettingsDlg(new SettingsDlg(this, &m_settings)),
		m_strLogFile(strLogFile),
		m_pStatusMsg(new QLabel(this)),
		m_pCoordQStatusMsg(new QLabel(this)),
		m_pCoordCursorStatusMsg(new QLabel(this)),
		m_sceneRecip(this),
		m_dlgRecipParam(this, &m_settings),
		m_dlgRealParam(this, &m_settings),
		m_pGotoDlg(new GotoDlg(this, &m_settings))
{
	const bool bSmallqVisible = 0;
	const bool bCoordAxesVisible = 1;
	const bool bBZVisible = 1;
	const bool bWSVisible = 1;
	const bool bAllPeaksVisible = 1;

	this->setupUi(this);
	this->setWindowTitle(s_strTitle.c_str());
	this->setFont(g_fontGen);
	this->setWindowIcon(load_icon("res/icons/takin.svg"));

	btnAtoms->setEnabled(g_bHasScatlens);

	if(m_settings.contains("main/geo"))
	{
		QByteArray geo = m_settings.value("main/geo").toByteArray();
		restoreGeometry(geo);
	}

	m_pStatusMsg->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	m_pCoordQStatusMsg->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	m_pCoordCursorStatusMsg->setFrameStyle(QFrame::Panel | QFrame::Sunken);

	// prevents window resizing by message length
	for(QLabel* pLabel : {m_pStatusMsg/*, m_pCoordQStatusMsg, m_pCoordCursorStatusMsg*/})
		pLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);

	QStatusBar *pStatusBar = new QStatusBar(this);
	pStatusBar->addWidget(m_pStatusMsg, 1);
	pStatusBar->addPermanentWidget(m_pCoordQStatusMsg, 0);
	pStatusBar->addPermanentWidget(m_pCoordCursorStatusMsg, 0);
	m_pCoordQStatusMsg->setMinimumWidth(350);
	m_pCoordQStatusMsg->setAlignment(Qt::AlignCenter);
	m_pCoordCursorStatusMsg->setMinimumWidth(325);
	m_pCoordCursorStatusMsg->setAlignment(Qt::AlignCenter);
	this->setStatusBar(pStatusBar);

	// --------------------------------------------------------------------------------

	m_vecEdits_real =
	{
		editA, editB, editC,
		editAlpha, editBeta, editGamma
	};
	m_vecEdits_recip =
	{
		editARecip, editBRecip, editCRecip,
		editAlphaRecip, editBetaRecip, editGammaRecip
	};

	m_vecEdits_plane =
	{
		editScatX0, editScatX1, editScatX2,
		editScatY0, editScatY1, editScatY2,

		editRealX0, editRealX1, editRealX2,
		editRealY0, editRealY1, editRealY2,
	};
	m_vecEdits_monoana =
	{
		editMonoD, editAnaD
	};

	m_vecCheckBoxesSenses = { checkSenseM, checkSenseS, checkSenseA };


	m_vecEditNames_real =
	{
		"sample/a", "sample/b", "sample/c",
		"sample/alpha", "sample/beta", "sample/gamma"
	};
	m_vecEditNames_recip =
	{
		"sample/a_recip", "sample/b_recip", "sample/c_recip",
		"sample/alpha_recip", "sample/beta_recip", "sample/gamma_recip"
	};
	m_vecEditNames_plane =
	{
		"plane/x0", "plane/x1", "plane/x2",
		"plane/y0", "plane/y1", "plane/y2",

		"plane/real_x0", "plane/real_x1", "plane/real_x2",
		"plane/real_y0", "plane/real_y1", "plane/real_y2",
	};
	m_vecEditNames_monoana =
	{
		"tas/mono_d", "tas/ana_d"
	};

	m_vecSpinBoxNamesSample = {"sample/phi", "sample/theta", "sample/psi"};
	m_vecCheckBoxNamesSenses = {"tas/sense_m", "tas/sense_s", "tas/sense_a"};


	// recip
	m_pviewRecip = new ScatteringTriangleView(groupRecip);
	groupRecip->addTab(m_pviewRecip, "Reciprocal Lattice");

	m_pviewProjRecip = new ProjLatticeView(groupRecip);
	groupRecip->addTab(m_pviewProjRecip, "Projection");

	m_pviewRecip->setScene(&m_sceneRecip);
	m_pviewProjRecip->setScene(&m_sceneProjRecip);

	if(m_settings.contains("main/recip_tab"))
		groupRecip->setCurrentIndex(m_settings.value("main/recip_tab").value<int>());
	else
		groupRecip->setCurrentIndex(0);

	// real
	m_pviewRealLattice = new LatticeView(groupReal);
	groupReal->addTab(m_pviewRealLattice, "Real Lattice");

	m_pviewReal = new TasLayoutView(groupReal);
	groupReal->addTab(m_pviewReal, "TAS Instrument");

	m_pviewTof = new TofLayoutView(groupReal);
	groupReal->addTab(m_pviewTof, "TOF Instrument");

	m_pviewRealLattice->setScene(&m_sceneRealLattice);
	m_pviewReal->setScene(&m_sceneReal);
	m_pviewTof->setScene(&m_sceneTof);

	if(m_settings.contains("main/real_tab"))
		groupReal->setCurrentIndex(m_settings.value("main/real_tab").value<int>());
	else
		groupReal->setCurrentIndex(1);


	QObject::connect(m_pSettingsDlg, &SettingsDlg::SettingsChanged, this, &TazDlg::SettingsChanged);

	QObject::connect(&m_sceneReal, &TasLayoutScene::nodeEvent, this, &TazDlg::RealNodeEvent);
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::nodeEvent, this, &TazDlg::RecipNodeEvent);
	QObject::connect(&m_sceneTof, &TofLayoutScene::nodeEvent, this, &TazDlg::TofNodeEvent);

	// TAS
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::triangleChanged, &m_sceneReal, &TasLayoutScene::triangleChanged);
	QObject::connect(&m_sceneReal, &TasLayoutScene::tasChanged, &m_sceneRecip, &ScatteringTriangleScene::tasChanged);
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged, &m_sceneReal, &TasLayoutScene::recipParamsChanged);

	// TOF
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::triangleChanged, &m_sceneTof, &TofLayoutScene::triangleChanged);
	QObject::connect(&m_sceneTof, &TofLayoutScene::tasChanged, &m_sceneRecip, &ScatteringTriangleScene::tasChanged);
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged, &m_sceneTof, &TofLayoutScene::recipParamsChanged);

	// connections between instruments
	QObject::connect(&m_sceneTof, &TofLayoutScene::tasChanged, &m_sceneReal, &TasLayoutScene::triangleChanged);
	QObject::connect(&m_sceneReal, &TasLayoutScene::tasChanged, &m_sceneTof, &TofLayoutScene::triangleChanged);

	// scale factor
	if(m_pviewRecip)
		QObject::connect(m_pviewRecip, &ScatteringTriangleView::scaleChanged, &m_sceneRecip, &ScatteringTriangleScene::scaleChanged);
	if(m_pviewProjRecip)
		QObject::connect(m_pviewProjRecip, &ProjLatticeView::scaleChanged, &m_sceneProjRecip, &ProjLatticeScene::scaleChanged);
	if(m_pviewRealLattice)
		QObject::connect(m_pviewRealLattice, &LatticeView::scaleChanged, &m_sceneRealLattice, &LatticeScene::scaleChanged);
	if(m_pviewReal)
		QObject::connect(m_pviewReal, &TasLayoutView::scaleChanged, &m_sceneReal, &TasLayoutScene::scaleChanged);
	if(m_pviewTof)
		QObject::connect(m_pviewTof, &TofLayoutView::scaleChanged, &m_sceneTof, &TofLayoutScene::scaleChanged);

	// parameter dialogs
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged, &m_dlgRecipParam, &RecipParamDlg::paramsChanged);
	QObject::connect(&m_sceneReal, &TasLayoutScene::paramsChanged, &m_dlgRealParam, &RealParamDlg::paramsChanged);

	// cursor position
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::coordsChanged, this, &TazDlg::RecipCoordsChanged);
	QObject::connect(&m_sceneProjRecip, &ProjLatticeScene::coordsChanged, this, &TazDlg::RecipCoordsChanged);
	QObject::connect(&m_sceneRealLattice, &LatticeScene::coordsChanged, this, &TazDlg::RealCoordsChanged);

	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::spurionInfo, this, &TazDlg::spurionInfo);

	QObject::connect(m_pGotoDlg, &GotoDlg::vars_changed, this, &TazDlg::VarsChanged);
	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged, m_pGotoDlg, &GotoDlg::RecipParamsChanged);

	QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged, this, &TazDlg::recipParamsChanged);


	for(QLineEdit* pEdit : m_vecEdits_monoana)
		QObject::connect(pEdit, &QLineEdit::textEdited, this, &TazDlg::UpdateDs);

	for(QLineEdit* pEdit : m_vecEdits_real)
	{
		QObject::connect(pEdit, &QLineEdit::textEdited, this, &TazDlg::CheckCrystalType);
		QObject::connect(pEdit, &QLineEdit::textEdited, this, &TazDlg::CalcPeaks);
	}

	for(QLineEdit* pEdit : m_vecEdits_plane)
		QObject::connect(pEdit, &QLineEdit::textEdited, this, &TazDlg::CalcPeaks);

	for(QLineEdit* pEdit : m_vecEdits_recip)
	{
		QObject::connect(pEdit, &QLineEdit::textEdited, this, &TazDlg::CheckCrystalType);
		QObject::connect(pEdit, &QLineEdit::textEdited, this, &TazDlg::CalcPeaksRecip);
	}

	QObject::connect(checkSenseM, &QCheckBox::stateChanged, this, &TazDlg::UpdateMonoSense);
	QObject::connect(checkSenseS, &QCheckBox::stateChanged, this, &TazDlg::UpdateSampleSense);
	QObject::connect(checkSenseA, &QCheckBox::stateChanged, this, &TazDlg::UpdateAnaSense);

	QObject::connect(editSpaceGroupsFilter, &QLineEdit::textEdited, this, &TazDlg::RepopulateSpaceGroups);
	QObject::connect(comboSpaceGroups, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &TazDlg::SetCrystalType);
	QObject::connect(comboSpaceGroups, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &TazDlg::CalcPeaks);
	QObject::connect(checkPowder, &QCheckBox::stateChanged, this, &TazDlg::CalcPeaks);

	QObject::connect(btnAtoms, &QPushButton::clicked, this, &TazDlg::ShowAtomsDlg);


	// --------------------------------------------------------------------------------
	// file menu
	QMenu *pMenuFile = new QMenu("File", this);

	QAction *pNew = new QAction("New", this);
	pNew->setIcon(load_icon("res/icons/document-new.svg"));
	pNew->setShortcut(QKeySequence::New);
	pMenuFile->addAction(pNew);

	pMenuFile->addSeparator();

	QAction *pLoad = new QAction("Open...", this);
	pLoad->setIcon(load_icon("res/icons/document-open.svg"));
	pLoad->setShortcut(QKeySequence::Open);
	pMenuFile->addAction(pLoad);

	m_pMenuRecent = new QMenu("Recently Opened", this);
	RecentFiles recent(&m_settings, "main/recent");
	recent.FillMenu(m_pMenuRecent, [this](const std::string& str){ LoadFile(str.c_str()); });
	pMenuFile->addMenu(m_pMenuRecent);

	QAction *pSave = new QAction("Save", this);
	pSave->setIcon(load_icon("res/icons/document-save.svg"));
	pSave->setShortcut(QKeySequence::Save);
	pMenuFile->addAction(pSave);

	QAction *pSaveAs = new QAction("Save as...", this);
	pSaveAs->setIcon(load_icon("res/icons/document-save-as.svg"));
	pSaveAs->setShortcut(QKeySequence::SaveAs);
	pMenuFile->addAction(pSaveAs);

	pMenuFile->addSeparator();

	QAction *pImport = new QAction("Import Data File...", this);
	pImport->setIcon(load_icon("res/icons/drive-harddisk.svg"));
	pMenuFile->addAction(pImport);

	m_pMenuRecentImport = new QMenu("Recently Imported Data Files", this);
	RecentFiles recentimport(&m_settings, "main/recent_import");
	recentimport.FillMenu(m_pMenuRecentImport, [this](const std::string& str){ ImportFile(str.c_str()); });
	pMenuFile->addMenu(m_pMenuRecentImport);


	QAction *pImportCIF = new QAction("Import CIF...", this);
	pImportCIF->setIcon(load_icon("res/icons/drive-harddisk.svg"));
	pMenuFile->addAction(pImportCIF);

	m_pMenuRecentImportCIF = new QMenu("Recently Imported CIFs", this);
	RecentFiles recentimportCIF(&m_settings, "main/recent_import_cif");
	recentimportCIF.FillMenu(m_pMenuRecentImportCIF, [this](const std::string& str){ ImportCIFFile(str.c_str()); });
	pMenuFile->addMenu(m_pMenuRecentImportCIF);


	pMenuFile->addSeparator();

	QAction *pSettings = new QAction("Preferences...", this);
	pSettings->setMenuRole(QAction::PreferencesRole);
	pSettings->setIcon(load_icon("res/icons/preferences-system.svg"));
	pMenuFile->addAction(pSettings);

	pMenuFile->addSeparator();

	QAction *pExit = new QAction("Quit Takin", this);
	pExit->setMenuRole(QAction::QuitRole);
	pExit->setIcon(load_icon("res/icons/system-log-out.svg"));
	pExit->setShortcut(QKeySequence::Quit);
	pMenuFile->addAction(pExit);


	// --------------------------------------------------------------------------------
	// recip menu
	m_pMenuViewRecip = new QMenu("Reciprocal Space", this);

	m_pGoto = new QAction("Go to Position...", this);
	m_pGoto->setIcon(load_icon("res/icons/goto.svg"));
	m_pMenuViewRecip->addAction(m_pGoto);

	QAction *pRecipParams = new QAction("Information...", this);
	m_pMenuViewRecip->addAction(pRecipParams);
	m_pMenuViewRecip->addSeparator();

	m_pCoordAxes = new QAction("Show Coordinate Axes", this);
	m_pCoordAxes->setCheckable(1);
	m_pCoordAxes->setChecked(bCoordAxesVisible);
	m_pMenuViewRecip->addAction(m_pCoordAxes);

	m_pSmallq = new QAction("Show Reduced Scattering Vector q", this);
	m_pSmallq->setIcon(load_icon("res/icons/q.svg"));
	m_pSmallq->setCheckable(1);
	m_pSmallq->setChecked(bSmallqVisible);
	m_pMenuViewRecip->addAction(m_pSmallq);

	m_pSnapSmallq = new QAction("Snap G to Bragg Peak", this);
	m_pSnapSmallq->setCheckable(1);
	m_pSnapSmallq->setChecked(m_sceneRecip.getSnapq());
	m_pMenuViewRecip->addAction(m_pSnapSmallq);

	QAction *pKeepAbsKiKf = new QAction("Keep |ki| and |kf| Fixed", this);
	pKeepAbsKiKf->setCheckable(1);
	pKeepAbsKiKf->setChecked(m_sceneRecip.getKeepAbsKiKf());
	m_pMenuViewRecip->addAction(pKeepAbsKiKf);

	m_pBZ = new QAction("Show First Brillouin Zone", this);
	m_pBZ->setIcon(load_icon("res/icons/brillouin.svg"));
	m_pBZ->setCheckable(1);
	m_pBZ->setChecked(bBZVisible);
	m_pMenuViewRecip->addAction(m_pBZ);

	m_pAllPeaks = new QAction("Show Forbidden Peaks", this);
	m_pAllPeaks->setCheckable(1);
	m_pAllPeaks->setChecked(bAllPeaksVisible);
	m_pMenuViewRecip->addAction(m_pAllPeaks);


	QMenu *pMenuEwald = new QMenu("Ewald Sphere", this);
	QActionGroup *pGroupEwald = new QActionGroup(this);
	m_pEwaldSphereNone = new QAction("Disabled", this);
	m_pEwaldSphereKi = new QAction("Around ki", this);
	m_pEwaldSphereKf = new QAction("Around kf", this);
	for(QAction* pAct : {m_pEwaldSphereNone, m_pEwaldSphereKi, m_pEwaldSphereKf})
	{
		pAct->setCheckable(1);
		pGroupEwald->addAction(pAct);
	}

	m_pEwaldSphereKf->setChecked(1);
	pMenuEwald->addActions(pGroupEwald->actions());
	m_pMenuViewRecip->addMenu(pMenuEwald);


	m_pMenuViewRecip->addSeparator();

	QMenu *pMenuProj = new QMenu("Projection", this);
	pMenuProj->setTearOffEnabled(1);
	QActionGroup *pGroupProj = new QActionGroup(this);

	m_pProjGnom = new QAction("Gnomonic", this);
	m_pProjStereo = new QAction("Stereographic", this);
	m_pProjPara = new QAction("Parallel", this);
	m_pProjPersp = new QAction("Perspectivic", this);

	for(QAction *pAct : {m_pProjGnom, m_pProjStereo, m_pProjPara, m_pProjPersp})
	{
		pAct->setCheckable(1);
		pGroupProj->addAction(pAct);
	}

	m_pProjStereo->setChecked(1);
	pMenuProj->addActions(pGroupProj->actions());

	m_pMenuViewRecip->addMenu(pMenuProj);
	m_pMenuViewRecip->addSeparator();

#if !defined NO_3D
	QAction *pView3DBZ = new QAction("3D Brillouin Zone...", this);
	pView3DBZ->setIcon(load_icon("res/icons/brillouin3d.svg"));
	m_pMenuViewRecip->addAction(pView3DBZ);

	QAction *pView3D = new QAction("3D Bragg Peaks...", this);
	//pView3D->setIcon(QIcon::fromTheme("applications-graphics"));
	m_pMenuViewRecip->addAction(pView3D);

	m_pMenuViewRecip->addSeparator();
#endif


	QAction *pRecipExport = new QAction("Export Lattice Graphics...", this);
	pRecipExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewRecip->addAction(pRecipExport);

	QAction *pProjExport = new QAction("Export Projection Graphics...", this);
	pProjExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewRecip->addAction(pProjExport);

	QAction *pBZ3DExport = new QAction("Export 3D Brillouin Zone Model...", this);
	pBZ3DExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewRecip->addAction(pBZ3DExport);

	QAction *pBZCutExport = new QAction("Export Brillouin Zone Cut...", this);
	pBZCutExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewRecip->addAction(pBZCutExport);


	// --------------------------------------------------------------------------------
	// real menu
	m_pMenuViewReal = new QMenu("Real Space", this);
	m_pMenuViewReal->addAction(m_pGoto);

	QAction *pRealParams = new QAction("Information...", this);
	m_pMenuViewReal->addAction(pRealParams);

	m_pMenuViewReal->addSeparator();

	m_pShowRealQDir = new QAction("Show Q Direction", this);
	m_pShowRealQDir->setCheckable(1);
	m_pShowRealQDir->setChecked(m_sceneReal.GetTasLayout()->GetRealQVisible());
	m_pMenuViewReal->addAction(m_pShowRealQDir);

	m_pWS = new QAction("Show Unit Cell", this);
	m_pWS->setIcon(load_icon("res/icons/brillouin.svg"));
	m_pWS->setCheckable(1);
	m_pWS->setChecked(bWSVisible);
	m_pMenuViewReal->addAction(m_pWS);

	m_pMenuViewReal->addSeparator();

	QAction *pDarkAngles = new QAction("Dark Angles...", this);
	m_pMenuViewReal->addAction(pDarkAngles);

	m_pMenuViewReal->addSeparator();

#if !defined NO_3D
	QAction *pView3DReal = new QAction("3D Unit Cell...", this);
	pView3DReal->setIcon(load_icon("res/icons/unitcell3d.svg"));
	m_pMenuViewReal->addAction(pView3DReal);

	m_pMenuViewReal->addSeparator();
#endif

	QAction *pRealLatticeExport = new QAction("Export Lattice Graphics...", this);
	pRealLatticeExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewReal->addAction(pRealLatticeExport);

	QAction *pRealExport = new QAction("Export TAS Layout...", this);
	pRealExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewReal->addAction(pRealExport);

	QAction *pTofExport = new QAction("Export TOF Layout...", this);
	pTofExport->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewReal->addAction(pTofExport);

	QAction *pExportUC = new QAction("Export 3D Unit Cell Model...", this);
	pExportUC->setIcon(load_icon("res/icons/image-x-generic.svg"));
	m_pMenuViewReal->addAction(pExportUC);


	// --------------------------------------------------------------------------------
	// resolution menu
	QMenu *pMenuReso = new QMenu("Resolution", this);

	QAction *pResoParams = new QAction("Parameters...", this);
	pResoParams->setIcon(load_icon("res/icons/accessories-calculator.svg"));
	pMenuReso->addAction(pResoParams);

	pMenuReso->addSeparator();

	QAction *pResoEllipses = new QAction("Ellipses...", this);
	pResoEllipses->setIcon(load_icon("res/icons/ellipses.svg"));
	pMenuReso->addAction(pResoEllipses);

#if !defined NO_3D
	QAction *pResoEllipses3D = new QAction("3D Ellipsoids...", this);
	pMenuReso->addAction(pResoEllipses3D);
#endif

	pMenuReso->addSeparator();

	QAction *pResoConv = new QAction("Convolution...", this);
	pResoConv->setIcon(load_icon("res/icons/utilities-system-monitor.svg"));
	pMenuReso->addAction(pResoConv);


	// --------------------------------------------------------------------------------
	// calc menu
	QMenu *pMenuCalc = new QMenu("Calculation", this);

	QAction *pNeutronProps = new QAction("Neutron Properties...", this);
	pNeutronProps->setIcon(load_icon("res/icons/x-office-spreadsheet-template.svg"));
	pMenuCalc->addAction(pNeutronProps);

	QAction *pCompProps = new QAction("Components...", this);
	//pCompProps->setIcon(load_icon("res/icons/x-office-spreadsheet-template.svg"));
	pMenuCalc->addAction(pCompProps);

	pMenuCalc->addSeparator();

	QAction *pFormfactor = nullptr;
	if(g_bHasFormfacts && g_bHasScatlens)
	{
		pFormfactor = new QAction("Elements && Form Factors...", this);
		pMenuCalc->addAction(pFormfactor);
	}

	QAction *pScatteringFactors = new QAction("Scattering Factors...", this);
	pMenuCalc->addAction(pScatteringFactors);

	pMenuCalc->addSeparator();

	QAction *pPowder = new QAction("Powder Lines...", this);
	pPowder->setIcon(load_icon("res/icons/weather-snow.svg"));
	pMenuCalc->addAction(pPowder);

	QAction *pPowderFit = new QAction("Powder Fitting...", this);
	pMenuCalc->addAction(pPowderFit);

	pMenuCalc->addSeparator();

	QAction *pDynPlane = new QAction("Kinematic Plane...", this);
	pMenuCalc->addAction(pDynPlane);

	QAction *pSpuri = new QAction("Spurious Scattering...", this);
	pMenuCalc->addAction(pSpuri);

	QAction *pElast = new QAction("Spurious Positions...", this);
	pMenuCalc->addAction(pElast);


#if !defined NO_NET
	// --------------------------------------------------------------------------------
	// network menu
	QMenu *pMenuNet = new QMenu("Network", this);

	QAction *pConn = new QAction("Connect to Instrument...", this);
	pConn->setIcon(load_icon("res/icons/network-transmit-receive.svg"));
	pMenuNet->addAction(pConn);

	QAction *pDisconn = new QAction("Disconnect", this);
	pDisconn->setIcon(load_icon("res/icons/network-offline.svg"));
	pMenuNet->addAction(pDisconn);

	pMenuNet->addSeparator();

	QAction *pNetScanMon = new QAction("Scan Monitor...", this);
	pMenuNet->addAction(pNetScanMon);

	QAction *pNetCache = new QAction("Network Cache...", this);
	pMenuNet->addAction(pNetCache);

	QAction *pNetRefresh = new QAction("Refresh", this);
	pNetRefresh->setIcon(load_icon("res/icons/view-refresh.svg"));
	pMenuNet->addSeparator();
	pMenuNet->addAction(pNetRefresh);
#endif


	// --------------------------------------------------------------------------------
	// tools menu
	QMenu *pMenuTools = new QMenu("Tools", this);

	// add menu entries for external tools
	std::string strTools = find_resource("res/conf/tools.xml");
	bool toolconfloaded = false;
	if(strTools != "")
	{
		tl::Prop<std::string> propTools;
		if(propTools.Load(strTools.c_str(), tl::PropType::XML))
		{
			// to avoid two separators in a row
			bool bJustAddedSeparator = false;

			// add all menu entries
			for(std::size_t entry = 0; true; ++entry)
			{
				std::ostringstream _xmlpath;
				_xmlpath << "tools/entry_" << entry;
				std::string xmlpath = _xmlpath.str();
				if(!propTools.Exists(xmlpath))
					break;

				std::string tooltype = propTools.Query<std::string>(xmlpath + "/type", "");

				if(tooltype == "separator" && !bJustAddedSeparator)
				{
					pMenuTools->addSeparator();
					bJustAddedSeparator = true;
				}
				// external tools
				else if(tooltype == "program")
				{
					std::string toolname = propTools.Query<std::string>(xmlpath + "/name", "");
					std::string toolprog = propTools.Query<std::string>(xmlpath + "/program", "");

					if(toolname == "")
					{
						tl::log_err("Invalid tool name");
						continue;
					}

					std::string toolbin = find_program_binary(toolprog);
					if(toolbin == "")
					{
						tl::log_warn("Tool binary \"", toolprog, "\" was not found.");
						continue;
					}

					// add menu item for external tool
					QAction *actionTool = new QAction(toolname.c_str(), this);
					pMenuTools->addAction(actionTool);
					bJustAddedSeparator = false;

					QObject::connect(actionTool, &QAction::triggered, [toolbin]()
					{
						// run exernal tool process
						tl::log_debug("Running process \"", toolbin, "\"...");

						tl::PipeProc<char> proc((toolbin + "&").c_str(), false);
						if(!proc.IsReady())
							tl::log_err("Process \"", toolbin, "\" could not be created.");
					});
				}
				// internal tools
				else if(tooltype == "program_internal")
				{
					std::string toolname = propTools.Query<std::string>(xmlpath + "/name", "");
					std::string toolprog = propTools.Query<std::string>(xmlpath + "/program", "");

					QAction *actionTool = new QAction(toolname.c_str(), this);

					if(toolprog == "takin_scanviewer")
					{
						QObject::connect(actionTool, &QAction::triggered, this, &TazDlg::ShowScanViewer);
					}
					else if(toolprog == "takin_scanpos")
					{
						QObject::connect(actionTool, &QAction::triggered, this, &TazDlg::ShowScanPos);
					}
					else if(toolprog == "takin_sgbrowser")
					{
						QObject::connect(actionTool, &QAction::triggered, this, &TazDlg::ShowSgListDlg);
					}
					else
					{
						tl::log_err("Unknown internal tool \"", toolprog, "\".");
						delete actionTool;
						actionTool = nullptr;
					}

					if(actionTool)
					{
						pMenuTools->addAction(actionTool);
						bJustAddedSeparator = false;
					}
				}
			}

			toolconfloaded = true;
		}
		else
		{
			tl::log_err("Cannot load tool configuration file \"", strTools, "\".");
		}
	}
	else
	{
		tl::log_err("No tool configuration file found.");
	}

	// default tools menu configuration, if nothing else given
	if(!toolconfloaded)
	{
		tl::log_warn("No tool configuration given, loading internal defaults.");

		QAction *pScanViewer = new QAction("Scan Viewer...", this);
		QObject::connect(pScanViewer, &QAction::triggered, this, &TazDlg::ShowScanViewer);
		pMenuTools->addAction(pScanViewer);

		QAction *pScanPos = new QAction("Scan Positions Plotter...", this);
		QObject::connect(pScanPos, &QAction::triggered, this, &TazDlg::ShowScanPos);
		pMenuTools->addAction(pScanPos);

		pMenuTools->addSeparator();

		QAction *pSgList = new QAction("Space Groups...", this);
		QObject::connect(pSgList, &QAction::triggered, this, &TazDlg::ShowSgListDlg);
		pMenuTools->addAction(pSgList);
	}


	// --------------------------------------------------------------------------------
	// help menu
	QMenu *pMenuHelp = new QMenu("Help", this);

	QAction *pHelp = new QAction("Show Help...", this);
	pHelp->setIcon(load_icon("res/icons/help-browser.svg"));
	pHelp->setShortcut(QKeySequence::HelpContents);
	pMenuHelp->addAction(pHelp);

	QAction *pDevelDoc = new QAction("Show Developer Help...", this);
	if(find_resource("doc/html/index.html", 0) == "")
		pDevelDoc->setEnabled(0);
	pDevelDoc->setIcon(load_icon("res/icons/help-browser.svg"));
	pMenuHelp->addAction(pDevelDoc);

	pMenuHelp->addSeparator();
	QAction *pLog = new QAction("Log...", this);
	pMenuHelp->addAction(pLog);

	pMenuHelp->addSeparator();
	QAction *pWebsite = new QAction("Show Takin Website...", this);
	pMenuHelp->addAction(pWebsite);
	QAction *pBugReport = new QAction("Report Bug...", this);
	pMenuHelp->addAction(pBugReport);

	pMenuHelp->addSeparator();

	QAction *pAboutQt = new QAction("About Qt...", this);
	pAboutQt->setMenuRole(QAction::AboutQtRole);
	pMenuHelp->addAction(pAboutQt);

	QAction *pExtLicenses = new QAction("Show 3rd Party Licenses...", this);
	pMenuHelp->addAction(pExtLicenses);

	pMenuHelp->addSeparator();

	QAction *pAbout = new QAction("About Takin...", this);
	pAbout->setMenuRole(QAction::AboutRole);
	pAbout->setIcon(load_icon("res/icons/dialog-information.svg"));
	pMenuHelp->addAction(pAbout);



	// --------------------------------------------------------------------------------
	QMenuBar *pMenuBar = new QMenuBar(this);
	pMenuBar->setNativeMenuBar(m_settings.value("main/native_dialogs", 1).toBool());

	pMenuBar->addMenu(pMenuFile);
	pMenuBar->addMenu(m_pMenuViewRecip);
	pMenuBar->addMenu(m_pMenuViewReal);
	pMenuBar->addMenu(pMenuReso);
	pMenuBar->addMenu(pMenuCalc);
	pMenuBar->addMenu(pMenuTools);
#if !defined NO_NET
	pMenuBar->addMenu(pMenuNet);
#endif
	pMenuBar->addMenu(pMenuHelp);


	QObject::connect(pNew, &QAction::triggered, this, &TazDlg::New);
	QObject::connect(pLoad, &QAction::triggered, this, static_cast<bool (TazDlg::*)()>(&TazDlg::Load));
	QObject::connect(pSave, &QAction::triggered, this, &TazDlg::Save);
	QObject::connect(pSaveAs, &QAction::triggered, this, &TazDlg::SaveAs);
	QObject::connect(pImport, &QAction::triggered, this, static_cast<bool (TazDlg::*)()>(&TazDlg::Import));
	QObject::connect(pImportCIF, &QAction::triggered, this, static_cast<bool (TazDlg::*)()>(&TazDlg::ImportCIF));

	QObject::connect(pPowderFit, &QAction::triggered, this, &TazDlg::ShowPowderFit);
	QObject::connect(pSettings, &QAction::triggered, this, &TazDlg::ShowSettingsDlg);
	QObject::connect(pExit, &QAction::triggered, this, &TazDlg::close);

	QObject::connect(m_pSmallq, &QAction::toggled, this, &TazDlg::EnableSmallq);
	QObject::connect(m_pCoordAxes, &QAction::toggled, this, &TazDlg::EnableCoordAxes);
	QObject::connect(m_pBZ, &QAction::toggled, this, &TazDlg::EnableBZ);
	QObject::connect(m_pWS, &QAction::toggled, this, &TazDlg::EnableWS);
	QObject::connect(m_pAllPeaks, &QAction::toggled, this, &TazDlg::ShowAllPeaks);

	for(QAction* pAct : {m_pEwaldSphereNone, m_pEwaldSphereKi, m_pEwaldSphereKf})
		QObject::connect(pAct, &QAction::triggered, this, &TazDlg::ShowEwaldSphere);

	QObject::connect(m_pShowRealQDir, &QAction::toggled, this, &TazDlg::EnableRealQDir);

	QObject::connect(m_pSnapSmallq, &QAction::toggled, &m_sceneRecip, &ScatteringTriangleScene::setSnapq);
	QObject::connect(pKeepAbsKiKf, &QAction::toggled, &m_sceneRecip, &ScatteringTriangleScene::setKeepAbsKiKf);

	QObject::connect(pRecipParams, &QAction::triggered, this, &TazDlg::ShowRecipParams);
	QObject::connect(pRealParams, &QAction::triggered, this, &TazDlg::ShowRealParams);

	for(QAction *pProj : {m_pProjGnom, m_pProjStereo, m_pProjPara, m_pProjPersp})
		QObject::connect(pProj, &QAction::triggered, this, &TazDlg::RecipProjChanged);

#if !defined NO_3D
	QObject::connect(pView3D, &QAction::triggered, this, &TazDlg::Show3D);
	QObject::connect(pView3DReal, &QAction::triggered, this, &TazDlg::Show3DReal);
	QObject::connect(pView3DBZ, &QAction::triggered, this, &TazDlg::Show3DBZ);
	QObject::connect(pResoEllipses3D, &QAction::triggered, this, &TazDlg::ShowResoEllipses3D);
#endif

	QObject::connect(pRecipExport, &QAction::triggered, this, &TazDlg::ExportRecip);
	QObject::connect(pBZ3DExport, &QAction::triggered, this, &TazDlg::ExportBZ3DModel);
	QObject::connect(pBZCutExport, &QAction::triggered, this, &TazDlg::ExportBZCut);
	QObject::connect(pProjExport, &QAction::triggered, this, &TazDlg::ExportProj);
	QObject::connect(pRealExport, &QAction::triggered, this, &TazDlg::ExportReal);
	QObject::connect(pTofExport, &QAction::triggered, this, &TazDlg::ExportTof);
	QObject::connect(pRealLatticeExport, &QAction::triggered, this, &TazDlg::ExportRealLattice);

	QObject::connect(pExportUC, &QAction::triggered, this, &TazDlg::ExportUCModel);

	QObject::connect(pResoParams, &QAction::triggered, this, &TazDlg::ShowResoParams);
	QObject::connect(pResoEllipses, &QAction::triggered, this, &TazDlg::ShowResoEllipses);
	QObject::connect(pResoConv, &QAction::triggered, this, &TazDlg::ShowResoConv);

	QObject::connect(pNeutronProps, &QAction::triggered, this, &TazDlg::ShowNeutronDlg);
	QObject::connect(pCompProps, &QAction::triggered, this, &TazDlg::ShowTofDlg);
	QObject::connect(m_pGoto, &QAction::triggered, this, &TazDlg::ShowGotoDlg);
	QObject::connect(pPowder, &QAction::triggered, this, &TazDlg::ShowPowderDlg);
	QObject::connect(pElast, &QAction::triggered, this, &TazDlg::ShowElasticDlg);
	QObject::connect(pSpuri, &QAction::triggered, this, &TazDlg::ShowSpurionDlg);
	QObject::connect(pScatteringFactors, &QAction::triggered, this, &TazDlg::ShowScatteringFactorsDlg);
	QObject::connect(pDynPlane, &QAction::triggered, this, &TazDlg::ShowDynPlaneDlg);
	QObject::connect(pDarkAngles, &QAction::triggered, this, &TazDlg::ShowDarkAnglesDlg);

#if !defined NO_NET
	QObject::connect(pConn, &QAction::triggered, this, &TazDlg::ShowConnectDlg);
	QObject::connect(pDisconn, &QAction::triggered, this, &TazDlg::Disconnect);
	QObject::connect(pNetRefresh, &QAction::triggered, this, &TazDlg::NetRefresh);
	QObject::connect(pNetCache, &QAction::triggered, this, &TazDlg::ShowNetCache);
	QObject::connect(pNetScanMon, &QAction::triggered, this, &TazDlg::ShowNetScanMonitor);
#endif

	if(pFormfactor)
		QObject::connect(pFormfactor, &QAction::triggered, this, &TazDlg::ShowFormfactorDlg);

	QObject::connect(pHelp, &QAction::triggered, this, &TazDlg::ShowHelp);
	QObject::connect(pDevelDoc, &QAction::triggered, this, &TazDlg::ShowDevelDoc);
	QObject::connect(pLog, &QAction::triggered, this, &TazDlg::ShowLog);
	QObject::connect(pWebsite, &QAction::triggered, this, &TazDlg::ShowWebsite);
	QObject::connect(pBugReport, &QAction::triggered, this, &TazDlg::ReportBug);
	QObject::connect(pAbout, &QAction::triggered, this, &TazDlg::ShowAbout);
	QObject::connect(pAboutQt, &QAction::triggered, qApp, &QApplication::aboutQt);
	QObject::connect(pExtLicenses, &QAction::triggered, this, &TazDlg::ShowExternalLicenses);


	setMenuBar(pMenuBar);
	// --------------------------------------------------------------------------------


	// --------------------------------------------------------------------------------
	// context menus
	if(m_pviewRecip) m_pviewRecip->setContextMenuPolicy(Qt::CustomContextMenu);
	if(m_pviewProjRecip) m_pviewProjRecip->setContextMenuPolicy(Qt::CustomContextMenu);
	if(m_pviewRealLattice) m_pviewRealLattice->setContextMenuPolicy(Qt::CustomContextMenu);
	if(m_pviewReal) m_pviewReal->setContextMenuPolicy(Qt::CustomContextMenu);
	if(m_pviewTof) m_pviewTof->setContextMenuPolicy(Qt::CustomContextMenu);

	if(m_pviewRecip)
		QObject::connect(m_pviewRecip, &ScatteringTriangleView::customContextMenuRequested, this, &TazDlg::RecipContextMenu);
	if(m_pviewProjRecip)
		QObject::connect(m_pviewProjRecip, &ProjLatticeView::customContextMenuRequested, this, &TazDlg::RecipContextMenu);
	if(m_pviewRealLattice)
		QObject::connect(m_pviewRealLattice, &LatticeView::customContextMenuRequested, this, &TazDlg::RealContextMenu);
	if(m_pviewReal)
		QObject::connect(m_pviewReal, &TasLayoutView::customContextMenuRequested, this, &TazDlg::RealContextMenu);
	if(m_pviewTof)
		QObject::connect(m_pviewTof, &TofLayoutView::customContextMenuRequested, this, &TazDlg::RealContextMenu);


	// --------------------------------------------------------------------------------
	// tool bars
	// file
	QToolBar *pFileTools = new QToolBar(this);
	pFileTools->setWindowTitle("File");
	pFileTools->addAction(pNew);
	pFileTools->addAction(pLoad);
	pFileTools->addAction(pImport);
	pFileTools->addAction(pSave);
	pFileTools->addAction(pSaveAs);
	addToolBar(pFileTools);

	// recip
	QToolBar *pRecipTools = new QToolBar(this);
	pRecipTools->setWindowTitle("Reciprocal Space");
	pRecipTools->addAction(m_pGoto);
	pRecipTools->addAction(m_pSmallq);
	pRecipTools->addAction(m_pBZ);
	//pRecipTools->addAction(m_pEwaldSphere);
	addToolBar(pRecipTools);

	// reso
	QToolBar *pResoTools = new QToolBar(this);
	pResoTools->setWindowTitle("Resolution");
	pResoTools->addAction(pResoParams);
	pResoTools->addAction(pResoEllipses);
	pResoTools->addAction(pResoConv);
	addToolBar(pResoTools);

	// calc
	QToolBar *pCalcTools = new QToolBar(this);
	pCalcTools->setWindowTitle("Calculation");
	pCalcTools->addAction(pNeutronProps);
	pCalcTools->addAction(pPowder);
	addToolBar(pCalcTools);

#if !defined NO_3D
	// 3d tools
	QToolBar *p3dTools = new QToolBar(this);
	p3dTools->setWindowTitle("3D Tools");
	p3dTools->addAction(pView3DBZ);
	p3dTools->addAction(pView3DReal);
	addToolBar(p3dTools);
#endif

#if !defined NO_NET
	// net
	QToolBar *pNetTools = new QToolBar(this);
	pNetTools->setWindowTitle("Network");
	pNetTools->addAction(pConn);
	pNetTools->addAction(pDisconn);
	pNetTools->addAction(pNetRefresh);
	addToolBar(pNetTools);
#endif

	// --------------------------------------------------------------------------------


	RepopulateSpaceGroups();

	unsigned int iMaxPeaks = m_settings.value("main/max_peaks", 10).toUInt();

	m_sceneRecip.GetTriangle()->SetMaxPeaks(iMaxPeaks);
	m_sceneRecip.GetTriangle()->SetPlaneDistTolerance(s_dPlaneDistTolerance);

	m_sceneProjRecip.GetLattice()->SetMaxPeaks(iMaxPeaks/2);

	m_sceneRealLattice.GetLattice()->SetMaxPeaks(iMaxPeaks);
	m_sceneRealLattice.GetLattice()->SetPlaneDistTolerance(s_dPlaneDistTolerance);

#if !defined NO_3D
	if(m_pRecip3d)
	{
		m_pRecip3d->SetMaxPeaks(t_real(iMaxPeaks/2));
		m_pRecip3d->SetPlaneDistTolerance(s_dPlaneDistTolerance);
	}
#endif

	m_bReady = true;
	UpdateDs();
	CalcPeaks();


	m_sceneRecip.GetTriangle()->SetqVisible(bSmallqVisible);
	m_sceneRecip.GetTriangle()->SetCoordAxesVisible(bCoordAxesVisible);
	m_sceneRecip.GetTriangle()->SetBZVisible(bBZVisible);
	m_sceneRecip.GetTriangle()->SetAllPeaksVisible(bAllPeaksVisible);
	m_sceneRecip.GetTriangle()->SetEwaldSphereVisible(EWALD_KF);
	m_sceneRealLattice.GetLattice()->SetWSVisible(bWSVisible);
	m_sceneRecip.emitUpdate();
	//m_sceneRecip.emitAllParams();

	setAcceptDrops(true);
}


TazDlg::~TazDlg()
{
	Disconnect();
	DeleteDialogs();

	// don't delete non-optional sub-modules in DeleteDialogs()
	if(m_pGotoDlg) { delete m_pGotoDlg; m_pGotoDlg = nullptr; }
	if(m_pSettingsDlg) { delete m_pSettingsDlg; m_pSettingsDlg = nullptr; }

	if(m_pviewRecip) { delete m_pviewRecip; m_pviewRecip = nullptr; }
	if(m_pviewProjRecip) { delete m_pviewProjRecip; m_pviewProjRecip = nullptr; }
	if(m_pviewRealLattice) { delete m_pviewRealLattice; m_pviewRealLattice = nullptr; }
	if(m_pviewReal) { delete m_pviewReal; m_pviewReal = nullptr; }
	if(m_pviewTof) { delete m_pviewTof; m_pviewTof = nullptr; }

	comboSpaceGroups->clear();
}


void TazDlg::DeleteDialogs()
{
	if(m_pAboutDlg) { delete m_pAboutDlg; m_pAboutDlg = nullptr; }
	if(m_pLogDlg) { delete m_pLogDlg; m_pLogDlg = nullptr; }
	if(m_pEllipseDlg) { delete m_pEllipseDlg; m_pEllipseDlg = nullptr; }
	if(m_pReso) { delete m_pReso; m_pReso = nullptr; }
	if(m_pConvoDlg) { delete m_pConvoDlg; m_pConvoDlg = nullptr; }
	if(m_pSpuri) { delete m_pSpuri; m_pSpuri = nullptr; }
	if(m_pNeutronDlg) { delete m_pNeutronDlg; m_pNeutronDlg = nullptr; }
	if(m_pTofDlg) { delete m_pTofDlg; m_pTofDlg = nullptr; }
	if(m_pScatteringFactorsDlg) { delete m_pScatteringFactorsDlg; m_pScatteringFactorsDlg = nullptr; }
	if(m_pDynPlaneDlg) { delete m_pDynPlaneDlg; m_pDynPlaneDlg = nullptr; }
	if(m_pScanViewer) { delete m_pScanViewer; m_pScanViewer = nullptr; }
	if(m_pScanPos) { delete m_pScanPos; m_pScanPos = nullptr; }
	if(m_pPowderFit) { delete m_pPowderFit; m_pPowderFit = nullptr; }
	if(m_pAtomsDlg) { delete m_pAtomsDlg; m_pAtomsDlg = nullptr; }
	if(m_pDarkAnglesDlg) { delete m_pDarkAnglesDlg; m_pDarkAnglesDlg = nullptr; }
	if(m_pPowderDlg) { delete m_pPowderDlg; m_pPowderDlg = nullptr; }
	if(m_pElasticDlg) { delete m_pElasticDlg; m_pElasticDlg = nullptr; }

#if !defined NO_3D
	if(m_pRecip3d) { delete m_pRecip3d; m_pRecip3d = nullptr; }
	if(m_pReal3d) { delete m_pReal3d; m_pReal3d = nullptr; }
	if(m_pBZ3d) { delete m_pBZ3d; m_pBZ3d = nullptr; }
	if(m_pEllipseDlg3D) { delete m_pEllipseDlg3D; m_pEllipseDlg3D = nullptr; }
#endif

#if !defined NO_NET
	if(m_pSrvDlg) { delete m_pSrvDlg; m_pSrvDlg = nullptr; }
	if(m_pScanMonDlg) { delete m_pScanMonDlg; m_pScanMonDlg = nullptr; }
	if(m_pNetCacheDlg) { delete m_pNetCacheDlg; m_pNetCacheDlg = nullptr; }
	if(m_pNetCache) { delete m_pNetCache; m_pNetCache = nullptr; }
#endif

	if(m_pSgListDlg) { delete m_pSgListDlg; m_pSgListDlg = nullptr; }
	if(m_pFormfactorDlg) { delete m_pFormfactorDlg; m_pFormfactorDlg = nullptr; }
}


void TazDlg::SettingsChanged()
{
	setFont(g_fontGen);

	m_sceneReal.update();
	m_sceneTof.update();
	m_sceneRealLattice.update();
	m_sceneRecip.update();
}


void TazDlg::keyPressEvent(QKeyEvent *pEvt)
{
	// x rotation
	if(pEvt->key() == Qt::Key_8)
		RotatePlane(0, tl::d2r<t_real>(-5.));
	else if(pEvt->key() == Qt::Key_2)
		RotatePlane(0, tl::d2r<t_real>(5.));

	// y rotation
	else if(pEvt->key() == Qt::Key_4)
		RotatePlane(1, tl::d2r<t_real>(-5.));
	else if(pEvt->key() == Qt::Key_6)
		RotatePlane(1, tl::d2r<t_real>(5.));

	// z rotation
	else if(pEvt->key() == Qt::Key_9)
		RotatePlane(2, tl::d2r<t_real>(-5.));
	else if(pEvt->key() == Qt::Key_7)
		RotatePlane(2, tl::d2r<t_real>(5.));

	QMainWindow::keyPressEvent(pEvt);
}


void TazDlg::keyReleaseEvent(QKeyEvent *pEvt)
{
	QMainWindow::keyReleaseEvent(pEvt);
}


void TazDlg::showEvent(QShowEvent *pEvt)
{
	QMainWindow::showEvent(pEvt);

	static bool bInitialShow = true;
	if(bInitialShow)
	{
		bInitialShow = false;

		if(m_pviewRecip) m_pviewRecip->centerOn(m_sceneRecip.GetTriangle()->GetGfxMid());
		if(m_pviewProjRecip) m_pviewProjRecip->centerOn(0.,0.);
		if(m_pviewReal) m_pviewReal->centerOn(m_sceneReal.GetTasLayout());
		if(m_pviewTof) m_pviewTof->centerOn(m_sceneTof.GetTofLayout());
		if(m_pviewRealLattice) m_pviewRealLattice->centerOn(0.,0.);
	}
}


void TazDlg::dragEnterEvent(QDragEnterEvent *pEvt)
{
	if(pEvt) pEvt->accept();
}


void TazDlg::dropEvent(QDropEvent *pEvt)
{
	if(!pEvt)
		return;
	const QMimeData* pMime = pEvt->mimeData();
	if(!pMime)
		return;

	std::string strFiles = pMime->text().toStdString();
	std::vector<std::string> vecFiles;
	tl::get_tokens<std::string, std::string>(strFiles, "\n", vecFiles);
	if(vecFiles.size() > 1)
		tl::log_warn("More than one file dropped, using first one.");

	if(vecFiles.size() >= 1)
	{
		std::string& strFile = vecFiles[0];
		tl::trim(strFile);

		const std::string strHead = "file://";
		if(algo::starts_with(strFile, strHead))
			algo::replace_head(strFile, strHead.length(), "");

		Load(strFile.c_str());
	}
}


void TazDlg::closeEvent(QCloseEvent* pEvt)
{
	m_bReady = false;

	m_settings.setValue("main/geo", saveGeometry());
	m_settings.setValue("main/recip_tab", groupRecip->currentIndex());
	m_settings.setValue("main/real_tab", groupReal->currentIndex());

	QMainWindow::closeEvent(pEvt);
}


void TazDlg::UpdateDs()
{
	t_real dMonoD = editMonoD->text().toDouble();
	t_real dAnaD = editAnaD->text().toDouble();

	m_sceneRecip.SetDs(dMonoD, dAnaD);

	ResoParams resoparams;
	resoparams.bMonoDChanged = true;
	resoparams.bAnaDChanged = true;
	resoparams.dMonoD = dMonoD;
	resoparams.dAnaD = dAnaD;

	if(m_pGotoDlg)
	{
		m_pGotoDlg->SetD(editMonoD->text().toDouble(), editAnaD->text().toDouble());
		m_pGotoDlg->CalcMonoAna();
		m_pGotoDlg->CalcSample();
	}

	if(m_pElasticDlg)
	{
		m_pElasticDlg->SetD(editMonoD->text().toDouble(), editAnaD->text().toDouble());
		m_pElasticDlg->CalcSpuriousPositions();
	}

	emit ResoParamsChanged(resoparams);
}


void TazDlg::UpdateSampleSense()
{
	const bool bSense = checkSenseS->isChecked();
	m_sceneReal.SetSampleSense(bSense);
	m_sceneTof.SetSampleSense(bSense);
	m_sceneRecip.SetSampleSense(bSense);

	if(m_pGotoDlg)
	{
		m_pGotoDlg->SetSampleSense(bSense);
		m_pGotoDlg->CalcSample();
	}

	if(m_pElasticDlg)
	{
		m_pElasticDlg->SetSampleSense(bSense);
		m_pElasticDlg->CalcSpuriousPositions();
	}

	m_dlgRealParam.SetSampleSense(bSense);

	ResoParams params;
	params.bSensesChanged[1] = true;
	params.bScatterSenses[1] = bSense;
	emit ResoParamsChanged(params);

	m_sceneRecip.emitUpdate();
}


void TazDlg::UpdateMonoSense()
{
	const bool bSense = checkSenseM->isChecked();
	m_sceneRecip.SetMonoSense(bSense);

	if(m_pGotoDlg)
	{
		m_pGotoDlg->SetMonoSense(bSense);
		m_pGotoDlg->CalcMonoAna();
	}

	if(m_pElasticDlg)
	{
		m_pElasticDlg->SetMonoSense(bSense);
		m_pElasticDlg->CalcSpuriousPositions();
	}

	ResoParams params;
	params.bSensesChanged[0] = true;
	params.bScatterSenses[0] = bSense;
	emit ResoParamsChanged(params);
}


void TazDlg::UpdateAnaSense()
{
	const bool bSense = checkSenseA->isChecked();
	m_sceneRecip.SetAnaSense(bSense);

	if(m_pGotoDlg)
	{
		m_pGotoDlg->SetAnaSense(bSense);
		m_pGotoDlg->CalcMonoAna();
	}

	if(m_pElasticDlg)
	{
		m_pElasticDlg->SetAnaSense(bSense);
		m_pElasticDlg->CalcSpuriousPositions();
	}

	ResoParams params;
	params.bSensesChanged[2] = true;
	params.bScatterSenses[2] = bSense;
	emit ResoParamsChanged(params);
}


void TazDlg::RecipNodeEvent(bool bStarted)
{
	// optimises reso dialog update policy
	if(m_pReso)
		m_pReso->SetUpdateOn(!bStarted, 1);
}


void TazDlg::RealNodeEvent(bool bStarted)
{
	// optimises reso dialog update policy
	if(m_pReso)
		m_pReso->SetUpdateOn(1, !bStarted);
}


void TazDlg::TofNodeEvent(bool bStarted)
{
	// optimises reso dialog update policy
	if(m_pReso)
		m_pReso->SetUpdateOn(1, !bStarted);
}


void TazDlg::RecipProjChanged()
{
	LatticeProj proj = LatticeProj::PARALLEL;
	if(m_pProjGnom->isChecked())
		proj = LatticeProj::GNOMONIC;
	else if(m_pProjStereo->isChecked())
		proj = LatticeProj::STEREOGRAPHIC;
	else if(m_pProjPara->isChecked())
		proj = LatticeProj::PARALLEL;
	else if(m_pProjPersp->isChecked())
		proj = LatticeProj::PERSPECTIVE;

	if(m_sceneProjRecip.GetLattice())
	{
		m_sceneProjRecip.GetLattice()->SetProjection(proj);
		m_sceneProjRecip.GetLattice()->CalcPeaks(m_latticecommon, true);
	}
}


// 3d stuff
#if !defined NO_3D

void TazDlg::Show3D()
{
	if(!m_pRecip3d)
	{
		m_pRecip3d = new Recip3DDlg(this, &m_settings);

		t_real dTol = s_dPlaneDistTolerance;
		m_pRecip3d->SetPlaneDistTolerance(dTol);

		// also track current Q position
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged,
			m_pRecip3d, &Recip3DDlg::RecipParamsChanged);
		m_sceneRecip.emitAllParams();
	}

	if(!m_pRecip3d->isVisible())
		m_pRecip3d->show();
	m_pRecip3d->activateWindow();

	m_pRecip3d->CalcPeaks(m_latticecommon);
	//CalcPeaks();
}


void TazDlg::Show3DReal()
{
	if(!m_pReal3d)
		m_pReal3d = new Real3DDlg(this, &m_settings);

	if(!m_pReal3d->isVisible())
		m_pReal3d->show();
	m_pReal3d->activateWindow();

	m_pReal3d->CalcPeaks(m_sceneRealLattice.GetLattice()->GetWS3D(), m_latticecommon);
	//CalcPeaks();
}


void TazDlg::Show3DBZ()
{
	if(!m_pBZ3d)
	{
		m_pBZ3d = new BZ3DDlg(this, &m_settings);

		// also track current q position
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged,
			m_pBZ3d, &BZ3DDlg::RecipParamsChanged);
		m_sceneRecip.emitAllParams();
	}

	if(!m_pBZ3d->isVisible())
		m_pBZ3d->show();
	m_pBZ3d->activateWindow();

	m_pBZ3d->RenderBZ(m_sceneRecip.GetTriangle()->GetBZ3D(),
		m_latticecommon,
		&m_sceneRecip.GetTriangle()->GetBZ3DPlaneVerts(),
		&m_sceneRecip.GetTriangle()->GetBZ3DSymmVerts());
	//CalcPeaks();
}


#else

void TazDlg::Show3D() {}
void TazDlg::Show3DReal() {}
void TazDlg::Show3DBZ() {}

#endif


void TazDlg::EnableSmallq(bool bEnable)
{
	m_sceneRecip.GetTriangle()->SetqVisible(bEnable);
}


void TazDlg::EnableCoordAxes(bool bEnable)
{
	m_sceneRecip.GetTriangle()->SetCoordAxesVisible(bEnable);
}


void TazDlg::EnableBZ(bool bEnable)
{
	m_sceneRecip.GetTriangle()->SetBZVisible(bEnable);
}


void TazDlg::EnableWS(bool bEnable)
{
	m_sceneRealLattice.GetLattice()->SetWSVisible(bEnable);
}


void TazDlg::ShowAllPeaks(bool bShowAllPeaks)
{
	m_sceneRecip.GetTriangle()->SetAllPeaksVisible(bShowAllPeaks);
	CalcPeaks();
}


void TazDlg::ShowEwaldSphere()
{
	EwaldSphere iEw = EWALD_NONE;
	if(m_pEwaldSphereNone->isChecked())
		iEw = EWALD_NONE;
	else if(m_pEwaldSphereKi->isChecked())
		iEw = EWALD_KI;
	else if(m_pEwaldSphereKf->isChecked())
		iEw = EWALD_KF;
	m_sceneRecip.GetTriangle()->SetEwaldSphereVisible(iEw);
}


void TazDlg::EnableRealQDir(bool bEnable)
{
	m_sceneReal.GetTasLayout()->SetRealQVisible(bEnable);
	m_sceneTof.GetTofLayout()->SetRealQVisible(bEnable);
}


// Q position
void TazDlg::recipParamsChanged(const RecipParams& params)
{
	t_real dQx = -params.Q_rlu[0], dQy = -params.Q_rlu[1], dQz = -params.Q_rlu[2];
	t_real dE = params.dE;

	tl::set_eps_0(dQx, g_dEps); tl::set_eps_0(dQy, g_dEps); tl::set_eps_0(dQz, g_dEps);
	tl::set_eps_0(dE, g_dEps);

	std::ostringstream ostrPos;
	ostrPos.precision(g_iPrecGfx);
	ostrPos << "Q = (" << dQx << ", " << dQy << ", " << dQz  << ") rlu";
	ostrPos << ", E = " << dE << " meV";

	ostrPos << ", BZ: ("
		<< params.G_rlu_accurate[0] << ", "
		<< params.G_rlu_accurate[1] << ", "
		<< params.G_rlu_accurate[2] << ")";

	m_pCoordQStatusMsg->setText(ostrPos.str().c_str());
}


// cursor position
void TazDlg::RecipCoordsChanged(t_real dh, t_real dk, t_real dl,
	bool bHasNearest, t_real dNearestH, t_real dNearestK, t_real dNearestL)
{
	tl::set_eps_0(dh, g_dEps); tl::set_eps_0(dk, g_dEps); tl::set_eps_0(dl, g_dEps);
	tl::set_eps_0(dNearestH, g_dEps); tl::set_eps_0(dNearestK, g_dEps); tl::set_eps_0(dNearestL, g_dEps);

	std::ostringstream ostrPos;
	ostrPos.precision(g_iPrecGfx);
	ostrPos << "Cur: (" << dh << ", " << dk << ", " << dl  << ") rlu";
	if(bHasNearest)
		ostrPos << ", BZ: ("
			<< dNearestH << ", " << dNearestK << ", " << dNearestL << ")";

	m_pCoordCursorStatusMsg->setText(ostrPos.str().c_str());
}


// cursor position
void TazDlg::RealCoordsChanged(t_real dh, t_real dk, t_real dl,
	bool bHasNearest, t_real dNearestH, t_real dNearestK, t_real dNearestL)
{
	tl::set_eps_0(dh, g_dEps); tl::set_eps_0(dk, g_dEps); tl::set_eps_0(dl, g_dEps);
	tl::set_eps_0(dNearestH, g_dEps); tl::set_eps_0(dNearestK, g_dEps); tl::set_eps_0(dNearestL, g_dEps);

	std::ostringstream ostrPos;
	ostrPos.precision(g_iPrecGfx);
	ostrPos << "Cur: (" << dh << ", " << dk << ", " << dl  << ") frac";
	if(bHasNearest)
		ostrPos << ", WS: ("
			<< dNearestH << ", " << dNearestK << ", " << dNearestL << ")";

	m_pCoordCursorStatusMsg->setText(ostrPos.str().c_str());
}



//--------------------------------------------------------------------------------
// parameter dialogs

void TazDlg::ShowRecipParams()
{
	focus_dlg(&m_dlgRecipParam);
}


void TazDlg::ShowRealParams()
{
	focus_dlg(&m_dlgRealParam);
}


//--------------------------------------------------------------------------------
// dialogs

void TazDlg::ShowNeutronDlg()
{
	if(!m_pNeutronDlg)
	{
		m_pNeutronDlg = new NeutronDlg(this, &m_settings);
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged,
						 m_pNeutronDlg, &NeutronDlg::paramsChanged);
		m_sceneRecip.emitAllParams();
	}

	focus_dlg(m_pNeutronDlg);
}


void TazDlg::ShowTofDlg()
{
	if(!m_pTofDlg)
		m_pTofDlg = new TOFDlg(this, &m_settings);

	focus_dlg(m_pTofDlg);
}


void TazDlg::InitGoto()
{
	if(!m_pGotoDlg)
	{
		m_pGotoDlg = new GotoDlg(this, &m_settings);
		m_pGotoDlg->SetD(editMonoD->text().toDouble(), editAnaD->text().toDouble());
		m_pGotoDlg->SetSenses(checkSenseM->isChecked(), checkSenseS->isChecked(), checkSenseA->isChecked());
		m_pGotoDlg->SetLattice(m_latticecommon.lattice);
		m_pGotoDlg->SetScatteringPlane(m_latticecommon.dir0RLU, m_latticecommon.dir1RLU);
	}
}


void TazDlg::ShowGotoDlg()
{
	InitGoto();
	focus_dlg(m_pGotoDlg);
}


void TazDlg::InitElasticDlg()
{
	if(!m_pElasticDlg)
	{
		m_pElasticDlg = new ElasticDlg(this, &m_settings);
		QObject::connect(m_pElasticDlg, &ElasticDlg::ChangedPosition, this, &TazDlg::VarsChanged);

		m_pElasticDlg->SetD(editMonoD->text().toDouble(), editAnaD->text().toDouble());
		m_pElasticDlg->SetSenses(checkSenseM->isChecked(), checkSenseS->isChecked(), checkSenseA->isChecked());
		m_pElasticDlg->SetLattice(m_latticecommon.lattice);
		m_pElasticDlg->SetScatteringPlane(m_latticecommon.dir0RLU, m_latticecommon.dir1RLU);
	}
}


void TazDlg::ShowElasticDlg()
{
	InitElasticDlg();
	focus_dlg(m_pElasticDlg);
}


void TazDlg::ShowPowderDlg()
{
	if(!m_pPowderDlg)
	{
		m_pPowderDlg = new PowderDlg(this, &m_settings);
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged,
						 m_pPowderDlg, &PowderDlg::paramsChanged);
		m_sceneRecip.emitAllParams();
	}

	focus_dlg(m_pPowderDlg);
}


void TazDlg::ShowSettingsDlg()
{
	if(!m_pSettingsDlg)
		m_pSettingsDlg = new SettingsDlg(this, &m_settings);

	focus_dlg(m_pSettingsDlg);
}


void TazDlg::ShowScatteringFactorsDlg()
{
	if(!m_pScatteringFactorsDlg)
		m_pScatteringFactorsDlg = new ScatteringFactorsDlg(this, &m_settings);

	focus_dlg(m_pScatteringFactorsDlg);
}


void TazDlg::ShowDynPlaneDlg()
{
	if(!m_pDynPlaneDlg)
	{
		m_pDynPlaneDlg = new DynPlaneDlg(this, &m_settings);
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged,
						 m_pDynPlaneDlg, &DynPlaneDlg::RecipParamsChanged);
		m_sceneRecip.emitAllParams();
	}

	focus_dlg(m_pDynPlaneDlg);
}


void TazDlg::ShowScanViewer()
{
	if(!m_pScanViewer)
		m_pScanViewer = new ScanViewerDlg(this, &m_settings);

	focus_dlg(m_pScanViewer);
}


void TazDlg::ShowScanPos()
{
	if(!m_pScanPos)
		m_pScanPos = new ScanPosDlg(this, &m_settings);

	focus_dlg(m_pScanPos);
}


void TazDlg::ShowPowderFit()
{
	if(!m_pPowderFit)
		m_pPowderFit = new PowderFitDlg(this, &m_settings);

	focus_dlg(m_pPowderFit);
}


//--------------------------------------------------------------------------------
// context menus

void TazDlg::RecipContextMenu(const QPoint& _pt)
{
	if(!m_pviewRecip) return;

	QPoint pt = this->m_pviewRecip->mapToGlobal(_pt);
	m_pMenuViewRecip->exec(pt);
}


void TazDlg::RealContextMenu(const QPoint& _pt)
{
	if(!m_pviewReal) return;

	QPoint pt = this->m_pviewReal->mapToGlobal(_pt);
	m_pMenuViewReal->exec(pt);
}



//--------------------------------------------------------------------------------
// obstacles

void TazDlg::InitDarkAngles()
{
	if(!m_pDarkAnglesDlg)
	{
		m_pDarkAnglesDlg = new DarkAnglesDlg(this, &m_settings);
		QObject::connect(m_pDarkAnglesDlg, &DarkAnglesDlg::ApplyDarkAngles,
			this, &TazDlg::ApplyDarkAngles);
	}
}


void TazDlg::ShowDarkAnglesDlg()
{
	InitDarkAngles();
	m_pDarkAnglesDlg->SetDarkAngles(m_vecDarkAngles);
	focus_dlg(m_pDarkAnglesDlg);
}


void TazDlg::ApplyDarkAngles(const std::vector<DarkAngle<t_real>>& vecAngles)
{
	m_vecDarkAngles = vecAngles;
	if(m_sceneReal.GetTasLayout())
		m_sceneReal.GetTasLayout()->SetDarkAngles(&m_vecDarkAngles);
}



//--------------------------------------------------------------------------------
// about, log & help dialogs

void TazDlg::ShowExternalLicenses()
{
	std::vector<std::string> license_dirs = find_resource_dirs("3rdparty_licenses");
	if(license_dirs.size() == 0)
	{
		QMessageBox::critical(this, "Error", "3rd party license directory could not be found.");
		return;
	}

	std::string file = "file:///" + fs::absolute(*license_dirs.begin()).string();
	QDesktopServices::openUrl(QUrl(file.c_str()));
}


void TazDlg::ShowAbout()
{
	if(!m_pAboutDlg)
		m_pAboutDlg = new AboutDlg(this, &m_settings);

	focus_dlg(m_pAboutDlg);
}


void TazDlg::ShowWebsite()
{
	QDesktopServices::openUrl(QUrl("http://www.ill.eu/takin"));
}


void TazDlg::ReportBug()
{
	QDesktopServices::openUrl(QUrl("https://github.com/ILLGrenoble/takin/issues"));
}


void TazDlg::ShowLog()
{
	if(!m_pLogDlg)
		m_pLogDlg = new LogDlg(this, &m_settings, m_strLogFile);

	focus_dlg(m_pLogDlg);
}


void TazDlg::ShowHelp()
{
	QDesktopServices::openUrl(QUrl("https://github.com/ILLGrenoble/takin/wiki"));
}


void TazDlg::ShowDevelDoc()
{
	std::string strHelpHtml = find_resource("doc/html/index.html");
	if(strHelpHtml != "")
	{
		std::string strFile = "file:///" + fs::absolute(strHelpHtml).string();
		if(QDesktopServices::openUrl(QUrl(strFile.c_str())))
			return;
	}

	QMessageBox::critical(this, "Error", "Documentation could not be displayed.");
}


#include "moc_taz.cpp"
