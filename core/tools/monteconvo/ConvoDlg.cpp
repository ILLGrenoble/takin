/**
 * monte carlo convolution tool
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015 - 2023
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

#include "ConvoDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/math/math.h"
#include "tlibs/math/rand.h"
#include "tlibs/file/file.h"

#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "libs/qt/recent.h"

#include <boost/predef.h>

#include <iostream>
#include <fstream>
#include <tuple>

#include <QFileDialog>
#include <QMessageBox>
#include <QMenu>


using t_real = t_real_reso;
const std::string ConvoDlg::s_strTitle = "Resolution Convolution";


ConvoDlg::ConvoDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent, Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::WindowMinMaxButtonsHint),
		m_pSett(pSett)
{
	if(m_pSett)
	{
		// make sure the maximum thread and process count is set up
		// (this might otherwise not be the case when starting monteconvo individually)
		if(m_pSett->contains("main/max_threads"))
			g_iMaxThreads = m_pSett->value("main/max_threads", g_iMaxThreads).toInt();
		if(m_pSett->contains("main/max_processes"))
			g_iMaxProcesses = m_pSett->value("main/max_processes", g_iMaxProcesses).toInt();
	}

	setupUi(this);
	setWindowTitle(s_strTitle.c_str());

	// -------------------------------------------------------------------------
	// default values
	editEpsRlu->setText(tl::var_to_str(EPS_RLU).c_str());
	editEpsPlaneDist->setText(tl::var_to_str(EPS_PLANE).c_str());
	// -------------------------------------------------------------------------

	// -------------------------------------------------------------------------
	// widgets
	m_vecSpinBoxes = { spinStartH, spinStartK, spinStartL, spinStartE,
		spinStopH, spinStopK, spinStopL, spinStopE,
		spinStopH2, spinStopK2, spinStopL2, spinStopE2,
		spinKfix,
		spinTolerance
	};

	m_vecSpinNames = {
		"monteconvo/h_from", "monteconvo/k_from", "monteconvo/l_from", "monteconvo/E_from",
		"monteconvo/h_to", "monteconvo/k_to", "monteconvo/l_to", "monteconvo/E_to",
		"monteconvo/h_to_2", "monteconvo/k_to_2", "monteconvo/l_to_2", "monteconvo/E_to_2",
		"monteconvo/kfix",
		"convofit/tolerance"
	};

	m_vecIntSpinBoxes = { spinNeutrons, spinSampleSteps,
		spinStepCnt,
		spinStrategy, spinMaxCalls
	};
	m_vecIntSpinNames = {
		"monteconvo/neutron_count", "monteconvo/sample_step_count",
		"monteconvo/step_count",
		"convofit/strategy", "convofit/max_calls"
	};

	m_vecEditBoxes = {
		editFilterCol, editFilterVal,
		editCrys, editRes, editSqw, editScan,
		editScale, editSlope, editOffs,
		editCounter, editMonitor,
		editTemp, editField,
		editAutosave,
		editEpsRlu, editEpsPlaneDist,
	};
	m_vecEditNames = {
		"monteconvo/filter_col", "monteconvo/filter_val",
		"monteconvo/crys", "monteconvo/instr",
		"monteconvo/sqw_conf", "monteconvo/scanfile",
		"monteconvo/S_scale", "monteconvo/S_slope", "monteconvo/S_offs",
		"convofit/counter", "convofit/monitor",
		"convofit/temp_override", "convofit/field_override",
		"monteconvo/autosave",
		"monteconvo/eps_rlu", "monteconvo/eps_plane_dist",
	};

	m_vecTextBoxes = { editSqwParams };
	m_vecTextNames = { "convofit/sqw_params" };

	m_vecComboBoxes = { comboAlgo,
		comboFixedK, comboFocMono, comboFocAna,
		comboFitter, comboAxis, comboAxis2,
		comboRnd
	};
	m_vecComboNames = { "monteconvo/algo_idx",
		"monteconvo/fixedk", "monteconvo/mono_foc", "monteconvo/ana_foc",
		"convofit/minimiser", "convofit/scanaxis", "convofit/scanaxis2",
		"convofit/recycle_neutrons"
	};

	m_vecCheckBoxes = { checkScan, check2dMap, checkNorm, checkFlip };
	m_vecCheckNames = { "monteconvo/has_scanfile", "monteconvo/scan_2d",
		"convofit/normalise", "convofit/flip_coords"
	};
	// -------------------------------------------------------------------------

	if(m_pSett)
	{
		QFont font;
		if(m_pSett->contains("main/font_gen") && font.fromString(m_pSett->value("main/font_gen", "").toString()))
			setFont(font);
	}

	btnStart->setIcon(load_icon("res/icons/media-playback-start.svg"));
	btnStartFit->setIcon(load_icon("res/icons/media-playback-start.svg"));
	btnStop->setIcon(load_icon("res/icons/media-playback-stop.svg"));

	/*
	 * curve 0,1	->	convolution
	 * curve 2	->	scan points
	 * curve 3-15	->	dispersion branches
	 */
	m_plotwrap.reset(new QwtPlotWrapper(plot, CONVO_MAX_CURVES, true));
	m_plotwrap->GetPlot()->setAxisTitle(QwtPlot::xBottom, "");
	m_plotwrap->GetPlot()->setAxisTitle(QwtPlot::yLeft, "S(Q,E) (a.u.)");

	m_plotwrap2d.reset(new QwtPlotWrapper(plot2d, 1, 0, 0, 1));
	m_plotwrap2d->GetPlot()->setAxisTitle(QwtPlot::yRight, "S(Q,E) (a.u.)");


	// --------------------------------------------------------------------
	// convolution lines
	QPen penCurve;
	penCurve.setColor(QColor(0,0,0x99));
	penCurve.setWidth(2);
	m_plotwrap->GetCurve(0)->setPen(penCurve);
	m_plotwrap->GetCurve(0)->setStyle(QwtPlotCurve::CurveStyle::Lines);
	m_plotwrap->GetCurve(0)->setTitle("S(Q,E)");

	// convolution points
	QPen penPoints;
	penPoints.setColor(QColor(0xff,0,0));
	penPoints.setWidth(4);
	m_plotwrap->GetCurve(1)->setPen(penPoints);
	m_plotwrap->GetCurve(1)->setStyle(QwtPlotCurve::CurveStyle::Dots);
	m_plotwrap->GetCurve(1)->setTitle("S(Q,E)");

	// scan data points
	QPen penScanPoints;
	penScanPoints.setColor(QColor(0x00,0x90,0x00));
	penScanPoints.setWidth(6);
	m_plotwrap->GetCurve(2)->setPen(penScanPoints);
	m_plotwrap->GetCurve(2)->setStyle(QwtPlotCurve::CurveStyle::Dots);
	m_plotwrap->GetCurve(2)->setTitle("S(Q,E)");
	m_plotwrap->GetCurve(2)->SetShowErrors(true);

	// dispersion branches
	for(int iCurve=CONVO_DISP_CURVE_START; iCurve<CONVO_MAX_CURVES; ++iCurve)
	{
		m_plotwrap->GetCurve(iCurve)->setPen(penCurve);
		m_plotwrap->GetCurve(iCurve)->setStyle(QwtPlotCurve::CurveStyle::Lines);
		m_plotwrap->GetCurve(iCurve)->setTitle("E(Q)");
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	QMenu *pMenuActions = new QMenu("Actions", this);

	QAction *pHK = new QAction("h <-> k", pMenuActions);
	QAction *pHL = new QAction("h <-> l", pMenuActions);
	QAction *pKL = new QAction("k <-> l", pMenuActions);
	pMenuActions->addAction(pHK);
	pMenuActions->addAction(pHL);
	pMenuActions->addAction(pKL);

	btnActions->setMenu(pMenuActions);
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// fill in resolution calculation methods
	comboAlgo->addItem("TAS: Cooper-Nathans (Pointlike)", static_cast<int>(ResoAlgo::CN));
	comboAlgo->addItem("TAS: Popovici (Pointlike)", static_cast<int>(ResoAlgo::POP_CN));
	comboAlgo->insertSeparator(2);
	comboAlgo->addItem("TAS: Popovici", static_cast<int>(ResoAlgo::POP));
	comboAlgo->addItem("TAS: Eckold-Sobolev", static_cast<int>(ResoAlgo::ECK));
	comboAlgo->addItem("TAS: Eckold-Sobolev (Extended)", static_cast<int>(ResoAlgo::ECK_EXT));
	comboAlgo->insertSeparator(6);
	comboAlgo->addItem("TOF: Violini", static_cast<int>(ResoAlgo::VIO));

	comboAlgo->setCurrentIndex(3);  // default: pop
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// fill sqw combo box
	load_sqw_plugins();
	auto vecSqwNames = get_sqw_names();
	for(const auto& tupSqw : vecSqwNames)
	{
		const std::string& strIdent = std::get<0>(tupSqw);
		const std::string& strName = std::get<1>(tupSqw);
		const std::string& strHelp = std::get<2>(tupSqw);

		// add module item to combo box
		comboSqw->addItem(strName.c_str(), strIdent.c_str());

		// add module help text
		m_help_texts.insert(std::make_pair(strIdent, strHelp));
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// menu bar
	m_pMenuBar = new QMenuBar(this);
	if(m_pSett)
		m_pMenuBar->setNativeMenuBar(m_pSett->value("main/native_dialogs", 1).toBool());


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
	RecentFiles recent(m_pSett, "monteconvo/recent");
	recent.FillMenu(m_pMenuRecent, [this](const std::string& str){ Load(str.c_str()); });
	pMenuFile->addMenu(m_pMenuRecent);

	QAction *pSave = new QAction("Save", this);
	pSave->setIcon(load_icon("res/icons/document-save.svg"));
	pSave->setShortcut(QKeySequence::Save);
	pMenuFile->addAction(pSave);

	QAction *pSaveAs = new QAction("Save As...", this);
	pSaveAs->setIcon(load_icon("res/icons/document-save-as.svg"));
	pSaveAs->setShortcut(QKeySequence::SaveAs);
	pMenuFile->addAction(pSaveAs);

	pMenuFile->addSeparator();

	QAction *pConvofit = new QAction("Export to Convofit...", this);
	pConvofit->setIcon(load_icon("res/icons/drive-harddisk.svg"));
	pMenuFile->addAction(pConvofit);

	pMenuFile->addSeparator();

	QAction *pExit = new QAction("Close Convo", this);
	//pExit->setMenuRole(QAction::QuitRole);
	pExit->setIcon(load_icon("res/icons/system-log-out.svg"));
	pExit->setShortcut(QKeySequence::Close);
	pMenuFile->addAction(pExit);


	// view menu
	QMenu *pMenuView = new QMenu("View", this);

	QAction *pActionLogY = new QAction("Toggle Logarithmic Scale", this);
	pMenuView->addAction(pActionLogY);


	// actions menu
	QMenu *pMenuConvoActions = new QMenu("Actions", this);

	QAction *pActionStart = new QAction("Start Convolution Simulation", this);
	pActionStart->setIcon(load_icon("res/icons/media-playback-start.svg"));
	pMenuConvoActions->addAction(pActionStart);

	QAction *pActionStartFit = new QAction("Start Convolution Fit", this);
	pActionStartFit->setIcon(load_icon("res/icons/media-playback-start.svg"));
	pMenuConvoActions->addAction(pActionStartFit);

	QAction *pActionDisp = new QAction("Plot Dispersion", this);
	pMenuConvoActions->addAction(pActionDisp);


	// results menu
	QMenu *pMenuPlots = new QMenu("Results", this);

	m_pLiveResults = new QAction("Live Results", this);
	m_pLiveResults->setCheckable(true);
	m_pLiveResults->setChecked(false);
	pMenuPlots->addAction(m_pLiveResults);

	m_pLivePlots = new QAction("Live Plots", this);
	m_pLivePlots->setCheckable(true);
	m_pLivePlots->setChecked(true);
	pMenuPlots->addAction(m_pLivePlots);

	pMenuPlots->addSeparator();

	QAction *pExportPlot = new QAction("Export Plot Data...", this);
	pMenuPlots->addAction(pExportPlot);

	QAction *pExportPlotGpl = new QAction("Export Plot to Gnuplot...", this);
	pMenuPlots->addAction(pExportPlotGpl);

	pMenuPlots->addSeparator();

	QAction *pExportPlot2d = new QAction("Export 2D Plot Data...", this);
	pMenuPlots->addAction(pExportPlot2d);

	QAction *pExportPlot2dGpl = new QAction("Export 2D Plot to Gnuplot...", this);
	pMenuPlots->addAction(pExportPlot2dGpl);

	pMenuPlots->addSeparator();

	QAction *pSaveResults = new QAction("Save Results...", this);
	pSaveResults->setIcon(load_icon("res/icons/document-save-as.svg"));
	pMenuPlots->addAction(pSaveResults);

	// help menu
	QMenu *pMenuHelp = new QMenu("Help", this);

	QAction *pAbout = new QAction("About...", this);
	pAbout->setMenuRole(QAction::AboutRole);
	pAbout->setIcon(load_icon("res/icons/dialog-information.svg"));
	pMenuHelp->addAction(pAbout);


	m_pMenuBar->addMenu(pMenuFile);
	m_pMenuBar->addMenu(pMenuView);
	m_pMenuBar->addMenu(pMenuConvoActions);
	m_pMenuBar->addMenu(pMenuPlots);
	m_pMenuBar->addMenu(pMenuHelp);


	connect(pExit, &QAction::triggered, this, &ConvoDlg::accept);
	connect(pNew, &QAction::triggered, this, &ConvoDlg::New);
	connect(pLoad, &QAction::triggered, this, static_cast<void (ConvoDlg::*)()>(&ConvoDlg::Load));
	connect(pSave, &QAction::triggered, this, static_cast<void (ConvoDlg::*)()>(&ConvoDlg::Save));
	connect(pSaveAs, &QAction::triggered, this, &ConvoDlg::SaveAs);
	connect(pConvofit, &QAction::triggered, this, &ConvoDlg::SaveConvofit);
	connect(pActionLogY, &QAction::triggered, m_plotwrap.get(), &QwtPlotWrapper::ToggleLogY);
	connect(pActionStart, &QAction::triggered, this, &ConvoDlg::Start);
	connect(pActionStartFit, &QAction::triggered, this, &ConvoDlg::StartFit);
	connect(pActionDisp, &QAction::triggered, this, &ConvoDlg::StartDisp);
	connect(pExportPlot, &QAction::triggered, m_plotwrap.get(), &QwtPlotWrapper::SavePlot);
	connect(pExportPlot2d, &QAction::triggered, m_plotwrap2d.get(), &QwtPlotWrapper::SavePlot);
	connect(pExportPlotGpl, &QAction::triggered, m_plotwrap.get(), &QwtPlotWrapper::ExportGpl);
	connect(pExportPlot2dGpl, &QAction::triggered, m_plotwrap2d.get(), &QwtPlotWrapper::ExportGpl);
	connect(pAbout, &QAction::triggered, this, &ConvoDlg::ShowAboutDlg);
	connect(pSaveResults, &QAction::triggered, [this]() { SaveResult(); });

	this->layout()->setMenuBar(m_pMenuBar);
	// --------------------------------------------------------------------


	m_pSqwParamDlg = new SqwParamDlg(this, m_pSett);
	connect(this, &ConvoDlg::SqwLoaded, m_pSqwParamDlg, &SqwParamDlg::SqwLoaded);
	connect(m_pSqwParamDlg, &SqwParamDlg::SqwParamsChanged, this, &ConvoDlg::SqwParamsChanged);

	m_pFavDlg = new FavDlg(this, m_pSett);
	connect(m_pFavDlg, &FavDlg::ChangePos, this, &ConvoDlg::ChangePos);

	connect(btnBrowseCrys, &QAbstractButton::clicked, this, &ConvoDlg::browseCrysFiles);
	connect(btnBrowseRes, &QAbstractButton::clicked, this, &ConvoDlg::browseResoFiles);
	connect(btnBrowseSqw, &QAbstractButton::clicked, this, &ConvoDlg::browseSqwFiles);
	connect(btnBrowseScan, &QAbstractButton::clicked, this, &ConvoDlg::browseScanFiles);
	connect(btnBrowseAutosave, &QAbstractButton::clicked, this, &ConvoDlg::browseAutosaveFile);
	connect(btnFav, &QAbstractButton::clicked, this, &ConvoDlg::ShowFavourites);
	connect(btnSqwParams, &QAbstractButton::clicked, this, &ConvoDlg::showSqwParamDlg);
	connect(btnSqwHelp, &QAbstractButton::clicked, this, &ConvoDlg::showSqwHelpDlg);
	connect(btnSaveResults, &QAbstractButton::clicked,
	[this]()
	{
		QString outfile = editAutosave->text();
		SaveResult(&outfile);
	});

	connect(comboSqw, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &ConvoDlg::SqwModelChanged);
	connect(editSqw, &QLineEdit::textChanged, this, &ConvoDlg::createSqwModel);
	connect(btnSqwReload, &QAbstractButton::clicked, [this]()
	{ createSqwModel(editSqw->text()); });

	connect(editScan, &QLineEdit::textChanged, this, &ConvoDlg::scanFileChanged);
	connect(editFilterCol, &QLineEdit::textChanged, this, [this]() -> void
	{ scanFileChanged(editScan->text()); });
	connect(editFilterVal, &QLineEdit::textChanged, this, [this]() -> void
	{ scanFileChanged(editScan->text()); });

	connect(editCounter, &QLineEdit::textChanged, this, [this]() -> void
	{ scanFileChanged(editScan->text()); });
	connect(editMonitor, &QLineEdit::textChanged, this, [this]() -> void
	{ scanFileChanged(editScan->text()); });

	connect(editEpsRlu, &QLineEdit::textChanged, this, [this]() -> void
	{ m_eps_rlu = tl::str_to_var<t_real>(editEpsRlu->text().toStdString()); });
	connect(editEpsPlaneDist, &QLineEdit::textChanged, this, [this]() -> void
	{ m_eps_plane = tl::str_to_var<t_real>(editEpsPlaneDist->text().toStdString()); });

	connect(editScale, &QLineEdit::textChanged, this, &ConvoDlg::scaleChanged);
	connect(editSlope, &QLineEdit::textChanged, this, &ConvoDlg::scaleChanged);
	connect(editOffs, &QLineEdit::textChanged, this, &ConvoDlg::scaleChanged);

	connect(btnStart, &QAbstractButton::clicked, this, &ConvoDlg::Start);
	connect(btnStartFit, &QAbstractButton::clicked, this, &ConvoDlg::StartFit);
	connect(btnStop, &QAbstractButton::clicked, this, &ConvoDlg::Stop);

	connect(checkScan, &QCheckBox::toggled, this, &ConvoDlg::scanCheckToggled);
	connect(checkFlip, &QCheckBox::toggled, this, &ConvoDlg::coordFlipToggled);

	connect(pHK, &QAction::triggered, this, &ConvoDlg::ChangeHK);
	connect(pHL, &QAction::triggered, this, &ConvoDlg::ChangeHL);
	connect(pKL, &QAction::triggered, this, &ConvoDlg::ChangeKL);

	for(QDoubleSpinBox* pSpin : {spinStartH, spinStartK, spinStartL, spinStartE,
		spinStopH, spinStopK, spinStopL, spinStopE})
	{
		connect(pSpin, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			this, &ConvoDlg::UpdateCurFavPos);
	}

	LoadSettings();


#if BOOST_OS_MACOS
	// check if system python is available
	if(!tl::dir_exists("/Library/Frameworks/Python.framework")
		&& find_resource_dirs("Frameworks/Python.framework", false).size()==0)
	{
		QMessageBox::information(this, "Python Module",
			"The <i>Python</i> S(Q,E) plugin module requires having the "
			"<a href=\"https://www.python.org/downloads/mac-osx/\">Python framework</a> installed.<br><br>"
			"Please also install the <i>numpy</i> and <i>scipy</i> packages using the following command:<br><br>"
			"<code>/Library/Frameworks/Python.framework/Versions/Current/bin/pip3 install numpy scipy</code>");
	}
#endif
}


ConvoDlg::~ConvoDlg()
{
	WaitForThread();

	if(m_pSqwParamDlg)
	{
		delete m_pSqwParamDlg;
		m_pSqwParamDlg = nullptr;
	}

	if(m_pFavDlg)
	{
		delete m_pFavDlg;
		m_pFavDlg = nullptr;
		}

	if(m_pMenuBar)
	{
		delete m_pMenuBar;
		m_pMenuBar = nullptr;
	}

	if(m_pSqw)
		m_pSqw.reset();

	unload_sqw_plugins();
}


void ConvoDlg::WaitForThread()
{
	if(!m_pth)
		return;

	if(m_pth->joinable())
		m_pth->join();

	delete m_pth;
	m_pth = nullptr;
}


void ConvoDlg::SqwModelChanged(int)
{
	if(!m_bAllowSqwReinit)
		return;

	editSqw->clear();
	createSqwModel("");
}


void ConvoDlg::createSqwModel(const QString& qstrFile)
{
	if(!m_bAllowSqwReinit)
		return;

	if(m_pSqw)
	{
		m_pSqw.reset();
		emit SqwLoaded(std::vector<SqwBase::t_var>{}, nullptr);
	}

	// identifier of the currently selected module
	std::string strSqwIdent = comboSqw->itemData(comboSqw->currentIndex()).toString().toStdString();
	if(strSqwIdent == "")
		return;

	std::string _strSqwFile = qstrFile.toStdString();
	tl::trim(_strSqwFile);
	std::string strSqwFile = find_file_in_global_paths(_strSqwFile);

	if(strSqwFile == "")
	{
		tl::log_warn("No S(Q,E) config file given.");
		//return;
	}

	m_pSqw.reset();
	m_pSqw = construct_sqw(strSqwIdent, strSqwFile);
	if(!m_pSqw)
	{
		QMessageBox::critical(this, "Error", "Unknown S(Q,E) model selected.");
		return;
	}

	if(m_pSqw && m_pSqw->IsOk())
	{
		emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
	}
	else
	{
		//QMessageBox::critical(this, "Error", "Could not create S(Q,E).");
		tl::log_err("Could not create S(Q,E).");
		return;
	}
}


void ConvoDlg::SqwParamsChanged(const std::vector<SqwBase::t_var>& vecVars,
	const std::vector<SqwBase::t_var_fit>* pvecVarsFit)
{
	if(!m_pSqw)
		return;

	m_pSqw->SetVars(vecVars);
	if(pvecVarsFit)
		m_pSqw->InitFitVars(*pvecVarsFit);

#ifndef NDEBUG
	// check: read parameters back in
	emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
#endif
}


/**
 * set a model parameter
 */
void ConvoDlg::SetSqwParam(const std::string& name, t_real_reso val)
{
	m_pSqw->SetVarIfAvail(name, tl::var_to_str(val));

	// read parameters back in to update paramters dialog
	emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
}


/**
 * set model parameters
 * [ ident, value, error ]
 */
void ConvoDlg::SetSqwParams(const std::vector<std::tuple<std::string, std::string, std::string>>& sqwparams)
{
	for(const auto& param : sqwparams)
	{
		m_pSqw->SetVarIfAvail(std::get<0>(param), std::get<1>(param));
		if(std::get<2>(param) != "")
			m_pSqw->SetErrIfAvail(std::get<0>(param), std::get<2>(param));
		//if(std::get<3>(param) != "")
		//	m_pSqw->SetRangeIfAvail(std::get<0>(param), std::get<3>(param));
	}

	// read parameters back in to update paramters dialog
	emit SqwLoaded(m_pSqw->GetVars(), &m_pSqw->GetFitVars());
}


/**
 * get the model parameters
 * [ ident, type, value, error, fit? ]
 */
ConvoDlg::t_sqwparams ConvoDlg::GetSqwParams(bool only_fitparams) const
{
	t_sqwparams params;

	std::vector<SqwBase::t_var> vars1 = m_pSqw->GetVars();
	std::vector<SqwBase::t_var_fit> vars2 = m_pSqw->GetFitVars();

	for(const SqwBase::t_var& var : vars1)
	{
		const std::string& strName = std::get<SQW_NAME>(var);
		std::string strErr, strRange;
		bool bFit = 0;

		// look for associated fit parameters: match with basic variable ident
		auto iterFit = std::find_if(vars2.begin(), vars2.end(),
			[&strName](const SqwBase::t_var_fit& varFit) -> bool
			{ return strName == std::get<0>(varFit); });

		if(iterFit != vars2.end())
		{
			strErr = std::get<1>(*iterFit);		// error
			strRange = std::get<3>(*iterFit);	// range
			bFit = std::get<2>(*iterFit);		// "is fit param" flag
		}

		if((only_fitparams && bFit) || !only_fitparams)
		{
			params.emplace_back(std::make_tuple(
				strName, std::get<SQW_TYPE>(var), std::get<SQW_VAL>(var),
				strErr, bFit, strRange));
		}
	}

	return params;
}

// -----------------------------------------------------------------------------


/**
 * clear plot curves
 */
void ConvoDlg::ClearPlot1D()
{
	static const std::vector<t_real> vecZero;
	if(!m_plotwrap)
		return;

	for(std::size_t iCurve=0; iCurve<CONVO_MAX_CURVES; ++iCurve)
		set_qwt_data<t_real_reso>()(*m_plotwrap, vecZero, vecZero, iCurve, false);
}


/**
 * start 1d or 2d convolutions
 */
void ConvoDlg::Start()
{
	if(check2dMap->isChecked())
		Start2D();
	else
		Start1D();
}


/**
 * stop running operations
 */
void ConvoDlg::Stop()
{
	m_atStop.store(true);
}


// -----------------------------------------------------------------------------


void ConvoDlg::ShowFavourites()
{
	focus_dlg(m_pFavDlg);
}

void ConvoDlg::UpdateCurFavPos()
{
	FavHklPos pos;
	pos.dhstart = spinStartH->value();
	pos.dkstart = spinStartK->value();
	pos.dlstart = spinStartL->value();
	pos.dEstart = spinStartE->value();
	pos.dhstop = spinStopH->value();
	pos.dkstop = spinStopK->value();
	pos.dlstop = spinStopL->value();
	pos.dEstop = spinStopE->value();

	m_pFavDlg->UpdateCurPos(pos);
}

void ConvoDlg::ChangePos(const struct FavHklPos& pos)
{
	spinStartH->setValue(pos.dhstart);
	spinStartK->setValue(pos.dkstart);
	spinStartL->setValue(pos.dlstart);
	spinStartE->setValue(pos.dEstart);
	spinStopH->setValue(pos.dhstop);
	spinStopK->setValue(pos.dkstop);
	spinStopL->setValue(pos.dlstop);
	spinStopE->setValue(pos.dEstop);
}


// -----------------------------------------------------------------------------


static void SwapSpin(QDoubleSpinBox* pS1, QDoubleSpinBox* pS2)
{
	double dVal = pS1->value();
	pS1->setValue(pS2->value());
	pS2->setValue(dVal);
}

void ConvoDlg::ChangeHK()
{
	SwapSpin(spinStartH, spinStartK);
	SwapSpin(spinStopH, spinStopK);
}

void ConvoDlg::ChangeHL()
{
	SwapSpin(spinStartH, spinStartL);
	SwapSpin(spinStopH, spinStopL);
}

void ConvoDlg::ChangeKL()
{
	SwapSpin(spinStartK, spinStartL);
	SwapSpin(spinStopK, spinStopL);
}


// -----------------------------------------------------------------------------


void ConvoDlg::scanCheckToggled(bool bChecked)
{
	if(bChecked)
		scanFileChanged(editScan->text());
}


void ConvoDlg::coordFlipToggled(bool bChecked)
{
	scanFileChanged(editScan->text());
}


void ConvoDlg::scanFileChanged(const QString& qstrFile)
{
	m_bUseScan = false;
	if(!checkScan->isChecked())
		return;

	Filter filter;
	if(editFilterCol->text() != "")
	{
		filter.colEquals = std::make_pair(
			editFilterCol->text().toStdString(),
			editFilterVal->text().toStdString());
	}

	m_scan = Scan();

	// optional counter and monitor overrides
	m_scan.strCntCol = editCounter->text().toStdString();
	m_scan.strMonCol = editMonitor->text().toStdString();

	bool allow_scan_merging = false;
	if(m_pSett)
		allow_scan_merging = m_pSett->value("main/allow_scan_merging", 0).toBool();

	if(!load_scan_file(qstrFile.toStdString(), m_scan,
		checkFlip->isChecked(), allow_scan_merging, filter))
	{
		tl::log_err("Cannot load scan(s).");
		return;
	}

	if(!m_scan.vecPoints.size())
	{
		tl::log_err("No points in scan(s).");
		return;
	}

	// get scan start and end coordinates
	comboFixedK->setCurrentIndex(m_scan.bKiFixed ? 0 : 1);
	spinKfix->setValue(m_scan.dKFix);

	spinStartH->setValue(m_scan.vecScanOrigin[0]);
	spinStartK->setValue(m_scan.vecScanOrigin[1]);
	spinStartL->setValue(m_scan.vecScanOrigin[2]);
	spinStartE->setValue(m_scan.vecScanOrigin[3]);

	spinStopH->setValue(m_scan.vecScanOrigin[0] + m_scan.vecScanDir[0]);
	spinStopK->setValue(m_scan.vecScanOrigin[1] + m_scan.vecScanDir[1]);
	spinStopL->setValue(m_scan.vecScanOrigin[2] + m_scan.vecScanDir[2]);
	spinStopE->setValue(m_scan.vecScanOrigin[3] + m_scan.vecScanDir[3]);

	m_bUseScan = true;
}


void ConvoDlg::scaleChanged()
{
	// get scan x axis
	bool bScanAxisFound = 0;
	int iScanAxisIdx = 0;
	std::string strScanVar = "";
	std::vector<std::vector<t_real>> vecAxes;
	std::tie(bScanAxisFound, iScanAxisIdx, strScanVar, vecAxes) = get_scan_axis<t_real>(
		true, comboAxis->currentIndex(), spinStepCnt->value(), m_eps_rlu,
		spinStartH->value(), spinStopH->value(), spinStartK->value(), spinStopK->value(),
		spinStartL->value(), spinStopL->value(), spinStartE->value(), spinStopE->value());
	if(!bScanAxisFound)
	{
		tl::log_err("No scan variable found.");
		return;
	}


	t_real dScale = tl::str_to_var<t_real>(editScale->text().toStdString());
	t_real dSlope = tl::str_to_var<t_real>(editSlope->text().toStdString());
	t_real dOffs = tl::str_to_var<t_real>(editOffs->text().toStdString());

	m_vecScaledS.resize(m_vecS.size());
	for(std::size_t i=0; i<m_vecS.size(); ++i)
	{
		const t_real dXVal = vecAxes[iScanAxisIdx][i];
		m_vecScaledS[i] = dScale*(m_vecS[i] + dSlope*dXVal) + dOffs;
		if(m_vecScaledS[i] < 0.)
			m_vecScaledS[i] = 0.;
	}

	set_qwt_data<t_real_reso>()(*m_plotwrap, m_vecQ, m_vecScaledS, 0, false);
	set_qwt_data<t_real_reso>()(*m_plotwrap, m_vecQ, m_vecScaledS, 1, false);

	m_plotwrap->GetPlot()->replot();
}


// -----------------------------------------------------------------------------


void ConvoDlg::showSqwParamDlg()
{
	focus_dlg(m_pSqwParamDlg);
}


void ConvoDlg::showSqwHelpDlg()
{
	// identifier of the currently selected S(Q,E) module
	std::string strSqwIdent = comboSqw->itemData(comboSqw->currentIndex()).toString().toStdString();
	if(strSqwIdent == "")
		return;

	// S(Q,E) module name
	QString strSqwLongName = comboSqw->itemText(comboSqw->currentIndex());

	auto iter =  m_help_texts.find(strSqwIdent);
	if(iter == m_help_texts.end())
	{
		QMessageBox::warning(this, strSqwLongName, "No help text was found.");
		return;
	}

	// no help text
	if(iter->second == "")
	{
		QMessageBox::warning(this, strSqwLongName, "No help text was defined.");
		return;
	}

	// display help
	QMessageBox::information(this, strSqwLongName, iter->second.c_str());
}


#include "libs/version.h"

void ConvoDlg::ShowAboutDlg()
{
	std::ostringstream ostrAbout;
	ostrAbout << "Takin/Convo version " << TAKIN_VER << ".\n";
	ostrAbout << "Written by Tobias Weber <tweber@ill.fr>,\n";
	ostrAbout << "2015 - 2024.\n";
	ostrAbout << "\n" << TAKIN_LICENSE("Takin/Convo");

	QMessageBox::about(this, "About Convo", ostrAbout.str().c_str());
}


#include "moc_ConvoDlg.cpp"
