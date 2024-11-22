/**
 * magnetic dynamics -- main gui setup and handler functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#include "magdyn.h"
#include "libs/loadcif.h"

#include <QtCore/QMimeData>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLabel>
#include <QtGui/QDesktopServices>

#include <iostream>



MagDynDlg::MagDynDlg(QWidget* pParent) : QDialog{pParent},
	m_sett{new QSettings{"takin", "magdyn", this}}
{
	// restore settings done from takin main settings dialog
	get_settings_from_takin_core();
	if(g_font != "")
	{
		QFont font = this->font();
		if(font.fromString(g_font))
			setFont(font);
	}

	InitSettingsDlg();

	// read settings that require a restart
	m_allow_ortho_spin = (g_allow_ortho_spin != 0);
	m_allow_general_J = (g_allow_general_J != 0);

	// create gui
	CreateMainWindow();
	CreateMenuBar();

	// create dialogs
	ShowInfoDlg(true);
	ShowNotesDlg(true);

	// create input panels
	CreateSitesPanel();
	CreateExchangeTermsPanel();
	CreateSamplePanel();
	CreateSampleEnvPanel();
	CreateVariablesPanel();

	// create output panels
	CreateDispersionPanel();
	CreateHamiltonPanel();
	CreateCoordinatesPanel();
	CreateExportPanel();


	// get space groups and symops
	auto spacegroups = get_sgs<t_mat_real>();
	m_SGops.clear();
	m_SGops.reserve(spacegroups.size());
	m_comboSG->clear();

	for(auto [sgnum, descr, ops] : spacegroups)
	{
		m_comboSG->addItem(descr.c_str(), m_comboSG->count());
		m_SGops.emplace_back(std::move(ops));
	}


	InitSettings();

	// restore settings
	if(m_sett)
	{
		// restore window size and position
		if(m_sett->contains("geo"))
			restoreGeometry(m_sett->value("geo").toByteArray());
		else
			resize(800, 600);

		if(m_sett->contains("recent_files"))
			m_recent.SetRecentFiles(m_sett->value("recent_files").toStringList());

		if(m_sett->contains("recent_struct_files"))
			m_recent_struct.SetRecentFiles(m_sett->value("recent_struct_files").toStringList());

		if(m_sett->contains("splitter"))
			m_split_inout->restoreState(m_sett->value("splitter").toByteArray());
	}

	setAcceptDrops(true);
}



MagDynDlg::~MagDynDlg()
{
	Clear();

	if(m_settings_dlg)
	{
		delete m_settings_dlg;
		m_settings_dlg = nullptr;
	}

	if(m_structplot_dlg)
	{
		delete m_structplot_dlg;
		m_structplot_dlg = nullptr;
	}

	if(m_table_import_dlg)
	{
		delete m_table_import_dlg;
		m_table_import_dlg = nullptr;
	}

	if(m_notes_dlg)
	{
		delete m_notes_dlg;
		m_notes_dlg = nullptr;
	}

	if(m_info_dlg)
	{
		delete m_info_dlg;
		m_info_dlg = nullptr;
	}
}



void MagDynDlg::CreateMainWindow()
{
	SetCurrentFile("");
	setSizeGripEnabled(true);

	m_tabs_in = new QTabWidget(this);
	m_tabs_out = new QTabWidget(this);

	// fixed status
	m_statusFixed = new QLabel(this);
	m_statusFixed->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_statusFixed->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Preferred);
	m_statusFixed->setFrameShape(QFrame::Panel);
	m_statusFixed->setFrameShadow(QFrame::Sunken);
	m_statusFixed->setText("Ready.");

	// expanding status
	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
	m_status->setFrameShape(QFrame::Panel);
	m_status->setFrameShadow(QFrame::Sunken);

	// progress bar
	m_progress = new QProgressBar(this);
	m_progress->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	// start/stop button
	m_btnStartStop = new QPushButton("Calculate", this);
	m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-start"));
	m_btnStartStop->setToolTip("Start calculation.");
	m_btnStartStop->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	// show structure
	QPushButton *btnShowStruct = new QPushButton("View Structure...", this);
	btnShowStruct->setIcon(QIcon::fromTheme("applications-graphics"));
	btnShowStruct->setToolTip("Show a 3D view of the magnetic sites and couplings.");
	btnShowStruct->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	// splitter for input and output tabs
	m_split_inout = new QSplitter(this);
	m_split_inout->setOrientation(Qt::Horizontal);
	m_split_inout->setChildrenCollapsible(true);
	m_split_inout->addWidget(m_tabs_in);
	m_split_inout->addWidget(m_tabs_out);

	// main grid
	m_maingrid = new QGridLayout(this);
	m_maingrid->setSpacing(4);
	m_maingrid->setContentsMargins(8, 8, 8, 8);
	m_maingrid->addWidget(m_split_inout, 0,0, 1,8);
	m_maingrid->addWidget(m_statusFixed, 1,0, 1,1);
	m_maingrid->addWidget(m_status, 1,1, 1,3);
	m_maingrid->addWidget(m_progress, 1,4, 1,2);
	m_maingrid->addWidget(m_btnStartStop, 1,6, 1,1);
	m_maingrid->addWidget(btnShowStruct, 1,7, 1,1);

	// signals
	connect(m_btnStartStop, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_startEnabled)
			this->CalcAll();
		else
			m_stopRequested = true;
	});

	connect(btnShowStruct, &QAbstractButton::clicked, this, &MagDynDlg::ShowStructPlotDlg);
}



/**
 * main menu
 */
void MagDynDlg::CreateMenuBar()
{
	m_menu = new QMenuBar(this);

	// file menu
	QMenu *menuFile = new QMenu("File", m_menu);
	QAction *acNew = new QAction("New", menuFile);
	QAction *acLoad = new QAction("Open...", menuFile);
	QAction *acImportStructure = new QAction("Import Structure...", menuFile);
	QAction *acSave = new QAction("Save", menuFile);
	QAction *acSaveAs = new QAction("Save As...", menuFile);
	QAction *acExit = new QAction("Quit", menuFile);

	// structure menu
	QMenu *menuStruct = new QMenu("Structure", m_menu);
	QAction *acStructSymIdx = new QAction("Assign Symmetry Indices", menuStruct);
	QAction *acStructImport = new QAction("Import From Table...", menuStruct);
	QAction *acStructExportSun = new QAction("Export To Sunny Code...");
	QAction *acStructExportSW = new QAction("Export To SpinW Code...");
	QAction *acStructNotes = new QAction("Notes...", menuStruct);
	QAction *acStructView = new QAction("View...", menuStruct);
	QAction *acGroundState = new QAction("Minimise Ground State...", menuStruct);

	// dispersion menu
	m_menuDisp = new QMenu("Dispersion", m_menu);
	m_plot_channels = new QAction("Plot Channels", m_menuDisp);
	m_plot_channels->setToolTip("Plot individual polarisation channels.");
	m_plot_channels->setCheckable(true);
	m_plot_channels->setChecked(false);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuDisp);
	QAction *acSaveFigure = new QAction("Save Figure...", m_menuDisp);
	QAction *acSaveDisp = new QAction("Save Data...", m_menuDisp);
	QAction *acSaveMultiDisp = new QAction("Save Data For All Qs...", m_menuDisp);
	QAction *acSaveDispScr = new QAction("Save Data As Script...", m_menuDisp);
	QAction *acSaveMultiDispScr = new QAction("Save Data As Script For All Qs...", m_menuDisp);

	// channels sub-menu
	m_menuChannels = new QMenu("Selected Channels", m_menuDisp);
	m_plot_channel[0] = new QAction("Channel xx", m_menuChannels);
	m_plot_channel[1] = new QAction("Channel yy", m_menuChannels);
	m_plot_channel[2] = new QAction("Channel zz", m_menuChannels);
	for(int i = 0; i < 3; ++i)
	{
		m_plot_channel[i]->setCheckable(true);
		m_plot_channel[i]->setChecked(true);
	}
	m_menuChannels->addAction(m_plot_channel[0]);
	m_menuChannels->addAction(m_plot_channel[1]);
	m_menuChannels->addAction(m_plot_channel[2]);
	m_menuChannels->setEnabled(m_plot_channels->isChecked());

	// weight plot sub-menu
	QMenu *menuWeights = new QMenu("Plot Weights", m_menuDisp);
	m_plot_weights_pointsize = new QAction("As Point Size", menuWeights);
	m_plot_weights_alpha = new QAction("As Colour Alpha", menuWeights);
	m_plot_weights_pointsize->setCheckable(true);
	m_plot_weights_pointsize->setChecked(true);
	m_plot_weights_alpha->setCheckable(true);
	m_plot_weights_alpha->setChecked(false);
	menuWeights->addAction(m_plot_weights_pointsize);
	menuWeights->addAction(m_plot_weights_alpha);

	// recent files menus
	m_menuOpenRecent = new QMenu("Open Recent", menuFile);
	m_menuImportStructRecent = new QMenu("Import Recent", menuFile);

	// recently opened files
	m_recent.SetRecentFilesMenu(m_menuOpenRecent);
	m_recent.SetMaxRecentFiles(g_maxnum_recents);
	m_recent.SetOpenFunc(&m_open_func);

	// recently imported structure files
	m_recent_struct.SetRecentFilesMenu(m_menuImportStructRecent);
	m_recent_struct.SetMaxRecentFiles(g_maxnum_recents);
	m_recent_struct.SetOpenFunc(&m_import_struct_func);

	// shortcuts
	acNew->setShortcut(QKeySequence::New);
	acLoad->setShortcut(QKeySequence::Open);
	acSave->setShortcut(QKeySequence::Save);
	acSaveAs->setShortcut(QKeySequence::SaveAs);
	acExit->setShortcut(QKeySequence::Quit);
	acExit->setMenuRole(QAction::QuitRole);

	// icons
	acNew->setIcon(QIcon::fromTheme("document-new"));
	acLoad->setIcon(QIcon::fromTheme("document-open"));
	acSave->setIcon(QIcon::fromTheme("document-save"));
	acSaveAs->setIcon(QIcon::fromTheme("document-save-as"));
	acExit->setIcon(QIcon::fromTheme("application-exit"));
	m_menuOpenRecent->setIcon(QIcon::fromTheme("document-open-recent"));
	acSaveFigure->setIcon(QIcon::fromTheme("image-x-generic"));
	acSaveDisp->setIcon(QIcon::fromTheme("text-x-generic"));
	acSaveMultiDisp->setIcon(QIcon::fromTheme("text-x-generic"));
	acSaveDispScr->setIcon(QIcon::fromTheme("text-x-script"));
	acSaveMultiDispScr->setIcon(QIcon::fromTheme("text-x-script"));
	acStructExportSun->setIcon(QIcon::fromTheme("weather-clear"));
	acStructExportSW->setIcon(QIcon::fromTheme("text-x-script"));
	acStructNotes->setIcon(QIcon::fromTheme("accessories-text-editor"));
	acStructView->setIcon(QIcon::fromTheme("applications-graphics"));

	// calculation menu
	QMenu *menuCalc = new QMenu("Calculation", m_menu);
	m_autocalc = new QAction("Automatically Calculate", menuCalc);
	m_autocalc->setToolTip("Automatically calculate the results.");
	m_autocalc->setCheckable(true);
	m_autocalc->setChecked(false);
	QAction *acCalc = new QAction("Start Calculation", menuCalc);
	acCalc->setIcon(QIcon::fromTheme("media-playback-start"));
	acCalc->setToolTip("Calculate all results.");
	m_use_dmi = new QAction("Use DMI", menuCalc);
	m_use_dmi->setToolTip("Enables the Dzyaloshinskij-Moriya interaction.");
	m_use_dmi->setCheckable(true);
	m_use_dmi->setChecked(true);

	if(m_allow_general_J)
	{
		m_use_genJ = new QAction("Use General J", menuCalc);
		m_use_genJ->setToolTip("Enables the general interaction matrix.");
		m_use_genJ->setCheckable(true);
		m_use_genJ->setChecked(true);
	}

	m_use_field = new QAction("Use External Field", menuCalc);
	m_use_field->setToolTip("Enables an external field.");
	m_use_field->setCheckable(true);
	m_use_field->setChecked(true);
	m_use_temperature = new QAction("Use Bose Factor", menuCalc);
	m_use_temperature->setToolTip("Enables the Bose factor.");
	m_use_temperature->setCheckable(true);
	m_use_temperature->setChecked(true);
	m_use_formfact = new QAction("Use Form Factor", menuCalc);
	m_use_formfact->setToolTip("Enables the magnetic form factor.");
	m_use_formfact->setCheckable(true);
	m_use_formfact->setChecked(false);
	m_use_weights = new QAction("Use Neutron Spectral Weights", menuCalc);
	m_use_weights->setToolTip("Enables calculation of the spin correlation function.");
	m_use_weights->setCheckable(true);
	m_use_weights->setChecked(true);
	m_use_projector = new QAction("Use Neutron Projector", menuCalc);
	m_use_projector->setToolTip("Enables the neutron orthogonal projector.");
	m_use_projector->setCheckable(true);
	m_use_projector->setChecked(true);
	m_unite_degeneracies = new QAction("Unite Degenerate Energies", menuCalc);
	m_unite_degeneracies->setToolTip("Unites the weight factors corresponding to degenerate eigenenergies.");
	m_unite_degeneracies->setCheckable(true);
	m_unite_degeneracies->setChecked(true);
	m_ignore_annihilation = new QAction("Ignore Magnon Annihilation", menuCalc);
	m_ignore_annihilation->setToolTip("Calculate only magnon creation..");
	m_ignore_annihilation->setCheckable(true);
	m_ignore_annihilation->setChecked(false);
	m_force_incommensurate = new QAction("Force Incommensurate", menuCalc);
	m_force_incommensurate->setToolTip("Enforce incommensurate calculation even for commensurate magnetic structures.");
	m_force_incommensurate->setCheckable(true);
	m_force_incommensurate->setChecked(false);

	// H components sub-menu
	QMenu *menuHamiltonians = new QMenu("Selected Hamiltonians", menuCalc);
	m_hamiltonian_comp[0] = new QAction("H(Q)", menuHamiltonians);
	m_hamiltonian_comp[1] = new QAction("H(Q + O)", menuHamiltonians);
	m_hamiltonian_comp[2] = new QAction("H(Q - O)", menuHamiltonians);
	for(int i = 0; i < 3; ++i)
	{
		m_hamiltonian_comp[i]->setCheckable(true);
		m_hamiltonian_comp[i]->setChecked(true);
	}
	menuHamiltonians->addAction(m_hamiltonian_comp[0]);
	menuHamiltonians->addAction(m_hamiltonian_comp[1]);
	menuHamiltonians->addAction(m_hamiltonian_comp[2]);

	// tools menu
	QMenu *menuTools = new QMenu("Tools", m_menu);
	QAction *acTopo = new QAction("Topology...", menuTools);
	QAction *acTrafoCalc = new QAction("Transformations...", menuTools);
	QAction *acPreferences = new QAction("Preferences...", menuTools);
	acTrafoCalc->setIcon(QIcon::fromTheme("accessories-calculator"));
	acPreferences->setIcon(QIcon::fromTheme("preferences-system"));

	acPreferences->setShortcut(QKeySequence::Preferences);
	acPreferences->setMenuRole(QAction::PreferencesRole);

	// help menu
	QMenu *menuHelp = new QMenu("Help", m_menu);
	QAction *acHelp = new QAction(
		QIcon::fromTheme("help-contents"),
		"Show Help...", menuHelp);
	QAction *acAboutQt = new QAction(
		QIcon::fromTheme("help-about"),
		"About Qt...", menuHelp);
	QAction *acAbout = new QAction(
		QIcon::fromTheme("help-about"),
		"About...", menuHelp);

	acAboutQt->setMenuRole(QAction::AboutQtRole);
	acAbout->setMenuRole(QAction::AboutRole);

	// actions
	menuFile->addAction(acNew);
	menuFile->addSeparator();
	menuFile->addAction(acLoad);
	menuFile->addMenu(m_menuOpenRecent);
	menuFile->addSeparator();
	menuFile->addAction(acSave);
	menuFile->addAction(acSaveAs);
	menuFile->addSeparator();
	menuFile->addAction(acImportStructure);
	menuFile->addMenu(m_menuImportStructRecent);
	menuFile->addSeparator();
	menuFile->addAction(acExit);

	menuStruct->addAction(acStructSymIdx);
	menuStruct->addSeparator();
	menuStruct->addAction(acStructNotes);
	menuStruct->addSeparator();
	menuStruct->addAction(acStructView);
#ifdef __TLIBS2_MAGDYN_USE_MINUIT__
	menuStruct->addAction(acGroundState);
#endif
	menuStruct->addSeparator();
	menuStruct->addAction(acStructImport);
	menuStruct->addAction(acStructExportSun);
	menuStruct->addAction(acStructExportSW);

	m_menuDisp->addAction(m_plot_channels);
	m_menuDisp->addMenu(m_menuChannels);
	m_menuDisp->addSeparator();
	m_menuDisp->addAction(acRescalePlot);
	m_menuDisp->addMenu(menuWeights);
	m_menuDisp->addSeparator();
	m_menuDisp->addAction(acSaveFigure);
	m_menuDisp->addSeparator();
	m_menuDisp->addAction(acSaveDisp);
	m_menuDisp->addAction(acSaveMultiDisp);
	m_menuDisp->addSeparator();
	m_menuDisp->addAction(acSaveDispScr);
	m_menuDisp->addAction(acSaveMultiDispScr);

	menuCalc->addAction(m_autocalc);
	menuCalc->addAction(acCalc);
	menuCalc->addSeparator();
	menuCalc->addAction(m_use_dmi);
	if(m_allow_general_J)
		menuCalc->addAction(m_use_genJ);
	menuCalc->addAction(m_use_field);
	menuCalc->addAction(m_use_temperature);
	menuCalc->addAction(m_use_formfact);
	menuCalc->addSeparator();
	menuCalc->addAction(m_use_weights);
	menuCalc->addAction(m_use_projector);
	menuCalc->addSeparator();
	menuCalc->addAction(m_unite_degeneracies);
	menuCalc->addAction(m_ignore_annihilation);
	menuCalc->addAction(m_force_incommensurate);
	menuCalc->addMenu(menuHamiltonians);

	menuTools->addAction(acTopo);
	menuTools->addAction(acTrafoCalc);
	menuTools->addSeparator();
	menuTools->addAction(acPreferences);

	menuHelp->addAction(acHelp);
	menuHelp->addSeparator();
	menuHelp->addAction(acAboutQt);
	menuHelp->addAction(acAbout);

	// signals
	connect(acNew, &QAction::triggered, this, &MagDynDlg::Clear);
	connect(acLoad, &QAction::triggered,
		this, static_cast<void (MagDynDlg::*)()>(&MagDynDlg::Load));
	connect(acImportStructure, &QAction::triggered,
		this, static_cast<void (MagDynDlg::*)()>(&MagDynDlg::ImportStructure));
	connect(acSave, &QAction::triggered,
		this, static_cast<void (MagDynDlg::*)()>(&MagDynDlg::Save));
	connect(acSaveAs, &QAction::triggered,
		this, static_cast<void (MagDynDlg::*)()>(&MagDynDlg::SaveAs));
	connect(acExit, &QAction::triggered, this, &QDialog::close);

	connect(acSaveFigure, &QAction::triggered, this, &MagDynDlg::SavePlotFigure);
	connect(acSaveDisp, &QAction::triggered,
		[this](){ this->SaveDispersion(false); });
	connect(acSaveMultiDisp, &QAction::triggered,
		[this](){ this->SaveMultiDispersion(false); });
	connect(acSaveDispScr, &QAction::triggered,
		[this](){ this->SaveDispersion(true); });
	connect(acSaveMultiDispScr, &QAction::triggered,
		[this](){ this->SaveMultiDispersion(true); });

	connect(acRescalePlot, &QAction::triggered, [this]()
	{
		if(!m_plot)
			return;

		m_plot->rescaleAxes();
		m_plot->replot();
	});

	auto calc_all = [this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcAll();
	};

	auto calc_all_dyn = [this]()
	{
		if(this->m_autocalc->isChecked())
		{
			this->CalcDispersion();
			this->CalcHamiltonian();
		}
	};

	connect(acStructNotes, &QAction::triggered, this, &MagDynDlg::ShowNotesDlg);
	connect(acStructSymIdx, &QAction::triggered, this, &MagDynDlg::CalcSymmetryIndices);
	connect(acStructView, &QAction::triggered, this, &MagDynDlg::ShowStructPlotDlg);
	connect(acGroundState, &QAction::triggered, this, &MagDynDlg::ShowGroundStateDlg);
	connect(acTopo, &QAction::triggered, this, &MagDynDlg::ShowTopologyDlg);
	connect(acStructImport, &QAction::triggered, this, &MagDynDlg::ShowTableImporter);
	connect(acStructExportSun, &QAction::triggered,
		this, static_cast<void (MagDynDlg::*)()>(&MagDynDlg::ExportToSunny));
	connect(acStructExportSW, &QAction::triggered,
		this, static_cast<void (MagDynDlg::*)()>(&MagDynDlg::ExportToSpinW));
	connect(m_use_dmi, &QAction::toggled, calc_all);
	if(m_allow_general_J)
		connect(m_use_genJ, &QAction::toggled, calc_all);
	connect(m_use_field, &QAction::toggled, calc_all);
	connect(m_use_temperature, &QAction::toggled, calc_all);
	connect(m_use_formfact, &QAction::toggled, calc_all);
	connect(m_use_weights, &QAction::toggled, calc_all_dyn);
	connect(m_use_projector, &QAction::toggled, calc_all_dyn);
	connect(m_unite_degeneracies, &QAction::toggled, calc_all_dyn);
	connect(m_ignore_annihilation, &QAction::toggled, calc_all_dyn);
	connect(m_force_incommensurate, &QAction::toggled, calc_all_dyn);
	connect(m_autocalc, &QAction::toggled, [this](bool checked)
	{
		if(checked)
			this->CalcAll();
	});

	connect(m_plot_channels, &QAction::toggled, [this](bool checked)
	{
		m_menuChannels->setEnabled(checked);
		this->PlotDispersion();
	});

	for(int i = 0; i < 3; ++i)
	{
		connect(m_hamiltonian_comp[i], &QAction::toggled, calc_all_dyn);

		connect(m_plot_channel[i], &QAction::toggled, [this](bool)
		{
			this->PlotDispersion();
		});
	}

	connect(m_plot_weights_pointsize, &QAction::toggled, [this](bool)
	{
		this->PlotDispersion();
	});
	connect(m_plot_weights_alpha, &QAction::toggled, [this](bool)
	{
		this->PlotDispersion();
	});

	connect(acCalc, &QAction::triggered, [this]()
	{
		this->CalcAll();
	});

	// show trafo dialog
	connect(acTrafoCalc, &QAction::triggered, [this]()
	{
		if(!m_trafos)
			m_trafos = new TrafoCalculator(this, m_sett);

		m_trafos->show();
		m_trafos->raise();
		m_trafos->activateWindow();
	});

	// show preferences dialog
	connect(acPreferences, &QAction::triggered, [this]()
	{
		if(!m_settings_dlg)
		{
			using t_SettingsDlg = SettingsDlg<g_settingsvariables.size(), &g_settingsvariables>;
			m_settings_dlg = new t_SettingsDlg(this, m_sett);

			dynamic_cast<t_SettingsDlg*>(m_settings_dlg)->AddChangedSettingsSlot([this]()
			{
				MagDynDlg::InitSettings();
			});
		}

		m_settings_dlg->show();
		m_settings_dlg->raise();
		m_settings_dlg->activateWindow();
	});

	// show info dialog
	connect(acHelp, &QAction::triggered, [this]()
	{
		QUrl url("https://github.com/ILLGrenoble/takin/wiki/Modelling-Magnetic-Structures");
		if(!QDesktopServices::openUrl(url))
			ShowError("Could not open the wiki.");
	});

	connect(acAboutQt, &QAction::triggered, []()
	{
		qApp->aboutQt();
	});

	// show info dialog
	connect(acAbout, &QAction::triggered, this, &MagDynDlg::ShowInfoDlg);

	// menu bar
	m_menu->addMenu(menuFile);
	m_menu->addMenu(menuStruct);
	m_menu->addMenu(m_menuDisp);
	m_menu->addMenu(menuCalc);
	m_menu->addMenu(menuTools);
	m_menu->addMenu(menuHelp);
	m_maingrid->setMenuBar(m_menu);
}



void MagDynDlg::mousePressEvent(QMouseEvent *evt)
{
	QDialog::mousePressEvent(evt);
}



/**
 * dialog is closing
 */
void MagDynDlg::closeEvent(QCloseEvent *)
{
	if(!m_sett)
		return;

	m_recent.TrimEntries();
	m_sett->setValue("recent_files", m_recent.GetRecentFiles());

	m_recent_struct.TrimEntries();
	m_sett->setValue("recent_struct_files", m_recent_struct.GetRecentFiles());

	m_sett->setValue("geo", saveGeometry());

	if(m_split_inout)
		m_sett->setValue("splitter", m_split_inout->saveState());
}



/**
 * a file is being dragged over the window
 */
void MagDynDlg::dragEnterEvent(QDragEnterEvent *evt)
{
	if(evt)
		evt->accept();
}



/**
 * a file is being dropped onto the window
 */
void MagDynDlg::dropEvent(QDropEvent *evt)
{
	const QMimeData *mime = evt->mimeData();
	if(!mime)
		return;

	for(const QUrl& url : mime->urls())
	{
		if(!url.isLocalFile())
			continue;

		Load(url.toLocalFile());
		evt->accept();
		break;
	}
}



/**
 * refresh and calculate everything
 */
void MagDynDlg::CalcAll()
{
	// calculate structure
	SyncToKernel();
	if(m_structplot_dlg)
		m_structplot_dlg->Sync();

	// calculate dynamics
	CalcDispersion();
	CalcHamiltonian();
}



/**
 * enable (or disable) GUI inputs after calculation threads have finished
 */
void MagDynDlg::EnableInput(bool enable)
{
	m_startEnabled = enable;

	if(enable)
	{
		m_tabs_in->setEnabled(true);
		m_tabs_out->setEnabled(true);
		m_menu->setEnabled(true);

		m_btnStartStop->setText("Calculate");
		m_btnStartStop->setToolTip("Start calculation.");
		m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_menu->setEnabled(false);
		m_tabs_out->setEnabled(false);
		m_tabs_in->setEnabled(false);

		m_btnStartStop->setText("Stop");
		m_btnStartStop->setToolTip("Stop calculation.");
		m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
}



void MagDynDlg::ShowError(const char* msg, bool critical) const
{
	MagDynDlg *dlg = const_cast<MagDynDlg*>(this);

	if(critical)
		QMessageBox::critical(dlg, windowTitle() + " -- Error", msg);
	else
		QMessageBox::warning(dlg, windowTitle() + " -- Warning", msg);
}
