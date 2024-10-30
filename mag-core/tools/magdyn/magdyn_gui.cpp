/**
 * magnetic dynamics -- gui setup
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

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>

#include <QtGui/QDesktopServices>


// instantiate settings dialog
template class SettingsDlg<g_settingsvariables.size(), &g_settingsvariables>;
using t_SettingsDlg = SettingsDlg<g_settingsvariables.size(), &g_settingsvariables>;



/**
 * initialise the static part of the settings dialog
 */
void MagDynDlg::InitSettingsDlg()
{
	// set-up common gui settings variables
	t_SettingsDlg::SetGuiTheme(&g_theme);
	t_SettingsDlg::SetGuiFont(&g_font);
	t_SettingsDlg::SetGuiUseNativeMenubar(&g_use_native_menubar);
	t_SettingsDlg::SetGuiUseNativeDialogs(&g_use_native_dialogs);

	// restore settings
	t_SettingsDlg::ReadSettings(m_sett);
}



/**
 * get changes from settings dialog
 */
void MagDynDlg::InitSettings()
{
	// calculator settings
	m_dyn.SetSilent(g_silent);
	m_dyn.SetEpsilon(g_eps);
	m_dyn.SetPrecision(g_prec);
	m_dyn.SetBoseCutoffEnergy(g_bose_cutoff);
	m_dyn.SetCholeskyMaxTries(g_cholesky_maxtries);
	m_dyn.SetCholeskyInc(g_cholesky_delta);

	m_recent.SetMaxRecentFiles(g_maxnum_recents);
	m_recent_struct.SetMaxRecentFiles(g_maxnum_recents);

	if(g_font != "")
	{
		QFont font = this->font();
		if(font.fromString(g_font))
			setFont(font);
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
 * allows the user to specify magnetic sites
 */
void MagDynDlg::CreateSitesPanel()
{
	m_sitespanel = new QWidget(this);

	m_sitestab = new QTableWidget(m_sitespanel);
	m_sitestab->setShowGrid(true);
	m_sitestab->setAlternatingRowColors(true);
	m_sitestab->setSortingEnabled(true);
	m_sitestab->setMouseTracking(true);
	m_sitestab->setSelectionBehavior(QTableWidget::SelectRows);
	m_sitestab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_sitestab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_sitestab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing() + 4);
	m_sitestab->verticalHeader()->setVisible(true);

	m_sitestab->setColumnCount(NUM_SITE_COLS);

	m_sitestab->setHorizontalHeaderItem(COL_SITE_NAME,
		new QTableWidgetItem{"Name"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_POS_X,
		new QTableWidgetItem{"x"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_POS_Y,
		new QTableWidgetItem{"y"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_POS_Z,
		new QTableWidgetItem{"z"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SYM_IDX,
		new QTableWidgetItem{"Sym. Idx."});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_X,
		new QTableWidgetItem{"Spin x"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_Y,
		new QTableWidgetItem{"Spin y"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_Z,
		new QTableWidgetItem{"Spin z"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_MAG,
		new QTableWidgetItem{"Spin |S|"});
	m_sitestab->setHorizontalHeaderItem(COL_SITE_RGB,
		new QTableWidgetItem{"Colour"});

	if(m_allow_ortho_spin)
	{
		m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_ORTHO_X,
			new QTableWidgetItem{"Spin ux"});
		m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_ORTHO_Y,
			new QTableWidgetItem{"Spin uy"});
		m_sitestab->setHorizontalHeaderItem(COL_SITE_SPIN_ORTHO_Z,
			new QTableWidgetItem{"Spin uz"});
	}
	else
	{
		m_sitestab->setColumnCount(NUM_SITE_COLS - 3);
	}

	m_sitestab->setColumnWidth(COL_SITE_NAME, 90);
	m_sitestab->setColumnWidth(COL_SITE_POS_X, 80);
	m_sitestab->setColumnWidth(COL_SITE_POS_Y, 80);
	m_sitestab->setColumnWidth(COL_SITE_POS_Z, 80);
	m_sitestab->setColumnWidth(COL_SITE_SYM_IDX, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_X, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_Y, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_Z, 80);
	m_sitestab->setColumnWidth(COL_SITE_SPIN_MAG, 80);
	m_sitestab->setColumnWidth(COL_SITE_RGB, 80);

	if(m_allow_ortho_spin)
	{
		m_sitestab->setColumnWidth(COL_SITE_SPIN_ORTHO_X, 80);
		m_sitestab->setColumnWidth(COL_SITE_SPIN_ORTHO_Y, 80);
		m_sitestab->setColumnWidth(COL_SITE_SPIN_ORTHO_Z, 80);
	}

	m_sitestab->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAdd = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_sitespanel);
	QPushButton *btnDel = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_sitespanel);
	QPushButton *btnUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_sitespanel);
	QPushButton *btnDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_sitespanel);

	btnAdd->setToolTip("Add a site.");
	btnDel->setToolTip("Delete selected site(s).");
	btnUp->setToolTip("Move selected site(s) up.");
	btnDown->setToolTip("Move selected site(s) down.");

	QPushButton *btnMirrorAtoms = new QPushButton("Mirror", m_sitespanel);
	QPushButton *btnShowNotes = new QPushButton(
		QIcon::fromTheme("accessories-text-editor"),
		"Notes...", m_sitespanel);
	QPushButton *btnGroundState = new QPushButton("Ground State...", m_sitespanel);
	btnMirrorAtoms->setToolTip("Flip the coordinates of the sites.");
	btnShowNotes->setToolTip("Add notes or comments describing the magnetic structure.");
	btnGroundState->setToolTip("Minimise ground state energy.");
#ifndef __TLIBS2_MAGDYN_USE_MINUIT__
	btnGroundState->setEnabled(false);
#endif

	// extend cell
	const char* idx_names[] = {"x = ", "y = ", "z = "};
	for(int cell_idx = 0; cell_idx < 3; ++cell_idx)
	{
		m_extCell[cell_idx] = new QSpinBox(m_sitespanel);
		m_extCell[cell_idx]->setMinimum(1);
		m_extCell[cell_idx]->setMaximum(99);
		m_extCell[cell_idx]->setValue(cell_idx == 2 ? 2 : 1);
		m_extCell[cell_idx]->setPrefix(idx_names[cell_idx]);
		m_extCell[cell_idx]->setToolTip("Order of supercell.");
		m_extCell[cell_idx]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	QPushButton *btnExtCell = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_sitespanel);
	btnExtCell->setToolTip("Extend the unit cell.");

	QPushButton *btnGenBySG = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_sitespanel);
	btnGenBySG->setToolTip("Create site positions from space group"
		" symmetry operators and existing positions.");

	for(QPushButton *btn : { btnAdd, btnDel, btnUp, btnDown, btnGenBySG, btnExtCell })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}


	auto grid = new QGridLayout(m_sitespanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_sitestab, y, 0, 1, 4);
	grid->addWidget(btnAdd, ++y, 0, 1, 1);
	grid->addWidget(btnDel, y, 1, 1, 1);
	grid->addWidget(btnUp, y, 2, 1, 1);
	grid->addWidget(btnDown, y++, 3, 1, 1);
	grid->addWidget(btnMirrorAtoms, y, 0, 1, 1);
	grid->addWidget(btnShowNotes, y, 2, 1, 1);
	grid->addWidget(btnGroundState, y++, 3, 1, 1);

	auto sep1 = new QFrame(m_sitespanel);
	sep1->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1, 1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1, 1);

	grid->addWidget(new QLabel("Extend Structure, Copying Existing Sites:",
		m_sitespanel), y++, 0, 1, 4);
	grid->addWidget(m_extCell[0], y, 0, 1, 1);
	grid->addWidget(m_extCell[1], y, 1, 1, 1);
	grid->addWidget(m_extCell[2], y, 2, 1, 1);
	grid->addWidget(btnExtCell, y++, 3, 1, 1);

	grid->addWidget(new QLabel("Create Symmetry-Equivalent Sites:",
		m_sitespanel), y, 0, 1, 3);
	grid->addWidget(btnGenBySG, y++, 3, 1, 1);

	// table CustomContextMenu
	QMenu *menuTableContext = new QMenu(m_sitestab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Site Before", this,
		[this]() { this->AddSiteTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Site After", this,
		[this]() { this->AddSiteTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Site", this,
		[this]() { this->AddSiteTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Site",this,
		[this]() { this->DelTabItem(m_sitestab); });


	// table CustomContextMenu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_sitestab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Site", this,
		[this]() { this->AddSiteTabItem(); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Site", this,
		[this]() { this->DelTabItem(m_sitestab); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddSiteTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_sitestab); });
	connect(btnUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_sitestab); });
	connect(btnDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_sitestab); });
	connect(btnGenBySG, &QAbstractButton::clicked,
		this, &MagDynDlg::GenerateSitesFromSG);
	connect(btnExtCell, &QAbstractButton::clicked,
		this, &MagDynDlg::ExtendStructure);

	connect(btnMirrorAtoms, &QAbstractButton::clicked, this, &MagDynDlg::MirrorAtoms);
	connect(btnShowNotes, &QAbstractButton::clicked, this, &MagDynDlg::ShowNotesDlg);
	connect(btnGroundState, &QAbstractButton::clicked, this, &MagDynDlg::ShowGroundStateDlg);

	connect(m_sitestab, &QTableWidget::itemSelectionChanged, [this]()
	{
		QList<QTableWidgetItem*> selected = m_sitestab->selectedItems();
		if(selected.size() == 0)
			return;

		const QTableWidgetItem* item = *selected.begin();
		m_sites_cursor_row = item->row();
		if(m_sites_cursor_row < 0 ||
			std::size_t(m_sites_cursor_row) >= m_dyn.GetMagneticSitesCount())
		{
			m_status->setText("");
			return;
		}

		const auto* site = GetSiteFromTableIndex(m_sites_cursor_row);
		if(!site)
		{
			m_status->setText("Invalid site selected.");
			return;
		}

		std::ostringstream ostr;
		ostr.precision(g_prec_gui);
		ostr << "Site " << site->name << ".";
		m_status->setText(ostr.str().c_str());
	});
	connect(m_sitestab, &QTableWidget::itemChanged,
		this, &MagDynDlg::SitesTableItemChanged);
	connect(m_sitestab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_sitestab, menuTableContext, menuTableContextNoItem, pt);
	});


	m_tabs_in->addTab(m_sitespanel, "Sites");
}



/**
 * allows the user to specify magnetic couplings between sites
 */
void MagDynDlg::CreateExchangeTermsPanel()
{
	m_termspanel = new QWidget(this);

	m_termstab = new QTableWidget(m_termspanel);
	m_termstab->setShowGrid(true);
	m_termstab->setAlternatingRowColors(true);
	m_termstab->setSortingEnabled(true);
	m_termstab->setMouseTracking(true);
	m_termstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_termstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_termstab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_termstab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing() + 4);
	m_termstab->verticalHeader()->setVisible(true);

	m_termstab->setColumnCount(NUM_XCH_COLS);
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_NAME, new QTableWidgetItem{"Name"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_ATOM1_IDX, new QTableWidgetItem{"Site 1"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_ATOM2_IDX, new QTableWidgetItem{"Site 2"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_DIST_X, new QTableWidgetItem{"Cell Δx"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_DIST_Y, new QTableWidgetItem{"Cell Δy"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_DIST_Z, new QTableWidgetItem{"Cell Δz"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_SYM_IDX, new QTableWidgetItem{"Sym. Idx."});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_INTERACTION, new QTableWidgetItem{"Exch. J"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_DMI_X, new QTableWidgetItem{"DMI x"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_DMI_Y, new QTableWidgetItem{"DMI y"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_DMI_Z, new QTableWidgetItem{"DMI z"});
	m_termstab->setHorizontalHeaderItem(
		COL_XCH_RGB, new QTableWidgetItem{"Colour"});

	if(m_allow_general_J)
	{
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_XX, new QTableWidgetItem{"J xx"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_XY, new QTableWidgetItem{"J xy"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_XZ, new QTableWidgetItem{"J xz"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_YX, new QTableWidgetItem{"J yx"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_YY, new QTableWidgetItem{"J yy"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_YZ, new QTableWidgetItem{"J yz"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_ZX, new QTableWidgetItem{"J zx"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_ZY, new QTableWidgetItem{"J zy"});
		m_termstab->setHorizontalHeaderItem(
			COL_XCH_GEN_ZZ, new QTableWidgetItem{"J zz"});
	}
	else
	{
		m_termstab->setColumnCount(NUM_XCH_COLS - 9);
	}

	m_termstab->setColumnWidth(COL_XCH_NAME, 90);
	m_termstab->setColumnWidth(COL_XCH_ATOM1_IDX, 80);
	m_termstab->setColumnWidth(COL_XCH_ATOM2_IDX, 80);
	m_termstab->setColumnWidth(COL_XCH_DIST_X, 80);
	m_termstab->setColumnWidth(COL_XCH_DIST_Y, 80);
	m_termstab->setColumnWidth(COL_XCH_DIST_Z, 80);
	m_termstab->setColumnWidth(COL_XCH_SYM_IDX, 80);
	m_termstab->setColumnWidth(COL_XCH_INTERACTION, 80);
	m_termstab->setColumnWidth(COL_XCH_DMI_X, 80);
	m_termstab->setColumnWidth(COL_XCH_DMI_Y, 80);
	m_termstab->setColumnWidth(COL_XCH_DMI_Z, 80);
	m_termstab->setColumnWidth(COL_XCH_RGB, 80);

	if(m_allow_general_J)
	{
		m_termstab->setColumnWidth(COL_XCH_GEN_XX, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_XY, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_XZ, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_YX, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_YY, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_YZ, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_ZX, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_ZY, 80);
		m_termstab->setColumnWidth(COL_XCH_GEN_ZZ, 80);
	}

	m_termstab->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAdd = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_termspanel);
	QPushButton *btnDel = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_termspanel);
	QPushButton *btnUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_termspanel);
	QPushButton *btnDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_termspanel);

	btnAdd->setToolTip("Add a coupling between two sites.");
	btnDel->setToolTip("Delete selected coupling(s).");
	btnUp->setToolTip("Move selected coupling(s) up.");
	btnDown->setToolTip("Move selected coupling(s) down.");

	for(QPushButton *btn : { btnAdd, btnDel, btnUp, btnDown })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}


	// couplings from distances
	m_maxdist = new QDoubleSpinBox(m_termspanel);
	m_maxdist->setDecimals(3);
	m_maxdist->setMinimum(0.001);
	m_maxdist->setMaximum(99.999);
	m_maxdist->setSingleStep(0.1);
	m_maxdist->setValue(5);
	m_maxdist->setPrefix("d = ");
	m_maxdist->setToolTip("Maximum distance between sites.");
	m_maxdist->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	m_maxSC = new QSpinBox(m_termspanel);
	m_maxSC->setMinimum(1);
	m_maxSC->setMaximum(99);
	m_maxSC->setValue(4);
	m_maxSC->setPrefix("order = ");
	m_maxSC->setToolTip("Maximum order of supercell to consider.");
	m_maxSC->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	m_maxcouplings = new QSpinBox(m_termspanel);
	m_maxcouplings->setMinimum(-1);
	m_maxcouplings->setMaximum(999);
	m_maxcouplings->setValue(100);
	m_maxcouplings->setPrefix("n = ");
	m_maxcouplings->setToolTip("Maximum number of couplings to generate (-1: no limit).");
	m_maxcouplings->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	QPushButton *btnGenByDist = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_termspanel);
	btnGenByDist->setToolTip("Create possible couplings by distances between sites.");
	btnGenByDist->setFocusPolicy(Qt::StrongFocus);
	btnGenByDist->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});


	// couplings from space group
	QPushButton *btnGenBySG = new QPushButton(
		QIcon::fromTheme("insert-object"),
		"Generate", m_termspanel);
	btnGenBySG->setToolTip("Create couplings from space group"
		" symmetry operators and existing couplings.");

	btnGenBySG->setFocusPolicy(Qt::StrongFocus);
	btnGenBySG->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	// ordering vector
	m_ordering[0] = new QDoubleSpinBox(m_termspanel);
	m_ordering[1] = new QDoubleSpinBox(m_termspanel);
	m_ordering[2] = new QDoubleSpinBox(m_termspanel);

	// normal axis
	m_normaxis[0] = new QDoubleSpinBox(m_termspanel);
	m_normaxis[1] = new QDoubleSpinBox(m_termspanel);
	m_normaxis[2] = new QDoubleSpinBox(m_termspanel);

	for(int i = 0; i < 3; ++i)
	{
		m_ordering[i]->setDecimals(4);
		m_ordering[i]->setMinimum(-9.9999);
		m_ordering[i]->setMaximum(+9.9999);
		m_ordering[i]->setSingleStep(0.01);
		m_ordering[i]->setValue(0.);
		m_ordering[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});

		m_normaxis[i]->setDecimals(4);
		m_normaxis[i]->setMinimum(-9.9999);
		m_normaxis[i]->setMaximum(+9.9999);
		m_normaxis[i]->setSingleStep(0.01);
		m_normaxis[i]->setValue(i==0 ? 1. : 0.);
		m_normaxis[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	m_ordering[0]->setPrefix("Oh = ");
	m_ordering[1]->setPrefix("Ok = ");
	m_ordering[2]->setPrefix("Ol = ");

	m_normaxis[0]->setPrefix("Nh = ");
	m_normaxis[1]->setPrefix("Nk = ");
	m_normaxis[2]->setPrefix("Nl = ");


	// grid
	auto grid = new QGridLayout(m_termspanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_termstab, y++,0,1,4);
	grid->addWidget(btnAdd, y,0,1,1);
	grid->addWidget(btnDel, y,1,1,1);
	grid->addWidget(btnUp, y,2,1,1);
	grid->addWidget(btnDown, y++,3,1,1);

	auto sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);
	auto sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Generate Possible Coupling Terms"
		" By Distance (\xe2\x84\xab):", m_termspanel), y++,0,1,4);
	grid->addWidget(m_maxdist, y,0,1,1);
	grid->addWidget(m_maxSC, y,1,1,1);
	grid->addWidget(m_maxcouplings, y,2,1,1);
	grid->addWidget(btnGenByDist, y++,3,1,1);
	grid->addWidget(new QLabel("Create Symmetry-Equivalent Couplings:",
		m_termspanel), y,0,1,3);
	grid->addWidget(btnGenBySG, y++,3,1,1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep2, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Ordering Vector:", m_termspanel), y,0,1,1);
	grid->addWidget(m_ordering[0], y,1,1,1);
	grid->addWidget(m_ordering[1], y,2,1,1);
	grid->addWidget(m_ordering[2], y++,3,1,1);
	grid->addWidget(new QLabel("Rotation Axis:", m_termspanel), y,0,1,1);
	grid->addWidget(m_normaxis[0], y,1,1,1);
	grid->addWidget(m_normaxis[1], y,2,1,1);
	grid->addWidget(m_normaxis[2], y++,3,1,1);

	// table CustomContextMenu
	QMenu *menuTableContext = new QMenu(m_termstab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Term Before", this,
		[this]() { this->AddTermTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Term After", this,
		[this]() { this->AddTermTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Term", this,
		[this]() { this->AddTermTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Term", this,
		[this]() { this->DelTabItem(m_termstab); });


	// table CustomContextMenu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_termstab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Term", this,
		[this]() { this->AddTermTabItem(); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Term", this,
		[this]() { this->DelTabItem(m_termstab); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddTermTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_termstab); });
	connect(btnUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_termstab); });
	connect(btnDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_termstab); });
	connect(btnGenByDist, &QAbstractButton::clicked,
		this, &MagDynDlg::GeneratePossibleCouplings);
	connect(btnGenBySG, &QAbstractButton::clicked,
		this, &MagDynDlg::GenerateCouplingsFromSG);

	connect(m_termstab, &QTableWidget::itemSelectionChanged, [this]()
	{
		QList<QTableWidgetItem*> selected = m_termstab->selectedItems();
		if(selected.size() == 0)
			return;

		const QTableWidgetItem* item = *selected.begin();
		m_terms_cursor_row = item->row();
		if(m_terms_cursor_row < 0 ||
			std::size_t(m_terms_cursor_row) >= m_dyn.GetExchangeTermsCount())
		{
			m_status->setText("");
			return;
		}

		const auto* term = GetTermFromTableIndex(m_terms_cursor_row);
		if(!term)
		{
			m_status->setText("Invalid coupling selected.");
			return;
		}

		std::ostringstream ostr;
		ostr.precision(g_prec_gui);
		ostr << "Coupling " << term->name
			<< ": length = " << term->length_calc << " \xe2\x84\xab.";
		m_status->setText(ostr.str().c_str());
	});
	connect(m_termstab, &QTableWidget::itemChanged,
		this, &MagDynDlg::TermsTableItemChanged);
	connect(m_termstab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_termstab, menuTableContext, menuTableContextNoItem, pt);
	});

	auto calc_all = [this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcAll();
	};

	for(int i = 0; i < 3; ++i)
	{
		connect(m_ordering[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);

		connect(m_normaxis[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);
	}


	m_tabs_in->addTab(m_termspanel, "Couplings");
}



/**
 * allows the user to specify the sample properties
 */
void MagDynDlg::CreateSamplePanel()
{
	m_samplepanel = new QWidget(this);

	// crystal lattice and angles
	static const char* latticestr[] = { "a = ", "b = ", "c = " };
	for(int i = 0; i < 3; ++i)
	{
		m_xtallattice[i] = new QDoubleSpinBox(m_samplepanel);
		m_xtallattice[i]->setDecimals(3);
		m_xtallattice[i]->setMinimum(0.001);
		m_xtallattice[i]->setMaximum(99.999);
		m_xtallattice[i]->setSingleStep(0.1);
		m_xtallattice[i]->setValue(5);
		m_xtallattice[i]->setPrefix(latticestr[i]);
		//m_xtallattice[i]->setSuffix(" \xe2\x84\xab");
		m_xtallattice[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	const char* anlesstr[] = { "α = ", "β = ", "γ = " };
	for(int i = 0; i < 3; ++i)
	{
		m_xtalangles[i] = new QDoubleSpinBox(m_samplepanel);
		m_xtalangles[i]->setDecimals(2);
		m_xtalangles[i]->setMinimum(0.01);
		m_xtalangles[i]->setMaximum(180.);
		m_xtalangles[i]->setSingleStep(0.1);
		m_xtalangles[i]->setValue(90);
		m_xtalangles[i]->setPrefix(anlesstr[i]);
		//m_xtalangles[i]->setSuffix("\xc2\xb0");
		m_xtalangles[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	// space groups
	m_comboSG = new QComboBox(m_samplepanel);
	m_comboSG->setFocusPolicy(Qt::StrongFocus);

	// scattering plane
	static const char* recipstr[] = { "h = ", "k = ", "l = " };
	for(int i = 0; i < 6; ++i)
	{
		m_scatteringplane[i] = new QDoubleSpinBox(m_samplepanel);
		m_scatteringplane[i]->setDecimals(3);
		m_scatteringplane[i]->setMinimum(-99.999);
		m_scatteringplane[i]->setMaximum(99.999);
		m_scatteringplane[i]->setSingleStep(1.);
		m_scatteringplane[i]->setValue(0);
		m_scatteringplane[i]->setPrefix(recipstr[i % 3]);
		//m_scatteringplane[i]->setSuffix(" rlu");
		m_scatteringplane[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}
	m_scatteringplane[0]->setValue(1);
	m_scatteringplane[4]->setValue(1);

	// form factor
	m_ffact = new QPlainTextEdit(m_samplepanel);


	auto grid = new QGridLayout(m_samplepanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;

	// crystal
	grid->addWidget(new QLabel("Crystal Definition", m_samplepanel), y++,0,1,4);
	grid->addWidget(new QLabel("Lattice (\xe2\x84\xab):", m_samplepanel), y,0,1,1);
	grid->addWidget(m_xtallattice[0], y,1,1,1);
	grid->addWidget(m_xtallattice[1], y,2,1,1);
	grid->addWidget(m_xtallattice[2], y++,3,1,1);
	grid->addWidget(new QLabel("Angles (\xc2\xb0):", m_samplepanel), y,0,1,1);
	grid->addWidget(m_xtalangles[0], y,1,1,1);
	grid->addWidget(m_xtalangles[1], y,2,1,1);
	grid->addWidget(m_xtalangles[2], y++,3,1,1);
	grid->addWidget(new QLabel("Space Group:", m_samplepanel), y,0,1,1);
	grid->addWidget(m_comboSG, y++,1,1,3);

	// separator
	auto sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	// scattering plane
	grid->addWidget(new QLabel("Scattering plane", m_samplepanel), y++,0,1,4);
	grid->addWidget(new QLabel("Vector 1 (rlu):", m_samplepanel), y,0,1,1);
	grid->addWidget(m_scatteringplane[0], y,1,1,1);
	grid->addWidget(m_scatteringplane[1], y,2,1,1);
	grid->addWidget(m_scatteringplane[2], y++,3,1,1);
	grid->addWidget(new QLabel("Vector 2 (rlu):", m_samplepanel), y,0,1,1);
	grid->addWidget(m_scatteringplane[3], y,1,1,1);
	grid->addWidget(m_scatteringplane[4], y,2,1,1);
	grid->addWidget(m_scatteringplane[5], y++,3,1,1);

	// separator
	auto sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep2, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	// form factor formula
	grid->addWidget(new QLabel("Magnetic Form Factor", m_samplepanel), y++,0,1,4);
	grid->addWidget(new QLabel("Enter Formula, f_M(Q) = ", m_samplepanel), y++,0,1,1);
	grid->addWidget(m_ffact, y++,0,1,4);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++,0, 1,1);


	auto calc_all = [this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcAll();
	};

	connect(m_ffact, &QPlainTextEdit::textChanged, calc_all);

	for(int i = 0; i < 2*3; ++i)
	{
		if(i < 3)
		{
			connect(m_xtallattice[i],
				static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
				calc_all);
			connect(m_xtalangles[i],
				static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
				calc_all);
		}

		connect(m_scatteringplane[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);
	}

	m_tabs_in->addTab(m_samplepanel, "Sample");
}



/**
 * lets the user define variables to be used for the J and DMI parameters
 */
void MagDynDlg::CreateVariablesPanel()
{
	m_varspanel = new QWidget(this);

	m_varstab = new QTableWidget(m_varspanel);
	m_varstab->setShowGrid(true);
	m_varstab->setAlternatingRowColors(true);
	m_varstab->setSortingEnabled(true);
	m_varstab->setMouseTracking(true);
	m_varstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_varstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_varstab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_varstab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing() + 4);
	m_varstab->verticalHeader()->setVisible(true);

	m_varstab->setColumnCount(NUM_VARS_COLS);
	m_varstab->setHorizontalHeaderItem(
		COL_VARS_NAME, new QTableWidgetItem{"Name"});
	m_varstab->setHorizontalHeaderItem(
		COL_VARS_VALUE_REAL, new QTableWidgetItem{"Value (Re)"});
	m_varstab->setHorizontalHeaderItem(
		COL_VARS_VALUE_IMAG, new QTableWidgetItem{"Value (Im)"});

	m_varstab->setColumnWidth(COL_VARS_NAME, 150);
	m_varstab->setColumnWidth(COL_VARS_VALUE_REAL, 150);
	m_varstab->setColumnWidth(COL_VARS_VALUE_IMAG, 150);
	m_varstab->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAdd = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_varspanel);
	QPushButton *btnDel = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_varspanel);
	QPushButton *btnUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_varspanel);
	QPushButton *btnDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_varspanel);
	QPushButton *btnReplace = new QPushButton(
		QIcon::fromTheme("edit-find-replace"),
		"Replace Values", m_varspanel);

	btnAdd->setToolTip("Add a variable.");
	btnDel->setToolTip("Delete selected variables(s).");
	btnUp->setToolTip("Move selected variable(s) up.");
	btnDown->setToolTip("Move selected variable(s) down.");
	btnReplace->setToolTip("Replace numeric values with variable names.");

	for(QPushButton *btn : { btnAdd, btnDel, btnUp, btnDown, btnReplace })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}


	// grid
	auto grid = new QGridLayout(m_varspanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_varstab, y++,0,1,4);
	grid->addWidget(btnAdd, y,0,1,1);
	grid->addWidget(btnDel, y,1,1,1);
	grid->addWidget(btnUp, y,2,1,1);
	grid->addWidget(btnDown, y++,3,1,1);
	grid->addWidget(btnReplace, y++,0,1,1);


	// table CustomContextMenu
	QMenu *menuTableContext = new QMenu(m_varstab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Variable Before", this,
		[this]() { this->AddVariableTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Variable After", this,
		[this]() { this->AddVariableTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Variable", this,
		[this]() { this->AddVariableTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Variable", this,
		[this]() { this->DelTabItem(m_varstab); });


	// table CustomContextMenu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_varstab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Variable", this,
		[this]() { this->AddVariableTabItem(); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Variable", this,
		[this]() { this->DelTabItem(m_varstab); });


	// signals
	connect(btnAdd, &QAbstractButton::clicked,
		[this]() { this->AddVariableTabItem(-1); });
	connect(btnDel, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_varstab); });
	connect(btnUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_varstab); });
	connect(btnDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_varstab); });
	connect(btnReplace, &QAbstractButton::clicked,
		this, &MagDynDlg::ReplaceValuesWithVariables);

	connect(m_varstab, &QTableWidget::itemSelectionChanged, [this]()
	{
		QList<QTableWidgetItem*> selected = m_varstab->selectedItems();
		if(selected.size() == 0)
			return;

		const QTableWidgetItem* item = *selected.begin();
		m_variables_cursor_row = item->row();
	});
	connect(m_varstab, &QTableWidget::itemChanged,
		this, &MagDynDlg::VariablesTableItemChanged);
	connect(m_varstab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_varstab, menuTableContext, menuTableContextNoItem, pt);
	});


	m_tabs_in->addTab(m_varspanel, "Variables");
}



/**
 * input for sample environment parameters (field, temperature)
 */
void MagDynDlg::CreateSampleEnvPanel()
{
	m_sampleenviropanel = new QWidget(this);

	// field magnitude
	m_field_mag = new QDoubleSpinBox(m_sampleenviropanel);
	m_field_mag->setDecimals(3);
	m_field_mag->setMinimum(0);
	m_field_mag->setMaximum(+99.999);
	m_field_mag->setSingleStep(0.1);
	m_field_mag->setValue(0.);
	m_field_mag->setPrefix("|B| = ");
	m_field_mag->setSuffix(" T");
	m_field_mag->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	// field direction
	m_field_dir[0] = new QDoubleSpinBox(m_sampleenviropanel);
	m_field_dir[1] = new QDoubleSpinBox(m_sampleenviropanel);
	m_field_dir[2] = new QDoubleSpinBox(m_sampleenviropanel);

	// align spins along field (field-polarised state)
	m_align_spins = new QCheckBox(
		"Align Spins Along Field Direction", m_sampleenviropanel);
	m_align_spins->setChecked(false);
	m_align_spins->setFocusPolicy(Qt::StrongFocus);

	// rotation axis
	m_rot_axis[0] = new QDoubleSpinBox(m_sampleenviropanel);
	m_rot_axis[1] = new QDoubleSpinBox(m_sampleenviropanel);
	m_rot_axis[2] = new QDoubleSpinBox(m_sampleenviropanel);

	// rotation angle
	m_rot_angle = new QDoubleSpinBox(m_sampleenviropanel);
	m_rot_angle->setDecimals(3);
	m_rot_angle->setMinimum(-360);
	m_rot_angle->setMaximum(+360);
	m_rot_angle->setSingleStep(0.1);
	m_rot_angle->setValue(90.);
	//m_rot_angle->setSuffix("\xc2\xb0");
	m_rot_angle->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	QPushButton *btn_rotate_ccw = new QPushButton(
		QIcon::fromTheme("object-rotate-left"),
		"Rotate CCW", m_sampleenviropanel);
	QPushButton *btn_rotate_cw = new QPushButton(
		QIcon::fromTheme("object-rotate-right"),
		"Rotate CW", m_sampleenviropanel);
	btn_rotate_ccw->setToolTip("Rotate the magnetic field in the counter-clockwise direction.");
	btn_rotate_cw->setToolTip("Rotate the magnetic field in the clockwise direction.");
	btn_rotate_ccw->setFocusPolicy(Qt::StrongFocus);
	btn_rotate_cw->setFocusPolicy(Qt::StrongFocus);


	// table with saved fields
	m_fieldstab = new QTableWidget(m_sampleenviropanel);
	m_fieldstab->setShowGrid(true);
	m_fieldstab->setAlternatingRowColors(true);
	m_fieldstab->setSortingEnabled(true);
	m_fieldstab->setMouseTracking(true);
	m_fieldstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_fieldstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_fieldstab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_fieldstab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing() + 4);
	m_fieldstab->verticalHeader()->setVisible(true);

	m_fieldstab->setColumnCount(NUM_FIELD_COLS);
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_H,
		new QTableWidgetItem{"Bh"});
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_K,
		new QTableWidgetItem{"Bk"});
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_L,
		new QTableWidgetItem{"Bl"});
	m_fieldstab->setHorizontalHeaderItem(COL_FIELD_MAG,
		new QTableWidgetItem{"|B|"});

	m_fieldstab->setColumnWidth(COL_FIELD_H, 150);
	m_fieldstab->setColumnWidth(COL_FIELD_K, 150);
	m_fieldstab->setColumnWidth(COL_FIELD_L, 150);
	m_fieldstab->setColumnWidth(COL_FIELD_MAG, 150);
	m_fieldstab->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAddField = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_sampleenviropanel);
	QPushButton *btnDelField = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_sampleenviropanel);
	QPushButton *btnFieldUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_sampleenviropanel);
	QPushButton *btnFieldDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_sampleenviropanel);

	btnAddField->setToolTip("Add a magnetic field.");
	btnDelField->setToolTip("Delete selected magnetic field(s).");
	btnFieldUp->setToolTip("Move selected magnetic field(s) up.");
	btnFieldDown->setToolTip("Move selected magnetic field(s) down.");

	QPushButton *btnSetField = new QPushButton("Set Field", m_sampleenviropanel);
	btnSetField->setToolTip("Set the selected field as the currently active one.");

	for(QPushButton *btn : { btnAddField, btnDelField, btnFieldUp, btnFieldDown })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}


	// table CustomContextMenu
	QMenu *menuTableContext = new QMenu(m_fieldstab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Field Before", this,
		[this]() { this->AddFieldTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Field After", this,
		[this]() { this->AddFieldTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Field", this,
		[this]() { this->AddFieldTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Field", this,
		[this]() { this->DelTabItem(m_fieldstab); });
	menuTableContext->addSeparator();
	menuTableContext->addAction(
		QIcon::fromTheme("go-home"),
		"Set As Current Field", this,
		[this]() { this->SetCurrentField(); });


	// table CustomContextMenu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_fieldstab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Field", this,
		[this]() { this->AddFieldTabItem(-1,
			m_field_dir[0]->value(),
			m_field_dir[1]->value(),
			m_field_dir[2]->value(),
			m_field_mag->value()); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Field", this,
		[this]() { this->DelTabItem(m_fieldstab); });


	// temperature
	m_temperature = new QDoubleSpinBox(m_sampleenviropanel);
	m_temperature->setDecimals(2);
	m_temperature->setMinimum(0);
	m_temperature->setMaximum(+999.99);
	m_temperature->setSingleStep(0.1);
	m_temperature->setValue(300.);
	m_temperature->setPrefix("T = ");
	m_temperature->setSuffix(" K");
	m_temperature->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});

	for(int i = 0; i < 3; ++i)
	{
		m_field_dir[i]->setDecimals(4);
		m_field_dir[i]->setMinimum(-99.9999);
		m_field_dir[i]->setMaximum(+99.9999);
		m_field_dir[i]->setSingleStep(0.1);
		m_field_dir[i]->setValue(i == 2 ? 1. : 0.);
		m_field_dir[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});

		m_rot_axis[i]->setDecimals(4);
		m_rot_axis[i]->setMinimum(-99.9999);
		m_rot_axis[i]->setMaximum(+99.9999);
		m_rot_axis[i]->setSingleStep(0.1);
		m_rot_axis[i]->setValue(i == 2 ? 1. : 0.);
		m_rot_axis[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	m_field_dir[0]->setPrefix("Bh = ");
	m_field_dir[1]->setPrefix("Bk = ");
	m_field_dir[2]->setPrefix("Bl = ");

	auto grid = new QGridLayout(m_sampleenviropanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(new QLabel("Magnetic Field:", m_sampleenviropanel), y++,0,1,2);
	grid->addWidget(new QLabel("Magnitude:", m_sampleenviropanel), y,0,1,1);
	grid->addWidget(m_field_mag, y++,1,1,1);
	grid->addWidget(new QLabel("Direction (rlu):", m_sampleenviropanel), y,0,1,1);
	grid->addWidget(m_field_dir[0], y,1,1,1);
	grid->addWidget(m_field_dir[1], y,2,1,1);
	grid->addWidget(m_field_dir[2], y++,3,1,1);
	grid->addWidget(m_align_spins, y++,0,1,4);

	auto sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);
	auto sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);
	auto sep3 = new QFrame(m_sampleenviropanel);
	sep3->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep1, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Rotate Magnetic Field:", m_sampleenviropanel), y++,0,1,2);
	grid->addWidget(new QLabel("Axis (rlu):", m_sampleenviropanel), y,0,1,1);
	grid->addWidget(m_rot_axis[0], y,1,1,1);
	grid->addWidget(m_rot_axis[1], y,2,1,1);
	grid->addWidget(m_rot_axis[2], y++,3,1,1);
	grid->addWidget(new QLabel("Angle (\xc2\xb0):", m_sampleenviropanel), y,0,1,1);
	grid->addWidget(m_rot_angle, y,1,1,1);
	grid->addWidget(btn_rotate_ccw, y,2,1,1);
	grid->addWidget(btn_rotate_cw, y++,3,1,1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep2, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Saved Fields:", m_sampleenviropanel), y++,0,1,4);
	grid->addWidget(m_fieldstab, y,0,1,4);
	grid->addWidget(btnAddField, ++y,0,1,1);
	grid->addWidget(btnDelField, y,1,1,1);
	grid->addWidget(btnFieldUp, y,2,1,1);
	grid->addWidget(btnFieldDown, y++,3,1,1);
	grid->addWidget(btnSetField, y++,3,1,1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep3, y++,0, 1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Temperature:", m_sampleenviropanel), y,0,1,1);
	grid->addWidget(m_temperature, y++,1,1,1);

	auto calc_all = [this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcAll();
	};

	// signals
	connect(m_field_mag,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		calc_all);

	for(int i = 0; i < 3; ++i)
	{
		connect(m_field_dir[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			calc_all);
	}

	connect(m_temperature,
		static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
		calc_all);

	connect(m_align_spins, &QCheckBox::toggled, calc_all);

	connect(btn_rotate_ccw, &QAbstractButton::clicked, [this]()
	{
		RotateField(true);
	});

	connect(btn_rotate_cw, &QAbstractButton::clicked, [this]()
	{
		RotateField(false);
	});

	connect(btnAddField, &QAbstractButton::clicked,
		[this]() { this->AddFieldTabItem(-1,
			m_field_dir[0]->value(),
			m_field_dir[1]->value(),
			m_field_dir[2]->value(),
			m_field_mag->value()); });
	connect(btnDelField, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_fieldstab); });
	connect(btnFieldUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_fieldstab); });
	connect(btnFieldDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_fieldstab); });

	connect(btnSetField, &QAbstractButton::clicked,
		[this]() { this->SetCurrentField(); });

	connect(m_fieldstab, &QTableWidget::itemSelectionChanged, [this]()
	{
		QList<QTableWidgetItem*> selected = m_fieldstab->selectedItems();
		if(selected.size() == 0)
			return;

		const QTableWidgetItem* item = *selected.begin();
		m_fields_cursor_row = item->row();
	});
	connect(m_fieldstab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_fieldstab, menuTableContext, menuTableContextNoItem, pt);
	});


	m_tabs_in->addTab(m_sampleenviropanel, "Environment");
}



/**
 * plots the dispersion relation for a given Q path
 */
void MagDynDlg::CreateDispersionPanel()
{
	const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	m_disppanel = new QWidget(this);

	// plotter
	m_plot = new QCustomPlot(m_disppanel);
	m_plot->xAxis->setLabel("Q (rlu)");
	m_plot->yAxis->setLabel("E (meV)");
	m_plot->setInteraction(QCP::iRangeDrag, true);
	m_plot->setInteraction(QCP::iRangeZoom, true);
	m_plot->setSelectionRectMode(QCP::srmZoom);
	m_plot->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	// start and stop coordinates
	m_Q_start[0] = new QDoubleSpinBox(m_disppanel);
	m_Q_start[1] = new QDoubleSpinBox(m_disppanel);
	m_Q_start[2] = new QDoubleSpinBox(m_disppanel);
	m_Q_end[0] = new QDoubleSpinBox(m_disppanel);
	m_Q_end[1] = new QDoubleSpinBox(m_disppanel);
	m_Q_end[2] = new QDoubleSpinBox(m_disppanel);

	m_Q_start[0]->setToolTip("Dispersion initial momentum transfer, h_i (rlu).");
	m_Q_start[1]->setToolTip("Dispersion initial momentum transfer, k_i (rlu).");
	m_Q_start[2]->setToolTip("Dispersion initial momentum transfer, l_i (rlu).");
	m_Q_end[0]->setToolTip("Dispersion final momentum transfer, h_f (rlu).");
	m_Q_end[1]->setToolTip("Dispersion final momentum transfer, k_f (rlu).");
	m_Q_end[2]->setToolTip("Dispersion final momentum transfer, l_f (rlu).");

	// number of Q points in the plot
	m_num_points = new QSpinBox(m_disppanel);
	m_num_points->setMinimum(1);
	m_num_points->setMaximum(99999);
	m_num_points->setValue(512);
	m_num_points->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_num_points->setToolTip("Number of Q points in the plot.");

	// scaling factor for weights
	for(auto** comp : {&m_weight_scale, &m_weight_min, &m_weight_max})
	{
		*comp = new QDoubleSpinBox(m_disppanel);
		(*comp)->setDecimals(4);
		(*comp)->setMinimum(0.);
		(*comp)->setMaximum(+9999.9999);
		(*comp)->setSingleStep(0.1);
		(*comp)->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	m_weight_scale->setValue(1.);
	m_weight_min->setValue(0.);
	m_weight_max->setValue(9999);
	m_weight_min->setMinimum(-1.);	// -1: disable clamping
	m_weight_max->setMinimum(-1.);	// -1: disable clamping
	m_weight_min->setToolTip("Minimum spectral weight for clamping.");
	m_weight_max->setToolTip("Maximum spectral weight for clamping.");
	m_weight_scale->setToolTip("Spectral weight scaling factor.");

	for(int i = 0; i < 3; ++i)
	{
		m_Q_start[i]->setDecimals(4);
		m_Q_start[i]->setMinimum(-99.9999);
		m_Q_start[i]->setMaximum(+99.9999);
		m_Q_start[i]->setSingleStep(0.01);
		m_Q_start[i]->setValue(0.);
		//m_Q_start[i]->setSuffix(" rlu");
		m_Q_start[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q_start[i]->setPrefix(hklPrefix[i]);

		m_Q_end[i]->setDecimals(4);
		m_Q_end[i]->setMinimum(-99.9999);
		m_Q_end[i]->setMaximum(+99.9999);
		m_Q_end[i]->setSingleStep(0.01);
		m_Q_end[i]->setValue(0.);
		//m_Q_end[i]->setSuffix(" rlu");
		m_Q_end[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q_end[i]->setPrefix(hklPrefix[i]);
	}

	m_Q_start[0]->setValue(-1.);
	m_Q_end[0]->setValue(+1.);

	auto grid = new QGridLayout(m_disppanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_plot, y++,0,1,4);
	grid->addWidget(new QLabel("Start Q (rlu):", m_disppanel), y,0,1,1);
	grid->addWidget(m_Q_start[0], y,1,1,1);
	grid->addWidget(m_Q_start[1], y,2,1,1);
	grid->addWidget(m_Q_start[2], y++,3,1,1);
	grid->addWidget(new QLabel("End Q (rlu):", m_disppanel), y,0,1,1);
	grid->addWidget(m_Q_end[0], y,1,1,1);
	grid->addWidget(m_Q_end[1], y,2,1,1);
	grid->addWidget(m_Q_end[2], y++,3,1,1);
	grid->addWidget(new QLabel("Q Count:", m_disppanel), y,0,1,1);
	grid->addWidget(m_num_points, y,1,1,1);
	grid->addWidget(new QLabel("Weight Scale:", m_disppanel), y,2,1,1);
	grid->addWidget(m_weight_scale, y++,3,1,1);
	grid->addWidget(new QLabel("Min. Weight:", m_disppanel), y,0,1,1);
	grid->addWidget(m_weight_min, y,1,1,1);
	grid->addWidget(new QLabel("Max. Weight:", m_disppanel), y,2,1,1);
	grid->addWidget(m_weight_max, y++,3,1,1);

	// signals
	for(int i = 0; i < 3; ++i)
	{
		connect(m_Q_start[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			if(this->m_autocalc->isChecked())
				this->CalcDispersion();
		});
		connect(m_Q_end[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			if(this->m_autocalc->isChecked())
				this->CalcDispersion();
		});
	}

	connect(m_num_points,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		[this]()
	{
		if(this->m_autocalc->isChecked())
			this->CalcDispersion();
	});

	for(auto* comp : {m_weight_scale, m_weight_min, m_weight_max})
	{
		connect(comp,
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			// update graph weights
			for(GraphWithWeights* graph : m_graphs)
				graph->SetWeightScale(m_weight_scale->value(), m_weight_min->value(), m_weight_max->value());
			if(m_plot)
				m_plot->replot();
		});
	}

	connect(m_plot, &QCustomPlot::mouseMove, this, &MagDynDlg::PlotMouseMove);
	connect(m_plot, &QCustomPlot::mousePress, this, &MagDynDlg::PlotMousePress);


	m_tabs_out->addTab(m_disppanel, "Dispersion");
}



/**
 * shows the hamilton operator for a given Q position
 */
void MagDynDlg::CreateHamiltonPanel()
{
	const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	m_hamiltonianpanel = new QWidget(this);

	// hamiltonian
	m_hamiltonian = new QTextEdit(m_hamiltonianpanel);
	m_hamiltonian->setReadOnly(true);
	m_hamiltonian->setWordWrapMode(QTextOption::NoWrap);
	m_hamiltonian->setLineWrapMode(QTextEdit::NoWrap);
	m_hamiltonian->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	// Q coordinates
	m_Q[0] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[1] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[2] = new QDoubleSpinBox(m_hamiltonianpanel);
	m_Q[0]->setToolTip("Momentum transfer component h (rlu).");
	m_Q[1]->setToolTip("Momentum transfer component k (rlu).");
	m_Q[2]->setToolTip("Momentum transfer component l (rlu).");

	for(int i = 0; i < 3; ++i)
	{
		m_Q[i]->setDecimals(4);
		m_Q[i]->setMinimum(-99.9999);
		m_Q[i]->setMaximum(+99.9999);
		m_Q[i]->setSingleStep(0.01);
		m_Q[i]->setValue(0.);
		m_Q[i]->setSuffix(" rlu");
		m_Q[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q[i]->setPrefix(hklPrefix[i]);
	}

	auto grid = new QGridLayout(m_hamiltonianpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(m_hamiltonian, y++,0,1,4);
	grid->addWidget(new QLabel("Q:", m_hamiltonianpanel), y,0,1,1);
	grid->addWidget(m_Q[0], y,1,1,1);
	grid->addWidget(m_Q[1], y,2,1,1);
	grid->addWidget(m_Q[2], y++,3,1,1);

	// signals
	for(int i = 0; i < 3; ++i)
	{
		connect(m_Q[i],
			static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
			[this]()
		{
			if(this->m_autocalc->isChecked())
				this->CalcHamiltonian();
		});
	}


	m_tabs_out->addTab(m_hamiltonianpanel, "Hamiltonian");
}



/**
 * panel for saved favourite Q positions and paths
 */
void MagDynDlg::CreateCoordinatesPanel()
{
	m_coordinatespanel = new QWidget(this);

	// table with saved fields
	m_coordinatestab = new QTableWidget(m_coordinatespanel);
	m_coordinatestab->setShowGrid(true);
	m_coordinatestab->setAlternatingRowColors(true);
	m_coordinatestab->setSortingEnabled(true);
	m_coordinatestab->setMouseTracking(true);
	m_coordinatestab->setSelectionBehavior(QTableWidget::SelectRows);
	m_coordinatestab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_coordinatestab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_coordinatestab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing() + 4);
	m_coordinatestab->verticalHeader()->setVisible(true);

	m_coordinatestab->setColumnCount(NUM_COORD_COLS);
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_NAME,
		new QTableWidgetItem{"Name"});
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_H,
		new QTableWidgetItem{"h"});
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_K,
		new QTableWidgetItem{"k"});
	m_coordinatestab->setHorizontalHeaderItem(COL_COORD_L,
		new QTableWidgetItem{"l"});

	m_coordinatestab->setColumnWidth(COL_COORD_NAME, 90);
	m_coordinatestab->setColumnWidth(COL_COORD_H, 90);
	m_coordinatestab->setColumnWidth(COL_COORD_K, 90);
	m_coordinatestab->setColumnWidth(COL_COORD_L, 90);
	m_coordinatestab->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	QPushButton *btnAddCoord = new QPushButton(
		QIcon::fromTheme("list-add"),
		"Add", m_coordinatespanel);
	QPushButton *btnDelCoord = new QPushButton(
		QIcon::fromTheme("list-remove"),
		"Delete", m_coordinatespanel);
	QPushButton *btnCoordUp = new QPushButton(
		QIcon::fromTheme("go-up"),
		"Up", m_coordinatespanel);
	QPushButton *btnCoordDown = new QPushButton(
		QIcon::fromTheme("go-down"),
		"Down", m_coordinatespanel);

	btnAddCoord->setToolTip("Add a Q coordinate.");
	btnDelCoord->setToolTip("Delete selected Q coordinate.");
	btnCoordUp->setToolTip("Move selected coordinate(s) up.");
	btnCoordDown->setToolTip("Move selected coordinate(s) down.");

	QPushButton *btnSaveMultiDisp = new QPushButton(
		QIcon::fromTheme("text-x-generic"),
		"Save Data...", m_coordinatespanel);
	QPushButton *btnSaveMultiDispScr = new QPushButton(
		QIcon::fromTheme("text-x-script"),
		"Save Script...", m_coordinatespanel);
	btnSaveMultiDisp->setToolTip("Calculate the dispersion paths and save them to a data file.");
	btnSaveMultiDispScr->setToolTip("Calculate the dispersion paths and save them to a script file.");

	QPushButton *btnSetDispersion = new QPushButton("To Dispersion", m_coordinatespanel);
	QPushButton *btnSetHamilton = new QPushButton("To Hamiltonian", m_coordinatespanel);
	btnSetDispersion->setToolTip("Calculate the dispersion relation for the currently selected Q path.");
	btnSetHamilton->setToolTip("Calculate the Hamiltonian for the currently selected Q coordinate.");

	for(QPushButton *btn : { btnAddCoord, btnDelCoord, btnCoordUp, btnCoordDown })
	{
		btn->setFocusPolicy(Qt::StrongFocus);
		btn->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}


	// table CustomContextMenu
	QMenu *menuTableContext = new QMenu(m_coordinatestab);
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coordinate Before", this,
		[this]() { this->AddCoordinateTabItem(-2); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coordinate After", this,
		[this]() { this->AddCoordinateTabItem(-3); });
	menuTableContext->addAction(
		QIcon::fromTheme("edit-copy"),
		"Clone Coordinate", this,
		[this]() { this->AddCoordinateTabItem(-4); });
	menuTableContext->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Coordinate", this,
		[this]() { this->DelTabItem(m_coordinatestab); });
	menuTableContext->addSeparator();
	menuTableContext->addAction(
		QIcon::fromTheme("go-home"),
		"Calculate Dispersion From This To Next Q", this,
		[this]() { this->SetCurrentCoordinate(0); });
	menuTableContext->addAction(
		QIcon::fromTheme("go-home"),
		"Calculate Hamiltonian For This Q", this,
		[this]() { this->SetCurrentCoordinate(1); });


	// table CustomContextMenu in case nothing is selected
	QMenu *menuTableContextNoItem = new QMenu(m_coordinatestab);
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-add"),
		"Add Coordinate", this,
		[this]() { this->AddCoordinateTabItem(-1, "", 1.,0.,0.); });
	menuTableContextNoItem->addAction(
		QIcon::fromTheme("list-remove"),
		"Delete Coordinate", this,
		[this]() { this->DelTabItem(m_coordinatestab); });


	auto grid = new QGridLayout(m_coordinatespanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(new QLabel("Saved Q Coordinates (rlu):", m_coordinatespanel), y++,0,1,4);
	grid->addWidget(m_coordinatestab, y,0,1,4);
	grid->addWidget(btnAddCoord, ++y,0,1,1);
	grid->addWidget(btnDelCoord, y,1,1,1);
	grid->addWidget(btnCoordUp, y,2,1,1);
	grid->addWidget(btnCoordDown, y++,3,1,1);
	grid->addWidget(btnSetDispersion, y,0,1,1);
	grid->addWidget(btnSetHamilton, y,1,1,1);
	grid->addWidget(btnSaveMultiDisp, y,2,1,1);
	grid->addWidget(btnSaveMultiDispScr, y++,3,1,1);

	// signals
	connect(btnAddCoord, &QAbstractButton::clicked,
		[this]() { this->AddCoordinateTabItem(-1, "", 1.,0.,0.); });
	connect(btnDelCoord, &QAbstractButton::clicked,
		[this]() { this->DelTabItem(m_coordinatestab); });
	connect(btnCoordUp, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemUp(m_coordinatestab); });
	connect(btnCoordDown, &QAbstractButton::clicked,
		[this]() { this->MoveTabItemDown(m_coordinatestab); });

	connect(btnSetDispersion, &QAbstractButton::clicked,
		[this]() { this->SetCurrentCoordinate(0); });
	connect(btnSetHamilton, &QAbstractButton::clicked,
		[this]() { this->SetCurrentCoordinate(1); });
	connect(btnSaveMultiDisp, &QAbstractButton::clicked,
		[this]() { this->SaveMultiDispersion(false); });
	connect(btnSaveMultiDispScr, &QAbstractButton::clicked,
		[this]() { this->SaveMultiDispersion(true); });

	connect(m_coordinatestab, &QTableWidget::itemSelectionChanged, [this]()
	{
		QList<QTableWidgetItem*> selected = m_coordinatestab->selectedItems();
		if(selected.size() == 0)
			return;

		const QTableWidgetItem* item = *selected.begin();
		m_coordinates_cursor_row = item->row();

	});
	connect(m_coordinatestab, &QTableWidget::customContextMenuRequested,
		[this, menuTableContext, menuTableContextNoItem](const QPoint& pt)
	{
		this->ShowTableContextMenu(m_coordinatestab, menuTableContext,
			menuTableContextNoItem, pt);
	});


	m_tabs_out->addTab(m_coordinatespanel, "Coordinates");
}



/**
 * exports data to different file types
 */
void MagDynDlg::CreateExportPanel()
{
	const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	m_exportpanel = new QWidget(this);

	// Q coordinates
	m_exportStartQ[0] = new QDoubleSpinBox(m_exportpanel);
	m_exportStartQ[1] = new QDoubleSpinBox(m_exportpanel);
	m_exportStartQ[2] = new QDoubleSpinBox(m_exportpanel);
	m_exportEndQ[0] = new QDoubleSpinBox(m_exportpanel);
	m_exportEndQ[1] = new QDoubleSpinBox(m_exportpanel);
	m_exportEndQ[2] = new QDoubleSpinBox(m_exportpanel);

	// number of grid points
	for(int i = 0; i < 3; ++i)
	{
		m_exportNumPoints[i] = new QSpinBox(m_exportpanel);
		m_exportNumPoints[i]->setMinimum(1);
		m_exportNumPoints[i]->setMaximum(99999);
		m_exportNumPoints[i]->setValue(128);
		m_exportNumPoints[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
	}

	// export
	m_exportFormat = new QComboBox(m_exportpanel);
	m_exportFormat->addItem("Takin Grid File", EXPORT_GRID);
#ifdef USE_HDF5
	m_exportFormat->addItem("HDF5 File", EXPORT_HDF5);
#endif
	m_exportFormat->addItem("Text File", EXPORT_TEXT);

	QPushButton *btn_export = new QPushButton(
		QIcon::fromTheme("document-save-as"),
		"Export...", m_exportpanel);
	btn_export->setFocusPolicy(Qt::StrongFocus);

	for(int i = 0; i < 3; ++i)
	{
		m_exportStartQ[i]->setDecimals(4);
		m_exportStartQ[i]->setMinimum(-99.9999);
		m_exportStartQ[i]->setMaximum(+99.9999);
		m_exportStartQ[i]->setSingleStep(0.01);
		m_exportStartQ[i]->setValue(-1.);
		m_exportStartQ[i]->setSuffix(" rlu");
		m_exportStartQ[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_exportStartQ[i]->setPrefix(hklPrefix[i]);

		m_exportEndQ[i]->setDecimals(4);
		m_exportEndQ[i]->setMinimum(-99.9999);
		m_exportEndQ[i]->setMaximum(+99.9999);
		m_exportEndQ[i]->setSingleStep(0.01);
		m_exportEndQ[i]->setValue(1.);
		m_exportEndQ[i]->setSuffix(" rlu");
		m_exportEndQ[i]->setSizePolicy(QSizePolicy{
			QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_exportEndQ[i]->setPrefix(hklPrefix[i]);
	}


	auto grid = new QGridLayout(m_exportpanel);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);

	int y = 0;
	grid->addWidget(new QLabel("Export Ranges:", m_exportpanel), y++,0,1,4);
	grid->addWidget(new QLabel("Start Q:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportStartQ[0], y,1,1,1);
	grid->addWidget(m_exportStartQ[1], y,2,1,1);
	grid->addWidget(m_exportStartQ[2], y++,3,1,1);
	grid->addWidget(new QLabel("End Q:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportEndQ[0], y,1,1,1);
	grid->addWidget(m_exportEndQ[1], y,2,1,1);
	grid->addWidget(m_exportEndQ[2], y++,3,1,1);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	auto sep1 = new QFrame(m_sampleenviropanel);
	sep1->setFrameStyle(QFrame::HLine);
	grid->addWidget(sep1, y++, 0,1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addWidget(new QLabel("Number of Grid Points per Q Direction:",
		m_exportpanel), y++,0,1,4);
	grid->addWidget(new QLabel("Points:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportNumPoints[0], y,1,1,1);
	grid->addWidget(m_exportNumPoints[1], y,2,1,1);
	grid->addWidget(m_exportNumPoints[2], y++,3,1,1);

	auto sep2 = new QFrame(m_sampleenviropanel);
	sep2->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep2, y++, 0,1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	QLabel* labelBoseInfo = new QLabel(QString(
		"Info: If this grid file is to be used in Takin's "
		"resolution convolution module (\"Model Source: Uniform Grid\"), "
		"please disable the Bose factor (\"Calculation\" -> \"Use Bose Factor\" [off]). "
		"The Bose factor is already managed by the convolution module."),
		m_exportpanel);
	labelBoseInfo->setWordWrap(true);
	grid->addWidget(labelBoseInfo, y++,0,1,4);

	auto sep3 = new QFrame(m_sampleenviropanel);
	sep3->setFrameStyle(QFrame::HLine);

	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);
	grid->addWidget(sep3, y++, 0,1,4);
	grid->addItem(new QSpacerItem(8, 8,
		QSizePolicy::Minimum, QSizePolicy::Fixed),
		y++,0, 1,1);

	grid->addItem(new QSpacerItem(16, 16,
		QSizePolicy::Minimum, QSizePolicy::Expanding),
		y++,0,1,4);

	grid->addWidget(new QLabel("Export Format:", m_exportpanel), y,0,1,1);
	grid->addWidget(m_exportFormat, y,1,1,1);
	grid->addWidget(btn_export, y++,3,1,1);

	// signals
	connect(btn_export, &QAbstractButton::clicked, this,
		static_cast<void (MagDynDlg::*)()>(&MagDynDlg::ExportSQE));


	m_tabs_out->addTab(m_exportpanel, "Export");
}



/**
 * notes dialog
 */
void MagDynDlg::ShowNotesDlg(bool only_create)
{
	if(!m_notes_dlg)
	{
		m_notes_dlg = new NotesDlg(this, m_sett);
		m_notes_dlg->setFont(this->font());
	}

	if(!only_create)
	{
		m_notes_dlg->show();
		m_notes_dlg->raise();
		m_notes_dlg->activateWindow();
	}
}



/**
 * about dialog
 */
void MagDynDlg::ShowInfoDlg(bool only_create)
{
	if(!m_info_dlg)
	{
		m_info_dlg = new InfoDlg(this, m_sett);
		m_info_dlg->setFont(this->font());
	}

	if(!only_create)
	{
		m_info_dlg->show();
		m_info_dlg->raise();
		m_info_dlg->activateWindow();
	}
}



/**
 * structure plotter dialog
 */
void MagDynDlg::ShowStructPlotDlg(bool only_create)
{
	if(!m_structplot_dlg)
	{
		m_structplot_dlg = new StructPlotDlg(this, m_sett, m_info_dlg);
		m_structplot_dlg->setFont(this->font());

		m_structplot_dlg->SetKernel(&m_dyn);
		m_structplot_dlg->SetTables(m_sitestab, m_termstab);

		QObject::connect(m_structplot_dlg, &StructPlotDlg::SelectSite,
			this, &MagDynDlg::SelectSite);
		QObject::connect(m_structplot_dlg, &StructPlotDlg::DeleteSite,
			this, &MagDynDlg::DeleteSite);
		QObject::connect(m_structplot_dlg, &StructPlotDlg::FlipSiteSpin,
			this, &MagDynDlg::FlipSiteSpin);

		QObject::connect(m_structplot_dlg, &StructPlotDlg::SelectTerm,
			this, &MagDynDlg::SelectTerm);
		QObject::connect(m_structplot_dlg, &StructPlotDlg::DeleteTerm,
			this, &MagDynDlg::DeleteTerm);
	}

	if(!only_create)
	{
		m_structplot_dlg->show();
		m_structplot_dlg->raise();
		m_structplot_dlg->activateWindow();
	}
}



/**
 * ground state minimiser dialog
 */
void MagDynDlg::ShowGroundStateDlg(bool only_create)
{
	if(!m_groundstate_dlg)
	{
		m_groundstate_dlg = new GroundStateDlg(this, m_sett);
		m_groundstate_dlg->setFont(this->font());

		m_groundstate_dlg->SetKernel(&m_dyn);
		QObject::connect(m_groundstate_dlg, &GroundStateDlg::SpinsUpdated,
			[this](const t_magdyn* dyn)
		{
			this->SetKernel(dyn, true, false, false);
		});
	}

	if(!only_create)
	{
		m_groundstate_dlg->show();
		m_groundstate_dlg->raise();
		m_groundstate_dlg->activateWindow();
	}
}



/**
 * main menu
 */
void MagDynDlg::CreateMenuBar()
{
	m_menu = new QMenuBar(this);

	// file menu
	auto menuFile = new QMenu("File", m_menu);
	auto acNew = new QAction("New", menuFile);
	auto acLoad = new QAction("Open...", menuFile);
	auto acImportStructure = new QAction("Import Structure...", menuFile);
	auto acSave = new QAction("Save", menuFile);
	auto acSaveAs = new QAction("Save As...", menuFile);
	auto acExit = new QAction("Quit", menuFile);

	// structure menu
	auto menuStruct = new QMenu("Structure", m_menu);
	auto acStructSymIdx = new QAction("Assign Symmetry Indices", menuStruct);
	auto acStructImport = new QAction("Import From Table...", menuStruct);
	auto acStructExportSun = new QAction("Export To Sunny Code...");
	auto acStructExportSW = new QAction("Export To SpinW Code...");
	auto acStructNotes = new QAction("Notes...", menuStruct);
	auto acStructView = new QAction("View...", menuStruct);
	auto acGroundState = new QAction("Minimise Ground State...", menuStruct);

	// dispersion menu
	m_menuDisp = new QMenu("Dispersion", m_menu);
	m_plot_channels = new QAction("Plot Channels", m_menuDisp);
	m_plot_channels->setToolTip("Plot individual polarisation channels.");
	m_plot_channels->setCheckable(true);
	m_plot_channels->setChecked(false);
	auto acRescalePlot = new QAction("Rescale Axes", m_menuDisp);
	auto acSaveFigure = new QAction("Save Figure...", m_menuDisp);
	auto acSaveDisp = new QAction("Save Data...", m_menuDisp);
	auto acSaveMultiDisp = new QAction("Save Data For All Qs...", m_menuDisp);
	auto acSaveDispScr = new QAction("Save Data As Script...", m_menuDisp);
	auto acSaveMultiDispScr = new QAction("Save Data As Script For All Qs...", m_menuDisp);

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
	auto menuCalc = new QMenu("Calculation", m_menu);
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
	auto menuTools = new QMenu("Tools", m_menu);
	auto acTrafoCalc = new QAction("Transformation Calculator...", menuTools);
	auto acPreferences = new QAction("Preferences...", menuTools);
	acTrafoCalc->setIcon(QIcon::fromTheme("accessories-calculator"));
	acPreferences->setIcon(QIcon::fromTheme("preferences-system"));

	acPreferences->setShortcut(QKeySequence::Preferences);
	acPreferences->setMenuRole(QAction::PreferencesRole);

	// help menu
	auto menuHelp = new QMenu("Help", m_menu);
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
			QMessageBox::critical(this, "Magnetic Dynamics", "Could not open the wiki.");
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
