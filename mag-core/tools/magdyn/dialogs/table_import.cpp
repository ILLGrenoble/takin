/**
 * import magnetic structure from a table
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
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

#include "table_import.h"
#include "defs.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QFrame>
#include <QtWidgets/QMenu>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QStyle>
#include <QtGui/QDesktopServices>

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
        #include <QtWidgets/QAction>
#else
        #include <QtGui/QAction>
#endif

#include "tlibs2/libs/str.h"



TableImportDlg::TableImportDlg(QWidget* parent, QSettings* sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Table Importer");
	setSizeGripEnabled(true);
	setFont(parent->font());

	// gui elements
	// --------------------------------------------------------------------------------
	QLabel *labelAtomIdx = new QLabel("Column Indices in Magnetic Sites Table:", this);
	QPushButton *btnSiteIndices = new QPushButton("Set Indices", this);
	QMenu *menuSiteIndices = new QMenu(this);
	btnSiteIndices->setMenu(menuSiteIndices);

	QAction *acSunSites = new QAction("Set Site Indices For Sunny", menuSiteIndices);
	acSunSites->setIcon(QIcon::fromTheme("weather-clear"));
	menuSiteIndices->addAction(acSunSites);
	menuSiteIndices->addSeparator();

	QAction *acSWAtomSites = new QAction("Set Site Indices For SpinW's \"Atom\" Table",
		menuSiteIndices);
	QAction *acSWMagSites = new QAction("Set Site Indices For SpinW's \"Mag\" Table",
		menuSiteIndices);
	menuSiteIndices->addAction(acSWAtomSites);
	menuSiteIndices->addAction(acSWMagSites);

	m_spinAtomName = new QSpinBox(this);
	m_spinAtomX = new QSpinBox(this);
	m_spinAtomY = new QSpinBox(this);
	m_spinAtomZ = new QSpinBox(this);
	m_spinAtomSX = new QSpinBox(this);
	m_spinAtomSY = new QSpinBox(this);
	m_spinAtomSZ = new QSpinBox(this);
	m_spinAtomSMag = new QSpinBox(this);

	m_spinAtomName->setPrefix("name = ");
	m_spinAtomX->setPrefix("x = ");
	m_spinAtomY->setPrefix("y = ");
	m_spinAtomZ->setPrefix("z = ");
	m_spinAtomSX->setPrefix("Sx = ");
	m_spinAtomSY->setPrefix("Sy = ");
	m_spinAtomSZ->setPrefix("Sz = ");
	m_spinAtomSMag->setPrefix("|S| = ");

	m_spinAtomName->setToolTip("Index the site's name.");
	m_spinAtomX->setToolTip("Index of the position vector's x component.");
	m_spinAtomY->setToolTip("Index of the position vector's y component.");
	m_spinAtomZ->setToolTip("Index of the position vector's z component.");
	m_spinAtomSX->setToolTip("Index of the spin vector's x component.");
	m_spinAtomSY->setToolTip("Index of the spin vector's y component.");
	m_spinAtomSZ->setToolTip("Index of the spin vector's z component.");
	m_spinAtomSMag->setToolTip("Index of the spin vector's magnitude.");

	for(QSpinBox* spin : {m_spinAtomName, m_spinAtomX, m_spinAtomY, m_spinAtomZ,
		m_spinAtomSX, m_spinAtomSY, m_spinAtomSZ, m_spinAtomSMag})
	{
		// -1 means not used
		spin->setMinimum(-1);
	}

	m_spinAtomName->setValue(0);
	m_spinAtomX->setValue(1);
	m_spinAtomY->setValue(2);
	m_spinAtomZ->setValue(3);
	m_spinAtomSX->setValue(4);
	m_spinAtomSY->setValue(5);
	m_spinAtomSZ->setValue(6);
	m_spinAtomSMag->setValue(7);

	QLabel *labelAtoms = new QLabel("Paste Magnetic Sites Table Here:", this);
	m_editAtoms = new QTextEdit(this);
	m_editAtoms->setLineWrapMode(QTextEdit::NoWrap);

	QFrame *sep1 = new QFrame(this);
	sep1->setFrameStyle(QFrame::HLine);
	// --------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------
	QLabel *labelCouplingIdx = new QLabel("Column Indices in Magnetic Couplings Table:", this);
	QPushButton *btnCouplingIndices = new QPushButton("Set Indices", this);
	QMenu *menuCouplingIndices = new QMenu(this);
	btnCouplingIndices->setMenu(menuCouplingIndices);

	QAction *acSunCouplings = new QAction("Set Coupling Indices For Sunny", menuCouplingIndices);
	acSunCouplings->setIcon(QIcon::fromTheme("weather-clear"));
	menuCouplingIndices->addAction(acSunCouplings);
	menuCouplingIndices->addSeparator();

	QAction *acSWCouplings = new QAction("Set Coupling Indices For SpinW's \"Bond\" Table",
		menuCouplingIndices);
	QAction *acSWSIA = new QAction("Set Coupling Indices For SpinW's \"Ion\" Table",
		menuCouplingIndices);
	menuCouplingIndices->addAction(acSWCouplings);
	menuCouplingIndices->addAction(acSWSIA);

	m_spinCouplingName = new QSpinBox(this);
	m_spinCouplingAtom1 = new QSpinBox(this);
	m_spinCouplingAtom2 = new QSpinBox(this);
	m_spinCouplingDX = new QSpinBox(this);
	m_spinCouplingDY = new QSpinBox(this);
	m_spinCouplingDZ = new QSpinBox(this);
	m_spinCouplingJ = new QSpinBox(this);
	m_spinCouplingDMIX = new QSpinBox(this);
	m_spinCouplingDMIY = new QSpinBox(this);
	m_spinCouplingDMIZ = new QSpinBox(this);
	m_spinCouplingJGen = new QSpinBox(this);

	m_spinCouplingName->setPrefix("name = ");
	m_spinCouplingAtom1->setPrefix("site1 = ");
	m_spinCouplingAtom2->setPrefix("site2 = ");
	m_spinCouplingDX->setPrefix("Δx = ");
	m_spinCouplingDY->setPrefix("Δy = ");
	m_spinCouplingDZ->setPrefix("Δz = ");
	m_spinCouplingJ->setPrefix("J = ");
	m_spinCouplingDMIX->setPrefix("DMIx = ");
	m_spinCouplingDMIY->setPrefix("DMIy = ");
	m_spinCouplingDMIZ->setPrefix("DMIz = ");
	m_spinCouplingJGen->setPrefix("Jgen = ");

	m_spinCouplingName->setToolTip("Index the coupling's name.");
	m_spinCouplingAtom1->setToolTip("Index of the first magnetic site in the coupling.");
	m_spinCouplingAtom2->setToolTip("Index of the second magnetic site in the coupling.");
	m_spinCouplingDX->setToolTip("Index unit cell vector's x component");
	m_spinCouplingDY->setToolTip("Index unit cell vector's y component");
	m_spinCouplingDZ->setToolTip("Index unit cell vector's z component");
	m_spinCouplingJ->setToolTip("Index of the exchange constant.");
	m_spinCouplingDMIX->setToolTip("Index of the DMI vector's x component.");
	m_spinCouplingDMIY->setToolTip("Index of the DMI vector's y component.");
	m_spinCouplingDMIZ->setToolTip("Index of the DMI vector's z component.");
	m_spinCouplingJGen->setToolTip("First index of the general coupling 3x3 matrix."
		"\nThe other 8 components are assumed to be in the subsequent columns.");

	for(QSpinBox* spin : {m_spinCouplingName, m_spinCouplingAtom1, m_spinCouplingAtom2,
		m_spinCouplingDX, m_spinCouplingDY, m_spinCouplingDZ, m_spinCouplingJ,
		m_spinCouplingDMIX, m_spinCouplingDMIY, m_spinCouplingDMIZ, m_spinCouplingJGen})
	{
		// -1 means not used
		spin->setMinimum(-1);
	}

	m_spinCouplingName->setValue(0);
	m_spinCouplingAtom1->setValue(1);
	m_spinCouplingAtom2->setValue(2);
	m_spinCouplingDX->setValue(3);
	m_spinCouplingDY->setValue(4);
	m_spinCouplingDZ->setValue(5);
	m_spinCouplingJ->setValue(6);
	m_spinCouplingDMIX->setValue(7);
	m_spinCouplingDMIY->setValue(8);
	m_spinCouplingDMIZ->setValue(9);
	m_spinCouplingJGen->setValue(-1);

	QLabel *labelCouplings = new QLabel("Paste Magnetic Couplings Table Here:", this);
	m_editCouplings = new QTextEdit(this);
	m_editCouplings->setLineWrapMode(QTextEdit::NoWrap);
	// --------------------------------------------------------------------------------

	QFrame *sep2 = new QFrame(this);
	sep2->setFrameStyle(QFrame::HLine);

	m_checkIndices1Based = new QCheckBox(this);
	m_checkUniteIncompleteTokens = new QCheckBox(this);
	m_checkIgnoreSymmetricCoupling = new QCheckBox(this);
	m_checkClearExisting = new QCheckBox(this);

	m_checkIndices1Based->setText("1-Based Indices");
	m_checkIndices1Based->setToolTip("Are the indices 1-based or 0-based?");
	m_checkIgnoreSymmetricCoupling->setText("Ignore Symmetric");
	m_checkIgnoreSymmetricCoupling->setToolTip("Ignore couplings with flipped site indices and inverted distance vectors.");
	m_checkUniteIncompleteTokens->setText("Unite Tokens");
	m_checkUniteIncompleteTokens->setToolTip("Unite incomplete tokens, e.g. bracket expressions.");
	m_checkClearExisting->setText("Clear Existing");
	m_checkClearExisting->setToolTip("Clears the existing sites or coupling before importing the new ones."
		"\nIf unchecked, new sites or couplings will be added to the end of the respective lists.");

	m_checkIndices1Based->setChecked(false);
	m_checkUniteIncompleteTokens->setChecked(true);
	m_checkIgnoreSymmetricCoupling->setChecked(false);
	m_checkClearExisting->setChecked(true);

	for(QCheckBox *box : { m_checkIndices1Based, m_checkUniteIncompleteTokens,
		m_checkIgnoreSymmetricCoupling, m_checkClearExisting })
		box->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

	QPushButton *btnImportAtoms = new QPushButton("Import Sites", this);
	QPushButton *btnImportCouplings = new QPushButton("Import Couplings", this);
	QPushButton *btnHelp = new QPushButton("Help", this);
	QPushButton *btnOK = new QPushButton("OK", this);

	btnHelp->setIcon(style()->standardIcon(QStyle::SP_DialogHelpButton));
	btnOK->setIcon(style()->standardIcon(QStyle::SP_DialogOkButton));

	for(QPushButton *btn : { btnImportAtoms, btnImportCouplings, btnHelp, btnOK })
		btn->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

	// grid
	QGridLayout* grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);
	int y = 0;
	grid->addWidget(labelAtomIdx, y, 0, 1, 3);
	grid->addWidget(btnSiteIndices, y++, 3, 1, 1);
	grid->addWidget(m_spinAtomName, y, 0, 1, 1);
	grid->addWidget(m_spinAtomX, y, 1, 1, 1);
	grid->addWidget(m_spinAtomY, y, 2, 1, 1);
	grid->addWidget(m_spinAtomZ, y++, 3, 1, 1);
	grid->addWidget(m_spinAtomSX, y, 0, 1, 1);
	grid->addWidget(m_spinAtomSY, y, 1, 1, 1);
	grid->addWidget(m_spinAtomSZ, y, 2, 1, 1);
	grid->addWidget(m_spinAtomSMag, y++, 3, 1, 1);
	grid->addWidget(labelAtoms, y++, 0, 1, 4);
	grid->addWidget(m_editAtoms, y++, 0, 1, 4);
	grid->addWidget(sep1, y++, 0, 1, 4);

	grid->addWidget(labelCouplingIdx, y, 0, 1, 3);
	grid->addWidget(btnCouplingIndices, y++, 3, 1, 1);
	grid->addWidget(m_spinCouplingName, y, 0, 1, 1);
	grid->addWidget(m_spinCouplingAtom1, y, 1, 1, 1);
	grid->addWidget(m_spinCouplingAtom2, y++, 2, 1, 1);
	grid->addWidget(m_spinCouplingDX, y, 0, 1, 1);
	grid->addWidget(m_spinCouplingDY, y, 1, 1, 1);
	grid->addWidget(m_spinCouplingDZ, y, 2, 1, 1);
	grid->addWidget(m_spinCouplingJ, y++, 3, 1, 1);
	grid->addWidget(m_spinCouplingDMIX, y, 0, 1, 1);
	grid->addWidget(m_spinCouplingDMIY, y, 1, 1, 1);
	grid->addWidget(m_spinCouplingDMIZ, y, 2, 1, 1);
	grid->addWidget(m_spinCouplingJGen, y++, 3, 1, 1);
	grid->addWidget(labelCouplings, y++, 0, 1, 4);
	grid->addWidget(m_editCouplings, y++, 0, 1, 4);
	grid->addWidget(sep2, y++, 0, 1, 4);
	grid->addWidget(m_checkIndices1Based, y, 0, 1, 1);
	grid->addWidget(m_checkUniteIncompleteTokens, y, 1, 1, 1);
	grid->addWidget(m_checkIgnoreSymmetricCoupling, y, 2, 1, 1);
	grid->addWidget(m_checkClearExisting, y++, 3, 1, 1);
	grid->addWidget(btnImportAtoms, y, 0, 1, 1);
	grid->addWidget(btnImportCouplings, y, 1, 1, 1);
	grid->addWidget(btnHelp, y, 2, 1, 1);
	grid->addWidget(btnOK, y++, 3, 1, 1);

	if(m_sett)
	{
		// restore window size and position
		if(m_sett->contains("tableimport/geo"))
			restoreGeometry(m_sett->value("tableimport/geo").toByteArray());
		else
			resize(500, 500);

		if(m_sett->contains("tableimport/idx_atom_name"))
			m_spinAtomName->setValue(m_sett->value("tableimport/idx_atom_name").toInt());
		if(m_sett->contains("tableimport/idx_atom_x"))
			m_spinAtomX->setValue(m_sett->value("tableimport/idx_atom_x").toInt());
		if(m_sett->contains("tableimport/idx_atom_y"))
			m_spinAtomY->setValue(m_sett->value("tableimport/idx_atom_y").toInt());
		if(m_sett->contains("tableimport/idx_atom_z"))
			m_spinAtomZ->setValue(m_sett->value("tableimport/idx_atom_z").toInt());
		if(m_sett->contains("tableimport/idx_atom_Sx"))
			m_spinAtomSX->setValue(m_sett->value("tableimport/idx_atom_Sx").toInt());
		if(m_sett->contains("tableimport/idx_atom_Sy"))
			m_spinAtomSY->setValue(m_sett->value("tableimport/idx_atom_Sy").toInt());
		if(m_sett->contains("tableimport/idx_atom_Sz"))
			m_spinAtomSZ->setValue(m_sett->value("tableimport/idx_atom_Sz").toInt());
		if(m_sett->contains("tableimport/idx_atom_Smag"))
			m_spinAtomSMag->setValue(m_sett->value("tableimport/idx_atom_Smag").toInt());
		if(m_sett->contains("tableimport/idx_coupling_name"))
			m_spinCouplingName->setValue(m_sett->value("tableimport/idx_coupling_name").toInt());
		if(m_sett->contains("tableimport/idx_coupling_atomidx_1"))
			m_spinCouplingAtom1->setValue(m_sett->value("tableimport/idx_coupling_atomidx_1").toInt());
		if(m_sett->contains("tableimport/idx_coupling_atomidx_2"))
			m_spinCouplingAtom2->setValue(m_sett->value("tableimport/idx_coupling_atomidx_2").toInt());
		if(m_sett->contains("tableimport/idx_coupling_Dx"))
			m_spinCouplingDX->setValue(m_sett->value("tableimport/idx_coupling_Dx").toInt());
		if(m_sett->contains("tableimport/idx_coupling_Dy"))
			m_spinCouplingDY->setValue(m_sett->value("tableimport/idx_coupling_Dy").toInt());
		if(m_sett->contains("tableimport/idx_coupling_Dz"))
			m_spinCouplingDZ->setValue(m_sett->value("tableimport/idx_coupling_Dz").toInt());
		if(m_sett->contains("tableimport/idx_coupling_J"))
			m_spinCouplingJ->setValue(m_sett->value("tableimport/idx_coupling_J").toInt());
		if(m_sett->contains("tableimport/idx_coupling_DMIx"))
			m_spinCouplingDMIX->setValue(m_sett->value("tableimport/idx_coupling_DMIx").toInt());
		if(m_sett->contains("tableimport/idx_coupling_DMIy"))
			m_spinCouplingDMIY->setValue(m_sett->value("tableimport/idx_coupling_DMIy").toInt());
		if(m_sett->contains("tableimport/idx_coupling_DMIz"))
			m_spinCouplingDMIZ->setValue(m_sett->value("tableimport/idx_coupling_DMIz").toInt());
		if(m_sett->contains("tableimport/idx_coupling_J_general"))
			m_spinCouplingJGen->setValue(m_sett->value("tableimport/idx_coupling_J_general").toInt());
		if(m_sett->contains("tableimport/coupling_indices_1based"))
			m_checkIndices1Based->setChecked(m_sett->value("tableimport/indices_1based").toBool());
		if(m_sett->contains("tableimport/unite_incomplete_tokens"))
			m_checkUniteIncompleteTokens->setChecked(m_sett->value("tableimport/unite_incomplete_tokens").toBool());
		if(m_sett->contains("tableimport/unite_incomplete_tokens"))
			m_checkIgnoreSymmetricCoupling->setChecked(m_sett->value("tableimport/ignore_symmetric_couplings").toBool());
		if(m_sett->contains("tableimport/clear_existing_items"))
			m_checkClearExisting->setChecked(m_sett->value("tableimport/clear_existing_items").toBool());
	}

	// connections
	connect(acSunSites, &QAction::triggered, [this]()
	{
		m_spinAtomName->setValue(0);
		m_spinAtomX->setValue(1);
		m_spinAtomY->setValue(2);
		m_spinAtomZ->setValue(3);
		m_spinAtomSX->setValue(4);
		m_spinAtomSY->setValue(5);
		m_spinAtomSZ->setValue(6);
		m_spinAtomSMag->setValue(7);

		m_checkIndices1Based->setChecked(true);
		m_checkUniteIncompleteTokens->setChecked(false);
		m_checkIgnoreSymmetricCoupling->setChecked(true);
		m_checkClearExisting->setChecked(true);
	});
	connect(acSWAtomSites, &QAction::triggered, [this]()
	{
		m_spinAtomName->setValue(0);
		m_spinAtomX->setValue(3);
		m_spinAtomY->setValue(4);
		m_spinAtomZ->setValue(5);
		m_spinAtomSX->setValue(-1);
		m_spinAtomSY->setValue(-1);
		m_spinAtomSZ->setValue(-1);
		m_spinAtomSMag->setValue(2);

		m_checkIndices1Based->setChecked(true);
		m_checkUniteIncompleteTokens->setChecked(true);
		m_checkIgnoreSymmetricCoupling->setChecked(false);
		m_checkClearExisting->setChecked(true);
	});
	connect(acSWMagSites, &QAction::triggered, [this]()
	{
		m_spinAtomName->setValue(1);
		m_spinAtomX->setValue(7);
		m_spinAtomY->setValue(8);
		m_spinAtomZ->setValue(9);
		m_spinAtomSX->setValue(4);
		m_spinAtomSY->setValue(5);
		m_spinAtomSZ->setValue(6);
		m_spinAtomSMag->setValue(3);

		m_checkIndices1Based->setChecked(true);
		m_checkUniteIncompleteTokens->setChecked(true);
		m_checkIgnoreSymmetricCoupling->setChecked(false);
		m_checkClearExisting->setChecked(true);
	});

	connect(acSunCouplings, &QAction::triggered, [this]()
	{
		m_spinCouplingName->setValue(0);
		m_spinCouplingAtom1->setValue(1);
		m_spinCouplingAtom2->setValue(2);
		m_spinCouplingDX->setValue(3);
		m_spinCouplingDY->setValue(4);
		m_spinCouplingDZ->setValue(5);
		m_spinCouplingJ->setValue(6);
		m_spinCouplingDMIX->setValue(7);
		m_spinCouplingDMIY->setValue(8);
		m_spinCouplingDMIZ->setValue(9);
		m_spinCouplingJGen->setValue(-1);

		m_checkIndices1Based->setChecked(true);
		m_checkUniteIncompleteTokens->setChecked(false);
		m_checkIgnoreSymmetricCoupling->setChecked(true);
		m_checkClearExisting->setChecked(true);
	});
	connect(acSWCouplings, &QAction::triggered, [this]()
	{
		m_spinCouplingName->setValue(-1);
		m_spinCouplingAtom1->setValue(10);
		m_spinCouplingAtom2->setValue(12);
		m_spinCouplingDX->setValue(2);
		m_spinCouplingDY->setValue(3);
		m_spinCouplingDZ->setValue(4);
		m_spinCouplingJ->setValue(16);
		m_spinCouplingDMIX->setValue(25);
		m_spinCouplingDMIY->setValue(26);
		m_spinCouplingDMIZ->setValue(27);
		m_spinCouplingJGen->setValue(-1);

		m_checkIndices1Based->setChecked(true);
		m_checkUniteIncompleteTokens->setChecked(true);
		m_checkIgnoreSymmetricCoupling->setChecked(false);
		m_checkClearExisting->setChecked(true);
	});
	connect(acSWSIA, &QAction::triggered, [this]()
	{
		m_spinCouplingName->setValue(-1);
		m_spinCouplingAtom1->setValue(1);
		m_spinCouplingAtom2->setValue(1);
		m_spinCouplingDX->setValue(-1);
		m_spinCouplingDY->setValue(-1);
		m_spinCouplingDZ->setValue(-1);
		m_spinCouplingJ->setValue(-1);
		m_spinCouplingDMIX->setValue(-1);
		m_spinCouplingDMIY->setValue(-1);
		m_spinCouplingDMIZ->setValue(-1);
		m_spinCouplingJGen->setValue(4);

		m_checkIndices1Based->setChecked(true);
		m_checkUniteIncompleteTokens->setChecked(true);
		m_checkIgnoreSymmetricCoupling->setChecked(false);
		m_checkClearExisting->setChecked(false);
	});

	connect(btnImportAtoms, &QAbstractButton::clicked, this, &TableImportDlg::ImportAtoms);
	connect(btnImportCouplings, &QAbstractButton::clicked, this, &TableImportDlg::ImportCouplings);
	connect(btnHelp, &QAbstractButton::clicked, this, &TableImportDlg::ShowHelp);
	connect(btnOK, &QAbstractButton::clicked, this, &TableImportDlg::accept);
}



TableImportDlg::~TableImportDlg()
{
}



/**
 * read in the magnetic sites from the table
 */
void TableImportDlg::ImportAtoms()
{
	// get the lines of the table
	std::string txt = m_editAtoms->toPlainText().toStdString();
	std::vector<std::string> lines;
	tl2::get_tokens<std::string, std::string>(txt, "\n", lines);

	const int idx_name = m_spinAtomName->value();
	const int idx_pos_x = m_spinAtomX->value();
	const int idx_pos_y = m_spinAtomY->value();
	const int idx_pos_z = m_spinAtomZ->value();
	const int idx_S_x = m_spinAtomSX->value();
	const int idx_S_y = m_spinAtomSY->value();
	const int idx_S_z = m_spinAtomSZ->value();
	const int idx_S_mag = m_spinAtomSMag->value();
	const bool unite_incomplete = m_checkUniteIncompleteTokens->isChecked();
	const bool clear_existing_sites = m_checkClearExisting->isChecked();

	std::vector<TableImportAtom> atompos_vec;
	atompos_vec.reserve(lines.size());

	for(const std::string& line : lines)
	{
		// get the column entries of the table
		std::vector<std::string> cols;
		tl2::get_tokens<std::string, std::string>(line, " \t", cols);
		if(unite_incomplete)
			cols = tl2::unite_incomplete_tokens<std::string>(cols, "([{", ")]}");

		TableImportAtom atompos;

		if(idx_name >= 0 && idx_name < int(cols.size()))
			atompos.name = cols[idx_name];

		if(idx_pos_x >= 0 && idx_pos_x < int(cols.size()))
			atompos.pos[0] = cols[idx_pos_x];
		if(idx_pos_y >= 0 && idx_pos_y < int(cols.size()))
			atompos.pos[1] = cols[idx_pos_y];
		if(idx_pos_z >= 0 && idx_pos_z < int(cols.size()))
			atompos.pos[2] = cols[idx_pos_z];

		if(idx_S_x >= 0 && idx_S_x < int(cols.size()))
			atompos.S[0] = cols[idx_S_x];
		if(idx_S_y >= 0 && idx_S_y < int(cols.size()))
			atompos.S[1] = cols[idx_S_y];
		if(idx_S_z >= 0 && idx_S_z < int(cols.size()))
			atompos.S[2] = cols[idx_S_z];

		if(idx_S_mag >= 0 && idx_S_mag < int(cols.size()))
			atompos.Smag = cols[idx_S_mag];

		atompos_vec.emplace_back(std::move(atompos));
	}

	emit SetAtomsSignal(atompos_vec, clear_existing_sites);
}



/**
 * read in the couplings from the table
 */
void TableImportDlg::ImportCouplings()
{
	// get the lines of the table
	std::string txt = m_editCouplings->toPlainText().toStdString();
	std::vector<std::string> lines;
	tl2::get_tokens<std::string, std::string>(txt, "\n", lines);

	const int idx_name = m_spinCouplingName->value();
	const int idx_atom1 = m_spinCouplingAtom1->value();
	const int idx_atom2 = m_spinCouplingAtom2->value();
	const int idx_dx = m_spinCouplingDX->value();
	const int idx_dy = m_spinCouplingDY->value();
	const int idx_dz = m_spinCouplingDZ->value();
	const int idx_J = m_spinCouplingJ->value();
	const int idx_dmix = m_spinCouplingDMIX->value();
	const int idx_dmiy = m_spinCouplingDMIY->value();
	const int idx_dmiz = m_spinCouplingDMIZ->value();
	const int idx_J_gen = m_spinCouplingJGen->value();
	const bool one_based = m_checkIndices1Based->isChecked();
	const bool unite_incomplete = m_checkUniteIncompleteTokens->isChecked();
	const bool ignore_symm = m_checkIgnoreSymmetricCoupling->isChecked();
	const bool clear_existing_couplings = m_checkClearExisting->isChecked();

	std::vector<TableImportCoupling> couplings;
	couplings.reserve(lines.size());

	for(const std::string& line : lines)
	{
		// get the column entries of the table
		std::vector<std::string> cols;
		tl2::get_tokens<std::string, std::string>(line, " \t", cols);
		if(unite_incomplete)
			cols = tl2::unite_incomplete_tokens<std::string>(cols);

		TableImportCoupling coupling;

		if(idx_name >= 0 && idx_name < int(cols.size()))
			coupling.name = cols[idx_name];
		if(idx_atom1 >= 0 && idx_atom1 < int(cols.size()))
		{
			coupling.atomidx1 = tl2::str_to_var<t_size>(cols[idx_atom1]);
			if(one_based)
				--*coupling.atomidx1;
		}
		if(idx_atom2 >= 0 && idx_atom2 < int(cols.size()))
		{
			coupling.atomidx2 = tl2::str_to_var<t_size>(cols[idx_atom2]);
			if(one_based)
				--*coupling.atomidx2;
		}
		if(idx_dx >= 0 && idx_dx < int(cols.size()))
			coupling.d[0] = cols[idx_dx];
		if(idx_dy >= 0 && idx_dy < int(cols.size()))
			coupling.d[1] = cols[idx_dy];
		if(idx_dz >= 0 && idx_dz < int(cols.size()))
			coupling.d[2] = cols[idx_dz];
		if(idx_J >= 0 && idx_J < int(cols.size()))
			coupling.J = cols[idx_J];
		if(idx_dmix >= 0 && idx_dmix < int(cols.size()))
			coupling.dmi[0] = cols[idx_dmix];
		if(idx_dmiy >= 0 && idx_dmiy < int(cols.size()))
			coupling.dmi[1] = cols[idx_dmiy];
		if(idx_dmiz >= 0 && idx_dmiz < int(cols.size()))
			coupling.dmi[2] = cols[idx_dmiz];
		if(idx_J_gen >= 0 && idx_J_gen + 8 < int(cols.size()))
		{
			for(std::size_t idx_j = 0; idx_j < 9; ++idx_j)
				coupling.Jgen[idx_j] = cols[idx_J_gen + idx_j];
		}

		if(ignore_symm && HasSymmetricCoupling(couplings, coupling))
			continue;

		couplings.emplace_back(std::move(coupling));
	}

	emit SetCouplingsSignal(couplings, clear_existing_couplings);
}



/**
 * does an equivalent coupling already exist in the couplings vector?
 */
bool TableImportDlg::HasSymmetricCoupling(
	const std::vector<TableImportCoupling>& couplings, const TableImportCoupling& coupling)
{
	// are the needed properties available?
	if(!coupling.atomidx1 || !coupling.atomidx2 ||
		!coupling.d[0].size() || !coupling.d[1].size() || !coupling.d[2].size())
		return false;

	for(const TableImportCoupling& c : couplings)
	{
		// are the needed properties available?
		if(!c.atomidx1 || !c.atomidx2 ||
			!c.d[0].size() || !c.d[1].size() || !c.d[2].size())
			continue;

		// TODO: this doesn't work if the strings in c.d[i] contain expressions
		if(*c.atomidx1 == *coupling.atomidx2 && *c.atomidx2 == *coupling.atomidx1
			&& tl2::equals<t_real>(tl2::str_to_var<t_real>(c.d[0]),
				-tl2::str_to_var<t_real>(coupling.d[0]), g_eps)
			&& tl2::equals<t_real>(tl2::str_to_var<t_real>(c.d[1]),
				-tl2::str_to_var<t_real>(coupling.d[1]), g_eps)
			&& tl2::equals<t_real>(tl2::str_to_var<t_real>(c.d[2]),
				-tl2::str_to_var<t_real>(coupling.d[2]), g_eps))
		{
			return true;
		}
	}

	return false;
}



/**
 * show the wiki's help page for importing structures
 */
void TableImportDlg::ShowHelp()
{
	QUrl url("https://github.com/ILLGrenoble/takin/wiki/Importing-Magnetic-Structures");
	if(!QDesktopServices::openUrl(url))
		QMessageBox::critical(this, windowTitle() + " -- Error", "Could not open the wiki.");
}



/**
 * dialog is closing
 */
void TableImportDlg::closeEvent(QCloseEvent *)
{
	if(!m_sett)
		return;

	m_sett->setValue("tableimport/geo", saveGeometry());

	m_sett->setValue("tableimport/idx_atom_name", m_spinAtomName->value());
	m_sett->setValue("tableimport/idx_atom_x", m_spinAtomX->value());
	m_sett->setValue("tableimport/idx_atom_y", m_spinAtomY->value());
	m_sett->setValue("tableimport/idx_atom_z", m_spinAtomZ->value());
	m_sett->setValue("tableimport/idx_atom_Sx", m_spinAtomSX->value());
	m_sett->setValue("tableimport/idx_atom_Sy", m_spinAtomSY->value());
	m_sett->setValue("tableimport/idx_atom_Sz", m_spinAtomSZ->value());
	m_sett->setValue("tableimport/idx_atom_Smag", m_spinAtomSMag->value());

	m_sett->setValue("tableimport/idx_coupling_name", m_spinCouplingName->value());
	m_sett->setValue("tableimport/idx_coupling_atomidx_1", m_spinCouplingAtom1->value());
	m_sett->setValue("tableimport/idx_coupling_atomidx_2", m_spinCouplingAtom2->value());
	m_sett->setValue("tableimport/idx_coupling_Dx", m_spinCouplingDX->value());
	m_sett->setValue("tableimport/idx_coupling_Dy", m_spinCouplingDY->value());
	m_sett->setValue("tableimport/idx_coupling_Dz", m_spinCouplingDZ->value());
	m_sett->setValue("tableimport/idx_coupling_J", m_spinCouplingJ->value());
	m_sett->setValue("tableimport/idx_coupling_DMIx", m_spinCouplingDMIX->value());
	m_sett->setValue("tableimport/idx_coupling_DMIy", m_spinCouplingDMIY->value());
	m_sett->setValue("tableimport/idx_coupling_DMIz", m_spinCouplingDMIZ->value());
	m_sett->setValue("tableimport/idx_coupling_J_general", m_spinCouplingJGen->value());

	m_sett->setValue("tableimport/indices_1based", m_checkIndices1Based->isChecked());
	m_sett->setValue("tableimport/unite_incomplete_tokens", m_checkUniteIncompleteTokens->isChecked());
	m_sett->setValue("tableimport/ignore_symmetric_couplings", m_checkIgnoreSymmetricCoupling->isChecked());
	m_sett->setValue("tableimport/clear_existing_items", m_checkClearExisting->isChecked());
}
