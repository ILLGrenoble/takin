/**
 * magnetic dynamics -- table management
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

#include <iostream>
#include <unordered_set>
#include <boost/scope_exit.hpp>


// types
using t_numitem = tl2::NumericTableWidgetItem<t_real>;
using t_sizeitem = tl2::NumericTableWidgetItem<t_size>;



/**
 * set a unique name for the given table item
 */
static void set_unique_tab_item_name(
	QTableWidget *tab, QTableWidgetItem *item,
	int name_col, const std::string& prefix)
{
	if(!item)
		return;

	// was the item renamed?
	if(tab->column(item) != name_col)
		return;

	int row = tab->row(item);
	if(row < 0 || row >= tab->rowCount())
		return;

	// check if the new name is unique and non-empty
	std::string new_name_base = item->text().toStdString();
	std::string new_name = new_name_base;
	if(tl2::trimmed(new_name) == "")
		new_name_base = new_name = prefix;

	// hash all other names
	std::unordered_set<std::string> used_names;
	for(int idx = 0; idx < tab->rowCount(); ++idx)
	{
		auto *name = tab->item(idx, name_col);
		if(!name || idx == row)
			continue;

		used_names.emplace(name->text().toStdString());
	}

	// possibly add a suffix to make the new name unique
	std::size_t ctr = 1;
	while(used_names.find(new_name) != used_names.end())
		new_name = new_name_base + "_" + tl2::var_to_str(ctr++);

	// rename the site if needed
	if(new_name != item->text().toStdString())
		item->setText(new_name.c_str());
}



/**
 * add an atom site
 */
void MagDynDlg::AddSiteTabItem(
	int row, const std::string& name, t_size sym_idx,
	const std::string& x, const std::string& y, const std::string& z,
	const std::string& sx, const std::string& sy, const std::string& sz,
	const std::string& S,
	const std::string& sox, const std::string& soy, const std::string& soz,
	const std::string& rgb)
{
	bool bclone = false;
	m_sitestab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_sitestab->blockSignals(false);
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	if(row == -1)	// append to end of table
		row = m_sitestab->rowCount();
	else if(row == -2 && m_sites_cursor_row >= 0)	// use row from member variable
		row = m_sites_cursor_row;
	else if(row == -3 && m_sites_cursor_row >= 0)	// use row from member variable +1
		row = m_sites_cursor_row + 1;
	else if(row == -4 && m_sites_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_sites_cursor_row + 1;
		bclone = true;
	}
	m_sitestab->setSortingEnabled(false);
	m_sitestab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < m_sitestab->columnCount(); ++thecol)
		{
			m_sitestab->setItem(row, thecol,
				m_sitestab->item(m_sites_cursor_row, thecol)->clone());
		}

		m_sitestab->item(row, COL_SITE_NAME)->setText(
			m_sitestab->item(row, COL_SITE_NAME)->text() + QString("_clone"));
	}
	else
	{
		m_sitestab->setItem(row, COL_SITE_NAME,
			new QTableWidgetItem(name.c_str()));
		m_sitestab->setItem(row, COL_SITE_POS_X, new t_numitem(x));
		m_sitestab->setItem(row, COL_SITE_POS_Y, new t_numitem(y));
		m_sitestab->setItem(row, COL_SITE_POS_Z, new t_numitem(z));
		m_sitestab->setItem(row, COL_SITE_SYM_IDX, new t_sizeitem(sym_idx));
		m_sitestab->setItem(row, COL_SITE_SPIN_X, new t_numitem(sx));
		m_sitestab->setItem(row, COL_SITE_SPIN_Y, new t_numitem(sy));
		m_sitestab->setItem(row, COL_SITE_SPIN_Z, new t_numitem(sz));
		m_sitestab->setItem(row, COL_SITE_SPIN_MAG, new t_numitem(S));
		m_sitestab->setItem(row, COL_SITE_RGB, new QTableWidgetItem(rgb.c_str()));
		if(m_allow_ortho_spin)
		{
			m_sitestab->setItem(row, COL_SITE_SPIN_ORTHO_X, new t_numitem(sox));
			m_sitestab->setItem(row, COL_SITE_SPIN_ORTHO_Y, new t_numitem(soy));
			m_sitestab->setItem(row, COL_SITE_SPIN_ORTHO_Z, new t_numitem(soz));
		}
	}

	set_unique_tab_item_name(m_sitestab, m_sitestab->item(row, COL_SITE_NAME),
		COL_SITE_NAME, "site");

	m_sitestab->scrollToItem(m_sitestab->item(row, 0));
	m_sitestab->setCurrentCell(row, 0);
	m_sitestab->setSortingEnabled(true);

	UpdateVerticalHeader(m_sitestab);
	SyncSiteComboBoxes();
}



/**
 * update the contents of all site selection combo boxes to match the sites table
 */
void MagDynDlg::SyncSiteComboBoxes()
{
	if(m_ignoreSitesCalc)
		return;

	// iterate couplings and update their site selection combo boxes
	for(int row = 0; row < m_termstab->rowCount(); ++row)
	{
		SitesComboBox *site_1 = reinterpret_cast<SitesComboBox*>(
			m_termstab->cellWidget(row, COL_XCH_ATOM1_IDX));
		SitesComboBox *site_2 = reinterpret_cast<SitesComboBox*>(
			m_termstab->cellWidget(row, COL_XCH_ATOM2_IDX));

		SyncSiteComboBox(site_1, site_1->currentText().toStdString());
		SyncSiteComboBox(site_2, site_2->currentText().toStdString());
	}
}



/**
 * update the contents of a site selection combo box to match the sites table
 */
void MagDynDlg::SyncSiteComboBox(SitesComboBox* combo, const std::string& selected_site)
{
	BOOST_SCOPE_EXIT(combo)
	{
		combo->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	combo->blockSignals(true);

	combo->clear();

	int selected_idx = -1;
	int selected_idx_alt = -1;  // alternate selection in case of renamed site

	// iterate sites to get their names
	std::unordered_set<std::string> seen_names;
	for(int row = 0; row < m_sitestab->rowCount(); ++row)
	{
		auto *name = m_sitestab->item(row, COL_SITE_NAME);
		if(!name)
			continue;

		std::string site_name = name->text().toStdString();
		if(seen_names.find(site_name) != seen_names.end())
			continue;
		seen_names.insert(site_name);

		// add the item
		combo->addItem(name->text());

		// check if this item has to be selected
		if(name->text().toStdString() == selected_site)
			selected_idx = row;

		// check if the selection corresponds to the site's previous name
		std::string old_name = name->data(DATA_ROLE_NAME).toString().toStdString();
		if(old_name != "" && old_name == selected_site)
			selected_idx_alt = row;
	}

	combo->addItem("<invalid>");

	if(selected_idx >= 0)
	{
		// select the site
		combo->setCurrentIndex(selected_idx);
	}
	else if(selected_idx_alt >= 0)
	{
		// use alternate selection in case of a renamed site
		combo->setCurrentIndex(selected_idx_alt);
	}
	else
	{
		// set selection to invalid
		combo->setCurrentIndex(combo->count() - 1);
	}
}



/**
 * create a combo box with the site names
 */
SitesComboBox* MagDynDlg::CreateSitesComboBox(const std::string& selected_site)
{
	/**
	 * filter out wheel events (used for the combo boxes in the table)
	 */
	struct WheelEventFilter : public QObject
	{
		WheelEventFilter(QObject* comp) : m_comp{comp} {}
		WheelEventFilter(const WheelEventFilter&) = delete;
		WheelEventFilter& operator=(const WheelEventFilter&) = delete;

		bool eventFilter(QObject *comp, QEvent *evt) override
		{
			// ignore wheel events on the given component
			return comp == m_comp && evt->type() == QEvent::Wheel;
		}

		QObject* m_comp = nullptr;
	};


	SitesComboBox* combo = new SitesComboBox();
	combo->setParent(m_termstab);
	combo->installEventFilter(new WheelEventFilter(combo));
	SyncSiteComboBox(combo, selected_site);

	// connections
	connect(combo, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
		[this]()
		{
			TermsTableItemChanged(nullptr);
		});

	return combo;
}



/**
 * add an exchange term
 */
void MagDynDlg::AddTermTabItem(
	int row, const std::string& name, t_size sym_idx,
	const std::string& atom_1, const std::string& atom_2,
	const std::string& dist_x, const std::string& dist_y, const std::string& dist_z,
	const std::string& J,
	const std::string& dmi_x, const std::string& dmi_y, const std::string& dmi_z,
	const std::string& gen_xx, const std::string& gen_xy, const std::string& gen_xz,
	const std::string& gen_yx, const std::string& gen_yy, const std::string& gen_yz,
	const std::string& gen_zx, const std::string& gen_zy, const std::string& gen_zz,
	const std::string& rgb)
{
	bool bclone = false;
	m_termstab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_termstab->blockSignals(false);
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	if(row == -1)	// append to end of table
		row = m_termstab->rowCount();
	else if(row == -2 && m_terms_cursor_row >= 0)	// use row from member variable
		row = m_terms_cursor_row;
	else if(row == -3 && m_terms_cursor_row >= 0)	// use row from member variable +1
		row = m_terms_cursor_row + 1;
	else if(row == -4 && m_terms_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_terms_cursor_row + 1;
		bclone = true;
	}

	m_termstab->setSortingEnabled(false);
	m_termstab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < m_termstab->columnCount(); ++thecol)
		{
			m_termstab->setItem(row, thecol,
				m_termstab->item(m_terms_cursor_row, thecol)->clone());

			// also clone site selection combo boxes
			if(thecol == COL_XCH_ATOM1_IDX || thecol == COL_XCH_ATOM2_IDX)
			{
				SitesComboBox *combo_old = reinterpret_cast<SitesComboBox*>(
					m_termstab->cellWidget(m_terms_cursor_row, thecol));
				SitesComboBox *combo = CreateSitesComboBox(
					combo_old->currentText().toStdString());
				m_termstab->setCellWidget(row, thecol, combo);
				m_termstab->setItem(row, thecol, combo);
			}
		}

		m_termstab->item(row, COL_XCH_NAME)->setText(
			m_termstab->item(row, COL_XCH_NAME)->text() + QString("_clone"));
	}
	else
	{
		m_termstab->setItem(row, COL_XCH_NAME,
			new QTableWidgetItem(name.c_str()));

		SitesComboBox* combo1 = CreateSitesComboBox(atom_1);
		SitesComboBox* combo2 = CreateSitesComboBox(atom_2);
		m_termstab->setCellWidget(row, COL_XCH_ATOM1_IDX, combo1);
		m_termstab->setCellWidget(row, COL_XCH_ATOM2_IDX, combo2);
		m_termstab->setItem(row, COL_XCH_ATOM1_IDX, combo1);
		m_termstab->setItem(row, COL_XCH_ATOM2_IDX, combo2);

		m_termstab->setItem(row, COL_XCH_DIST_X, new t_numitem(dist_x));
		m_termstab->setItem(row, COL_XCH_DIST_Y, new t_numitem(dist_y));
		m_termstab->setItem(row, COL_XCH_DIST_Z, new t_numitem(dist_z));
		m_termstab->setItem(row, COL_XCH_SYM_IDX, new t_sizeitem(sym_idx));
		m_termstab->setItem(row, COL_XCH_INTERACTION, new t_numitem(J));
		m_termstab->setItem(row, COL_XCH_DMI_X, new t_numitem(dmi_x));
		m_termstab->setItem(row, COL_XCH_DMI_Y, new t_numitem(dmi_y));
		m_termstab->setItem(row, COL_XCH_DMI_Z, new t_numitem(dmi_z));
		m_termstab->setItem(row, COL_XCH_RGB, new QTableWidgetItem(rgb.c_str()));
		if(m_allow_general_J)
		{
			m_termstab->setItem(row, COL_XCH_GEN_XX, new t_numitem(gen_xx));
			m_termstab->setItem(row, COL_XCH_GEN_XY, new t_numitem(gen_xy));
			m_termstab->setItem(row, COL_XCH_GEN_XZ, new t_numitem(gen_xz));
			m_termstab->setItem(row, COL_XCH_GEN_YX, new t_numitem(gen_yx));
			m_termstab->setItem(row, COL_XCH_GEN_YY, new t_numitem(gen_yy));
			m_termstab->setItem(row, COL_XCH_GEN_YZ, new t_numitem(gen_yz));
			m_termstab->setItem(row, COL_XCH_GEN_ZX, new t_numitem(gen_zx));
			m_termstab->setItem(row, COL_XCH_GEN_ZY, new t_numitem(gen_zy));
			m_termstab->setItem(row, COL_XCH_GEN_ZZ, new t_numitem(gen_zz));
		}
	}

	set_unique_tab_item_name(m_termstab, m_termstab->item(row, COL_XCH_NAME),
		COL_XCH_NAME, "coupling");

	m_termstab->scrollToItem(m_termstab->item(row, 0));
	m_termstab->setCurrentCell(row, 0);
	m_termstab->setSortingEnabled(true);

	UpdateVerticalHeader(m_termstab);
}



/**
 * add a variable
 */
void MagDynDlg::AddVariableTabItem(int row, const std::string& name, const t_cplx& value)
{
	bool bclone = false;
	m_varstab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_varstab->blockSignals(false);
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	if(row == -1)	// append to end of table
		row = m_varstab->rowCount();
	else if(row == -2 && m_variables_cursor_row >= 0)	// use row from member variable
		row = m_variables_cursor_row;
	else if(row == -3 && m_variables_cursor_row >= 0)	// use row from member variable +1
		row = m_variables_cursor_row + 1;
	else if(row == -4 && m_variables_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_variables_cursor_row + 1;
		bclone = true;
	}

	m_varstab->setSortingEnabled(false);
	m_varstab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < m_varstab->columnCount(); ++thecol)
		{
			m_varstab->setItem(row, thecol,
				m_varstab->item(m_variables_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_varstab->setItem(row, COL_VARS_NAME, new QTableWidgetItem(name.c_str()));
		m_varstab->setItem(row, COL_VARS_VALUE_REAL, new t_numitem(value.real()));
		m_varstab->setItem(row, COL_VARS_VALUE_IMAG, new t_numitem(value.imag()));
	}

	set_unique_tab_item_name(m_varstab, m_varstab->item(row, COL_VARS_NAME),
		COL_VARS_NAME, "var");

	m_varstab->scrollToItem(m_varstab->item(row, 0));
	m_varstab->setCurrentCell(row, 0);
	m_varstab->setSortingEnabled(true);

	UpdateVerticalHeader(m_varstab);
}



/**
 * add a magnetic field
 */
void MagDynDlg::AddFieldTabItem(int row,
	t_real Bh, t_real Bk, t_real Bl,
	t_real Bmag)
{
	bool bclone = false;

	if(row == -1)	// append to end of table
		row = m_fieldstab->rowCount();
	else if(row == -2 && m_fields_cursor_row >= 0)	// use row from member variable
		row = m_fields_cursor_row;
	else if(row == -3 && m_fields_cursor_row >= 0)	// use row from member variable +1
		row = m_fields_cursor_row + 1;
	else if(row == -4 && m_fields_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_fields_cursor_row + 1;
		bclone = true;
	}

	m_fieldstab->setSortingEnabled(false);
	m_fieldstab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < m_fieldstab->columnCount(); ++thecol)
		{
			m_fieldstab->setItem(row, thecol,
				m_fieldstab->item(m_fields_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_fieldstab->setItem(row, COL_FIELD_H, new t_numitem(Bh));
		m_fieldstab->setItem(row, COL_FIELD_K, new t_numitem(Bk));
		m_fieldstab->setItem(row, COL_FIELD_L, new t_numitem(Bl));
		m_fieldstab->setItem(row, COL_FIELD_MAG, new t_numitem(Bmag));
	}

	m_fieldstab->scrollToItem(m_fieldstab->item(row, 0));
	m_fieldstab->setCurrentCell(row, 0);
	m_fieldstab->setSortingEnabled(true);

	UpdateVerticalHeader(m_fieldstab);
}



/**
 * add a coordinate path
 */
void MagDynDlg::AddCoordinateTabItem(int row, const std::string& name,
	t_real h, t_real k, t_real l)
{
	bool bclone = false;

	if(row == -1)	// append to end of table
		row = m_coordinatestab->rowCount();
	else if(row == -2 && m_coordinates_cursor_row >= 0)	// use row from member variable
		row = m_coordinates_cursor_row;
	else if(row == -3 && m_coordinates_cursor_row >= 0)	// use row from member variable +1
		row = m_coordinates_cursor_row + 1;
	else if(row == -4 && m_coordinates_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_coordinates_cursor_row + 1;
		bclone = true;
	}

	m_coordinatestab->setSortingEnabled(false);
	m_coordinatestab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < m_coordinatestab->columnCount(); ++thecol)
		{
			m_coordinatestab->setItem(row, thecol,
				m_coordinatestab->item(m_coordinates_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_coordinatestab->setItem(row, COL_COORD_NAME,
			new QTableWidgetItem(name.c_str()));
		m_coordinatestab->setItem(row, COL_COORD_H, new t_numitem(h));
		m_coordinatestab->setItem(row, COL_COORD_K, new t_numitem(k));
		m_coordinatestab->setItem(row, COL_COORD_L, new t_numitem(l));
	}

	m_coordinatestab->scrollToItem(m_coordinatestab->item(row, 0));
	m_coordinatestab->setCurrentCell(row, 0);
	m_coordinatestab->setSortingEnabled(true);

	UpdateVerticalHeader(m_coordinatestab);
}



/**
 * delete table widget items
 */
void MagDynDlg::DelTabItem(QTableWidget *pTab, int begin, int end)
{
	bool needs_recalc = true;
	if(pTab == m_fieldstab || pTab == m_coordinatestab)
		needs_recalc = false;

	if(needs_recalc)
		pTab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_, pTab, needs_recalc)
	{
		if(needs_recalc)
		{
			pTab->blockSignals(false);
			if(this_->m_autocalc->isChecked())
				this_->CalcAll();
		}
	} BOOST_SCOPE_EXIT_END

	// if nothing is selected, clear all items
	if(begin == -1 || pTab->selectedItems().count() == 0)
	{
		pTab->clearContents();
		pTab->setRowCount(0);
	}
	else if(begin == -2)	// clear selected
	{
		for(int row : GetSelectedRows(pTab, true))
			pTab->removeRow(row);
	}
	else if(begin >= 0 && end >= 0)		// clear given range
	{
		for(int row=end-1; row>=begin; --row)
			pTab->removeRow(row);
	}

	UpdateVerticalHeader(pTab);
	if(pTab == m_sitestab)
		SyncSiteComboBoxes();
}



void MagDynDlg::MoveTabItemUp(QTableWidget *pTab)
{
	bool needs_recalc = true;
	if(pTab == m_fieldstab || pTab == m_coordinatestab)
		needs_recalc = false;

	if(needs_recalc)
		pTab->blockSignals(true);
	pTab->setSortingEnabled(false);
	BOOST_SCOPE_EXIT(this_, pTab, needs_recalc)
	{
		if(needs_recalc)
		{
			pTab->blockSignals(false);
			if(this_->m_autocalc->isChecked())
				this_->CalcAll();
		}
	} BOOST_SCOPE_EXIT_END

	auto selected = GetSelectedRows(pTab, false);
	for(int row : selected)
	{
		if(row == 0)
			continue;

		auto *item = pTab->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		pTab->insertRow(row-1);
		for(int col = 0; col < pTab->columnCount(); ++col)
		{
			pTab->setItem(row-1, col, pTab->item(row+1, col)->clone());

			// also clone site selection combo boxes
			if(pTab == m_termstab && (col == COL_XCH_ATOM1_IDX || col == COL_XCH_ATOM2_IDX))
			{
				SitesComboBox *combo_old = reinterpret_cast<SitesComboBox*>(
					pTab->cellWidget(row+1, col));
				SitesComboBox *combo = CreateSitesComboBox(
					combo_old->currentText().toStdString());
				pTab->setCellWidget(row-1, col, combo);
				pTab->setItem(row-1, col, combo);
			}
		}
		pTab->removeRow(row+1);
	}

	for(int row = 0; row < pTab->rowCount(); ++row)
	{
		if(auto *item = pTab->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row+1) != selected.end())
		{
			for(int col=0; col<pTab->columnCount(); ++col)
				pTab->item(row, col)->setSelected(true);
		}
	}

	UpdateVerticalHeader(pTab);
}



void MagDynDlg::MoveTabItemDown(QTableWidget *pTab)
{
	bool needs_recalc = true;
	if(pTab == m_fieldstab || pTab == m_coordinatestab)
		needs_recalc = false;

	if(needs_recalc)
		pTab->blockSignals(true);
	pTab->setSortingEnabled(false);
	BOOST_SCOPE_EXIT(this_, pTab, needs_recalc)
	{
		if(needs_recalc)
		{
			pTab->blockSignals(false);
			if(this_->m_autocalc->isChecked())
				this_->CalcAll();
		}
	} BOOST_SCOPE_EXIT_END

	auto selected = GetSelectedRows(pTab, true);
	for(int row : selected)
	{
		if(row == pTab->rowCount()-1)
			continue;

		auto *item = pTab->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		pTab->insertRow(row+2);
		for(int col = 0; col < pTab->columnCount(); ++col)
		{
			pTab->setItem(row+2, col, pTab->item(row, col)->clone());

			// also clone site selection combo boxes
			if(pTab == m_termstab && (col == COL_XCH_ATOM1_IDX || col == COL_XCH_ATOM2_IDX))
			{
				SitesComboBox *combo_old = reinterpret_cast<SitesComboBox*>(
					pTab->cellWidget(row, col));
				SitesComboBox *combo = CreateSitesComboBox(
					combo_old->currentText().toStdString());
				pTab->setCellWidget(row+2, col, combo);
				pTab->setItem(row+2, col, combo);
			}
		}
		pTab->removeRow(row);
	}

	for(int row = 0; row < pTab->rowCount(); ++row)
	{
		if(auto *item = pTab->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row-1) != selected.end())
		{
			for(int col=0; col<pTab->columnCount(); ++col)
				pTab->item(row, col)->setSelected(true);
		}
	}

	UpdateVerticalHeader(pTab);
}



/**
 * insert a vertical header column showing the row index
 */
void MagDynDlg::UpdateVerticalHeader(QTableWidget *pTab)
{
	for(int row = 0; row < pTab->rowCount(); ++row)
	{
		QTableWidgetItem *item = pTab->verticalHeaderItem(row);
		if(!item)
			item = new QTableWidgetItem{};
		item->setText(QString::number(row + 1));
		pTab->setVerticalHeaderItem(row, item);
	}
}



std::vector<int> MagDynDlg::GetSelectedRows(QTableWidget *pTab, bool sort_reversed) const
{
	std::vector<int> vec;
	vec.reserve(pTab->selectedItems().size());

	for(int row=0; row<pTab->rowCount(); ++row)
	{
		if(auto *item = pTab->item(row, 0); item && item->isSelected())
			vec.push_back(row);
	}

	if(sort_reversed)
	{
		std::stable_sort(vec.begin(), vec.end(), [](int row1, int row2)
		{ return row1 > row2; });
	}

	return vec;
}



/**
 * sites table item contents changed
 */
void MagDynDlg::SitesTableItemChanged(QTableWidgetItem *item)
{
	if(!item)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_sitestab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_sitestab->blockSignals(true);

	// was the site renamed?
	if(m_sitestab->column(item) == COL_SITE_NAME)
	{
		int row = m_sitestab->row(item);
		if(row >= 0 && row < m_sitestab->rowCount() && row < (int)m_dyn.GetMagneticSitesCount())
		{
			// get the previous name of the site
			std::string old_name = m_dyn.GetMagneticSite(row).name;
			// and set it as temporary, alternate site name
			item->setData(DATA_ROLE_NAME, QString(old_name.c_str()));

			set_unique_tab_item_name(m_sitestab, item, COL_SITE_NAME, "site");
		}
	}

	SyncSiteComboBoxes();

	// clear alternate name
	item->setData(DATA_ROLE_NAME, QVariant());

	if(m_autocalc->isChecked())
		CalcAll();
}



/**
 * coupling table item contents changed
 */
void MagDynDlg::TermsTableItemChanged(QTableWidgetItem *item)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_termstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_termstab->blockSignals(true);

	set_unique_tab_item_name(m_termstab, item, COL_XCH_NAME, "coupling");

	if(m_autocalc->isChecked())
		CalcAll();
}



/**
 * variable table item contents changed
 */
void MagDynDlg::VariablesTableItemChanged(QTableWidgetItem *item)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_varstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_varstab->blockSignals(true);

	set_unique_tab_item_name(m_varstab, item, COL_VARS_NAME, "var");

	if(m_autocalc->isChecked())
		CalcAll();
}



/**
 * show the popup menu for the tables
 */
void MagDynDlg::ShowTableContextMenu(
	QTableWidget *pTab, QMenu *pMenu, QMenu *pMenuNoItem, const QPoint& ptLocal)
{
	QPoint pt = ptLocal;
	// transform the point to widget coordinates first if it has a viewport
	if(pTab->viewport())
		pt = pTab->viewport()->mapToParent(pt);

	// transform the point to global coordinates
	auto ptGlob = pTab->mapToGlobal(pt);

	if(const auto* item = pTab->itemAt(ptLocal); item)
		pMenu->popup(ptGlob);
	else
		pMenuNoItem->popup(ptGlob);
}



/**
 * set the selection to the given site in the corresponding table
 */
void MagDynDlg::SelectSite(const std::string& site)
{
	if(t_size idx = m_dyn.GetMagneticSiteIndex(site);
	   idx < m_dyn.GetMagneticSitesCount())
	{
		// select current site in table
		m_tabs_in->setCurrentWidget(m_sitespanel);
		m_sitestab->setCurrentCell(idx, 0);
	}
}



/**
 * delete the given site from the corresponding table
 */
void MagDynDlg::DeleteSite(const std::string& site)
{
	if(t_size idx = m_dyn.GetMagneticSiteIndex(site);
		idx < m_dyn.GetMagneticSitesCount())
	{
		DelTabItem(m_sitestab, idx, idx + 1);
	}
}



/**
 * invert the site's spin
 */
void MagDynDlg::FlipSiteSpin(const std::string& site)
{
	if(t_size idx = m_dyn.GetMagneticSiteIndex(site);
		idx < m_dyn.GetMagneticSitesCount())
	{
		auto *spin_x = static_cast<t_numitem*>(m_sitestab->item(idx, COL_SITE_SPIN_X));
		auto *spin_y = static_cast<t_numitem*>(m_sitestab->item(idx, COL_SITE_SPIN_Y));
		auto *spin_z = static_cast<t_numitem*>(m_sitestab->item(idx, COL_SITE_SPIN_Z));

		spin_x->SetValue(-spin_x->GetValue());
		spin_y->SetValue(-spin_y->GetValue());
		spin_z->SetValue(-spin_z->GetValue());
	}
}



/**
 * set the selection to the given term in the corresponding table
 */
void MagDynDlg::SelectTerm(const std::string& term)
{
	if(t_size idx = m_dyn.GetExchangeTermIndex(term);
	   idx < m_dyn.GetExchangeTermsCount())
	{
		// select current term in table
		m_tabs_in->setCurrentWidget(m_termspanel);
		m_termstab->setCurrentCell(idx, 0);
	}
}



/**
 * delete the given term from the corresponding table
 */
void MagDynDlg::DeleteTerm(const std::string& term)
{
	if(t_size idx = m_dyn.GetExchangeTermIndex(term);
		idx < m_dyn.GetExchangeTermsCount())
	{
		DelTabItem(m_termstab, idx, idx + 1);
	}
}



/**
 * get all entries from the variables table
 */
std::vector<t_magdyn::Variable> MagDynDlg::GetVariables() const
{
	std::vector<t_magdyn::Variable> vars;
	vars.reserve(m_varstab->rowCount());

	// iterate all variables
	for(int row = 0; row < m_varstab->rowCount(); ++row)
	{
		auto *name = m_varstab->item(row, COL_VARS_NAME);
		auto *val_re = static_cast<t_numitem*>(m_varstab->item(row, COL_VARS_VALUE_REAL));
		auto *val_im = static_cast<t_numitem*>(m_varstab->item(row, COL_VARS_VALUE_IMAG));

		if(!name || !val_re || !val_im)
		{
			std::cerr << "Invalid entry in variables table row "
				<< row << "." << std::endl;
			continue;
		}

		t_magdyn::Variable var;
		var.name = name->text().toStdString();
		var.value = val_re->GetValue() + val_im->GetValue() * t_cplx(0, 1);

		vars.emplace_back(std::move(var));
	}

	return vars;
}



/**
 * replace numeric values in the tables with variable names
 */
t_size MagDynDlg::ReplaceValuesWithVariables()
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		this_->m_ignoreSitesCalc = false;
	} BOOST_SCOPE_EXIT_END

	m_ignoreCalc = true;
	m_ignoreSitesCalc = true;

	t_size num_replacements = 0;
	for(const t_magdyn::Variable& var : GetVariables())
		num_replacements += ReplaceValueWithVariable(var.name, var.value);

	m_status->setText(QString("Replaced %1 values.").arg(num_replacements));
	return num_replacements;
}



/**
 * replace numeric values in the tables with a variable name
 */
t_size MagDynDlg::ReplaceValueWithVariable(const std::string& var, const t_cplx& val)
{
	// replace the item's numeric value if it's equal to the variable
	auto replace = [&var, &val](t_numitem* item) -> bool
	{
		bool val_ok = false;
		if(bool equ = tl2::equals<t_cplx>(item->GetValue(&val_ok), val, g_eps);
			equ && val_ok)
		{
			// variable is equal to table value
			item->setText(var.c_str());
			return true;
		}
		else if(bool equ = tl2::equals<t_cplx>(item->GetValue(&val_ok), -val, g_eps);
			equ && val_ok)
		{
			// variable is equal to negative table value
			item->setText(("-" + var).c_str());
			return true;
		}

		return false;
	};


	// iterate lines in the exchange terms table and replace the values
	t_size num_replacements = 0;
	for(int row = 0; row < m_termstab->rowCount(); ++row)
	{
		t_numitem *xch = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_INTERACTION));
		t_numitem *dmi_x = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_DMI_X));
		t_numitem *dmi_y = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_DMI_Y));
		t_numitem *dmi_z = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_DMI_Z));

		t_numitem *gen_xx = nullptr, *gen_xy = nullptr, *gen_xz = nullptr;
		t_numitem *gen_yx = nullptr, *gen_yy = nullptr, *gen_yz = nullptr;
		t_numitem *gen_zx = nullptr, *gen_zy = nullptr, *gen_zz = nullptr;
		if(m_allow_general_J)
		{
			gen_xx = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_XX));
			gen_xy = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_XY));
			gen_xz = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_XZ));
			gen_yx = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_YX));
			gen_yy = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_YY));
			gen_yz = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_YZ));
			gen_zx = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_ZX));
			gen_zy = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_ZY));
			gen_zz = static_cast<t_numitem*>(m_termstab->item(row, COL_XCH_GEN_ZZ));
		}

		for(t_numitem *item : { xch, dmi_x, dmi_y, dmi_z,
			gen_xx, gen_xy, gen_xz,
			gen_yx, gen_yy, gen_yz,
			gen_zx, gen_zy, gen_zz })
		{
			if(!item)
				continue;

			if(replace(item))
				++num_replacements;
		}
	}

	return num_replacements;
}

