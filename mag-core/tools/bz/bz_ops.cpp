/**
 * brillouin zone tool
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#include "bz.h"

#include <QtWidgets/QTabWidget>
#include <QtWidgets/QMessageBox>

#include <iostream>
#include <tuple>

#include <boost/scope_exit.hpp>

#include "libs/symops.h"
#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/qt/helper.h"
#include "tlibs2/libs/qt/numerictablewidgetitem.h"


using namespace tl2_ops;


void BZDlg::AddSymOpTabItem(int row, const t_mat& op)
{
	m_symOpIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_symOpIgnoreChanges = false;
		this_->CalcBZ(true);
	} BOOST_SCOPE_EXIT_END


	bool bclone = false;

	if(row == -1)	// append to end of table
		row = m_symops->rowCount();
	else if(row == -2 && m_symOpCursorRow >= 0)	// use row from member variable
		row = m_symOpCursorRow;
	else if(row == -3 && m_symOpCursorRow >= 0)	// use row from member variable +1
		row = m_symOpCursorRow + 1;
	else if(row == -4 && m_symOpCursorRow >= 0)	// use row from member variable +1
	{
		row = m_symOpCursorRow + 1;
		bclone = true;
	}

	m_symops->setSortingEnabled(false);
	m_symops->insertRow(row);

	if(bclone)
	{
		for(int thecol=0; thecol<NUM_SYMOP_COLS; ++thecol)
			m_symops->setItem(row, thecol,
				m_symops->item(m_symOpCursorRow, thecol)->clone());
	}
	else
	{
		std::string prop = get_op_properties<t_mat, t_vec>(op, g_eps);
		t_real det = tl2::det(op);
		//t_real det_rot = tl2::det(tl2::submat(op, 3, 3));

		m_symops->setItem(row, COL_OP,
			new QTableWidgetItem(op_to_str<t_mat>(op, g_prec, g_eps).c_str()));
		m_symops->setItem(row, COL_PROP,
			new QTableWidgetItem(prop.c_str()));
		m_symops->setItem(row, COL_DET,
			new tl2::NumericTableWidgetItem<t_real>(det));
	}

	m_symops->scrollToItem(m_symops->item(row, 0));
	m_symops->setCurrentCell(row, 0);

	m_symops->setSortingEnabled(true);
}


void BZDlg::DelSymOpTabItem(int begin, int end)
{
	m_symOpIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_symOpIgnoreChanges = false;
		this_->CalcBZ(true);
	} BOOST_SCOPE_EXIT_END

	// if nothing is selected, clear all items
	if(begin == -1 || m_symops->selectedItems().count() == 0)
	{
		m_symops->clearContents();
		m_symops->setRowCount(0);
	}
	else if(begin == -2)	// clear selected
	{
		for(int row : GetSelectedSymOpRows(true))
			m_symops->removeRow(row);
	}
	else if(begin >= 0 && end >= 0)		// clear given range
	{
		for(int row=end-1; row>=begin; --row)
			m_symops->removeRow(row);
	}
}


void BZDlg::MoveSymOpTabItemUp()
{
	m_symOpIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_symOpIgnoreChanges = false;
	} BOOST_SCOPE_EXIT_END

	m_symops->setSortingEnabled(false);

	auto selected = GetSelectedSymOpRows(false);
	for(int row : selected)
	{
		if(row == 0)
			continue;

		auto *item = m_symops->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		m_symops->insertRow(row-1);
		for(int col = 0; col < m_symops->columnCount(); ++col)
			m_symops->setItem(row-1, col, m_symops->item(row+1, col)->clone());
		m_symops->removeRow(row+1);
	}

	for(int row = 0; row < m_symops->rowCount(); ++row)
	{
		if(auto *item = m_symops->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row+1)
				!= selected.end())
		{
			for(int col = 0; col < m_symops->columnCount(); ++col)
				m_symops->item(row, col)->setSelected(true);
		}
	}
}


void BZDlg::MoveSymOpTabItemDown()
{
	m_symOpIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_symOpIgnoreChanges = false;
	} BOOST_SCOPE_EXIT_END

	m_symops->setSortingEnabled(false);

	auto selected = GetSelectedSymOpRows(true);
	for(int row : selected)
	{
		if(row == m_symops->rowCount()-1)
			continue;

		auto *item = m_symops->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		m_symops->insertRow(row+2);
		for(int col = 0; col < m_symops->columnCount(); ++col)
			m_symops->setItem(row+2, col, m_symops->item(row, col)->clone());
		m_symops->removeRow(row);
	}

	for(int row = 0; row < m_symops->rowCount(); ++row)
	{
		if(auto *item = m_symops->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row-1)
				!= selected.end())
		{
			for(int col=0; col<m_symops->columnCount(); ++col)
				m_symops->item(row, col)->setSelected(true);
		}
	}
}


std::vector<int> BZDlg::GetSelectedSymOpRows(bool sort_reversed) const
{
	std::vector<int> vec;
	vec.reserve(m_symops->selectedItems().size());

	for(int row = 0; row < m_symops->rowCount(); ++row)
	{
		if(auto *item = m_symops->item(row, 0); item && item->isSelected())
			vec.push_back(row);
	}

	if(sort_reversed)
	{
		std::stable_sort(vec.begin(), vec.end(), [](int row1, int row2)
		{
			return row1 > row2;
		});
	}

	return vec;
}


/**
 * item contents changed
 */
void BZDlg::SymOpTableItemChanged(QTableWidgetItem *item)
{
	// update properties
	if(item->column() == COL_OP)
	{
		t_mat op = str_to_op<t_mat>(item->text().toStdString());
		std::string prop = get_op_properties<t_mat, t_vec>(op, g_eps);
		t_real det = tl2::det(op);
		//t_real det_rot = tl2::det(tl2::submat(op, 3, 3));

		if(QTableWidgetItem *itemProp = m_symops->item(item->row(), COL_PROP); itemProp)
			itemProp->setText(prop.c_str());
		else
			m_symops->setItem(item->row(), COL_PROP, new QTableWidgetItem(prop.c_str()));

		if(tl2::NumericTableWidgetItem<t_real> *itemDet =
			static_cast<tl2::NumericTableWidgetItem<t_real>*>(m_symops->item(item->row(), COL_DET)); itemDet)
			itemDet->SetValue(det);
		else
			m_symops->setItem(item->row(), COL_DET, new tl2::NumericTableWidgetItem<t_real>(det));
	}

	if(!m_symOpIgnoreChanges)
		CalcBZ(true);
}


void BZDlg::ShowSymOpTableContextMenu(const QPoint& _pt)
{
        QPoint pt = _pt;
        // transform the point to widget coordinates first if it has a viewport
        if(m_symops->viewport())
                pt = m_symops->viewport()->mapToParent(pt);

        // transform the point to global coordinates
	auto ptGlob = m_symops->mapToGlobal(pt);

	if(const auto* item = m_symops->itemAt(_pt); item)
	{
		m_symOpCursorRow = item->row();
		m_symOpContextMenu->popup(ptGlob);
	}
	else
	{
		m_symOpContextMenuNoItem->popup(ptGlob);
	}
}


/**
 * get symmetry operations from selected space group
 */
void BZDlg::GetSymOpsFromSG()
{
	m_ignoreCalc = 1;

	try
	{
		// symops of current space group
		auto sgidx = m_comboSG->itemData(m_comboSG->currentIndex()).toInt();
		if(sgidx < 0 || std::size_t(sgidx) >= m_sg_ops.size())
		{
			QMessageBox::critical(this, "Space Group Conversion",
				"Invalid space group selected.");
			m_ignoreCalc = 0;
			return;
		}

		// remove original symops
		DelSymOpTabItem(-1);

		// add symops
		for(const auto& op : m_sg_ops[sgidx])
		{
			AddSymOpTabItem(-1, op);
		}
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Space Group Conversion", ex.what());
	}


	m_ignoreCalc = 0;
	CalcB(true);
}


/**
 * reads symmetry operations from table
 */
std::vector<t_mat> BZDlg::GetSymOps(bool only_centring) const
{
	std::vector<t_mat> ops;

	for(int row = 0; row < m_symops->rowCount(); ++row)
	{
		auto *op_item = m_symops->item(row, COL_OP);
		if(!op_item)
		{
			std::cerr << "Invalid entry in symop table row " << row << "." << std::endl;
			continue;
		}

		t_mat op = str_to_op<t_mat>(op_item->text().toStdString());

		bool add_op = true;
		if(only_centring)
			add_op = tl2::hom_is_centring<t_mat>(op, g_eps);

		if(add_op)
			ops.emplace_back(std::move(op));
	}

	return ops;
}
