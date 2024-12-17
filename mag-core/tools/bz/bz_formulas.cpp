/**
 * brillouin zone tool -- formulas table
 * @author Tobias Weber <tweber@ill.fr>
 * @date June-2022
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

#include <iostream>
#include <tuple>

#include <boost/scope_exit.hpp>

#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/qt/helper.h"

using namespace tl2_ops;


void BZDlg::AddFormulaTabItem(int row, const std::string& formula)
{
	m_formulaIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_formulaIgnoreChanges = false;
		this_->CalcFormulas();
	} BOOST_SCOPE_EXIT_END

	bool bclone = false;

	if(row == -1)	// append to end of table
	{
		row = m_formulas->rowCount();
	}
	else if(row == -2 && m_formulaCursorRow >= 0)	// use row from member variable
	{
		row = m_formulaCursorRow;
	}
	else if(row == -3 && m_formulaCursorRow >= 0)	// use row from member variable +1
	{
		row = m_formulaCursorRow + 1;
	}
	else if(row == -4 && m_formulaCursorRow >= 0)	// use row from member variable +1
	{
		row = m_formulaCursorRow + 1;
		bclone = true;
	}

	//bool sorting = m_formulas->isSortingEnabled();
	m_formulas->setSortingEnabled(false);
	m_formulas->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < NUM_FORMULAS_COLS; ++thecol)
		{
			m_formulas->setItem(row, thecol,
				m_formulas->item(m_formulaCursorRow, thecol)->clone());
		}
	}
	else
	{
		m_formulas->setItem(row, COL_FORMULA,
			new QTableWidgetItem(formula.c_str()));
	}

	m_formulas->scrollToItem(m_formulas->item(row, 0));
	m_formulas->setCurrentCell(row, 0);

	m_formulas->setSortingEnabled(/*sorting*/ true);
}


void BZDlg::DelFormulaTabItem(int begin, int end)
{
	m_formulaIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_formulaIgnoreChanges = false;
		this_->CalcFormulas();
	} BOOST_SCOPE_EXIT_END

	// if nothing is selected, clear all items
	if(begin == -1 || m_formulas->selectedItems().count() == 0)
	{
		m_formulas->clearContents();
		m_formulas->setRowCount(0);
	}
	else if(begin == -2)	// clear selected
	{
		for(int row : GetSelectedFormulaRows(true))
			m_formulas->removeRow(row);
	}
	else if(begin >= 0 && end >= 0)		// clear given range
	{
		for(int row=end-1; row>=begin; --row)
			m_formulas->removeRow(row);
	}
}


void BZDlg::MoveFormulaTabItemUp()
{
	m_formulaIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_formulaIgnoreChanges = false;
	} BOOST_SCOPE_EXIT_END

	m_formulas->setSortingEnabled(false);

	auto selected = GetSelectedFormulaRows(false);
	for(int row : selected)
	{
		if(row == 0)
			continue;

		auto *item = m_formulas->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		m_formulas->insertRow(row-1);
		for(int col = 0; col < m_formulas->columnCount(); ++col)
		{
			m_formulas->setItem(row - 1, col,
				m_formulas->item(row + 1, col)->clone());
		}
		m_formulas->removeRow(row + 1);
	}

	for(int row = 0; row < m_formulas->rowCount(); ++row)
	{
		if(auto *item = m_formulas->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row + 1)
				!= selected.end())
		{
			for(int col = 0; col < m_formulas->columnCount(); ++col)
				m_formulas->item(row, col)->setSelected(true);
		}
	}
}


void BZDlg::MoveFormulaTabItemDown()
{
	m_formulaIgnoreChanges = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_formulaIgnoreChanges = false;
	} BOOST_SCOPE_EXIT_END

	m_formulas->setSortingEnabled(false);

	auto selected = GetSelectedFormulaRows(true);
	for(int row : selected)
	{
		if(row == m_formulas->rowCount() - 1)
			continue;

		auto *item = m_formulas->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		m_formulas->insertRow(row+2);
		for(int col = 0; col < m_formulas->columnCount(); ++col)
			m_formulas->setItem(row+2, col, m_formulas->item(row, col)->clone());
		m_formulas->removeRow(row);
	}

	for(int row=0; row<m_formulas->rowCount(); ++row)
	{
		if(auto *item = m_formulas->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row - 1)
				!= selected.end())
		{
			for(int col = 0; col < m_formulas->columnCount(); ++col)
				m_formulas->item(row, col)->setSelected(true);
		}
	}
}


std::vector<int> BZDlg::GetSelectedFormulaRows(bool sort_reversed) const
{
	std::vector<int> vec;
	vec.reserve(m_formulas->selectedItems().size());

	for(int row = 0; row < m_formulas->rowCount(); ++row)
	{
		if(auto *item = m_formulas->item(row, 0); item && item->isSelected())
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
void BZDlg::FormulaTableItemChanged([[maybe_unused]] QTableWidgetItem *item)
{
	if(!m_formulaIgnoreChanges)
		CalcFormulas();
}


void BZDlg::ShowFormulaTableContextMenu(const QPoint& _pt)
{
        QPoint pt = _pt;
        // transform the point to widget coordinates first if it has a viewport
        if(m_formulas->viewport())
                pt = m_formulas->viewport()->mapToParent(pt);

        // transform the point to global coordinates
	auto ptGlob = m_formulas->mapToGlobal(pt);

	if(const auto* item = m_formulas->itemAt(_pt); item)
	{
		m_formulaCursorRow = item->row();
		m_formulasContextMenu->popup(ptGlob);
	}
	else
	{
		m_formulasContextMenuNoItem->popup(ptGlob);
	}
}


/**
 * reads all formulas from the table
 */
std::vector<std::string> BZDlg::GetFormulas() const
{
	std::vector<std::string> formulas;

	for(int row = 0; row < m_formulas->rowCount(); ++row)
	{
		auto *op_item = m_formulas->item(row, COL_FORMULA);
		if(!op_item)
		{
			std::cerr << "Invalid entry in formula table row "
				<< row << "." << std::endl;
			continue;
		}

		formulas.emplace_back(op_item->text().toStdString());
	}

	return formulas;
}
