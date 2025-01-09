/**
 * numeric table widget item
 * @author Tobias Weber <tweber@ill.fr>
 * @date Dec-2018
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "misc" (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
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

#ifndef __NUM_TABWIDGETITEM_H__
#define __NUM_TABWIDGETITEM_H__


#include <QtWidgets/QTableWidget>
#include "../str.h"


namespace tl2 {


template<class T = double>
class NumericTableWidgetItem : public QTableWidgetItem
{
public:
	/**
	 * construct a numeric table widget item from a numeric value
	 */
	explicit NumericTableWidgetItem(const T& val, std::streamsize prec = 6)
		: QTableWidgetItem(var_to_str(val, prec).c_str()),
		  m_prec{prec}
	{}


	/**
	 * construct a numeric table widget item from a string value
	 */
	explicit NumericTableWidgetItem(const QString& val, std::streamsize prec = 6)
		: QTableWidgetItem(val), m_prec{prec}
	{}


	/**
	 * construct a numeric table widget item from a string value
	 */
	explicit NumericTableWidgetItem(const std::string& val, std::streamsize prec = 6)
		: QTableWidgetItem(QString(val.c_str())), m_prec{prec}
	{}


	/**
	 * compare items
	 */
	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		std::pair<bool, T> res1 = eval_expr<std::string, T>(text().toStdString());
		std::pair<bool, T> res2 = eval_expr<std::string, T>(item.text().toStdString());

		// if the expressions cannot be parsed, order lexicographically
		if(!res1.first || !res2.first)
			return QTableWidgetItem::operator<(item);

		return res1.second < res2.second;
	}


	/**
	 * get a cloned copy of this table item widget
	 */
	virtual QTableWidgetItem* clone() const override
	{
		auto item = new (std::nothrow)
			NumericTableWidgetItem<T>(GetValue(), m_prec);
		if(item)
			*static_cast<QTableWidgetItem*>(item) = *this;

		return item;
	};


	/**
	 * get the numerical value
	 */
	T GetValue(bool *valid = nullptr) const
	{
		std::string str = text().toStdString();
		std::pair<bool, T> pairResult = eval_expr<std::string, T>(str);

		if(valid)
			*valid = pairResult.first;

		// invalid result?
		if(!pairResult.first)
			return T{};

		return pairResult.second;
	}


	/**
	 * set the numerical value
	 */
	void SetValue(T val)
	{
		QString str = var_to_str(val, m_prec).c_str();
		setText(str);
	}


private:
	// precision
	std::streamsize m_prec{6};
};


}
#endif
