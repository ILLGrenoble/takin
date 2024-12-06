/**
 * qt helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2016
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

#ifndef __QT_HELPER_H__
#define __QT_HELPER_H__

#include <vector>
#include <string>
#include <type_traits>

#include <QTableWidget>
#include <QTableWidgetItem>
#include <QDialog>

#include "tlibs/string/string.h"


extern bool save_table(const char* pcFile, const QTableWidget* pTable);


// a table widget item with a stored numeric value
template<typename T=double>
class QTableWidgetItemWrapper : public QTableWidgetItem
{
protected:
	T m_tVal = T(0);
	unsigned int m_iPrec = std::numeric_limits<T>::digits10;

public:
	QTableWidgetItemWrapper() : QTableWidgetItem(), m_tVal()
	{}

	// same value and item text
	QTableWidgetItemWrapper(T tVal)
		: QTableWidgetItem(tl::var_to_str<T>(tVal, std::numeric_limits<T>::digits10).c_str()),
		m_tVal(tVal)
	{}

	// one value, but different text
	QTableWidgetItemWrapper(T tVal, const std::string& strText)
		: QTableWidgetItem(strText.c_str()),
		m_tVal(tVal)
	{}

	virtual ~QTableWidgetItemWrapper() = default;

	unsigned int GetPrec() const
	{
		return m_iPrec;
	}

	void SetPrec(unsigned int iPrec)
	{
		m_iPrec = iPrec;
	}

	T GetValue() const
	{
		return m_tVal;
	}

	// same value and item text
	void SetValue(T val)
	{
		m_tVal = val;
		this->setText(tl::var_to_str<T>(m_tVal, m_iPrec).c_str());
	}

	// one value, but different text
	void SetValue(T val, const std::string& str)
	{
		m_tVal = val;
		this->setText(str.c_str());
	}

	void SetValueKeepText(T val)
	{
		m_tVal = val;
	}

	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		const QTableWidgetItemWrapper<T>* pItem =
			dynamic_cast<const QTableWidgetItemWrapper<T>*>(&item);
		if(!pItem)
			return false;

		return this->GetValue() < pItem->GetValue();
	}
};


// ----------------------------------------------------------------------------


enum class QtStdPath
{
	HOME,
	FONTS
};

extern std::vector<std::string> get_qt_std_path(QtStdPath path);


// ----------------------------------------------------------------------------


/**
 * focus a dialog
 */
extern void focus_dlg(QDialog* pDlg);


#endif
