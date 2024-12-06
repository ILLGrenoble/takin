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

#include "qthelper.h"
#include "tlibs/math/math.h"
#include "tlibs/string/string.h"
#include "tlibs/helper/misc.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>

#include <QFileDialog>
#include <QMessageBox>
#include <QMouseEvent>


bool save_table(const char* pcFile, const QTableWidget* pTable)
{
	std::ofstream ofstr(pcFile);
	if(!ofstr)
		return false;

	const int iNumCols = pTable->columnCount();
	const int iNumRows = pTable->rowCount();

	// item lengths
	std::unique_ptr<int[]> ptrMaxTxtLen(new int[iNumCols]);
	for(int iCol = 0; iCol < iNumCols; ++iCol)
		ptrMaxTxtLen[iCol] = 0;

	for(int iCol = 0; iCol < iNumCols; ++iCol)
	{
		const QTableWidgetItem *pItem = pTable->horizontalHeaderItem(iCol);
		ptrMaxTxtLen[iCol] = std::max(pItem ? pItem->text().length() : 0, ptrMaxTxtLen[iCol]);
	}

	for(int iRow = 0; iRow < iNumRows; ++iRow)
	{
		for(int iCol = 0; iCol < iNumCols; ++iCol)
		{
			const QTableWidgetItem *pItem = pTable->item(iRow, iCol);
			ptrMaxTxtLen[iCol] = std::max(pItem ? pItem->text().length() : 0, ptrMaxTxtLen[iCol]);
		}
	}


	// write items
	for(int iCol = 0; iCol < iNumCols; ++iCol)
	{
		const QTableWidgetItem *pItem = pTable->horizontalHeaderItem(iCol);
		ofstr << std::setw(ptrMaxTxtLen[iCol]+4) << (pItem ? pItem->text().toStdString() : "");
	}
	ofstr << "\n";

	for(int iRow = 0; iRow < iNumRows; ++iRow)
	{
		for(int iCol = 0; iCol < iNumCols; ++iCol)
		{
			const QTableWidgetItem *pItem = pTable->item(iRow, iCol);
			ofstr << std::setw(ptrMaxTxtLen[iCol]+4) << (pItem ? pItem->text().toStdString() : "");
		}

		ofstr << "\n";
	}

	return true;
}


// ----------------------------------------------------------------------------


#include <QStandardPaths>

std::vector<std::string> get_qt_std_path(QtStdPath path)
{
	QStandardPaths::StandardLocation iLoc;
	switch(path)
	{
		case QtStdPath::FONTS:
			iLoc = QStandardPaths::FontsLocation;
			break;
		case QtStdPath::HOME:
		default:
			iLoc = QStandardPaths::HomeLocation;
			break;
	}

	QStringList lst = QStandardPaths::standardLocations(iLoc);

	std::vector<std::string> vecPaths;
	for(int iStr=0; iStr<lst.size(); ++iStr)
		vecPaths.push_back(lst.at(iStr).toStdString());
	return vecPaths;
}


// ----------------------------------------------------------------------------


void focus_dlg(QDialog* pDlg)
{
	if(!pDlg)
		return;

	pDlg->show();
	pDlg->raise();
	pDlg->activateWindow();
}
