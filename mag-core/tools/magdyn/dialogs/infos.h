/**
 * magnon dynamics -- info dialog
 * @author Tobias Weber <tweber@ill.fr>
 * @date june-2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

#ifndef __MAGDYN_INFOS_H__
#define __MAGDYN_INFOS_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>



class InfoDlg : public QDialog
{ Q_OBJECT
public:
	InfoDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~InfoDlg() = default;

	InfoDlg(const InfoDlg&) = delete;
	const InfoDlg& operator=(const InfoDlg&) = delete;

	void SetGlDeviceInfos(const std::string& ver, const std::string& shader_ver,
		const std::string& vendor, const std::string& renderer);


private:
	QSettings *m_sett{};
	QLabel *m_labelGlInfos[4]{nullptr, nullptr, nullptr, nullptr};


protected slots:
	virtual void accept() override;
};


#endif
