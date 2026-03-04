/**
 * dark angles dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jun-2017, 28-jul-2022
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __DARK_ANGLES_H__
#define __DARK_ANGLES_H__

#include <QDialog>
#include <QSettings>

#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "tlibs/file/prop.h"

#include "ui/ui_darkangles.h"


template<class T = double>
struct DarkAngle
{
	T dAngleStart;
	T dAngleEnd;
	T dAngleOffs;

	int iCentreOn;    // 0=mono, 1=sample, 2=ana
	int iRelativeTo;  // 0=crystal angle, 1=in axis, 2=out axis
};


class DarkAnglesDlg : public QDialog, Ui::DarkAnglesDlg
{ Q_OBJECT
protected:
	QSettings *m_pSettings = nullptr;

protected:
	virtual void closeEvent(QCloseEvent*) override;
	void SendApplyDarkAngles();

	// list
	void ClearList();
	void AddAnglesToList(const std::vector<DarkAngle<t_real_glob>>& angles);
	void SetAnglesFromList(QListWidgetItem* item);

protected slots:
	void ButtonBoxClicked(QAbstractButton* pBtn);
	void RemoveAngle();
	void AddAngle();

	// list
	void AddAnglesToList();
	void RemAnglesFromList();
	void LoadList();
	void SaveList();
	void ListItemSelected();
	void ListItemDoubleClicked(QListWidgetItem*);

public:
	DarkAnglesDlg(QWidget* pParent = nullptr, QSettings *pSettings = nullptr);
	virtual ~DarkAnglesDlg();

	void SetDarkAngles(const std::vector<DarkAngle<t_real_glob>>& vecAngle);
	std::vector<DarkAngle<t_real_glob>> GetDarkAngles() const;

	void Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);
	void Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot);

signals:
	void ApplyDarkAngles(const std::vector<DarkAngle<t_real_glob>>& vecAngle);
};

#endif
