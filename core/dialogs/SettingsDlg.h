/**
 * settings
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 5-dec-2014
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __TAZ_SETTINGS_H__
#define __TAZ_SETTINGS_H__

#include <QDialog>
#include <QSettings>

#include <tuple>
#include <vector>
#include <string>

#include "ui/ui_settings.h"


class SettingsDlg : public QDialog, Ui::SettingsDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = 0;

		// key, default, lineedit
		typedef std::tuple<std::string, std::string, QLineEdit*> t_tupEdit;
		std::vector<t_tupEdit> m_vecEdits;

		// checkboxes
		typedef std::tuple<std::string, bool, QCheckBox*> t_tupCheck;
		std::vector<t_tupCheck> m_vecChecks;

		// spins
		typedef std::tuple<std::string, int, QSpinBox*> t_tupSpin;
		std::vector<t_tupSpin> m_vecSpins;

		// combos
		typedef std::tuple<std::string, int, QComboBox*> t_tupCombo;
		std::vector<t_tupCombo> m_vecCombos;

		// tables
		std::vector<int> m_elli_rows, m_ello_rows;

	public:
		SettingsDlg(QWidget* pParent = nullptr, QSettings* pSett = nullptr);
		virtual ~SettingsDlg();

	signals:
		void SettingsChanged() const;

	protected:
		bool LoadSettings(QSettings* = nullptr);
		bool SaveSettings(QSettings* = nullptr) const;

		void SetDefaults(bool bOverwrite = false);
		void SetGlobals(QSettings* = nullptr) const;

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;

	protected slots:
		void ButtonBoxClicked(QAbstractButton* pBtn);

		void SelectGLFont();
		void SelectGfxFont();
		void SelectGenFont();

		void SelectGplTool();
};

#endif
