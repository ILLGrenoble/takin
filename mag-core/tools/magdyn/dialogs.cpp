/**
 * magnetic dynamics -- dialog setup
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


// instantiate settings dialog
template class SettingsDlg<g_settingsvariables.size(), &g_settingsvariables>;
using t_SettingsDlg = SettingsDlg<g_settingsvariables.size(), &g_settingsvariables>;



/**
 * initialise the static part of the settings dialog
 */
void MagDynDlg::InitSettingsDlg()
{
	// set-up common gui settings variables
	t_SettingsDlg::SetGuiTheme(&g_theme);
	t_SettingsDlg::SetGuiFont(&g_font);
	t_SettingsDlg::SetGuiUseNativeMenubar(&g_use_native_menubar);
	t_SettingsDlg::SetGuiUseNativeDialogs(&g_use_native_dialogs);

	// restore settings
	t_SettingsDlg::ReadSettings(m_sett);
}



/**
 * get changes from settings dialog
 */
void MagDynDlg::InitSettings()
{
	// calculator settings
	m_dyn.SetSilent(g_silent);
	m_dyn.SetPerformChecks(g_checks);
	m_dyn.SetEpsilon(g_eps);
	m_dyn.SetPrecision(g_prec);
	m_dyn.SetBoseCutoffEnergy(g_bose_cutoff);
	m_dyn.SetCholeskyMaxTries(g_cholesky_maxtries);
	m_dyn.SetCholeskyInc(g_cholesky_delta);

	m_recent.SetMaxRecentFiles(g_maxnum_recents);
	m_recent_struct.SetMaxRecentFiles(g_maxnum_recents);

	if(g_font != "")
	{
		QFont font = this->font();
		if(font.fromString(g_font))
			setFont(font);
	}
}



/**
 * settings dialog
 */
void MagDynDlg::ShowSettingsDlg()
{
	if(!m_settings_dlg)
	{
		m_settings_dlg = new t_SettingsDlg(this, m_sett);

		dynamic_cast<t_SettingsDlg*>(m_settings_dlg)->AddChangedSettingsSlot([this]()
		{
			MagDynDlg::InitSettings();
		});
	}

	m_settings_dlg->show();
	m_settings_dlg->raise();
	m_settings_dlg->activateWindow();
}



/**
 * notes dialog
 */
void MagDynDlg::ShowNotesDlg(bool only_create)
{
	if(!m_notes_dlg)
	{
		m_notes_dlg = new NotesDlg(this, m_sett);
		m_notes_dlg->setFont(this->font());
	}

	if(!only_create)
	{
		m_notes_dlg->show();
		m_notes_dlg->raise();
		m_notes_dlg->activateWindow();
	}
}



/**
 * about dialog
 */
void MagDynDlg::ShowInfoDlg(bool only_create)
{
	if(!m_info_dlg)
	{
		m_info_dlg = new InfoDlg(this, m_sett);
		m_info_dlg->setFont(this->font());
	}

	if(!only_create)
	{
		m_info_dlg->show();
		m_info_dlg->raise();
		m_info_dlg->activateWindow();
	}
}



/**
 * structure plotter dialog
 */
void MagDynDlg::ShowStructPlotDlg(bool only_create)
{
	if(!m_structplot_dlg)
	{
		m_structplot_dlg = new StructPlotDlg(this, m_sett);
		m_structplot_dlg->setFont(this->font());

		m_structplot_dlg->SetKernel(&m_dyn);
		m_structplot_dlg->SetTables(m_sitestab, m_termstab);

		connect(m_structplot_dlg, &StructPlotDlg::SelectSite,
			this, &MagDynDlg::SelectSite);
		connect(m_structplot_dlg, &StructPlotDlg::DeleteSite,
			this, &MagDynDlg::DeleteSite);
		connect(m_structplot_dlg, &StructPlotDlg::FlipSiteSpin,
			this, &MagDynDlg::FlipSiteSpin);

		connect(m_structplot_dlg, &StructPlotDlg::SelectTerm,
			this, &MagDynDlg::SelectTerm);
		connect(m_structplot_dlg, &StructPlotDlg::DeleteTerm,
			this, &MagDynDlg::DeleteTerm);

		connect(m_structplot_dlg, &StructPlotDlg::GlDeviceInfos,
			m_info_dlg, &InfoDlg::SetGlDeviceInfos);
	}

	if(!only_create)
	{
		m_structplot_dlg->show();
		m_structplot_dlg->raise();
		m_structplot_dlg->activateWindow();
	}
}



/**
 * ground state minimiser dialog
 */
void MagDynDlg::ShowGroundStateDlg(bool only_create)
{
	if(!m_groundstate_dlg)
	{
		m_groundstate_dlg = new GroundStateDlg(this, m_sett);
		m_groundstate_dlg->setFont(this->font());

		m_groundstate_dlg->SetKernel(&m_dyn);

		connect(m_groundstate_dlg, &GroundStateDlg::SpinsUpdated,
			[this](const t_magdyn* dyn)
		{
			this->SetKernel(dyn, true, false, false);
		});
	}

	if(!only_create)
	{
		m_groundstate_dlg->show();
		m_groundstate_dlg->raise();
		m_groundstate_dlg->activateWindow();
	}
}



/**
 * topology dialog
 */
void MagDynDlg::ShowTopologyDlg(bool only_create)
{
	if(!m_topo_dlg)
	{
		m_topo_dlg = new TopologyDlg(this, m_sett);
		m_topo_dlg->setFont(this->font());

		m_topo_dlg->SetKernel(&m_dyn);

		// set Q position
		auto [Q_start, Q_end] = GetDispersionQ();
		m_topo_dlg->SetDispersionQ(Q_start, Q_end);
	}

	if(!only_create)
	{
		m_topo_dlg->show();
		m_topo_dlg->raise();
		m_topo_dlg->activateWindow();
	}
}



/**
 * 3d dispersion dialog
 */
void MagDynDlg::ShowDispersion3DDlg(bool only_create)
{
	if(!m_disp3d_dlg)
	{
		m_disp3d_dlg = new Dispersion3DDlg(this, m_sett);
		m_disp3d_dlg->setFont(this->font());

		m_disp3d_dlg->SetKernel(&m_dyn);

		// set Q position
		auto [Q_start, Q_end] = GetDispersionQ();
		m_disp3d_dlg->SetDispersionQ(Q_start, Q_end);

		connect(m_disp3d_dlg, &Dispersion3DDlg::GlDeviceInfos,
			m_info_dlg, &InfoDlg::SetGlDeviceInfos);
	}

	if(!only_create)
	{
		m_disp3d_dlg->show();
		m_disp3d_dlg->raise();
		m_disp3d_dlg->activateWindow();
	}
}
