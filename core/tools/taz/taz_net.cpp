/**
 * TAS tool (server stuff)
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2014
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

#include "taz.h"

#include <QStatusBar>
#include <QMessageBox>

#define DEFAULT_MSG_TIMEOUT 4000


void TazDlg::ShowConnectDlg()
{
	if(!m_pSrvDlg)
	{
		m_pSrvDlg = new SrvDlg(this, &m_settings);
		QObject::connect(m_pSrvDlg, &SrvDlg::ConnectTo, this, &TazDlg::ConnectTo);
	}

	focus_dlg(m_pSrvDlg);
}


void TazDlg::ConnectTo(ControlSystem control_sys,
	const QString& _strHost, const QString& _strPort,
	const QString& _strUser, const QString& _strPass)
{
	Disconnect();

	std::string strHost = _strHost.toStdString();
	std::string strPort = _strPort.toStdString();
	std::string strUser = _strUser.toStdString();
	std::string strPass = _strPass.toStdString();

	if(control_sys == ControlSystem::NICOS)
	{
		tl::log_info("Trying to establish NICOS connection to ", strHost, ":", strPort, "...");
		m_pNetCache = new NicosCache(&m_settings);
	}
	else if(control_sys == ControlSystem::SICS)
	{
		tl::log_info("Trying to establish SICS connection to ", strHost, ":", strPort, "...");
		m_pNetCache = new SicsCache(&m_settings);
	}
	else
	{
		tl::log_err("Unknown instrument control system selected.");
		return;
	}


	QObject::connect(m_pNetCache, &NetCache::vars_changed, this, &TazDlg::VarsChanged);
	QObject::connect(m_pNetCache, &NetCache::connected, this, &TazDlg::Connected);
	QObject::connect(m_pNetCache, &NetCache::disconnected, this, &TazDlg::Disconnected);

	if(!m_pNetCacheDlg)
		m_pNetCacheDlg = new NetCacheDlg(this, &m_settings);
	if(!m_pScanMonDlg)
		m_pScanMonDlg = new ScanMonDlg(this, &m_settings);

	m_pNetCacheDlg->ClearAll();
	m_pScanMonDlg->ClearPlot();

	QObject::connect(m_pNetCache, &NetCache::updated_cache_value, m_pNetCacheDlg, &NetCacheDlg::UpdateValue);
	QObject::connect(m_pNetCache, &NetCache::updated_cache_value, m_pScanMonDlg, &ScanMonDlg::UpdateValue);


	// no manual node movement
	if(m_sceneReal.GetTasLayout()) m_sceneReal.GetTasLayout()->AllowMouseMove(0);
	if(m_sceneTof.GetTofLayout()) m_sceneTof.GetTofLayout()->AllowMouseMove(0);
	if(m_sceneRecip.GetTriangle()) m_sceneRecip.GetTriangle()->AllowMouseMove(0);

	m_pNetCache->connect(strHost, strPort, strUser, strPass);
}


void TazDlg::Disconnect()
{
	if(m_pNetCache)
	{
		m_pNetCache->disconnect();

		QObject::disconnect(m_pNetCache, &NetCache::vars_changed, this, &TazDlg::VarsChanged);
		QObject::disconnect(m_pNetCache, &NetCache::connected, this, &TazDlg::Connected);
		QObject::disconnect(m_pNetCache, &NetCache::disconnected, this, &TazDlg::Disconnected);

		QObject::disconnect(m_pNetCache, &NetCache::updated_cache_value, m_pNetCacheDlg, &NetCacheDlg::UpdateValue);

		delete m_pNetCache;
		m_pNetCache = nullptr;
	}

	// re-enable manual node movement
	if(m_sceneReal.GetTasLayout()) m_sceneReal.GetTasLayout()->AllowMouseMove(1);
	if(m_sceneTof.GetTofLayout()) m_sceneTof.GetTofLayout()->AllowMouseMove(1);
	if(m_sceneRecip.GetTriangle()) m_sceneRecip.GetTriangle()->AllowMouseMove(1);

	statusBar()->showMessage("Disconnected.", DEFAULT_MSG_TIMEOUT);
}


void TazDlg::ShowNetCache()
{
	if(!m_pNetCacheDlg)
		m_pNetCacheDlg = new NetCacheDlg(this, &m_settings);

	focus_dlg(m_pNetCacheDlg);
}


void TazDlg::ShowNetScanMonitor()
{
	if(!m_pScanMonDlg)
		m_pScanMonDlg = new ScanMonDlg(this, &m_settings);

	focus_dlg(m_pScanMonDlg);
}


void TazDlg::NetRefresh()
{
	if(m_pNetCache)
		m_pNetCache->refresh();
	else
		QMessageBox::warning(this, "Warning", "Not connected to an instrument server.");
}


void TazDlg::Connected(const QString& strHost, const QString& strSrv)
{
	m_strCurFile = "";
	m_vecAtoms.clear();

	setWindowTitle((s_strTitle + " - ").c_str() + strHost + ":" + strSrv);
	statusBar()->showMessage("Connected to " + strHost + " on port " + strSrv + ".", DEFAULT_MSG_TIMEOUT);
}


void TazDlg::Disconnected()
{
	setWindowTitle((s_strTitle).c_str());
}
