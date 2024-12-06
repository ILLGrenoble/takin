/**
 * recent files helper
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2015
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

#ifndef __RECENT_FILES_H__
#define __RECENT_FILES_H__

#include <QSettings>
#include <QString>
#include <QStringList>
#include <QMenu>
#include <QAction>
#include <QFile>

#include <algorithm>
#include <iostream>



/**
 * conversion to QStringList
 */
template<class t_cont, class t_str = typename t_cont::value_type>
QStringList cont_to_qstrlist(const t_cont& cont)
{
	QStringList lst;
	for(const t_str& str : cont)
		lst.push_back(str.c_str());
	return lst;
}


/**
 * conversion from QStringList
 */
template<class t_cont, class t_str = typename t_cont::value_type>
t_cont qstrlist_to_cont(const QStringList& lstStr)
{
	t_cont cont;
	for(int i = 0; i < lstStr.size(); ++i)
		cont.push_back(lstStr[i].toStdString());
	return cont;
}


/**
 * conversion from QStringList with a predicate
 */
template<class t_cont, class t_str = typename t_cont::value_type, class t_predicate>
t_cont qstrlist_to_cont(const QStringList& lstStr, t_predicate func)
{
	t_cont cont;
	for(int i = 0; i < lstStr.size(); ++i)
	{
		if(!func(lstStr[i]))
			continue;
		cont.push_back(lstStr[i].toStdString());
	}
	return cont;
}


class RecentFiles
{
	private:
		QSettings *m_pSettings = nullptr;
		std::string m_strKey;

		typedef std::list<std::string> t_lstFiles;
		t_lstFiles m_lstFiles;

		std::size_t m_iMaxFiles = 16;
		bool m_bOnlyExistingFiles = true;


	public:
		RecentFiles(QSettings *pSett, const char* pcKey)
			: m_pSettings(pSett), m_strKey(pcKey)
		{
			LoadList();
		}

		virtual ~RecentFiles() = default;


		void SetOnlyExistingFiles(bool b)
		{
			m_bOnlyExistingFiles = b;
		}


		void LoadList()
		{
			QString strKey(m_strKey.c_str());
			QStringList lstStr = m_pSettings->value(strKey).value<QStringList>();
			lstStr.removeDuplicates();

			m_lstFiles = qstrlist_to_cont<t_lstFiles>(lstStr,
				[this](const QString& file) -> bool
			{
				if(m_bOnlyExistingFiles && !QFile::exists(file))
					return false;
				return true;
			});

			if(m_lstFiles.size() > m_iMaxFiles)
				m_lstFiles.resize(m_iMaxFiles);
		}


		void SaveList()
		{
			QString strKey(m_strKey.c_str());
			QStringList lstStr = cont_to_qstrlist<t_lstFiles>(m_lstFiles);
			m_pSettings->setValue(strKey, lstStr);
		}


		bool HasFile(const char* pcFile) const
		{
			t_lstFiles::const_iterator iter = std::find(m_lstFiles.begin(), m_lstFiles.end(), std::string(pcFile));
			return (iter != m_lstFiles.end());
		}


		void AddFile(const char* pcFile)
		{
			if(m_bOnlyExistingFiles && !QFile::exists(pcFile))
				return;

			RemoveFile(pcFile);

			m_lstFiles.push_front(std::string(pcFile));
			if(m_lstFiles.size() > m_iMaxFiles)
				m_lstFiles.resize(m_iMaxFiles);
		}


		void RemoveFile(const char* pcFile)
		{
			t_lstFiles::iterator iter = std::remove(m_lstFiles.begin(), m_lstFiles.end(), std::string(pcFile));
			if(iter != m_lstFiles.end())
				m_lstFiles.resize(m_lstFiles.size()-1);
		}


		void SetMaxSize(std::size_t iMaxFiles)
		{
			m_iMaxFiles = iMaxFiles;
			if(m_lstFiles.size() > iMaxFiles)
				m_lstFiles.resize(iMaxFiles);
		}


		template<class t_func>
		void FillMenu(QMenu* pMenu, t_func func)
		{
			// clear old actions
			QList<QAction*> lstActions = pMenu->actions();
			for(int iAction=0; iAction<lstActions.size(); ++iAction)
			{
				delete lstActions[iAction];
				lstActions[iAction] = nullptr;
			}
			pMenu->clear();

			// setup new actions
			for(const t_lstFiles::value_type& str : m_lstFiles)
			{
				if(m_bOnlyExistingFiles && !QFile::exists(str.c_str()))
					continue;

				QAction *pAction = new QAction(pMenu);
				QObject::connect(pAction, &QAction::triggered, [str, func]() -> void
				{
					func(str.c_str());
				});

				pAction->setText(str.c_str());
				pMenu->addAction(pAction);
			}
		}
};

#endif
