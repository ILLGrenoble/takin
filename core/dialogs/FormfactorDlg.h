/**
 * Form Factor & Scattering Length Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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

#ifndef __XTL_FF_DLG_H__
#define __XTL_FF_DLG_H__

#include <QDialog>
#include <QSettings>

#include <vector>
#include <memory>

#include "ui/ui_formfactors.h"

#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "libs/globals.h"


class FormfactorDlg : public QDialog, Ui::FormFactorDlg
{ Q_OBJECT
private:
	QListWidgetItem *m_pCurAtom = nullptr;
	QListWidgetItem *m_pCurMagAtom = nullptr;


protected:
	QSettings *m_pSettings = nullptr;

	// form factors
	std::vector<t_real_glob> m_vecQ, m_vecFF;
	std::unique_ptr<QwtPlotWrapper> m_plotwrap;

	// mag form factors
	std::vector<t_real_glob> m_vecQ_m, m_vecFF_m;
	std::unique_ptr<QwtPlotWrapper> m_plotwrap_m;

	// scattering lengths
	std::vector<t_real_glob> m_vecElem, m_vecSc;
	std::unique_ptr<QwtPlotWrapper> m_plotwrapSc;

	// elements
	std::vector<t_real_glob> m_vecElemX, m_vecElemY;
	std::unique_ptr<QwtPlotWrapper> m_plotwrapElem;


protected:
	virtual void closeEvent(QCloseEvent* pEvt) override;
	virtual void accept() override;

	void SetupAtoms();
	void SetupMagAtoms();

	t_real_glob GetFormFact(t_real_glob dQ) const;
	t_real_glob GetMagFormFact(t_real_glob dQ) const;

protected slots:
	void SearchAtom(const QString& qstr);
	void AtomSelected(QListWidgetItem *pItem, QListWidgetItem *pItemPrev);
	void RefreshAtom();

	void SearchMagAtom(const QString& qstr);
	void MagAtomSelected(QListWidgetItem *pItem, QListWidgetItem *pItemPrev);
	void RefreshMagAtom();
	void CalcTermSymbol(const QString& qstr);
	void Calcg();

	void SearchSLAtom(const QString& qstr);
	void SearchElement(const QString& qstr);

	void PlotScatteringLengths();
	void PlotElements();

	void SetupScatteringLengths();
	void SetupElements();

	void cursorMoved(const QPointF& pt);

public:
	FormfactorDlg(QWidget* pParent = nullptr, QSettings *pSettings = nullptr);
	virtual ~FormfactorDlg();
};

#endif
