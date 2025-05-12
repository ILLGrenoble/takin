/**
 * S(Q, E) parameters dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2015
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

#ifndef __SQW_DLG_H__
#define __SQW_DLG_H__

#include <QDialog>
#include <QSettings>

#include "ui/ui_sqwparams.h"
#include "sqwbase.h"



class SqwParamDlg : public QDialog, Ui::SqwParamDlg
{ Q_OBJECT
protected:
	QSettings *m_pSett = nullptr;

protected:
	void SaveSqwParams();
	virtual void showEvent(QShowEvent *pEvt) override;

public:
	SqwParamDlg(QWidget* pParent=nullptr, QSettings* pSett=nullptr);
	virtual ~SqwParamDlg();

public slots:
	void SqwLoaded(const std::vector<SqwBase::t_var>&, const std::vector<SqwBase::t_var_fit>*);

protected slots:
	void ButtonBoxClicked(QAbstractButton *pBtn);

signals:
	void SqwParamsChanged(const std::vector<SqwBase::t_var>&, const std::vector<SqwBase::t_var_fit>*);
};

#endif
