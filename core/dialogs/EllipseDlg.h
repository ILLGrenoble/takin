/**
 * Ellipse Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013 - 2024
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __RESO_ELLI_DLG__
#define __RESO_ELLI_DLG__

#include <QDialog>
#include <QSettings>
#include <QtWidgets>

#include <vector>
#include <memory>

#include "tlibs/math/linalg.h"
#include "tools/res/ellipse.h"
#include "tools/res/defs.h"
#include "ui/ui_ellipses.h"
#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"


struct EllipseDlgParams
{
	const ublas::matrix<t_real_reso>* reso = nullptr;
	const ublas::vector<t_real_reso>* reso_v = nullptr;
	t_real_reso reso_s = 0;
	const ublas::vector<t_real_reso>* Q_avg = nullptr;

	const ublas::matrix<t_real_reso>* resoHKL = nullptr;
	const ublas::vector<t_real_reso>* reso_vHKL = nullptr;
	const ublas::vector<t_real_reso>* Q_avgHKL = nullptr;

	const ublas::matrix<t_real_reso>* resoOrient = nullptr;
	const ublas::vector<t_real_reso>* reso_vOrient = nullptr;
	const ublas::vector<t_real_reso>* Q_avgOrient = nullptr;

	const std::vector<ublas::vector<t_real_reso>>* vecMC_direct = nullptr;
	const std::vector<ublas::vector<t_real_reso>>* vecMC_HKL = nullptr;

	ResoAlgo algo = ResoAlgo::UNKNOWN;
};


class EllipseDlg : public QDialog, Ui::EllipseDlg
{ Q_OBJECT
	protected:
		// plot windows ready?
		bool m_bReady = false;

		std::vector<std::unique_ptr<QwtPlotWrapper>> m_vecplotwrap;

		std::vector<struct Ellipse2d<t_real_reso>> m_elliProj;
		std::vector<struct Ellipse2d<t_real_reso>> m_elliSlice;

		std::vector<std::vector<t_real_reso>> m_vecXCurvePoints;
		std::vector<std::vector<t_real_reso>> m_vecYCurvePoints;

		std::vector<std::vector<t_real_reso>> m_vecMCXCurvePoints;
		std::vector<std::vector<t_real_reso>> m_vecMCYCurvePoints;

		QSettings *m_pSettings = nullptr;

	protected:
		// pointers to original parameters
		EllipseDlgParams m_params;

		// center ellipses around (0, 0)
		bool m_bCenterOn0 = true;

		// copied or modified parameters
		ublas::matrix<t_real_reso> m_reso, m_resoHKL, m_resoOrient;
		ublas::vector<t_real_reso> m_reso_v = ublas::zero_vector<t_real_reso>(4),
			m_reso_vHKL = ublas::zero_vector<t_real_reso>(4),
			m_reso_vOrient = ublas::zero_vector<t_real_reso>(4);
		t_real_reso m_reso_s = 0;
		ublas::vector<t_real_reso> m_Q_avg, m_Q_avgHKL, m_Q_avgOrient;

		ResoAlgo m_algo = ResoAlgo::UNKNOWN;

	public:
		EllipseDlg(QWidget* pParent=nullptr, QSettings* pSett=nullptr, Qt::WindowFlags fl=Qt::Tool);
		virtual ~EllipseDlg();

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void closeEvent(QCloseEvent *pEvt) override;
		virtual void accept() override;

	protected slots:
		void cursorMoved(const QPointF& pt);
		void SaveEllipses();

	public slots:
		void SetParams(const EllipseDlgParams& params);
		void Calc();
		void SetCenterOn0(bool bCenter);

	public:
		void SetTitle(const char* pcTitle);
};

#endif
