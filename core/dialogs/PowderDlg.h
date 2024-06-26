/**
 * Powder Line Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013, 2-dec-2014
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

#ifndef __POWDER_DLG_H__
#define __POWDER_DLG_H__

#include <QDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QSettings>
#include "ui/ui_powder.h"

#include <map>
#include <string>
#include <vector>
#include <memory>

#include "libs/spacegroups/spacegroup.h"
#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "tlibs/file/prop.h"
#include "AtomsDlg.h"
#include "RecipParamDlg.h"


#define POWDER_MAX_CURVES 64


struct PowderLine
{
	int h, k, l;

	t_real_glob dAngle;
	std::string strAngle;

	t_real_glob dQ;
	std::string strQ;

	std::string strPeaks;

	unsigned int iMult;
	t_real_glob dFn, dFx;	// neutron/xray structure factors
	t_real_glob dIn, dIx;	// neutron/xray intensities
};


class PowderDlg : public QDialog, Ui::PowderDlg
{ Q_OBJECT
	protected:
		// neutron powder lines
		std::vector<t_real_glob> m_vecTT, m_vecInt;
		std::unique_ptr<QwtPlotWrapper> m_plotwrapN;

		// x-ray powder lines
		std::vector<t_real_glob> m_vecTTx, m_vecIntx;
		std::unique_ptr<QwtPlotWrapper> m_plotwrapX;

		// powder angles vs ki plot
		std::vector<std::vector<t_real_glob>> m_vecKis;
		std::vector<std::vector<t_real_glob>> m_vecAngles;
		std::unique_ptr<QwtPlotWrapper> m_plotwrapAnglesKi;

	protected:
		bool m_bDontCalc = 1;
		QSettings *m_pSettings = 0;

		xtl::CrystalSystem m_crystalsys = xtl::CrystalSystem::CRYS_NOT_SET;
		const xtl::SpaceGroups<t_real_glob>::t_mapSpaceGroups* m_pmapSpaceGroups = nullptr;

		AtomsDlg *m_pAtomsDlg = nullptr;
		std::vector<xtl::AtomPos<t_real_glob>> m_vecAtoms;

		t_real_glob m_dExtKi = 0.;
		t_real_glob m_dExtKf = 0.;

	public:
		PowderDlg(QWidget* pParent=0, QSettings* pSett=0);
		virtual ~PowderDlg();

	protected:
		void PlotPowderLines(const std::vector<const PowderLine*>& vecLines);
		void ClearPlots();

	protected slots:
		void CalcPeaks();

		void CheckCrystalType();
		void SpaceGroupChanged();
		void RepopulateSpaceGroups();

		void SaveTable();
		void SavePowder();
		void LoadPowder();

		void ShowAtomDlg();
		void ApplyAtoms(const std::vector<xtl::AtomPos<t_real_glob>>&);

		void cursorMoved(const QPointF& pt);

		virtual void dragEnterEvent(QDragEnterEvent *pEvt) override;
		virtual void dropEvent(QDropEvent *pEvt) override;

		void SetExtKi();
		void SetExtKf();

	public slots:
		void paramsChanged(const RecipParams& parms);

	protected:
		virtual void showEvent(QShowEvent *pEvt) override;
		virtual void accept() override;

		const xtl::SpaceGroup<t_real_glob>* GetCurSpaceGroup() const;

		void Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);
		void Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot);
};

#endif
