/**
 * scan viewer
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2020
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

#ifndef __TAZ_SCANVIEWER_H__
#define __TAZ_SCANVIEWER_H__

#include "ui/ui_scanviewer.h"

#include <QDialog>
#include <QSettings>
#include <QFileSystemWatcher>
#include <QKeyEvent>
#include <string>
#include <vector>
#include <memory>

#include "tlibs/file/loadinstr.h"
#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "libs/globals.h"
#include "FitParamDlg.h"



/**
 * normalise detector (y) by monitor (m) counts
 * y_new = y / m
 * dy_new = 1/m dy - y/m^2 dm
 */
template<class t_real = t_real_glob>
std::pair<t_real, t_real> norm_cnts_to_mon(t_real y, t_real dy, t_real m, t_real dm)
{
	t_real val = y / m;
	t_real err = std::sqrt(std::pow(dy/m, 2.) + std::pow(dm*y/(m*m), 2.));

	return std::make_pair(val, err);
}



class ScanViewerDlg : public QDialog, Ui::ScanViewerDlg
{ Q_OBJECT
public:
	ScanViewerDlg(QWidget* pParent = nullptr, QSettings *core_settings = nullptr);
	virtual ~ScanViewerDlg();

	void SelectDir(const QString& path);


private:
	QSettings m_settings, *m_core_settings{nullptr};

	std::unique_ptr<QFileSystemWatcher> m_pWatcher;
	std::string m_strCurDir, m_strCurFile;
	std::string m_strSelectedKey;
	std::vector<std::string> m_vecExts;

	bool m_bDoUpdate = false;
	tl::FileInstrBase<t_real_glob> *m_pInstr = nullptr;
	std::vector<t_real_glob> m_vecX, m_vecY, m_vecYErr;
	std::vector<t_real_glob> m_vecFitX, m_vecFitY;
	std::unique_ptr<QwtPlotWrapper> m_plotwrap;
	std::string m_strX, m_strY, m_strMon, m_strCmd;

	FitParamDlg *m_pFitParamDlg = nullptr;


private:
	void SetAbout();


protected:
	void ClearPlot();
	void PlotScan();
	void ShowProps();

	int HasRecentPath(const QString& strPath);

	virtual void closeEvent(QCloseEvent* pEvt) override;
	virtual void keyPressEvent(QKeyEvent* pEvt) override;

	template<std::size_t iNumArgs, class t_func>
	bool Fit(t_func&& func,
		const std::vector<std::string>& vecParamNames,
		std::vector<t_real_glob>& vecVals,
		std::vector<t_real_glob>& vecErrs,
		const std::vector<bool>& vecFixed);

	void ShowRawFiles(const std::vector<std::string>& files);


protected slots:
	void GenerateExternal(int iLang=0);

	void UpdateFileList();
	void FileSelected();
	void PropSelected(QTableWidgetItem *pItem, QTableWidgetItem *pItemPrev);
	void SelectDir();
	void ChangedPath();
	void DirWasModified();
	void SearchProps(const QString&);

	void XAxisSelected(int);
	void YAxisSelected(int);
	void MonAxisSelected(int);
	void NormaliseStateChanged(int iState);
	void StartOrSkipChanged(int);

	void ShowFitParams();
	void FitGauss();
	void FitLorentz();
	void FitVoigt();
	void FitLine();
	void FitParabola();
	void FitSine();

	void CalcPol();
};


#endif
