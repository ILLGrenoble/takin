/**
 * monte carlo convolution tool
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

#ifndef __MCONVO_GUI_H__
#define __MCONVO_GUI_H__

#include <QDialog>
#include <QSettings>
#include <QMenuBar>
#include <QMenu>

#include <thread>
#include <atomic>
#include <memory>
#include <map>
#include <unordered_map>

#include "ui/ui_monteconvo.h"

#include "libs/qt/qthelper.h"
#include "libs/qt/qwthelper.h"
#include "tlibs/file/prop.h"

#include "sqwfactory.h"
#include "monteconvo_common.h"

#include "tools/res/defs.h"

#include "dialogs/FavDlg.h"
#include "SqwParamDlg.h"
#include "TASReso.h"


#define CONVO_MAX_CURVES          32
#define CONVO_DISP_CURVE_START    3



class ConvoDlg : public QDialog, Ui::ConvoDlg
{ Q_OBJECT
protected:
	std::thread *m_pth = nullptr;
	std::atomic<bool> m_atStop;

	QSettings *m_pSett = nullptr;
	QMenuBar *m_pMenuBar = nullptr;
	SqwParamDlg *m_pSqwParamDlg = nullptr;
	FavDlg *m_pFavDlg = nullptr;

	bool m_bAllowSqwReinit = true;
	std::shared_ptr<SqwBase> m_pSqw;
	std::vector<t_real_reso> m_vecQ, m_vecS, m_vecScaledS;
	std::vector<std::vector<t_real_reso>> m_vecvecQ, m_vecvecE, m_vecvecW;
	std::unique_ptr<QwtPlotWrapper> m_plotwrap, m_plotwrap2d;

	bool m_bUseScan = false;
	Scan m_scan;

	static const std::string s_strTitle;
	std::string m_strLastFile;

	// distance in 1/A to be considered inside the scattering plane
	t_real_reso m_eps_plane = EPS_PLANE;
	t_real_reso m_eps_rlu = EPS_RLU;

	t_real_reso m_chi2 = 0.;

	// maps module identifiers to the corresponding module help texts
	std::unordered_map<std::string, std::string> m_help_texts;


protected:
	std::vector<QDoubleSpinBox*> m_vecSpinBoxes;
	std::vector<QSpinBox*> m_vecIntSpinBoxes;
	std::vector<QLineEdit*> m_vecEditBoxes;
	std::vector<QPlainTextEdit*> m_vecTextBoxes;
	std::vector<QComboBox*> m_vecComboBoxes;
	std::vector<QCheckBox*> m_vecCheckBoxes;

	std::vector<std::string> m_vecSpinNames, m_vecIntSpinNames,
		m_vecEditNames, m_vecTextNames,
		m_vecComboNames, m_vecCheckNames;

	QAction *m_pLiveResults = nullptr, *m_pLivePlots = nullptr;

	// recent files
	QMenu *m_pMenuRecent = nullptr;


protected:
	void LoadSettings();
	virtual void showEvent(QShowEvent *pEvt) override;

	void ClearPlot1D();
	void Start1D();
	void Start2D();


public:
	void Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot);
	void Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);

	t_real_reso GetChi2() const { return m_chi2; }

	void SetSqwParam(const std::string& name, t_real_reso val);

	// [ ident, value, error ]
	void SetSqwParams(const std::vector<std::tuple<
		std::string, std::string, std::string>>& sqwparams);


	// [ ident, type, value, error, fit?, range ]
	using t_sqwparams = std::vector<std::tuple<
		std::string, std::string, std::string, std::string, bool, std::string>>;
	t_sqwparams GetSqwParams(bool only_fitparams) const;


	void StartSim1D(bool bForceDeferred = false, unsigned int seed = tl::get_rand_seed());

	bool StopRequested() const
	{
		return m_atStop.load();
	}


protected slots:
	void showSqwParamDlg();
	void showSqwHelpDlg();

	void browseCrysFiles();
	void browseResoFiles();
	void browseSqwFiles();
	void browseScanFiles();
	void browseAutosaveFile();

	void SqwModelChanged(int);
	void createSqwModel(const QString& qstrFile);
	void SqwParamsChanged(const std::vector<SqwBase::t_var>&,
		const std::vector<SqwBase::t_var_fit>*);

	void scanFileChanged(const QString& qstrFile);
	void scanCheckToggled(bool);
	void coordFlipToggled(bool);
	void scaleChanged();

	void SaveResult(const QString* outfile = nullptr);

	void Start();       // convolution
	void StartFit();    // convolution fit
	void StartDisp();   // plot dispersion
	void Stop();        // stop convo

	void ChangeHK();
	void ChangeHL();
	void ChangeKL();

	void ShowFavourites();
	void UpdateCurFavPos();
	void ChangePos(const struct FavHklPos& pos);

	virtual void accept() override;

	void New();
	void Load(const QString&);
	void Load();
	void Save();
	void SaveAs();
	void SaveConvofit();

	void ShowAboutDlg();


public:
	ConvoDlg(QWidget* pParent=nullptr, QSettings* pSett=nullptr);
	virtual ~ConvoDlg();


signals:
	void SqwLoaded(const std::vector<SqwBase::t_var>&, const std::vector<SqwBase::t_var_fit>*);
};


#endif
