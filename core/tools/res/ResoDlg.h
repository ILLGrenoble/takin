/**
 * resolution calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013 - 2016
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

#ifndef __RESO_DLG_H__
#define __RESO_DLG_H__

#include <QSettings>
#include <QDialog>
#include <QLabel>

#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include "ui/ui_reso.h"
#include "ellipse.h"

#include "algos/cn.h"
#include "algos/pop.h"
#include "algos/eck.h"
#include "algos/eck_ext.h"
#include "algos/vio.h"
#include "algos/vio_ext.h"
#include "algos/simple.h"

#include "tlibs/math/linalg.h"
#include "tlibs/file/prop.h"
#ifndef NO_3D
	#include "libs/plotgl.h"
#endif
#include "dialogs/RecipParamDlg.h"
#include "dialogs/RealParamDlg.h"
#include "dialogs/EllipseDlg.h"
#include "dialogs/TOFDlg.h"


// parameters that are not already in RealParams or RecipParams
struct ResoParams
{
	bool bMonoDChanged = 0, bAnaDChanged = 0;
	bool bSensesChanged[3] = {0,0,0};

	t_real_reso dMonoD, dAnaD;
	bool bScatterSenses[3];
};


struct SampleParams
{
	t_real_reso dLattice[3];
	t_real_reso dAngles[3];
	t_real_reso dPlane1[3], dPlane2[3];
};


class ResoDlg : public QDialog, Ui::ResoDlg
{Q_OBJECT
private:
	std::string m_strCurDir;
	std::unordered_map<std::string, std::shared_ptr<ReflCurve<t_real_reso>>> m_mapRefl;
	std::shared_ptr<ReflCurve<t_real_reso>> load_cache_refl(const std::string& strFile);

	void DebugOutput();

protected:
	std::vector<QDoubleSpinBox*> m_vecSpinBoxes;
	std::vector<QString> m_vecSpinNames;

	std::vector<QSpinBox*> m_vecIntSpinBoxes;
	std::vector<QString> m_vecIntSpinNames;

	std::vector<QCheckBox*> m_vecCheckBoxes;
	std::vector<QString> m_vecCheckNames;

	std::vector<QLineEdit*> m_vecPosEditBoxes;
	std::vector<QString> m_vecPosEditNames;

	std::vector<QLineEdit*> m_vecEditBoxes;
	std::vector<QString> m_vecEditNames;

	std::vector<QRadioButton*> m_vecRadioPlus;
	std::vector<QRadioButton*> m_vecRadioMinus;
	std::vector<QString> m_vecRadioNames;

	std::vector<QComboBox*> m_vecComboBoxes;
	std::vector<QString> m_vecComboNames;

	void WriteLastConfig();
	void ReadLastConfig();


	// -------------------------------------------------------------------------
	ublas::vector<t_real_reso> m_vecOrient1, m_vecOrient2;
	ublas::matrix<t_real_reso> m_matU, m_matB, m_matUinv, m_matBinv;
	ublas::matrix<t_real_reso> m_matUrlu, m_matUinvrlu;
	ublas::matrix<t_real_reso> m_matUB, m_matUBinv;
	bool m_bHasUB = 0;
	t_real_reso m_dAngleQVec0 = 0.;
	// -------------------------------------------------------------------------


	EckParams m_tasparams;
	VioExtParams m_tofparams;
	SimpleResoParams m_simpleparams;

	ResoResults m_res;
	ublas::matrix<t_real_reso> m_resoHKL, m_resoOrient;
	ublas::vector<t_real_reso> m_reso_vHKL, m_reso_vOrient;
	ublas::vector<t_real_reso> m_Q_avgHKL, m_Q_avgOrient;

	std::vector<ublas::vector<t_real_reso>> m_vecMC_direct;
	std::vector<ublas::vector<t_real_reso>> m_vecMC_HKL;


	bool m_bDontCalc;
	bool m_bEll4dCurrent = 0;
	Ellipsoid4d<t_real_reso> m_ell4d;

	QSettings* m_pSettings = 0;
	bool m_bUpdateOnRealEvent = 1;
	bool m_bUpdateOnRecipEvent = 1;

	std::unique_ptr<TOFDlg> m_pTOFDlg;


	ResoAlgo GetSelectedAlgo() const;
	void SetSelectedAlgo(ResoAlgo algo);

public:
	ResoDlg(QWidget* pParent, QSettings* pSettings=0);
	virtual ~ResoDlg();

	void EmitResults();

protected slots:
	void Calc();
	void AlgoChanged();

	void ShowTOFCalcDlg();

	void SaveRes();
	void LoadRes();

	void LoadMonoRefl();
	void LoadAnaEffic();

	void ButtonBoxClicked(QAbstractButton*);
	void hideEvent (QHideEvent *event);
	void showEvent(QShowEvent *event);

	void checkAutoCalcElli4dChanged();
	void CalcElli4d();
	void MCGenerate();

	void RefreshQEPos();

protected:
	void setupAlgos();
	void RefreshSimCmd();

public slots:
	void ResoParamsChanged(const ResoParams& params);
	void RecipParamsChanged(const RecipParams& parms);
	void RealParamsChanged(const RealParams& parms);
	void SampleParamsChanged(const SampleParams& parms);

public:
	void Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot);
	void Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);

	void SetUpdateOn(bool bRealEvent, bool bRecipEvent)
	{
		m_bUpdateOnRealEvent = bRealEvent;
		m_bUpdateOnRecipEvent = bRecipEvent;
	}

signals:
	void ResoResultsSig(const EllipseDlgParams& params);
};

#endif
