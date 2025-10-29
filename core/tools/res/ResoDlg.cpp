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

#include "ResoDlg.h"
#include <iostream>
#include <fstream>
#include <iomanip>

#include "tlibs/string/string.h"
#include "tlibs/helper/flags.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/helper/misc.h"
#include "tlibs/math/math.h"
#include "tlibs/math/geo.h"
#include "tlibs/math/rand.h"
#include "tlibs/phys/lattice.h"

#include "libs/globals.h"
#include "libs/qt/qthelper.h"
#include "ellipse.h"
#include "mc.h"

#include <QPainter>
#include <QFileDialog>
#include <QMessageBox>
#include <QGridLayout>


using t_mat = ublas::matrix<t_real_reso>;
using t_vec = ublas::vector<t_real_reso>;

static const auto angs = tl::get_one_angstrom<t_real_reso>();
static const auto rads = tl::get_one_radian<t_real_reso>();
static const auto meV = tl::get_one_meV<t_real_reso>();
static const auto cm = tl::get_one_centimeter<t_real_reso>();
static const auto meters = tl::get_one_meter<t_real_reso>();
static const auto sec = tl::get_one_second<t_real_reso>();



ResoDlg::ResoDlg(QWidget *pParent, QSettings* pSettings)
	: QDialog(pParent), m_bDontCalc(true), m_pSettings(pSettings)
{
	setupUi(this);
	spinMCSample->setEnabled(false);      // TODO
	spinMCSampleLive->setEnabled(false);  // TODO

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);
	}

	btnSave->setIcon(load_icon("res/icons/document-save.svg"));
	btnLoad->setIcon(load_icon("res/icons/document-open.svg"));

	setupAlgos();
	connect(comboAlgo, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &ResoDlg::AlgoChanged);
	comboAlgo->setCurrentIndex(3);  // default: pop

	groupGuide->setChecked(false);

	// -------------------------------------------------------------------------
	// widgets
	m_vecSpinBoxes = {spinMonod, spinMonoMosaic, spinAnad,
		spinAnaMosaic, spinSampleMosaic,
		spinHCollMono, spinHCollBSample,
		spinHCollASample, spinHCollAna, spinVCollMono,
		spinVCollBSample, spinVCollASample, spinVCollAna,
		spinMonoRefl, spinAnaEffic,

		spinMonoW, spinMonoH, spinMonoThick, spinMonoCurvH, spinMonoCurvV,
		spinSampleW_Q, spinSampleW_perpQ, spinSampleH,
		spinAnaW, spinAnaH, spinAnaThick, spinAnaCurvH, spinAnaCurvV,
		spinSrcW, spinSrcH,
		spinGuideDivH, spinGuideDivV,
		spinDetW, spinDetH,
		spinDistMonoSample, spinDistSampleAna, spinDistAnaDet,
		spinDistVSrcMono, spinDistHSrcMono,

		spinMonitorW, spinMonitorH, spinMonitorThick,
		spinDistMonoMonitor,
		spinScatterKfAngle,

		spinMonoMosaicV, spinSampleMosaicV, spinAnaMosaicV,
		spinSamplePosX, spinSamplePosY, spinSamplePosZ,

		spinDistTofPulseMono, spinDistTofMonoSample, spinDistTofSampleDet,
		spinDistTofPulseMonoSig, spinDistTofMonoSampleSig, spinDistTofSampleDetSig,
		spinTofPulseSig, spinTofMonoSig, spinTofDetSig,
		spinTof2thI, spinTofphI, spinTofphF,
		spinTof2thISig, spinTof2thFSig, spinTofphISig, spinTofphFSig,

		spinSigKi, spinSigKi_perp, spinSigKi_z,
		spinSigKf, spinSigKf_perp, spinSigKf_z,
	};

	m_vecSpinNames = {"reso/mono_d", "reso/mono_mosaic", "reso/ana_d",
		"reso/ana_mosaic", "reso/sample_mosaic",
		"reso/h_coll_mono", "reso/h_coll_before_sample",
		"reso/h_coll_after_sample", "reso/h_coll_ana", "reso/v_coll_mono",
		"reso/v_coll_before_sample", "reso/v_coll_after_sample", "reso/v_coll_ana",
		"reso/mono_refl", "reso/ana_effic",

		"reso/pop_mono_w", "reso/pop_mono_h", "reso/pop_mono_thick", "reso/pop_mono_curvh", "reso/pop_mono_curvv",
		"reso/pop_sample_wq", "reso/pop_sample_wperpq", "reso/pop_sample_h",
		"reso/pop_ana_w", "reso/pop_ana_h", "reso/pop_ana_thick", "reso/pop_ana_curvh", "reso/pop_ana_curvv",
		"reso/pop_src_w", "reso/pop_src_h",
		"reso/pop_guide_divh", "reso/pop_guide_divv",
		"reso/pop_det_w", "reso/pop_det_h",
		"reso/pop_dist_mono_sample", "reso/pop_dist_sample_ana", "reso/pop_dist_ana_det",
		"reso/pop_dist_vsrc_mono", "reso/pop_dist_hsrc_mono",

		"reso/pop_monitor_w", "reso/pop_monitor_h", "reso/pop_monitor_thick",
		"reso/pop_dist_mono_monitor",
		"reso/scatter_kf_angle",

		"reso/eck_mono_mosaic_v", "reso/eck_sample_mosaic_v", "reso/eck_ana_mosaic_v",
		"reso/eck_sample_pos_x", "reso/eck_sample_pos_y", "reso/eck_sample_pos_z",

		"reso/viol_dist_pulse_mono", "reso/viol_dist_mono_sample", "reso/viol_dist_sample_det",
		"reso/viol_dist_pulse_mono_sig", "reso/viol_dist_mono_sample_sig", "reso/viol_dist_sample_det_sig",
		"reso/viol_time_pulse_sig", "reso/viol_time_mono_sig", "reso/viol_time_det_sig",
		"reso/viol_angle_tt_i", "reso/viol_angle_ph_i", "reso/viol_angle_ph_f",
		"reso/viol_angle_tt_i_sig", "reso/viol_angle_tt_f_sig", "reso/viol_angle_ph_i_sig", "reso/viol_angle_ph_f_sig",

		"reso/simple_sig_ki", "reso/simple_sig_ki_perp", "reso/simple_sig_ki_z",
		"reso/simple_sig_kf", "reso/simple_sig_kf_perp", "reso/simple_sig_kf_z",
	};

	m_vecIntSpinBoxes = { spinMCNeutronsLive, spinMCSampleLive };
	m_vecIntSpinNames = { "reso/mc_live_neutrons", "reso/mc_live_sample_neutrons" };

	m_vecEditBoxes = {editMonoRefl, editAnaEffic};
	m_vecEditNames = {"reso/mono_refl_file", "reso/ana_effic_file"};

	m_vecPosEditBoxes = {editE, editQ, editKi, editKf};
	m_vecPosEditNames = {"reso/E", "reso/Q", "reso/ki", "reso/kf"};

	m_vecCheckBoxes = {checkUseAltR0, checkUseKi3, checkUseKf3,
		checkUseKfKi, checkUseKi, checkUseMonitor, checkUseSampleVol,
		checkUseResVol};
	m_vecCheckNames = {"reso/use_alt_R0", "reso/use_ki3", "reso/use_kf3",
		"reso/use_kfki", "reso/use_monki", "reso/use_mon", "reso/use_samplevol",
		"reso/use_resvol"};

	m_vecRadioPlus = {radioMonoScatterPlus, radioAnaScatterPlus, radioSampleScatterPlus,
		radioSampleCub, radioSrcRect, radioDetRect, radioMonitorRect,
		radioTofDetSph};
	m_vecRadioMinus = {radioMonoScatterMinus, radioAnaScatterMinus, radioSampleScatterMinus,
		radioSampleCyl, radioSrcCirc, radioDetCirc, radioMonitorCirc,
		radioTofDetCyl};
	m_vecRadioNames = {"reso/mono_scatter_sense", "reso/ana_scatter_sense", "reso/sample_scatter_sense",
		"reso/pop_sample_cuboid", "reso/pop_source_rect", "reso/pop_det_rect", "reso/pop_monitor_rect",
		"reso/viol_det_sph"};

	m_vecComboBoxes = {/*comboAlgo,*/
		comboAnaHori, comboAnaVert,
		comboMonoHori, comboMonoVert};
	m_vecComboNames = {/*"reso/algo",*/
		"reso/pop_ana_use_curvh", "reso/pop_ana_use_curvv",
		"reso/pop_mono_use_curvh", "reso/pop_mono_use_curvv"};
	// -------------------------------------------------------------------------

	ReadLastConfig();

	for(QDoubleSpinBox* pSpinBox : m_vecSpinBoxes)
		connect(pSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &ResoDlg::Calc);
	for(QSpinBox* pSpinBox : m_vecIntSpinBoxes)
		connect(pSpinBox, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, &ResoDlg::Calc);
	for(QLineEdit* pEditBox : m_vecEditBoxes)
		connect(pEditBox, &QLineEdit::textChanged, this, &ResoDlg::Calc);
	for(QLineEdit* pEditBox : m_vecPosEditBoxes)
		connect(pEditBox, &QLineEdit::textEdited, this, &ResoDlg::RefreshQEPos);
	for(QRadioButton* pRadio : m_vecRadioPlus)
		connect(pRadio, &QRadioButton::toggled, this, &ResoDlg::Calc);
	for(QCheckBox* pCheck : m_vecCheckBoxes)
		connect(pCheck, &QCheckBox::toggled, this, &ResoDlg::Calc);
	for(QComboBox* pCombo : m_vecComboBoxes)
		connect(pCombo, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &ResoDlg::Calc);

	connect(comboAlgo, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, &ResoDlg::Calc);

	connect(groupGuide, &QGroupBox::toggled, this, &ResoDlg::Calc);
#if QT_VERSION >= QT_VERSION_CHECK(6, 7, 0)
	connect(checkElli4dAutoCalc, &QCheckBox::checkStateChanged, this, &ResoDlg::checkAutoCalcElli4dChanged);
#else
	connect(checkElli4dAutoCalc, &QCheckBox::stateChanged, this, &ResoDlg::checkAutoCalcElli4dChanged);
#endif
	connect(btnCalcElli4d, &QPushButton::clicked, this, &ResoDlg::CalcElli4d);
	connect(btnMCGenerate, &QPushButton::clicked, this, &ResoDlg::MCGenerate);
	connect(buttonBox, &QDialogButtonBox::clicked, this, &ResoDlg::ButtonBoxClicked);
	connect(btnSave, &QPushButton::clicked, this, &ResoDlg::SaveRes);
	connect(btnLoad, &QPushButton::clicked, this, &ResoDlg::LoadRes);
	connect(btnTOFCalc, &QPushButton::clicked, this, &ResoDlg::ShowTOFCalcDlg);
	connect(btnMonoRefl, &QToolButton::clicked, this, &ResoDlg::LoadMonoRefl);
	connect(btnAnaEffic, &QToolButton::clicked, this, &ResoDlg::LoadAnaEffic);

	m_bDontCalc = false;
	RefreshQEPos();
	//Calc();
}



ResoDlg::~ResoDlg()
{ }



void ResoDlg::setupAlgos()
{
	comboAlgo->addItem("TAS: Cooper-Nathans (Pointlike)", static_cast<int>(ResoAlgo::CN));
	comboAlgo->addItem("TAS: Popovici (Pointlike)", static_cast<int>(ResoAlgo::POP_CN));
	comboAlgo->insertSeparator(2);
	comboAlgo->addItem("TAS: Popovici", static_cast<int>(ResoAlgo::POP));
	comboAlgo->addItem("TAS: Eckold-Sobolev", static_cast<int>(ResoAlgo::ECK));
	comboAlgo->addItem("TAS: Eckold-Sobolev (Extended)", static_cast<int>(ResoAlgo::ECK_EXT));
	comboAlgo->insertSeparator(6);
	comboAlgo->addItem("TOF: Violini", static_cast<int>(ResoAlgo::VIO));
	comboAlgo->addItem("TOF: Violini (Extended)", static_cast<int>(ResoAlgo::VIO_EXT));
	comboAlgo->insertSeparator(8);
	comboAlgo->addItem("Simple", static_cast<int>(ResoAlgo::SIMPLE));
}



void ResoDlg::RefreshQEPos()
{
	try
	{
		tl::t_wavenumber_si<t_real_reso> Q = t_real_reso(editQ->text().toDouble()) / angs;
		tl::t_wavenumber_si<t_real_reso> ki = t_real_reso(editKi->text().toDouble()) / angs;
		tl::t_wavenumber_si<t_real_reso> kf = t_real_reso(editKf->text().toDouble()) / angs;
		//t_real_reso dE = editE->text().toDouble();
		tl::t_energy_si<t_real_reso> E = tl::get_energy_transfer(ki, kf);

		tl::t_angle_si<t_real_reso> kiQ = tl::get_angle_ki_Q(ki, kf, Q, true, false);
		tl::t_angle_si<t_real_reso> kfQ = tl::get_angle_kf_Q(ki, kf, Q, true, true);
		tl::t_angle_si<t_real_reso> twotheta = tl::get_sample_twotheta(ki, kf, Q, true);

		const t_real_reso dMono = spinMonod->value();
		const t_real_reso dAna = spinAnad->value();

		m_simpleparams.ki = m_tofparams.ki = m_tasparams.ki = ki;
		m_simpleparams.kf = m_tofparams.kf = m_tasparams.kf = kf;
		m_simpleparams.E = m_tofparams.E = m_tasparams.E = E;
		m_simpleparams.Q = m_tofparams.Q = m_tasparams.Q = Q;

		m_simpleparams.twotheta = m_tofparams.twotheta = m_tasparams.twotheta = twotheta;
		m_simpleparams.angle_ki_Q = m_tofparams.angle_ki_Q = m_tasparams.angle_ki_Q = kiQ;
		m_simpleparams.angle_kf_Q = m_tofparams.angle_kf_Q = m_tasparams.angle_kf_Q = kfQ;

		m_tasparams.thetam = t_real_reso(0.5) * tl::get_mono_twotheta(ki, dMono*angs, true);
		m_tasparams.thetaa = t_real_reso(0.5) * tl::get_mono_twotheta(kf, dAna*angs, true);

#ifndef NDEBUG
		tl::log_debug("Manually changed parameters: ",
			"ki=",ki, ", kf=", kf, ", Q=",Q, ", E=", E,
			", tt=", twotheta, ", kiQ=", kiQ, ", kfQ=", kfQ, ".");
#endif
		Calc();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
		labelStatus->setText(QString("<font color='red'>Error: ") + ex.what() + QString("</font>"));
	}
}



/**
 * loads a reflectivity curve and caches it in a map
 */
std::shared_ptr<ReflCurve<t_real_reso>> ResoDlg::load_cache_refl(const std::string& strFile)
{
	std::shared_ptr<ReflCurve<t_real_reso>> pRefl = nullptr;
	if(strFile == "")
		return pRefl;

	std::vector<std::string> vecRelDirs = { m_strCurDir, "~", "." };
	const std::vector<std::string>& vecGlobalPaths = get_global_paths();
	for(const std::string& strGlobalPath : vecGlobalPaths)
		vecRelDirs.push_back(strGlobalPath);

	auto iter = m_mapRefl.find(strFile);
	if(iter == m_mapRefl.end())
	{ // no yet cached -> load curve
		pRefl = std::make_shared<ReflCurve<t_real_reso>>(strFile, &vecRelDirs);
		if(pRefl && *pRefl)
			m_mapRefl[strFile] = pRefl;
	}
	else
	{ // curve available in cache
		pRefl = iter->second;
	}

	return pRefl;
}



/**
 * calculates the resolution
 */
void ResoDlg::Calc()
{
	try
	{
		m_bEll4dCurrent = false;
		if(m_bDontCalc)
			return;

		EckParams &tas = m_tasparams;
		VioExtParams &tof = m_tofparams;
		SimpleResoParams &simple = m_simpleparams;

		ResoResults &res = m_res;

		// cn parameters
		tas.mono_d = t_real_reso(spinMonod->value()) * angs;
		tas.ana_d = t_real_reso(spinAnad->value()) * angs;
		tas.mono_mosaic = t_real_reso(tl::m2r(spinMonoMosaic->value())) * rads;
		tas.sample_mosaic = t_real_reso(tl::m2r(spinSampleMosaic->value())) * rads;
		tas.ana_mosaic = t_real_reso(tl::m2r(spinAnaMosaic->value())) * rads;
		tas.mono_mosaic_v = t_real_reso(tl::m2r(spinMonoMosaicV->value())) * rads;
		tas.sample_mosaic_v = t_real_reso(tl::m2r(spinSampleMosaicV->value())) * rads;
		tas.ana_mosaic_v = t_real_reso(tl::m2r(spinAnaMosaicV->value())) * rads;

		tas.dmono_sense = (radioMonoScatterPlus->isChecked() ? +1. : -1.);
		tas.dana_sense = (radioAnaScatterPlus->isChecked() ? +1. : -1.);
		tas.dsample_sense = (radioSampleScatterPlus->isChecked() ? +1. : -1.);
		//if(spinQ->value() < 0.)
		//	tas.dsample_sense = -tas.dsample_sense;

		tas.coll_h_pre_mono = t_real_reso(tl::m2r(spinHCollMono->value())) * rads;
		tas.coll_h_pre_sample = t_real_reso(tl::m2r(spinHCollBSample->value())) * rads;
		tas.coll_h_post_sample = t_real_reso(tl::m2r(spinHCollASample->value())) * rads;
		tas.coll_h_post_ana = t_real_reso(tl::m2r(spinHCollAna->value())) * rads;

		tas.coll_v_pre_mono = t_real_reso(tl::m2r(spinVCollMono->value())) * rads;
		tas.coll_v_pre_sample = t_real_reso(tl::m2r(spinVCollBSample->value())) * rads;
		tas.coll_v_post_sample = t_real_reso(tl::m2r(spinVCollASample->value())) * rads;
		tas.coll_v_post_ana = t_real_reso(tl::m2r(spinVCollAna->value())) * rads;

		tas.dmono_refl = spinMonoRefl->value();
		tas.dana_effic = spinAnaEffic->value();
		std::string strMonoRefl = editMonoRefl->text().toStdString();
		std::string strAnaEffic = editAnaEffic->text().toStdString();
		if(strMonoRefl != "")
		{
			tas.mono_refl_curve = load_cache_refl(strMonoRefl);
			if(!tas.mono_refl_curve)
				tl::log_err("Cannot load mono reflectivity file \"", strMonoRefl, "\".");
		}
		if(strAnaEffic != "")
		{
			tas.ana_effic_curve = load_cache_refl(strAnaEffic);
			if(!tas.ana_effic_curve)
				tl::log_err("Cannot load ana reflectivity file \"", strAnaEffic, "\".");
		}


		if(checkUseAltR0->isChecked())
			tas.flags |= CALC_ALT_R0;
		else
			tas.flags &= ~CALC_ALT_R0;
		if(checkUseKi3->isChecked())
			tas.flags |= CALC_KI3;
		else
			tas.flags &= ~CALC_KI3;
		if(checkUseKf3->isChecked())
			tas.flags |= CALC_KF3;
		else
			tas.flags &= ~CALC_KF3;
		if(checkUseKfKi->isChecked())
			tas.flags |= CALC_KFKI;
		else
			tas.flags &= ~CALC_KFKI;
		if(checkUseKi->isChecked())
			tas.flags |= CALC_MONKI;
		else
			tas.flags &= ~CALC_MONKI;
		if(checkUseMonitor->isChecked())
			tas.flags |= CALC_MON;
		else
			tas.flags &= ~CALC_MON;
		if(checkUseSampleVol->isChecked())
			tas.flags |= NORM_TO_SAMPLE;
		else
			tas.flags &= ~NORM_TO_SAMPLE;
		if(checkUseResVol->isChecked())
			tas.flags |= NORM_TO_RESVOL;
		else
			tas.flags &= ~NORM_TO_RESVOL;


		// Position
		/*tas.ki = t_real_reso(editKi->text().toDouble()) / angs;
		tas.kf = t_real_reso(editKf->text().toDouble()) / angs;
		tas.Q = t_real_reso(editQ->text().toDouble()) / angs;
		tas.E = t_real_reso(editE->text().toDouble()) * meV;
		//tas.E = tl::get_energy_transfer(tas.ki, tas.kf);*/


		// pop parameters
		tas.mono_w = t_real_reso(spinMonoW->value()) * cm;
		tas.mono_h = t_real_reso(spinMonoH->value()) * cm;
		tas.mono_thick = t_real_reso(spinMonoThick->value()) * cm;
		tas.mono_curvh = t_real_reso(spinMonoCurvH->value()) * cm;
		tas.mono_curvv = t_real_reso(spinMonoCurvV->value()) * cm;
		tas.bMonoIsCurvedH = tas.bMonoIsCurvedV = false;
		tas.bMonoIsOptimallyCurvedH = tas.bMonoIsOptimallyCurvedV = false;

		spinMonoCurvH->setEnabled(comboMonoHori->currentIndex() == 2);
		spinMonoCurvV->setEnabled(comboMonoVert->currentIndex() == 2);

		if(comboMonoHori->currentIndex() == 2)
			tas.bMonoIsCurvedH = 1;
		else if(comboMonoHori->currentIndex() == 1)
			tas.bMonoIsCurvedH = tas.bMonoIsOptimallyCurvedH = 1;

		if(comboMonoVert->currentIndex() == 2)
			tas.bMonoIsCurvedV = 1;
		else if(comboMonoVert->currentIndex() == 1)
			tas.bMonoIsCurvedV = tas.bMonoIsOptimallyCurvedV = 1;

		tas.ana_w = t_real_reso(spinAnaW->value()) *cm;
		tas.ana_h = t_real_reso(spinAnaH->value()) * cm;
		tas.ana_thick = t_real_reso(spinAnaThick->value()) * cm;
		tas.ana_curvh = t_real_reso(spinAnaCurvH->value()) * cm;
		tas.ana_curvv = t_real_reso(spinAnaCurvV->value()) * cm;
		tas.bAnaIsCurvedH = tas.bAnaIsCurvedV = 0;
		tas.bAnaIsOptimallyCurvedH = tas.bAnaIsOptimallyCurvedV = 0;

		spinAnaCurvH->setEnabled(comboAnaHori->currentIndex() == 2);
		spinAnaCurvV->setEnabled(comboAnaVert->currentIndex() == 2);

		if(comboAnaHori->currentIndex() == 2)
			tas.bAnaIsCurvedH = 1;
		else if(comboAnaHori->currentIndex()==1)
			tas.bAnaIsCurvedH = tas.bAnaIsOptimallyCurvedH = 1;

		if(comboAnaVert->currentIndex() == 2)
			tas.bAnaIsCurvedV = 1;
		else if(comboAnaVert->currentIndex()==1)
			tas.bAnaIsCurvedV = tas.bAnaIsOptimallyCurvedV = 1;

		tas.bSampleCub = radioSampleCub->isChecked();
		tas.sample_w_q = t_real_reso(spinSampleW_Q->value()) * cm;
		tas.sample_w_perpq = t_real_reso(spinSampleW_perpQ->value()) * cm;
		tas.sample_h = t_real_reso(spinSampleH->value()) * cm;

		tas.bSrcRect = radioSrcRect->isChecked();
		tas.src_w = t_real_reso(spinSrcW->value()) * cm;
		tas.src_h = t_real_reso(spinSrcH->value()) * cm;

		tas.bDetRect = radioDetRect->isChecked();
		tas.det_w = t_real_reso(spinDetW->value()) * cm;
		tas.det_h = t_real_reso(spinDetH->value()) * cm;

		tas.bGuide = groupGuide->isChecked();
		tas.guide_div_h = t_real_reso(tl::m2r(spinGuideDivH->value())) * rads;
		tas.guide_div_v = t_real_reso(tl::m2r(spinGuideDivV->value())) * rads;

		tas.dist_mono_sample = t_real_reso(spinDistMonoSample->value()) * cm;
		tas.dist_sample_ana = t_real_reso(spinDistSampleAna->value()) * cm;
		tas.dist_ana_det = t_real_reso(spinDistAnaDet->value()) * cm;
		tas.dist_vsrc_mono = t_real_reso(spinDistVSrcMono->value()) * cm;
		tas.dist_hsrc_mono = t_real_reso(spinDistHSrcMono->value()) * cm;

		tas.bMonitorRect = radioMonitorRect->isChecked();
		tas.monitor_w = t_real_reso(spinMonitorW->value()) * cm;
		tas.monitor_h = t_real_reso(spinMonitorH->value()) * cm;
		tas.monitor_thick = t_real_reso(spinMonitorThick->value()) * cm;
		tas.dist_mono_monitor = t_real_reso(spinDistMonoMonitor->value()) * cm;


		// eck parameters
		tas.pos_x = t_real_reso(spinSamplePosX->value()) * cm;
		tas.pos_y = t_real_reso(spinSamplePosY->value()) * cm;
		tas.pos_z = t_real_reso(spinSamplePosZ->value()) * cm;
		tas.angle_kf = tl::d2r(t_real_reso(spinScatterKfAngle->value())) * rads;

		// TODO
		tas.mono_numtiles_h = 1;
		tas.mono_numtiles_v = 1;
		tas.ana_numtiles_h = 1;
		tas.ana_numtiles_v = 1;


		// TOF parameters
		tof.len_pulse_mono = t_real_reso(spinDistTofPulseMono->value()) * cm;
		tof.len_mono_sample = t_real_reso(spinDistTofMonoSample->value()) * cm;
		tof.len_sample_det = t_real_reso(spinDistTofSampleDet->value()) * cm;

		tof.sig_len_pulse_mono = t_real_reso(spinDistTofPulseMonoSig->value()) * cm;
		tof.sig_len_mono_sample = t_real_reso(spinDistTofMonoSampleSig->value()) * cm;
		tof.sig_len_sample_det = t_real_reso(spinDistTofSampleDetSig->value()) * cm;

		tof.sig_pulse = t_real_reso(spinTofPulseSig->value() * 1e-6) * sec;
		tof.sig_mono = t_real_reso(spinTofMonoSig->value() * 1e-6) * sec;
		tof.sig_det = t_real_reso(spinTofDetSig->value() * 1e-6) * sec;

		tof.twotheta_i = tl::d2r(t_real_reso(spinTof2thI->value())) * rads;
		tof.angle_outplane_i = tl::d2r(t_real_reso(spinTofphI->value())) * rads;
		tof.angle_outplane_f = tl::d2r(t_real_reso(spinTofphF->value())) * rads;

		tof.sig_twotheta_i = tl::d2r(t_real_reso(spinTof2thISig->value())) * rads;
		tof.sig_twotheta_f = tl::d2r(t_real_reso(spinTof2thFSig->value())) * rads;
		tof.sig_outplane_i = tl::d2r(t_real_reso(spinTofphISig->value())) * rads;
		tof.sig_outplane_f = tl::d2r(t_real_reso(spinTofphFSig->value())) * rads;

		tof.det_shape = radioTofDetSph->isChecked() ? TofDetShape::SPH : TofDetShape::CYL;


		// parameters for simple resolution model
		simple.sig_ki = t_real_reso(spinSigKi->value()) / angs;
		simple.sig_kf = t_real_reso(spinSigKf->value()) / angs;
		simple.sig_ki_perp = t_real_reso(spinSigKi_perp->value()) / angs;
		simple.sig_kf_perp = t_real_reso(spinSigKf_perp->value()) / angs;
		simple.sig_ki_z = t_real_reso(spinSigKi_z->value()) / angs;
		simple.sig_kf_z = t_real_reso(spinSigKf_z->value()) / angs;


		// pre-calculate optimal curvature parameters to show in the gui
		tl::t_length_si<t_real_reso> mono_curvh =
			tl::foc_curv(tas.dist_hsrc_mono, tas.dist_mono_sample, tas.ki, tas.mono_d, false);
		tl::t_length_si<t_real_reso> mono_curvv =
			tl::foc_curv(tas.dist_vsrc_mono, tas.dist_mono_sample, tas.ki, tas.mono_d, true);
		tl::t_length_si<t_real_reso> ana_curvh =
			tl::foc_curv(tas.dist_sample_ana, tas.dist_ana_det, tas.kf, tas.ana_d, false);
		tl::t_length_si<t_real_reso> ana_curvv =
			tl::foc_curv(tas.dist_sample_ana, tas.dist_ana_det, tas.kf, tas.ana_d, true);

		std::stringstream ostrCurv;
		ostrCurv.precision(g_iPrecGfx);
		ostrCurv << "Opt. curvatures: "
			<< "Mono.-H.: " << t_real_reso(mono_curvh / cm) << " cm, "
			<< "Mono.-V.: " << t_real_reso(mono_curvv / cm) << " cm, "
			<< "Ana.-H.: " << t_real_reso(ana_curvh / cm) << " cm, "
			<< "Ana.-V.: " << t_real_reso(ana_curvv / cm) << " cm.";
		labelTASGeoInfo->setText(ostrCurv.str().c_str());


		// calculation
		switch(ResoDlg::GetSelectedAlgo())
		{
			case ResoAlgo::CN: res = calc_cn(tas); break;
			case ResoAlgo::POP_CN: res = calc_pop_cn(tas); break;
			case ResoAlgo::POP: res = calc_pop(tas); break;
			case ResoAlgo::ECK: res = calc_eck(tas); break;
			case ResoAlgo::ECK_EXT: res = calc_eck_ext(tas); break;
			case ResoAlgo::VIO: res = calc_vio(tof); break;
			case ResoAlgo::VIO_EXT: res = calc_vio_ext(tof); break;
			case ResoAlgo::SIMPLE: res = calc_simplereso(simple); break;
			default: tl::log_err("Unknown resolution algorithm selected."); return;
		}

		editE->setText(tl::var_to_str(t_real_reso(tas.E/meV), g_iPrec).c_str());
		//if(m_pInstDlg) m_pInstDlg->SetParams(tas, res);
		//if(m_pScatterDlg) m_pScatterDlg->SetParams(tas, res);

		if(res.bOk)
		{
			const std::string strAA_1 = tl::get_spec_char_utf8("AA")
				+ tl::get_spec_char_utf8("sup-")
				+ tl::get_spec_char_utf8("sup1");
			const std::string strAA_3 = tl::get_spec_char_utf8("AA")
				+ tl::get_spec_char_utf8("sup-")
				+ tl::get_spec_char_utf8("sup3");

#ifndef NDEBUG
			// check against ELASTIC approximation for perp. slope from (Shirane 2002), p. 268
			// valid for small mosaicities
			t_real_reso dEoverQperp = tl::co::hbar*tl::co::hbar*tas.ki / tl::co::m_n
				* units::cos(tas.twotheta/2.)
				* (1. + units::tan(units::abs(tas.twotheta/2.))
				* units::tan(units::abs(tas.twotheta/2.) - units::abs(tas.thetam)))
					/ meV / angs;

			tl::log_info("E/Q_perp (approximation for ki=kf) = ", dEoverQperp, " meV*A");
			tl::log_info("E/Q_perp (2nd approximation for ki=kf) = ", t_real_reso(4.*tas.ki * angs), " meV*A");
#endif

			if(checkElli4dAutoCalc->isChecked())
			{
				CalcElli4d();
				m_bEll4dCurrent = true;
			}

			if(groupSim->isChecked())
				RefreshSimCmd();

			// calculate rlu quadric if a sample is defined
			if(m_bHasUB)
			{
				std::tie(m_resoHKL, m_reso_vHKL, m_Q_avgHKL) =
					conv_lab_to_rlu<t_mat, t_vec, t_real_reso>
						(m_dAngleQVec0, m_matUB, m_matUBinv,
						res.reso, res.reso_v, res.Q_avg);
				std::tie(m_resoOrient, m_reso_vOrient, m_Q_avgOrient) =
					conv_lab_to_rlu_orient<t_mat, t_vec, t_real_reso>
						(m_dAngleQVec0, m_matUB, m_matUBinv,
						m_matUrlu, m_matUinvrlu,
						res.reso, res.reso_v, res.Q_avg);
			}

			// print results
			std::ostringstream ostrRes;

			//ostrRes << std::scientific;
			ostrRes.precision(g_iPrec);
			ostrRes << "<html><body>\n";

			ostrRes << "<p><b>Correction Factors:</b>\n";
			ostrRes << "\t<ul><li>Resolution Volume: " << res.dResVol << " meV " << strAA_3 << "</li>\n";
			ostrRes << "\t<li>R0: " << res.dR0 << "</li></ul></p>\n\n";

			// --------------------------------------------------------------------------------
			// Bragg widths
			ostrRes << "<p><b>Coherent (Bragg) FWHMs:</b>\n";
			ostrRes << "\t<ul><li>Q_para: " << res.dBraggFWHMs[0] << " " << strAA_1 << "</li>\n";
			ostrRes << "\t<li>Q_ortho: " << res.dBraggFWHMs[1] << " " << strAA_1 << "</li>\n";
			ostrRes << "\t<li>Q_z: " << res.dBraggFWHMs[2] << " " << strAA_1 << "</li>\n";
			ostrRes << "\t<li>E: " << res.dBraggFWHMs[3] << " meV</li>\n";

			static const char* pcHkl[] = { "h", "k", "l" };
			if(m_bHasUB)
			{
				const std::vector<t_real_reso> vecFwhms = calc_bragg_fwhms(m_resoHKL);

				for(unsigned iHkl=0; iHkl<3; ++iHkl)
				{
					ostrRes << "\t<li>" << pcHkl[iHkl] << ": "
						<< vecFwhms[iHkl] << " rlu</li>\n";
				}
			}

			ostrRes << "</ul></p>\n";
			// --------------------------------------------------------------------------------

			// --------------------------------------------------------------------------------
			// Vanadium widths
			auto dVanadiumFWHMs = calc_vanadium_fwhms(res.reso);
			ostrRes << "<p><b>Incoherent (Vanadium) FWHMs:</b>\n";
			ostrRes << "\t<ul><li>Q_para: " << dVanadiumFWHMs[0] << " " << strAA_1 << "</li>\n";
			ostrRes << "\t<li>Q_ortho: " << dVanadiumFWHMs[1] << " " << strAA_1 << "</li>\n";
			ostrRes << "\t<li>Q_z: " << dVanadiumFWHMs[2] << " " << strAA_1 << "</li>\n";
			ostrRes << "\t<li>E: " << dVanadiumFWHMs[3] << " meV</li>\n";

			if(m_bHasUB)
			{
				const std::vector<t_real_reso> vecFwhms = calc_vanadium_fwhms(m_resoHKL);

				for(unsigned iHkl=0; iHkl<3; ++iHkl)
				{
					ostrRes << "\t<li>" << pcHkl[iHkl] << ": "
						<< vecFwhms[iHkl] << " rlu</li>\n";
				}
			}

			ostrRes << "</ul></p>\n";
			// --------------------------------------------------------------------------------

			ostrRes << "<p><b>Resolution Matrix (Q_para, Q_ortho, Q_z, E) in 1/A, meV and using Gaussian sigmas:</b>\n\n";
			ostrRes << "<blockquote><table border=\"0\" width=\"75%\">\n";
			for(std::size_t i=0; i<res.reso.size1(); ++i)
			{
				ostrRes << "<tr>\n";
				for(std::size_t j=0; j<res.reso.size2(); ++j)
				{
					t_real_reso dVal = res.reso(i,j);
					tl::set_eps_0(dVal, g_dEps);

					ostrRes << "<td>" << std::setw(g_iPrec*2) << dVal << "</td>";
				}
				ostrRes << "</tr>\n";

				if(i!=res.reso.size1()-1)
					ostrRes << "\n";
			}
			ostrRes << "</table></blockquote></p>\n";

			ostrRes << "<p><b>Resolution Vector in 1/A, meV:</b> ";
			for(std::size_t iVec=0; iVec<res.reso_v.size(); ++iVec)
			{
				ostrRes << res.reso_v[iVec];
				if(iVec != res.reso_v.size()-1)
					ostrRes << ", ";
			}
			ostrRes << "<br>\n";
			ostrRes << "<b>Resolution Scalar</b>: " << res.reso_s << "</p>\n";


			{
				tl::Quadric<t_real_reso> quadr(res.reso, res.reso_v, res.reso_s);
				int rank,rankext, pos_evals,neg_evals,zero_evals, pos_evalsext,neg_evalsext,zero_evalsext;
				std::tie(rank,rankext, pos_evals,neg_evals,zero_evals, pos_evalsext,neg_evalsext,zero_evalsext)
					= quadr.ClassifyQuadric(g_dEps);
				ostrRes << "<p><b>Resolution Matrix Rank and Signature:</b> ";
				ostrRes << rank << ", (" << pos_evals << ", " << neg_evals << ", " << zero_evals << ")";
				ostrRes << "<br>\n";
				ostrRes << "<b>Extended Resolution Matrix Rank and Signature:</b> ";
				ostrRes << rankext << ", (" << pos_evalsext << ", " << neg_evalsext << ", " << zero_evalsext << ")";
				ostrRes << "</p>\n";
			}


			if(m_bHasUB)
			{
				ostrRes << "<p><b>Resolution Matrix (h, k, l, E) in rlu, meV:</b>\n\n";
				ostrRes << "<blockquote><table border=\"0\" width=\"75%\">\n";
				for(std::size_t i=0; i<m_resoHKL.size1(); ++i)
				{
					ostrRes << "<tr>\n";
					for(std::size_t j=0; j<m_resoHKL.size2(); ++j)
					{
						t_real_reso dVal = m_resoHKL(i,j);
						tl::set_eps_0(dVal, g_dEps);
						ostrRes << "<td>" << std::setw(g_iPrec*2) << dVal << "</td>";
					}
					ostrRes << "</tr>\n";

					if(i != m_resoHKL.size1()-1)
						ostrRes << "\n";
				}
				ostrRes << "</table></blockquote></p>\n";

				ostrRes << "<p><b>Resolution Vector in rlu, meV:</b> ";
				for(std::size_t iVec=0; iVec<m_reso_vHKL.size(); ++iVec)
				{
					ostrRes << m_reso_vHKL[iVec];
					if(iVec != m_reso_vHKL.size()-1)
						ostrRes << ", ";
				}
				ostrRes << "</p>\n";
				//ostrRes << "<p><b>Resolution Scalar</b>: " << res.reso_s << "</p>\n";

				tl::Quadric<t_real_reso> quadr(m_resoHKL, m_reso_vHKL, res.reso_s);
				int rank,rankext, pos_evals,neg_evals,zero_evals, pos_evalsext,neg_evalsext,zero_evalsext;
				std::tie(rank,rankext, pos_evals,neg_evals,zero_evals, pos_evalsext,neg_evalsext,zero_evalsext)
					= quadr.ClassifyQuadric(g_dEps);
				ostrRes << "<p><b>Resolution Matrix (rlu) Rank and Signature:</b> ";
				ostrRes << rank << ", (" << pos_evals << ", " << neg_evals << ", " << zero_evals << ")";
				ostrRes << "<br>\n";
				ostrRes << "<b>Extended Resolution Matrix (rlu) Rank and Signature:</b> ";
				ostrRes << rankext << ", (" << pos_evalsext << ", " << neg_evalsext << ", " << zero_evalsext << ")";
				ostrRes << "</p>\n";
			}


			ostrRes << "</body></html>";

			editResults->setHtml(QString::fromUtf8(ostrRes.str().c_str()));
			labelStatus->setText("Calculation successful.");



			// generate live MC neutrons
			const std::size_t iNumMC = spinMCNeutronsLive->value();
			if(iNumMC)
			{
				McNeutronOpts<t_mat> opts;
				opts.bCenter = 0;
				opts.matU = m_matU;
				opts.matB = m_matB;
				opts.matUB = m_matUB;
				opts.matUinv = m_matUinv;
				opts.matBinv = m_matBinv;
				opts.matUBinv = m_matUBinv;

				t_mat* pMats[] = {&opts.matU, &opts.matB, &opts.matUB,
					&opts.matUinv, &opts.matBinv, &opts.matUBinv};

				for(t_mat *pMat : pMats)
				{
					pMat->resize(4,4,1);

					for(int i0=0; i0<3; ++i0)
						(*pMat)(i0,3) = (*pMat)(3,i0) = 0.;
					(*pMat)(3,3) = 1.;
				}

				opts.dAngleQVec0 = m_dAngleQVec0;

				if(m_bHasUB)
				{
					// rlu system
					opts.coords = McNeutronCoords::RLU;
					if(m_vecMC_HKL.size() != iNumMC)
						m_vecMC_HKL.resize(iNumMC);
					mc_neutrons<t_vec>(m_ell4d, iNumMC, opts, m_vecMC_HKL.begin());
				}
				else
					m_vecMC_HKL.clear();

				// Qpara, Qperp system
				opts.coords = McNeutronCoords::DIRECT;
				if(m_vecMC_direct.size() != iNumMC)
					m_vecMC_direct.resize(iNumMC);
				mc_neutrons<t_vec>(m_ell4d, iNumMC, opts, m_vecMC_direct.begin());
			}
			else
			{
				m_vecMC_direct.clear();
				m_vecMC_HKL.clear();
			}

			EmitResults();
		}
		else
		{
			labelStatus->setText(QString("<font color='red'>Error: ")
				+ res.strErr.c_str() + QString("</font>"));
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err("Cannot calculate resolution: ", ex.what(), ".");

		labelStatus->setText(QString("<font color='red'>Error: ")
			+ ex.what() + QString("</font>"));
	}
}



void ResoDlg::SetSelectedAlgo(ResoAlgo algo)
{
	for(int iItem = 0; iItem < comboAlgo->count(); ++iItem)
	{
		QVariant varAlgo = comboAlgo->itemData(iItem);
		if(algo == static_cast<ResoAlgo>(varAlgo.toInt()))
		{
			comboAlgo->setCurrentIndex(iItem);
			return;
		}
	}

	tl::log_err("Unknown resolution algorithm set, index: ", static_cast<int>(algo), ".");
}



ResoAlgo ResoDlg::GetSelectedAlgo() const
{
	ResoAlgo algoSel = ResoAlgo::UNKNOWN;
	QVariant varAlgo = comboAlgo->itemData(comboAlgo->currentIndex());
    if(varAlgo == QVariant::Invalid)
		tl::log_err("Unknown resolution algorithm selected, index: ", static_cast<int>(algoSel), ".");
	else
		algoSel = static_cast<ResoAlgo>(varAlgo.toInt());
	return algoSel;
}



void ResoDlg::EmitResults()
{
	ResoAlgo algoSel = ResoDlg::GetSelectedAlgo();
	EllipseDlgParams params;

	params.reso = &m_res.reso;
	params.reso_v = &m_res.reso_v;
	params.reso_s = m_res.reso_s;
	params.Q_avg = &m_res.Q_avg;

	params.resoHKL = &m_resoHKL;
	params.reso_vHKL = &m_reso_vHKL;
	params.Q_avgHKL = &m_Q_avgHKL;

	params.resoOrient = &m_resoOrient;
	params.reso_vOrient = &m_reso_vOrient;
	params.Q_avgOrient = &m_Q_avgOrient;

	params.vecMC_direct = &m_vecMC_direct;
	params.vecMC_HKL = &m_vecMC_HKL;

	params.algo = algoSel;

	emit ResoResultsSig(params);
}



void ResoDlg::ResoParamsChanged(const ResoParams& params)
{
	//tl::log_debug("reso params changed, recalc: ", !m_bDontCalc);

	bool bOldDontCalc = m_bDontCalc;
	m_bDontCalc = true;

	if(params.bSensesChanged[0])
		params.bScatterSenses[0] ? radioMonoScatterPlus->setChecked(1) : radioMonoScatterMinus->setChecked(1);
	if(params.bSensesChanged[1])
		params.bScatterSenses[1] ? radioSampleScatterPlus->setChecked(1) : radioSampleScatterMinus->setChecked(1);
	if(params.bSensesChanged[2])
		params.bScatterSenses[2] ? radioAnaScatterPlus->setChecked(1) : radioAnaScatterMinus->setChecked(1);

	if(params.bMonoDChanged) spinMonod->setValue(params.dMonoD);
	if(params.bAnaDChanged) spinAnad->setValue(params.dAnaD);

	m_bDontCalc = bOldDontCalc;

	// need to recalculate the angles in case the d-spacings have changed
	if(params.bMonoDChanged || params.bAnaDChanged)
		RefreshQEPos();
	Calc();
}



void ResoDlg::RecipParamsChanged(const RecipParams& parms)
{
	//tl::log_debug("recip params changed");

	bool bOldDontCalc = m_bDontCalc;
	m_bDontCalc = true;

	try
	{
		m_simpleparams.twotheta = m_tofparams.twotheta = m_tasparams.twotheta =
			t_real_reso(parms.d2Theta) * rads;

		m_simpleparams.ki = m_tofparams.ki = m_tasparams.ki = t_real_reso(parms.dki) / angs;
		m_simpleparams.kf = m_tofparams.kf = m_tasparams.kf = t_real_reso(parms.dkf) / angs;
		m_simpleparams.E = m_tofparams.E = m_tasparams.E = t_real_reso(parms.dE) * meV;

		t_vec vecHKL = -tl::make_vec({parms.Q_rlu[0], parms.Q_rlu[1], parms.Q_rlu[2]});
		t_real_reso dQ = parms.dQ;

		if(m_bHasUB)
		{
			if(m_matUB.size1() != vecHKL.size())
				vecHKL.resize(m_matUB.size1(), true);
			t_vec vecQ = ublas::prod(m_matUB, vecHKL);
			vecQ.resize(2,1);
			m_dAngleQVec0 = -tl::vec_angle(vecQ);
			dQ = ublas::norm_2(vecQ);
		}

		m_simpleparams.Q = m_tofparams.Q = m_tasparams.Q = dQ / angs;

		m_simpleparams.angle_ki_Q = m_tofparams.angle_ki_Q = m_tasparams.angle_ki_Q =
			tl::get_angle_ki_Q(m_tasparams.ki, m_tasparams.kf, m_tasparams.Q, true, false);
		m_simpleparams.angle_kf_Q = m_tofparams.angle_kf_Q = m_tasparams.angle_kf_Q =
			tl::get_angle_kf_Q(m_tasparams.ki, m_tasparams.kf, m_tasparams.Q, true, true);

		editQ->setText(tl::var_to_str(dQ, g_iPrec).c_str());
		editE->setText(tl::var_to_str(parms.dE, g_iPrec).c_str());
		editKi->setText(tl::var_to_str(parms.dki, g_iPrec).c_str());
		editKf->setText(tl::var_to_str(parms.dkf, g_iPrec).c_str());
	}
	catch(const std::exception& ex)
	{
		tl::log_err("Cannot set reciprocal parameters for resolution: ", ex.what(), ".");
	}

	m_bDontCalc = bOldDontCalc;
	if(m_bUpdateOnRecipEvent)
		Calc();
}



void ResoDlg::RealParamsChanged(const RealParams& parms)
{
	//tl::log_debug("real params changed");

	bool bOldDontCalc = m_bDontCalc;
	m_bDontCalc = true;

	m_tasparams.thetam = units::abs(t_real_reso(parms.dMonoT) * rads);
	m_tasparams.thetaa = units::abs(t_real_reso(parms.dAnaT) * rads);

	m_simpleparams.twotheta = m_tofparams.twotheta = m_tasparams.twotheta =
		t_real_reso(parms.dSampleTT) * rads;

	m_bDontCalc = bOldDontCalc;
	if(m_bUpdateOnRealEvent)
		Calc();
}



void ResoDlg::SampleParamsChanged(const SampleParams& parms)
{
	try
	{
		//tl::log_debug("sample params changed");

		tl::Lattice<t_real_reso> lattice(parms.dLattice[0],parms.dLattice[1],parms.dLattice[2],
			parms.dAngles[0],parms.dAngles[1],parms.dAngles[2]);

		m_vecOrient1 = tl::make_vec<t_vec>({parms.dPlane1[0], parms.dPlane1[1], parms.dPlane1[2]});
		m_vecOrient2 = tl::make_vec<t_vec>({parms.dPlane2[0], parms.dPlane2[1], parms.dPlane2[2]});
		//m_vecOrient1 /= ublas::norm_2(m_vecOrient1);
		//m_vecOrient2 /= ublas::norm_2(m_vecOrient2);

		m_matB = tl::get_B(lattice, 1);
		m_matU = tl::get_U(m_vecOrient1, m_vecOrient2, &m_matB);
		m_matUrlu = tl::get_U(m_vecOrient1, m_vecOrient2);
		m_matUB = ublas::prod(m_matU, m_matB);

		bool bHasB = tl::inverse(m_matB, m_matBinv);
		bool bHasU = tl::inverse(m_matU, m_matUinv);
		bool bHasUrlu = tl::inverse(m_matUrlu, m_matUinvrlu);
		m_matUBinv = ublas::prod(m_matBinv, m_matUinv);

		for(auto* pmat : {&m_matB, &m_matU, &m_matUB, &m_matUBinv, &m_matUrlu, &m_matUinvrlu})
		{
			pmat->resize(4,4,1);
			(*pmat)(3,0) = (*pmat)(3,1) = (*pmat)(3,2) = 0.;
			(*pmat)(0,3) = (*pmat)(1,3) = (*pmat)(2,3) = 0.;
			(*pmat)(3,3) = 1.;
		}

		m_bHasUB = bHasB && bHasU && bHasUrlu;
	}
	catch(const std::exception& ex)
	{
		m_bHasUB = false;
		tl::log_err("Cannot set sample parameters for resolution: ", ex.what(), ".");
	}
}



// --------------------------------------------------------------------------------
// Monte-Carlo stuff

void ResoDlg::checkAutoCalcElli4dChanged()
{
	if(checkElli4dAutoCalc->isChecked() && !m_bEll4dCurrent)
		CalcElli4d();
}



void ResoDlg::CalcElli4d()
{
	m_ell4d = calc_res_ellipsoid4d<t_real_reso>(
		m_res.reso, m_res.reso_v, m_res.reso_s, m_res.Q_avg);

	std::ostringstream ostrElli;
	ostrElli << "<html><body>\n";

	ostrElli << "<p><b>Ellipsoid volume:</b> " << m_ell4d.vol << "</p>\n\n";

	ostrElli << "<p><b>Ellipsoid offsets:</b>\n"
		<< "\t<ul><li>Qx = " << m_ell4d.x_offs << "</li>\n"
		<< "\t<li>Qy = " << m_ell4d.y_offs << "</li>\n"
		<< "\t<li>Qz = " << m_ell4d.z_offs << "</li>\n"
		<< "\t<li>E = " << m_ell4d.w_offs << "</li></ul></p>\n\n";

	ostrElli << "<p><b>Ellipsoid HWHMs (unsorted):</b>\n"
		<< "\t<ul><li>" << m_ell4d.x_hwhm << "</li>\n"
		<< "\t<li>" << m_ell4d.y_hwhm << "</li>\n"
		<< "\t<li>" << m_ell4d.z_hwhm << "</li>\n"
		<< "\t<li>" << m_ell4d.w_hwhm << "</li></ul></p>\n\n";

	ostrElli << "</body></html>\n";

	editElli->setHtml(QString::fromUtf8(ostrElli.str().c_str()));
}



void ResoDlg::MCGenerate()
{
	if(!m_bEll4dCurrent)
		CalcElli4d();

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strLastDir = m_pSettings ? m_pSettings->value("reso/mc_dir", "~").toString() : "~";
	QString _strFile = QFileDialog::getSaveFileName(this, "Save MC neutron data...",
		strLastDir, "Data files (*.dat *.DAT);;All files (*.*)", nullptr, fileopt);
	if(_strFile == "")
		return;

	std::string strFile = _strFile.toStdString();

	const int iNeutrons = spinMCNeutrons->value();
	const bool bCenter = checkMCCenter->isChecked();

	std::ofstream ofstr(strFile);
	if(!ofstr.is_open())
	{
		QMessageBox::critical(this, "Error", "Cannot open file.");
		return;
	}

	std::vector<t_vec> vecNeutrons;
	McNeutronOpts<t_mat> opts;
	opts.bCenter = bCenter;
	opts.coords = McNeutronCoords(comboMCCoords->currentIndex());
	opts.matU = m_matU;
	opts.matB = m_matB;
	opts.matUB = m_matUB;
	opts.matUinv = m_matUinv;
	opts.matBinv = m_matBinv;
	opts.matUBinv = m_matUBinv;

	t_mat* pMats[] = {&opts.matU, &opts.matB, &opts.matUB,
		&opts.matUinv, &opts.matBinv, &opts.matUBinv};

	for(t_mat *pMat : pMats)
	{
		pMat->resize(4,4,1);

		for(int i0=0; i0<3; ++i0)
			(*pMat)(i0,3) = (*pMat)(3,i0) = 0.;
		(*pMat)(3,3) = 1.;
	}


	opts.dAngleQVec0 = m_dAngleQVec0;
	vecNeutrons.resize(iNeutrons);
	mc_neutrons<t_vec>(m_ell4d, iNeutrons, opts, vecNeutrons.begin());


	ofstr.precision(g_iPrec);

	if(opts.coords == McNeutronCoords::DIRECT)
	{
		ofstr << "# coord_sys: direct\n";
		ofstr << "# " << std::setw(std::max<int>(g_iPrec*2-2, 4)) << m_ell4d.x_lab << " "
			<< std::setw(g_iPrec*2) << m_ell4d.y_lab << " "
			<< std::setw(g_iPrec*2) << m_ell4d.z_lab << " "
			<< std::setw(g_iPrec*2) << m_ell4d.w_lab << " \n";
	}
	else if(opts.coords == McNeutronCoords::ANGS)
	{
		ofstr << "# coord_sys: angstrom\n";
		ofstr << "# " << std::setw(std::max<int>(g_iPrec*2-2, 4)) << "Qx (1/A) "
			<< std::setw(g_iPrec*2) << "Qy (1/A) "
			<< std::setw(g_iPrec*2) << "Qz (1/A) "
			<< std::setw(g_iPrec*2) << "E (meV) " << "\n";
	}
	else if(opts.coords == McNeutronCoords::RLU)
	{
		ofstr << "# coord_sys: rlu\n";
		ofstr << "# " << std::setw(std::max<int>(g_iPrec*2-2, 4)) << "h (rlu) "
			<< std::setw(g_iPrec*2) << "k (rlu) "
			<< std::setw(g_iPrec*2) << "l (rlu) "
			<< std::setw(g_iPrec*2) << "E (meV) " << "\n";
	}
	else
	{
		ofstr << "# coord_sys: unknown\n";
	}


	for(const t_vec& vecNeutron : vecNeutrons)
	{
		for(unsigned i = 0; i < 4; ++i)
			ofstr << std::setw(g_iPrec*2) << vecNeutron[i] << " ";
		ofstr << "\n";
	}

	if(m_pSettings)
		m_pSettings->setValue("reso/mc_dir", QString(tl::get_dir(strFile).c_str()));
}



// --------------------------------------------------------------------------------



void ResoDlg::AlgoChanged()
{
	std::string strAlgo = "<html><body>\n";

	switch(GetSelectedAlgo())
	{
		case ResoAlgo::CN:
		case ResoAlgo::POP_CN:
		{
			tabWidget->setTabEnabled(0,1);
			tabWidget->setTabEnabled(1,0);
			tabWidget->setTabEnabled(2,0);
			tabWidget->setTabEnabled(3,0);
			tabWidget->setTabEnabled(4,0);

			strAlgo = "<b>M. J. Cooper and <br>R. Nathans</b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1107/S0365110X67002816>"
				"Acta Cryst. 23, <br>pp. 357-367</a>,<br>\n";
			strAlgo += "1967.";

			strAlgo += "<br><br><b>P. W. Mitchell <i>et al.</i></b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1107/S0108767384000325>"
				"Acta Cryst. A 40(2), <br>pp. 152-160</a>,<br>\n";
			strAlgo += "1984.";
			break;
		}
		case ResoAlgo::POP:
		{
			tabWidget->setTabEnabled(0,1);
			tabWidget->setTabEnabled(1,1);
			tabWidget->setTabEnabled(2,1);
			tabWidget->setTabEnabled(3,0);
			tabWidget->setTabEnabled(4,0);

			strAlgo = "<b>M. Popovici</b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1107/S0567739475001088>"
				"Acta Cryst. A 31, <br>pp. 507-513</a>,<br>\n";
			strAlgo += "1975.";
			break;
		}
		case ResoAlgo::ECK:
		{
			tabWidget->setTabEnabled(0,1);
			tabWidget->setTabEnabled(1,1);
			tabWidget->setTabEnabled(2,1);
			tabWidget->setTabEnabled(3,0);
			tabWidget->setTabEnabled(4,0);

			strAlgo = "<b>G. Eckold and <br>O. Sobolev</b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1016/j.nima.2014.03.019>"
				"NIM A 752, <br>pp. 54-64</a>,<br>\n";
			strAlgo += "2014.";

			strAlgo += "<br><br><b>G. Eckold</b>,<br>\n";
			strAlgo += "personal communication,<br>\n";
			strAlgo += "2020.";
			break;
		}
		case ResoAlgo::ECK_EXT:
		{
			tabWidget->setTabEnabled(0,1);
			tabWidget->setTabEnabled(1,1);
			tabWidget->setTabEnabled(2,1);
			tabWidget->setTabEnabled(3,0);
			tabWidget->setTabEnabled(4,0);

			strAlgo = "<b>G. Eckold and <br>O. Sobolev</b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1016/j.nima.2014.03.019>"
				"NIM A 752, <br>pp. 54-64</a>,<br>\n";
			strAlgo += "2014.";

			strAlgo += "<br><br><b>G. Eckold</b>,<br>\n";
			strAlgo += "personal communication,<br>\n";
			strAlgo += "2020.";

			strAlgo += "<br><br><b>M. Enderle</b>,<br>\n";
			strAlgo += "personal communication,<br>\n";
			strAlgo += "2025.";
			break;
		}
		case ResoAlgo::VIO:
		{
			tabWidget->setTabEnabled(0,0);
			tabWidget->setTabEnabled(1,0);
			tabWidget->setTabEnabled(2,0);
			tabWidget->setTabEnabled(3,1);
			tabWidget->setTabEnabled(4,0);

			strAlgo = "<b>N. Violini <i>et al.</i></b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1016/j.nima.2013.10.042>"
				"NIM A 736, <br>pp. 31-39</a>,<br>\n";
			strAlgo += "2014.";
			break;
		}
		case ResoAlgo::VIO_EXT:
		{
			tabWidget->setTabEnabled(0,0);
			tabWidget->setTabEnabled(1,0);
			tabWidget->setTabEnabled(2,0);
			tabWidget->setTabEnabled(3,1);
			tabWidget->setTabEnabled(4,0);

			strAlgo = "<b>N. Violini <i>et al.</i></b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1016/j.nima.2013.10.042>"
				"NIM A 736, <br>pp. 31-39</a>,<br>\n";
			strAlgo += "2014.";

			strAlgo += "<br><br><b>V. Mecoli</b>,<br>\n";
			strAlgo += "personal communication,<br>\n";
			strAlgo += "2025.";

			break;
		}
		case ResoAlgo::SIMPLE:
		{
			tabWidget->setTabEnabled(0,0);
			tabWidget->setTabEnabled(1,0);
			tabWidget->setTabEnabled(2,0);
			tabWidget->setTabEnabled(3,0);
			tabWidget->setTabEnabled(4,1);

			strAlgo = "<b>Simple</b><br>\n";
			break;
		}
		default:
		{
			strAlgo += "<i>unknown</i>";
			break;
		}
	}

	strAlgo += "\n</body></html>\n";
	labelAlgoRef->setText(strAlgo.c_str());
	labelAlgoRef->setOpenExternalLinks(1);
}



/**
 * quick hack to scan a variable
 */
void ResoDlg::DebugOutput()
{
	std::ostream &ostr = std::cout;

	for(t_real_reso val : tl::linspace<t_real_reso, t_real_reso, std::vector>(1., 500., 128))
	{
		spinMonoCurvH->setValue(val);
		//spinMonoCurvV->setValue(val);
		//spinAnaCurvH->setValue(val);
		//spinAnaCurvV->setValue(val);

		//Calc();
		ostr << val << "\t" << m_res.dR0 << "\t" << m_res.dResVol << std::endl;
	}
}



// --------------------------------------------------------------------------------



void ResoDlg::ShowTOFCalcDlg()
{
	if(!m_pTOFDlg)
		m_pTOFDlg.reset(new TOFDlg(this, m_pSettings));

	focus_dlg(m_pTOFDlg.get());
}



// --------------------------------------------------------------------------------



void ResoDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
	   buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		//DebugOutput();
		WriteLastConfig();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		QDialog::accept();
	}
}



void ResoDlg::hideEvent(QHideEvent *event)
{
	if(m_pSettings)
		m_pSettings->setValue("reso/wnd_geo", saveGeometry());
}



void ResoDlg::showEvent(QShowEvent *event)
{
	if(m_pSettings)
		restoreGeometry(m_pSettings->value("reso/wnd_geo").toByteArray());
}



#include "moc_ResoDlg.cpp"
