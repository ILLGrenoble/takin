/**
 * resolution calculation
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013 - 2016
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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
static const auto sec = tl::get_one_second<t_real_reso>();
static const t_real_reso pi = tl::get_pi<t_real_reso>();
static const t_real_reso fwhm2sig = 1./tl::get_SIGMA2FWHM<t_real_reso>();



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
	m_vecSpinBoxes = {
		// cn
		spinMonod, spinMonoMosaic, spinAnad,
		spinAnaMosaic, spinSampleMosaic,
		spinHCollMono, spinHCollBSample,
		spinHCollASample, spinHCollAna, spinVCollMono,
		spinVCollBSample, spinVCollASample, spinVCollAna,
		spinMonoRefl, spinAnaEffic,

		// pop
		spinMonoW, spinMonoH, spinMonoThick, spinMonoCurvH, spinMonoCurvV,
		spinSampleW_Q, spinSampleW_perpQ, spinSampleH,
		spinAnaW, spinAnaH, spinAnaThick, spinAnaCurvH, spinAnaCurvV,
		spinSrcW, spinSrcH,
		spinGuideDivH, spinGuideDivV,
		spinDetW, spinDetH,
		spinDistMonoSample, spinDistSampleAna, spinDistAnaDet,
		spinDistVSrcMono, spinDistHSrcMono,

		// pop
		spinMonitorW, spinMonitorH, spinMonitorThick,
		spinDistMonoMonitor,
		spinScatterKfAngle,

		// eck
		spinMonoMosaicV, spinSampleMosaicV, spinAnaMosaicV,
		spinSamplePosX, spinSamplePosY, spinSamplePosZ,

		// vio
		spinDistTofPulseMono, spinDistTofMonoSample, spinDistTofSampleDet,
		spinDistTofPulseMonoSig, spinDistTofMonoSampleSig, spinDistTofSampleDetSig,
		spinTofPulseSig, spinTofMonoSig, spinTofDetSig,
		spinTof2thI, spinTofphI, spinTofphF,
		spinTof2thISig, spinTof2thFSig, spinTofphISig, spinTofphFSig,

		// vio_ext
		spinDistTof2PulseGuide, spinDistTof2MonoGuide, spinDistTof2GuideSample, spinDistTof2SampleDet,
		spinTof2PulseWin, spinTof2PulseBeam, spinTof2PulseWidth, spinTof2PulseRPM,
		spinTof2MonoWin, spinTof2MonoBeam, spinTof2MonoWidth, spinTof2MonoRPM,
		spinDistTof2GuideWidth, spinTof2GuideHeight, spinTof2GuideCoating,
		spinTof2DetTubeWidth, spinTof2DetHeight, spinTof2DetZ,
		spinTof2SampleWidth, spinTof2SampleHeight,

		// simple
		spinSigKi, spinSigKi_perp, spinSigKi_z,
		spinSigKf, spinSigKf_perp, spinSigKf_z,
	};

	m_vecSpinNames = {
		// cn
		"reso/mono_d", "reso/mono_mosaic", "reso/ana_d",
		"reso/ana_mosaic", "reso/sample_mosaic",
		"reso/h_coll_mono", "reso/h_coll_before_sample",
		"reso/h_coll_after_sample", "reso/h_coll_ana", "reso/v_coll_mono",
		"reso/v_coll_before_sample", "reso/v_coll_after_sample", "reso/v_coll_ana",
		"reso/mono_refl", "reso/ana_effic",

		// pop
		"reso/pop_mono_w", "reso/pop_mono_h", "reso/pop_mono_thick", "reso/pop_mono_curvh", "reso/pop_mono_curvv",
		"reso/pop_sample_wq", "reso/pop_sample_wperpq", "reso/pop_sample_h",
		"reso/pop_ana_w", "reso/pop_ana_h", "reso/pop_ana_thick", "reso/pop_ana_curvh", "reso/pop_ana_curvv",
		"reso/pop_src_w", "reso/pop_src_h",
		"reso/pop_guide_divh", "reso/pop_guide_divv",
		"reso/pop_det_w", "reso/pop_det_h",
		"reso/pop_dist_mono_sample", "reso/pop_dist_sample_ana", "reso/pop_dist_ana_det",
		"reso/pop_dist_vsrc_mono", "reso/pop_dist_hsrc_mono",

		// pop
		"reso/pop_monitor_w", "reso/pop_monitor_h", "reso/pop_monitor_thick",
		"reso/pop_dist_mono_monitor",
		"reso/scatter_kf_angle",

		// eck
		"reso/eck_mono_mosaic_v", "reso/eck_sample_mosaic_v", "reso/eck_ana_mosaic_v",
		"reso/eck_sample_pos_x", "reso/eck_sample_pos_y", "reso/eck_sample_pos_z",

		// vio
		"reso/viol_dist_pulse_mono", "reso/viol_dist_mono_sample", "reso/viol_dist_sample_det",
		"reso/viol_dist_pulse_mono_sig", "reso/viol_dist_mono_sample_sig", "reso/viol_dist_sample_det_sig",
		"reso/viol_time_pulse_sig", "reso/viol_time_mono_sig", "reso/viol_time_det_sig",
		"reso/viol_angle_tt_i", "reso/viol_angle_ph_i", "reso/viol_angle_ph_f",
		"reso/viol_angle_tt_i_sig", "reso/viol_angle_tt_f_sig", "reso/viol_angle_ph_i_sig", "reso/viol_angle_ph_f_sig",

		// vio_ext
		"reso/vio_ext_dist_pulse_guide", "reso/vio_ext_dist_mono_guide", "reso/vio_ext_dist_guide_sample", "reso/vio_ext_dist_sample_det",
		"reso/vio_ext_pulse_chopper_win", "reso/vio_ext_pulse_chopper_beam", "reso/vio_ext_pulse_chopper_width", "reso/vio_ext_pulse_chopper_rpm",
		"reso/vio_ext_mono_chopper_win", "reso/vio_ext_mono_chopper_beam", "reso/vio_ext_mono_chopper_width", "reso/vio_ext_mono_chopper_rpm",
		"reso/vio_ext_guide_width", "reso/vio_ext_guide_height", "reso/vio_ext_guide_coating",
		"reso/vio_ext_det_tube_width", "reso/vio_ext_det_height", "reso/vio_ext_det_z",
		"reso/vio_ext_sample_width", "reso/vio_ext_sample_height",

		// simple
		"reso/simple_sig_ki", "reso/simple_sig_ki_perp", "reso/simple_sig_ki_z",
		"reso/simple_sig_kf", "reso/simple_sig_kf_perp", "reso/simple_sig_kf_z",
	};

	m_vecIntSpinBoxes = { spinMCNeutronsLive, spinMCSampleLive, spinTof2MCLengths };
	m_vecIntSpinNames = { "reso/mc_live_neutrons", "reso/mc_live_sample_neutrons", "reso/vio_ext_mc_lengths" };

	m_vecEditBoxes = {editMonoRefl, editAnaEffic};
	m_vecEditNames = {"reso/mono_refl_file", "reso/ana_effic_file"};

	m_vecPosEditBoxes = {editE, editQ, editKi, editKf};
	m_vecPosEditNames = {"reso/E", "reso/Q", "reso/ki", "reso/kf"};

	m_vecCheckBoxes = {
		checkUseAltR0, checkUseKi3, checkUseKf3,
		checkUseKfKi, checkUseKi, checkUseMonitor, checkUseSampleVol,
		checkUseResVol,

		checkTof2PulseCounterRot, checkTof2MonoCounterRot,
	};
	m_vecCheckNames = {
		"reso/use_alt_R0", "reso/use_ki3", "reso/use_kf3",
		"reso/use_kfki", "reso/use_monki", "reso/use_mon", "reso/use_samplevol",
		"reso/use_resvol",

		"reso/vio_ext_pulse_chopper_counterrot", "reso/vio_ext_mono_chopper_counterrot",
	};

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
	connect(btnTOFSet, &QPushButton::clicked, this, &ResoDlg::ConvertVioExtParams);
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



// --------------------------------------------------------------------------------



void ResoDlg::AlgoChanged()
{
	std::string strAlgo = "<html><body>\n";

	switch(GetSelectedAlgo())
	{
		case ResoAlgo::CN:
		case ResoAlgo::POP_CN:
		{
			tabWidget->setTabEnabled(0, 1);
			tabWidget->setTabEnabled(1, 0);
			tabWidget->setTabEnabled(2, 0);
			tabWidget->setTabEnabled(3, 0);
			tabWidget->setTabEnabled(4, 0);
			tabWidget->setTabEnabled(5, 0);

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
			tabWidget->setTabEnabled(0, 1);
			tabWidget->setTabEnabled(1, 1);
			tabWidget->setTabEnabled(2, 1);
			tabWidget->setTabEnabled(3, 0);
			tabWidget->setTabEnabled(4, 0);
			tabWidget->setTabEnabled(5, 0);

			strAlgo = "<b>M. Popovici</b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1107/S0567739475001088>"
				"Acta Cryst. A 31, <br>pp. 507-513</a>,<br>\n";
			strAlgo += "1975.";
			break;
		}
		case ResoAlgo::ECK:
		{
			tabWidget->setTabEnabled(0, 1);
			tabWidget->setTabEnabled(1, 1);
			tabWidget->setTabEnabled(2, 1);
			tabWidget->setTabEnabled(3, 0);
			tabWidget->setTabEnabled(4, 0);
			tabWidget->setTabEnabled(5, 0);

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
			tabWidget->setTabEnabled(0, 1);
			tabWidget->setTabEnabled(1, 1);
			tabWidget->setTabEnabled(2, 1);
			tabWidget->setTabEnabled(3, 0);
			tabWidget->setTabEnabled(4, 0);
			tabWidget->setTabEnabled(5, 0);

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
			tabWidget->setTabEnabled(0, 0);
			tabWidget->setTabEnabled(1, 0);
			tabWidget->setTabEnabled(2, 0);
			tabWidget->setTabEnabled(3, 1);
			tabWidget->setTabEnabled(4, 0);
			tabWidget->setTabEnabled(5, 0);

			strAlgo = "<b>N. Violini <i>et al.</i></b>,<br>\n";
			strAlgo += "<a href=http://dx.doi.org/10.1016/j.nima.2013.10.042>"
				"NIM A 736, <br>pp. 31-39</a>,<br>\n";
			strAlgo += "2014.";
			break;
		}
		case ResoAlgo::VIO_EXT:
		{
			tabWidget->setTabEnabled(0, 0);
			tabWidget->setTabEnabled(1, 0);
			tabWidget->setTabEnabled(2, 0);
			tabWidget->setTabEnabled(3, 0);
			tabWidget->setTabEnabled(4, 1);
			tabWidget->setTabEnabled(5, 0);

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
			tabWidget->setTabEnabled(0, 0);
			tabWidget->setTabEnabled(1, 0);
			tabWidget->setTabEnabled(2, 0);
			tabWidget->setTabEnabled(3, 0);
			tabWidget->setTabEnabled(4, 0);
			tabWidget->setTabEnabled(5, 1);

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



/**
 * converts "vio_ext" parameters to "vio" parameters
 * @references based on conversion formulas by V. Mecoli (personal communication, 19-nov-2025)
 */
void ResoDlg::ConvertVioExtParams()
{
	const auto ki = editKi->text().toDouble() / angs;
	const auto kf = editKf->text().toDouble() / angs;
	const t_real_reso mu = 1e-6;
	const t_real_reso vi = tl::k2v(ki) * mu*sec / cm;
	const t_real_reso vf = tl::k2v(kf) * mu*sec / cm;
	const t_real_reso lam = tl::k2lam(ki) / angs;
	const t_real_reso cr_mult_pulse = checkTof2PulseCounterRot->isChecked() ? 2. : 1.;
	const t_real_reso cr_mult_mono = checkTof2MonoCounterRot->isChecked() ? 2. : 1.;
	const t_real_reso Nb = 9.41e-4;  // sld of nickel, see: J. Phys.: Conf. Ser. 528 012005 (2014)

	// distances in cm
	t_real_reso dist_sample_det = spinDistTof2SampleDet->value();
	spinDistTofPulseMono->setValue(spinDistTof2PulseGuide->value() - spinDistTof2MonoGuide->value());
	spinDistTofMonoSample->setValue(spinDistTof2MonoGuide->value() + spinDistTof2GuideSample->value());
	spinDistTofSampleDet->setValue(dist_sample_det);

	// time bases in mus
	const t_real_reso c0 = 1./2./std::sqrt(3.);
	t_real_reso det_tube_w = spinTof2DetTubeWidth->value();
	spinTofPulseSig->setValue(c0 * spinTof2PulseWin->value() / (cr_mult_pulse*spinTof2PulseRPM->value() * 6.*mu) * fwhm2sig);
	spinTofMonoSig->setValue(c0 * spinTof2MonoWin->value() / (cr_mult_mono*spinTof2MonoRPM->value() * 6.*mu) * fwhm2sig);
	spinTofDetSig->setValue(det_tube_w / vf * fwhm2sig);

	// path length deviations in cm
	spinDistTofMonoSampleSig->setValue(vi*spinTofMonoSig->value() - spinTof2MonoWidth->value() * fwhm2sig);
	spinDistTofPulseMonoSig->setValue(std::sqrt(
		  std::pow(vi*spinTofPulseSig->value() - spinTof2PulseWidth->value(), 2.)
		+ std::pow(spinDistTofMonoSampleSig->value(), 2.)) * fwhm2sig);
	spinDistTofSampleDetSig->setValue(det_tube_w * fwhm2sig);

	// angular acceptances in deg
	t_real_reso sample_w = spinTof2SampleWidth->value();
	t_real_reso thetacrit = spinTof2GuideCoating->value()*std::asin(lam/10.*std::sqrt(Nb/pi));
	spinTof2thISig->setValue(tl::r2d(thetacrit));
	spinTofphISig->setValue(tl::r2d(thetacrit));
	spinTof2thFSig->setValue(tl::r2d(2.*det_tube_w/sample_w/sample_w
		* (dist_sample_det - std::sqrt(dist_sample_det*dist_sample_det - sample_w*sample_w))));
	spinTofphFSig->setValue(tl::r2d(std::atan2(spinTof2DetHeight->value()/100., dist_sample_det)));

	radioTofDetCyl->setChecked(true);
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
