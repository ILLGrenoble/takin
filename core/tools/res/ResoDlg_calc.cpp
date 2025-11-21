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


using t_mat = ublas::matrix<t_real_reso>;
using t_vec = ublas::vector<t_real_reso>;

static const auto angs = tl::get_one_angstrom<t_real_reso>();
static const auto rads = tl::get_one_radian<t_real_reso>();
static const auto meV = tl::get_one_meV<t_real_reso>();
static const auto cm = tl::get_one_centimeter<t_real_reso>();
static const auto meters = tl::get_one_meter<t_real_reso>();
static const auto sec = tl::get_one_second<t_real_reso>();



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


		// alternate TOF parameters
		tof.len_ch_pulse_guide = t_real_reso(spinDistTof2PulseGuide->value()) * cm;
		tof.len_ch_mono_guide = t_real_reso(spinDistTof2MonoGuide->value()) * cm;
		tof.len_guide_sample = t_real_reso(spinDistTof2GuideSample->value()) * cm;
		tof.len_sample_det2 = t_real_reso(spinDistTof2SampleDet->value()) * cm;  // TODO: unite with len_sample_det

		tof.ch_pulse_angle_win = tl::d2r(t_real_reso(spinTof2PulseWin->value())) * rads;
		tof.ch_pulse_angle_beam = tl::d2r(t_real_reso(spinTof2PulseBeam->value())) * rads;
		tof.ch_pulse_width = t_real_reso(spinTof2PulseWidth->value()) * cm;
		tof.ch_pulse_rpm = t_real_reso(spinTof2PulseRPM->value());
		tof.ch_pulse_counterrot = checkTof2PulseCounterRot->isChecked();

		tof.ch_mono_angle_win = tl::d2r(t_real_reso(spinTof2MonoWin->value())) * rads;
		tof.ch_mono_angle_beam = tl::d2r(t_real_reso(spinTof2MonoBeam->value())) * rads;
		tof.ch_mono_width = t_real_reso(spinTof2MonoWidth->value()) * cm;
		tof.ch_mono_rpm = t_real_reso(spinTof2MonoRPM->value());
		tof.ch_mono_counterrot = checkTof2MonoCounterRot->isChecked();

		tof.endguide_width = t_real_reso(spinDistTof2GuideWidth->value()) * cm;
		tof.endguide_height = t_real_reso(spinTof2GuideHeight->value()) * cm;
		tof.endguide_coating = t_real_reso(spinTof2GuideCoating->value());

		tof.det_tube_width = t_real_reso(spinTof2DetTubeWidth->value()) * cm;
		tof.det_height = t_real_reso(spinTof2DetHeight->value()) * cm;
		tof.det_z = t_real_reso(spinTof2DetZ->value()) * cm;

		tof.sample_width = t_real_reso(spinTof2SampleWidth->value()) * cm;
		tof.sample_height = t_real_reso(spinTof2SampleHeight->value()) * cm;

		tof.mc_lengths = spinTof2MCLengths->value();


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
