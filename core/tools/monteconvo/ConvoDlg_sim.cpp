/**
 * monte carlo convolution tool -> convolution simulations
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2015, 2016
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

#include "ConvoDlg.h"

#include "tlibs/time/chrono.h"
#include "tlibs/time/stopwatch.h"
#include "tlibs/helper/thread.h"
#include "tlibs/math/stat.h"


using t_real = t_real_reso;
using t_stopwatch = tl::Stopwatch<t_real>;


/**
 * create 1d convolution
 */
void ConvoDlg::Start1D()
{
	StartSim1D(false, tl::get_rand_seed());
}


/**
 * create 1d convolution
 */
void ConvoDlg::StartSim1D(bool bForceDeferred, unsigned int seed)
{
	m_atStop.store(false);
	ClearPlot1D();

	bool bUseScan = m_bUseScan && checkScan->isChecked();
	t_real dScale = tl::str_to_var<t_real>(editScale->text().toStdString());
	t_real dSlope = tl::str_to_var<t_real>(editSlope->text().toStdString());
	t_real dOffs = tl::str_to_var<t_real>(editOffs->text().toStdString());

	bool bRecycleNeutrons = checkRnd->isChecked();
	bool bFlipCoords = checkFlip->isChecked();
	bool bLiveResults = m_pLiveResults->isChecked();
	bool bLivePlots = m_pLivePlots->isChecked();
	std::string strAutosave = editAutosave->text().toStdString();

	btnStart->setEnabled(false);
	btnStartFit->setEnabled(false);
	tabSettings->setEnabled(false);
	tabOptions->setEnabled(false);
	m_pMenuBar->setEnabled(false);
	if(m_pSqwParamDlg)
		m_pSqwParamDlg->setEnabled(false);
	editScale->setEnabled(false);
	editSlope->setEnabled(false);
	editOffs->setEnabled(false);
	btnStop->setEnabled(true);
	tabWidget->setCurrentWidget(tabPlot);

	Qt::ConnectionType connty = bForceDeferred
		? Qt::ConnectionType::DirectConnection
		: Qt::ConnectionType::BlockingQueuedConnection;

	std::function<void()> fkt = [this, connty, bForceDeferred, bUseScan, bFlipCoords,
		seed, bRecycleNeutrons, dScale, dSlope, dOffs, bLiveResults, bLivePlots, strAutosave]
	{
		std::function<void()> fktEnableButtons = [this]
		{
			QMetaObject::invokeMethod(btnStop, "setEnabled", Q_ARG(bool, false));
			QMetaObject::invokeMethod(tabSettings, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(tabOptions, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(m_pMenuBar, "setEnabled", Q_ARG(bool, true));
			if(m_pSqwParamDlg)
				QMetaObject::invokeMethod(m_pSqwParamDlg, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editScale, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editSlope, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editOffs, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(btnStart, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(btnStartFit, "setEnabled", Q_ARG(bool, true));
		};

		t_stopwatch watch;
		watch.start();

		const unsigned int iNumNeutrons = spinNeutrons->value();
		const unsigned int iNumSampleSteps = spinSampleSteps->value();
		const unsigned int iNumSteps = spinStepCnt->value();

		bool bScanAxisFound = 0;
		int iScanAxisIdx = 0;
		std::string strScanVar = "";
		std::vector<std::vector<t_real>> vecAxes;
		std::tie(bScanAxisFound, iScanAxisIdx, strScanVar, vecAxes) = get_scan_axis<t_real>(
			true, comboAxis->currentIndex(), spinStepCnt->value(), m_eps_rlu,
			t_real(spinStartH->value()), t_real(spinStopH->value()), t_real(spinStartK->value()), t_real(spinStopK->value()),
			t_real(spinStartL->value()), t_real(spinStopL->value()), t_real(spinStartE->value()), t_real(spinStopE->value()));
		if(!bScanAxisFound)
		{
			//QMessageBox::critical(this, "Error", "No scan variable found.");
			tl::log_err("No scan variable found.");
			fktEnableButtons();
			return;
		}

		const std::vector<t_real> *pVecScanX = &vecAxes[iScanAxisIdx];
		const std::vector<t_real>& vecH = vecAxes[0];
		const std::vector<t_real>& vecK = vecAxes[1];
		const std::vector<t_real>& vecL = vecAxes[2];
		const std::vector<t_real>& vecE = vecAxes[3];


		QMetaObject::invokeMethod(m_plotwrap.get(), "setAxisTitle",
			Q_ARG(int, QwtPlot::yLeft),
			Q_ARG(const QString&, QString("S(Q,E) (a.u.)")));
		QMetaObject::invokeMethod(m_plotwrap.get(), "setAxisTitle",
			Q_ARG(int, QwtPlot::xBottom),
			Q_ARG(const QString&, QString(strScanVar.c_str())));


		// -------------------------------------------------------------------------
		// Load reso file
		TASReso reso;
		reso.SetPlaneDistTolerance(m_eps_plane);

		std::string _strResoFile = editRes->text().toStdString();
		tl::trim(_strResoFile);
		const std::string strResoFile = find_file_in_global_paths(_strResoFile);

		tl::log_debug("Loading resolution from \"", strResoFile, "\".");
		if(strResoFile == "" || !reso.LoadRes(strResoFile.c_str()))
		{
			//QMessageBox::critical(this, "Error", "Could not load resolution file.");
			fktEnableButtons();
			return;
		}
		// -------------------------------------------------------------------------


		if(bUseScan)	// get crystal definition from scan file
		{
			ublas::vector<t_real> vec1 =
				tl::make_vec({m_scan.plane.vec1[0], m_scan.plane.vec1[1], m_scan.plane.vec1[2]});
			ublas::vector<t_real> vec2 =
				tl::make_vec({m_scan.plane.vec2[0], m_scan.plane.vec2[1], m_scan.plane.vec2[2]});

			reso.SetLattice(m_scan.sample.a, m_scan.sample.b, m_scan.sample.c,
				m_scan.sample.alpha, m_scan.sample.beta, m_scan.sample.gamma,
				vec1, vec2);
		}
		else	// use crystal config file
		{
			// -------------------------------------------------------------------------
			// Load lattice
			std::string _strLatticeFile = editCrys->text().toStdString();
			tl::trim(_strLatticeFile);
			const std::string strLatticeFile = find_file_in_global_paths(_strLatticeFile);

			tl::log_debug("Loading crystal from \"", strLatticeFile, "\".");
			if(strLatticeFile == "" || !reso.LoadLattice(strLatticeFile.c_str(), bFlipCoords))
			{
				//QMessageBox::critical(this, "Error", "Could not load crystal file.");
				fktEnableButtons();
				return;
			}
			// -------------------------------------------------------------------------
		}

		reso.SetAlgo(ResoAlgo(comboAlgo->currentIndex()+1));
		reso.SetKiFix(comboFixedK->currentIndex()==0);
		reso.SetKFix(spinKfix->value());
		reso.SetOptimalFocus(get_reso_focus(comboFocMono->currentIndex(), comboFocAna->currentIndex()));


		if(m_pSqw == nullptr || !m_pSqw->IsOk())
		{
			//QMessageBox::critical(this, "Error", "No valid S(Q,E) model loaded.");
			fktEnableButtons();
			return;
		}



		// meta data
		std::ostringstream ostrOut;
		ostrOut.precision(g_iPrec);

		ostrOut << "#\n";
		write_takin_metadata(ostrOut);
		ostrOut << "# MC neutrons: " << iNumNeutrons << "\n";
		ostrOut << "# MC sample steps: " << iNumSampleSteps << "\n";
		ostrOut << "# Scale: " << dScale << "\n";
		ostrOut << "# Slope: " << dSlope << "\n";
		ostrOut << "# Offset: " << dOffs << "\n";
		if(m_strLastFile != "")
			ostrOut << "# File: " << m_strLastFile << "\n";
		if(editScan->text() != "")
			ostrOut << "# Scan file: " << editScan->text().toStdString() << "\n";
		ostrOut << "#\n";

		ostrOut << std::left << std::setw(g_iPrec*2) << "# h" << " "
			<< std::left << std::setw(g_iPrec*2) << "k" << " "
			<< std::left << std::setw(g_iPrec*2) << "l" << " "
			<< std::left << std::setw(g_iPrec*2) << "E" << " "
			<< std::left << std::setw(g_iPrec*2) << "S(Q, E)" << " "
			<< std::left << std::setw(g_iPrec*2) << "S_scaled(Q, E)"
			<< "\n";

		QMetaObject::invokeMethod(editStartTime, "setText",
			Q_ARG(const QString&, QString(watch.GetStartTimeStr().c_str())));

		QMetaObject::invokeMethod(progress, "setMaximum", Q_ARG(int, iNumSteps));
		QMetaObject::invokeMethod(progress, "setValue", Q_ARG(int, 0));

		QMetaObject::invokeMethod(textResult, "clear", connty);


		m_vecQ.clear();
		m_vecS.clear();
		m_vecScaledS.clear();

		m_vecQ.reserve(iNumSteps);
		m_vecS.reserve(iNumSteps);
		m_vecScaledS.reserve(iNumSteps);

		unsigned int iNumThreads = bForceDeferred ? 0 : get_max_threads();
		tl::log_debug("Calculating using ", iNumThreads, (iNumThreads == 1 ? " thread." : " threads."));

		// function to be called before each thread
		auto th_start_func = [seed, bRecycleNeutrons, bForceDeferred]
		{
			// TODO: init random seeds for non-deferred, threaded simulation
			if(bRecycleNeutrons && bForceDeferred)
				tl::init_rand_seed(seed);
			else
				tl::init_rand();
		};

		// call the start function directly in non-threaded mode
		if(bForceDeferred)
			th_start_func();

		tl::ThreadPool<std::pair<bool, t_real>(), decltype(th_start_func)> tp(iNumThreads, &th_start_func);
		auto& lstFuts = tp.GetResults();

		for(unsigned int iStep = 0; iStep < iNumSteps; ++iStep)
		{
			t_real dCurH = vecH[iStep];
			t_real dCurK = vecK[iStep];
			t_real dCurL = vecL[iStep];
			t_real dCurE = vecE[iStep];

			tp.AddTask(
			[&reso, dCurH, dCurK, dCurL, dCurE, iNumNeutrons, iNumSampleSteps, this]()
				-> std::pair<bool, t_real>
			{
				if(this->StopRequested())
					return std::pair<bool, t_real>(false, 0.);

				t_real dS = 0.;
				t_real dhklE_mean[4] = { 0., 0., 0., 0. };

				if(iNumNeutrons == 0)
				{	// if no neutrons are given, just plot the unconvoluted S(Q,E)
					dS += (*m_pSqw)(dCurH, dCurK, dCurL, dCurE);
					dS += m_pSqw->GetBackground(dCurH, dCurK, dCurL, dCurE);
				}
				else
				{	// convolution
					TASReso localreso = reso;
					localreso.SetRandomSamplePos(iNumSampleSteps);
					std::vector<ublas::vector<t_real>> vecNeutrons;

					try
					{
						if(!localreso.SetHKLE(dCurH, dCurK, dCurL, dCurE))
						{
							std::ostringstream ostrErr;
							ostrErr << "Invalid crystal position: (" <<
								dCurH << " " << dCurK << " " << dCurL << ") rlu, "
								<< dCurE << " meV.";
							throw tl::Err(ostrErr.str().c_str());
						}
					}
					catch(const std::exception& ex)
					{
						tl::log_err(ex.what());
						return std::pair<bool, t_real>(false, 0.);
					}

					Ellipsoid4d<t_real> elli =
						localreso.GenerateMC_deferred(iNumNeutrons, vecNeutrons);

					for(const ublas::vector<t_real>& vecHKLE : vecNeutrons)
					{
						if(this->StopRequested())
							return std::pair<bool, t_real>(false, 0.);

						dS += (*m_pSqw)(vecHKLE[0], vecHKLE[1], vecHKLE[2], vecHKLE[3]);

						for(int i = 0; i < 4; ++i)
							dhklE_mean[i] += vecHKLE[i];
					}

					// normalise to mc neutron count
					dS /= t_real(iNumNeutrons*iNumSampleSteps);

					// add background
					dS += m_pSqw->GetBackground(dCurH, dCurK, dCurL, dCurE);

					for(int i=0; i<4; ++i)
						dhklE_mean[i] /= t_real(iNumNeutrons*iNumSampleSteps);

					dS *= localreso.GetResoResults().dR0 * localreso.GetR0Scale();
				}
				return std::pair<bool, t_real>(true, dS);
			});
		}

		tp.Start();
		auto iterTask = tp.GetTasks().begin();
		unsigned int iStep = 0;
		for(auto &fut : lstFuts)
		{
			if(this->StopRequested())
				break;

			// deferred (in main thread), eval this task manually
			if(iNumThreads == 0)
			{
				(*iterTask)();
				++iterTask;
			}

			std::pair<bool, t_real> pairS = fut.get();
			if(!pairS.first) break;
			t_real dS = pairS.second;
			if(tl::is_nan_or_inf(dS))
			{
				dS = t_real(0);
				tl::log_warn("S(Q,E) is invalid.");
			}

			const t_real dXVal = (*pVecScanX)[iStep];
			t_real dSScale = dScale*(dS + dSlope*dXVal) + dOffs;
			if(dSScale < 0.)
				dSScale = 0.;

			ostrOut << std::left << std::setw(g_iPrec*2) << vecH[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecK[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecL[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecE[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << dS << " "
				<< std::left << std::setw(g_iPrec*2) << dSScale
				<< "\n";

			m_vecQ.push_back(dXVal);
			m_vecS.push_back(dS);
			m_vecScaledS.push_back(dSScale);

			static const std::vector<t_real> vecNull;
			bool bIsLastStep = (iStep == lstFuts.size()-1);

			if(bLivePlots || bIsLastStep)
			{
				set_qwt_data<t_real>()(*m_plotwrap, m_vecQ, m_vecScaledS, 0, false);
				set_qwt_data<t_real>()(*m_plotwrap, m_vecQ, m_vecScaledS, 1, false);
				if(bUseScan)
					set_qwt_data<t_real>()(*m_plotwrap, m_scan.vechklE[iScanAxisIdx], m_scan.vecCts, 2, false, &m_scan.vecCtsErr);
				else
					set_qwt_data<t_real>()(*m_plotwrap, vecNull, vecNull, 2, false);

				if(bIsLastStep)
				{
					t_real_qwt dLeft = std::numeric_limits<t_real_qwt>::max();
					t_real_qwt dRight = -dLeft;
					t_real_qwt dTop = dRight;
					t_real_qwt dBottom = dLeft;

					auto minmaxQ = std::minmax_element(m_vecQ.begin(), m_vecQ.end());
					auto minmaxS = std::minmax_element(m_vecScaledS.begin(), m_vecScaledS.end());

					if(minmaxQ.first != m_vecQ.end())
					{
						dLeft = *minmaxQ.first;
						dRight = *minmaxQ.second;
					}
					if(minmaxS.first != m_vecScaledS.end())
					{
						dBottom = *minmaxS.first;
						dTop = *minmaxS.second;
					}

					if(bUseScan && m_scan.vechklE[iScanAxisIdx].size() && m_scan.vecCts.size())
					{
						auto minmaxX = std::minmax_element(m_scan.vechklE[iScanAxisIdx].begin(), m_scan.vechklE[iScanAxisIdx].end());
						auto minmaxY = std::minmax_element(m_scan.vecCts.begin(), m_scan.vecCts.end());

						dLeft = std::min<t_real_qwt>(dLeft, *minmaxX.first);
						dRight = std::max<t_real_qwt>(dRight, *minmaxX.second);
						dBottom = std::min<t_real_qwt>(dBottom, *minmaxY.first);
						dTop = std::max<t_real_qwt>(dTop, *minmaxY.second);
					}

					set_zoomer_base(m_plotwrap->GetZoomer(),
						dLeft, dRight, dTop, dBottom,
						!bForceDeferred, m_plotwrap.get());
				}
				QMetaObject::invokeMethod(m_plotwrap.get(), "doUpdate", connty);
			}

			if(bLiveResults || bIsLastStep)
			{
				if(bIsLastStep)
					ostrOut << "# ------------------------- EOF -------------------------\n";

				QMetaObject::invokeMethod(textResult, "setPlainText", connty,
					Q_ARG(const QString&, QString(ostrOut.str().c_str())));

				// autosave output
				if(strAutosave != "")
				{
					std::ofstream ofstrAutosave(strAutosave);
					ofstrAutosave << ostrOut.str() << std::endl;
				}
			}

			QMetaObject::invokeMethod(progress, "setValue", Q_ARG(int, iStep+1));
			QMetaObject::invokeMethod(editStopTime, "setText",
				Q_ARG(const QString&, QString(watch.GetEstStopTimeStr(t_real(iStep+1)/t_real(iNumSteps)).c_str())));
			++iStep;
		}


		// approximate chi^2
		if(bUseScan && m_pSqw)
		{
			const std::size_t iNumScanPts = m_scan.vecPoints.size();
			std::vector<t_real> vecSFuncY;
			vecSFuncY.reserve(iNumScanPts);

			for(std::size_t iScanPt = 0; iScanPt < iNumScanPts; ++iScanPt)
			{
				const ScanPoint& pt = m_scan.vecPoints[iScanPt];
				t_real E = pt.E / tl::one_meV;
				ublas::vector<t_real> vecScanHKLE = tl::make_vec({ pt.h, pt.k, pt.l, E });


				// find point on S(Q,E) curve closest to scan point
				std::size_t iMinIdx = 0;
				t_real dMinDist = std::numeric_limits<t_real>::max();
				for(std::size_t iStep = 0; iStep < iNumSteps; ++iStep)
				{
					ublas::vector<t_real> vecCurveHKLE =
						tl::make_vec({ vecH[iStep], vecK[iStep], vecL[iStep], vecE[iStep] });

					t_real dDist = ublas::norm_2(vecCurveHKLE - vecScanHKLE);
					if(dDist < dMinDist)
					{
						dMinDist = dDist;
						iMinIdx = iStep;
					}
				}

				// add the scaled S value from the closest point
				vecSFuncY.push_back(m_vecScaledS[iMinIdx]);
			}

			m_chi2 = tl::chi2_direct<t_real>(iNumScanPts,
				vecSFuncY.data(), m_scan.vecCts.data(), m_scan.vecCtsErr.data());
			tl::log_info("chi^2 = ", m_chi2);

			if(strAutosave != "")
			{
				std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
				ofstrAutosave << "# chi^2: " << m_chi2 << std::endl;
			}
		}


		// output elapsed time
		watch.stop();
		QMetaObject::invokeMethod(editStopTime, "setText",
			Q_ARG(const QString&, QString(watch.GetStopTimeStr().c_str())));

		if(strAutosave != "")
		{
			std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
			ofstrAutosave << "# Simulation start time: " << watch.GetStartTimeStr() << "\n";
			ofstrAutosave << "# Simulation stop time: " << watch.GetStopTimeStr() << std::endl;
		}

		fktEnableButtons();
	};


	if(bForceDeferred)
	{
		fkt();
	}
	else
	{
		if(m_pth)
		{
			if(m_pth->joinable())
				m_pth->join();
			delete m_pth;
		}
		m_pth = new std::thread(std::move(fkt));
	}
}



/**
 * create 2d convolution
 */
void ConvoDlg::Start2D()
{
	m_atStop.store(false);

	bool bFlipCoords = checkFlip->isChecked();
	bool bLiveResults = m_pLiveResults->isChecked();
	bool bLivePlots = m_pLivePlots->isChecked();
	std::string strAutosave = editAutosave->text().toStdString();

	btnStart->setEnabled(false);
	btnStartFit->setEnabled(false);
	tabSettings->setEnabled(false);
	tabOptions->setEnabled(false);
	m_pMenuBar->setEnabled(false);
	if(m_pSqwParamDlg)
		m_pSqwParamDlg->setEnabled(false);
	editScale->setEnabled(false);
	editSlope->setEnabled(false);
	editOffs->setEnabled(false);
	btnStop->setEnabled(true);
	tabWidget->setCurrentWidget(tabPlot2d);

	bool bForceDeferred = false;
	Qt::ConnectionType connty = bForceDeferred
		? Qt::ConnectionType::DirectConnection
		: Qt::ConnectionType::BlockingQueuedConnection;

	std::function<void()> fkt = [this, connty, bFlipCoords, bForceDeferred,
		bLiveResults, bLivePlots, strAutosave]
	{
		std::function<void()> fktEnableButtons = [this]
		{
			QMetaObject::invokeMethod(btnStop, "setEnabled", Q_ARG(bool, false));
			QMetaObject::invokeMethod(tabSettings, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(tabOptions, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(m_pMenuBar, "setEnabled", Q_ARG(bool, true));
			if(m_pSqwParamDlg)
				QMetaObject::invokeMethod(m_pSqwParamDlg, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editScale, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editSlope, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editOffs, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(btnStart, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(btnStartFit, "setEnabled", Q_ARG(bool, true));
		};

		t_stopwatch watch;
		watch.start();

		const unsigned int iNumNeutrons = spinNeutrons->value();
		const unsigned int iNumSampleSteps = spinSampleSteps->value();

		const unsigned int iNumSteps = std::sqrt(spinStepCnt->value());
		const t_real dStartHKL[] =
		{
			t_real(spinStartH->value()), t_real(spinStartK->value()),
			t_real(spinStartL->value()), t_real(spinStartE->value())
		};
		const t_real dDeltaHKL1[] =
		{
			t_real(spinStopH->value() - spinStartH->value()) / t_real(iNumSteps),
			t_real(spinStopK->value() - spinStartK->value()) / t_real(iNumSteps),
			t_real(spinStopL->value() - spinStartL->value()) / t_real(iNumSteps),
			t_real(spinStopE->value() - spinStartE->value()) / t_real(iNumSteps)
		};
		const t_real dDeltaHKL2[] =
		{
			t_real(spinStopH2->value() - spinStartH->value()) / t_real(iNumSteps),
			t_real(spinStopK2->value() - spinStartK->value()) / t_real(iNumSteps),
			t_real(spinStopL2->value() - spinStartL->value()) / t_real(iNumSteps),
			t_real(spinStopE2->value() - spinStartE->value()) / t_real(iNumSteps)
		};


		// -------------------------------------------------------------------------
		// find axis labels and ranges
		const int iScanAxis1 = comboAxis->currentIndex();
		const int iScanAxis2 = comboAxis2->currentIndex();

		std::string strScanVar1 = "";
		t_real dStart1{}, dStop1{};
		if(iScanAxis1 == 1 || (iScanAxis1 == 0 && !tl::float_equal(spinStartH->value(), spinStopH->value(), m_eps_rlu)))
		{
			strScanVar1 = "h (rlu)";
			dStart1 = spinStartH->value();
			dStop1 = spinStopH->value();
		}
		else if(iScanAxis1 == 2 || (iScanAxis1 == 0 && !tl::float_equal(spinStartK->value(), spinStopK->value(), m_eps_rlu)))
		{
			strScanVar1 = "k (rlu)";
			dStart1 = spinStartK->value();
			dStop1 = spinStopK->value();
		}
		else if(iScanAxis1 == 3 || (iScanAxis1 == 0 && !tl::float_equal(spinStartL->value(), spinStopL->value(), m_eps_rlu)))
		{
			strScanVar1 = "l (rlu)";
			dStart1 = spinStartL->value();
			dStop1 = spinStopL->value();
		}
		else if(iScanAxis1 == 4 || (iScanAxis1 == 0 && !tl::float_equal(spinStartE->value(), spinStopE->value(), m_eps_rlu)))
		{
			strScanVar1 = "E (meV)";
			dStart1 = spinStartE->value();
			dStop1 = spinStopE->value();
		}

		std::string strScanVar2 = "";
		t_real dStart2{}, dStop2{};
		if(iScanAxis2 == 1 || (iScanAxis2 == 0 && !tl::float_equal(spinStartH->value(), spinStopH2->value(), m_eps_rlu)))
		{
			strScanVar2 = "h (rlu)";
			dStart2 = spinStartH->value();
			dStop2 = spinStopH2->value();
		}
		else if(iScanAxis2 == 2 || (iScanAxis2 == 0 && !tl::float_equal(spinStartK->value(), spinStopK2->value(), m_eps_rlu)))
		{
			strScanVar2 = "k (rlu)";
			dStart2 = spinStartK->value();
			dStop2 = spinStopK2->value();
		}
		else if(iScanAxis2 == 3 || (iScanAxis2 == 0 && !tl::float_equal(spinStartL->value(), spinStopL2->value(), m_eps_rlu)))
		{
			strScanVar2 = "l (rlu)";
			dStart2 = spinStartL->value();
			dStop2 = spinStopL2->value();
		}
		else if(iScanAxis2 == 4 || (iScanAxis2 == 0 && !tl::float_equal(spinStartE->value(), spinStopE2->value(), m_eps_rlu)))
		{
			strScanVar2 = "E (meV)";
			dStart2 = spinStartE->value();
			dStop2 = spinStopE2->value();
		}

		QMetaObject::invokeMethod(m_plotwrap2d.get(), "setAxisTitle",
			Q_ARG(int, QwtPlot::xBottom),
			Q_ARG(const QString&, QString(strScanVar1.c_str())));

		QMetaObject::invokeMethod(m_plotwrap2d.get(), "setAxisTitle",
			Q_ARG(int, QwtPlot::yLeft),
			Q_ARG(const QString&, QString(strScanVar2.c_str())));
		// -------------------------------------------------------------------------



		// -------------------------------------------------------------------------
		// Load reso file
		TASReso reso;
		reso.SetPlaneDistTolerance(m_eps_plane);

		std::string _strResoFile = editRes->text().toStdString();
		tl::trim(_strResoFile);
		const std::string strResoFile = find_file_in_global_paths(_strResoFile);

		tl::log_debug("Loading resolution from \"", strResoFile, "\".");
		if(strResoFile == "" || !reso.LoadRes(strResoFile.c_str()))
		{
			//QMessageBox::critical(this, "Error", "Could not load resolution file.");
			fktEnableButtons();
			return;
		}
		// -------------------------------------------------------------------------


		// -------------------------------------------------------------------------
		// Load lattice
		std::string _strLatticeFile = editCrys->text().toStdString();
		tl::trim(_strLatticeFile);
		const std::string strLatticeFile = find_file_in_global_paths(_strLatticeFile);

		tl::log_debug("Loading crystal from \"", strLatticeFile, "\".");
		if(strLatticeFile == "" || !reso.LoadLattice(strLatticeFile.c_str(), bFlipCoords))
		{
			//QMessageBox::critical(this, "Error", "Could not load crystal file.");
			fktEnableButtons();
			return;
		}
		// -------------------------------------------------------------------------


		reso.SetAlgo(ResoAlgo(comboAlgo->currentIndex()+1));
		reso.SetKiFix(comboFixedK->currentIndex()==0);
		reso.SetKFix(spinKfix->value());
		reso.SetOptimalFocus(get_reso_focus(comboFocMono->currentIndex(), comboFocAna->currentIndex()));

		if(m_pSqw == nullptr || !m_pSqw->IsOk())
		{
			//QMessageBox::critical(this, "Error", "No valid S(Q,E) model loaded.");
			fktEnableButtons();
			return;
		}


		std::ostringstream ostrOut;
		ostrOut.precision(g_iPrec);
		ostrOut << "#\n";
		write_takin_metadata(ostrOut);
		ostrOut << "# MC neutrons: " << iNumNeutrons << "\n";
		ostrOut << "# MC sample steps: " << iNumSampleSteps << "\n";
		//ostrOut << "# Scale: " << dScale << "\n";
		//ostrOut << "# Slope: " << dSlope << "\n";
		//ostrOut << "# Offset: " << dOffs << "\n";
		if(m_strLastFile != "")
			ostrOut << "# File: " << m_strLastFile << "\n";
		ostrOut << "#\n";
		ostrOut << std::left << std::setw(g_iPrec*2) << "# h" << " "
			<< std::left << std::setw(g_iPrec*2) << "k" << " "
			<< std::left << std::setw(g_iPrec*2) << "l" << " "
			<< std::left << std::setw(g_iPrec*2) << "E" << " "
			<< std::left << std::setw(g_iPrec*2) << "S(Q,E)" << "\n";

		QMetaObject::invokeMethod(editStartTime2d, "setText",
			Q_ARG(const QString&, QString(watch.GetStartTimeStr().c_str())));

		QMetaObject::invokeMethod(progress, "setMaximum", Q_ARG(int, iNumSteps*iNumSteps));
		QMetaObject::invokeMethod(progress, "setValue", Q_ARG(int, 0));

		QMetaObject::invokeMethod(textResult, "clear", connty);

		// raster width & height
		m_plotwrap2d->GetRaster()->Init(iNumSteps, iNumSteps);
		m_plotwrap2d->GetRaster()->SetXRange(dStart1, dStop1);
		m_plotwrap2d->GetRaster()->SetYRange(dStart2, dStop2);
		set_zoomer_base(m_plotwrap2d->GetZoomer(),
			m_plotwrap2d->GetRaster()->GetXMin(), m_plotwrap2d->GetRaster()->GetXMax(),
			m_plotwrap2d->GetRaster()->GetYMax(), m_plotwrap2d->GetRaster()->GetYMin(),
			!bForceDeferred, m_plotwrap2d.get());

		std::vector<t_real> vecH; vecH.reserve(iNumSteps*iNumSteps);
		std::vector<t_real> vecK; vecK.reserve(iNumSteps*iNumSteps);
		std::vector<t_real> vecL; vecL.reserve(iNumSteps*iNumSteps);
		std::vector<t_real> vecE; vecE.reserve(iNumSteps*iNumSteps);

		for(unsigned int iStepY=0; iStepY<iNumSteps; ++iStepY)
		{
			for(unsigned int iStepX=0; iStepX<iNumSteps; ++iStepX)
			{
				vecH.push_back(dStartHKL[0] + dDeltaHKL2[0]*t_real(iStepY) + dDeltaHKL1[0]*t_real(iStepX));
				vecK.push_back(dStartHKL[1] + dDeltaHKL2[1]*t_real(iStepY) + dDeltaHKL1[1]*t_real(iStepX));
				vecL.push_back(dStartHKL[2] + dDeltaHKL2[2]*t_real(iStepY) + dDeltaHKL1[2]*t_real(iStepX));
				vecE.push_back(dStartHKL[3] + dDeltaHKL2[3]*t_real(iStepY) + dDeltaHKL1[3]*t_real(iStepX));
			}
		}

		unsigned int iNumThreads = bForceDeferred ? 0 : get_max_threads();
		tl::log_debug("Calculating using ", iNumThreads, (iNumThreads == 1 ? " thread." : " threads."));

		void (*pThStartFunc)() = []{ tl::init_rand(); };
		tl::ThreadPool<std::pair<bool, t_real>()> tp(iNumThreads, pThStartFunc);
		auto& lstFuts = tp.GetResults();

		for(unsigned int iStep = 0; iStep < iNumSteps*iNumSteps; ++iStep)
		{
			t_real dCurH = vecH[iStep];
			t_real dCurK = vecK[iStep];
			t_real dCurL = vecL[iStep];
			t_real dCurE = vecE[iStep];

			tp.AddTask(
			[&reso, dCurH, dCurK, dCurL, dCurE, iNumNeutrons, iNumSampleSteps, this]()
				-> std::pair<bool, t_real>
			{
				if(this->StopRequested())
					return std::pair<bool, t_real>(false, 0.);

				t_real dS = 0.;
				t_real dhklE_mean[4] = {0., 0., 0., 0.};

				if(iNumNeutrons == 0)
				{	// if no neutrons are given, just plot the unconvoluted S(Q,E)
					dS += (*m_pSqw)(dCurH, dCurK, dCurL, dCurE);
					dS += m_pSqw->GetBackground(dCurH, dCurK, dCurL, dCurE);
				}
				else
				{	// convolution
					TASReso localreso = reso;
					localreso.SetRandomSamplePos(iNumSampleSteps);
					std::vector<ublas::vector<t_real>> vecNeutrons;

					try
					{
						if(!localreso.SetHKLE(dCurH, dCurK, dCurL, dCurE))
						{
							std::ostringstream ostrErr;
							ostrErr << "Invalid crystal position: ("
								<< dCurH << " " << dCurK << " " << dCurL << ") rlu, "
								<< dCurE << " meV.";
							throw tl::Err(ostrErr.str().c_str());
						}
					}
					catch(const std::exception& ex)
					{
						tl::log_err(ex.what());
						return std::pair<bool, t_real>(false, 0.);
					}

					Ellipsoid4d<t_real> elli =
						localreso.GenerateMC_deferred(iNumNeutrons, vecNeutrons);

					for(const ublas::vector<t_real>& vecHKLE : vecNeutrons)
					{
						if(this->StopRequested())
							return std::pair<bool, t_real>(false, 0.);

						dS += (*m_pSqw)(vecHKLE[0], vecHKLE[1], vecHKLE[2], vecHKLE[3]);

						for(int i = 0; i < 4; ++i)
							dhklE_mean[i] += vecHKLE[i];
					}

					dS /= t_real(iNumNeutrons*iNumSampleSteps);
					for(int i = 0; i < 4; ++i)
						dhklE_mean[i] /= t_real(iNumNeutrons*iNumSampleSteps);

					dS *= localreso.GetResoResults().dR0 * localreso.GetR0Scale();
				}
				return std::pair<bool, t_real>(true, dS);
			});
		}

		tp.Start();
		auto iterTask = tp.GetTasks().begin();
		unsigned int iStep = 0;
		for(auto &fut : lstFuts)
		{
			if(this->StopRequested())
				break;

			// deferred (in main thread), eval this task manually
			if(iNumThreads == 0)
			{
				(*iterTask)();
				++iterTask;
			}

			std::pair<bool, t_real> pairS = fut.get();
			if(!pairS.first) break;
			t_real dS = pairS.second;
			if(tl::is_nan_or_inf(dS))
			{
				dS = t_real(0);
				tl::log_warn("S(Q,E) is invalid.");
			}

			ostrOut << std::left << std::setw(g_iPrec*2) << vecH[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecK[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecL[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecE[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << dS << "\n";

			m_plotwrap2d->GetRaster()->SetPixel(iStep%iNumSteps, iStep/iNumSteps, t_real_qwt(dS));

			bool bIsLastStep = (iStep == lstFuts.size()-1);

			if(bLivePlots || bIsLastStep)
			{
				m_plotwrap2d->GetRaster()->SetZRange();

				QMetaObject::invokeMethod(m_plotwrap2d.get(), "scaleColorBar", connty);
				QMetaObject::invokeMethod(m_plotwrap2d.get(), "doUpdate", connty);
			}

			if(bLiveResults || bIsLastStep)
			{
				if(bIsLastStep)
					ostrOut << "# ------------------------- EOF -------------------------\n";
				QMetaObject::invokeMethod(textResult, "setPlainText", connty,
					Q_ARG(const QString&, QString(ostrOut.str().c_str())));

				// autosave output
				if(strAutosave != "")
				{
					std::ofstream ofstrAutosave(strAutosave);
					ofstrAutosave << ostrOut.str() << std::endl;
				}
			}

			QMetaObject::invokeMethod(progress, "setValue", Q_ARG(int, iStep+1));
			QMetaObject::invokeMethod(editStopTime2d, "setText",
				Q_ARG(const QString&, QString(
					watch.GetEstStopTimeStr(t_real(iStep + 1)/t_real(iNumSteps*iNumSteps)).c_str())));
			++iStep;
		}

		// output elapsed time
		watch.stop();
		QMetaObject::invokeMethod(editStopTime2d, "setText",
			Q_ARG(const QString&, QString(watch.GetStopTimeStr().c_str())));

		if(strAutosave != "")
		{
			std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
			ofstrAutosave << "# Simulation start time: " << watch.GetStartTimeStr() << "\n";
			ofstrAutosave << "# Simulation stop time: " << watch.GetStopTimeStr() << std::endl;
		}

		fktEnableButtons();
	};


	if(bForceDeferred)
	{
		fkt();
	}
	else
	{
		if(m_pth)
		{
			if(m_pth->joinable())
				m_pth->join();
			delete m_pth;
		}
		m_pth = new std::thread(std::move(fkt));
	}
}



/**
 * start dispersion plot
 */
void ConvoDlg::StartDisp()
{
	m_atStop.store(false);
	ClearPlot1D();

	bool bLiveResults = m_pLiveResults->isChecked();
	bool bLivePlots = m_pLivePlots->isChecked();
	std::string strAutosave = editAutosave->text().toStdString();

	btnStart->setEnabled(false);
	btnStartFit->setEnabled(false);
	tabSettings->setEnabled(false);
	tabOptions->setEnabled(false);
	m_pMenuBar->setEnabled(false);
	if(m_pSqwParamDlg)
		m_pSqwParamDlg->setEnabled(false);
	editScale->setEnabled(false);
	editSlope->setEnabled(false);
	editOffs->setEnabled(false);
	btnStop->setEnabled(true);
	tabWidget->setCurrentWidget(tabPlot);

	bool bForceDeferred = false;
	Qt::ConnectionType connty = bForceDeferred
		? Qt::ConnectionType::DirectConnection
		: Qt::ConnectionType::BlockingQueuedConnection;

	std::function<void()> fkt = [this, connty, bForceDeferred, bLiveResults, bLivePlots, strAutosave]
	{
		std::function<void()> fktEnableButtons = [this]
		{
			QMetaObject::invokeMethod(btnStop, "setEnabled", Q_ARG(bool, false));
			QMetaObject::invokeMethod(tabSettings, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(tabOptions, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(m_pMenuBar, "setEnabled", Q_ARG(bool, true));
			if(m_pSqwParamDlg)
				QMetaObject::invokeMethod(m_pSqwParamDlg, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editScale, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editSlope, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(editOffs, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(btnStart, "setEnabled", Q_ARG(bool, true));
			QMetaObject::invokeMethod(btnStartFit, "setEnabled", Q_ARG(bool, true));
		};

		t_stopwatch watch;
		watch.start();

		const unsigned int iNumSteps = spinStepCnt->value();

		bool bScanAxisFound = 0;
		int iScanAxisIdx = 0;
		std::string strScanVar = "";
		std::vector<std::vector<t_real>> vecAxes;
		std::tie(bScanAxisFound, iScanAxisIdx, strScanVar, vecAxes) = get_scan_axis<t_real>(
			false, comboAxis->currentIndex(), spinStepCnt->value(), m_eps_rlu,
			t_real(spinStartH->value()), t_real(spinStopH->value()), t_real(spinStartK->value()), t_real(spinStopK->value()),
			t_real(spinStartL->value()), t_real(spinStopL->value()), t_real(spinStartE->value()), t_real(spinStopE->value()));
		if(!bScanAxisFound)
		{
			//QMessageBox::critical(this, "Error", "No scan variable found.");
			tl::log_err("No scan variable found.");
			fktEnableButtons();
			return;
		}

		const std::vector<t_real> *pVecScanX = &vecAxes[iScanAxisIdx];
		const std::vector<t_real>& vecH = vecAxes[0];
		const std::vector<t_real>& vecK = vecAxes[1];
		const std::vector<t_real>& vecL = vecAxes[2];


		QMetaObject::invokeMethod(m_plotwrap.get(), "setAxisTitle",
			Q_ARG(int, QwtPlot::yLeft),
			Q_ARG(const QString&, QString("E(Q) (meV)")));
		QMetaObject::invokeMethod(m_plotwrap.get(), "setAxisTitle",
			Q_ARG(int, QwtPlot::xBottom),
			Q_ARG(const QString&, QString(strScanVar.c_str())));


		if(m_pSqw == nullptr || !m_pSqw->IsOk())
		{
			//QMessageBox::critical(this, "Error", "No valid S(Q,E) model loaded.");
			fktEnableButtons();
			return;
		}



		std::ostringstream ostrOut;
		ostrOut.precision(g_iPrec);
		ostrOut << "#\n";
		write_takin_metadata(ostrOut);
		ostrOut << "# Format: h k l E1 w1 E2 w2 ... En wn\n";
		ostrOut << "#\n";

		QMetaObject::invokeMethod(editStartTime, "setText",
			Q_ARG(const QString&, QString(watch.GetStartTimeStr().c_str())));

		QMetaObject::invokeMethod(progress, "setMaximum", Q_ARG(int, iNumSteps));
		QMetaObject::invokeMethod(progress, "setValue", Q_ARG(int, 0));

		QMetaObject::invokeMethod(textResult, "clear", connty);


		m_vecvecQ.clear();
		m_vecvecE.clear();
		m_vecvecW.clear();

		unsigned int iNumThreads = bForceDeferred ? 0 : get_max_threads();
		tl::log_debug("Calculating using ", iNumThreads, (iNumThreads == 1 ? " thread." : " threads."));

		tl::ThreadPool<std::tuple<bool, std::vector<t_real>, std::vector<t_real>>()>
			tp(iNumThreads);
		auto& lstFuts = tp.GetResults();

		for(unsigned int iStep = 0; iStep < iNumSteps; ++iStep)
		{
			t_real dCurH = vecH[iStep];
			t_real dCurK = vecK[iStep];
			t_real dCurL = vecL[iStep];

			tp.AddTask([dCurH, dCurK, dCurL, this]() ->
			std::tuple<bool, std::vector<t_real>, std::vector<t_real>>
			{
				if(this->StopRequested())
					return std::make_tuple(false, std::vector<t_real>(), std::vector<t_real>());

				std::vector<t_real> vecE, vecW;
				std::tie(vecE, vecW) = m_pSqw->disp(dCurH, dCurK, dCurL);
				return std::tuple<bool, std::vector<t_real>, std::vector<t_real>>
					(true, vecE, vecW);
			});
		}

		tp.Start();
		auto iterTask = tp.GetTasks().begin();
		unsigned int iStep = 0;
		for(auto &fut : lstFuts)
		{
			if(this->StopRequested()) break;

			// deferred (in main thread), eval this task manually
			if(iNumThreads == 0)
			{
				(*iterTask)();
				++iterTask;
			}

			auto tupEW = fut.get();
			if(!std::get<0>(tupEW)) break;

			ostrOut << std::left << std::setw(g_iPrec*2) << vecH[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecK[iStep] << " "
				<< std::left << std::setw(g_iPrec*2) << vecL[iStep] << " ";
			for(std::size_t iE=0; iE<std::get<1>(tupEW).size(); ++iE)
			{
				ostrOut << std::left << std::setw(g_iPrec*2) << std::get<1>(tupEW)[iE] << " ";
				ostrOut << std::left << std::setw(g_iPrec*2) << std::get<2>(tupEW)[iE] << " ";
			}
			ostrOut << "\n";


			// store dispersion branches as separate curves
			if(std::get<1>(tupEW).size() > m_vecvecE.size())
			{
				m_vecvecQ.resize(std::get<1>(tupEW).size());
				m_vecvecE.resize(std::get<1>(tupEW).size());
				m_vecvecW.resize(std::get<2>(tupEW).size());
			}
			for(std::size_t iBranch=0; iBranch<std::get<1>(tupEW).size(); ++iBranch)
			{
				m_vecvecQ[iBranch].push_back((*pVecScanX)[iStep]);
				m_vecvecE[iBranch].push_back(std::get<1>(tupEW)[iBranch]);
				m_vecvecW[iBranch].push_back(std::get<2>(tupEW)[iBranch]);
			}


			bool bIsLastStep = (iStep == lstFuts.size()-1);

			if(bLivePlots || bIsLastStep)
			{
				for(std::size_t iBranch=0; iBranch<m_vecvecE.size() && iBranch+CONVO_DISP_CURVE_START<CONVO_MAX_CURVES; ++iBranch)
				{
					set_qwt_data<t_real>()(*m_plotwrap, m_vecvecQ[iBranch], m_vecvecE[iBranch], CONVO_DISP_CURVE_START+iBranch, false);
				}

				if(bIsLastStep)
					set_zoomer_base(m_plotwrap->GetZoomer(),
					tl::container2_cast<t_real_qwt, t_real, std::vector>()(m_vecvecQ),
					tl::container2_cast<t_real_qwt, t_real, std::vector>()(m_vecvecE),
					!bForceDeferred, m_plotwrap.get());
				QMetaObject::invokeMethod(m_plotwrap.get(), "doUpdate", connty);
			}

			if(bLiveResults || bIsLastStep)
			{
				if(bIsLastStep)
					ostrOut << "# ------------------------- EOF -------------------------\n";

				QMetaObject::invokeMethod(textResult, "setPlainText", connty,
					Q_ARG(const QString&, QString(ostrOut.str().c_str())));

				// autosave output
				if(strAutosave != "")
				{
					std::ofstream ofstrAutosave(strAutosave);
					ofstrAutosave << ostrOut.str() << std::endl;
				}
			}

			QMetaObject::invokeMethod(progress, "setValue", Q_ARG(int, iStep+1));
			QMetaObject::invokeMethod(editStopTime, "setText",
				Q_ARG(const QString&, QString(
					watch.GetEstStopTimeStr(t_real(iStep + 1)/t_real(iNumSteps)).c_str())));
			++iStep;
		}

		// output elapsed time
		watch.stop();
		QMetaObject::invokeMethod(editStopTime, "setText",
			Q_ARG(const QString&, QString(watch.GetStopTimeStr().c_str())));

		if(strAutosave != "")
		{
			std::ofstream ofstrAutosave(strAutosave, std::ios_base::app);
			ofstrAutosave << "# calculation start time = " << watch.GetStartTimeStr() << "\n";
			ofstrAutosave << "# calculation stop time = " << watch.GetStopTimeStr() << std::endl;
		}

		fktEnableButtons();
	};


	if(bForceDeferred)
	{
		fkt();
	}
	else
	{
		if(m_pth)
		{
			if(m_pth->joinable())
				m_pth->join();
			delete m_pth;
		}
		m_pth = new std::thread(std::move(fkt));
	}
}
