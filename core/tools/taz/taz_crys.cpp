/**
 * TAS tool (crystal stuff)
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2014 - 2016
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

#include "taz.h"
#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/helper/exception.h"
#include <boost/algorithm/string.hpp>


using t_real = t_real_glob;


static inline QString dtoqstr(t_real dVal, unsigned int iPrec=8)
{
	std::string str = tl::var_to_str(dVal, iPrec);
	return QString(str.c_str());
}


std::ostream& operator<<(std::ostream& ostr, const tl::Lattice<t_real>& lat)
{
	ostr << "a = " << lat.GetA();
	ostr << ", b = " << lat.GetB();
	ostr << ", c = " << lat.GetC();
	ostr << ", alpha = " << lat.GetAlpha();
	ostr << ", beta = " << lat.GetBeta();
	ostr << ", gamma = " << lat.GetGamma();
	return ostr;
}


void TazDlg::emitSampleParams()
{
	t_real a = editA->text().toDouble();
	t_real b = editB->text().toDouble();
	t_real c = editC->text().toDouble();

	t_real alpha = tl::d2r(editAlpha->text().toDouble());
	t_real beta = tl::d2r(editBeta->text().toDouble());
	t_real gamma = tl::d2r(editGamma->text().toDouble());

	t_real dX0 = editScatX0->text().toDouble();
	t_real dX1 = editScatX1->text().toDouble();
	t_real dX2 = editScatX2->text().toDouble();

	t_real dY0 = editScatY0->text().toDouble();
	t_real dY1 = editScatY1->text().toDouble();
	t_real dY2 = editScatY2->text().toDouble();

	SampleParams sampleparams;
	sampleparams.dAngles[0] = alpha; sampleparams.dAngles[1] = beta; sampleparams.dAngles[2] = gamma;
	sampleparams.dLattice[0] = a; sampleparams.dLattice[1] = b; sampleparams.dLattice[2] = c;
	sampleparams.dPlane1[0] = dX0; sampleparams.dPlane1[1] = dX1; sampleparams.dPlane1[2] = dX2;
	sampleparams.dPlane2[0] = dY0; sampleparams.dPlane2[1] = dY1; sampleparams.dPlane2[2] = dY2;
	emit SampleParamsChanged(sampleparams);
}


void TazDlg::SetCrystalType()
{
	m_crystalsys = xtl::CrystalSystem::CRYS_NOT_SET;

	xtl::SpaceGroup<t_real> *pSpaceGroup = 0;
	int iSpaceGroupIdx = comboSpaceGroups->currentIndex();
	if(iSpaceGroupIdx != 0)
		pSpaceGroup = (xtl::SpaceGroup<t_real>*)comboSpaceGroups->itemData(iSpaceGroupIdx).value<void*>();
	if(pSpaceGroup)
		m_crystalsys = pSpaceGroup->GetCrystalSystem();

	CheckCrystalType();
}


void TazDlg::CheckCrystalType()
{
	if(m_settings.value("main/ignore_xtal_restrictions", 0).toBool())
	{
		for(QLineEdit* pEdit : {editA, editB, editC, editAlpha, editBeta, editGamma,
			editARecip, editBRecip, editCRecip,
			editAlphaRecip, editBetaRecip, editGammaRecip})
			pEdit->setEnabled(true);
		return;
	}

	set_crystal_system_edits(m_crystalsys, editA, editB, editC,
		editAlpha, editBeta, editGamma,
		editARecip, editBRecip, editCRecip,
		editAlphaRecip, editBetaRecip, editGammaRecip);
}


void TazDlg::CalcPeaksRecip()
{
	if(!m_bReady) return;

	try
	{
		t_real a = editARecip->text().toDouble();
		t_real b = editBRecip->text().toDouble();
		t_real c = editCRecip->text().toDouble();

		t_real alpha = tl::d2r(editAlphaRecip->text().toDouble());
		t_real beta = tl::d2r(editBetaRecip->text().toDouble());
		t_real gamma = tl::d2r(editGammaRecip->text().toDouble());

		tl::Lattice<t_real> lattice(a,b,c, alpha,beta,gamma);
		tl::Lattice<t_real> recip = lattice.GetRecip();

		editA->setText(dtoqstr(recip.GetA(), g_iPrec));
		editB->setText(dtoqstr(recip.GetB(), g_iPrec));
		editC->setText(dtoqstr(recip.GetC(), g_iPrec));
		editAlpha->setText(dtoqstr(tl::r2d(recip.GetAlpha()), g_iPrec));
		editBeta->setText(dtoqstr(tl::r2d(recip.GetBeta()), g_iPrec));
		editGamma->setText(dtoqstr(tl::r2d(recip.GetGamma()), g_iPrec));

		m_bUpdateRecipEdits = 0;
		CalcPeaks();
		m_bUpdateRecipEdits = 1;
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


void TazDlg::CalcPeaks()
{
	if(!m_bReady || !m_sceneRecip.GetTriangle() || !m_sceneRealLattice.GetLattice())
		return;

	try
	{
		const bool bPowder = checkPowder->isChecked();

		// lattice
		const t_real a = editA->text().toDouble();
		const t_real b = editB->text().toDouble();
		const t_real c = editC->text().toDouble();

		const t_real alpha = tl::d2r(editAlpha->text().toDouble());
		const t_real beta = tl::d2r(editBeta->text().toDouble());
		const t_real gamma = tl::d2r(editGamma->text().toDouble());

#ifndef NDEBUG
		tl::log_debug("lattice consts = ", a, ", ", b, ", ", c);
		tl::log_debug("lattice angles = ", alpha, ", ", beta, ", ", gamma);
#endif

		if(a<=0. || b<=0. || c<=0. || alpha<=0. || beta<=0. || gamma<=0.)
			throw tl::Err("Invalid lattice definition.");

		tl::Lattice<t_real> lattice(a,b,c, alpha,beta,gamma);
		tl::Lattice<t_real> recip = lattice.GetRecip();


		//----------------------------------------------------------------------
		// scattering plane
		t_real dX0 = editScatX0->text().toDouble();
		t_real dX1 = editScatX1->text().toDouble();
		t_real dX2 = editScatX2->text().toDouble();
		ublas::vector<t_real> vecPlaneXRLU = tl::make_vec({dX0, dX1, dX2});

		t_real dY0 = editScatY0->text().toDouble();
		t_real dY1 = editScatY1->text().toDouble();
		t_real dY2 = editScatY2->text().toDouble();
		ublas::vector<t_real> vecPlaneYRLU = tl::make_vec({dY0, dY1, dY2});

		// display up vector
		ublas::matrix<t_real> matGCov = recip.GetMetricCov();
		ublas::vector<t_real> vecPlaneZRLU =
			tl::cross_prod_contra(matGCov, vecPlaneXRLU, vecPlaneYRLU, false);
		vecPlaneZRLU /= tl::veclen(vecPlaneZRLU);
		tl::set_eps_0(vecPlaneZRLU, g_dEpsGfx);

		editScatZ0->setText(tl::var_to_str(vecPlaneZRLU[0], g_iPrecGfx).c_str());
		editScatZ1->setText(tl::var_to_str(vecPlaneZRLU[1], g_iPrecGfx).c_str());
		editScatZ2->setText(tl::var_to_str(vecPlaneZRLU[2], g_iPrecGfx).c_str());

		ublas::vector<t_real> vecX0 = ublas::zero_vector<t_real>(3);
		tl::Plane<t_real> planeRLU(vecX0, vecPlaneXRLU, vecPlaneYRLU);
		if(!planeRLU.IsValid())
		{
			tl::log_err("Invalid scattering plane.");
			return;
		}
		//----------------------------------------------------------------------


		//----------------------------------------------------------------------
		// view plane for real lattice
		t_real dX0R = editRealX0->text().toDouble();
		t_real dX1R = editRealX1->text().toDouble();
		t_real dX2R = editRealX2->text().toDouble();
		ublas::vector<t_real> vecPlaneXR = tl::make_vec({dX0R, dX1R, dX2R});

		t_real dY0R = editRealY0->text().toDouble();
		t_real dY1R = editRealY1->text().toDouble();
		t_real dY2R = editRealY2->text().toDouble();
		ublas::vector<t_real> vecPlaneYR = tl::make_vec({dY0R, dY1R, dY2R});

		// display up vector
		ublas::matrix<t_real> matGCovR = lattice.GetMetricCov();
		ublas::vector<t_real> vecPlaneZR =
			tl::cross_prod_contra(matGCovR, vecPlaneXR, vecPlaneYR, false);
		vecPlaneZR /= tl::veclen(vecPlaneZR);
		tl::set_eps_0(vecPlaneZR, g_dEpsGfx);

		editRealZ0->setText(tl::var_to_str(vecPlaneZR[0], g_iPrecGfx).c_str());
		editRealZ1->setText(tl::var_to_str(vecPlaneZR[1], g_iPrecGfx).c_str());
		editRealZ2->setText(tl::var_to_str(vecPlaneZR[2], g_iPrecGfx).c_str());

		ublas::vector<t_real> vecX0R = ublas::zero_vector<t_real>(3);
		tl::Plane<t_real> planeRealFrac(vecX0R, vecPlaneXR, vecPlaneYR);
		if(!planeRealFrac.IsValid())
		{
			tl::log_err("Invalid view plane for real lattice.");
			return;
		}
		//----------------------------------------------------------------------


		if(m_pGotoDlg)
		{
			m_pGotoDlg->SetLattice(lattice);
			m_pGotoDlg->SetScatteringPlane(
				tl::make_vec({dX0, dX1, dX2}),
				tl::make_vec({dY0, dY1, dY2}));
			m_pGotoDlg->CalcSample();
		}

		if(m_pElasticDlg)
		{
			m_pElasticDlg->SetLattice(lattice);
			m_pElasticDlg->SetScatteringPlane(
				tl::make_vec({dX0, dX1, dX2}),
				tl::make_vec({dY0, dY1, dY2}));
			m_pElasticDlg->CalcSpuriousPositions();
		}

		emitSampleParams();


		if(m_bUpdateRecipEdits)
		{
			editARecip->setText(dtoqstr(recip.GetA(), g_iPrec));
			editBRecip->setText(dtoqstr(recip.GetB(), g_iPrec));
			editCRecip->setText(dtoqstr(recip.GetC(), g_iPrec));
			editAlphaRecip->setText(dtoqstr(tl::r2d(recip.GetAlpha()), g_iPrec));
			editBetaRecip->setText(dtoqstr(tl::r2d(recip.GetBeta()), g_iPrec));
			editGammaRecip->setText(dtoqstr(tl::r2d(recip.GetGamma()), g_iPrec));
		}

		const std::wstring& strAA = tl::get_spec_char_utf16("AA");
		const std::wstring& strMinus = tl::get_spec_char_utf16("sup-");
		const std::wstring& strThree = tl::get_spec_char_utf16("sup3");

		t_real dVol = lattice.GetVol();
		t_real dVol_recip = recip.GetVol() /*/ (2.*M_PI*2.*M_PI*2.*M_PI)*/;

		std::wostringstream ostrSample;
		ostrSample.precision(g_iPrecGfx);
		ostrSample << "Unit Cell Volume: ";
		ostrSample << "Real: " << dVol << " " << strAA << strThree;
		ostrSample << ", Recip.: " << dVol_recip << " " << strAA << strMinus << strThree;
		m_pStatusMsg->setText(QString::fromWCharArray(ostrSample.str().c_str()));


		std::string strCryTy = "<not set>";
		xtl::SpaceGroup<t_real>* pSpaceGroup = nullptr;
		int iSpaceGroupIdx = comboSpaceGroups->currentIndex();
		if(iSpaceGroupIdx != 0)
			pSpaceGroup = (xtl::SpaceGroup<t_real>*)
				comboSpaceGroups->itemData(iSpaceGroupIdx).value<void*>();

		if(pSpaceGroup)
			strCryTy = pSpaceGroup->GetCrystalSystemName();

		editCrystalSystem->setText(strCryTy.c_str());


		m_latticecommon = xtl::LatticeCommon<t_real_glob>();
		if(m_latticecommon.Calc(lattice, recip, planeRLU, planeRealFrac, pSpaceGroup, &m_vecAtoms))
		{
			m_sceneRecip.GetTriangle()->CalcPeaks(m_latticecommon, bPowder);
			if(m_sceneRecip.getSnapq())
				m_sceneRecip.GetTriangle()->SnapToNearestPeak(
					m_sceneRecip.GetTriangle()->GetNodeGq());
			m_sceneRecip.emitUpdate();

			m_sceneProjRecip.GetLattice()->CalcPeaks(m_latticecommon, true);
			m_sceneRealLattice.GetLattice()->CalcPeaks(m_latticecommon);

#ifndef NO_3D
			if(m_pRecip3d)
				m_pRecip3d->CalcPeaks(m_latticecommon);
			if(m_pReal3d)
				m_pReal3d->CalcPeaks(m_sceneRealLattice.GetLattice()->GetWS3D(),
					m_latticecommon);
			if(m_pBZ3d)
				m_pBZ3d->RenderBZ(m_sceneRecip.GetTriangle()->GetBZ3D(),
					m_latticecommon,
					&m_sceneRecip.GetTriangle()->GetBZ3DPlaneVerts(),
					&m_sceneRecip.GetTriangle()->GetBZ3DSymmVerts());
#endif
		}
		else
		{
			tl::log_err("Lattice calculations failed.");
		}

		m_dlgRealParam.CrystalChanged(m_latticecommon);
	}
	catch(const std::exception& ex)
	{
		m_sceneRecip.GetTriangle()->ClearPeaks();
		tl::log_err(ex.what());
	}
}


void TazDlg::VarsChanged(const CrystalOptions& crys, const TriangleOptions& triag)
{
	// update crystal
	bool bRecalcPeaks = false;

	if(crys.bChangedSampleName)
		editDescr->setText(tl::trimmed(crys.strSampleName).c_str());

	if(crys.bChangedLattice)
	{
		QString qstr0 = QString::number(crys.dLattice[0]);
		QString qstr1 = QString::number(crys.dLattice[1]);
		QString qstr2 = QString::number(crys.dLattice[2]);

		if(qstr0!=editA->text() || qstr1!=editB->text() || qstr2!=editC->text())
		{
			editA->setText(qstr0);
			editB->setText(qstr1);
			editC->setText(qstr2);

			bRecalcPeaks = true;
		}
	}

	if(crys.bChangedLatticeAngles)
	{
		QString qstr0 = QString::number(crys.dLatticeAngles[0]);
		QString qstr1 = QString::number(crys.dLatticeAngles[1]);
		QString qstr2 = QString::number(crys.dLatticeAngles[2]);

		if(qstr0!=editAlpha->text() || qstr1!=editBeta->text() || qstr2!=editGamma->text())
		{
			editAlpha->setText(qstr0);
			editBeta->setText(qstr1);
			editGamma->setText(qstr2);

			bRecalcPeaks = true;
		}
	}

	if(crys.bChangedSpacegroup)
	{
		editSpaceGroupsFilter->clear();
		int iSGIdx = comboSpaceGroups->findText(crys.strSpacegroup.c_str(), Qt::MatchContains /*| Qt::MatchCaseSensitive*/);
		if(iSGIdx >= 0)
			comboSpaceGroups->setCurrentIndex(iSGIdx);
		else
			comboSpaceGroups->setCurrentIndex(0);

		bRecalcPeaks = true;
	}

	if(crys.bChangedPlane1)
	{
		QString qstr0 = QString::number(crys.dPlane1[0]);
		QString qstr1 = QString::number(crys.dPlane1[1]);
		QString qstr2 = QString::number(crys.dPlane1[2]);

		if(qstr0!=editScatX0->text() || qstr1!=editScatX1->text() || qstr2!=editScatX2->text())
		{
			editScatX0->setText(qstr0);
			editScatX1->setText(qstr1);
			editScatX2->setText(qstr2);

			bRecalcPeaks = true;
		}
	}
	if(crys.bChangedPlane2)
	{
		QString qstr0 = QString::number(crys.dPlane2[0]);
		QString qstr1 = QString::number(crys.dPlane2[1]);
		QString qstr2 = QString::number(crys.dPlane2[2]);

		if(qstr0!=editScatY0->text() || qstr1!=editScatY1->text() || qstr2!=editScatY2->text())
		{
			editScatY0->setText(qstr0);
			editScatY1->setText(qstr1);
			editScatY2->setText(qstr2);

			bRecalcPeaks = true;
		}
	}

	if(bRecalcPeaks)
	{
		try
		{
			CalcPeaks();
		}
		catch(const std::exception& ex)
		{
			tl::log_err("Cannot recalculate peaks: ", ex.what(), ".");
		}
	}


	// ------------------------------------------------------------------------


	// update scattering triangle
	if(triag.bChangedMonoD)
	{
		QString qstr = QString::number(triag.dMonoD);
		if(qstr != editMonoD->text())
		{
			editMonoD->setText(qstr);
			UpdateDs();
		}
	}

	if(triag.bChangedAnaD)
	{
		QString qstr = QString::number(triag.dAnaD);
		if(qstr != editAnaD->text())
		{
			editAnaD->setText(qstr);
			UpdateDs();
		}
	}


	// negative sign already handled in scattering sense
	if(triag.bChangedTwoTheta && !checkSenseS->isChecked())
		const_cast<TriangleOptions&>(triag).dTwoTheta = -triag.dTwoTheta;

	if(triag.IsAnythingChanged() /*|| crys.IsAnythingChanged()*/)
	{
		m_sceneReal.triangleChanged(triag);
		m_sceneReal.emitUpdate(triag);

		UpdateMonoSense();
		UpdateAnaSense();
		UpdateSampleSense();
	}

	if(triag.bChangedAngleKiVec0)
	{
		m_sceneRecip.tasChanged(triag);
		m_sceneRecip.emitUpdate();
	}
}


void TazDlg::RotatePlane(unsigned iAxis, t_real dAngle)
{
	try
	{
		m_bReady = false;

		t_real dX0 = editScatX0->text().toDouble();
		t_real dX1 = editScatX1->text().toDouble();
		t_real dX2 = editScatX2->text().toDouble();
		ublas::vector<t_real> vecX = tl::make_vec({dX0, dX1, dX2});

		t_real dY0 = editScatY0->text().toDouble();
		t_real dY1 = editScatY1->text().toDouble();
		t_real dY2 = editScatY2->text().toDouble();
		ublas::vector<t_real> vecY = tl::make_vec({dY0, dY1, dY2});

		ublas::vector<t_real> vecZ = tl::cross_3(vecX, vecY);

		std::vector<ublas::vector<t_real>> vecOrth =
			tl::gram_schmidt<ublas::vector<t_real>>
				({vecX, vecY, vecZ}, 1);

		ublas::matrix<t_real> matRot =
			tl::rotation_matrix(vecOrth[iAxis], dAngle);
		vecX = ublas::prod(matRot, vecOrth[0]);
		vecY = ublas::prod(matRot, vecOrth[1]);

		tl::set_eps_0(vecX, g_dEps);
		tl::set_eps_0(vecY, g_dEps);

		editScatX0->setText(tl::var_to_str(vecX[0], g_iPrec).c_str());
		editScatX1->setText(tl::var_to_str(vecX[1], g_iPrec).c_str());
		editScatX2->setText(tl::var_to_str(vecX[2], g_iPrec).c_str());
		editScatY0->setText(tl::var_to_str(vecY[0], g_iPrec).c_str());
		editScatY1->setText(tl::var_to_str(vecY[1], g_iPrec).c_str());
		editScatY2->setText(tl::var_to_str(vecY[2], g_iPrec).c_str());

		m_bReady = true;
		CalcPeaks();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


void TazDlg::RepopulateSpaceGroups()
{
	std::shared_ptr<const xtl::SpaceGroups<t_real>> sgs = xtl::SpaceGroups<t_real>::GetInstance();
	const xtl::SpaceGroups<t_real>::t_mapSpaceGroups* pmapSpaceGroups = sgs->get_space_groups();
	if(!pmapSpaceGroups)
		return;

	QString strCurSG = comboSpaceGroups->currentText();
	comboSpaceGroups->setCurrentIndex(0);

	for(int iCnt=comboSpaceGroups->count()-1; iCnt>0; --iCnt)
		comboSpaceGroups->removeItem(iCnt);

	// filter string
	std::string strFilter = editSpaceGroupsFilter->text().toStdString();
	// look for space group number instead of name?
	bool bWantsNr = tl::str_is_digits<decltype(strFilter)>(strFilter);

	// number, string containing number + name, pointer
	using t_sgitem = std::tuple<unsigned int, std::string, void*>;
	std::vector<t_sgitem> vecSGs;

	// get all space groups and apply filter
	for(const xtl::SpaceGroups<t_real>::t_mapSpaceGroups::value_type& pair : *pmapSpaceGroups)
	{
		// space group number and name
		const std::string strSGNr = tl::var_to_str(pair.second.GetNr());
		const std::string& strName = pair.second.GetName();

		// apply user filter
		typedef const boost::iterator_range<std::string::const_iterator> t_striterrange;
		if(strFilter!="")
		{
			// sg number or name?
			const std::string *pStrToFilter = bWantsNr ? &strSGNr : &strName;

			if(!boost::ifind_first(t_striterrange(pStrToFilter->begin(), pStrToFilter->end()),
				t_striterrange(strFilter.begin(), strFilter.end())))
				continue;
		}

		vecSGs.push_back(std::make_tuple(pair.second.GetNr(),
			"(" + strSGNr + ") " + strName,
			(void*)&pair.second));
	}

	// sort by space group number
	std::sort(vecSGs.begin(), vecSGs.end(),
		[](const t_sgitem& sg1, const t_sgitem& sg2)
		{
			return std::get<0>(sg1) < std::get<0>(sg2);
		});

	for(std::size_t iSG=0; iSG<vecSGs.size(); ++iSG)
	{
		// insert space groups
		comboSpaceGroups->insertItem(comboSpaceGroups->count(),
			std::get<1>(vecSGs[iSG]).c_str(), QVariant::fromValue(std::get<2>(vecSGs[iSG])));
	}


	int iSGIdx = comboSpaceGroups->findText(strCurSG);
	if(iSGIdx > 0)
	{	// select previously selected space group
		comboSpaceGroups->setCurrentIndex(iSGIdx);
	}
	else if(vecSGs.size() == 1)
	{	// if old sg is not compatible with filter -> select new sg, if unique
		int iNewSGIdx = comboSpaceGroups->findText(std::get<1>(vecSGs[0]).c_str());
		if(iNewSGIdx > 0)
			comboSpaceGroups->setCurrentIndex(iNewSGIdx);
	}
	else
	{	// not set
		comboSpaceGroups->setCurrentIndex(0);
	}
}


//--------------------------------------------------------------------------------
// spurion stuff

void TazDlg::ShowSpurionDlg()
{
	if(!m_pSpuri)
	{
		m_pSpuri = new SpurionDlg(this, &m_settings);
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged,
			m_pSpuri, &SpurionDlg::paramsChanged);

		m_sceneRecip.emitAllParams();
	}

	focus_dlg(m_pSpuri);
}


void TazDlg::spurionInfo(const tl::ElasticSpurion& spuri,
	const std::vector<tl::InelasticSpurion<t_real>>& vecInelCKI,
	const std::vector<tl::InelasticSpurion<t_real>>& vecInelCKF)
{
	std::ostringstream ostrMsg;
	if(spuri.bAType || spuri.bMType || vecInelCKI.size() || vecInelCKF.size())
		ostrMsg << "Warning: ";

	if(spuri.bAType || spuri.bMType)
	{
		ostrMsg << "Spurious elastic scattering of type ";
		if(spuri.bAType && spuri.bMType)
		{
			ostrMsg << "A and M";
			ostrMsg << (spuri.bAKfSmallerKi ? " (kf<ki)" : " (kf>ki)");
		}
		else if(spuri.bAType)
		{
			ostrMsg << "A";
			ostrMsg << (spuri.bAKfSmallerKi ? " (kf<ki)" : " (kf>ki)");
		}
		else if(spuri.bMType)
		{
			ostrMsg << "M";
			ostrMsg << (spuri.bMKfSmallerKi ? " (kf<ki)" : " (kf>ki)");
		}
		ostrMsg << ".";

		if(vecInelCKI.size() || vecInelCKF.size())
			ostrMsg << " ";
	}

	const std::string& strDelta = tl::get_spec_char_utf8("Delta");

	if(vecInelCKI.size())
	{
		ostrMsg << "Spurious inelastic CKI scattering at "
			<< strDelta << "E = ";
		for(unsigned int iInel=0; iInel<vecInelCKI.size(); ++iInel)
		{
			ostrMsg << vecInelCKI[iInel].dE_meV << " meV";
			if(iInel != vecInelCKI.size()-1)
				ostrMsg << ", ";
		}
		ostrMsg << ".";

		if(vecInelCKF.size())
			ostrMsg << " ";
	}

	if(vecInelCKF.size())
	{
		ostrMsg << "Spurious inelastic CKF scattering at "
			<< strDelta << "E = ";
		for(unsigned int iInel=0; iInel<vecInelCKF.size(); ++iInel)
		{
			ostrMsg << vecInelCKF[iInel].dE_meV << " meV";
			if(iInel != vecInelCKF.size()-1)
				ostrMsg << ", ";
		}
		ostrMsg << ".";
	}

	m_pStatusMsg->setText(QString::fromUtf8(ostrMsg.str().c_str(), ostrMsg.str().size()));
}



//--------------------------------------------------------------------------------
// reso stuff
void TazDlg::InitReso()
{
	if(!m_pReso)
	{
		m_pReso = new ResoDlg(this, &m_settings);

		QObject::connect(this, &TazDlg::ResoParamsChanged, m_pReso, &ResoDlg::ResoParamsChanged);
		QObject::connect(&m_sceneRecip, &ScatteringTriangleScene::paramsChanged, m_pReso, &ResoDlg::RecipParamsChanged);
		QObject::connect(&m_sceneReal, &TasLayoutScene::paramsChanged, m_pReso, &ResoDlg::RealParamsChanged);
		QObject::connect(this, &TazDlg::SampleParamsChanged, m_pReso, &ResoDlg::SampleParamsChanged);

		UpdateDs();
		UpdateMonoSense();
		UpdateSampleSense();
		UpdateAnaSense();

		emitSampleParams();
		m_sceneRecip.emitAllParams();
		m_sceneReal.emitAllParams();
	}
}


void TazDlg::ShowResoParams()
{
	InitReso();
	focus_dlg(m_pReso);
}


void TazDlg::ShowResoEllipses()
{
	InitReso();

	if(!m_pEllipseDlg)
	{
		m_pEllipseDlg = new EllipseDlg(this, &m_settings);
		QObject::connect(m_pReso, &ResoDlg::ResoResultsSig, m_pEllipseDlg, &EllipseDlg::SetParams);

		m_pReso->EmitResults();
	}

	focus_dlg(m_pEllipseDlg);
}


void TazDlg::InitResoConv()
{
	if(!m_pConvoDlg)
		m_pConvoDlg = new ConvoDlg(this, &m_settings);
}


void TazDlg::ShowResoConv()
{
	InitResoConv();
	focus_dlg(m_pConvoDlg);
}


#ifndef NO_3D
void TazDlg::ShowResoEllipses3D()
{
	InitReso();

	if(!m_pEllipseDlg3D)
	{
		m_pEllipseDlg3D = new EllipseDlg3D(this, &m_settings);
		QObject::connect(m_pReso, &ResoDlg::ResoResultsSig, m_pEllipseDlg3D, &EllipseDlg3D::SetParams);

		m_pReso->EmitResults();
	}

	focus_dlg(m_pEllipseDlg3D);
}

#else
void TazDlg::ShowResoEllipses3D() {}
#endif



//--------------------------------------------------------------------------------
// spacegroups & formfactors dialogs

void TazDlg::ShowFormfactorDlg()
{
	if(!m_pFormfactorDlg)
		m_pFormfactorDlg = new FormfactorDlg(this, &m_settings);

	focus_dlg(m_pFormfactorDlg);
}


void TazDlg::ShowAtomsDlg(bool bOnlyCreate)
{
	if(!m_pAtomsDlg)
	{
		m_pAtomsDlg = new AtomsDlg(this, &m_settings);
		QObject::connect(m_pAtomsDlg, &AtomsDlg::ApplyAtoms, this, &TazDlg::ApplyAtoms);
	}

	m_pAtomsDlg->SetAtoms(m_vecAtoms);

	if(!bOnlyCreate)
		focus_dlg(m_pAtomsDlg);
}


void TazDlg::ApplyAtoms(const std::vector<xtl::AtomPos<t_real>>& vecAtoms)
{
	try
	{
		m_vecAtoms = vecAtoms;
		CalcPeaks();
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


void TazDlg::ShowSgListDlg()
{
	if(!m_pSgListDlg)
		m_pSgListDlg = new SgListDlg(this);

	focus_dlg(m_pSgListDlg);
}
