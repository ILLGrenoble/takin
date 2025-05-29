/**
 * Goto Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 15-oct-2014
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

#include "GotoDlg.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"

#include <QFileDialog>
#include <QMessageBox>


using t_real = t_real_glob;
static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();


GotoDlg::GotoDlg(QWidget* pParent, QSettings* pSett) : QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);
	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 2);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") &&
			font.fromString(m_pSettings->value("main/font_gen", "").toString()))
		{
			setFont(font);
		}
	}

	btnAdd->setIcon(load_icon("res/icons/list-add.svg"));
	btnDel->setIcon(load_icon("res/icons/list-remove.svg"));
	btnSave->setIcon(load_icon("res/icons/document-save.svg"));
	btnLoad->setIcon(load_icon("res/icons/document-open.svg"));

	std::vector<QLineEdit*> vecObjs {editH, editK, editL};
	std::vector<QLineEdit*> vecAngles {edit2ThetaM, editThetaM,
		edit2ThetaA, editThetaA, edit2ThetaS, editThetaS};

	QObject::connect(buttonBox, &QDialogButtonBox::clicked, this, &GotoDlg::ButtonBoxClicked);

	QObject::connect(btnAdd, &QAbstractButton::clicked,
		this, static_cast<void(GotoDlg::*)()>(&GotoDlg::AddPosToList));
	QObject::connect(btnDel, &QAbstractButton::clicked,
		this, &GotoDlg::RemPosFromList);
	QObject::connect(btnLoad, &QAbstractButton::clicked, this, &GotoDlg::LoadList);
	QObject::connect(btnSave, &QAbstractButton::clicked, this, &GotoDlg::SaveList);
	QObject::connect(listSeq, &QListWidget::itemSelectionChanged,
		this, &GotoDlg::ListItemSelected);
	QObject::connect(listSeq, &QListWidget::itemDoubleClicked,
		this, &GotoDlg::ListItemDoubleClicked);

	QObject::connect(editKi, &QLineEdit::textEdited, this, &GotoDlg::EditedKiKf);
	QObject::connect(editKf, &QLineEdit::textEdited, this, &GotoDlg::EditedKiKf);
	QObject::connect(editE, &QLineEdit::textEdited, this, &GotoDlg::EditedE);
	QObject::connect(btnGetPos, &QAbstractButton::clicked, this, &GotoDlg::GetCurPos);

	for(QLineEdit* pObj : vecObjs)
		QObject::connect(pObj, &QLineEdit::textEdited, this, &GotoDlg::CalcSample);
	for(QLineEdit* pObj : vecAngles)
		QObject::connect(pObj, &QLineEdit::textEdited, this, &GotoDlg::EditedAngles);

	if(m_pSettings)
	{
		if(m_pSettings->contains("goto_pos/geo"))
			restoreGeometry(m_pSettings->value("goto_pos/geo").toByteArray());
		if(m_pSettings->contains("goto_pos/ki_fix"))
		{
			bool bKi = m_pSettings->value("goto_pos/ki_fix").toBool();
			if(bKi)
				radioFixedKi->setChecked(true);
			else
				radioFixedKf->setChecked(true);
		}
	}
}


GotoDlg::~GotoDlg()
{
	ClearList();
}


void GotoDlg::CalcSample()
{
	m_bSampleOk = false;

	std::vector<QLineEdit*> vecEdits{editThetaS, edit2ThetaS};
	for(QLineEdit* pEdit : vecEdits)
		pEdit->setText("");

	t_real dH = tl::str_to_var_parse<t_real>(editH->text().toStdString());
	t_real dK = tl::str_to_var_parse<t_real>(editK->text().toStdString());
	t_real dL = tl::str_to_var_parse<t_real>(editL->text().toStdString());
	t_real dKi = tl::str_to_var_parse<t_real>(editKi->text().toStdString());
	t_real dKf = tl::str_to_var_parse<t_real>(editKf->text().toStdString());

	ublas::vector<t_real> vecQ;
	bool bFailed = false;
	try
	{
		tl::get_tas_angles(m_lattice,
			m_vec1, m_vec2,
			dKi, dKf,
			dH, dK, dL,
			m_bSenseS,
			&m_dSampleTheta, &m_dSample2Theta,
			&vecQ);

		if(tl::is_nan_or_inf<t_real>(m_dSample2Theta))
			throw tl::Err("Invalid sample 2theta angle.");
		if(tl::is_nan_or_inf<t_real>(m_dSampleTheta))
			throw tl::Err("Invalid sample theta angle.");
	}
	catch(const std::exception& ex)
	{
		editThetaS->setText("invalid");
		edit2ThetaS->setText("invalid");

		labelStatus->setText((std::string("Error (sample): ") + ex.what()).c_str());
		bFailed = true;
	}

	tl::set_eps_0(vecQ, g_dEps);
	tl::set_eps_0(m_dSample2Theta, g_dEps);
	tl::set_eps_0(m_dSampleTheta, g_dEps);

	const std::wstring strAA = tl::get_spec_char_utf16("AA") +
		tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");

	std::wostringstream ostrStatus;
	ostrStatus << "Q = " << vecQ << " " << strAA << ", |Q| = "
		<< ublas::norm_2(vecQ) << " " << strAA;
	labelQ->setText(QString::fromWCharArray(ostrStatus.str().c_str()));

#ifndef NDEBUG
	try
	{
		ublas::vector<t_real> vecQ0 = vecQ;
		vecQ0.resize(2, true);
		t_real dAngleQVec0 = tl::r2d(tl::vec_angle(vecQ0));
		tl::log_info("Angle Q Orient0: ", dAngleQVec0, " deg.");
	}
	catch(const std::exception& ex)
	{}
#endif

	if(bFailed) return;

	editThetaS->setText(tl::var_to_str(tl::r2d(m_dSampleTheta), g_iPrec).c_str());
	edit2ThetaS->setText(tl::var_to_str(tl::r2d(m_dSample2Theta), g_iPrec).c_str());

	m_bSampleOk = true;

	if(m_bMonoAnaOk && m_bSampleOk)
		labelStatus->setText("Position OK.");
}


void GotoDlg::CalcMonoAna()
{
	m_bMonoAnaOk = false;

	std::vector<QLineEdit*> vecEdits{editThetaM, edit2ThetaM, editThetaA, edit2ThetaA};
	for(QLineEdit* pEdit : vecEdits)
		pEdit->setText("");

	t_real dKi = tl::str_to_var_parse<t_real>(editKi->text().toStdString());
	t_real dKf = tl::str_to_var_parse<t_real>(editKf->text().toStdString());

	tl::t_wavenumber_si<t_real> ki = dKi / angs;
	tl::t_wavenumber_si<t_real> kf = dKf / angs;

	bool bMonoOk = false;
	bool bAnaOk = false;

	try
	{
		m_dMono2Theta = tl::get_mono_twotheta(ki, m_dMono*angs, m_bSenseM) / rads;
		if(tl::is_nan_or_inf<t_real>(m_dMono2Theta))
			throw tl::Err("Invalid monochromator angle.");

		tl::set_eps_0(m_dMono2Theta, g_dEps);
		bMonoOk = true;
	}
	catch(const std::exception& ex)
	{
		editThetaM->setText("invalid");
		edit2ThetaM->setText("invalid");

		//tl::log_err("mono: ", ex.what());
		labelStatus->setText((std::string("Error (mono)): ") + ex.what()).c_str());
	}

	try
	{
		m_dAna2Theta = tl::get_mono_twotheta(kf, m_dAna*angs, m_bSenseA) / rads;
		if(tl::is_nan_or_inf<t_real>(m_dAna2Theta))
			throw tl::Err("Invalid analysator angle.");

		tl::set_eps_0(m_dAna2Theta, g_dEps);
		bAnaOk = true;
	}
	catch(const std::exception& ex)
	{
		editThetaA->setText("invalid");
		edit2ThetaA->setText("invalid");

		//tl::log_err("ana: ", ex.what());
		labelStatus->setText((std::string("Error (ana): ") + ex.what()).c_str());
	}


	if(bMonoOk)
	{
		t_real dTMono = m_dMono2Theta / 2.;

		edit2ThetaM->setText(tl::var_to_str(tl::r2d(m_dMono2Theta), g_iPrec).c_str());
		editThetaM->setText(tl::var_to_str(tl::r2d(dTMono), g_iPrec).c_str());
	}
	if(bAnaOk)
	{
		t_real dTAna = m_dAna2Theta / 2.;

		edit2ThetaA->setText(tl::var_to_str(tl::r2d(m_dAna2Theta), g_iPrec).c_str());
		editThetaA->setText(tl::var_to_str(tl::r2d(dTAna), g_iPrec).c_str());
	}

	if(bMonoOk && bAnaOk)
		m_bMonoAnaOk = true;
	if(m_bMonoAnaOk && m_bSampleOk)
		labelStatus->setText("Position OK.");
}


void GotoDlg::EditedKiKf()
{
	t_real dKi = tl::str_to_var_parse<t_real>(editKi->text().toStdString());
	t_real dKf = tl::str_to_var_parse<t_real>(editKf->text().toStdString());

	tl::t_energy_si<t_real> Ei = tl::k2E(dKi / angs);
	tl::t_energy_si<t_real> Ef = tl::k2E(dKf / angs);

	tl::t_energy_si<t_real> E = Ei - Ef;
	t_real dE = E/meV;
	tl::set_eps_0(dE, g_dEps);

	std::string strE = tl::var_to_str<t_real>(dE, g_iPrec);
	editE->setText(strE.c_str());

	CalcMonoAna();
	CalcSample();
}


void GotoDlg::EditedE()
{
	t_real dE = tl::str_to_var_parse<t_real>(editE->text().toStdString());
	tl::t_energy_si<t_real> E = dE * meV;

	bool bImag = false;
	tl::t_wavenumber_si<t_real> k_E = tl::E2k(E, bImag);
	t_real dSign = 1.;
	if(bImag) dSign = -1.;

	if(radioFixedKi->isChecked())
	{
		t_real dKi = tl::str_to_var_parse<t_real>(editKi->text().toStdString());

		tl::t_wavenumber_si<t_real> ki = dKi / angs;
		tl::t_wavenumber_si<t_real> kf =
			tl::my_units_sqrt<tl::t_wavenumber_si<t_real>>(ki*ki - dSign*k_E*k_E);

		t_real dKf = kf*angs;
		tl::set_eps_0(dKf, g_dEps);

		std::string strKf = tl::var_to_str<t_real>(dKf, g_iPrec);
		editKf->setText(strKf.c_str());
	}
	else
	{
		t_real dKf = tl::str_to_var_parse<t_real>(editKf->text().toStdString());

		tl::t_wavenumber_si<t_real> kf = dKf / angs;
		tl::t_wavenumber_si<t_real> ki =
			tl::my_units_sqrt<tl::t_wavenumber_si<t_real>>(kf*kf + dSign*k_E*k_E);

		t_real dKi = ki*angs;
		tl::set_eps_0(dKi, g_dEps);

		std::string strKi = tl::var_to_str<t_real>(dKi, g_iPrec);
		editKi->setText(strKi.c_str());
	}

	CalcMonoAna();
	CalcSample();
}


/**
 * calc. tas angles -> hkl
 */
void GotoDlg::EditedAngles()
{
	t_real th_m = tl::d2r(tl::str_to_var_parse<t_real>(edit2ThetaM->text().toStdString())/2.);
	tl::set_eps_0(th_m, g_dEps);
	editThetaM->setText(tl::var_to_str<t_real>(tl::r2d(th_m), g_iPrec).c_str());

	t_real th_a = tl::d2r(tl::str_to_var_parse<t_real>(edit2ThetaA->text().toStdString())/2.);
	tl::set_eps_0(th_a, g_dEps);
	editThetaA->setText(tl::var_to_str<t_real>(tl::r2d(th_a), g_iPrec).c_str());

	t_real th_s = tl::d2r(tl::str_to_var_parse<t_real>(editThetaS->text().toStdString()));
	t_real tt_s = tl::d2r(tl::str_to_var_parse<t_real>(edit2ThetaS->text().toStdString()));

	t_real h,k,l;
	t_real dKi, dKf, dE;
	ublas::vector<t_real> vecQ;
	bool bFailed = false;
	try
	{
		tl::get_hkl_from_tas_angles<t_real>(m_lattice,
			m_vec1, m_vec2,
			m_dMono, m_dAna,
			th_m, th_a, th_s, tt_s,
			m_bSenseM, m_bSenseA, m_bSenseS,
			&h, &k, &l,
			&dKi, &dKf, &dE, 0,
			&vecQ);

		if(tl::is_nan_or_inf<t_real>(h) || tl::is_nan_or_inf<t_real>(k) || tl::is_nan_or_inf<t_real>(l))
			throw tl::Err("Invalid hkl.");
	}
	catch(const std::exception& ex)
	{
		labelStatus->setText((std::string("Error (hkl): ") + ex.what()).c_str());
		bFailed = true;
	}

	const std::wstring strAA = tl::get_spec_char_utf16("AA") +
		tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");
	std::wostringstream ostrStatus;
	ostrStatus << "Q = " << vecQ << " " << strAA
		<< ", |Q| = " << ublas::norm_2(vecQ) << " " << strAA;
	labelQ->setText(QString::fromWCharArray(ostrStatus.str().c_str()));

	if(bFailed) return;

	for(t_real* d : {&h,&k,&l, &dKi,&dKf,&dE})
		tl::set_eps_0(*d, g_dEps);

	editH->setText(tl::var_to_str<t_real>(h, g_iPrec).c_str());
	editK->setText(tl::var_to_str<t_real>(k, g_iPrec).c_str());
	editL->setText(tl::var_to_str<t_real>(l, g_iPrec).c_str());

	editKi->setText(tl::var_to_str<t_real>(dKi, g_iPrec).c_str());
	editKf->setText(tl::var_to_str<t_real>(dKf, g_iPrec).c_str());
	editE->setText(tl::var_to_str<t_real>(dE, g_iPrec).c_str());

	m_dMono2Theta = th_m*2.;
	m_dAna2Theta = th_a*2.;
	m_dSample2Theta = tt_s;
	m_dSampleTheta = th_s;
	m_bMonoAnaOk = true;
	m_bSampleOk = true;

	if(m_bMonoAnaOk && m_bSampleOk)
		labelStatus->setText("Position OK.");
}


void GotoDlg::GetCurPos()
{
	if(!m_bHasParamsRecip)
		return;

	tl::set_eps_0(m_paramsRecip.dki, g_dEps);
	tl::set_eps_0(m_paramsRecip.dkf, g_dEps);
	tl::set_eps_0(m_paramsRecip.Q_rlu[0], g_dEps);
	tl::set_eps_0(m_paramsRecip.Q_rlu[1], g_dEps);
	tl::set_eps_0(m_paramsRecip.Q_rlu[2], g_dEps);

	editKi->setText(tl::var_to_str(m_paramsRecip.dki, g_iPrec).c_str());
	editKf->setText(tl::var_to_str(m_paramsRecip.dkf, g_iPrec).c_str());

	editH->setText(tl::var_to_str(-m_paramsRecip.Q_rlu[0], g_iPrec).c_str());
	editK->setText(tl::var_to_str(-m_paramsRecip.Q_rlu[1], g_iPrec).c_str());
	editL->setText(tl::var_to_str(-m_paramsRecip.Q_rlu[2], g_iPrec).c_str());

	EditedKiKf();
	CalcMonoAna();
	CalcSample();
}


void GotoDlg::RecipParamsChanged(const RecipParams& params)
{
	m_bHasParamsRecip = true;
	m_paramsRecip = params;
}


void GotoDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ApplyRole ||
	   buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		ApplyCurPos();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}

	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
		{
			m_pSettings->setValue("goto_pos/geo", saveGeometry());
			m_pSettings->setValue("goto_pos/ki_fix", radioFixedKi->isChecked());
		}

		QDialog::accept();
	}
}


void GotoDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


//------------------------------------------------------------------------------


/**
 * instrument position stored in the list
 */
struct HklPos
{
	t_real dh, dk, dl;
	t_real dki, dkf;
	t_real dE;
};


bool GotoDlg::ApplyCurPos()
{
	if(!m_bMonoAnaOk || !m_bSampleOk)
		return false;

	CrystalOptions crys;
	TriangleOptions triag;

	triag.bChangedMonoTwoTheta = true;
	triag.dMonoTwoTheta = m_dMono2Theta;

	triag.bChangedAnaTwoTheta = true;
	triag.dAnaTwoTheta = m_dAna2Theta;

	triag.bChangedTheta = true;
	triag.dTheta = m_dSampleTheta;

	triag.bChangedAngleKiVec0 = true;
	triag.dAngleKiVec0 = tl::get_pi<t_real>()/2. - m_dSampleTheta;

	triag.bChangedTwoTheta = true;
	triag.dTwoTheta = m_dSample2Theta;

	emit vars_changed(crys, triag);
	return true;
}


bool GotoDlg::GotoPos(QListWidgetItem* pItem, bool bApply)
{
	if(!pItem) return false;
	HklPos* pPos = (HklPos*)pItem->data(Qt::UserRole).value<void*>();
	if(!pPos) return false;

	editH->setText(tl::var_to_str(pPos->dh, g_iPrec).c_str());
	editK->setText(tl::var_to_str(pPos->dk, g_iPrec).c_str());
	editL->setText(tl::var_to_str(pPos->dl, g_iPrec).c_str());
	editKi->setText(tl::var_to_str(pPos->dki, g_iPrec).c_str());
	editKf->setText(tl::var_to_str(pPos->dkf, g_iPrec).c_str());

	EditedKiKf();
	CalcSample();

	if(bApply)
		return ApplyCurPos();
	return m_bMonoAnaOk && m_bSampleOk;
}


bool GotoDlg::GotoPos(unsigned int iItem)
{
	if(int(iItem) >= listSeq->count())
		return false;
	return GotoPos(listSeq->item(int(iItem)), 1);
}


void GotoDlg::ListItemSelected()
{
	QListWidgetItem *pItem = listSeq->currentItem();
	GotoPos(pItem, 0);
}


void GotoDlg::ListItemDoubleClicked(QListWidgetItem* pItem)
{
	bool bOk = GotoPos(pItem, 1);
	if(!bOk)
		QMessageBox::critical(this, "Error", "Invalid position.");
}


void GotoDlg::AddPosToList(t_real dh, t_real dk, t_real dl, t_real dki, t_real dkf)
{
	HklPos *pPos = new HklPos;

	pPos->dh = dh;
	pPos->dk = dk;
	pPos->dl = dl;
	pPos->dki = dki;
	pPos->dkf = dkf;
	pPos->dE = (tl::k2E(pPos->dki/angs) - tl::k2E(pPos->dkf/angs))/tl::meV;

	tl::set_eps_0(pPos->dh, g_dEps);
	tl::set_eps_0(pPos->dk, g_dEps);
	tl::set_eps_0(pPos->dl, g_dEps);
	tl::set_eps_0(pPos->dki, g_dEps);
	tl::set_eps_0(pPos->dkf, g_dEps);
	tl::set_eps_0(pPos->dE, g_dEps);

	const std::wstring strAA = tl::get_spec_char_utf16("AA") +
		tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");

	std::wostringstream ostrHKL;
	ostrHKL.precision(g_iPrecGfx);
	ostrHKL << "(" << pPos->dh << ", " << pPos->dk << ", " << pPos->dl << ") rlu\n";
	ostrHKL << "ki = " << pPos->dki << " " << strAA;
	ostrHKL << ", kf = " << pPos->dkf << " " << strAA << "\n";
	ostrHKL << "E = " << pPos->dE << " meV";

	QString qstr = QString::fromWCharArray(ostrHKL.str().c_str());
	QListWidgetItem* pItem = new QListWidgetItem(qstr, listSeq);
	pItem->setData(Qt::UserRole, QVariant::fromValue<void*>(pPos));
}


void GotoDlg::AddPosToList()
{
	t_real dh = tl::str_to_var_parse<t_real>(editH->text().toStdString());
	t_real dk = tl::str_to_var_parse<t_real>(editK->text().toStdString());
	t_real dl = tl::str_to_var_parse<t_real>(editL->text().toStdString());
	t_real dki = tl::str_to_var_parse<t_real>(editKi->text().toStdString());
	t_real dkf = tl::str_to_var_parse<t_real>(editKf->text().toStdString());

	AddPosToList(dh, dk, dl, dki, dkf);
}


void GotoDlg::RemPosFromList()
{
	QListWidgetItem *pItem = listSeq->currentItem();
	if(pItem)
	{
		HklPos* pPos = (HklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(pPos)
			delete pPos;
		delete pItem;
	}
}


void GotoDlg::ClearList()
{
	while(listSeq->count())
	{
		QListWidgetItem *pItem = listSeq->item(0);
		if(!pItem)
			break;

		HklPos* pPos = (HklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(pPos)
			delete pPos;
		delete pItem;
	}
}


void GotoDlg::LoadList()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = "~";
	if(m_pSettings)
		strDirLast = m_pSettings->value("goto_pos/last_dir", "~").toString();
	QString qstrFile = QFileDialog::getOpenFileName(this,
		"Load Positions", strDirLast,
		"TAZ files (*.taz *.TAZ)", nullptr,
		fileopt);
	if(qstrFile == "")
		return;


	std::string strFile = qstrFile.toStdString();
	std::string strDir = tl::get_dir(strFile);

	tl::Prop<std::string> xml;
	if(!xml.Load(strFile.c_str(), tl::PropType::XML))
	{
		QMessageBox::critical(this, "Error", "Could not load positions.");
		return;
	}

	Load(xml, strXmlRoot);
	if(m_pSettings)
		m_pSettings->setValue("goto_pos/last_dir", QString(strDir.c_str()));
}


void GotoDlg::SaveList()
{
	const std::string strXmlRoot("taz/");

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_pSettings && !m_pSettings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = "~";
	if(m_pSettings)
		m_pSettings->value("goto_pos/last_dir", "~").toString();
	QString qstrFile = QFileDialog::getSaveFileName(this,
		"Save Positions", strDirLast,
		"TAZ files (*.taz *.TAZ)", nullptr,
		fileopt);

	if(qstrFile == "")
		return;

	std::string strFile = qstrFile.toStdString();
	std::string strDir = tl::get_dir(strFile);
	if(tl::get_fileext(strFile,1) != "taz")
		strFile += ".taz";

	std::map<std::string, std::string> mapConf;
	Save(mapConf, strXmlRoot);

	tl::Prop<std::string> xml;
	xml.Add(mapConf);
	bool bOk = xml.Save(strFile.c_str(), tl::PropType::XML);
	if(!bOk)
		QMessageBox::critical(this, "Error", "Could not save positions.");

	if(bOk && m_pSettings)
		m_pSettings->setValue("goto_pos/last_dir", QString(strDir.c_str()));
}


/**
 * save configuration
 */
void GotoDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	mapConf[strXmlRoot + "goto_pos/h"] = editH->text().toStdString();
	mapConf[strXmlRoot + "goto_pos/k"] = editK->text().toStdString();
	mapConf[strXmlRoot + "goto_pos/l"] = editL->text().toStdString();
	mapConf[strXmlRoot + "goto_pos/ki"] = editKi->text().toStdString();
	mapConf[strXmlRoot + "goto_pos/kf"] = editKf->text().toStdString();
	mapConf[strXmlRoot + "goto_pos/cki"] = radioFixedKi->isChecked()?"1":"0";

	// favlist
	for(int iItem = 0; iItem < listSeq->count(); ++iItem)
	{
		const QListWidgetItem *pItem = listSeq->item(iItem);
		if(!pItem)
			continue;
		const HklPos* pPos = (HklPos*)pItem->data(Qt::UserRole).value<void*>();
		if(!pPos)
			continue;

		std::ostringstream ostrItemBase;
		ostrItemBase << "goto_favlist/pos_" << iItem << "/";
		std::string strItemBase = ostrItemBase.str();

		mapConf[strXmlRoot + strItemBase + "h"] = tl::var_to_str(pPos->dh);
		mapConf[strXmlRoot + strItemBase + "k"] = tl::var_to_str(pPos->dk);
		mapConf[strXmlRoot + strItemBase + "l"] = tl::var_to_str(pPos->dl);
		mapConf[strXmlRoot + strItemBase + "ki"] = tl::var_to_str(pPos->dki);
		mapConf[strXmlRoot + strItemBase + "kf"] = tl::var_to_str(pPos->dkf);
	}
}


/**
 * load configuration
 */
void GotoDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	bool bOk = false;

	editH->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "goto_pos/h", 1., &bOk), g_iPrec).c_str());
	editK->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "goto_pos/k", 0., &bOk), g_iPrec).c_str());
	editL->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "goto_pos/l", 0., &bOk), g_iPrec).c_str());
	editKi->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "goto_pos/ki", 1.4, &bOk), g_iPrec).c_str());
	editKf->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "goto_pos/kf", 1.4, &bOk), g_iPrec).c_str());
	radioFixedKi->setChecked(xml.Query<bool>(strXmlRoot + "goto_pos/cki", false, &bOk));

	// favlist
	ClearList();
	std::size_t iItem = 0;
	while(true)
	{
		std::ostringstream ostrItemBase;
		ostrItemBase << "goto_favlist/pos_" << iItem << "/";
		std::string strItemBase = ostrItemBase.str();

		// last position reached?
		if(!xml.Exists(strXmlRoot + strItemBase))
			break;

		t_real dh = xml.Query<t_real>(strXmlRoot + strItemBase + "h", 0., &bOk);
		t_real dk = xml.Query<t_real>(strXmlRoot + strItemBase + "k", 0., &bOk);
		t_real dl = xml.Query<t_real>(strXmlRoot + strItemBase + "l", 0., &bOk);
		t_real dki = xml.Query<t_real>(strXmlRoot + strItemBase + "ki", 0., &bOk);
		t_real dkf = xml.Query<t_real>(strXmlRoot + strItemBase + "kf", 0., &bOk);

		AddPosToList(dh, dk, dl, dki, dkf);
		++iItem;
	}

	EditedKiKf();
	CalcSample();
}


#include "moc_GotoDlg.cpp"
