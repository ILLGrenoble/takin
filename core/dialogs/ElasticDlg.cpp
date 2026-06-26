/**
 * Elastic (and Inelastic) Positions Dialog
 * @author Tobias Weber <tweber@ill.fr>
 * @date 5-jan-2024
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

#include "ElasticDlg.h"

#include "dialogs/FilePreviewDlg.h"

#include "tlibs/phys/neutrons.h"
#include "tlibs/file/loadinstr.h"
#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"

#include <boost/scope_exit.hpp>

#include <QMessageBox>


// position table columns
#define POSTAB_H     0
#define POSTAB_K     1
#define POSTAB_L     2
#define POSTAB_KI    3
#define POSTAB_KF    4

#define POSTAB_COLS  5


using t_real = t_real_glob;
using t_mat = ublas::matrix<t_real>;
using t_vec = ublas::vector<t_real>;

static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();


ElasticDlg::ElasticDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") &&
			font.fromString(m_pSettings->value("main/font_gen", "").toString()))
		{
			setFont(font);
		}

		if(m_pSettings->contains("elastic_pos/geo"))
			restoreGeometry(m_pSettings->value("elastic_pos/geo").toByteArray());
	}

	btnAddPosition->setIcon(load_icon("res/icons/list-add.svg"));
	btnDelPosition->setIcon(load_icon("res/icons/list-remove.svg"));

	// connections
	QObject::connect(btnAddPosition, &QAbstractButton::clicked,
		this, static_cast<void (ElasticDlg::*)()>(&ElasticDlg::AddPosition));
	QObject::connect(btnDelPosition, &QAbstractButton::clicked,
		this, &ElasticDlg::DelPosition);
	QObject::connect(btnGenPositions, &QAbstractButton::clicked,
		this, &ElasticDlg::GeneratePositions);
	QObject::connect(btnLoadScan, &QAbstractButton::clicked,
		this, &ElasticDlg::ImportPositions);
	QObject::connect(btnSync, &QAbstractButton::toggled,
		this, &ElasticDlg::SyncToggled);
	QObject::connect(btnGotoInel, &QAbstractButton::clicked,
		[this]() { GotoPosition(ElasticDlgGoto::INEL_ORIG); });
	QObject::connect(btnGotoElastKfKi, &QAbstractButton::clicked,
		[this]() { GotoPosition(ElasticDlgGoto::ELAST_KFKI); });
	QObject::connect(btnGotoElastKiKf, &QAbstractButton::clicked,
		[this]() { GotoPosition(ElasticDlgGoto::ELAST_KIKF); });
	QObject::connect(btnGotoInelnKfmKi, &QAbstractButton::clicked,
		[this]() { GotoPosition(ElasticDlgGoto::INEL_nKFmKI); });
	QObject::connect(btnGotoInelnKimKf, &QAbstractButton::clicked,
		[this]() { GotoPosition(ElasticDlgGoto::INEL_nKImKF); });
	QObject::connect(tablePositions, &QTableWidget::itemChanged,
		this, &ElasticDlg::CalcSpuriousPositions);
	QObject::connect(spinN, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
	[this]()
	{
		SetNM();
		CalcSpuriousPositions();
	});
	QObject::connect(spinM, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
	[this]()
	{
		SetNM();
		CalcSpuriousPositions();
	});
	QObject::connect(buttonBox, &QDialogButtonBox::clicked,
		this, &ElasticDlg::ButtonBoxClicked);

	SetNM();
}


ElasticDlg::~ElasticDlg()
{
}


/**
 * set the n and m variables and labels*
 */
void ElasticDlg::SetNM()
{
	m_n = spinN->value();
	m_m = spinM->value();
	tl::set_eps_0(m_n, g_dEps);
	tl::set_eps_0(m_m, g_dEps);

	std::ostringstream ostrnKfmKi, ostrnKimKf;
	ostrnKfmKi.precision(g_iPrec);
	ostrnKimKf.precision(g_iPrec);
	ostrnKfmKi << m_n << "*kf' = m*" << m_m << "*ki";
	ostrnKimKf << m_n << "*ki' = m*" << m_m << "*kf";

	btnGotoInelnKfmKi->setText(ostrnKfmKi.str().c_str());
	btnGotoInelnKimKf->setText(ostrnKimKf.str().c_str());
}


void ElasticDlg::SyncToggled(bool sync)
{
	m_bSyncWithMainWindow = sync;

	CalcSpuriousPositions();
}


/**
 * calculate the elastic and inelastic positions corresponding
 * to the given original inelastic ones
 */
void ElasticDlg::CalcSpuriousPositions()
{
	if(!m_bAllowCalculation)
		return;

	// clear old cached positions
	for(std::vector<ElasticDlgPos>* pos : { &m_positions_inel,
		&m_positions_elast_kfki, &m_positions_elast_kikf,
		&m_positions_inel_nkimkf, &m_positions_inel_nkfmki })
	{
		pos->clear();
		pos->reserve(tablePositions->rowCount());
	}

	t_mat matB = tl::get_B(GetLattice(), 1);
	t_mat matU = tl::get_U(GetScatteringPlaneVec1(), GetScatteringPlaneVec2(), &matB);
	t_mat matUB = ublas::prod(matU, matB);

	t_mat matBinv, matUinv;
	bool B_ok = tl::inverse(matB, matBinv);
	bool U_ok = tl::inverse(matU, matUinv);
	t_mat matUBinv = ublas::prod(matBinv, matUinv);

	const std::wstring strAA = tl::get_spec_char_utf16("AA") +
		tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");

	std::wostringstream ostrResultsKfKi, ostrResultsKiKf;
	std::wostringstream ostrResultsnKimKf, ostrResultsnKfmKi;
	ostrResultsKfKi.precision(g_iPrec);
	ostrResultsKiKf.precision(g_iPrec);
	ostrResultsnKimKf.precision(g_iPrec);
	ostrResultsnKfmKi.precision(g_iPrec);

	ostrResultsKfKi << "<b>Elastic Positions Q' Corresponding to (Q, E) with kf' := ki:</b>";
	ostrResultsKfKi << "<center><table border=\"1\" cellpadding=\"0\" width=\"95%\">";
	ostrResultsKfKi << "<tr>";
	ostrResultsKfKi << "<th><b>No.</b></th>";
	if(B_ok && U_ok)
		ostrResultsKfKi << "<th><b>Q (rlu)</b></th>";
	else
		ostrResultsKfKi << "<th><b>Q (" << strAA << ")</b></th>";
	ostrResultsKfKi << "<th><b>|Q| (" << strAA << ")</b></th>";
	ostrResultsKfKi << "<th><b>E (meV)</b></th>";
	if(B_ok && U_ok)
		ostrResultsKfKi << "<th><b>Q' (rlu)</b></th>";
	else
		ostrResultsKfKi << "<th><b>Q' (" << strAA << ")</b></th>";
	ostrResultsKfKi << "<th><b>|Q'| (" << strAA << ")</b></th>";
	ostrResultsKfKi << "</tr>";

	ostrResultsKiKf << "<b>Elastic Positions Q' Corresponding to (Q, E) with ki' := kf:</b>";
	ostrResultsKiKf << "<center><table border=\"1\" cellpadding=\"0\" width=\"95%\">";
	ostrResultsKiKf << "<tr>";
	ostrResultsKiKf << "<th><b>No.</b></th>";
	if(B_ok && U_ok)
		ostrResultsKiKf << "<th><b>Q (rlu)</b></th>";
	else
		ostrResultsKiKf << "<th><b>Q (" << strAA << ")</b></th>";
	ostrResultsKiKf << "<th><b>|Q| (" << strAA << ")</b></th>";
	ostrResultsKiKf << "<th><b>E (meV)</b></th>";
	if(B_ok && U_ok)
		ostrResultsKiKf << "<th><b>Q' (rlu)</b></th>";
	else
		ostrResultsKiKf << "<th><b>Q' (" << strAA << ")</b></th>";
	ostrResultsKiKf << "<th><b>|Q'| (" << strAA << ")</b></th>";
	ostrResultsKiKf << "</tr>";

	ostrResultsnKimKf << "<b>Inelastic Positions Q' Corresponding to (Q, E) with "
		<< m_n << " * ki' := " << m_m << " * kf:</b>";
	ostrResultsnKimKf << "<center><table border=\"1\" cellpadding=\"0\" width=\"95%\">";
	ostrResultsnKimKf << "<tr>";
	ostrResultsnKimKf << "<th><b>No.</b></th>";
	if(B_ok && U_ok)
		ostrResultsnKimKf << "<th><b>Q (rlu)</b></th>";
	else
		ostrResultsnKimKf << "<th><b>Q (" << strAA << ")</b></th>";
	ostrResultsnKimKf << "<th><b>|Q| (" << strAA << ")</b></th>";
	ostrResultsnKimKf << "<th><b>E (meV)</b></th>";
	if(B_ok && U_ok)
		ostrResultsnKimKf << "<th><b>Q' (rlu)</b></th>";
	else
		ostrResultsnKimKf << "<th><b>Q' (" << strAA << ")</b></th>";
	ostrResultsnKimKf << "<th><b>|Q'| (" << strAA << ")</b></th>";
	ostrResultsnKimKf << "<th><b>E' (meV)</b></th>";
	ostrResultsnKimKf << "</tr>";

	ostrResultsnKfmKi << "<b>Inelastic Positions Q' Corresponding to (Q, E) with "
		<< m_n << " * kf' := " << m_m << " * ki:</b>";
	ostrResultsnKfmKi << "<center><table border=\"1\" cellpadding=\"0\" width=\"95%\">";
	ostrResultsnKfmKi << "<tr>";
	ostrResultsnKfmKi << "<th><b>No.</b></th>";
	if(B_ok && U_ok)
		ostrResultsnKfmKi << "<th><b>Q (rlu)</b></th>";
	else
		ostrResultsnKfmKi << "<th><b>Q (" << strAA << ")</b></th>";
	ostrResultsnKfmKi << "<th><b>|Q| (" << strAA << ")</b></th>";
	ostrResultsnKfmKi << "<th><b>E (meV)</b></th>";
	if(B_ok && U_ok)
		ostrResultsnKfmKi << "<th><b>Q' (rlu)</b></th>";
	else
		ostrResultsnKfmKi << "<th><b>Q' (" << strAA << ")</b></th>";
	ostrResultsnKfmKi << "<th><b>|Q'| (" << strAA << ")</b></th>";
	ostrResultsnKfmKi << "<th><b>E' (meV)</b></th>";
	ostrResultsnKfmKi << "</tr>";

	spinPosIdx->setMinimum(1);
	spinPosIdx->setMaximum(tablePositions->rowCount());

	// iterate positions
	for(int row = 0; row < tablePositions->rowCount(); ++row)
	{
		try
		{
			QTableWidgetItem *item_h = tablePositions->item(row, POSTAB_H);
			QTableWidgetItem *item_k = tablePositions->item(row, POSTAB_K);
			QTableWidgetItem *item_l = tablePositions->item(row, POSTAB_L);
			QTableWidgetItem *item_ki = tablePositions->item(row, POSTAB_KI);
			QTableWidgetItem *item_kf = tablePositions->item(row, POSTAB_KF);

			if(!item_h || !item_k || !item_l || !item_ki || !item_kf)
				throw tl::Err("Invalid hkl, ki or kf.");

			// get coordinates
			ElasticDlgPos pos_inel;
			pos_inel.h = tl::str_to_var_parse<t_real>(item_h->text().toStdString());
			pos_inel.k = tl::str_to_var_parse<t_real>(item_k->text().toStdString());
			pos_inel.l = tl::str_to_var_parse<t_real>(item_l->text().toStdString());
			pos_inel.ki = tl::str_to_var_parse<t_real>(item_ki->text().toStdString());
			pos_inel.kf = tl::str_to_var_parse<t_real>(item_kf->text().toStdString());

			pos_inel.mono_2theta = tl::get_mono_twotheta(pos_inel.ki / angs, GetMonoD()*angs, GetMonoSense()) / rads;
			pos_inel.ana_2theta = tl::get_mono_twotheta(pos_inel.kf / angs, GetAnaD()*angs, GetAnaSense()) / rads;

			t_vec vecQ;
			if(tl::is_nan_or_inf<t_real>(pos_inel.mono_2theta) || tl::is_nan_or_inf<t_real>(pos_inel.ana_2theta))
				throw tl::Err("Invalid monochromator or analyser angle.");

			// get angles corresponding to given inelastic position
			tl::get_tas_angles(GetLattice(),
				GetScatteringPlaneVec1(), GetScatteringPlaneVec2(),
				pos_inel.ki, pos_inel.kf, pos_inel.h, pos_inel.k, pos_inel.l,
				GetSampleSense(),
				&pos_inel.sample_theta, &pos_inel.sample_2theta,
				&vecQ);

			if(tl::is_nan_or_inf<t_real>(pos_inel.sample_2theta)
				|| tl::is_nan_or_inf<t_real>(pos_inel.sample_theta))
				throw tl::Err("Invalid sample 2theta angle.");

			t_real dE = (tl::k2E(pos_inel.ki / angs) - tl::k2E(pos_inel.kf / angs)) / meV;
			t_vec vecQrlu = tl::prod_mv(matUBinv, vecQ);

			tl::set_eps_0(dE, g_dEps);
			tl::set_eps_0(vecQ, g_dEps);
			tl::set_eps_0(vecQrlu, g_dEps);
			tl::set_eps_0(pos_inel.sample_2theta, g_dEps);
			tl::set_eps_0(pos_inel.sample_theta, g_dEps);


			try
			{
				// get corresponding elastic analyser angle for kf' = ki
				ElasticDlgPos pos_elast = pos_inel;
				pos_elast.E = 0;
				pos_elast.kf = pos_inel.ki;
				pos_elast.ana_2theta = tl::get_mono_twotheta(pos_elast.kf/angs, GetAnaD()*angs, GetAnaSense()) / rads;

				// find Q position for kf' = ki elastic position
				t_vec vecQ1;
				tl::get_hkl_from_tas_angles<t_real>(GetLattice(),
					GetScatteringPlaneVec1(), GetScatteringPlaneVec2(),
					GetMonoD(), GetAnaD(),
					pos_elast.mono_2theta*0.5, pos_elast.ana_2theta*0.5,
					pos_elast.sample_theta, pos_elast.sample_2theta,
					GetMonoSense(), GetAnaSense(), GetSampleSense(),
					&pos_elast.h, &pos_elast.k, &pos_elast.l,
					&pos_elast.ki, &pos_elast.kf, &pos_elast.E, 0,
					&vecQ1);

				if(tl::is_nan_or_inf<t_real>(pos_elast.h)
					|| tl::is_nan_or_inf<t_real>(pos_elast.k)
					|| tl::is_nan_or_inf<t_real>(pos_elast.l))
					throw tl::Err("Invalid h' k' l'.");

				t_vec vecQ1rlu = tl::prod_mv(matUBinv, vecQ1);

				tl::set_eps_0(vecQ1, g_dEps);
				tl::set_eps_0(vecQ1rlu, g_dEps);

				for(t_real* d : { &pos_elast.h, &pos_elast.k, &pos_elast.l,
					&pos_elast.ki, &pos_elast.kf, &pos_elast.E })
					tl::set_eps_0(*d, g_dEps);


				// print inelastic position
				ostrResultsKfKi << "<tr>";
				ostrResultsKfKi << "<td>" << row+1 << "</td>";
				if(B_ok && U_ok)
					ostrResultsKfKi << "<td>" << vecQrlu[0] << ", " << vecQrlu[1] << ", " << vecQrlu[2] << "</td>";
				else
					ostrResultsKfKi << "<td>" << vecQ[0] << ", " << vecQ[1] << ", " << vecQ[2] << "</td>";
				ostrResultsKfKi << "<td>" << tl::veclen(vecQ) << "</td>";
				ostrResultsKfKi << "<td>" << dE << "</td>";

				// print results for elastic kf' = ki position
				if(B_ok && U_ok)
					ostrResultsKfKi << "<td>" << vecQ1rlu[0] << ", " << vecQ1rlu[1] << ", " << vecQ1rlu[2] << "</td>";
				else
					ostrResultsKfKi << "<td>" << vecQ1[0] << ", " << vecQ1[1] << ", " << vecQ1[2] << "</td>";
				ostrResultsKfKi << "<td>" << tl::veclen(vecQ1) << "</td>";
				ostrResultsKfKi << "</tr>";

				m_positions_elast_kfki.emplace_back(std::move(pos_elast));
			}
			catch(const std::exception& ex)
			{
				ostrResultsKfKi << "<tr>";
				ostrResultsKfKi << "<td>" << row+1 << "</td>";
				ostrResultsKfKi << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				ostrResultsKfKi << "</tr>";
			}


			try
			{
				// get corresponding elastic monochromator angle for ki' = kf
				ElasticDlgPos pos_elast = pos_inel;
				pos_elast.E = 0;
				pos_elast.ki = pos_inel.kf;
				pos_elast.mono_2theta = tl::get_mono_twotheta(pos_elast.ki/angs, GetMonoD()*angs, GetMonoSense()) / rads;

				// find Q position for ki' = kf elastic position
				t_vec vecQ2;
				tl::get_hkl_from_tas_angles<t_real>(GetLattice(),
					GetScatteringPlaneVec1(), GetScatteringPlaneVec2(),
					GetMonoD(), GetAnaD(),
					pos_elast.mono_2theta*0.5, pos_elast.ana_2theta*0.5,
					pos_elast.sample_theta, pos_elast.sample_2theta,
					GetMonoSense(), GetAnaSense(), GetSampleSense(),
					&pos_elast.h, &pos_elast.k, &pos_elast.l,
					&pos_elast.ki, &pos_elast.kf, &pos_elast.E, 0,
					&vecQ2);

				if(tl::is_nan_or_inf<t_real>(pos_elast.h)
					|| tl::is_nan_or_inf<t_real>(pos_elast.k)
					|| tl::is_nan_or_inf<t_real>(pos_elast.l))
					throw tl::Err("Invalid h' k' l'.");

				t_vec vecQ2rlu = tl::prod_mv(matUBinv, vecQ2);

				tl::set_eps_0(vecQ2, g_dEps);
				tl::set_eps_0(vecQ2rlu, g_dEps);

				for(t_real* d : { &pos_elast.h, &pos_elast.k, &pos_elast.l,
					&pos_elast.ki, &pos_elast.kf, &pos_elast.E })
					tl::set_eps_0(*d, g_dEps);


				// print inelastic position
				ostrResultsKiKf << "<tr>";
				ostrResultsKiKf << "<td>" << row+1 << "</td>";
				if(B_ok && U_ok)
					ostrResultsKiKf << "<td>" << vecQrlu[0] << ", " << vecQrlu[1] << ", " << vecQrlu[2] << "</td>";
				else
					ostrResultsKiKf << "<td>" << vecQ[0] << ", " << vecQ[1] << ", " << vecQ[2] << "</td>";
				ostrResultsKiKf << "<td>" << tl::veclen(vecQ) << "</td>";
				ostrResultsKiKf << "<td>" << dE << "</td>";

				// print results for elastic ki' = kf position
				if(B_ok && U_ok)
					ostrResultsKiKf << "<td>" << vecQ2rlu[0] << ", " << vecQ2rlu[1] << ", " << vecQ2rlu[2] << "</td>";
				else
					ostrResultsKiKf << "<td>" << vecQ2[0] << ", " << vecQ2[1] << ", " << vecQ2[2] << "</td>";
				ostrResultsKiKf << "<td>" << tl::veclen(vecQ2) << "</td>";
				ostrResultsKiKf << "</tr>";

				m_positions_elast_kikf.emplace_back(std::move(pos_elast));
			}
			catch(const std::exception& ex)
			{
				ostrResultsKiKf << "<tr>";
				ostrResultsKiKf << "<td>" << row+1 << "</td>";
				ostrResultsKiKf << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				ostrResultsKiKf << "</tr>";
			}


			try
			{
				// get corresponding monochromator angle for n*ki' = m*kf
				ElasticDlgPos pos = pos_inel;
				pos.E = 0;
				pos.ki = m_m/m_n * pos_inel.kf;
				pos.mono_2theta = tl::get_mono_twotheta(pos.ki/angs, GetMonoD()*angs, GetMonoSense()) / rads;

				// find Q position for n*ki' = m*kf position
				t_vec vecQ2;
				tl::get_hkl_from_tas_angles<t_real>(GetLattice(),
					GetScatteringPlaneVec1(), GetScatteringPlaneVec2(),
					GetMonoD(), GetAnaD(),
					pos.mono_2theta*0.5, pos.ana_2theta*0.5,
					pos.sample_theta, pos.sample_2theta,
					GetMonoSense(), GetAnaSense(), GetSampleSense(),
					&pos.h, &pos.k, &pos.l,
					&pos.ki, &pos.kf, &pos.E, 0,
					&vecQ2);

				if(tl::is_nan_or_inf<t_real>(pos.h)
					|| tl::is_nan_or_inf<t_real>(pos.k)
					|| tl::is_nan_or_inf<t_real>(pos.l))
					throw tl::Err("Invalid h' k' l'.");

				t_vec vecQ2rlu = tl::prod_mv(matUBinv, vecQ2);

				tl::set_eps_0(vecQ2, g_dEps);
				tl::set_eps_0(vecQ2rlu, g_dEps);

				for(t_real* d : { &pos.h, &pos.k, &pos.l,
					&pos.ki, &pos.kf, &pos.E })
					tl::set_eps_0(*d, g_dEps);


				// print inelastic position
				ostrResultsnKimKf << "<tr>";
				ostrResultsnKimKf << "<td>" << row+1 << "</td>";
				if(B_ok && U_ok)
					ostrResultsnKimKf << "<td>" << vecQrlu[0] << ", " << vecQrlu[1] << ", " << vecQrlu[2] << "</td>";
				else
					ostrResultsnKimKf << "<td>" << vecQ[0] << ", " << vecQ[1] << ", " << vecQ[2] << "</td>";
				ostrResultsnKimKf << "<td>" << tl::veclen(vecQ) << "</td>";
				ostrResultsnKimKf << "<td>" << dE << "</td>";

				// print results for n*ki' = m*kf position
				if(B_ok && U_ok)
					ostrResultsnKimKf << "<td>" << vecQ2rlu[0] << ", " << vecQ2rlu[1] << ", " << vecQ2rlu[2] << "</td>";
				else
					ostrResultsnKimKf << "<td>" << vecQ2[0] << ", " << vecQ2[1] << ", " << vecQ2[2] << "</td>";
				ostrResultsnKimKf << "<td>" << tl::veclen(vecQ2) << "</td>";
				ostrResultsnKimKf << "<td>" << pos.E << "</td>";
				ostrResultsnKimKf << "</tr>";

				m_positions_inel_nkimkf.emplace_back(std::move(pos));
			}
			catch(const std::exception& ex)
			{
				ostrResultsnKimKf << "<tr>";
				ostrResultsnKimKf << "<td>" << row+1 << "</td>";
				ostrResultsnKimKf << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				ostrResultsnKimKf << "</tr>";
			}


			try
			{
				// get corresponding analyser angle for n*kf' = m*ki
				ElasticDlgPos pos = pos_inel;
				pos.E = 0;
				pos.kf = m_m/m_n * pos_inel.ki;
				pos.ana_2theta = tl::get_mono_twotheta(pos.kf/angs, GetAnaD()*angs, GetAnaSense()) / rads;

				// find Q position for n*kf' = m*ki position
				t_vec vecQ1;
				tl::get_hkl_from_tas_angles<t_real>(GetLattice(),
					GetScatteringPlaneVec1(), GetScatteringPlaneVec2(),
					GetMonoD(), GetAnaD(),
					pos.mono_2theta*0.5, pos.ana_2theta*0.5,
					pos.sample_theta, pos.sample_2theta,
					GetMonoSense(), GetAnaSense(), GetSampleSense(),
					&pos.h, &pos.k, &pos.l,
					&pos.ki, &pos.kf, &pos.E, 0,
					&vecQ1);

				if(tl::is_nan_or_inf<t_real>(pos.h)
					|| tl::is_nan_or_inf<t_real>(pos.k)
					|| tl::is_nan_or_inf<t_real>(pos.l))
					throw tl::Err("Invalid h' k' l'.");

				t_vec vecQ1rlu = tl::prod_mv(matUBinv, vecQ1);

				tl::set_eps_0(vecQ1, g_dEps);
				tl::set_eps_0(vecQ1rlu, g_dEps);

				for(t_real* d : { &pos.h, &pos.k, &pos.l,
					&pos.ki, &pos.kf, &pos.E })
					tl::set_eps_0(*d, g_dEps);


				// print inelastic position
				ostrResultsnKfmKi << "<tr>";
				ostrResultsnKfmKi << "<td>" << row+1 << "</td>";
				if(B_ok && U_ok)
					ostrResultsnKfmKi << "<td>" << vecQrlu[0] << ", " << vecQrlu[1] << ", " << vecQrlu[2] << "</td>";
				else
					ostrResultsnKfmKi << "<td>" << vecQ[0] << ", " << vecQ[1] << ", " << vecQ[2] << "</td>";
				ostrResultsnKfmKi << "<td>" << tl::veclen(vecQ) << "</td>";
				ostrResultsnKfmKi << "<td>" << dE << "</td>";

				// print results for n*kf' = m*ki position
				if(B_ok && U_ok)
					ostrResultsnKfmKi << "<td>" << vecQ1rlu[0] << ", " << vecQ1rlu[1] << ", " << vecQ1rlu[2] << "</td>";
				else
					ostrResultsnKfmKi << "<td>" << vecQ1[0] << ", " << vecQ1[1] << ", " << vecQ1[2] << "</td>";
				ostrResultsnKfmKi << "<td>" << tl::veclen(vecQ1) << "</td>";
				ostrResultsnKfmKi << "<td>" << pos.E << "</td>";
				ostrResultsnKfmKi << "</tr>";

				m_positions_inel_nkfmki.emplace_back(std::move(pos));
			}
			catch(const std::exception& ex)
			{
				ostrResultsnKfmKi << "<tr>";
				ostrResultsnKfmKi << "<td>" << row+1 << "</td>";
				ostrResultsnKfmKi << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				ostrResultsnKfmKi << "</tr>";
			}


			m_positions_inel.emplace_back(std::move(pos_inel));
		}
		catch(const std::exception& ex)
		{
			for(std::wostream* ostr : { &ostrResultsKfKi, &ostrResultsKiKf,
				&ostrResultsnKimKf, &ostrResultsnKfmKi })
			{
				(*ostr) << "<tr>";
				(*ostr) << "<td>" << row+1 << "</td>";
				(*ostr) << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				(*ostr) << "</tr>";
			}
		}
	}

	ostrResultsKfKi << "</table></center>";
	ostrResultsKiKf << "</table></center>";
	ostrResultsnKimKf << "</table></center>";
	ostrResultsnKfmKi << "</table></center>";

	std::wostringstream ostrResults;
	ostrResults << "<html><body>";
	ostrResults << "<p>" << ostrResultsKfKi.str() << "</p>";
	ostrResults << "<br>";
	ostrResults << "<p>" << ostrResultsKiKf.str() << "</p>";
	ostrResults << "<br>";
	ostrResults << "<p>" << ostrResultsnKimKf.str() << "</p>";
	ostrResults << "<br>";
	ostrResults << "<p>" << ostrResultsnKfmKi.str() << "</p>";
	ostrResults << "</body></html>";
	textResults->setHtml(QString::fromWCharArray(ostrResults.str().c_str()));
}


/**
 * set the scattering triangle to the given inelastic or elastic position
 */
void ElasticDlg::GotoPosition(ElasticDlgGoto which)
{
	const std::vector<ElasticDlgPos>* positions = nullptr;
	switch(which)
	{
		case ElasticDlgGoto::ELAST_KFKI:
			positions = &m_positions_elast_kfki;
			break;
		case ElasticDlgGoto::ELAST_KIKF:
			positions = &m_positions_elast_kikf;
			break;
		case ElasticDlgGoto::INEL_nKFmKI:
			positions = &m_positions_inel_nkfmki;
			break;
		case ElasticDlgGoto::INEL_nKImKF:
			positions = &m_positions_inel_nkimkf;
			break;
		case ElasticDlgGoto::INEL_ORIG:
		default:
			positions = &m_positions_inel;
			break;
	}

	if(!m_bSyncWithMainWindow)
	{
		QMessageBox::critical(this, "Error",
			"Not synchronised with main sample / instrument parameters. Please click on \"Sync\" to do so.");
		return;
	}

	int pos_idx = spinPosIdx->value() - 1;
	if(pos_idx >= (int)positions->size() || pos_idx < 0)
	{
		QMessageBox::critical(this, "Error", "Invalid position index.");
		return;
	}

	const ElasticDlgPos& pos = (*positions)[pos_idx];

	CrystalOptions crys;
	TriangleOptions triag;

	triag.bChangedMonoTwoTheta = true;
	triag.dMonoTwoTheta = pos.mono_2theta;
	triag.bChangedAnaTwoTheta = true;
	triag.dAnaTwoTheta = pos.ana_2theta;
	triag.bChangedTheta = true;
	triag.dTheta = pos.sample_theta;
	triag.bChangedAngleKiVec0 = true;
	triag.dAngleKiVec0 = tl::get_pi<t_real>()/2. - pos.sample_theta;
	triag.bChangedTwoTheta = true;
	triag.dTwoTheta = pos.sample_2theta;

	emit ChangedPosition(crys, triag);
}



// ----------------------------------------------------------------------------
// positions table
// ----------------------------------------------------------------------------
/**
 * add an empty inelastic position to the table
 */
void ElasticDlg::AddPosition()
{
	// add a row
	int row = tablePositions->rowCount();
	tablePositions->setRowCount(row + 1);

	// insert items with default values
	for(int col=0; col<POSTAB_COLS; ++col)
		tablePositions->setItem(row, col, new QTableWidgetItem("0"));
}


/**
 * add a given inelastic position to the table
 */
void ElasticDlg::AddPosition(t_real h, t_real k, t_real l, t_real ki, t_real kf)
{
	AddPosition();

	tablePositions->item(tablePositions->rowCount() - 1, POSTAB_H)->setText(
		tl::var_to_str(h, g_iPrec).c_str());
	tablePositions->item(tablePositions->rowCount() - 1, POSTAB_K)->setText(
		tl::var_to_str(k, g_iPrec).c_str());
	tablePositions->item(tablePositions->rowCount() - 1, POSTAB_L)->setText(
		tl::var_to_str(l, g_iPrec).c_str());
	tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KI)->setText(
		tl::var_to_str(ki, g_iPrec).c_str());
	tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KF)->setText(
		tl::var_to_str(kf, g_iPrec).c_str());
}


/**
 * delete the selected inelastic positions from the table
 */
void ElasticDlg::DelPosition()
{
	int row = tablePositions->currentRow();
	if(row >= 0)
	{
		tablePositions->removeRow(row);
	}
	else
	{
		tablePositions->clearContents();
		tablePositions->setRowCount(0);
	}

	CalcSpuriousPositions();
}


/**
 * create the position generation dialog
 */
void ElasticDlg::InitGeneratePositionsDlg()
{
	if(!m_pGenPosDlg)
	{
		m_pGenPosDlg = new GenPosDlg(this, m_pSettings);
		QObject::connect(m_pGenPosDlg, &GenPosDlg::GeneratedPositions,
			this, &ElasticDlg::GeneratedPositions);
	}
}


/**
 * show dialog to generate positions
 */
void ElasticDlg::GeneratePositions()
{
	InitGeneratePositionsDlg();
	focus_dlg(m_pGenPosDlg);
}


/**
 * add received generated positions
 */
void ElasticDlg::GeneratedPositions(const std::vector<ScanPosition>& positions)
{
	m_bAllowCalculation = false;
	BOOST_SCOPE_EXIT(this_, &m_bAllowCalculation)
	{
		// recalculate
		m_bAllowCalculation = true;
		this_->CalcSpuriousPositions();
	} BOOST_SCOPE_EXIT_END

	// clear old positions
	tablePositions->clearContents();
	tablePositions->setRowCount(0);

	for(const ScanPosition& pos : positions)
	{
		AddPosition();

		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_H)->setText(
			tl::var_to_str(pos.h, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_K)->setText(
			tl::var_to_str(pos.k, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_L)->setText(
			tl::var_to_str(pos.l, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KI)->setText(
			tl::var_to_str(pos.ki, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KF)->setText(
			tl::var_to_str(pos.kf, g_iPrec).c_str());
	}
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// importing positions from scan files
// ----------------------------------------------------------------------------
/**
 * import positions from a scan file
 */
void ElasticDlg::ImportPositions()
{
	std::vector<std::string> vecFiles = GetFiles();
	if(!vecFiles.size() || vecFiles[0]=="")
		return;

	m_bAllowCalculation = false;
	BOOST_SCOPE_EXIT(this_, &m_bAllowCalculation)
	{
		// recalculate
		m_bAllowCalculation = true;
		this_->CalcSpuriousPositions();
	} BOOST_SCOPE_EXIT_END

	// clear old positions
	tablePositions->clearContents();
	tablePositions->setRowCount(0);
	btnSync->setChecked(false);

	// iterate scan files
	bool first_file = true;
	for(const std::string& strFile : vecFiles)
	{
		std::unique_ptr<tl::FileInstrBase<t_real>> pScan(
			tl::FileInstrBase<t_real>::LoadInstr(strFile.c_str()));
		if(!pScan)
		{
			tl::log_err("Invalid scan file: \"", strFile, "\".");
			continue;
		}

		if(first_file)
		{
			// load external crystal and instrument configuration from first file
			std::array<t_real, 3> latt = pScan->GetSampleLattice();
			std::array<t_real, 3> angles = pScan->GetSampleAngles();
			m_cfg_ext.m_lattice = tl::Lattice<t_real>(
				latt[0], latt[1], latt[2],
				angles[0], angles[1], angles[2]);

			std::array<bool, 3> senses = pScan->GetScatterSenses();
			m_cfg_ext.m_bSenses[0] = senses[0];
			m_cfg_ext.m_bSenses[1] = senses[1];
			m_cfg_ext.m_bSenses[2] = senses[2];

			std::array<t_real, 3> vec1 = pScan->GetScatterPlane0();
			std::array<t_real, 3> vec2 = pScan->GetScatterPlane1();
			m_cfg_ext.m_vec1 = tl::make_vec<t_vec>({ vec1[0], vec1[1], vec1[2] });
			m_cfg_ext.m_vec2 = tl::make_vec<t_vec>({ vec2[0], vec2[1], vec2[2] });

			std::array<t_real, 2> ds = pScan->GetMonoAnaD();
			m_cfg_ext.m_dMono = ds[0];
			m_cfg_ext.m_dAna = ds[1];

			first_file = false;
		}

		// iterate scan positions and add them to the table
		for(std::size_t idx = 0; idx < pScan->GetScanCount(); ++idx)
		{
			auto pos = pScan->GetScanHKLKiKf(idx);
			AddPosition(std::get<0>(pos), std::get<1>(pos), std::get<2>(pos),
				std::get<3>(pos), std::get<4>(pos));
		}
	}
}


/**
 * opens a file selection dialog to select scan files
 */
std::vector<std::string> ElasticDlg::GetFiles()
{
	std::vector<std::string> vecFiles;

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	bool bShowPreview = true;
	QString strDirLast;
	if(m_pSettings)
	{
		bShowPreview = m_pSettings->value("main/dlg_previews", true).toBool();
		strDirLast = m_pSettings->value("elastic_pos/last_import_dir", ".").toString();
		if(!m_pSettings->value("main/native_dialogs", true).toBool())
			fileopt = QFileDialog::DontUseNativeDialog;
	}

	QFileDialog *pdlg = nullptr;
	if(bShowPreview)
	{
		pdlg = new FilePreviewDlg(this, "Import Data File...", m_pSettings);
		pdlg->setOptions(QFileDialog::DontUseNativeDialog);
	}
	else
	{
		pdlg = new QFileDialog(this, "Import Data File...");
		pdlg->setOptions(fileopt);
	}
	std::unique_ptr<QFileDialog> ptrdlg(pdlg);

	pdlg->setDirectory(strDirLast);
	pdlg->setFileMode(QFileDialog::ExistingFiles);
	pdlg->setViewMode(QFileDialog::Detail);
#if !defined NO_IOSTR
	QString strFilter = "Data files (*.dat *.scn *.DAT *.SCN *.ng0 *.NG0 *.log *.LOG *.scn.gz *.SCN.GZ *.dat.gz *.DAT.GZ *.ng0.gz *.NG0.GZ *.log.gz *.LOG.GZ *.scn.bz2 *.SCN.BZ2 *.dat.bz2 *.DAT.BZ2 *.ng0.bz2 *.NG0.BZ2 *.log.bz2 *.LOG.BZ2);;All files (*.*)";
#else
	QString strFilter = "Data files (*.dat *.scn *.DAT *.SCN *.NG0 *.ng0 *.log *.LOG);;All files (*.*)";
#endif
	pdlg->setNameFilter(strFilter);
	if(!pdlg->exec())
		return vecFiles;
	if(!pdlg->selectedFiles().size())
		return vecFiles;

	vecFiles.reserve(pdlg->selectedFiles().size());
	for(int iFile=0; iFile<pdlg->selectedFiles().size(); ++iFile)
		vecFiles.push_back(pdlg->selectedFiles()[iFile].toStdString());

	if(m_pSettings && vecFiles[0] != "")
	{
		std::string strDir = tl::get_dir(vecFiles[0]);
		m_pSettings->setValue("elastic_pos/last_import_dir", QString(strDir.c_str()));
	}

	return vecFiles;
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// loading / saving
// ----------------------------------------------------------------------------
/**
 * save the configuration
 */
void ElasticDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	// save positions
	for(int row = 0; row < tablePositions->rowCount(); ++row)
	{
		QTableWidgetItem *item_h = tablePositions->item(row, POSTAB_H);
		QTableWidgetItem *item_k = tablePositions->item(row, POSTAB_K);
		QTableWidgetItem *item_l = tablePositions->item(row, POSTAB_L);
		QTableWidgetItem *item_ki = tablePositions->item(row, POSTAB_KI);
		QTableWidgetItem *item_kf = tablePositions->item(row, POSTAB_KF);

		if(!item_h || !item_k || !item_l || !item_ki || !item_kf)
			continue;

		std::ostringstream ostrPos;
		ostrPos << "elastic_pos/pos_" << row << "/";
		std::string strPos = ostrPos.str();

		mapConf[strXmlRoot + strPos + "h"] = item_h->text().toStdString();
		mapConf[strXmlRoot + strPos + "k"] = item_k->text().toStdString();
		mapConf[strXmlRoot + strPos + "l"] = item_l->text().toStdString();
		mapConf[strXmlRoot + strPos + "ki"] = item_ki->text().toStdString();
		mapConf[strXmlRoot + strPos + "kf"] = item_kf->text().toStdString();
	}

	mapConf[strXmlRoot + "elastic_pos/inel_n"] = tl::var_to_str(spinN->value(), g_iPrec).c_str();
	mapConf[strXmlRoot + "elastic_pos/inel_m"] = tl::var_to_str(spinM->value(), g_iPrec).c_str();
	mapConf[strXmlRoot + "elastic_pos/pos_idx"] = tl::var_to_str(spinPosIdx->value()).c_str();

	if(m_pGenPosDlg)
		m_pGenPosDlg->Save(mapConf, strXmlRoot + "elastic_pos/");
}


/**
 * load the configuration
 */
void ElasticDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	m_bAllowCalculation = false;
	BOOST_SCOPE_EXIT(this_, &m_bAllowCalculation)
	{
		// recalculate
		m_bAllowCalculation = true;
		this_->CalcSpuriousPositions();
	} BOOST_SCOPE_EXIT_END

	bool bOk = false;

	// clear old positions
	tablePositions->clearContents();
	tablePositions->setRowCount(0);
	btnSync->setChecked(true);

	// load positions
	std::size_t row = 0;
	while(true)
	{
		std::ostringstream ostrPos;
		ostrPos << "elastic_pos/pos_" << row << "/";
		std::string strPos = ostrPos.str();

		// last position reached?
		if(!xml.Exists(strXmlRoot + strPos))
			break;

		t_real dh = xml.Query<t_real>(strXmlRoot + strPos + "h", 0., &bOk);
		t_real dk = xml.Query<t_real>(strXmlRoot + strPos + "k", 0., &bOk);
		t_real dl = xml.Query<t_real>(strXmlRoot + strPos + "l", 0., &bOk);
		t_real dki = xml.Query<t_real>(strXmlRoot + strPos + "ki", 0., &bOk);
		t_real dkf = xml.Query<t_real>(strXmlRoot + strPos + "kf", 0., &bOk);

		AddPosition(dh, dk, dl, dki, dkf);
		++row;
	}

	spinN->setValue(xml.Query<double>(strXmlRoot + "elastic_pos/inel_n", 1., &bOk));
	spinM->setValue(xml.Query<double>(strXmlRoot + "elastic_pos/inel_m", 2., &bOk));
	spinPosIdx->setValue(xml.Query<int>(strXmlRoot + "elastic_pos/pos_idx", 1, &bOk));

	// restore the settings from the position generation dialog
	if(xml.Exists(strXmlRoot + "elastic_pos/gen_pos"))
	{
		InitGeneratePositionsDlg();
		m_pGenPosDlg->Load(xml, strXmlRoot + "elastic_pos/");
	}
}
// ----------------------------------------------------------------------------


void ElasticDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("elastic_pos/geo", saveGeometry());

		QDialog::accept();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}
}


void ElasticDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


// ----------------------------------------------------------------------------
// setter and getter
// ----------------------------------------------------------------------------
void ElasticDlg::SetLattice(const tl::Lattice<t_real>& lattice)
{
	m_cfg.m_lattice = lattice;
}


void ElasticDlg::SetScatteringPlane(const t_vec& vec1, const t_vec& vec2)
{
	m_cfg.m_vec1 = vec1;
	m_cfg.m_vec2 = vec2;
}


void ElasticDlg::SetD(t_real dMono, t_real dAna)
{
	m_cfg.m_dMono = dMono;
	m_cfg.m_dAna = dAna;
}


void ElasticDlg::SetMonoSense(bool bSense)
{
	m_cfg.m_bSenses[0] = bSense;
}


void ElasticDlg::SetSampleSense(bool bSense)
{
	m_cfg.m_bSenses[1] = bSense;
}


void ElasticDlg::SetAnaSense(bool bSense)
{
	m_cfg.m_bSenses[2] = bSense;
}


void ElasticDlg::SetSenses(bool bM, bool bS, bool bA)
{
	m_cfg.m_bSenses[0] = bM;
	m_cfg.m_bSenses[1] = bS;
	m_cfg.m_bSenses[2] = bA;
}


const tl::Lattice<t_real>& ElasticDlg::GetLattice() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_lattice;
}


const t_vec& ElasticDlg::GetScatteringPlaneVec1() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_vec1;
}


const t_vec& ElasticDlg::GetScatteringPlaneVec2() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_vec2;
}


t_real ElasticDlg::GetMonoD() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_dMono;
}


t_real ElasticDlg::GetAnaD() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_dAna;
}


bool ElasticDlg::GetMonoSense() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_bSenses[0];
}


bool ElasticDlg::GetSampleSense() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_bSenses[1];
}


bool ElasticDlg::GetAnaSense() const
{
	const ElasticDlgCfg *cfg = m_bSyncWithMainWindow ? &m_cfg : &m_cfg_ext;
	return cfg->m_bSenses[2];
}
// ----------------------------------------------------------------------------


#include "moc_ElasticDlg.cpp"
