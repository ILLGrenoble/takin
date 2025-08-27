/**
 * Generate Scan Positions
 * @author Tobias Weber <tweber@ill.fr>
 * @date 8-jan-2024
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

#include "GenPosDlg.h"

#include "tlibs/phys/neutrons.h"
#include "tlibs/string/string.h"


using t_real = t_real_glob;

static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();


GenPosDlg::GenPosDlg(QWidget* pParent, QSettings* pSett)
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

		if(m_pSettings->contains("gen_pos/geo"))
			restoreGeometry(m_pSettings->value("gen_pos/geo").toByteArray());
	}

	// connections
	QObject::connect(btnGenerate, &QAbstractButton::clicked,
		this, &GenPosDlg::GeneratePositions);
	QObject::connect(buttonBox, &QDialogButtonBox::clicked,
		this, &GenPosDlg::ButtonBoxClicked);
#if QT_VERSION >= QT_VERSION_CHECK(6, 7, 0)
	QObject::connect(checkCKf, &QCheckBox::checkStateChanged, [this](Qt::CheckState state)
#else
	QObject::connect(checkCKf, &QCheckBox::stateChanged, [this](int state)
#endif
	{
		if(state == Qt::Unchecked)
			labelFixedK->setText("ki (Å⁻¹):");
		else
			labelFixedK->setText("kf (Å⁻¹):");
	});
}


GenPosDlg::~GenPosDlg()
{
}


/**
 * generate the scan positions and send them to the listeners
 */
void GenPosDlg::GeneratePositions()
{
	t_real hi = tl::str_to_var<t_real>(editHi->text().toStdString());
	t_real ki = tl::str_to_var<t_real>(editKi->text().toStdString());
	t_real li = tl::str_to_var<t_real>(editLi->text().toStdString());
	t_real Ei = tl::str_to_var<t_real>(editEi->text().toStdString());
	t_real hf = tl::str_to_var<t_real>(editHf->text().toStdString());
	t_real kf = tl::str_to_var<t_real>(editKf->text().toStdString());
	t_real lf = tl::str_to_var<t_real>(editLf->text().toStdString());
	t_real Ef = tl::str_to_var<t_real>(editEf->text().toStdString());

	t_real kfix = tl::str_to_var<t_real>(editKfix->text().toStdString());
	bool kf_fixed = checkCKf->isChecked();

	std::vector<ScanPosition> positions;
	int steps = spinSteps->value();
	positions.reserve(steps);

	for(int step = 0; step < steps; ++step)
	{
		ScanPosition pos;

		pos.h = tl::lerp(hi, hf, t_real(step)/t_real(steps-1));
		pos.k = tl::lerp(ki, kf, t_real(step)/t_real(steps-1));
		pos.l = tl::lerp(li, lf, t_real(step)/t_real(steps-1));
		pos.E = tl::lerp(Ei, Ef, t_real(step)/t_real(steps-1));

		if(kf_fixed)
		{
			pos.kf = kfix;
			pos.ki = tl::get_other_k(pos.E*meV, kfix/angs, !kf_fixed) * angs;
		}
		else
		{
			pos.ki = kfix;
			pos.kf = tl::get_other_k(pos.E*meV, kfix/angs, !kf_fixed) * angs;
		}

		positions.emplace_back(std::move(pos));
	}

	emit GeneratedPositions(positions);
}


/**
 * save configuration
 */
void GenPosDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	mapConf[strXmlRoot + "gen_pos/hi"] = editHi->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/ki"] = editKi->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/li"] = editLi->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/Ei"] = editEi->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/hf"] = editHf->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/kf"] = editKf->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/lf"] = editLf->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/Ef"] = editEf->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/ckf"] = checkCKf->isChecked() ? "1" : "0";
	mapConf[strXmlRoot + "gen_pos/kfix"] = editKfix->text().toStdString();
	mapConf[strXmlRoot + "gen_pos/steps"] = tl::var_to_str(spinSteps->value()).c_str();
}


/**
 * load configuration
 */
void GenPosDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	bool bOk;
	editHi->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/hi", 1., &bOk), g_iPrec).c_str());
	editKi->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/ki", 0., &bOk), g_iPrec).c_str());
	editLi->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/li", 0., &bOk), g_iPrec).c_str());
	editEi->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/Ei", 0, &bOk), g_iPrec).c_str());
	editHf->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/hf", 1., &bOk), g_iPrec).c_str());
	editKf->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/kf", 0., &bOk), g_iPrec).c_str());
	editLf->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/lf", 0., &bOk), g_iPrec).c_str());
	editEf->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/Ef", 1, &bOk), g_iPrec).c_str());
	checkCKf->setChecked(xml.Query<bool>(strXmlRoot + "gen_pos/ckf", true, &bOk));
	editKfix->setText(tl::var_to_str(xml.Query<t_real>(strXmlRoot + "gen_pos/kfix", 1.4, &bOk), g_iPrec).c_str());
	spinSteps->setValue(xml.Query<int>(strXmlRoot + "gen_pos/steps", 16, &bOk));
}


void GenPosDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("gen_pos/geo", saveGeometry());

		QDialog::accept();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}
}


#include "moc_GenPosDlg.cpp"
