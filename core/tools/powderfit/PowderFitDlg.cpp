/**
 * Powder Fit Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date Mar-2018
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

#include <QMessageBox>
#include <QFileDialog>
#include <QPushButton>

#include <fstream>

#include "PowderFitDlg.h"
#include "tlibs/string/string.h"
#include "tlibs/phys/powder.h"
#include "tlibs/string/spec_char.h"
#include "libs/globals.h"

using t_real = t_real_glob;


PowderFitDlg::PowderFitDlg(QWidget* pParent, QSettings *pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	setupUi(this);

	buttonBox->addButton(new QPushButton("Plot..."), QDialogButtonBox::ActionRole);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("powderfit/geo"))
			restoreGeometry(m_pSettings->value("powderfit/geo").toByteArray());

		// load recent values
		if(m_pSettings->contains("powderfit/d"))
			editD->setText(m_pSettings->value("powderfit/d").toString());
		if(m_pSettings->contains("powderfit/ki"))
			editKi->setText(m_pSettings->value("powderfit/ki").toString());
		if(m_pSettings->contains("powderfit/Gs"))
			editGs->setText(m_pSettings->value("powderfit/Gs").toString());
		if(m_pSettings->contains("powderfit/tt"))
			editAngles->setText(m_pSettings->value("powderfit/tt").toString());
		if(m_pSettings->contains("powderfit/dtt"))
			editAngleErrs->setText(m_pSettings->value("powderfit/dtt").toString());
	}

	QObject::connect(buttonBox, &QDialogButtonBox::clicked, this, &PowderFitDlg::ButtonBoxClicked);

	for(QLineEdit* pEdit : {editD, editKi, editGs, editAngles, editAngleErrs})
		QObject::connect(pEdit, &QLineEdit::textChanged, this, &PowderFitDlg::Calc);

	Calc();
}


/**
 * saves the plot script
 */
void PowderFitDlg::SavePlot()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	QString strDirLast;
	if(m_pSettings)
	{
		if(!m_pSettings->value("main/native_dialogs", 1).toBool())
			fileopt = QFileDialog::DontUseNativeDialog;
		strDirLast = m_pSettings->value("powderfit/last_save_dir", ".").toString();
	}

	std::string strFile = QFileDialog::getSaveFileName(this,
		"Save Plot Script", strDirLast, "Plot files (*.gpl *.GPL)",
		nullptr, fileopt).toStdString();

	if(strFile != "")
	{
		std::ofstream ofstr(strFile);
		if(!ofstr)
		{
			QMessageBox::critical(this, "Error", "Cannot save plot script.");
			return;
		}

		std::string strScr = editScript->toPlainText().toStdString();
		ofstr << strScr;
		ofstr.close();

		if(m_pSettings)
		{
			std::string strDir = tl::get_dir(strFile);
			m_pSettings->setValue("powderfit/last_save_dir", strDir.c_str());
		}
	}
}


/**
 * fits the powder lines
 */
void PowderFitDlg::Calc()
{
	try
	{
		t_real dMono = tl::str_to_var<t_real>(editD->text().toStdString());
		t_real dKi = tl::str_to_var<t_real>(editKi->text().toStdString());

		std::string strGs = editGs->text().toStdString();
		std::string strTTs = editAngles->text().toStdString();
		std::string strTTErrs = editAngleErrs->text().toStdString();

		std::vector<t_real> vecGs, vecTTs, vecTTErrs;
		tl::get_tokens<t_real, std::string>(strGs, ",;", vecGs);
		tl::get_tokens<t_real, std::string>(strTTs, ",;", vecTTs);
		tl::get_tokens<t_real, std::string>(strTTErrs, ",;", vecTTErrs);

		for(t_real& dTT : vecTTs) dTT = tl::d2r(dTT);
		for(t_real& dTTErr : vecTTErrs) dTTErr = tl::d2r(dTTErr);

		if(vecGs.size()!=vecTTs.size() || vecTTs.size()!=vecTTErrs.size())
		{
			editResults->setPlainText("Error: Same number of input values required!");
			return;
		}


		std::vector<t_real> vecRes = { dKi, 0., 0., 0. };
		std::vector<t_real> vecResErr = { dKi*0.25, 0.5, vecRes[2]*0.25, vecRes[3]*0.25,};

		bool bFitOk = 0;
		try
		{
			bFitOk = tl::powder_align<t_real>(dMono, vecGs, vecTTs, vecTTErrs, vecRes, vecResErr);
		}
		catch(const std::exception& ex)
		{
			tl::log_err(ex.what());
		}


		// write results
		{
			std::ostringstream ostrResults;
			ostrResults.precision(g_iPrec);

			std::string strAA_1 = tl::get_spec_char_utf8("AA") +
				tl::get_spec_char_utf8("sup-") + tl::get_spec_char_utf8("sup1");
			std::string strDeg = tl::get_spec_char_utf8("deg");
			std::string strTT = std::string("2") + tl::get_spec_char_utf8("theta");
			std::string strDTT = tl::get_spec_char_utf8("Delta") + strTT;

			ostrResults << "Fit converged: " << std::boolalpha << bFitOk << "\n";
			ostrResults << "ki = " << vecRes[0] << " " << strAA_1 << "\n";
			ostrResults << strTT << "m = " << tl::r2d(vecRes[2]) << strDeg << "\n";
			ostrResults << strDTT << "m = " << tl::r2d(vecRes[3]) << strDeg << "\n";
			ostrResults << strDTT << "s = " << tl::r2d(vecRes[1]) << strDeg << "\n";

			editResults->setPlainText(QString::fromUtf8(ostrResults.str().c_str()));
		}


		// create plot script
		{
			std::ostringstream ostrPlot;
			ostrPlot.precision(g_iPrecGfx);

			ostrPlot << "set xlabel \"Wavenumber G (1/A)\"\n";
			ostrPlot << "set ylabel \"Scattering angle (deg)\"\n";
			ostrPlot << "unset key\n";

			ostrPlot << "\nk = " << vecRes[0] << "\n";
			ostrPlot << "dtt = " << tl::r2d(vecRes[1]) << "\n";


			ostrPlot << "\nplot \\"
				<< "\n\t 2.*asin(x/(2*k))/pi*180. + dtt, \\"
				<< "\n\t\"-\" u 1:2:3 w yerrorbars pt 7\n";

			for(std::size_t i=0; i<vecGs.size(); ++i)
			{
				ostrPlot << vecGs[i] << "\t" << tl::r2d(vecTTs[i])
					<< "\t" << tl::r2d(vecTTErrs[i]) << "\n";
			}
			ostrPlot << "e\n";

			editScript->setPlainText(ostrPlot.str().c_str());
		}
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


void PowderFitDlg::Plot()
{
	try
	{
		std::string strVer = get_gpltool_version();
		tl::log_info("Invoking ", strVer, ".");

		if(!m_pPlotProc)
			m_pPlotProc.reset(new tl::PipeProc<char>(("\"" + g_strGplTool + "\" -p 2>/dev/null 1>/dev/null").c_str(), 1));

		if(strVer == "" || !m_pPlotProc || !m_pPlotProc->IsReady())
		{
			QMessageBox::critical(this, "Error", "Gnuplot cannot be invoked.");
			return;
		}

		(*m_pPlotProc) << editScript->toPlainText().toStdString();
		m_pPlotProc->flush();
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Critical Error", ex.what());
	}
}


void PowderFitDlg::accept()
{
	// save recent values
	if(m_pSettings)
	{
		m_pSettings->setValue("powderfit/d", editD->text());
		m_pSettings->setValue("powderfit/ki", editKi->text());
		m_pSettings->setValue("powderfit/Gs", editGs->text());
		m_pSettings->setValue("powderfit/tt", editAngles->text());
		m_pSettings->setValue("powderfit/dtt", editAngleErrs->text());
	}
}


void PowderFitDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->standardButton(pBtn) == QDialogButtonBox::Save)
	{
		SavePlot();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::ActionRole)
	{
		Plot();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("powderfit/geo", saveGeometry());
		QDialog::accept();
	}
}


#include "moc_PowderFitDlg.cpp"
