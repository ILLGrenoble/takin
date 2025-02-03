/**
 * Scan Position Plotter Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date Feb-2017
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

#include <QMessageBox>
#include <QFileDialog>

#include "dialogs/FilePreviewDlg.h"
#include "ScanPosDlg.h"
#include "scanpos.h"
#include "tlibs/string/string.h"

using t_real = t_real_glob;
using t_vec = tl::ublas::vector<t_real>;
using t_mat = tl::ublas::matrix<t_real>;


ScanPosDlg::ScanPosDlg(QWidget* pParent, QSettings *pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	setupUi(this);

	buttonBox->addButton(new QPushButton("Plot..."), QDialogButtonBox::ActionRole);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("scanpos/geo"))
			restoreGeometry(m_pSettings->value("scanpos/geo").toByteArray());
	}

	QObject::connect(buttonBox, &QDialogButtonBox::clicked, this, &ScanPosDlg::ButtonBoxClicked);
	QObject::connect(editBragg, &QLineEdit::textChanged, this, &ScanPosDlg::UpdatePlot);
	QObject::connect(editPlane, &QLineEdit::textChanged, this, &ScanPosDlg::UpdatePlot);
	QObject::connect(editPos, &QPlainTextEdit::textChanged, this, &ScanPosDlg::UpdatePlot);
	QObject::connect(checkFlip, &QCheckBox::toggled, this, &ScanPosDlg::UpdatePlot);
	QObject::connect(check3D, &QCheckBox::toggled, this, &ScanPosDlg::UpdatePlot);
	QObject::connect(btnFromScan, &QPushButton::clicked, this, &ScanPosDlg::LoadPlaneFromFile);
	QObject::connect(btnAddFromScan, &QPushButton::clicked, this, &ScanPosDlg::LoadPosFromFile);
}


/**
 * opens a file selection dialog to select scan files
 */
std::vector<std::string> ScanPosDlg::GetFiles(bool bMultiple)
{
	std::vector<std::string> vecFiles;

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	bool bShowPreview = true;
	QString strDirLast;
	if(m_pSettings)
	{
		bShowPreview = m_pSettings->value("main/dlg_previews", true).toBool();
		strDirLast = m_pSettings->value("scanpos/last_import_dir", ".").toString();
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
	pdlg->setFileMode(bMultiple ? QFileDialog::ExistingFiles : QFileDialog::ExistingFile);
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
		m_pSettings->setValue("scanpos/last_import_dir", QString(strDir.c_str()));
	}

	return vecFiles;
}


/**
 * gets the scattering plane from a scan file
 */
void ScanPosDlg::LoadPlaneFromFile()
{
	std::vector<std::string> vecFiles = GetFiles(false);
	if(!vecFiles.size() || vecFiles[0]=="")
		return;

	std::unique_ptr<tl::FileInstrBase<t_real>> ptrScan(
		tl::FileInstrBase<t_real>::LoadInstr(vecFiles[0].c_str()));
	if(!ptrScan)
	{
		tl::log_err("Invalid scan file: \"", vecFiles[0], "\".");
		return;
	}

	auto arrVec0 = ptrScan->GetScatterPlane0();
	t_vec vec0 = tl::make_vec<t_vec>({arrVec0[0],arrVec0[1],arrVec0[2]});
	auto arrVec1 = ptrScan->GetScatterPlane1();
	t_vec vec1 = tl::make_vec<t_vec>({arrVec1[0],arrVec1[1],arrVec1[2]});

	if(checkFlip->isChecked()) vec1 = -vec1;
	if(vec0.size() != 3) vec0.resize(3);
	if(vec1.size() != 3) vec1.resize(3);

	std::ostringstream ostrPlane;
	ostrPlane.precision(g_iPrec);
	ostrPlane << vec0[0] << ", " << vec0[1] << ", " << vec0[2] << "; "
		<< vec1[0] << ", " << vec1[1] << ", " << vec1[2];
	editPlane->setText(ostrPlane.str().c_str());


	// first scan position
	auto vecPos = get_coord<t_vec, t_real>(vec0, vec1, *ptrScan);
	if(vecPos.first.size() != 3)
	{
		tl::log_err("Invalid scan position for file \"", vecFiles[0], "\".");
		return;
	}

	// guess Bragg peak
	t_vec vecBraggHKL = tl::make_vec<t_vec>(
		{std::round(vecPos.first[0]), std::round(vecPos.first[1]), std::round(vecPos.first[2])});
	std::ostringstream ostrBragg;
	ostrBragg.precision(g_iPrec);
	ostrBragg << vecBraggHKL[0] << ", " << vecBraggHKL[1] << ", " << vecBraggHKL[2];
	editBragg->setText(ostrBragg.str().c_str());
}


/**
 * inserts scan position from data files
 */
void ScanPosDlg::LoadPosFromFile()
{
	std::vector<std::string> vecFiles = GetFiles(true);
	if(!vecFiles.size() || vecFiles[0]=="")
		return;

	t_vec vec0, vec1; // dummies
	vec0 = vec1 = tl::ublas::zero_vector<t_real>(3);

	std::ostringstream ostrPos;
	ostrPos.precision(g_iPrec);

	for(const std::string& strFile : vecFiles)
	{
		std::unique_ptr<tl::FileInstrBase<t_real>> ptrScanOther(
			tl::FileInstrBase<t_real>::LoadInstr(strFile.c_str()));
		if(!ptrScanOther)
		{
			tl::log_err("Invalid scan file: \"", strFile, "\".");
			continue;
		}

		// scan positions
		auto vecPosOther = get_coord<t_vec, t_real>(vec0, vec1, *ptrScanOther);
		if(vecPosOther.first.size() != 3)
		{
			tl::log_err("Invalid scan position for file \"", strFile, "\".");
			continue;
		}

		ostrPos << vecPosOther.first[0] << ", " << vecPosOther.first[1] << ", " << vecPosOther.first[2] << "; ";
	}

	QString strText = editPos->toPlainText();
	strText += ostrPos.str().c_str();
	editPos->setPlainText(strText);
}


/**
 * saves the plot script
 */
void ScanPosDlg::SavePlot()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	QString strDirLast;
	if(m_pSettings)
	{
		if(!m_pSettings->value("main/native_dialogs", 1).toBool())
			fileopt = QFileDialog::DontUseNativeDialog;
		strDirLast = m_pSettings->value("scanpos/last_save_dir", ".").toString();
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
			m_pSettings->setValue("scanpos/last_save_dir", strDir.c_str());
		}
	}
}


/**
 * generates the plot script
 */
void ScanPosDlg::UpdatePlot()
{
	try
	{
		std::string strBragg = editBragg->text().toStdString();
		std::string strPlane = editPlane->text().toStdString();
		std::string strAllPos = editPos->toPlainText().toStdString();
		bool bFlip = checkFlip->isChecked();
		bool b3D = check3D->isChecked();


		// Bragg reflection
		std::vector<t_real> _vecBragg;
		tl::parse_tokens<t_real, std::string>(strBragg, ",", _vecBragg);
		t_vec vecBragg = tl::convert_vec_full<t_real, t_real, std::vector, tl::ublas::vector>(_vecBragg);
		if(vecBragg.size() != 3) vecBragg.resize(3);


		// Scattering plane
		auto pairPlane = tl::split_first<std::string>(strPlane, ";", 1, 0);
		std::vector<t_real> _vec1, _vec2;
		tl::parse_tokens<t_real, std::string>(pairPlane.first, ",", _vec1);
		tl::parse_tokens<t_real, std::string>(pairPlane.second, ",", _vec2);
		t_vec vec1 = tl::convert_vec_full<t_real, t_real, std::vector, tl::ublas::vector>(_vec1);
		t_vec vec2 = tl::convert_vec_full<t_real, t_real, std::vector, tl::ublas::vector>(_vec2);
		if(vec1.size() != 3) vec1.resize(3);
		if(vec2.size() != 3) vec2.resize(3);


		// Scan positions
		std::vector<std::string> vecPos;
		std::vector<t_vec> vecAllHKL, vecAllPos;
		tl::get_tokens<std::string, std::string>(strAllPos, ";", vecPos);
		for(std::string& strPos : vecPos)
		{
			tl::trim(strPos);
			if(strPos == "") continue;

			std::vector<t_real> _vecHKL;
			tl::parse_tokens<t_real, std::string>(strPos, ",", _vecHKL);
			t_vec vecHKL = tl::convert_vec_full<t_real, t_real, std::vector, tl::ublas::vector>(_vecHKL);
			if(vecHKL.size() != 3) vecHKL.resize(3);

			vecAllPos.push_back(get_plane_coord<t_vec, t_real>(vec1, vec2, vecHKL));
			vecAllHKL.push_back(std::move(vecHKL));
		}


		// Create plot script
		std::ostringstream ostrPlot;
		make_plot<t_vec, t_real>(ostrPlot, vec1, vec2, vecBragg, vecAllHKL, vecAllPos, bFlip, b3D);
		editScript->setPlainText(ostrPlot.str().c_str());
	}
	catch(const std::exception& ex)
	{
		tl::log_err(ex.what());
	}
}


void ScanPosDlg::Plot()
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


void ScanPosDlg::accept()
{}


void ScanPosDlg::ButtonBoxClicked(QAbstractButton* pBtn)
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
			m_pSettings->setValue("scanpos/geo", saveGeometry());
		QDialog::accept();
	}
}


#include "moc_ScanPosDlg.cpp"
