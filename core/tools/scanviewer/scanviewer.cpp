/**
 * scan viewer
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2024
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

#include "scanviewer.h"

#include <QFileDialog>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QMessageBox>

#include <iostream>
#include <set>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <iterator>

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/filesystem.hpp>

#include "exporters.h"

#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/stat.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/file/file.h"
#include "tlibs/log/log.h"
#include "tlibs/helper/misc.h"

#ifdef USE_WIDE_STR
	#define T_STR std::wstring
#else
	#define T_STR std::string
#endif


using t_real = t_real_glob;
using t_vec = tl::ublas::vector<t_real>;

namespace fs = boost::filesystem;


ScanViewerDlg::ScanViewerDlg(QWidget* pParent, QSettings* core_settings)
	: QDialog(pParent, Qt::WindowTitleHint|Qt::WindowCloseButtonHint|Qt::WindowMinMaxButtonsHint),
		m_settings("takin", "scanviewer", this),
		m_core_settings{core_settings},
		m_vecExts({ ".dat", ".DAT", ".scn", ".SCN", ".ng0", ".NG0", ".log", ".LOG", ".nxs", ".NXS", ".hdf", ".HDF", "" }),
		m_pFitParamDlg(new FitParamDlg(this, &m_settings))
{
	this->setupUi(this);
	this->setFocusPolicy(Qt::StrongFocus);
	SetAbout();

	// also load takin's core settings if not explicitly given
	if(!m_core_settings)
		m_core_settings = new QSettings("takin", "core", this);

	QFont font;
	if(m_core_settings && m_core_settings->contains("main/font_gen") &&
		font.fromString(m_core_settings->value("main/font_gen", "").toString()))
		setFont(font);

	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 2);

	SetupPlotter(2);

	// -------------------------------------------------------------------------
	// property map stuff
	tableProps->setColumnCount(2);
	tableProps->setColumnWidth(0, 150);
	tableProps->setColumnWidth(1, 350);
	//tableProps->sortByColumn(0);

	tableProps->setHorizontalHeaderItem(0, new QTableWidgetItem("Property"));
	tableProps->setHorizontalHeaderItem(1, new QTableWidgetItem("Value"));

	tableProps->verticalHeader()->setVisible(false);
	tableProps->verticalHeader()->setDefaultSectionSize(tableProps->verticalHeader()->minimumSectionSize()+4);
	// -------------------------------------------------------------------------

	ScanViewerDlg *pThis = this;
	QObject::connect(comboPath, &QComboBox::editTextChanged, pThis, &ScanViewerDlg::ChangedPath);
	QObject::connect(listFiles, &QListWidget::itemSelectionChanged, pThis, &ScanViewerDlg::FileSelected);
	QObject::connect(editSearch, &QLineEdit::textEdited, pThis, &ScanViewerDlg::SearchProps);
	QObject::connect(btnBrowse, &QToolButton::clicked, pThis, static_cast<void (ScanViewerDlg::*)()>(&ScanViewerDlg::SelectDir));
	QObject::connect(btnRefresh, &QToolButton::clicked, pThis, static_cast<void (ScanViewerDlg::*)()>(&ScanViewerDlg::DirWasModified));
	for(QLineEdit* pEdit : {editPolVec1, editPolVec2, editPolCur1, editPolCur2})
		QObject::connect(pEdit, &QLineEdit::textEdited, pThis, &ScanViewerDlg::CalcPol);

	QObject::connect(btnParam, &QToolButton::clicked, pThis, &ScanViewerDlg::ShowFitParams);
	QObject::connect(btnGauss, &QToolButton::clicked, pThis, &ScanViewerDlg::FitGauss);
	QObject::connect(btnLorentz, &QToolButton::clicked, pThis, &ScanViewerDlg::FitLorentz);
	QObject::connect(btnVoigt, &QToolButton::clicked, pThis, &ScanViewerDlg::FitVoigt);
	QObject::connect(btnLine, &QToolButton::clicked, pThis, &ScanViewerDlg::FitLine);
	QObject::connect(btnParabola, &QToolButton::clicked, pThis, &ScanViewerDlg::FitParabola);
	QObject::connect(btnSine, &QToolButton::clicked, pThis, &ScanViewerDlg::FitSine);

	QObject::connect(comboX, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::XAxisSelected);
	QObject::connect(comboY, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::YAxisSelected);
	QObject::connect(comboMon, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::MonAxisSelected);
	QObject::connect(checkNorm, static_cast<void (QCheckBox::*)(int)>(&QCheckBox::stateChanged), pThis, &ScanViewerDlg::NormaliseStateChanged);
	QObject::connect(checkMerge, static_cast<void (QCheckBox::*)(int)>(&QCheckBox::stateChanged), pThis, &ScanViewerDlg::MergeStateChanged);
	QObject::connect(spinStart, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), pThis, &ScanViewerDlg::StartOrSkipChanged);
	QObject::connect(spinStop, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), pThis, &ScanViewerDlg::StartOrSkipChanged);
	QObject::connect(spinSkip, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), pThis, &ScanViewerDlg::StartOrSkipChanged);
	QObject::connect(tableProps, &QTableWidget::currentItemChanged, pThis, &ScanViewerDlg::PropSelected);
	QObject::connect(comboExport, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), pThis, &ScanViewerDlg::GenerateExternal);


	// fill recent paths combobox
	QStringList lstDirs = m_settings.value("recent_dirs").toStringList();
	for(const QString& strDir : lstDirs)
	{
		fs::path dir(strDir.toStdString());
		if(*tl::wstr_to_str(dir.native()).rbegin() != fs::path::preferred_separator)
			dir /= T_STR{} + fs::path::preferred_separator;

		if(HasRecentPath(strDir) < 0)
			comboPath->addItem(tl::wstr_to_str(dir.native()).c_str());
	}

	// last selected dir
	QString strDir = m_settings.value("last_dir", tl::wstr_to_str(fs::current_path().native()).c_str()).toString();

	int idx = HasRecentPath(strDir);
	if(idx < 0)
	{
		fs::path dir(strDir.toStdString());
		if(*tl::wstr_to_str(dir.native()).rbegin() != fs::path::preferred_separator)
			dir /= T_STR{} + fs::path::preferred_separator;

		comboPath->addItem(tl::wstr_to_str(dir.native()).c_str());
		idx = comboPath->findText(strDir);
	}

	comboPath->setCurrentIndex(idx);
	//comboPath->setEditText(strDir);


	if(m_settings.contains("norm"))
		checkNorm->setChecked(m_settings.value("norm").toBool());
	if(m_settings.contains("merge_mode"))
		checkMerge->setChecked(m_settings.value("merge_mode").toBool());
	if(m_settings.contains("pol/vec1"))
		editPolVec1->setText(m_settings.value("pol/vec1").toString());
	if(m_settings.contains("pol/vec2"))
		editPolVec2->setText(m_settings.value("pol/vec2").toString());
	if(m_settings.contains("pol/cur1"))
		editPolCur1->setText(m_settings.value("pol/cur1").toString());
	if(m_settings.contains("pol/cur2"))
		editPolCur2->setText(m_settings.value("pol/cur2").toString());

	m_bDoUpdate = true;
	ChangedPath();

#ifndef HAS_COMPLEX_ERF
	btnVoigt->setEnabled(false);
#endif

	if(m_settings.contains("geo"))
		restoreGeometry(m_settings.value("geo").toByteArray());
	if(m_settings.contains("splitter"))
		splitter->restoreState(m_settings.value("splitter").toByteArray());
}


ScanViewerDlg::~ScanViewerDlg()
{
	ClearPlot();
	tableProps->setRowCount(0);
	if(m_pFitParamDlg) { delete m_pFitParamDlg; m_pFitParamDlg = nullptr; }
}


void ScanViewerDlg::SetupPlotter(unsigned int numCurves)
{
	auto get_colour = [numCurves](unsigned int curve, bool is_curve) -> QColor
	{
		if(numCurves == 2)
		{
			// plotting single data set only
			if(is_curve)  // curve
				return QColor(0, 0, 0x99);
			else          // data points
				return QColor(0xff, 0, 0);
		}

		// plotting multiple data sets
		t_real lerpval = t_real(curve) / t_real(numCurves - 2);
		t_real r = tl::lerp<t_real, t_real>(1., 0., lerpval);
		t_real g = 0.;
		t_real b = tl::lerp<t_real, t_real>(0., 1., lerpval);

		return QColor(int(r * 255.), int(g * 255.), int(b * 255.));
	};

	// already set up with this number of curves?
	if(m_plotwrap && numCurves == m_plotwrap->GetNumCurves())
		return;

	QColor colorBck(240, 240, 240, 255);
	plot->setCanvasBackground(colorBck);

	bool show_legend = (numCurves > 2);  // only when plotting multiple scans
	m_plotwrap.reset(new QwtPlotWrapper(plot, numCurves, true, false, false, show_legend));

	for(unsigned int curve = 0; curve < numCurves; curve += 2)
	{
		// even indices are curves
		QPen penCurve;
		penCurve.setColor(get_colour(curve, true));
		penCurve.setWidth(2);
		m_plotwrap->GetCurve(curve)->setPen(penCurve);
		m_plotwrap->GetCurve(curve)->setStyle(QwtPlotCurve::CurveStyle::Lines);
		m_plotwrap->GetCurve(curve)->setTitle("Curve");
		m_plotwrap->GetCurve(curve)->setItemAttribute(QwtPlotCurve::Legend, false);

		// odd indices are scan points
		QPen penPoints;
		penPoints.setColor(get_colour(curve, false));
		penPoints.setWidth(4);
		m_plotwrap->GetCurve(curve + 1)->setPen(penPoints);
		m_plotwrap->GetCurve(curve + 1)->setStyle(QwtPlotCurve::CurveStyle::Dots);
		m_plotwrap->GetCurve(curve + 1)->setTitle("Data");
		m_plotwrap->GetCurve(curve + 1)->setItemAttribute(QwtPlotCurve::Legend, false);
	}

	m_plotwrap->GetPlot()->updateLegend();
}


void ScanViewerDlg::SetAbout()
{
	labelVersion->setText("Version " TAKIN_VER ".");
	labelWritten->setText("Written by Tobias Weber <tweber@ill.fr>.");
	labelYears->setText("Years: 2015 - 2025.");

	std::string strCC = "Built";
#ifdef BOOST_PLATFORM
	strCC += " for " + std::string(BOOST_PLATFORM);
#endif
	strCC += " using " + std::string(BOOST_COMPILER);
#ifdef __cplusplus
	strCC += " (standard: " + tl::var_to_str(__cplusplus) + ")";
#endif
#ifdef BOOST_STDLIB
	strCC += " with " + std::string(BOOST_STDLIB);
#endif
	strCC += " on " + std::string(__DATE__) + ", " + std::string(__TIME__);
	strCC += ".";
	labelCC->setText(strCC.c_str());
}


void ScanViewerDlg::closeEvent(QCloseEvent* pEvt)
{
	// save settings
	m_settings.setValue("norm", checkNorm->isChecked());
	m_settings.setValue("merge_mode", checkMerge->isChecked());
	m_settings.setValue("pol/vec1", editPolVec1->text());
	m_settings.setValue("pol/vec2", editPolVec2->text());
	m_settings.setValue("pol/cur1", editPolCur1->text());
	m_settings.setValue("pol/cur2", editPolCur2->text());
	m_settings.setValue("last_dir", QString(m_strCurDir.c_str()));
	m_settings.setValue("geo", saveGeometry());
	m_settings.setValue("splitter", splitter->saveState());

	// save recent directory list
	QStringList lstDirs;
	for(int iDir = 0; iDir < comboPath->count(); ++iDir)
	{
		QString qstrPath = comboPath->itemText(iDir);

		// clean up recent path list
		fs::path dir(qstrPath.toStdString());
		if(!fs::exists(dir) || !fs::is_directory(dir) || fs::is_empty(dir))
			continue;

		lstDirs << qstrPath;
	}
	m_settings.setValue("recent_dirs", lstDirs);

	QDialog::closeEvent(pEvt);
}


void ScanViewerDlg::keyPressEvent(QKeyEvent* pEvt)
{
	if(pEvt->key() == Qt::Key_R)
	{
		tl::log_debug("Refreshing file list...");
		ChangedPath();
	}

	QDialog::keyPressEvent(pEvt);
}


void ScanViewerDlg::ClearPlot()
{
	for(tl::FileInstrBase<t_real_glob>* instr : m_instrs)
	{
		if(!instr)
			continue;
		delete instr;
	}

	m_instrs.clear();

	m_vecX.clear();
	m_vecY.clear();
	m_vecYErr.clear();
	m_vecFitX.clear();
	m_vecFitY.clear();

	// set dummy data
	static const std::vector<t_real> dummy;
	for(unsigned int curve = 0; curve < m_plotwrap->GetNumCurves(); ++curve)
		set_qwt_data<t_real>()(*m_plotwrap, dummy, dummy, curve, false);

	m_strX = m_strY = m_strCmd = "";
	plot->setAxisTitle(QwtPlot::xBottom, "");
	plot->setAxisTitle(QwtPlot::yLeft, "");
	plot->setTitle("");

	auto edits = { editA, editB, editC,
		editAlpha, editBeta, editGamma,
		editPlaneX0, editPlaneX1, editPlaneX2,
		editPlaneY0, editPlaneY1, editPlaneY2,
		editTitle, editSample,
		editUser, editContact,
		editKfix, editTimestamp };
	for(auto* pEdit : edits)
		pEdit->setText("");

	comboX->clear();
	comboY->clear();
	comboMon->clear();
	textExportedFile->clear();
	textRawFile->clear();
	spinStart->setValue(0);
	spinStop->setValue(0);
	spinSkip->setValue(0);

	m_plotwrap->GetPlot()->updateLegend();
	m_plotwrap->GetPlot()->replot();
}


/**
 * check if given path is already in the "recent paths" list
 */
int ScanViewerDlg::HasRecentPath(const QString& strPath)
{
	fs::path dir(strPath.toStdString());
	if(*tl::wstr_to_str(dir.native()).rbegin() != fs::path::preferred_separator)
		dir /= T_STR{} + fs::path::preferred_separator;

	for(int iDir = 0; iDir < comboPath->count(); ++iDir)
	{
		QString strOtherPath = comboPath->itemText(iDir);
		fs::path dirOther(strOtherPath.toStdString());
		if(*tl::wstr_to_str(dirOther.native()).rbegin() != fs::path::preferred_separator)
			dirOther /= T_STR{} + fs::path::preferred_separator;

		if(dir == dirOther)
			return iDir;
	}

	return -1;
}


/**
 * select a new directory to browse
 */
void ScanViewerDlg::SelectDir(const QString& path)
{
	if(!tl::dir_exists(path.toStdString().c_str()))
		return;

	int idx = HasRecentPath(path);
	if(idx < 0)
	{
		fs::path dir(path.toStdString());
		if(*tl::wstr_to_str(dir.native()).rbegin() != fs::path::preferred_separator)
			dir /= T_STR{} + fs::path::preferred_separator;

		comboPath->addItem(tl::wstr_to_str(dir.native()).c_str());
		idx = comboPath->findText(tl::wstr_to_str(dir.native()).c_str());
	}

	comboPath->setCurrentIndex(idx);
	//comboPath->setEditText(path);
	ChangedPath();
}


/**
 * new scan directory selected
 */
void ScanViewerDlg::SelectDir()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(m_core_settings && !m_core_settings->value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strCurDir = (m_strCurDir=="" ? "~" : m_strCurDir.c_str());
	QString strDir = QFileDialog::getExistingDirectory(this, "Select directory",
		strCurDir, QFileDialog::ShowDirsOnly | fileopt);
	if(strDir != "")
		SelectDir(strDir);
}


void ScanViewerDlg::XAxisSelected(int) { PlotScan(); }
void ScanViewerDlg::YAxisSelected(int) { PlotScan(); }
void ScanViewerDlg::MonAxisSelected(int) { PlotScan(); }
void ScanViewerDlg::NormaliseStateChanged(int iState) { PlotScan(); }
void ScanViewerDlg::MergeStateChanged(int iState) { FileSelected(); }
void ScanViewerDlg::StartOrSkipChanged(int) { PlotScan(); }


/**
 * new file selected
 */
void ScanViewerDlg::FileSelected()
{
	// all selected items
	QList<QListWidgetItem*> lstSelected = listFiles->selectedItems();
	if(lstSelected.size() == 0)
		return;

	// clear previous files
	ClearPlot();

	// get first selected file
	m_strCurFile = lstSelected.first()->data(Qt::UserRole).toString().toStdString();

	// get the rest of the selected files
	std::vector<std::string> vecSelectedFiles, vecSelectedFilesRest;
	for(const QListWidgetItem *pLstItem : lstSelected)
	{
		if(!pLstItem)
			continue;

		std::string selectedFile = pLstItem->data(Qt::UserRole).toString().toStdString();
		vecSelectedFiles.push_back(selectedFile);

		// ignore first file
		if(pLstItem == lstSelected.first())
			continue;

		vecSelectedFilesRest.push_back(selectedFile);
	}

	ShowRawFiles(vecSelectedFiles);

	bool allow_incompatible_merging = false;
	if(m_core_settings)
		allow_incompatible_merging = m_core_settings->value("main/allow_scan_merging", 0).toBool();

	if(checkMerge->isChecked())
	{
		m_instrs.resize(1);
		m_instrs[0] = tl::FileInstrBase<t_real>::LoadInstr(m_strCurFile.c_str());
		if(!m_instrs[0])
			return;

		// merge with other selected files
		for(const std::string& strOtherFile : vecSelectedFilesRest)
		{
			std::unique_ptr<tl::FileInstrBase<t_real>> pToMerge(
				tl::FileInstrBase<t_real>::LoadInstr(strOtherFile.c_str()));
			if(!pToMerge)
				continue;

			m_instrs[0]->MergeWith(pToMerge.get(), allow_incompatible_merging);
		}
	}
	else
	{
		m_instrs.resize(vecSelectedFilesRest.size() + 1);

		std::size_t instr_idx = 0;
		m_instrs[instr_idx] = tl::FileInstrBase<t_real>::LoadInstr(m_strCurFile.c_str());
		if(m_instrs[instr_idx])
			++instr_idx;

		for(const std::string& strOtherFile : vecSelectedFilesRest)
		{
			m_instrs[instr_idx] = tl::FileInstrBase<t_real>::LoadInstr(strOtherFile.c_str());
			if(!m_instrs[instr_idx])
				continue;

			// check if files are compatible (i.e. would be mergeable)
			if(m_instrs[instr_idx] && !allow_incompatible_merging && instr_idx > 0
				&& !m_instrs[0]->IsCompatible(m_instrs[instr_idx]))
			{
				delete m_instrs[instr_idx];
				m_instrs[instr_idx] = nullptr;
				continue;
			}

			++instr_idx;
		}
	}

	std::vector<std::string> vecScanVars = m_instrs[0]->GetScannedVars();
	std::string strCntVar = m_instrs[0]->GetCountVar();
	std::string strMonVar = m_instrs[0]->GetMonVar();
	//tl::log_info("Count var: ", strCntVar, ", mon var: ", strMonVar);

	const std::wstring strPM = tl::get_spec_char_utf16("pm");  // +-

	m_bDoUpdate = false;
	int iIdxX = -1, iIdxY = -1, iIdxMon = -1, iCurIdx = 0;
	int iAlternateX = 0;
	const tl::FileInstrBase<t_real>::t_vecColNames& vecColNames = m_instrs[0]->GetColNames();
	for(const tl::FileInstrBase<t_real>::t_vecColNames::value_type& strCol : vecColNames)
	{
		const tl::FileInstrBase<t_real>::t_vecVals& vecCol = m_instrs[0]->GetCol(strCol);

		t_real dMean = tl::mean_value(vecCol);
		t_real dStd = tl::std_dev(vecCol);
		bool bStdZero = tl::float_equal(dStd, t_real(0), g_dEpsGfx);

		std::wstring _strCol = tl::str_to_wstr(strCol);
		_strCol += bStdZero ? L" (value: " : L" (mean: ";
		_strCol += tl::var_to_str<t_real, std::wstring>(dMean, g_iPrecGfx);
		if(!bStdZero)
		{
			_strCol += L" " + strPM + L" ";
			_strCol += tl::var_to_str<t_real, std::wstring>(dStd, g_iPrecGfx);
		}
		_strCol += L")";

		comboX->addItem(QString::fromWCharArray(_strCol.c_str()), QString(strCol.c_str()));
		comboY->addItem(QString::fromWCharArray(_strCol.c_str()), QString(strCol.c_str()));
		comboMon->addItem(QString::fromWCharArray(_strCol.c_str()), QString(strCol.c_str()));

		std::string strFirstScanVar = vecScanVars.size() ? tl::str_to_lower(vecScanVars[0]) : "";
		std::string strColLower = tl::str_to_lower(strCol);

		if(vecScanVars.size())
		{
			if(strFirstScanVar == strColLower)
				iIdxX = iCurIdx;
			else if(strFirstScanVar.substr(0, strCol.length()) == strColLower)
				iAlternateX = iCurIdx;
			// sometimes the scanned variable is named "QH", but the data column "H"
			else if(strFirstScanVar.substr(1) == strColLower)
				iAlternateX = iCurIdx;
		}

		// count and monitor variables
		if(tl::str_to_lower(strCntVar) == strColLower)
			iIdxY = iCurIdx;
		if(tl::str_to_lower(strMonVar) == strColLower)
			iIdxMon = iCurIdx;

		++iCurIdx;
	}

	if(iIdxX < 0 && iAlternateX >= 0)
		iIdxX = iAlternateX;
	comboX->setCurrentIndex(iIdxX);
	comboY->setCurrentIndex(iIdxY);
	comboMon->setCurrentIndex(iIdxMon);

	CalcPol();

	int iNumPol = static_cast<int>(m_instrs[0]->NumPolChannels()) - 1;
	if(iNumPol < 0)
		iNumPol = 0;
	spinSkip->setValue(iNumPol);

	m_bDoUpdate = true;
	ShowMetaData();
	ShowProps();
	PlotScan();
}


/**
 * highlights a scan property field
 */
void ScanViewerDlg::SearchProps(const QString& qstr)
{
	QList<QTableWidgetItem*> lstItems = tableProps->findItems(qstr, Qt::MatchContains);
	if(lstItems.size())
		tableProps->setCurrentItem(lstItems[0]);
}


void ScanViewerDlg::PlotScan()
{
	if(!m_bDoUpdate)
		return;

	SetupPlotter(m_instrs.size() * 2);  // curves and points

	bool bNormalise = checkNorm->isChecked();

	m_strX = comboX->itemData(comboX->currentIndex(), Qt::UserRole).toString().toStdString();
	m_strY = comboY->itemData(comboY->currentIndex(), Qt::UserRole).toString().toStdString();
	m_strMon = comboMon->itemData(comboMon->currentIndex(), Qt::UserRole).toString().toStdString();
	if(m_strX == "" || m_strY == "" || m_strMon == "")
		return;

	const int iStartIdx = spinStart->value();
	const int iEndSkip = spinStop->value();
	const int iSkipRows = spinSkip->value();

	m_vecX.resize(m_instrs.size());
	m_vecY.resize(m_instrs.size());
	m_vecYErr.resize(m_instrs.size());

	for(std::size_t instr_idx = 0; instr_idx < m_instrs.size(); ++instr_idx)
	{
		const tl::FileInstrBase<t_real_glob> *instr = m_instrs[instr_idx];
		if(!instr)
			continue;

		std::vector<t_real>& vecX = m_vecX[instr_idx];
		std::vector<t_real>& vecY = m_vecY[instr_idx];
		std::vector<t_real>& vecYErr = m_vecYErr[instr_idx];

		// get the data vectors
		vecX = instr->GetCol(m_strX.c_str());
		vecY = instr->GetCol(m_strY.c_str());
		std::vector<t_real> vecMon = instr->GetCol(m_strMon.c_str());

		// TODO: sort the data vectors
		//tl::sort_3(vecX.begin(), vecX.end(), vecY.begin(), vecMon.begin());

		bool bYIsACountVar = (m_strY == instr->GetCountVar() || m_strY == instr->GetMonVar());
		m_plotwrap->GetCurve(instr_idx*2 + 1)->SetShowErrors(bYIsACountVar);

		// see if there's a corresponding error column for the selected counter or monitor
		std::string ctr_err_col, mon_err_col;

		// get counter error if defined
		if(m_strY == instr->GetCountVar())
			ctr_err_col = instr->GetCountErr();
		else if(m_strY == instr->GetMonVar())
			ctr_err_col = instr->GetMonErr();

		if(ctr_err_col != "")
		{
			// use given error column
			vecYErr = instr->GetCol(ctr_err_col);
		}
		else
		{
			vecYErr.clear();
			vecYErr.reserve(vecY.size());

			// calculate error
			for(std::size_t iY = 0; iY < vecY.size(); ++iY)
			{
				t_real err = tl::float_equal(vecY[iY], t_real(0), g_dEps) ? t_real(1) : std::sqrt(std::abs(vecY[iY]));
				vecYErr.push_back(err);
			}
		}

		// get monitor error if defined
		std::vector<t_real> vecMonErr;
		if(m_strMon == instr->GetCountVar())
			mon_err_col = instr->GetCountErr();
		else if(m_strMon == instr->GetMonVar())
			mon_err_col = instr->GetMonErr();

		if(mon_err_col != "")
		{
			// use given error column
			vecMonErr = instr->GetCol(mon_err_col);
		}
		else
		{
			vecMonErr.reserve(vecY.size());

			// calculate error
			for(std::size_t iY = 0; iY < vecMon.size(); ++iY)
			{
				t_real err = tl::float_equal(vecMon[iY], t_real(0), g_dEps) ? t_real(1) : std::sqrt(std::abs(vecMon[iY]));
				vecMonErr.push_back(err);
			}
		}


		// remove points from start
		if(iStartIdx != 0)
		{
			if(std::size_t(iStartIdx) >= vecX.size())
				vecX.clear();
			else
				vecX.erase(vecX.begin(), vecX.begin()+iStartIdx);

			if(std::size_t(iStartIdx) >= vecY.size())
			{
				vecY.clear();
				vecMon.clear();
				vecYErr.clear();
				vecMonErr.clear();
			}
			else
			{
				vecY.erase(vecY.begin(), vecY.begin() + iStartIdx);
				vecMon.erase(vecMon.begin(), vecMon.begin() + iStartIdx);
				vecYErr.erase(vecYErr.begin(), vecYErr.begin() + iStartIdx);
				vecMonErr.erase(vecMonErr.begin(), vecMonErr.begin() + iStartIdx);
			}
		}

		// remove points from end
		if(iEndSkip != 0)
		{
			if(std::size_t(iEndSkip) >= vecX.size())
				vecX.clear();
			else
				vecX.erase(vecX.end() - iEndSkip, vecX.end());

			if(std::size_t(iEndSkip) >= vecY.size())
			{
				vecY.clear();
				vecMon.clear();
				vecYErr.clear();
				vecMonErr.clear();
			}
			else
			{
				vecY.erase(vecY.end() - iEndSkip, vecY.end());
				vecMon.erase(vecMon.end() - iEndSkip, vecMon.end());
				vecYErr.erase(vecYErr.end() - iEndSkip, vecYErr.end());
				vecMonErr.erase(vecMonErr.end() - iEndSkip, vecMonErr.end());
			}
		}

		// interleave rows
		if(iSkipRows != 0)
		{
			std::vector<t_real> vecXNew, vecYNew, vecMonNew, vecYErrNew, vecMonErrNew;

			for(std::size_t iRow = 0; iRow < std::min(vecX.size(), vecY.size()); ++iRow)
			{
				vecXNew.push_back(vecX[iRow]);
				vecYNew.push_back(vecY[iRow]);
				vecMonNew.push_back(vecMon[iRow]);
				vecYErrNew.push_back(vecYErr[iRow]);
				vecMonErrNew.push_back(vecMonErr[iRow]);

				iRow += iSkipRows;
			}

			vecX = std::move(vecXNew);
			vecY = std::move(vecYNew);
			vecMon = std::move(vecMonNew);
			vecYErr = std::move(vecYErrNew);
			vecMonErr = std::move(vecMonErrNew);
		}


		// errors
		if(vecMon.size() != vecY.size() || vecMonErr.size() != vecYErr.size())
		{
			bNormalise = false;
			tl::log_err("Counter and monitor data count do not match, cannot normalise.");
		}

		// normalise to monitor?
		if(bNormalise)
		{
			for(std::size_t iY = 0; iY < vecY.size(); ++iY)
			{
				if(tl::float_equal(vecMon[iY], t_real(0), g_dEps))
				{
					tl::log_warn("Monitor counter is zero for point ", iY + 1, ".");

					vecY[iY] = 0.;
					vecYErr[iY] = 1.;
				}
				else
				{
					std::tie(vecY[iY], vecYErr[iY]) = norm_cnts_to_mon(
						vecY[iY], vecYErr[iY], vecMon[iY], vecMonErr[iY]);
				}
			}
		}


		// show fit (for fist scan file)
		if(m_vecFitX.size() && instr_idx == 0)
			set_qwt_data<t_real>()(*m_plotwrap, m_vecFitX, m_vecFitY, instr_idx*2, false);
		else
			set_qwt_data<t_real>()(*m_plotwrap, vecX, vecY, instr_idx*2, false);
		set_qwt_data<t_real>()(*m_plotwrap, vecX, vecY, instr_idx*2 + 1, false, &vecYErr);

		// legend
		m_plotwrap->GetCurve(instr_idx*2 + 1)->setTitle(instr->GetScanNumber().c_str());
		m_plotwrap->GetCurve(instr_idx*2 + 0)->setItemAttribute(QwtPlotCurve::Legend, false);
		m_plotwrap->GetCurve(instr_idx*2 + 1)->setItemAttribute(QwtPlotCurve::Legend, m_instrs.size() > 1);
	}

	// labels
	QString strY = m_strY.c_str();
	if(bNormalise)
	{
		strY += " / ";
		strY += m_strMon.c_str();
	}
	plot->setAxisTitle(QwtPlot::xBottom, m_strX.c_str());
	plot->setAxisTitle(QwtPlot::yLeft, strY);
	plot->setTitle(m_strCmd.c_str());

	// replot
	set_zoomer_base(m_plotwrap->GetZoomer(), m_vecX, m_vecY, false, m_plotwrap.get(), true, &m_vecYErr);
	m_plotwrap->GetPlot()->updateLegend();
	m_plotwrap->GetPlot()->replot();

	GenerateExternal(comboExport->currentIndex());
}


/**
 * shown metadata from the first selected scan file (if multiple are selected)
 */
void ScanViewerDlg::ShowMetaData()
{
	if(!m_instrs.size())
		return;

	const tl::FileInstrBase<t_real_glob> *instr = m_instrs[0];
	if(!instr || !m_bDoUpdate)
		return;

	std::array<t_real, 3> arrLatt = instr->GetSampleLattice();
	std::array<t_real, 3> arrAng = instr->GetSampleAngles();
	std::array<t_real, 3> arrPlaneX = instr->GetScatterPlane0();
	std::array<t_real, 3> arrPlaneY = instr->GetScatterPlane1();

	// scattering plane normal
	t_vec plane_1 = tl::make_vec<t_vec>({ arrPlaneX[0],  arrPlaneX[1], arrPlaneX[2] });
	t_vec plane_2 = tl::make_vec<t_vec>({ arrPlaneY[0],  arrPlaneY[1], arrPlaneY[2] });
	t_vec plane_n = tl::cross_3(plane_1, plane_2);

	editA->setText(tl::var_to_str(arrLatt[0], g_iPrec).c_str());
	editB->setText(tl::var_to_str(arrLatt[1], g_iPrec).c_str());
	editC->setText(tl::var_to_str(arrLatt[2], g_iPrec).c_str());
	editAlpha->setText(tl::var_to_str(tl::r2d(arrAng[0]), g_iPrec).c_str());
	editBeta->setText(tl::var_to_str(tl::r2d(arrAng[1]), g_iPrec).c_str());
	editGamma->setText(tl::var_to_str(tl::r2d(arrAng[2]), g_iPrec).c_str());

	editPlaneX0->setText(tl::var_to_str(arrPlaneX[0], g_iPrec).c_str());
	editPlaneX1->setText(tl::var_to_str(arrPlaneX[1], g_iPrec).c_str());
	editPlaneX2->setText(tl::var_to_str(arrPlaneX[2], g_iPrec).c_str());
	editPlaneY0->setText(tl::var_to_str(arrPlaneY[0], g_iPrec).c_str());
	editPlaneY1->setText(tl::var_to_str(arrPlaneY[1], g_iPrec).c_str());
	editPlaneY2->setText(tl::var_to_str(arrPlaneY[2], g_iPrec).c_str());
	editPlaneZ0->setText(tl::var_to_str(plane_n[0], g_iPrec).c_str());
	editPlaneZ1->setText(tl::var_to_str(plane_n[1], g_iPrec).c_str());
	editPlaneZ2->setText(tl::var_to_str(plane_n[2], g_iPrec).c_str());

	labelKfix->setText(instr->IsKiFixed()
		? QString::fromWCharArray(L"ki (1/\x212b):")
		: QString::fromWCharArray(L"kf (1/\x212b):"));
	editKfix->setText(tl::var_to_str(instr->GetKFix()).c_str());

	editTitle->setText(instr->GetTitle().c_str());
	editProposal->setText(instr->GetProposal().c_str());
	editSample->setText(instr->GetSampleName().c_str());
	editUser->setText(instr->GetUser().c_str());
	editContact->setText(instr->GetLocalContact().c_str());
	editInstrument->setText(instr->GetInstrument().c_str());
	editTimestamp->setText(instr->GetTimestamp().c_str());

	m_strCmd = instr->GetScanCommand();
}


/**
 * convert to external plotter format
 * TODO: include all plot curves, not only the first
 */
void ScanViewerDlg::GenerateExternal(int iLang)
{
	textExportedFile->clear();
	if(!m_vecX.size() || !m_vecY.size() || !m_vecYErr.size())
		return;

	const std::vector<t_real>& vecX = m_vecX[0];
	const std::vector<t_real>& vecY = m_vecY[0];
	const std::vector<t_real>& vecYErr = m_vecYErr[0];

	std::string strSrc;

	if(iLang == 0)	// gnuplot
	{
		strSrc = export_scan_to_gnuplot<std::vector<t_real>>(
			vecX, vecY, vecYErr, m_strX, m_strY, m_strCmd, m_strCurFile);
	}
	else if(iLang == 1)	// root
	{
		strSrc = export_scan_to_root<std::vector<t_real>>(
			vecX, vecY, vecYErr, m_strX, m_strY, m_strCmd, m_strCurFile);
	}
	else if(iLang == 2)	// python
	{
		strSrc = export_scan_to_python<std::vector<t_real>>(
			vecX, vecY, vecYErr, m_strX, m_strY, m_strCmd, m_strCurFile);
	}
	else if(iLang == 3) // julia
	{
		strSrc = export_scan_to_julia<std::vector<t_real>>(
			vecX, vecY, vecYErr, m_strX, m_strY, m_strCmd, m_strCurFile);
	}
	else if(iLang == 4) // hermelin
	{
		strSrc = export_scan_to_hermelin<std::vector<t_real>>(
			vecX, vecY, vecYErr, m_strX, m_strY, m_strCmd, m_strCurFile);
	}
	else
	{
		tl::log_err("Unknown external language.");
	}

	textExportedFile->setText(strSrc.c_str());
}


/**
 * show raw scan files
 */
void ScanViewerDlg::ShowRawFiles(const std::vector<std::string>& files)
{
	QString rawFiles;

	for(const std::string& file : files)
	{
		std::size_t size = tl::get_file_size(file);
		std::ifstream ifstr(file);
		if(!ifstr)
			continue;

		auto ch_ptr = std::unique_ptr<char[]>(new char[size + 1]);
		ch_ptr[size] = 0;
		ifstr.read(ch_ptr.get(), size);

		// check if the file is of non-binary type
		bool is_printable = std::all_of(ch_ptr.get(), ch_ptr.get() + size, [](char ch) -> bool
		{
			bool is_bin = std::iscntrl(ch) != 0 && std::isspace(ch) == 0;
			//if(is_bin)
			//	std::cerr << "Non-printable character: 0x" << std::hex << int((unsigned char)ch) << std::endl;
			return !is_bin;
		});

		if(is_printable)
		{
			//rawFiles += ("# File: " + file + "\n").c_str();
			rawFiles += ch_ptr.get();
			rawFiles += "\n";
		}
		else
		{
			std::string file_ext = tl::get_fileext(file);
			if(file_ext == "nxs" || file_ext == "hdf")
			{
				tl::FileInstrBase<t_real> *instr = tl::FileInstrBase<t_real>::LoadInstr(file.c_str());
				if(!instr)
				{
					tl::log_err("Cannot load binary file \"", file, "\".");
					continue;
				}

				std::ostringstream ostr;
				ostr.precision(g_iPrec);

				if(!instr->Save(ostr))
				{
					tl::log_err("Cannot convert binary file \"", file, "\".");
					continue;
				}

				rawFiles += ostr.str().c_str();
			}
			else
			{
				rawFiles += "<unknown binary file>";
			}
		}
	}

	textRawFile->setText(rawFiles);
}


/**
 * save selected property key for later
 */
void ScanViewerDlg::PropSelected(QTableWidgetItem *pItem, QTableWidgetItem *pItemPrev)
{
	if(!pItem)
		m_strSelectedKey = "";

	for(int iItem = 0; iItem < tableProps->rowCount(); ++iItem)
	{
		const QTableWidgetItem *pKey = tableProps->item(iItem, 0);
		const QTableWidgetItem *pVal = tableProps->item(iItem, 1);

		if(pKey == pItem || pVal == pItem)
		{
			m_strSelectedKey = pKey->text().toStdString();
			break;
		}
	}
}


/**
 * save selected property key for later
 */
void ScanViewerDlg::ShowProps()
{
	if(!m_instrs.size())
		return;

	const tl::FileInstrBase<t_real_glob> *instr = m_instrs[0];
	if(!instr || !m_bDoUpdate)
		return;

	const tl::FileInstrBase<t_real>::t_mapParams& params = instr->GetAllParams();
	tableProps->setRowCount(params.size());

	const bool bSort = tableProps->isSortingEnabled();
	tableProps->setSortingEnabled(0);
	unsigned int iItem = 0;
	for(const tl::FileInstrBase<t_real>::t_mapParams::value_type& pair : params)
	{
		QTableWidgetItem *pItemKey = tableProps->item(iItem, 0);
		if(!pItemKey)
		{
			pItemKey = new QTableWidgetItem();
			tableProps->setItem(iItem, 0, pItemKey);
		}

		QTableWidgetItem* pItemVal = tableProps->item(iItem, 1);
		if(!pItemVal)
		{
			pItemVal = new QTableWidgetItem();
			tableProps->setItem(iItem, 1, pItemVal);
		}

		pItemKey->setText(pair.first.c_str());
		pItemVal->setText(pair.second.c_str());

		++iItem;
	}

	tableProps->setSortingEnabled(bSort);
	//tableProps->sortItems(0, Qt::AscendingOrder);


	// retain previous selection
	bool bHasSelection = false;
	for(int iItem = 0; iItem < tableProps->rowCount(); ++iItem)
	{
		const QTableWidgetItem *pItem = tableProps->item(iItem, 0);
		if(!pItem) continue;

		if(pItem->text().toStdString() == m_strSelectedKey)
		{
			tableProps->selectRow(iItem);
			bHasSelection = true;
			break;
		}
	}

	if(!bHasSelection)
		tableProps->selectRow(0);
}


/**
 * entered new directory
 */
void ScanViewerDlg::ChangedPath()
{
	listFiles->clear();
	ClearPlot();
	tableProps->setRowCount(0);

	m_strCurDir = "";
	std::vector<std::string> dirs;
	tl::get_tokens<std::string, std::string, std::vector<std::string>>(comboPath->currentText().toStdString(), ";", dirs);

	// watch directory for changes
	m_pWatcher.reset(new QFileSystemWatcher(this));
	QObject::connect(m_pWatcher.get(), &QFileSystemWatcher::directoryChanged,
		 this, &ScanViewerDlg::DirWasModified);

	for(const std::string& curDir : dirs)
	{
		if(!tl::dir_exists(curDir.c_str()))
			continue;

		if(m_strCurDir != "")
			m_strCurDir += ";";
		m_strCurDir += tl::wstr_to_str(fs::path(curDir).native());
		tl::trim(m_strCurDir);
		std::size_t len = m_strCurDir.length();
		if(len > 0 && *(m_strCurDir.begin() + len - 1) != fs::path::preferred_separator)
			m_strCurDir += fs::path::preferred_separator;

		// watch directory for changes
		m_pWatcher->addPath(curDir.c_str());
	}

	UpdateFileList();
}


/**
 * the current directory has been modified externally
 */
void ScanViewerDlg::DirWasModified()
{
	// get currently selected item
	QString strTxt;
	const QListWidgetItem *pCur = listFiles->currentItem();
	if(pCur)
		strTxt = pCur->text();

	UpdateFileList();

	// re-select previously selected item
	if(pCur)
	{
		QList<QListWidgetItem*> lstItems = listFiles->findItems(strTxt, Qt::MatchExactly);
		if(lstItems.size())
			listFiles->setCurrentItem(*lstItems.begin(), QItemSelectionModel::SelectCurrent);
	}
}


/**
 * re-populate file list
 */
void ScanViewerDlg::UpdateFileList()
{
	listFiles->clear();

	std::vector<std::string> dirs;
	tl::get_tokens<std::string, std::string, std::vector<std::string>>(m_strCurDir, ";", dirs);

	for(const std::string& curDir : dirs)
	{
		try
		{
			fs::path dir(curDir);
			fs::directory_iterator dir_begin(dir), dir_end;

			std::set<fs::path> lst;
			std::copy_if(dir_begin, dir_end, std::insert_iterator<decltype(lst)>(lst, lst.end()),
				[this](const fs::path& p) -> bool
				{
					// ignore non-existing files and directories
					if(!tl::file_exists(p.string().c_str()))
						return false;

					std::string strExt = tl::wstr_to_str(p.extension().native());
					if(strExt == ".bz2" || strExt == ".gz" || strExt == ".z")
						strExt = "." + tl::wstr_to_str(tl::get_fileext2(p.filename().native()));

					// allow everything if no extensions are defined
					if(this->m_vecExts.size() == 0)
						return true;

					// see if extension is in list
					return std::find(this->m_vecExts.begin(), this->m_vecExts.end(),
						strExt) != this->m_vecExts.end();
				});

			for(const fs::path& d : lst)
			{
				QListWidgetItem *item = new QListWidgetItem(tl::wstr_to_str(d.filename().native()).c_str());
				item->setData(Qt::UserRole, QString(tl::wstr_to_str(d.string()).c_str()));
				listFiles->addItem(item);
			}
		}
		catch(const std::exception& ex)
		{}
	}
}


#include "moc_scanviewer.cpp"
