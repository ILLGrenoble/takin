/**
 * scan viewer
 * @author Tobias Weber <tweber@ill.fr>
 * @date mar-2015 - 2025
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

#include <QTableWidget>
#include <QTableWidgetItem>

#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/filesystem.hpp>

#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/stat.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/string/spec_char.h"
#include "tlibs/file/file.h"
#include "tlibs/log/log.h"
#include "tlibs/helper/misc.h"

#include "libs/version.h"

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

	font.setStyleHint(QFont::Monospace);
	textRawFile->setFont(font);

	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 2);

	SetupColumnAliases();
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
	tableProps->verticalHeader()->setDefaultSectionSize(tableProps->verticalHeader()->minimumSectionSize() + 4);
	// -------------------------------------------------------------------------

	ScanViewerDlg *pThis = this;
	QObject::connect(comboPath, &QComboBox::editTextChanged, pThis, &ScanViewerDlg::ChangedPath);
	QObject::connect(listFiles, &QListWidget::itemSelectionChanged, pThis, &ScanViewerDlg::FileSelected);
	QObject::connect(editSearch, &QLineEdit::textEdited, pThis, &ScanViewerDlg::SearchProps);
	QObject::connect(editFileExts, &QLineEdit::textEdited, pThis, static_cast<void (ScanViewerDlg::*)()>(&ScanViewerDlg::DirWasModified));
	QObject::connect(btnResetExts, &QToolButton::clicked, pThis, static_cast<void (ScanViewerDlg::*)()>(&ScanViewerDlg::ResetFileExtensions));
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
#if QT_VERSION >= QT_VERSION_CHECK(6, 7, 0)
	QObject::connect(checkNorm, static_cast<void (QCheckBox::*)(Qt::CheckState)>(&QCheckBox::checkStateChanged), pThis, &ScanViewerDlg::NormaliseStateChanged);
	QObject::connect(checkMerge, static_cast<void (QCheckBox::*)(Qt::CheckState)>(&QCheckBox::checkStateChanged), pThis, &ScanViewerDlg::MergeStateChanged);
#else
	QObject::connect(checkNorm, static_cast<void (QCheckBox::*)(int)>(&QCheckBox::stateChanged), pThis, &ScanViewerDlg::NormaliseStateChanged);
	QObject::connect(checkMerge, static_cast<void (QCheckBox::*)(int)>(&QCheckBox::stateChanged), pThis, &ScanViewerDlg::MergeStateChanged);
#endif
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
	if(m_settings.contains("file_exts"))
		editFileExts->setText(m_settings.value("file_exts").toString());

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
	m_settings.setValue("file_exts", editFileExts->text());
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


void ScanViewerDlg::XAxisSelected(int) { PlotScan(); }
void ScanViewerDlg::YAxisSelected(int) { PlotScan(); }
void ScanViewerDlg::MonAxisSelected(int) { PlotScan(); }
void ScanViewerDlg::NormaliseStateChanged(int iState) { PlotScan(); }
void ScanViewerDlg::MergeStateChanged(int iState) { FileSelected(); }
void ScanViewerDlg::StartOrSkipChanged(int) { PlotScan(); }


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
 * highlights a scan property field
 */
void ScanViewerDlg::SearchProps(const QString& qstr)
{
	QList<QTableWidgetItem*> lstItems = tableProps->findItems(qstr, Qt::MatchContains);
	if(lstItems.size())
		tableProps->setCurrentItem(lstItems[0]);
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


#include "moc_scanviewer.cpp"
