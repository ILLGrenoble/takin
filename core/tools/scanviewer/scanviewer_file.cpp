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

#include <QFileDialog>

#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <algorithm>
#include <iterator>

#include <boost/filesystem.hpp>

#include "exporters.h"

#include "tlibs/math/math.h"
#include "tlibs/math/linalg.h"
#include "tlibs/math/stat.h"
#include "tlibs/phys/neutrons.h"
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


/**
 * get a data column from a scan file
 */
const std::vector<t_real>& ScanViewerDlg::GetCol(const tl::FileInstrBase<t_real_glob>* instr,
	const std::string& colname) const
{
	static const std::vector<t_real> nullcol;
	if(!instr)
		return nullcol;

	const std::size_t invalid_idx = instr->GetColNames().size();

	// try given column name
	std::size_t idx = 0;
	const std::vector<t_real>& col = instr->GetCol(colname, &idx);
	if(idx != invalid_idx)
		return col;

	// try column name aliases
	auto iter = m_col_aliases.find(colname);
	if(iter != m_col_aliases.end())
	{
		for(const std::string& alias : iter->second)
		{
			const std::vector<t_real>& aliascol = instr->GetCol(colname, &idx);
			if(idx != invalid_idx)
				return aliascol;
		}
	}

	return nullcol;
}



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
		std::string file_ext = tl::get_fileext(file);

		// see if it's an nxs file; if so, convert it to text
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
			rawFiles += "\n";
			continue;
		}

		// else if it's a text file, display it in raw format
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
			rawFiles += "<unknown binary file>\n";
		}
	}

	textRawFile->setText(rawFiles);
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


/**
 * load configuration file with alternate names for data columns
 */
void ScanViewerDlg::SetupColumnAliases()
{
	m_col_aliases.clear();

	std::string aliasfile = find_resource("column_aliases.cfg", false);
	if(aliasfile == "")
		return;

	std::ifstream ifstr(aliasfile);
	if(!ifstr)
	{
		tl::log_err("Could not open column alias configuration file \"", aliasfile, "\".");
		return;
	}

	std::size_t num_aliases = 0;
	while(!ifstr.eof())
	{
		std::string line;
		std::getline(ifstr, line);

		std::pair<std::string, std::string> keyval = tl::split_first<std::string>(line, "=:", true);
		if(keyval.first == "" || keyval.second == "")
			continue;

		std::vector<std::string> vals;
		tl::get_tokens<std::string, std::string, std::vector<std::string>>(keyval.second, ",;", vals);
		if(vals.size() == 0)
			continue;

		m_col_aliases.emplace(std::make_pair(keyval.first, vals));
		++num_aliases;
	}

	tl::log_info("Found ", num_aliases, " column alias configurations in file \"", aliasfile, "\".");
}
