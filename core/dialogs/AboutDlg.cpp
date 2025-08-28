/**
 * About Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date nov-2015
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

#include "AboutDlg.h"

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <qwt_global.h>
#include "tlibs/version.h"

#include "tlibs/string/string.h"
#include "libs/formfactors/formfact.h"
#include "libs/spacegroups/spacegroup.h"
#include "libs/globals.h"
#include "libs/version.h"
#include <sstream>

#ifndef NO_QHULL
	#include <Qhull.h>
#endif


#define __STR__(VAR) #VAR				// string of macro
#define __STR_DEREF__(VAR) __STR__(VAR)	// expand macro


AboutDlg::AboutDlg(QWidget* pParent, QSettings *pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	setupUi(this);
	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") && font.fromString(m_pSettings->value("main/font_gen", "").toString()))
			setFont(font);

		if(m_pSettings->contains("about/geo"))
			restoreGeometry(m_pSettings->value("about/geo").toByteArray());
	}

	labelVersion->setText("Version " TAKIN_VER ".");
	labelWritten->setText("Written by Tobias Weber <tweber@ill.fr>.");
	labelYears->setText("2014 - 2017 for Technische Universität München (TUM), Garching, Germany"
		";\n2017 - 2025 for Institut Laue-Langevin (ILL), Grenoble, France.");

	labelRepo->setText("Source repository: <a href=\"https://github.com/ILLGrenoble/takin\">https://github.com/ILLGrenoble/takin</a>.");

	labelDesc->setText("Papers about Takin can be found here, <b>please cite</b> the version closest to the one you have used:"
		"<ul>"
		"<li>Version 2.5: <a href=\"https://doi.org/10.1016/j.softx.2023.101471\">DOI: 10.1016/j.softx.2023.101471</a>,</li>"
		"<li>Version 2.0: <a href=\"https://doi.org/10.1016/j.softx.2021.100667\">DOI: 10.1016/j.softx.2021.100667</a>,</li>"
		"<li>Version 1.5: <a href=\"https://doi.org/10.1016/j.softx.2017.06.002\">DOI: 10.1016/j.softx.2017.06.002</a>,</li>"
		"<li>Version 1.0: <a href=\"https://doi.org/10.1016/j.softx.2016.06.002\">DOI: 10.1016/j.softx.2016.06.002</a>.</li>"
		"</ul>");
	labelDesc->setOpenExternalLinks(1);
	labelLicense->setOpenExternalLinks(1);

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
	strCC += ".";
	labelCC->setText(strCC.c_str());
	labelBuildDate->setText(QString("Build date: ") +
		QString(__DATE__) + ", " + QString(__TIME__) + ".");


	// -------------------------------------------------------------------------
	// Libraries

	std::ostringstream ostrLibs;
	ostrLibs << "<html><body>";
	ostrLibs << "<dl>";

	ostrLibs << "<dt>Uses Qt version " << QT_VERSION_STR << ".</dt>";
	ostrLibs << "<dd><a href=\"http://qt-project.org\">http://qt-project.org</a><br></dd>";

	ostrLibs << "<dt>Uses Qwt version " << QWT_VERSION_STR << ".</dt>";
	ostrLibs << "<dd><a href=\"http://qwt.sourceforge.net\">http://qwt.sourceforge.net</a><br></dd>";

	std::string strBoost = BOOST_LIB_VERSION;
	tl::find_all_and_replace<std::string>(strBoost, "_", ".");
	ostrLibs << "<dt>Uses Boost version " << strBoost << ".</dt>";
	ostrLibs << "<dd><a href=\"http://www.boost.org\">http://www.boost.org</a><br></dd>";

	ostrLibs << "<dt>Uses tlibs version " << TLIBS_VERSION << ".</dt>";
	ostrLibs << "<dd><a href=\"https://github.com/ILLGrenoble/takin/tree/master/tlibs\">https://github.com/ILLGrenoble/takin/tree/master/tlibs</a><br></dd>";

#ifndef NO_LAPACK
	ostrLibs << "<dt>Uses Lapack/e version 3.</dt>";
	ostrLibs << "<dd><a href=\"http://www.netlib.org/lapack\">http://www.netlib.org/lapack</a><br></dd>";
#endif

	ostrLibs << "<dt>Uses Minuit version 2.</dt>";
	ostrLibs << "<dd><a href=\"https://root.cern.ch\">https://root.cern.ch</a><br></dd>";

	ostrLibs << "<dt>Uses Clipper crystallography library.</dt>";
	ostrLibs << "<dd><a href=\"http://www.ysbl.york.ac.uk/~cowtan/clipper\">http://www.ysbl.york.ac.uk/~cowtan/clipper</a><br></dd>";

#ifdef HAS_COMPLEX_ERF
	ostrLibs << "<dt>Uses Faddeeva library.</dt>";
	ostrLibs << "<dd><a href=\"http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package\">http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package</a><br></dd>";
#endif

#ifndef NO_QHULL
	ostrLibs << "<dt>Uses QHull version " << qh_version << ". " << "</dt>";
	ostrLibs << "<dd><a href=\"http://www.qhull.org\">http://www.qhull.org</a><br></dd>";
#endif

	ostrLibs << "<dt>Uses resolution algorithms ported from Rescal version 5.</dt>";
	ostrLibs << "<dd><a href=\"http://www.ill.eu/en/html/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab\">http://www.ill.eu/en/html/instruments-support/computing-for-science/cs-software/all-software/matlab-ill/rescal-for-matlab</a><br></dd>";

	ostrLibs << "<dt>Uses Tango icons.</dt>";
	ostrLibs << "<dd><a href=\"http://tango.freedesktop.org\">http://tango.freedesktop.org</a><br></dd>";

	ostrLibs << "<dt>Uses DejaVu fonts.</dt>";
	ostrLibs << "<dd><a href=\"https://dejavu-fonts.github.io/\">https://dejavu-fonts.github.io</a><br></dd>";

	ostrLibs << "</dl>";
	//ostrLibs << "<p>See the LICENSES file in the Takin root directory.</p>";
	ostrLibs << "</body></html>";

	editLibs->setText(ostrLibs.str().c_str());



	// -------------------------------------------------------------------------
	// Features

	std::ostringstream ostrFeat;
	ostrFeat << "Feature flags: ";

#if !defined NO_NET
	ostrFeat << "Network ";
#endif

#if defined USE_PLUGINS
	ostrFeat << "Plugins ";
#endif

#if !defined NO_3D
	ostrFeat << "GL ";
#endif

#if !defined NO_LAPACK
	ostrFeat << "Lapack ";
#endif

#if defined HAS_COMPLEX_ERF
	ostrFeat << "Faddeeva ";
#endif

#if !defined NO_QHULL
	ostrFeat << "QHull ";
#endif

#if !defined NO_IOSTR
	ostrFeat << "Boost.Iostr ";
#endif

	ostrFeat << ".";
	labelFeatures->setText(ostrFeat.str().c_str());



	// -------------------------------------------------------------------------
	// Tables

	std::ostringstream ostrConst;
	ostrConst << "<html><body>";
	ostrConst << "<dl>";

	ostrConst << "<dt>Physical constants from Boost.Units.</dt>";
	ostrConst << "<dd><a href=\"http://www.boost.org/doc/libs/release/libs/units/\">http://www.boost.org/doc/libs/release/libs/units/</a><br></dd>";

	std::shared_ptr<const xtl::SpaceGroups<t_real_glob>> sgs = xtl::SpaceGroups<t_real_glob>::GetInstance();
	ostrConst << "<dt>" << sgs->get_sgsource(0) <<"</dt>";
	ostrConst << "<dd><a href=\"" << sgs->get_sgsource(1) << "\">" << sgs->get_sgsource(1) << "</a><br></dd>";

	std::shared_ptr<const xtl::FormfactList<t_real_glob>> ff = xtl::FormfactList<t_real_glob>::GetInstance();
	std::shared_ptr<const xtl::MagFormfactList<t_real_glob>> mff = xtl::MagFormfactList<t_real_glob>::GetInstance();
	std::shared_ptr<const xtl::ScatlenList<t_real_glob>> sl = xtl::ScatlenList<t_real_glob>::GetInstance();
	std::shared_ptr<const xtl::PeriodicSystem<t_real_glob>> pt = xtl::PeriodicSystem<t_real_glob>::GetInstance();

	if(g_bHasFormfacts)
	{
		ostrConst << "<dt>" << ff->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << ff->GetSourceUrl() << "\">" << ff->GetSourceUrl() << "</a><br></dd>";
	}
	if(g_bHasMagFormfacts)
	{
		ostrConst << "<dt>" << mff->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << mff->GetSourceUrl() << "\">" << mff->GetSourceUrl() << "</a><br></dd>";
	}
	if(g_bHasScatlens)
	{
		ostrConst << "<dt>" << sl->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << sl->GetSourceUrl() << "\">" << sl->GetSourceUrl() << "</a><br></dd>";
	}
	if(g_bHasElements)
	{
		ostrConst << "<dt>" << pt->GetSource() << "</dt>";
		ostrConst << "<dd><a href=\"" << pt->GetSourceUrl() << "\">" << pt->GetSourceUrl() << "</a><br></dd>";
	}

	ostrConst << "</dl>";
	ostrConst << "</body></html>";
	editTables->setText(ostrConst.str().c_str());



	std::string strLicensesFile = find_resource("LICENSES.txt");
	std::ifstream ifstrLicenses(strLicensesFile);
	std::string strLicenses;
	while(ifstrLicenses)
	{
		std::string strLic;
		std::getline(ifstrLicenses, strLic);
		strLicenses += strLic + "\n";
	}
	editAllLicenses->setPlainText(strLicenses.c_str());



	std::string strLitFile = find_resource("LITERATURE.txt");
	std::ifstream ifstrLit(strLitFile);
	std::string strLit;
	while(ifstrLit)
	{
		std::string _strLit;
		std::getline(ifstrLit, _strLit);
		strLit += _strLit + "\n";
	}
	editLiterature->setPlainText(strLit.c_str());
}


void AboutDlg::accept()
{
	if(m_pSettings)
		m_pSettings->setValue("about/geo", saveGeometry());

	QDialog::accept();
}


#include "moc_AboutDlg.cpp"
