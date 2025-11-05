/**
 * xmonteconvo
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2015
 * @copyright GPLv2
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

#include <clocale>
#include <QLocale>
#include <QApplication>
#include <QDir>
#include <QSettings>
#include <QMessageBox>

#include "tlibs/math/rand.h"
#include "tlibs/file/prop.h"
#include "libs/version.h"
#include "libs/globals.h"
#include "ConvoDlg.h"


int main(int argc, char** argv)
{
	tl::init_rand();

	QApplication app(argc, argv);
	app.setApplicationName("Takin/Monteconvo");
	app.setApplicationVersion(TAKIN_VER);

	// set up resource paths
	std::string strHome = QDir::homePath().toStdString() + "/.takin";
	std::string strApp = app.applicationDirPath().toStdString();
	add_resource_path(strHome, 0);
	add_resource_path(strApp);

	// set locale
	std::setlocale(LC_ALL, "C");
	std::locale::global(std::locale::classic());
	QLocale::setDefault(QLocale::English);

	QSettings m_settings("takin", "core");
	ConvoDlg dlg(nullptr, &m_settings);
	dlg.setWindowFlags(Qt::Window);

	// find arguments
	int iNextValidOption = -1;
	for(int iArg = 1; iArg < argc; ++iArg)
	{
		std::string strArg = argv[iArg];

		// max. number of threads
		if(strArg == "-t" && iArg+1 < argc)
		{
			std::string strVal = argv[iArg+1];

			g_iMaxThreads = tl::str_to_var<int>(strVal);
			++iArg;
			iNextValidOption = iArg+1;
		}
	}

	// load a given file
	if(argc > 1 && iNextValidOption < argc)
	{
		tl::log_info("Loading \"", argv[iNextValidOption], "\".");

		const std::string strXmlRoot("taz/");
		tl::Prop<std::string> xml;
		if(xml.Load(argv[iNextValidOption], tl::PropType::XML))
			dlg.Load(xml, strXmlRoot);
		else
			QMessageBox::critical(nullptr, "Error", "Could not load convolution file.");
	}

	// run
	dlg.show();
	return app.exec();
}
