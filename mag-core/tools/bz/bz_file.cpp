/**
 * brillouin zone tool
 * @author Tobias Weber <tweber@ill.fr>
 * @date Maz-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include "bz.h"
#include "bz_conf.h"

#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMessageBox>
#include <QtSvg/QSvgGenerator>

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;

#include <boost/scope_exit.hpp>

#include "libs/loadcif.h"
#include "tlibs2/libs/algos.h"

using namespace tl2_ops;


void BZDlg::NewFile()
{
	m_ignoreCalc = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		this_->CalcB(true);
	} BOOST_SCOPE_EXIT_END

	// clear old tables
	DelSymOpTabItem(-1);
	DelFormulaTabItem(-1);

	// set some defaults
	m_comboSG->setCurrentIndex(0);
	m_editA->setValue(5.);
	m_editB->setValue(5.);
	m_editC->setValue(5.);
	m_editAlpha->setValue(90.);
	m_editBeta->setValue(90.);
	m_editGamma->setValue(90.);

	m_cutX->setValue(1);
	m_cutY->setValue(0);
	m_cutZ->setValue(0);
	m_cutNX->setValue(0);
	m_cutNY->setValue(0);
	m_cutNZ->setValue(1);
	m_cutD->setValue(0);
	m_BZDrawOrder->setValue(4);
	m_BZCalcOrder->setValue(4);
}


bool BZDlg::Load(const QString& filename, bool use_stdin)
{
	m_ignoreCalc = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		this_->CalcB(true);
	} BOOST_SCOPE_EXIT_END

	try
	{
		BZConfig cfg = load_bz_config(filename.toStdString(), use_stdin);

		// clear old items
		DelSymOpTabItem(-1);

		// settings
		if(cfg.xtal_a)
			m_editA->setValue(*cfg.xtal_a);
		if(cfg.xtal_b)
			m_editB->setValue(*cfg.xtal_b);
		if(cfg.xtal_c)
			m_editC->setValue(*cfg.xtal_c);
		if(cfg.xtal_alpha)
			m_editAlpha->setValue(*cfg.xtal_alpha);
		if(cfg.xtal_beta)
			m_editBeta->setValue(*cfg.xtal_beta);
		if(cfg.xtal_gamma)
			m_editGamma->setValue(*cfg.xtal_gamma);
		if(cfg.order)
			m_BZCalcOrder->setValue(*cfg.order);
		if(cfg.cut_order)
			m_BZDrawOrder->setValue(*cfg.cut_order);
		if(cfg.cut_x)
			m_cutX->setValue(*cfg.cut_x);
		if(cfg.cut_y)
			m_cutY->setValue(*cfg.cut_y);
		if(cfg.cut_z)
			m_cutZ->setValue(*cfg.cut_z);
		if(cfg.cut_nx)
			m_cutNX->setValue(*cfg.cut_nx);
		if(cfg.cut_ny)
			m_cutNY->setValue(*cfg.cut_ny);
		if(cfg.cut_nz)
			m_cutNZ->setValue(*cfg.cut_nz);
		if(cfg.cut_d)
			m_cutD->setValue(*cfg.cut_d);
		if(cfg.sg_idx)
			m_comboSG->setCurrentIndex(*cfg.sg_idx);

		// symops
		for(const t_mat& symop : cfg.symops)
			AddSymOpTabItem(-1, symop);

		// formulas
		for(const std::string& formula : cfg.formulas)
			AddFormulaTabItem(-1, formula);
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Brillouin Zones", ex.what());
		m_ignoreCalc = false;
		return false;
	}


	return true;
}


bool BZDlg::Save(const QString& filename)
{
	pt::ptree node;

	// meta infos
	const char* user = std::getenv("USER");
	if(!user) user = "";

	node.put<std::string>("bz.meta.info", "bz_tool");
	node.put<std::string>("bz.meta.date", tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()));
	node.put<std::string>("bz.meta.user", user);
	node.put<std::string>("bz.meta.url", "https://github.com/ILLGrenoble/takin");
	node.put<std::string>("bz.meta.doi", "https://doi.org/10.5281/zenodo.4117437");
	node.put<std::string>("bz.meta.doi_tlibs", "https://doi.org/10.5281/zenodo.5717779");

	// lattice
	node.put<t_real>("bz.xtal.a", m_editA->value());
	node.put<t_real>("bz.xtal.b", m_editB->value());
	node.put<t_real>("bz.xtal.c", m_editC->value());
	node.put<t_real>("bz.xtal.alpha", m_editAlpha->value());
	node.put<t_real>("bz.xtal.beta", m_editBeta->value());
	node.put<t_real>("bz.xtal.gamma", m_editGamma->value());
	node.put<int>("bz.order", m_BZCalcOrder->value());
	node.put<int>("bz.cut.order", m_BZDrawOrder->value());
	node.put<t_real>("bz.cut.x", m_cutX->value());
	node.put<t_real>("bz.cut.y", m_cutY->value());
	node.put<t_real>("bz.cut.z", m_cutZ->value());
	node.put<t_real>("bz.cut.nx", m_cutNX->value());
	node.put<t_real>("bz.cut.ny", m_cutNY->value());
	node.put<t_real>("bz.cut.nz", m_cutNZ->value());
	node.put<t_real>("bz.cut.d", m_cutD->value());
	node.put<int>("bz.sg_idx", m_comboSG->currentIndex());

	// symop list
	pt::ptree symops;
	for(int row = 0; row < m_symops->rowCount(); ++row)
	{
		std::string opstr = m_symops->item(row, COL_OP)->text().toStdString();
		algo::replace_all(opstr, "\n", " ");

		pt::ptree itemNode;
		itemNode.put<std::string>("op", opstr);
		symops.insert(symops.end(), itemNode.begin(), itemNode.end());
	}
	node.add_child("bz.symops", symops);

	// formula list
	pt::ptree formulas;
	for(int row = 0; row < m_formulas->rowCount(); ++row)
	{
		std::string opstr = m_formulas->item(row, COL_FORMULA)->text().toStdString();
		algo::replace_all(opstr, "\n", " ");

		pt::ptree itemNode;
		itemNode.put<std::string>("expr", opstr);
		formulas.insert(formulas.end(), itemNode.begin(), itemNode.end());
	}
	node.add_child("bz.formulas", formulas);

	std::ofstream ofstr{filename.toStdString()};
	if(!ofstr)
	{
		QMessageBox::critical(this, "Brillouin Zones", "Cannot open file for writing.");
		return false;
	}
	ofstr.precision(g_prec);
	pt::write_xml(ofstr, node, pt::xml_writer_make_settings('\t', 1, std::string{"utf-8"}));

	return true;
}


void BZDlg::Load()
{
	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getOpenFileName(
		this, "Load File", dirLast, "XML Files (*.xml *.XML)");
	if(filename=="" || !QFile::exists(filename))
		return;

	if(Load(filename))
	{
		m_sett->setValue("dir", QFileInfo(filename).path());
		m_recent.AddRecentFile(filename);
	}
}


void BZDlg::Save()
{
	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save File", dirLast, "XML Files (*.xml *.XML)");
	if(filename=="")
		return;

	if(Save(filename))
	{
		m_sett->setValue("dir", QFileInfo(filename).path());
		m_recent.AddRecentFile(filename);
	}
}


/**
 * load a CIF
 */
void BZDlg::ImportCIF()
{
	m_ignoreCalc = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		this_->CalcB(true);
	} BOOST_SCOPE_EXIT_END

	try
	{
		QString dirLast = m_sett->value("dir_cif", "").toString();
		QString filename = QFileDialog::getOpenFileName(
			this, "Import CIF", dirLast, "CIF Files (*.cif *.CIF)");
		if(filename=="" || !QFile::exists(filename))
			return;
		m_sett->setValue("dir_cif", QFileInfo(filename).path());

		auto [errstr, atoms, generatedatoms, atomnames, lattice, symops] =
			load_cif<t_vec, t_mat>(filename.toStdString(), g_eps);
		if(errstr != "")
		{
			QMessageBox::critical(this, "CIF Importer", errstr.c_str());
			return;
		}


		// clear old symops
		DelSymOpTabItem(-1);

		// lattice
		m_editA->setValue(lattice.a);
		m_editB->setValue(lattice.b);
		m_editC->setValue(lattice.c);
		m_editAlpha->setValue(lattice.alpha);
		m_editBeta->setValue(lattice.beta);
		m_editGamma->setValue(lattice.gamma);


		// symops
		for(std::size_t opnum = 0; opnum < symops.size(); ++opnum)
			AddSymOpTabItem(-1, symops[opnum]);
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Brillouin Zones", ex.what());
	}
}


void BZDlg::SaveCutSVG()
{
	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save SVG File", dirLast, "SVG Files (*.svg *.SVG)");
	if(filename=="")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	QSvgGenerator svg;
	svg.setSize(QSize{800, 800});
	svg.setViewBox(QRect{0, 0, 800, 800});
	svg.setFileName(filename);
	svg.setTitle("Brillouin Zone Cut");
	svg.setDescription("Created with Takin (https://doi.org/10.5281/zenodo.4117437).");

	QPainter painter;
	painter.begin(&svg);
	m_bzscene->render(&painter);
	painter.end();
}
