/**
 * magnetic dynamics -- loading and saving of the configuration
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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

#include <boost/scope_exit.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include "magdyn.h"

#include <QtCore/QString>

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "tlibs2/libs/str.h"


// precision
extern int g_prec;



void MagDynDlg::Clear()
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		//this_->SyncToKernel();

		if(this_->m_structplot_dlg)
			this_->m_structplot_dlg->Sync();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// clear old tables
	DelTabItem(m_sitestab, -1);
	DelTabItem(m_termstab, -1);
	DelTabItem(m_varstab, -1);
	DelTabItem(m_fieldstab, -1);
	DelTabItem(m_coordinatestab, -1);

	ClearDispersion(true);
	m_hamiltonian->clear();
	m_dyn.Clear();

	SetCurrentFile("");

	// reset some defaults
	m_comboSG->setCurrentIndex(0);

	m_ordering[0]->setValue(0.);
	m_ordering[1]->setValue(0.);
	m_ordering[2]->setValue(0.);

	m_normaxis[0]->setValue(1.);
	m_normaxis[1]->setValue(0.);
	m_normaxis[2]->setValue(0.);

	m_weight_scale->setValue(1.);
	m_weight_min->setValue(0.);
	m_weight_max->setValue(9999.);

	if(m_groundstate_dlg)
		m_groundstate_dlg->SyncFromKernel();
	if(m_notes_dlg)
		m_notes_dlg->ClearNotes();

	// reset some options
	for(int i = 0; i < 3; ++i)
	{
		m_hamiltonian_comp[i]->setChecked(true);
		m_plot_channel[i]->setChecked(true);
	}

	m_statusFixed->setText("Ready.");
	m_status->setText("");
}



/**
 * set the currently open file and the corresponding window title
 */
void MagDynDlg::SetCurrentFile(const QString& filename)
{
	m_recent.SetCurFile(filename);

	QString title = "Magnetic Dynamics";
	if(filename != "")
		title += " - " + filename;
	setWindowTitle(title);
}



/**
 * set the currently open file and its directory
 */
void MagDynDlg::SetCurrentFileAndDir(const QString& filename)
{
	if(filename == "")
		return;

	m_sett->setValue("dir", QFileInfo(filename).path());
	m_recent.AddRecentFile(filename);
	SetCurrentFile(filename);
}



// --------------------------------------------------------------------------------
/**
 * load magnetic structure configuration
 */
void MagDynDlg::Load()
{
	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getOpenFileName(
		this, "Load File", dirLast, "Magnetic Dynamics Files (*.magdyn *.xml)");
	if(filename == "" || !QFile::exists(filename))
		return;

	Clear();

	if(Load(filename))
		SetCurrentFileAndDir(filename);
}



/**
 * load magnetic structure configuration
 */
bool MagDynDlg::Load(const QString& filename, bool calc_dynamics)
{
	try
	{
		BOOST_SCOPE_EXIT(this_, calc_dynamics)
		{
			this_->m_ignoreCalc = false;
			if(this_->m_autocalc->isChecked())
			{
				if(calc_dynamics)
					this_->CalcAll();
				else
					this_->SyncToKernel();
			}
		} BOOST_SCOPE_EXIT_END
		m_ignoreCalc = true;

		// properties tree
		pt::ptree node;

		// load from file
		std::ifstream ifstr{filename.toStdString()};
		pt::read_xml(ifstr, node);

		// check signature
		if(auto optInfo = node.get_optional<std::string>("magdyn.meta.info");
			!optInfo || !(*optInfo==std::string{"magdyn_tool"}))
		{
			QMessageBox::critical(this, "Magnetic Dynamics", "Unrecognised file format.");
			return false;
		}

		// read in comment
		if(auto optNotes = node.get_optional<std::string>("magdyn.meta.notes"); optNotes)
			m_notes_dlg->SetNotes(*optNotes);

		const pt::ptree &magdyn = node.get_child("magdyn");

		// settings
		if(auto optVal = magdyn.get_optional<t_real>("config.h_start"))
			m_Q_start[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.k_start"))
			m_Q_start[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.l_start"))
			m_Q_start[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.h_end"))
			m_Q_end[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.k_end"))
			m_Q_end[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.l_end"))
			m_Q_end[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.h"))
			m_Q[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.k"))
			m_Q[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.l"))
			m_Q[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_size>("config.num_Q_points"))
			m_num_points->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.weight_scale"))
			m_weight_scale->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.weight_min"))
			m_weight_min->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.weight_max"))
			m_weight_max->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.plot_channels"))
			m_plot_channels->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.plot_weight_as_pointsize"))
			m_plot_weights_pointsize->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.plot_weight_as_alpha"))
			m_plot_weights_alpha->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.auto_calc"))
			m_autocalc->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_DMI"))
			m_use_dmi->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_field"))
			m_use_field->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_temperature"))
			m_use_temperature->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_magffact"))
			m_use_formfact->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_weights"))
			m_use_weights->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.unite_degeneracies"))
			m_unite_degeneracies->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.ignore_annihilation"))
			m_ignore_annihilation->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.force_incommensurate"))
			m_force_incommensurate->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.calc_H"))
			m_hamiltonian_comp[0]->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.calc_Hp"))
			m_hamiltonian_comp[1]->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.calc_Hm"))
			m_hamiltonian_comp[2]->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_projector"))
			m_use_projector->setChecked(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.field_axis_h"))
			m_rot_axis[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.field_axis_k"))
			m_rot_axis[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.field_axis_l"))
			m_rot_axis[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.field_angle"))
			m_rot_angle->setValue(*optVal);
		bool found_sg = false;
		if(auto optVal = magdyn.get_optional<std::string>("config.spacegroup"))
		{
			// set space group by name
			int idx = m_comboSG->findText(optVal->c_str(), Qt::MatchContains);
			if(idx >= 0)
			{
				m_comboSG->setCurrentIndex(idx);
				found_sg = true;
			}
		}
		if(!found_sg)
		{
			if(auto optVal = magdyn.get_optional<int>("config.spacegroup_index"))
				m_comboSG->setCurrentIndex(*optVal);
		}
		if(auto optVal = magdyn.get_optional<t_real>("config.export_start_h"))
			m_exportStartQ[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.export_start_k"))
			m_exportStartQ[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.export_start_l"))
			m_exportStartQ[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.export_end_h"))
			m_exportEndQ[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.export_end_k"))
			m_exportEndQ[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.export_end_l"))
			m_exportEndQ[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_size>("config.export_num_points_1"))
			m_exportNumPoints[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_size>("config.export_num_points_2"))
			m_exportNumPoints[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_size>("config.export_num_points_3"))
			m_exportNumPoints[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<t_real>("config.couplings_max_dist"))
			m_maxdist->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<int>("config.couplings_max_supercell"))
			m_maxSC->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<int>("config.couplings_max_count"))
			m_maxcouplings->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<int>("config.sites_extcell_x"))
			m_extCell[0]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<int>("config.sites_extcell_y"))
			m_extCell[1]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<int>("config.sites_extcell_z"))
			m_extCell[2]->setValue(*optVal);
		if(auto optVal = magdyn.get_optional<bool>("config.use_genJ"))
		{
			if(!m_allow_general_J && *optVal)
			{
				QMessageBox::warning(this, "Magnetic Structure",
					"This file requires support for general exchange matrices J, "
					"please activate them in the preferences.");
			}
			else if(m_allow_general_J)
			{
				m_use_genJ->setChecked(*optVal);
			}
		}

		if(!m_dyn.Load(magdyn))
		{
			QMessageBox::critical(this, "Magnetic Dynamics", "Cannot load magdyn file.");
			return false;
		}

		// external field
		m_field_dir[0]->setValue(m_dyn.GetExternalField().dir[0]);
		m_field_dir[1]->setValue(m_dyn.GetExternalField().dir[1]);
		m_field_dir[2]->setValue(m_dyn.GetExternalField().dir[2]);
		m_field_mag->setValue(m_dyn.GetExternalField().mag);
		m_align_spins->setChecked(m_dyn.GetExternalField().align_spins);
		if(!m_use_field->isChecked())
			m_dyn.ClearExternalField();

		// ordering vector
		const t_vec_real& ordering = m_dyn.GetOrderingWavevector();
		if(ordering.size() == 3)
		{
			m_ordering[0]->setValue(ordering[0]);
			m_ordering[1]->setValue(ordering[1]);
			m_ordering[2]->setValue(ordering[2]);
		}

		// normal axis
		const t_vec_real& norm = m_dyn.GetRotationAxis();
		if(norm.size() == 3)
		{
			m_normaxis[0]->setValue(norm[0]);
			m_normaxis[1]->setValue(norm[1]);
			m_normaxis[2]->setValue(norm[2]);
		}

		// temperature
		t_real temp = m_dyn.GetTemperature();
		if(temp >= 0.)
			m_temperature->setValue(temp);
		if(!m_use_temperature->isChecked())
			m_dyn.SetTemperature(-1.);

		// form factor
		std::string ffact = m_dyn.GetMagneticFormFactor();
		if(ffact != "")
			m_ffact->setPlainText(ffact.c_str());
		if(!m_use_formfact->isChecked())
			m_dyn.SetMagneticFormFactor("");

		// clear old tables
		DelTabItem(m_sitestab, -1);
		DelTabItem(m_termstab, -1);
		DelTabItem(m_varstab, -1);
		DelTabItem(m_fieldstab, -1);
		DelTabItem(m_coordinatestab, -1);

		// variables
		for(const auto& var : m_dyn.GetVariables())
			AddVariableTabItem(-1, var.name, var.value);

		// crystal lattice
		const auto& xtal = m_dyn.GetCrystalLattice();
		m_xtallattice[0]->setValue(xtal[0]);
		m_xtallattice[1]->setValue(xtal[1]);
		m_xtallattice[2]->setValue(xtal[2]);
		m_xtalangles[0]->setValue(tl2::r2d<t_real>(xtal[3]));
		m_xtalangles[1]->setValue(tl2::r2d<t_real>(xtal[4]));
		m_xtalangles[2]->setValue(tl2::r2d<t_real>(xtal[5]));

		// scattering plane
		const auto* plane = m_dyn.GetScatteringPlane();
		m_scatteringplane[0]->setValue(plane[0][0]);
		m_scatteringplane[1]->setValue(plane[0][1]);
		m_scatteringplane[2]->setValue(plane[0][2]);
		m_scatteringplane[3]->setValue(plane[1][0]);
		m_scatteringplane[4]->setValue(plane[1][1]);
		m_scatteringplane[5]->setValue(plane[1][2]);

		// sync magnetic sites and additional entries
		SyncSitesFromKernel(magdyn.get_child_optional("atom_sites"));

		// sync exchange terms and additional entries
		SyncTermsFromKernel(magdyn.get_child_optional("exchange_terms"));

		// saved fields
		if(auto vars = magdyn.get_child_optional("saved_fields"); vars)
		{
			for(const auto &var : *vars)
			{
				t_real Bh = var.second.get<t_real>("direction_h", 0.);
				t_real Bk = var.second.get<t_real>("direction_k", 0.);
				t_real Bl = var.second.get<t_real>("direction_l", 0.);
				t_real Bmag = var.second.get<t_real>("magnitude", 0.);

				AddFieldTabItem(-1, Bh, Bk, Bl, Bmag);
			}
		}

		// saved coordinates
		if(auto vars = magdyn.get_child_optional("saved_coordinates"); vars)
		{
			for(const auto &var : *vars)
			{
				std::string name = var.second.get<std::string>("name", "");
				t_real h = var.second.get<t_real>("h", 0.);
				t_real k = var.second.get<t_real>("k", 0.);
				t_real l = var.second.get<t_real>("l", 0.);

				AddCoordinateTabItem(-1, name, h, k, l);
			}
		}

		if(m_groundstate_dlg)
			m_groundstate_dlg->SyncFromKernel();
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Magnetic Dynamics", ex.what());
		return false;
	}

	return true;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
/**
 * import magnetic structure configuration
 */
void MagDynDlg::ImportStructure()
{
	QString dirLast = m_sett->value("dir_struct", "").toString();
	QString filename = QFileDialog::getOpenFileName(
		this, "Import File", dirLast, "Magnetic Structure Files (*.xml)");
	if(filename == "" || !QFile::exists(filename))
		return;

	Clear();

	if(ImportStructure(filename))
	{
		m_sett->setValue("dir_struct", QFileInfo(filename).path());
		m_recent_struct.AddRecentFile(filename);
	}
}



/**
 * import magnetic structure configuration
 */
bool MagDynDlg::ImportStructure(const QString& filename)
{
	try
	{
		BOOST_SCOPE_EXIT(this_)
		{
			this_->m_ignoreCalc = false;
			if(this_->m_autocalc->isChecked())
			{
				this_->SyncToKernel();
			}
		} BOOST_SCOPE_EXIT_END
		m_ignoreCalc = true;

		// properties tree
		pt::ptree node;

		// load from file
		std::ifstream ifstr{filename.toStdString()};
		pt::read_xml(ifstr, node);

		// check signature
		if(auto optInfo = node.get_optional<std::string>("sfact.meta.info");
		   !optInfo || !(*optInfo==std::string{"magsfact_tool"} || *optInfo==std::string{"sfact_tool"}))
		{
			QMessageBox::critical(this, "Magnetic Structure", "Unrecognised structure file format.");
			return false;
		}

		const auto &sfact = node.get_child("sfact");

		if(auto optVal = sfact.get_optional<t_real>("xtal.a"))
			m_xtallattice[0]->setValue(*optVal);
		if(auto optVal = sfact.get_optional<t_real>("xtal.b"))
			m_xtallattice[1]->setValue(*optVal);
		if(auto optVal = sfact.get_optional<t_real>("xtal.c"))
			m_xtallattice[2]->setValue(*optVal);
		if(auto optVal = sfact.get_optional<t_real>("xtal.alpha"))
			m_xtalangles[0]->setValue(*optVal);
		if(auto optVal = sfact.get_optional<t_real>("xtal.beta"))
			m_xtalangles[1]->setValue(*optVal);
		if(auto optVal = sfact.get_optional<t_real>("xtal.gamma"))
			m_xtalangles[2]->setValue(*optVal);
		if(auto optVal = sfact.get_optional<int>("sg_idx"))
			m_comboSG->setCurrentIndex(*optVal);

		// spin structure
		if(auto nuclei = sfact.get_child_optional("nuclei"); nuclei)
		{
			for(const auto &nucl : *nuclei)
			{
				std::string name = nucl.second.get<std::string>("name", "n/a");
				std::string x = nucl.second.get<std::string>("x", "0");
				std::string y = nucl.second.get<std::string>("y", "0");
				std::string z = nucl.second.get<std::string>("z", "0");
				std::string M_mag = nucl.second.get<std::string>("M_mag", "1");
				std::string ReMx = nucl.second.get<std::string>("ReMx", "0");
				std::string ReMy = nucl.second.get<std::string>("ReMy", "0");
				std::string ReMz = nucl.second.get<std::string>("ReMz", "1");
				std::string rgb = nucl.second.get<std::string>("col", "auto");

				AddSiteTabItem(-1, name, 0,
					x, y, z,
				   ReMx, ReMy, ReMz, M_mag,
				   "auto", "auto", "auto",
				   rgb);
			}
		}

		// propagation vectors
		if(auto propvecs = sfact.get_child_optional("propvecs"); propvecs)
		{
			if(propvecs->size() >= 1)
			{
				// use first propagation vector
				t_real x = propvecs->begin()->second.get<t_real>("x", 0.);
				t_real y = propvecs->begin()->second.get<t_real>("y", 0.);
				t_real z = propvecs->begin()->second.get<t_real>("z", 0.);

				m_ordering[0]->setValue(x);
				m_ordering[1]->setValue(y);
				m_ordering[2]->setValue(z);
			}

			if(propvecs->size() > 1)
			{
				QMessageBox::warning(this, "Magnetic Structure",
					"Only one propagation vector is supported.");
			}
		}
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Magnetic Structure", ex.what());
		return false;
	}

	return true;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
/**
 * save current magnetic structure configuration
 */
void MagDynDlg::Save()
{
	const QString& curFile = m_recent.GetCurFile();
	if(curFile == "")
		SaveAs();
	else
		Save(curFile);
}



/**
 * save current magnetic structure configuration
 */
void MagDynDlg::SaveAs()
{
	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save File", dirLast, "Magnetic Dynamics Files (*.magdyn)");
	if(filename == "")
		return;

	if(Save(filename))
		SetCurrentFileAndDir(filename);
}



/**
 * save current magnetic structure configuration
 */
bool MagDynDlg::Save(const QString& filename)
{
	using t_item = tl2::NumericTableWidgetItem<t_real>;

	try
	{
		SyncToKernel();

		// properties tree
		pt::ptree magdyn;

		const char* user = std::getenv("USER");
		if(!user)
			user = "";

		magdyn.put<std::string>("meta.user", user);
		magdyn.put<std::string>("meta.url", "https://github.com/ILLGrenoble/takin");
		magdyn.put<std::string>("meta.doi", "https://doi.org/10.5281/zenodo.4117437");

		// save user comment as utf8 to avoid collisions with possible xml tags
		magdyn.put<std::string>("meta.notes", m_notes_dlg->GetNotes());

		// settings
		magdyn.put<t_real>("config.h_start", m_Q_start[0]->value());
		magdyn.put<t_real>("config.k_start", m_Q_start[1]->value());
		magdyn.put<t_real>("config.l_start", m_Q_start[2]->value());
		magdyn.put<t_real>("config.h_end", m_Q_end[0]->value());
		magdyn.put<t_real>("config.k_end", m_Q_end[1]->value());
		magdyn.put<t_real>("config.l_end", m_Q_end[2]->value());
		magdyn.put<t_real>("config.h", m_Q[0]->value());
		magdyn.put<t_real>("config.k", m_Q[1]->value());
		magdyn.put<t_real>("config.l", m_Q[2]->value());
		magdyn.put<t_size>("config.num_Q_points", m_num_points->value());
		magdyn.put<t_real>("config.weight_scale", m_weight_scale->value());
		magdyn.put<t_real>("config.weight_min", m_weight_min->value());
		magdyn.put<t_real>("config.weight_max", m_weight_max->value());
		magdyn.put<bool>("config.plot_channels", m_plot_channels->isChecked());
		magdyn.put<bool>("config.plot_weight_as_pointsize", m_plot_weights_pointsize->isChecked());
		magdyn.put<bool>("config.plot_weight_as_alpha", m_plot_weights_alpha->isChecked());
		magdyn.put<bool>("config.auto_calc", m_autocalc->isChecked());
		magdyn.put<bool>("config.use_DMI", m_use_dmi->isChecked());
		magdyn.put<bool>("config.use_genJ", m_allow_general_J && m_use_genJ->isChecked());
		magdyn.put<bool>("config.use_field", m_use_field->isChecked());
		magdyn.put<bool>("config.use_temperature", m_use_temperature->isChecked());
		magdyn.put<bool>("config.use_magffact", m_use_formfact->isChecked());
		magdyn.put<bool>("config.use_weights", m_use_weights->isChecked());
		magdyn.put<bool>("config.unite_degeneracies", m_unite_degeneracies->isChecked());
		magdyn.put<bool>("config.ignore_annihilation", m_ignore_annihilation->isChecked());
		magdyn.put<bool>("config.force_incommensurate", m_force_incommensurate->isChecked());
		magdyn.put<bool>("config.calc_H", m_hamiltonian_comp[0]->isChecked());
		magdyn.put<bool>("config.calc_Hp", m_hamiltonian_comp[1]->isChecked());
		magdyn.put<bool>("config.calc_Hm", m_hamiltonian_comp[2]->isChecked());
		magdyn.put<bool>("config.use_projector", m_use_projector->isChecked());
		magdyn.put<t_real>("config.field_axis_h", m_rot_axis[0]->value());
		magdyn.put<t_real>("config.field_axis_k", m_rot_axis[1]->value());
		magdyn.put<t_real>("config.field_axis_l", m_rot_axis[2]->value());
		magdyn.put<t_real>("config.field_angle", m_rot_angle->value());
		magdyn.put<std::string>("config.spacegroup", m_comboSG->currentText().toStdString());
		magdyn.put<t_real>("config.spacegroup_index", m_comboSG->currentIndex());
		magdyn.put<t_real>("config.export_start_h", m_exportStartQ[0]->value());
		magdyn.put<t_real>("config.export_start_k", m_exportStartQ[1]->value());
		magdyn.put<t_real>("config.export_start_l", m_exportStartQ[2]->value());
		magdyn.put<t_real>("config.export_end_h", m_exportEndQ[0]->value());
		magdyn.put<t_real>("config.export_end_k", m_exportEndQ[1]->value());
		magdyn.put<t_real>("config.export_end_l", m_exportEndQ[2]->value());
		magdyn.put<t_size>("config.export_num_points_1", m_exportNumPoints[0]->value());
		magdyn.put<t_size>("config.export_num_points_2", m_exportNumPoints[1]->value());
		magdyn.put<t_size>("config.export_num_points_3", m_exportNumPoints[2]->value());
		magdyn.put<t_real>("config.couplings_max_dist", m_maxdist->value());
		magdyn.put<int>("config.couplings_max_supercell", m_maxSC->value());
		magdyn.put<int>("config.couplings_max_count", m_maxcouplings->value());
		magdyn.put<int>("config.sites_extcell_x", m_extCell[0]->value());
		magdyn.put<int>("config.sites_extcell_y", m_extCell[1]->value());
		magdyn.put<int>("config.sites_extcell_z", m_extCell[2]->value());

		// save magnon calculator configuration
		if(!m_dyn.Save(magdyn))
		{
			QMessageBox::critical(this, "Magnetic Dynamics",
				"Cannot save magdyn file.");
			return false;
		}


		// saved fields
		for(int field_row = 0; field_row < m_fieldstab->rowCount(); ++field_row)
		{
			const auto* Bh = static_cast<t_item*>(m_fieldstab->item(field_row, COL_FIELD_H));
			const auto* Bk = static_cast<t_item*>(m_fieldstab->item(field_row, COL_FIELD_K));
			const auto* Bl = static_cast<t_item*>(m_fieldstab->item(field_row, COL_FIELD_L));
			const auto* Bmag = static_cast<t_item*>(m_fieldstab->item(field_row, COL_FIELD_MAG));

			boost::property_tree::ptree itemNode;
			itemNode.put<t_real>("direction_h", Bh ? Bh->GetValue() : 0.);
			itemNode.put<t_real>("direction_k", Bk ? Bk->GetValue() : 0.);
			itemNode.put<t_real>("direction_l", Bl ? Bl->GetValue() : 0.);
			itemNode.put<t_real>("magnitude", Bmag ? Bmag->GetValue() : 0.);

			magdyn.add_child("saved_fields.field", itemNode);
		}

		// get site entries for putting additional infos
		if(auto sites = magdyn.get_child_optional("atom_sites"); sites)
		{
			auto siteiter = (*sites).begin();

			for(std::size_t site_idx = 0; site_idx < std::size_t(m_sitestab->rowCount()); ++site_idx)
			{
				// set additional data from site entry
				if(site_idx >= sites->size())
					break;

				// write colour
				std::string rgb = m_sitestab->item(site_idx, COL_SITE_RGB)->text().toStdString();
				siteiter->second.put<std::string>("colour", rgb);

				std::advance(siteiter, 1);
			}
		}

		// get exchange terms entries for putting additional infos
		if(auto terms = magdyn.get_child_optional("exchange_terms"); terms)
		{
			auto termiter = (*terms).begin();

			for(std::size_t term_idx = 0; term_idx < std::size_t(m_termstab->rowCount()); ++term_idx)
			{
				// set additional data from exchange term entry
				if(term_idx >= terms->size())
					break;

				// write colour
				std::string rgb = m_termstab->item(term_idx, COL_XCH_RGB)->text().toStdString();
				termiter->second.put<std::string>("colour", rgb);

				std::advance(termiter, 1);
			}
		}

		// saved coordinates
		for(int coord_row = 0; coord_row < m_coordinatestab->rowCount(); ++coord_row)
		{
			std::string name = m_coordinatestab->item(coord_row, COL_COORD_NAME)->text().toStdString();
			const auto* h = static_cast<t_item*>(m_coordinatestab->item(coord_row, COL_COORD_H));
			const auto* k = static_cast<t_item*>(m_coordinatestab->item(coord_row, COL_COORD_K));
			const auto* l = static_cast<t_item*>(m_coordinatestab->item(coord_row, COL_COORD_L));

			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", name);
			itemNode.put<t_real>("h", h ? h->GetValue() : 0.);
			itemNode.put<t_real>("k", k ? k->GetValue() : 0.);
			itemNode.put<t_real>("l", l ? l->GetValue() : 0.);

			magdyn.add_child("saved_coordinates.coordinate", itemNode);
		}


		pt::ptree node;
		node.put_child("magdyn", magdyn);

		// save to file
		std::ofstream ofstr{filename.toStdString()};
		if(!ofstr)
		{
			QMessageBox::critical(this, "Magnetic Dynamics",
				"Cannot open file for writing.");
			return false;
		}

		ofstr.precision(g_prec);
		pt::write_xml(ofstr, node,
			pt::xml_writer_make_settings('\t', 1, std::string{"utf-8"}));
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Magnetic Dynamics", ex.what());
		return false;
	}

	return true;
}
// --------------------------------------------------------------------------------
