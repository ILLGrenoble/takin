/**
 * magnetic dynamics -- minimise ground state energy
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

#include "ground_state.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QTableWidgetItem>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QMessageBox>

#include <string>
#include <sstream>
#include <unordered_set>

#include "tlibs2/libs/qt/numerictablewidgetitem.h"
using t_numitem = tl2::NumericTableWidgetItem<t_real>;



/**
 * columns of the sites table
 */
enum : int
{
	COL_SPIN_NAME = 0,             // label
	COL_SPIN_PHI, COL_SPIN_THETA,  // spin direction
	COL_SPIN_U, COL_SPIN_V,        // spin direction
	COL_SPIN_U_FIXED,              // don't minimise u?
	COL_SPIN_V_FIXED,              // don't minimise v?

	NUM_SPIN_COLS
};



/**
 * shows the 3d view of the magnetic structure
 */
GroundStateDlg::GroundStateDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Ground State Energy Minimisation");
	setSizeGripEnabled(true);

	m_spinstab = new QTableWidget(this);
	m_spinstab->setShowGrid(true);
	m_spinstab->setAlternatingRowColors(true);
	m_spinstab->setSortingEnabled(true);
	m_spinstab->setMouseTracking(true);
	m_spinstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_spinstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_spinstab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	m_spinstab->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_spinstab->verticalHeader()->setVisible(true);

	m_spinstab->setColumnCount(NUM_SPIN_COLS);

	m_spinstab->setHorizontalHeaderItem(COL_SPIN_NAME, new QTableWidgetItem{"Name"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_PHI, new QTableWidgetItem{"Spin φ"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_THETA, new QTableWidgetItem{"Spin θ"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_U, new QTableWidgetItem{"Spin u"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_V, new QTableWidgetItem{"Spin v"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_U_FIXED, new QTableWidgetItem{"Fix u"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_V_FIXED, new QTableWidgetItem{"Fix v"});

	m_spinstab->setColumnWidth(COL_SPIN_NAME, 125);
	m_spinstab->setColumnWidth(COL_SPIN_PHI, 90);
	m_spinstab->setColumnWidth(COL_SPIN_THETA, 90);
	m_spinstab->setColumnWidth(COL_SPIN_U, 90);
	m_spinstab->setColumnWidth(COL_SPIN_V, 90);
	m_spinstab->setColumnWidth(COL_SPIN_U_FIXED, 50);
	m_spinstab->setColumnWidth(COL_SPIN_V_FIXED, 50);

	m_btnFromKernel = new QPushButton("Get Spins", this);
	m_btnToKernel = new QPushButton("Set Spins", this);
	m_btnMinimise = new QPushButton("Minimise", this);
	QPushButton *btnOk = new QPushButton("OK", this);

	m_btnFromKernel->setToolTip("Fetch spins from main dialog.");
	m_btnToKernel->setToolTip("Send spins back to main dialog.");

	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	int y = 0;
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);
	grid->addWidget(m_spinstab, y++, 0, 1, 4);
	grid->addWidget(m_btnFromKernel, y, 0, 1, 1);
	grid->addWidget(m_btnToKernel, y, 1, 1, 1);
	grid->addWidget(m_btnMinimise, y, 2, 1, 1);
	grid->addWidget(btnOk, y++, 3, 1, 1);
	grid->addWidget(m_status, y, 0, 1, 4);

	if(m_sett && m_sett->contains("ground_state/geo"))
		restoreGeometry(m_sett->value("ground_state/geo").toByteArray());
	else
		resize(640, 480);

	connect(m_spinstab, &QTableWidget::itemChanged, this, &GroundStateDlg::SpinsTableItemChanged);
	connect(btnOk, &QAbstractButton::clicked, this, &QDialog::accept);
	connect(m_btnMinimise, &QAbstractButton::clicked, this, &GroundStateDlg::Minimise);
	connect(m_btnToKernel, &QAbstractButton::clicked, this, &GroundStateDlg::SyncToKernel);
	connect(m_btnFromKernel, &QAbstractButton::clicked, [this]()
	{
		GroundStateDlg::SyncFromKernel();
	});

	EnableMinimisation(true);
}



GroundStateDlg::~GroundStateDlg()
{
	m_stop_request = true;

	if(m_thread)
	{
		m_thread->join();
		m_thread = nullptr;
	}
}



/**
 * a spin property was edited in the table
 */
void GroundStateDlg::SpinsTableItemChanged(QTableWidgetItem *item)
{
	if(!item)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_spinstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_spinstab->blockSignals(true);

	const int row = m_spinstab->row(item);
	const int col = m_spinstab->column(item);
	const int num_rows = m_spinstab->rowCount();
	const int num_cols = m_spinstab->columnCount();

	if(row < 0 || col < 0 || row >= num_rows || col >= num_cols)
		return;

	switch(col)
	{
		// phi was edited
		case COL_SPIN_PHI:
		{
			// get phi and theta
			auto *item_phi = static_cast<t_numitem*>(item);
			auto *item_theta = static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_THETA));

			// calculate new u and v
			const auto [ u, v ] = tl2::sph_to_uv<t_real>(
				tl2::d2r(item_phi->GetValue()), tl2::d2r(item_theta->GetValue()));

			// set new u and v
			static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_U))->SetValue(u);
			static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_V))->SetValue(v);

			break;
		}

		// theta was edited
		case COL_SPIN_THETA:
		{
			// get phi and theta
			auto *item_theta = static_cast<t_numitem*>(item);
			auto *item_phi = static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_PHI));

			// calculate new u and v
			const auto [ u, v ] = tl2::sph_to_uv<t_real>(
				tl2::d2r(item_phi->GetValue()), tl2::d2r(item_theta->GetValue()));

			// set new u and v
			static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_U))->SetValue(u);
			static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_V))->SetValue(v);

			break;
		}

		// u was edited
		case COL_SPIN_U:
		{
			// get u and v
			auto *item_u = static_cast<t_numitem*>(item);
			auto *item_v = static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_V));

			// calculate new phi and theta
			const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(
				item_u->GetValue(), item_v->GetValue());

			// set new phi and theta
			static_cast<t_numitem*>(
				m_spinstab->item(row, COL_SPIN_PHI))->SetValue(tl2::r2d(phi));
			static_cast<t_numitem*>(
				m_spinstab->item(row, COL_SPIN_THETA))->SetValue(tl2::r2d(theta));

			break;
		}

		// v was edited
		case COL_SPIN_V:
		{
			// get u and v
			auto *item_v = static_cast<t_numitem*>(item);
			auto *item_u = static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_U));

			// calculate new phi and theta
			const auto [ phi, theta ] = tl2::uv_to_sph<t_real>(
				item_u->GetValue(), item_v->GetValue());

			// set new phi and theta
			static_cast<t_numitem*>(
				m_spinstab->item(row, COL_SPIN_PHI))->SetValue(tl2::r2d(phi));
			static_cast<t_numitem*>(
				m_spinstab->item(row, COL_SPIN_THETA))->SetValue(tl2::r2d(theta));

			break;
		}
	}

	UpdateSpinFromTable(row);
	CalcGroundStateEnergy();
}



/**
 * set the kernel's spins to the ones given in the table
 */
void GroundStateDlg::UpdateSpinsFromTable()
{
	if(!m_dyn)
		return;

	for(int row = 0; row < m_spinstab->rowCount(); ++row)
		UpdateSpinFromTable(row);
}



/**
 * set the kernel's spins to the ones given in the table
 */
void GroundStateDlg::UpdateSpinFromTable(int row)
{
	if(!m_dyn)
		return;

	std::string site_name = m_spinstab->item(row, COL_SPIN_NAME)->text().toStdString();

	// get phi and theta
	auto *item_phi = static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_PHI));
	auto *item_theta = static_cast<t_numitem*>(m_spinstab->item(row, COL_SPIN_THETA));

	if(!item_phi || !item_theta)
		return;

	// calculate spin (x, y, z)
	const auto [ x, y, z ] = tl2::sph_to_cart<t_real>(
		1., tl2::d2r(item_phi->GetValue()), tl2::d2r(item_theta->GetValue()));

	// set new spin (x, y, z)
	if(t_size idx = m_dyn->GetMagneticSiteIndex(site_name);
		idx < m_dyn->GetMagneticSitesCount())
	{
		t_site& site = m_dyn->GetMagneticSites()[idx];

		site.spin_dir[0] = tl2::var_to_str(x, g_prec);
		site.spin_dir[1] = tl2::var_to_str(y, g_prec);
		site.spin_dir[2] = tl2::var_to_str(z, g_prec);
		site.spin_dir_calc = tl2::create<t_vec_real>({ x, y, z });

		m_dyn->CalcMagneticSite(site);
	}
}



/**
 * set a pointer to the main magdyn kernel and synchronise the spins
 */
void GroundStateDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn_kern = dyn;
	if(!m_dyn)
		SyncFromKernel();
}



/**
 * get spin configuration from kernel
 */
void GroundStateDlg::SyncFromKernel(const t_magdyn *dyn,
	const std::unordered_set<std::string> *fixed_spins)
{
	// if no kernel is given, use the main one
	if(!dyn)
		dyn = m_dyn_kern;
	if(!dyn)
		return;

	// copy the data (if it's not the same structure anyway)
	if(dyn != &*m_dyn)
		m_dyn = *dyn;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_spinstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_spinstab->blockSignals(true);

	m_spinstab->setSortingEnabled(false);
	m_spinstab->clearContents();
	m_spinstab->setRowCount(0);

	for(const t_site& site : m_dyn->GetMagneticSites())
	{
		const t_vec_real& S = site.spin_dir_calc;
		const auto [ rho, phi, theta ] =  tl2::cart_to_sph<t_real>(S[0], S[1], S[2]);
		const auto [ u, v ] = tl2::sph_to_uv<t_real>(phi, theta);

		int row = m_spinstab->rowCount();
		m_spinstab->insertRow(row);

		QTableWidgetItem *item_name = new QTableWidgetItem(site.name.c_str());
		QTableWidgetItem *item_phi = new t_numitem(tl2::r2d(phi));
		QTableWidgetItem *item_theta = new t_numitem(tl2::r2d(theta));
		QTableWidgetItem *item_u = new t_numitem(u);
		QTableWidgetItem *item_v = new t_numitem(v);

		// write-protect site identifier
		item_name->setFlags(item_name->flags() & ~Qt::ItemIsEditable);

		QCheckBox* u_fixed = new QCheckBox(this);
		QCheckBox* v_fixed = new QCheckBox(this);
		u_fixed->setChecked(true);
		v_fixed->setChecked(false);

		m_spinstab->setItem(row, COL_SPIN_NAME, item_name);
		m_spinstab->setItem(row, COL_SPIN_PHI, item_phi);
		m_spinstab->setItem(row, COL_SPIN_THETA, item_theta);
		m_spinstab->setItem(row, COL_SPIN_U, item_u);
		m_spinstab->setItem(row, COL_SPIN_V, item_v);
		m_spinstab->setCellWidget(row, COL_SPIN_U_FIXED, u_fixed);
		m_spinstab->setCellWidget(row, COL_SPIN_V_FIXED, v_fixed);

		// keep fixed spins
		if(fixed_spins)
		{
			u_fixed->setChecked(fixed_spins->find(site.name + "_phi") != fixed_spins->end());
			v_fixed->setChecked(fixed_spins->find(site.name + "_theta") != fixed_spins->end());
		}
	}

	m_spinstab->setSortingEnabled(true);
	CalcGroundStateEnergy();
}



/**
 * send spin configuration to kernel
 */
void GroundStateDlg::SyncToKernel()
{
	if(!m_dyn)
		return;

	emit SpinsUpdated(&*m_dyn);
}



void GroundStateDlg::EnableMinimisation(bool enable)
{
	if(enable)
	{
		m_btnMinimise->setText("Minimise");
		m_btnFromKernel->setEnabled(true);
		m_btnToKernel->setEnabled(true);
		m_spinstab->setEnabled(true);

		m_btnMinimise->setToolTip("Start minimisation of ground state energy.");
	}
	else
	{
		m_btnMinimise->setText("Stop");
		m_btnFromKernel->setEnabled(false);
		m_btnToKernel->setEnabled(false);
		m_spinstab->setEnabled(false);

		m_btnMinimise->setToolTip("Stop minimisation of ground state energy.");
	}
}



/**
 * minimise the ground state energy
 */
void GroundStateDlg::Minimise()
{
	if(!m_dyn)
		return;

	if(m_running)
	{
		m_stop_request = true;
		m_status->setText("Stopping calculation.");
		return;
	}

	if(m_thread)
	{
		m_thread->join();
		m_thread = nullptr;
	}

	m_status->setText("Calculating ground state.");
	m_running = true;
	EnableMinimisation(false);

	m_thread = std::make_unique<std::thread>([this]()
	{
		BOOST_SCOPE_EXIT(this_)
		{
			QMetaObject::invokeMethod(this_, [this_]()
			{
				this_->CalcGroundStateEnergy();
				this_->EnableMinimisation(true);
			});
			this_->m_running = false;
		} BOOST_SCOPE_EXIT_END

		m_stop_request = false;
		bool cancelled = false;

		// set fixed spin parameters
		std::unordered_set<std::string> fixed_spins;

		for(int row = 0; row < m_spinstab->rowCount(); ++row)
		{
			QCheckBox *u_fixed = reinterpret_cast<QCheckBox*>(
				m_spinstab->cellWidget(row, COL_SPIN_U_FIXED));
			if(u_fixed && u_fixed->isChecked())
			{
				std::string fixed_name = m_spinstab->item(
					row, COL_SPIN_NAME)->text().toStdString() + "_phi";
				fixed_spins.emplace(std::move(fixed_name));
			}

			QCheckBox *v_fixed = reinterpret_cast<QCheckBox*>(
				m_spinstab->cellWidget(row, COL_SPIN_V_FIXED));
			if(v_fixed && v_fixed->isChecked())
			{
				std::string fixed_name = m_spinstab->item(
					row, COL_SPIN_NAME)->text().toStdString() + "_theta";
				fixed_spins.emplace(std::move(fixed_name));
			}
		}

		try
		{
			// minimise
			if(!m_dyn->CalcGroundState(&fixed_spins, true, &m_stop_request))
			{
				QMetaObject::invokeMethod(this, [this]()
				{
					this->ShowError("Could not find ground state.");
				});
			}
		}
		catch(const tl2::StopRequestException&)
		{
			cancelled = true;
		}

		if(!cancelled)
		{
			// set new parameters
			QMetaObject::invokeMethod(this, [this, fixed_spins]()
			{
				this->SyncFromKernel(&*m_dyn, &fixed_spins);
			});
		}
	});
}



void GroundStateDlg::CalcGroundStateEnergy()
{
	if(!m_dyn)
		return;

	t_real E0 = m_dyn->CalcGroundStateEnergy();
	std::ostringstream ostr;
	ostr.precision(g_prec_gui);
	ostr << "Ground state energy: E0 = " << E0 << ".";

	m_status->setText(ostr.str().c_str());
}



void GroundStateDlg::ShowError(const char* msg)
{
	QMessageBox::critical(this, "Magnetic Dynamics", msg);
}



/**
 * dialog is closing
 */
void GroundStateDlg::accept()
{
	if(m_sett)
		m_sett->setValue("ground_state/geo", saveGeometry());

	QDialog::accept();
}
