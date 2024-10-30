/**
 * magnetic dynamics
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

#ifndef __MAG_DYN_GUI_H__
#define __MAG_DYN_GUI_H__

#include <QtCore/QSettings>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTableWidgetItem>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMenu>
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	#include <QtWidgets/QAction>
#else
	#include <QtGui/QAction>
#endif

#include <qcustomplot.h>

#include <vector>
#include <unordered_map>
#include <optional>

#include <boost/property_tree/ptree.hpp>

#include "tlibs2/libs/magdyn.h"
#include "tlibs2/libs/qt/numerictablewidgetitem.h"
#include "tlibs2/libs/qt/recent.h"

#include "gui_defs.h"
#include "graph.h"
#include "table_import.h"
#include "struct_plot.h"
#include "ground_state.h"
#include "trafos.h"
#include "notes.h"
#include "infos.h"



/**
 * combo box showing the magnetic sites and sorting according to their index
 */
struct SitesComboBox : public QComboBox, QTableWidgetItem
{
	SitesComboBox() = default;
	virtual ~SitesComboBox() = default;


	virtual bool operator<(const QTableWidgetItem& item) const override
	{
		const SitesComboBox* combo = dynamic_cast<const SitesComboBox*>(&item);
		if(!combo)
			return true;

		return currentIndex() < combo->currentIndex();
	}
};



/**
 * magnon calculation dialog
 */
class MagDynDlg : public QDialog
{ Q_OBJECT
public:
	MagDynDlg(QWidget* pParent = nullptr);
	virtual ~MagDynDlg();

	MagDynDlg(const MagDynDlg&) = delete;
	const MagDynDlg& operator=(const MagDynDlg&) = delete;


protected:
	QSettings *m_sett{};
	QMenuBar *m_menu{};
	QSplitter *m_split_inout{};
	QLabel *m_statusFixed{}, *m_status{};
	QProgressBar *m_progress{};
	QPushButton* m_btnStartStop{};

	QAction *m_autocalc{};
	QAction *m_use_dmi{}, *m_use_genJ{};
	QAction *m_use_field{}, *m_use_temperature{};
	QAction *m_use_formfact{};
	QAction *m_use_weights{}, *m_use_projector{};
	QAction *m_unite_degeneracies{};
	QAction *m_ignore_annihilation{};
	QAction *m_force_incommensurate{};
	QAction *m_plot_channels{};
	QMenu *m_menuChannels{};
	QAction *m_plot_channel[3]{};
	QAction *m_plot_weights_pointsize{};
	QAction *m_plot_weights_alpha{};
	QAction *m_hamiltonian_comp[3]{};

	// recently opened files
	tl2::RecentFiles m_recent{}, m_recent_struct{};
	QMenu *m_menuOpenRecent{}, *m_menuImportStructRecent{};

	// function to call for the recent file->open menu items
	std::function<bool(const QString& filename)> m_open_func
		= [this](const QString& filename) -> bool
	{
		this->Clear();
		this->SetCurrentFile(filename);
		return this->Load(filename);
	};

	// function to call for the recent file->import menu items
	std::function<bool(const QString& filename)> m_import_struct_func
		= [this](const QString& filename) -> bool
	{
		this->Clear();
		return this->ImportStructure(filename);
	};

	QGridLayout *m_maingrid{};

	// tabs
	QTabWidget *m_tabs_in{}, *m_tabs_out{};

	// panels
	QWidget *m_sitespanel{}, *m_termspanel{};
	QWidget *m_samplepanel{}, *m_sampleenviropanel{};
	QWidget *m_disppanel{}, *m_hamiltonianpanel{};
	QWidget *m_coordinatespanel{}, *m_exportpanel{};
	QWidget *m_varspanel{};

	// sites
	QTableWidget *m_sitestab{};
	QSpinBox *m_extCell[3]{nullptr, nullptr, nullptr};

	// terms, ordering vector, and rotation axis
	QTableWidget *m_termstab{};
	QDoubleSpinBox *m_maxdist{};
	QSpinBox *m_maxSC{};
	QSpinBox *m_maxcouplings{};
	QDoubleSpinBox *m_ordering[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_normaxis[3]{nullptr, nullptr, nullptr};

	// sample
	QDoubleSpinBox *m_xtallattice[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_xtalangles[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_scatteringplane[2*3]{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
	std::vector<std::vector<t_mat_real>> m_SGops{};
	QComboBox *m_comboSG{};
	QPlainTextEdit *m_ffact{};    // magnetic form factor formula

	// variables
	QTableWidget *m_varstab{};

	// dispersion
	QCustomPlot *m_plot{};
	std::vector<GraphWithWeights*> m_graphs{};
	QMenu *m_menuDisp{};
	QDoubleSpinBox *m_Q_start[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_Q_end[3]{nullptr, nullptr, nullptr};
	QSpinBox *m_num_points{};
	QDoubleSpinBox *m_weight_scale{}, *m_weight_min{}, *m_weight_max{};

	// hamiltonian
	QTextEdit *m_hamiltonian{};
	QDoubleSpinBox *m_Q[3]{nullptr, nullptr, nullptr};

	// external magnetic field
	QDoubleSpinBox *m_field_dir[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_field_mag{};
	QCheckBox *m_align_spins{};
	QDoubleSpinBox *m_rot_axis[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_rot_angle{};
	QTableWidget *m_fieldstab{};

	// temperature
	QDoubleSpinBox *m_temperature{};

	// coordinates
	QTableWidget *m_coordinatestab{};

	// export
	QDoubleSpinBox *m_exportStartQ[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_exportEndQ[3]{nullptr, nullptr, nullptr};
	QSpinBox *m_exportNumPoints[3]{nullptr, nullptr, nullptr};
	QComboBox *m_exportFormat{nullptr};

	// magnon dynamics calculator
	t_magdyn m_dyn{};

	// settings dialog
	QDialog *m_settings_dlg{};

	// dialogs
	TableImportDlg *m_table_import_dlg{};  // table import dialog
	NotesDlg *m_notes_dlg{};               // notes dialog
	TrafoCalculator *m_trafos{};           // trafo calculator
	InfoDlg *m_info_dlg{};                 // info dialog
	StructPlotDlg *m_structplot_dlg{};     // magnetic structure plotter
	GroundStateDlg *m_groundstate_dlg{};   // ground state minimiser


protected:
	// set up gui
	void CreateMainWindow();
	void CreateMenuBar();

	// set up dialogs
	void ShowInfoDlg(bool only_create = false);
	void ShowNotesDlg(bool only_create = false);
	void ShowStructPlotDlg(bool only_create = false);
	void ShowGroundStateDlg(bool only_create = false);
	void InitSettingsDlg();
	void InitSettings();

	// set up input panels
	void CreateSitesPanel();
	void CreateExchangeTermsPanel();
	void CreateSamplePanel();
	void CreateVariablesPanel();
	void CreateSampleEnvPanel();

	// set up output panels
	void CreateDispersionPanel();
	void CreateHamiltonPanel();
	void CreateExportPanel();
	void CreateCoordinatesPanel();

	// general table operations
	void MoveTabItemUp(QTableWidget *pTab);
	void MoveTabItemDown(QTableWidget *pTab);
	void ShowTableContextMenu(QTableWidget *pTab,
		QMenu *pMenu, QMenu *pMenuNoItem, const QPoint& pt);
	std::vector<int> GetSelectedRows(QTableWidget *pTab,
		bool sort_reversed = false) const;

	// add a site to the table
	void AddSiteTabItem(int row = -1,
		const std::string& name = "", t_size sym_idx = 0,
		const std::string& x = "0", const std::string& y = "0", const std::string& z = "0",
		const std::string& sx = "0",
		const std::string& sy = "0",
		const std::string& sz = "1",
		const std::string&  S = "1",
		const std::string& sox = "auto",
		const std::string& soy = "auto",
		const std::string& soz = "auto",
		const std::string& rgb = "auto");

	// add a coupling to the table
	void AddTermTabItem(int row = -1,
		const std::string& name = "", t_size sym_idx = 0,
		const std::string& atom_1 = "", const std::string& atom_2 = "",
		const std::string& dist_x = "0", const std::string& dist_y = "0", const std::string& dist_z = "0",
		const std::string& J = "0",
		const std::string& dmi_x = "0", const std::string& dmi_y = "0", const std::string& dmi_z = "0",
		const std::string& gen_xx = "0", const std::string& gen_xy = "0", const std::string& gen_xz = "0",
		const std::string& gen_yx = "0", const std::string& gen_yy = "0", const std::string& gen_yz = "0",
		const std::string& gen_zx = "0", const std::string& gen_zy = "0", const std::string& gen_zz = "0",
		const std::string& rgb = "#00bf00");

	void SyncSiteComboBoxes();
	void SyncSiteComboBox(SitesComboBox* combo, const std::string& selected_site);
	SitesComboBox* CreateSitesComboBox(const std::string& selected_site);

	// add a variable to the table
	void AddVariableTabItem(int row = -1,
		const std::string& name = "",
		const t_cplx& var = t_cplx{0, 0});

	// add a coordinate to the table
	void AddCoordinateTabItem(int row = -1,
		const std::string& name = "",
		t_real h = 0., t_real k = 0., t_real l = 1.);

	// add a magnetic field to the table
	void AddFieldTabItem(int row = -1,
		t_real Bh = 0., t_real Bk = 0., t_real Bl = 1.,
		t_real Bmag = 1.);

	void SetCurrentField();
	void SetCurrentCoordinate(int which = 0);

	std::vector<t_magdyn::Variable> GetVariables() const;
	t_size ReplaceValueWithVariable(const std::string& var, const t_cplx& val);
	t_size ReplaceValuesWithVariables();

	void DelTabItem(QTableWidget *pTab, int begin=-2, int end=-2);
	void UpdateVerticalHeader(QTableWidget *pTab);

	void SitesTableItemChanged(QTableWidgetItem *item);
	void TermsTableItemChanged(QTableWidgetItem *item);
	void VariablesTableItemChanged(QTableWidgetItem *item);

	void ClearDispersion(bool replot = false);
	void Clear();

	void Load();
	void Save();
	void SaveAs();

	void ImportStructure();
	void ExportToSunny();
	void ExportToSpinW();
	void ExportSQE();

	void SavePlotFigure();
	void SaveDispersion(bool as_scr = false);
	void SaveMultiDispersion(bool as_scr = false);

	void MirrorAtoms();
	void RotateField(bool ccw = true);
	void GenerateSitesFromSG();
	void GenerateCouplingsFromSG();
	void GeneratePossibleCouplings();
	void ExtendStructure();
	const std::vector<t_mat_real>& GetSymOpsForCurrentSG(bool show_err = true) const;

	// transfer sites from the kernel
	void SyncSitesFromKernel(boost::optional<const boost::property_tree::ptree&> extra_infos = boost::none);
	void SyncTermsFromKernel(boost::optional<const boost::property_tree::ptree&> extra_infos = boost::none);
	void SyncSymmetryIndicesFromKernel();
	void SyncToKernel();         // transfer all data to the kernel
	void CalcSymmetryIndices();  // assign symmetry groups to sites and couplings
	void CalcAll();              // syncs sites and terms and calculates all dynamics
	void SetKernel(const t_magdyn* dyn, bool sync_sites = true, bool sync_terms = true, bool sync_idx = true);

	void PlotDispersion();

	void PlotMouseMove(QMouseEvent* evt);
	void PlotMousePress(QMouseEvent* evt);

	virtual void mousePressEvent(QMouseEvent *evt) override;
	virtual void closeEvent(QCloseEvent *evt) override;
	virtual void dragEnterEvent(QDragEnterEvent *evt) override;
	virtual void dropEvent(QDropEvent *evt) override;

	// table importer
	void ShowTableImporter();
	void ImportAtoms(const std::vector<TableImportAtom>&, bool clear_existing = true);
	void ImportCouplings(const std::vector<TableImportCoupling>&, bool clear_existing = true);

	// disable/enable gui input for threaded operations
	void EnableInput();
	void DisableInput();

	// get the site corresponding to the given table index
	const t_magdyn::MagneticSite* GetSiteFromTableIndex(int idx) const;

	// get the coupling corresponding to the given table index
	const t_magdyn::ExchangeTerm* GetTermFromTableIndex(int idx) const;


public:
	bool Load(const QString& filename, bool calc_dynamics = true);
	bool Save(const QString& filename);

	bool ImportStructure(const QString& filename);
	bool ExportToSunny(const QString& filename);
	bool ExportToSpinW(const QString& filename);
	bool ExportSQE(const QString& filename);

	void SetCurrentFileAndDir(const QString& filename);
	void SetCurrentFile(const QString& filename);

	void SetNumQPoints(t_size num_Q_pts);
	void SetCoordinates(const t_vec_real& Qi, const t_vec_real& Qf, bool calc_dynamics = true);

	void CalcDispersion();
	void CalcHamiltonian();


private:
	int m_sites_cursor_row { -1 };
	int m_terms_cursor_row { -1 };
	int m_variables_cursor_row { -1 };
	int m_fields_cursor_row { -1 };
	int m_coordinates_cursor_row { -1 };

	bool m_ignoreCalc { false };
	bool m_ignoreSitesCalc { false };
	bool m_stopRequested { false };
	bool m_startEnabled { true };     // "calc" behaves as start or stop button?

	// optional features
	bool m_allow_ortho_spin { false };
	bool m_allow_general_J { false };

	// data for dispersion plot
	QVector<qreal> m_qs_data{}, m_Es_data{}, m_ws_data{};
	QVector<qreal> m_qs_data_channel[3]{}, m_Es_data_channel[3]{}, m_ws_data_channel[3]{};
	t_size m_Q_idx{};                 // plot x axis
	t_real m_Q_min{}, m_Q_max{};      // plot x axis range


public slots:
	void SelectSite(const std::string& site);
	void DeleteSite(const std::string& site);
	void FlipSiteSpin(const std::string& site);

	void SelectTerm(const std::string& term);
	void DeleteTerm(const std::string& term);
};


#endif
