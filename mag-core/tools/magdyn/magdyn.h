/**
 * magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
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

#include "tlibs2/libs/magdyn.h"
#include "tlibs2/libs/qt/numerictablewidgetitem.h"
#include "tlibs2/libs/qt/recent.h"
#include "tlibs2/libs/qt/glplot.h"

#include "defs.h"
#include "graph.h"
#include "table_import.h"
#include "trafos.h"

using namespace tl2_mag;


// let the user explicitly specify the orthogonal spin
//#define MAGDYN_ALLOW_SPIN_ORTHO_SETTABLE

// make the general interaction matrix accessible
#define MAGDYN_ALLOW_GENERAL_J


// magnon calculation core
using t_magdyn = MagDyn<t_mat, t_vec, t_mat_real, t_vec_real, t_cplx, t_real, t_size>;



/**
 * export file types
 */
enum : int
{
	EXPORT_HDF5 = 0,
	EXPORT_GRID = 1,
	EXPORT_TEXT = 2,
};


/**
 * columns of the sites table
 */
enum : int
{
	COL_SITE_NAME = 0,
	COL_SITE_POS_X, COL_SITE_POS_Y, COL_SITE_POS_Z,
	COL_SITE_SPIN_X, COL_SITE_SPIN_Y, COL_SITE_SPIN_Z,
	COL_SITE_SPIN_MAG,
	COL_SITE_RGB,
#ifdef MAGDYN_ALLOW_SPIN_ORTHO_SETTABLE
	COL_SITE_SPIN_ORTHO_X, COL_SITE_SPIN_ORTHO_Y, COL_SITE_SPIN_ORTHO_Z,
#endif

	NUM_SITE_COLS
};


/**
 * columns of the exchange terms table
 */
enum : int
{
	COL_XCH_NAME = 0,
	COL_XCH_ATOM1_IDX, COL_XCH_ATOM2_IDX,
	COL_XCH_DIST_X, COL_XCH_DIST_Y, COL_XCH_DIST_Z,
	COL_XCH_INTERACTION,                             // isotropic exchange interaction
	COL_XCH_DMI_X, COL_XCH_DMI_Y, COL_XCH_DMI_Z,     // antisymmetric DMI
#ifdef MAGDYN_ALLOW_GENERAL_J
	COL_XCH_GEN_XX, COL_XCH_GEN_XY, COL_XCH_GEN_XZ,  //
	COL_XCH_GEN_YX, COL_XCH_GEN_YY, COL_XCH_GEN_YZ,  // general interaction matrix
	COL_XCH_GEN_ZX, COL_XCH_GEN_ZY, COL_XCH_GEN_ZZ,  //
#endif
	COL_XCH_RGB,

	NUM_XCH_COLS
};


/**
 * columns of the variables table
 */
enum : int
{
	COL_VARS_NAME = 0,
	COL_VARS_VALUE_REAL,
	COL_VARS_VALUE_IMAG,

	NUM_VARS_COLS
};


/**
 * columns of the table with saved fields
 */
enum : int
{
	COL_FIELD_H = 0, COL_FIELD_K, COL_FIELD_L,
	COL_FIELD_MAG,

	NUM_FIELD_COLS
};


/**
 * columns of the table with saved Q coordinates
 */
enum : int
{
	COL_COORD_HI = 0, COL_COORD_KI, COL_COORD_LI,
	COL_COORD_HF, COL_COORD_KF, COL_COORD_LF,

	NUM_COORD_COLS
};


/**
 * infos for magnetic sites
 */
struct AtomSiteInfo
{
	const t_magdyn::MagneticSite* site = nullptr;
};


/**
 * infos for exchange term
 */
struct ExchangeTermInfo
{
	const t_magdyn::ExchangeTerm* term = nullptr;
};


/**
 * magnon calculation dialog
 */
class MagDynDlg : public QDialog
{
public:
	MagDynDlg(QWidget* pParent = nullptr);
	virtual ~MagDynDlg();

	MagDynDlg(const MagDynDlg&) = delete;
	const MagDynDlg& operator=(const MagDynDlg&) = delete;


protected:
	QSettings *m_sett{};
	QMenuBar *m_menu{};
	QSplitter *m_split_inout{};
	QLabel *m_status{};
	QProgressBar *m_progress{};
	QPushButton* m_btnStart{};

	QAction *m_autocalc{};
	QAction *m_use_dmi{};
#ifdef MAGDYN_ALLOW_GENERAL_J
	QAction *m_use_genJ{};
#endif
	QAction *m_use_field{};
	QAction *m_use_temperature{};
	QAction *m_use_weights{};
	QAction *m_use_projector{};
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
	tl2::RecentFiles m_recent{};
	QMenu *m_menuOpenRecent{};
	// function to call for the recent file menu items
	std::function<bool(const QString& filename)> m_open_func
		= [this](const QString& filename) -> bool
	{
		this->Clear();
		this->SetCurrentFile(filename);
		return this->Load(filename);
	};

	QGridLayout *m_maingrid{};

	// tabs
	QTabWidget *m_tabs_in{}, *m_tabs_out{};

	// panels
	QWidget *m_sitespanel{};
	QWidget *m_termspanel{};
	QWidget *m_sampleenviropanel{};
	QWidget *m_varspanel{};
	QWidget *m_notespanel{};
	QWidget *m_disppanel{};
	QWidget *m_hamiltonianpanel{};
	QWidget *m_coordinatespanel{};
	QWidget *m_exportpanel{};

	// sites
	QTableWidget *m_sitestab{};
	QComboBox *m_comboSG{};
	std::vector<std::vector<t_mat_real>> m_SGops{};
	QDoubleSpinBox *m_xtallattice[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_xtalangles[3]{nullptr, nullptr, nullptr};

	// terms, ordering vector, and rotation axis
	QTableWidget *m_termstab{};
	QDoubleSpinBox *m_maxdist{};
	QSpinBox *m_maxSC{};
	QSpinBox *m_maxcouplings{};
	QComboBox *m_comboSG2{};  // copy of m_comboSG
	QDoubleSpinBox *m_ordering[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_normaxis[3]{nullptr, nullptr, nullptr};

	// variables
	QTableWidget *m_varstab{};

	// dispersion
	QCustomPlot *m_plot{};
	std::vector<GraphWithWeights*> m_graphs{};
	QMenu *m_menuDisp{};
	QDoubleSpinBox *m_q_start[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_q_end[3]{nullptr, nullptr, nullptr};
	QSpinBox *m_num_points{};
	QDoubleSpinBox *m_weight_scale{}, *m_weight_min{}, *m_weight_max{};

	// hamiltonian
	QTextEdit *m_hamiltonian{};
	QDoubleSpinBox *m_q[3]{nullptr, nullptr, nullptr};

	// external magnetic field
	QDoubleSpinBox *m_field_dir[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_field_mag{};
	QCheckBox *m_align_spins{};
	QDoubleSpinBox *m_rot_axis[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_rot_angle{};
	QTableWidget *m_fieldstab{};

	// temperature
	QDoubleSpinBox *m_temperature{};

	// notes
	QPlainTextEdit *m_notes{};

	// coordinates
	QTableWidget *m_coordinatestab{};

	// export
	QDoubleSpinBox *m_exportStartQ[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_exportEndQ[3]{nullptr, nullptr, nullptr};
	QSpinBox *m_exportNumPoints[3]{nullptr, nullptr, nullptr};
	QComboBox *m_exportFormat{nullptr};

	// magnon dynamics calculator
	t_magdyn m_dyn{};

	// trafo calculator
	TrafoCalculator *m_trafos{};

	// settings dialog
	QDialog *m_settings_dlg{};

	// structure plotter
	QDialog *m_structplot_dlg{};
	QLabel *m_structplot_status{};
	QCheckBox *m_structplot_coordcross{};
	QCheckBox *m_structplot_labels{};
	QMenu *m_structplot_context{};
	QLabel *m_labelGlInfos[4]{nullptr, nullptr, nullptr, nullptr};
	tl2::GlPlot *m_structplot{};
	std::unordered_map<std::size_t, AtomSiteInfo> m_structplot_atoms{};
	std::unordered_map<std::size_t, ExchangeTermInfo> m_structplot_terms{};
	std::optional<std::size_t> m_structplot_cur_obj{};
	std::optional<std::size_t> m_structplot_cur_atom{};
	std::optional<std::size_t> m_structplot_cur_term{};

	TableImportDlg *m_table_import_dlg{};  // table import dialog
	QDialog *m_info_dlg{};                 // info dialog


protected:
	// set up gui
	void CreateMainWindow();
	void CreateInfoDlg();
	void CreateMenuBar();

	void InitSettingsDlg();
	void InitSettings();

	// set up input panels
	void CreateSitesPanel();
	void CreateExchangeTermsPanel();
	void CreateVariablesPanel();
	void CreateSampleEnvPanel();
	void CreateNotesPanel();

	// set up output panels
	void CreateDispersionPanel();
	void CreateHamiltonPanel();
	void CreateExportPanel();
	void CreateCoordinatesPanel();

	// general table operations
	void MoveTabItemUp(QTableWidget *pTab);
	void MoveTabItemDown(QTableWidget *pTab);
	void ShowTableContextMenu(
		QTableWidget *pTab, QMenu *pMenu, QMenu *pMenuNoItem, const QPoint& pt);
	std::vector<int> GetSelectedRows(
		QTableWidget *pTab, bool sort_reversed = false) const;

	// add a site to the table
	void AddSiteTabItem(int row = -1,
		const std::string& name = "",
		t_real x = 0., t_real y = 0., t_real z = 0.,
		const std::string& sx = "0",
		const std::string& sy = "0",
		const std::string& sz = "1",
		t_real S = 1.,
		const std::string& sox = "auto",
		const std::string& soy = "auto",
		const std::string& soz = "auto",
		const std::string& rgb = "auto");

	// add a coupling to the table
	void AddTermTabItem(int row = -1,
		const std::string& name = "",
		t_size atom_1 = 0, t_size atom_2 = 0,
		t_real dist_x = 0., t_real dist_y = 0., t_real dist_z = 0.,
		const std::string& J = "0",
		const std::string& dmi_x = "0", const std::string& dmi_y = "0", const std::string& dmi_z = "0",
		const std::string& gen_xx = "0", const std::string& gen_xy = "0", const std::string& gen_xz = "0",
		const std::string& gen_yx = "0", const std::string& gen_yy = "0", const std::string& gen_yz = "0",
		const std::string& gen_zx = "0", const std::string& gen_zy = "0", const std::string& gen_zz = "0",
		const std::string& rgb = "#0x00bf00");

	// add a variable to the table
	void AddVariableTabItem(int row = -1,
		const std::string& name = "var",
		const t_cplx& var = t_cplx{0, 0});

	// add a coordinate to the table
	void AddCoordinateTabItem(int row = -1,
		t_real hi = 0., t_real ki = 0., t_real li = 1.,
		t_real hf = 0., t_real kf = 0., t_real lf = 1.);

	// add a magnetic field to the table
	void AddFieldTabItem(int row = -1,
		t_real Bh = 0., t_real Bk = 0., t_real Bl = 1.,
		t_real Bmag = 1.);

	void SetCurrentField();
	void SetCurrentCoordinate(int which = 0);

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
	void ExportSQE();

	void SavePlotFigure();
	void SaveDispersion();

	void MirrorAtoms();
	void RotateField(bool ccw = true);
	void GenerateSitesFromSG();
	void GenerateCouplingsFromSG();
	void GeneratePossibleCouplings();

	std::optional<t_size> GetTermAtomIndex(int row, int num) const;
	void SyncSitesAndTerms();
	void CalcAll();          // syncs sites and terms and calculates all dynamics

	void PlotDispersion();

	void PlotMouseMove(QMouseEvent* evt);
	void PlotMousePress(QMouseEvent* evt);

	virtual void mousePressEvent(QMouseEvent *evt) override;
	virtual void closeEvent(QCloseEvent *evt) override;
	virtual void dragEnterEvent(QDragEnterEvent *evt) override;
	virtual void dropEvent(QDropEvent *evt) override;

	// table importer
	void ShowTableImporter();
	void ImportAtoms(const std::vector<TableImportAtom>&);
	void ImportCouplings(const std::vector<TableImportCoupling>&);

	// structure plotter
	void ShowStructurePlot();
	void StructPlotMouseClick(bool left, bool mid, bool right);
	void StructPlotMouseDown(bool left, bool mid, bool right);
	void StructPlotMouseUp(bool left, bool mid, bool right);
	void StructPlotPickerIntersection(
		const t_vec3_gl* pos, std::size_t objIdx,
		const t_vec3_gl* posSphere);
	void StructPlotAfterGLInitialisation();
	void StructPlotSync();
	void StructPlotDelete();
	void StructPlotShowCoordCross(bool show);
	void StructPlotShowLabels(bool show);
	void StructPlotCentreCamera();

	// disable/enable gui input for threaded operations
	void EnableInput();
	void DisableInput();


public:
	bool Load(const QString& filename, bool calc_dynamics = true);
	bool Save(const QString& filename);

	bool ExportSQE(const QString& filename);

	void SetCurrentFileAndDir(const QString& filename);
	void SetCurrentFile(const QString& filename);

	void SetCoordinates(const t_vec_real& Qi, const t_vec_real& Qf, bool calc_dynamics = true);

	void CalcDispersion();
	void CalcHamiltonian();


private:
	int m_sites_cursor_row = -1;
	int m_terms_cursor_row = -1;
	int m_variables_cursor_row = -1;
	int m_fields_cursor_row = -1;
	int m_coordinates_cursor_row = -1;

	// site/coupling counter for newly inserted items
	std::size_t m_curSiteCtr = 0;
	std::size_t m_curCouplingCtr = 0;

	bool m_ignoreTableChanges = true;
	bool m_ignoreCalc = false;
	bool m_stopRequested = false;

	// data for dispersion plot
	QVector<qreal> m_qs_data{}, m_Es_data{}, m_ws_data{};
	QVector<qreal> m_qs_data_channel[3]{}, m_Es_data_channel[3]{}, m_ws_data_channel[3]{};
	t_size m_Q_idx{};                 // plot x axis
	t_real m_Q_start{}, m_Q_end{};    // plot x axis range

	// reference object handles
	std::size_t m_structplot_sphere = 0;
	std::size_t m_structplot_arrow = 0;
	std::size_t m_structplot_cyl = 0;
};


#endif
