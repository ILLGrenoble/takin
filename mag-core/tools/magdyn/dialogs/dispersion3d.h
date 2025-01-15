/**
 * magnetic dynamics -- 3d dispersion plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date January 2025
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2025  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MAG_DYN_DISP3D_DLG_H__
#define __MAG_DYN_DISP3D_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>

#include <qcustomplot.h>
#include <vector>
#include <array>
#include <tuple>
#include <unordered_map>

#include "tlibs2/libs/qt/glplot.h"
#include "gui_defs.h"



/**
 * 3d dispersion dialog
 */
class Dispersion3DDlg : public QDialog
{ Q_OBJECT
protected:
	using t_data_Q = std::tuple<t_vec_real /*0: Q*/, t_real /*1: E*/, t_real /*2: S*/,
		t_size /*3: Q_idx_1*/, t_size /*4: Q_idx_2*/,
		t_size /*5: degeneracy*/, bool /*6: valid*/>;
	using t_data_Qs = std::vector<t_data_Q>;
	using t_data_bands = std::vector<t_data_Qs>;


public:
	Dispersion3DDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~Dispersion3DDlg();

	Dispersion3DDlg(const Dispersion3DDlg&) = delete;
	Dispersion3DDlg& operator=(const Dispersion3DDlg&) = delete;

	// set kernel and Q path from the main window
	void SetKernel(const t_magdyn* dyn);
	void SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend);


protected:
	virtual void accept() override;

	// calculation functions
	void EnableCalculation(bool enable = true);
	void Calculate();
	void Plot(bool clear_settings = true);

	// calculation helper functions
	t_real GetMeanEnergy(const t_data_Qs& data) const;
	t_real GetMeanEnergy(t_size band_idx) const;
	std::tuple<t_vec_real, t_vec_real, t_vec_real> GetQVectors() const;
	std::array<int, 3> GetBranchColour(t_size branch_idx, t_size num_branches) const;
	void FromMainQ();
	void ShowError(const char* msg);

	// band table functions
	void ClearBands();
	void AddBand(const std::string& name, const QColor& colour, bool enabled = true);
	bool IsBandEnabled(t_size idx) const;

	// export data
	void WriteHeader(std::ostream& ostr) const;
	void SaveData();
	void SaveScript();

	// ------------------------------------------------------------------------
	// plotter interface
	void AfterPlotGLInitialisation();
	void PlotPickerIntersection(const t_vec3_gl* pos,
		std::size_t objIdx, std::size_t triagIdx,
		const t_vec3_gl* posSphere);

	void PlotCameraHasUpdated();
	void CentrePlotCamera();
	void CentrePlotCameraOnObject();

	void PlotMouseClick(bool left, bool mid, bool right);
	void PlotMouseDown(bool left, bool mid, bool right);
	void PlotMouseUp(bool left, bool mid, bool right);

	void SetPlotCoordinateSystem(int which);
	void ShowPlotCoordCross(bool show);
	void ShowPlotLabels(bool show);
	void SetPlotPerspectiveProjection(bool proj);
	void SetPlotCameraRotation(t_real_gl phi, t_real_gl theta);
	// ------------------------------------------------------------------------


private:
	t_size m_Q_count_1{}, m_Q_count_2{}; // number of Q points along the two directions
	t_data_bands m_data{};               // data for all energy bands
	t_real m_max_E{};                    // maximum band energy

	// ------------------------------------------------------------------------
	// from main dialog
	const t_magdyn *m_dyn{};             // main calculation kernel
	t_vec_real m_Qstart{}, m_Qend{};     // Qs from the main window
	QSettings *m_sett{};                 // program settings
	// ------------------------------------------------------------------------

	// dispersion
	tl2::GlPlot *m_dispplot{};           // 3d plotter
	std::optional<std::size_t> m_cur_obj{};
	std::unordered_map<std::size_t /*plot object*/, t_size /*band index*/> m_band_objs{};
	t_vec_gl m_cam_centre{tl2::zero<t_vec_gl>(3)};

	QSplitter *m_split_plot{};
	QTableWidget *m_table_bands{};       // table listing the magnon bands
	QDoubleSpinBox *m_Q_origin[3]{}, *m_Q_dir1[3]{}, *m_Q_dir2[3]{};
	QSpinBox *m_num_Q_points[2]{};
	QCheckBox *m_only_pos_E{};           // ignore magnon annihilation
	QMenu *m_context{};                  // general plot context menu
	QMenu *m_context_band{};             // context menu for magnon bands

	// plot options
	QDoubleSpinBox *m_Q_scale1{}, *m_Q_scale2{}, *m_E_scale{};
	QDoubleSpinBox *m_cam_phi{}, *m_cam_theta{};
	QCheckBox *m_perspective{};

	QPushButton *m_btn_start_stop{};     // start/stop calculation
	QProgressBar *m_progress{};          // progress bar
	QLabel *m_status{};                  // status bar

	bool m_calc_enabled{};               // enable calculations
	bool m_stop_requested{};             // stop running calculations


signals:
	void GlDeviceInfos(const std::string& ver, const std::string& shader_ver,
		const std::string& vendor, const std::string& renderer);
};


#endif
