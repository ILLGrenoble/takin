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
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QLabel>

#include <qcustomplot.h>
#include <vector>

#include "tlibs2/libs/qt/glplot.h"
#include "gui_defs.h"



/**
 * 3d dispersion dialog
 */
class Dispersion3DDlg : public QDialog
{ Q_OBJECT
public:
	Dispersion3DDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~Dispersion3DDlg();

	Dispersion3DDlg(const Dispersion3DDlg&) = delete;
	Dispersion3DDlg& operator=(const Dispersion3DDlg&) = delete;

	// set kernel from main window
	void SetKernel(const t_magdyn* dyn);


protected:
	virtual void accept() override;

	void ShowError(const char* msg);
	void EnableCalculation(bool enable = true);
	void Calculate();
	void Plot();

	// ------------------------------------------------------------------------
	// plotter interface
	void AfterPlotGLInitialisation();
	void PlotPickerIntersection(const t_vec3_gl* pos,
		std::size_t objIdx, std::size_t triagIdx,
		const t_vec3_gl* posSphere);

	void PlotCameraHasUpdated();

	void PlotMouseClick(bool left, bool mid, bool right);
	void PlotMouseDown(bool left, bool mid, bool right);
	void PlotMouseUp(bool left, bool mid, bool right);
	// ------------------------------------------------------------------------


private:
	// ------------------------------------------------------------------------
	// from main dialog
	const t_magdyn *m_dyn{};            // main calculation kernel
	QSettings *m_sett{};                // program settings
	// ------------------------------------------------------------------------

	tl2::GlPlot *m_dispplot{};          // 3d plotter

	QDoubleSpinBox *m_Q_origin[3]{}, *m_Q_dir1[3]{}, *m_Q_dir2[3]{};
	QSpinBox *m_num_Q_points[2]{};
	QPushButton *m_btn_start_stop{};    // start/stop calculation
	QProgressBar *m_progress{};         // progress bar
	QLabel *m_status{};                 // status bar

	bool m_calc_enabled{};              // enable calculations
	bool m_stop_requested{};            // stop running calculations


signals:
	void GlDeviceInfos(const std::string& ver, const std::string& shader_ver,
		const std::string& vendor, const std::string& renderer);
};


#endif
