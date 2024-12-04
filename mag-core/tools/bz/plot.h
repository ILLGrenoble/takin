/**
 * brillouin zone tool -- 3d plot
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
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

#ifndef __BZTOOL_PLOTTER_H__
#define __BZTOOL_PLOTTER_H__

#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QCheckBox>
#include <QtCore/QSettings>

#include <memory>
#include <vector>
#include <string>

#include "globals.h"

#include "tlibs2/libs/qt/glplot.h"


class BZPlotDlg : public QDialog
{ Q_OBJECT
public:
	BZPlotDlg(QWidget* pParent = nullptr, QSettings *sett = nullptr,
		QLabel **infos = nullptr);
	~BZPlotDlg() = default;

	void Clear();

	void SetABTrafo(const t_mat& crystA, const t_mat& crystB);
	void AddVoronoiVertex(const t_vec& pos);
	void AddBraggPeak(const t_vec& pos);
	void AddTriangles(const std::vector<t_vec>& vecs);
	void SetPlane(const t_vec& norm, t_real d);


protected:
	void SetStatusMsg(const std::string& msg);
	void ShowCoordCross(bool show);
	void ShowLabels(bool show);
	void ShowPlane(bool show);
	void ShowQVertices(bool show);

	void PlotMouseDown(bool left, bool mid, bool right);
	void PlotMouseUp(bool left, bool mid, bool right);
	void PickerIntersection(const t_vec3_gl* pos,
		std::size_t objIdx, const t_vec3_gl* posSphere);
	void AfterGLInitialisation();

	virtual void closeEvent(QCloseEvent *evt) override;


private:
	t_mat m_crystA = tl2::unit<t_mat>(3);   // crystal A matrix
	t_mat m_crystB = tl2::unit<t_mat>(3);   // crystal B matrix

	QSettings *m_sett = nullptr;

	std::shared_ptr<tl2::GlPlot> m_plot;
	std::size_t m_sphere = 0, m_plane = 0;

	QLabel **m_labelGlInfos = nullptr;
	QLabel *m_status = nullptr;
	QCheckBox *m_show_coordcross = nullptr;
	QCheckBox *m_show_labels = nullptr;
	QCheckBox *m_show_plane = nullptr;
	QCheckBox *m_show_Qs = nullptr;

	long m_curPickedObj = -1;                 // current 3d bz object
	std::vector<std::size_t> m_objsBragg{};   // Bragg peak plot objects
	std::vector<std::size_t> m_objsVoronoi{}; // Voronoi vertex plot objects
	std::vector<std::size_t> m_objsBZ{};      // BZ triangle plot objects


signals:
	void NeedRecalc();
};


#endif
