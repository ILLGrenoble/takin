/**
 * atom dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date Dec-2019
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __MOLDYN_GUI_H__
#define __MOLDYN_GUI_H__

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QLabel>
#include <QtWidgets/QDoubleSpinBox>
#include <QtCore/QSettings>

#include <vector>
#include <tuple>

#include "tlibs2/libs/qt/glplot.h"
#include "tlibs2/libs/maths.h"
#include "libs/defs.h"

#include "moldyn-loader.h"


using t_vec = tl2::vec<t_real, std::vector>;
using t_mat = tl2::mat<t_real, std::vector>;

using t_real_gl = tl2::t_real_gl;
using t_vec2_gl = tl2::t_vec2_gl;
using t_vec3_gl = tl2::t_vec3_gl;
using t_vec_gl = tl2::t_vec_gl;
using t_mat_gl = tl2::t_mat_gl;


/**
 * atom indices for convex hull calculation
 */
struct Hull
{
	std::vector<std::tuple<std::size_t, std::size_t>> vertices;
	std::optional<std::size_t> plotObj;
};


class MolDynDlg : public QMainWindow
{
public:
	MolDynDlg(QWidget* pParent = nullptr);
	~MolDynDlg() = default;

protected:
	std::size_t Add3DAtom(const t_vec& vec, const t_vec& col, t_real scale, const std::string& typelabel, int atomindex=-1);
	void Change3DAtom(std::size_t obj, const t_vec* vec, const t_vec* col=nullptr, const t_real *scale=nullptr,
		const std::string *typelabel=nullptr, int atomindex=-1);

	void SetStatusMsg(const std::string& msg);
	void UpdateAtomsStatusMsg();

	void New();
	void Load();
	void SaveAs();

	void PlotMouseDown(bool left, bool mid, bool right);
	void PlotMouseUp(bool left, bool mid, bool right);
	void PlotMouseClick(bool left, bool mid, bool right);
	void PickerIntersection(const t_vec3_gl* pos,
		std::size_t objIdx, std::size_t triagIdx,
		const t_vec3_gl* posSphere);
	void AfterGLInitialisation();
	void GLInitialisationFailed();

	std::tuple<bool, std::size_t, std::size_t, std::size_t> GetAtomIndexFromHandle(std::size_t handle) const;
	std::tuple<std::string, int> SplitDataString(const std::string&) const;

	void CalculateDistanceBetweenAtoms();
	void CalculatePositionsOfAtoms();
	void CalculateDeltaDistancesOfAtoms();
	void CalculateConvexHullOfAtoms();

	void CalculateConvexHulls();

	void SliderValueChanged(int val);

	void SelectAll();
	void SelectNone();

	std::vector<std::tuple<std::size_t, std::size_t>> GetSelectedAtoms();

    void DeleteSelectedAtoms();
    void OnlyKeepSelectedAtoms();
	void SelectAtomsOfSameType();
	void DeleteAtomUnderCursor();
	void DeleteAllAtomsOfSameType();
	void KeepAtomsOfSameType();

	virtual void closeEvent(QCloseEvent *evt) override;
	virtual void keyPressEvent(QKeyEvent *evt) override;


protected:
	MolDyn<t_real, t_vec> m_mol;
	t_mat m_crystA = tl2::unit<t_mat>(3);
	t_mat m_crystB = tl2::unit<t_mat>(3);

	QSettings *m_sett = nullptr;
	QMenuBar *m_menu = nullptr;

	QStatusBar *m_status = nullptr;
	QLabel *m_statusCurAtom = nullptr;
	QLabel *m_statusAtoms = nullptr;

	QSlider *m_sliderFrame = nullptr;
	QDoubleSpinBox *m_spinScale = nullptr;
	QMenu *m_atomContextMenu = nullptr;

	tl2::GlPlot *m_plot = nullptr;
	std::size_t m_sphere = 0;
	std::vector<std::size_t> m_sphereHandles;

	std::vector<Hull> m_hulls;


private:
	long m_curPickedObj = -1;
	bool m_ignoreChanges = 1;
};


#endif
