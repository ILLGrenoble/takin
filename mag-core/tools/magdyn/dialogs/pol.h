/**
 * Calculation of polarisation vector
 * @author Tobias Weber <tweber@ill.fr>
 * @date Oct-2018
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "magtools" project
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
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

#ifndef __TAKIN_POL__
#define __TAKIN_POL__

#include <QtCore/QSettings>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>

#include <vector>
#include <optional>

#include "tlibs2/libs/qt/glplot.h"
#include "libs/defs.h"


class PolDlg : public QDialog
{ Q_OBJECT
public:
	// create UI
	PolDlg(QWidget* pParent, QSettings *sett = nullptr);

	PolDlg() = delete;
	PolDlg(const PolDlg&) = delete;
	PolDlg& operator=(const PolDlg&) = delete;

	// calculate final polarisation vector
	void CalcPol();


protected:
	virtual void closeEvent(QCloseEvent *) override;


protected slots:
	//called after the plotter has initialised
	void AfterGLInitialisation();

	// get the length of a vector
	t_real GetArrowLen(std::size_t objIdx) const;

	// called when the mouse hovers over an object
	void PickerIntersection(const tl2::t_vec3_gl* pos,
		std::size_t objIdx, std::size_t /*triagIdx*/,
		const tl2::t_vec3_gl* posSphere);

	// mouse button pressed
	void MouseDown(bool left, bool mid, bool right);

	// mouse button released
	void MouseUp(bool left, bool mid, bool right);


private:
	QSettings *m_sett{};

	std::shared_ptr<tl2::GlPlot> m_plot{std::make_shared<tl2::GlPlot>(this)};

	QLineEdit* m_editNRe = new QLineEdit("0", this);
	QLineEdit* m_editNIm = new QLineEdit("0", this);

	QLineEdit* m_editMPerpReX = new QLineEdit("0", this);
	QLineEdit* m_editMPerpReY = new QLineEdit("1", this);
	QLineEdit* m_editMPerpReZ = new QLineEdit("0", this);
	QLineEdit* m_editMPerpImX = new QLineEdit("0", this);
	QLineEdit* m_editMPerpImY = new QLineEdit("0", this);
	QLineEdit* m_editMPerpImZ = new QLineEdit("0", this);

	QLineEdit* m_editPiX = new QLineEdit("1", this);
	QLineEdit* m_editPiY = new QLineEdit("0", this);
	QLineEdit* m_editPiZ = new QLineEdit("0", this);
	QLineEdit* m_editPfX = new QLineEdit(this);
	QLineEdit* m_editPfY = new QLineEdit(this);
	QLineEdit* m_editPfZ = new QLineEdit(this);

	QLabel* m_labelStatus = new QLabel(this);

	// 3d object handles
	std::size_t m_arrow_pi = 0;
	std::size_t m_arrow_pf = 0;
	std::size_t m_arrow_M_Re = 0;
	std::size_t m_arrow_M_Im = 0;

	bool m_3dobjsReady = false;
	bool m_mouseDown[3] = { false, false, false };
	std::optional<std::size_t> m_curPickedObj{};
	std::optional<std::size_t> m_curDraggedObj{};


signals:
	void GlDeviceInfos(const std::string& ver, const std::string& shader_ver,
		const std::string& vendor, const std::string& renderer);
};


#endif
