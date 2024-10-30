/**
 * magnetic dynamics -- magnetic structure plotting
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
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

#ifndef __MAG_DYN_STRUCTPLOT_H__
#define __MAG_DYN_STRUCTPLOT_H__

#include <QtCore/QSettings>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QDialog>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>

#include "tlibs2/libs/qt/glplot.h"

#include "gui_defs.h"
#include "infos.h"



/**
 * infos for magnetic sites
 */
struct MagneticSiteInfo
{
	std::string name{};
	t_vec_real uc_pos{};  // position in unit cell
	t_vec_real sc_pos{};  // position in super cell
};



/**
 * infos for exchange term
 */
struct ExchangeTermInfo
{
	std::string name{};
	t_real length{};
};



/**
 * magnon calculation dialog
 */
class StructPlotDlg : public QDialog
{ Q_OBJECT
public:
	StructPlotDlg(QWidget *parent, QSettings *sett = nullptr, InfoDlg *info = nullptr);
	~StructPlotDlg() = default;
	StructPlotDlg(const StructPlotDlg&) = delete;
	StructPlotDlg& operator=(const StructPlotDlg&) = delete;

	void Sync();
	void SetKernel(const t_magdyn* dyn) { m_dyn = dyn; }
	void SetTables(QTableWidget *sites, QTableWidget *terms)
	{
		m_sitestab = sites;
		m_termstab = terms;
	}


protected:
	void MouseClick(bool left, bool mid, bool right);
	void MouseDown(bool left, bool mid, bool right);
	void MouseUp(bool left, bool mid, bool right);

	void AfterGLInitialisation();

	void PickerIntersection(
		const t_vec3_gl* pos, std::size_t objIdx,
		const t_vec3_gl* posSphere);

	void DeleteItem();
	void FlipSpin();

	void ShowCoordCross(bool show);
	void ShowLabels(bool show);
	void SetPerspectiveProjection(bool proj);
	void SetCameraRotation(t_real_gl phi, t_real_gl theta);
	void SetCoordinateSystem(int which);

	void CameraHasUpdated();
	void CentreCamera();
	void CentreCameraOnObject();

	void HighlightSite(const std::string& name);
	void HighlightTerm(const std::string& name);

	virtual void accept() override;


private:
	const t_magdyn *m_dyn{};

	// connections to main dialog
	QSettings *m_sett{};
	InfoDlg *m_info_dlg{};
	QTableWidget *m_sitestab{}, *m_termstab{};

	QLabel *m_status{};
	QCheckBox *m_coordcross{}, *m_labels{}, *m_perspective{};
	QComboBox *m_coordsys{};
	QDoubleSpinBox *m_cam_phi{}, *m_cam_theta{};
	QMenu *m_context{}, *m_context_site{}, *m_context_term{};

	tl2::GlPlot *m_structplot{};
	std::unordered_map<std::size_t, MagneticSiteInfo> m_sites{};
	std::unordered_map<std::size_t, ExchangeTermInfo> m_terms{};
	std::optional<std::size_t> m_cur_obj{};
	std::optional<std::string> m_cur_site{};
	std::optional<std::string> m_cur_term{};
	t_vec_gl m_centre{tl2::zero<t_vec_gl>(3)};

	// reference object handles
	std::size_t m_sphere{};
	std::size_t m_arrow{};
	std::size_t m_cyl{};


signals:
	void SelectSite(const std::string& site);
	void DeleteSite(const std::string& site);
	void FlipSiteSpin(const std::string& site);

	void SelectTerm(const std::string& term);
	void DeleteTerm(const std::string& term);
};


#endif
