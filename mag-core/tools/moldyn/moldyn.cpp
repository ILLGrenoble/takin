/**
 * atom dynamics tool
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

#include "moldyn.h"

#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QProgressDialog>

#include <iostream>
#include <tuple>
#include <memory>

#include "tlibs2/libs/algos.h"
#include "tlibs2/libs/qt/helper.h"

using namespace tl2_ops;

constexpr t_real g_eps = 1e-6;
constexpr int g_prec = 6;


#define PROG_NAME "Molecular Dynamics Tool"
#define PROG_VER "0.2"

#define DATA_SEP "#|#"	// data separation string



// ----------------------------------------------------------------------------
/**
 * File dialog with options
 */
class MolDynFileDlg : public QFileDialog
{
	public:
		MolDynFileDlg(QWidget *parent, const QString& title, const QString& dir, const QString& filter)
			: QFileDialog(parent, title, dir, filter)
		{
			// options panel with frame skip
			QLabel *labelFrameSkip = new QLabel("Frame Skip: ", this);
			m_spinFrameSkip = new QSpinBox(this);
			m_spinFrameSkip->setValue(10);
			m_spinFrameSkip->setSingleStep(1);
			m_spinFrameSkip->setRange(0, 9999999);

			labelFrameSkip->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
			m_spinFrameSkip->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});

			QWidget *pPanel = new QWidget();
			auto pPanelGrid = new QGridLayout(pPanel);
			pPanelGrid->setSpacing(2);
			pPanelGrid->setContentsMargins(4,4,4,4);

			pPanelGrid->addWidget(labelFrameSkip, 0,0,1,1);
			pPanelGrid->addWidget(m_spinFrameSkip, 0,1,1,1);

			// add the options panel to the layout
			setOptions(QFileDialog::DontUseNativeDialog);
			QGridLayout *pGrid = reinterpret_cast<QGridLayout*>(layout());
			if(pGrid)
				pGrid->addWidget(pPanel, pGrid->rowCount(), 0, 1, pGrid->columnCount());
		}


		int GetFrameSkip() const
		{
			return m_spinFrameSkip->value();
		}


	private:
		QSpinBox *m_spinFrameSkip = nullptr;
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
MolDynDlg::MolDynDlg(QWidget* pParent) : QMainWindow{pParent},
	m_sett{new QSettings{"takin", "moldyn", this}}
{
	setWindowTitle(QString(PROG_NAME) + QString(" Version ") + QString(PROG_VER));
	this->setObjectName("moldyn");

	m_status = new QStatusBar(this);
	m_statusCurAtom = new QLabel(m_status);
	m_statusAtoms = new QLabel(m_status);
	m_status->addWidget(m_statusCurAtom);
	m_status->addPermanentWidget(m_statusAtoms);
	this->setStatusBar(m_status);


	QWidget *pMainPanel = new QWidget();
	auto pMainGrid = new QGridLayout(pMainPanel);
	pMainGrid->setSpacing(2);
	pMainGrid->setContentsMargins(4,4,4,4);
	this->setCentralWidget(pMainPanel);


	// menu bar
	{
		m_menu = new QMenuBar(this);
		m_menu->setNativeMenuBar(m_sett ? m_sett->value("native_gui", false).toBool() : false);


		// File
		auto menuFile = new QMenu("File", m_menu);

		auto acNew = new QAction("New", menuFile);
		auto acLoad = new QAction("Load...", menuFile);
		auto acSaveAs = new QAction("Save As...", menuFile);
		auto acExit = new QAction("Quit", menuFile);

		menuFile->addAction(acNew);
		menuFile->addSeparator();
		menuFile->addAction(acLoad);
		menuFile->addAction(acSaveAs);
		menuFile->addSeparator();
		menuFile->addAction(acExit);

		connect(acNew, &QAction::triggered, this, &MolDynDlg::New);
		connect(acLoad, &QAction::triggered, this, &MolDynDlg::Load);
		connect(acSaveAs, &QAction::triggered, this, &MolDynDlg::SaveAs);
		connect(acExit, &QAction::triggered, this, &QDialog::close);


		// Edit
		auto menuEdit = new QMenu("Edit", m_menu);
		auto acSelectAll = new QAction("Select All", menuEdit);
		auto acSelectNone = new QAction("Select None", menuEdit);

		menuEdit->addAction(acSelectAll);
		menuEdit->addAction(acSelectNone);

		connect(acSelectAll, &QAction::triggered, this, &MolDynDlg::SelectAll);
		connect(acSelectNone, &QAction::triggered, this, &MolDynDlg::SelectNone);


		// Calculate
		auto menuCalc = new QMenu("Calculate", m_menu);
		auto acCalcDist = new QAction("Distance Between Selected Atoms...", menuEdit);
		auto acCalcPos = new QAction("Positions Of Selected Atoms...", menuEdit);
		auto acCalcDeltaDist = new QAction("Distances to Initial Position of Selected Atoms...", menuEdit);
		auto acCalcHull = new QAction("Convex Hull of Selected Atoms", menuEdit);

		menuCalc->addAction(acCalcDist);
		menuCalc->addAction(acCalcPos);
		menuCalc->addAction(acCalcDeltaDist);
#ifdef USE_QHULL
		menuCalc->addSeparator();
		menuCalc->addAction(acCalcHull);
#endif

		connect(acCalcDist, &QAction::triggered, this, &MolDynDlg::CalculateDistanceBetweenAtoms);
		connect(acCalcPos, &QAction::triggered, this, &MolDynDlg::CalculatePositionsOfAtoms);
		connect(acCalcDeltaDist, &QAction::triggered, this, &MolDynDlg::CalculateDeltaDistancesOfAtoms);
		connect(acCalcHull, &QAction::triggered, this, &MolDynDlg::CalculateConvexHullOfAtoms);


		// Help
		auto menuHelp = new QMenu("Help", m_menu);
		auto acHelpInfo = new QAction("Infos...", menuHelp);

		menuHelp->addAction(acHelpInfo);

		connect(acHelpInfo, &QAction::triggered, this, [this]()
		{
			QString strHelp;

			strHelp = QString{PROG_NAME} + " v" + QString{PROG_VER} + ", part of the Takin 2 package.\n";
			strHelp += "Written by Tobias Weber <tweber@ill.fr>\nin December 2019.\n\n";

			strHelp += "This program is free software: you can redistribute it and/or modify "
				"it under the terms of the GNU General Public License as published by "
				"the Free Software Foundation, version 3 of the License.\n\n"
				"This program is distributed in the hope that it will be useful, "
				"but WITHOUT ANY WARRANTY; without even the implied warranty of "
				"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the "
				"GNU General Public License for more details.\n\n"
				"You should have received a copy of the GNU General Public License "
				"along with this program.  If not, see <http://www.gnu.org/licenses/>.";

			QMessageBox::information(this, PROG_NAME, strHelp);
		});


		m_menu->addMenu(menuFile);
		m_menu->addMenu(menuEdit);
		m_menu->addMenu(menuCalc);
		m_menu->addMenu(menuHelp);
		this->setMenuBar(m_menu);
	}


	// context menus
	{
		m_atomContextMenu = new QMenu(this);
		m_atomContextMenu->setTitle("Atoms");
		m_atomContextMenu->addAction("Delete Atom", this, &MolDynDlg::DeleteAtomUnderCursor);
        m_atomContextMenu->addSeparator();
        m_atomContextMenu->addAction("Delete Selected Atoms", this, &MolDynDlg::DeleteSelectedAtoms);
        m_atomContextMenu->addAction("Delete All But Selected Atoms", this, &MolDynDlg::OnlyKeepSelectedAtoms);
        m_atomContextMenu->addSeparator();
		m_atomContextMenu->addAction("Delete All Atoms Of Selected Type", this, &MolDynDlg::DeleteAllAtomsOfSameType);
		m_atomContextMenu->addAction("Only Keep Atoms Of Selected Type", this, &MolDynDlg::KeepAtomsOfSameType);
		m_atomContextMenu->addSeparator();
		m_atomContextMenu->addAction("Select All Atoms of Same Type", this, &MolDynDlg::SelectAtomsOfSameType);
	}


	// plot widget
	{
		m_plot = new tl2::GlPlot(this);
		m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

		m_plot->GetRenderer()->EnablePicker(1);
		m_plot->GetRenderer()->SetLight(0, tl2::create<t_vec3_gl>({ 5, 5, 5 }));
		m_plot->GetRenderer()->SetLight(1, tl2::create<t_vec3_gl>({ -5, -5, -5 }));
		m_plot->GetRenderer()->SetCoordMax(1.);
		m_plot->GetRenderer()->GetCamera().SetDist(5.);
		m_plot->GetRenderer()->GetCamera().UpdateTransformation();

		connect(m_plot, &tl2::GlPlot::AfterGLInitialisation, this, &MolDynDlg::AfterGLInitialisation);
		connect(m_plot, &tl2::GlPlot::GLInitialisationFailed, this, &MolDynDlg::GLInitialisationFailed);
		connect(m_plot->GetRenderer(), &tl2::GlPlotRenderer::PickerIntersection, this, &MolDynDlg::PickerIntersection);
		connect(m_plot, &tl2::GlPlot::MouseDown, this, &MolDynDlg::PlotMouseDown);
		connect(m_plot, &tl2::GlPlot::MouseUp, this, &MolDynDlg::PlotMouseUp);
		connect(m_plot, &tl2::GlPlot::MouseClick, this, &MolDynDlg::PlotMouseClick);

		//this->setCentralWidget(m_plot);
		pMainGrid->addWidget(m_plot, 0,0,1,9);
	}


	// controls
	{
		auto labCoordSys = new QLabel("Coordinates:", this);
		auto labFrames = new QLabel("Frames:", this);
		auto labScale = new QLabel("Scale:", this);
		labCoordSys->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
		labFrames->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});
		labScale->setSizePolicy(QSizePolicy{QSizePolicy::Fixed, QSizePolicy::Fixed});

		auto comboCoordSys = new QComboBox(this);
		comboCoordSys->addItem("Fractional Units (rlu)");
		comboCoordSys->addItem("Lab Units (\xe2\x84\xab)");
		comboCoordSys->setFocusPolicy(Qt::StrongFocus);

		m_spinScale = new QDoubleSpinBox(this);
		m_spinScale->setDecimals(4);
		m_spinScale->setRange(1e-4, 1e4);
		m_spinScale->setSingleStep(0.1);
		m_spinScale->setValue(0.4);
		m_spinScale->setFocusPolicy(Qt::StrongFocus);

		m_sliderFrame = new QSlider(Qt::Horizontal, this);
		m_sliderFrame->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Minimum});
		m_sliderFrame->setMinimum(0);
		m_sliderFrame->setSingleStep(1);
		m_sliderFrame->setPageStep(10);
		m_sliderFrame->setTracking(1);
		m_sliderFrame->setFocusPolicy(Qt::StrongFocus);

		connect(m_sliderFrame, &QSlider::valueChanged, this, &MolDynDlg::SliderValueChanged);
		connect(comboCoordSys, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), this, [this](int val)
		{
			if(this->m_plot)
				this->m_plot->GetRenderer()->SetCoordSys(val);
		});
		connect(m_spinScale, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](double val)
		{
			if(!this->m_plot) return;

			// hack to trigger update
			SliderValueChanged(m_sliderFrame->value());
		});

		pMainGrid->addWidget(labCoordSys, 1,0,1,1);
		pMainGrid->addWidget(comboCoordSys, 1,1,1,1);
		pMainGrid->addWidget(labScale, 1,2,1,1);
		pMainGrid->addWidget(m_spinScale, 1,3,1,1);
		pMainGrid->addWidget(labFrames, 1,4,1,1);
		pMainGrid->addWidget(m_sliderFrame, 1,5,1,4);
	}


	// restore window size and position
	if(m_sett && m_sett->contains("geo"))
		restoreGeometry(m_sett->value("geo").toByteArray());
	else
		resize(600, 500);

	m_ignoreChanges = 0;
}



// ----------------------------------------------------------------------------
/**
 * add an atom
 */
std::size_t MolDynDlg::Add3DAtom(const t_vec& vec, const t_vec& col, t_real scale, const std::string& typelabel, int atomindex)
{
	auto obj = m_plot->GetRenderer()->AddLinkedObject(m_sphere, 0,0,0, col[0],col[1],col[2],1);
	Change3DAtom(obj, &vec, &col, &scale, &typelabel, atomindex);
	return obj;
}


/**
 * change an atom
 */
void MolDynDlg::Change3DAtom(std::size_t obj, const t_vec *vec, const t_vec *col, const t_real *scale,
	const std::string *label, int atomindex)
{
	if(vec)
	{
		t_mat_gl mat = tl2::hom_translation<t_mat_gl>(
			t_real_gl((*vec)[0]), t_real_gl((*vec)[1]), t_real_gl((*vec)[2]));
		if(scale) mat *= tl2::hom_scaling<t_mat_gl>(
			t_real_gl(*scale), t_real_gl(*scale), t_real_gl(*scale));
		m_plot->GetRenderer()->SetObjectMatrix(obj, mat);
	}

	if(col) m_plot->GetRenderer()->SetObjectCol(obj, (*col)[0], (*col)[1], (*col)[2], 1);
	if(label)
	{
		std::ostringstream ostrData;
		ostrData << *label << DATA_SEP << atomindex;

		m_plot->GetRenderer()->SetObjectLabel(obj, *label);
		m_plot->GetRenderer()->SetObjectDataString(obj, ostrData.str());
	}
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------

std::vector<std::tuple<std::size_t, std::size_t>> MolDynDlg::GetSelectedAtoms()
{
	// get selected atoms
	std::vector<std::tuple<std::size_t, std::size_t>> objs;

	for(const auto& obj : m_sphereHandles)
	{
		// continue if object isn't selected
		if(!m_plot->GetRenderer()->GetObjectHighlight(obj))
			continue;

		//const auto [typelabel, _atomSubTypeIdx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(obj));

		// get indices for selected atoms
		const auto [bOk, atomTypeIdx, atomSubTypeIdx, sphereIdx] = GetAtomIndexFromHandle(obj);
		if(!bOk /*|| _atomSubTypeIdx != atomSubTypeIdx*/)
		{
			QMessageBox::critical(this, PROG_NAME, "Atom handle not found, data is corrupted.");
			return std::vector<std::tuple<std::size_t, std::size_t>>{};
		}

		objs.emplace_back(std::make_tuple(atomTypeIdx, atomSubTypeIdx));
	}

	return objs;
}


/**
 * calculate the distance between selected atoms
 */
void MolDynDlg::CalculateDistanceBetweenAtoms()
{
	try
	{
		// get selected atoms
		std::vector<std::tuple<std::size_t, std::size_t>> objs = GetSelectedAtoms();

		if(objs.size() <= 1)
		{
			QMessageBox::critical(this, PROG_NAME, "At least two atoms have to be selected.");
			return;
		}


		// create file
		QString dirLast = m_sett->value("dir", "").toString();
		QString filename = QFileDialog::getSaveFileName(this, "Save File", dirLast, "Data File (*.dat)");
		if(filename == "")
			return;

		std::ofstream ofstr(filename.toStdString());
		if(!ofstr)
		{
			QMessageBox::critical(this, PROG_NAME, "Cannot open file.");
			return;
		}

		ofstr.precision(g_prec);
		m_sett->setValue("dir", QFileInfo(filename).path());


		// get coordinates of first atom
		auto [firstObjTypeIdx, firstObjSubTypeIdx] = objs[0];
		auto firstObjCoords = m_mol.GetAtomCoords(firstObjTypeIdx, firstObjSubTypeIdx);


		// output data file header infos
		ofstr << "#\n";
		ofstr << "# Column 1: Frame\n";
		ofstr << "# Columns 2, 3, ...: Distances to first atom (\xe2\x84\xab)\n";
		ofstr << "# Atoms: ";

		for(std::size_t objIdx=0; objIdx<objs.size(); ++objIdx)
		{
			auto [objTypeIdx, objSubTypeIdx] = objs[objIdx];
			ofstr << m_mol.GetAtomName(objTypeIdx) << "#" << (objSubTypeIdx+1);
			if(objIdx < objs.size()-1)
				ofstr << ", ";
		}
		ofstr << "\n#\n";


		// progress dialog
		std::shared_ptr<QProgressDialog> dlgProgress = std::make_shared<QProgressDialog>(
			"Calculating...", "Cancel", 1, objs.size(), this);
		dlgProgress->setWindowModality(Qt::WindowModal);

		// get distances to other selected atoms
		for(std::size_t objIdx=1; objIdx<objs.size(); ++objIdx)
		{
			auto [objTypeIdx, objSubTypeIdx] = objs[objIdx];
			const auto objCoords = m_mol.GetAtomCoords(objTypeIdx, objSubTypeIdx);

			for(std::size_t frameidx=0; frameidx<m_mol.GetFrameCount(); ++frameidx)
			{
				t_real dist = tl2::get_dist_uc(m_crystA, firstObjCoords[frameidx], objCoords[frameidx]);

				ofstr
					<< std::left << std::setw(g_prec*1.5) << frameidx << " "
					<< std::left << std::setw(g_prec*1.5) << dist << "\n";
			}

			dlgProgress->setValue(objIdx+1);
			if(dlgProgress->wasCanceled())
			{
				ofstr << "\n# WARNING: Calculation aborted by user.\n";
				break;
			}
		}
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, PROG_NAME, ex.what());
	}
}


/**
 * calculate positions of selected atoms
 */
void MolDynDlg::CalculatePositionsOfAtoms()
{
	try
	{
		// get selected atoms
		std::vector<std::tuple<std::size_t, std::size_t>> objs = GetSelectedAtoms();

		if(objs.size() <= 0)
		{
			QMessageBox::critical(this, PROG_NAME, "At least one atom has to be selected.");
			return;
		}


		// create file
		QString dirLast = m_sett->value("dir", "").toString();
		QString filename = QFileDialog::getSaveFileName(this, "Save File", dirLast, "Data File (*.dat)");
		if(filename == "")
			return;

		std::ofstream ofstr(filename.toStdString());
		if(!ofstr)
		{
			QMessageBox::critical(this, PROG_NAME, "Cannot open file.");
			return;
		}

		ofstr.precision(g_prec);
		m_sett->setValue("dir", QFileInfo(filename).path());


		// output data file header infos
		ofstr << "#\n";
		ofstr << "# Column 1: Frame\n";
		ofstr << "# Columns 2, 3, 4: x, y, z position of atom  (\xe2\x84\xab)\n";
		ofstr << "# Atoms: ";

		for(std::size_t objIdx=0; objIdx<objs.size(); ++objIdx)
		{
			auto [objTypeIdx, objSubTypeIdx] = objs[objIdx];
			ofstr << m_mol.GetAtomName(objTypeIdx) << "#" << (objSubTypeIdx+1);
			if(objIdx < objs.size()-1)
				ofstr << ", ";
		}
		ofstr << "\n#\n";


		// progress dialog
		std::shared_ptr<QProgressDialog> dlgProgress = std::make_shared<QProgressDialog>(
			"Calculating...", "Cancel", 0, m_mol.GetFrameCount(), this);
		dlgProgress->setWindowModality(Qt::WindowModal);

		// iterate all selected atoms
		for(std::size_t frameidx=0; frameidx<m_mol.GetFrameCount(); ++frameidx)
		{
			ofstr << std::left << std::setw(g_prec*1.5) << frameidx << " ";

			for(std::size_t objIdx=0; objIdx<objs.size(); ++objIdx)
			{
				auto [objTypeIdx, objSubTypeIdx] = objs[objIdx];
				const t_vec& coords = m_mol.GetAtomCoords(objTypeIdx, objSubTypeIdx, frameidx);

				ofstr
					<< std::left << std::setw(g_prec*1.5) << coords[0] << " "
					<< std::left << std::setw(g_prec*1.5) << coords[1] << " "
					<< std::left << std::setw(g_prec*1.5) << coords[2] << "  ";
			}

			ofstr << "\n";

			dlgProgress->setValue(frameidx+1);
			if(dlgProgress->wasCanceled())
			{
				ofstr << "\n# WARNING: Calculation aborted by user.\n";
				break;
			}
		}
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, PROG_NAME, ex.what());
	}
}


/**
 * calculate distance to initial positions of selected atoms
 */
void MolDynDlg::CalculateDeltaDistancesOfAtoms()
{
	try
	{
		// get selected atoms
		std::vector<std::tuple<std::size_t, std::size_t>> objs = GetSelectedAtoms();

		if(objs.size() <= 0)
		{
			QMessageBox::critical(this, PROG_NAME, "At least one atom has to be selected.");
			return;
		}


		// create file
		QString dirLast = m_sett->value("dir", "").toString();
		QString filename = QFileDialog::getSaveFileName(this, "Save File", dirLast, "Data File (*.dat)");
		if(filename == "")
			return;

		std::ofstream ofstr(filename.toStdString());
		if(!ofstr)
		{
			QMessageBox::critical(this, PROG_NAME, "Cannot open file.");
			return;
		}

		ofstr.precision(g_prec);
		m_sett->setValue("dir", QFileInfo(filename).path());


		// output data file header infos
		ofstr << "#\n";
		ofstr << "# Column 1: Frame\n";
		ofstr << "# Column 2, 3, ...: Distance deltas (\xe2\x84\xab)\n";
		ofstr << "# Atoms: ";

		for(std::size_t objIdx=0; objIdx<objs.size(); ++objIdx)
		{
			auto [objTypeIdx, objSubTypeIdx] = objs[objIdx];
			ofstr << m_mol.GetAtomName(objTypeIdx) << "#" << (objSubTypeIdx+1);
			if(objIdx < objs.size()-1)
				ofstr << ", ";
		}
		ofstr << "\n#\n";


		// progress dialog
		std::shared_ptr<QProgressDialog> dlgProgress = std::make_shared<QProgressDialog>(
			"Calculating...", "Cancel", 0, m_mol.GetFrameCount(), this);
		dlgProgress->setWindowModality(Qt::WindowModal);

		// iterate all selected atoms
		for(std::size_t frameidx=0; frameidx<m_mol.GetFrameCount(); ++frameidx)
		{
			ofstr << std::left << std::setw(g_prec*1.5) << frameidx << " ";

			for(std::size_t objIdx=0; objIdx<objs.size(); ++objIdx)
			{
				auto [objTypeIdx, objSubTypeIdx] = objs[objIdx];
			 	const t_vec& coords = m_mol.GetAtomCoords(objTypeIdx, objSubTypeIdx, frameidx);
				const t_vec& coordsInitial = m_mol.GetAtomCoords(objTypeIdx, objSubTypeIdx, 0);

				t_real dist = tl2::get_dist_uc(m_crystA, coords, coordsInitial);

				ofstr << std::left << std::setw(g_prec*1.5) << dist << " ";
			}

			ofstr << "\n";

			dlgProgress->setValue(frameidx+1);
			if(dlgProgress->wasCanceled())
			{
				ofstr << "\n# WARNING: Calculation aborted by user.\n";
				break;
			}
		}
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, PROG_NAME, ex.what());
	}
}


/**
 * calculate the convex hull of selected atoms
 */
void MolDynDlg::CalculateConvexHullOfAtoms()
{
	// get selected atoms
	Hull hull;
	hull.vertices = GetSelectedAtoms();

	if(hull.vertices.size() <= 3)
	{
		QMessageBox::critical(this, PROG_NAME, "At least four atoms have to be selected.");
		return;
	}

	// mark indices for hull calculation
	m_hulls.emplace_back(std::move(hull));
	CalculateConvexHulls();
}


/**
 * calculate all convex hulls for the atom indices given in m_hulls
 */
void MolDynDlg::CalculateConvexHulls()
{
#ifdef USE_QHULL
	std::size_t frameidx = m_sliderFrame->value();

	for(auto& hull : m_hulls)
	{
		// remove old plot object
		if(hull.plotObj)
		{
			m_plot->GetRenderer()->RemoveObject(*hull.plotObj);
			hull.plotObj.reset();
		}

		std::vector<t_vec> vertices;

		for(const auto [objTypeIdx, objSubTypeIdx] : hull.vertices)
		{
			const t_vec& coords = m_mol.GetAtomCoords(objTypeIdx, objSubTypeIdx, frameidx);
			vertices.push_back(coords);
		}


		auto [polys, normals, dists] = tl2_qh::get_convexhull<t_vec>(vertices);


		std::vector<t_vec3_gl> glvertices;
		std::vector<t_vec3_gl> glnormals;

		for(std::size_t polyidx=0; polyidx<polys.size(); ++polyidx)
		{
			const auto& poly = polys[polyidx];
			const auto& normal = normals[polyidx];

			for(const auto& vert : poly)
				glvertices.emplace_back(tl2::convert<t_vec3_gl, t_vec>(vert));

			glnormals.emplace_back(tl2::convert<t_vec3_gl, t_vec>(normal));
		}

		hull.plotObj = m_plot->GetRenderer()->AddTriangleObject(glvertices, glnormals, 0,0,1,0.5);

		// TODO
	}

#else
	QMessageBox::critical(this, "Error", "Calculation of convex hull is disabled.");
#endif
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
void MolDynDlg::New()
{
	m_mol.Clear();

	for(auto& hull : m_hulls)
	{
		if(hull.plotObj)
			m_plot->GetRenderer()->RemoveObject(*hull.plotObj);
	}

	m_hulls.clear();

	for(const auto& obj : m_sphereHandles)
		m_plot->GetRenderer()->RemoveObject(obj);

	m_sphereHandles.clear();
	m_sliderFrame->setValue(0);

	m_plot->update();
}


void MolDynDlg::Load()
{
	if(!m_plot) return;

	try
	{
		QString dirLast = m_sett->value("dir", "").toString();
		auto filedlg = std::make_shared<MolDynFileDlg>(this, "Load File", dirLast, "Molecular Dynamics File (*)");
		if(!filedlg->exec())
			return;
		auto files = filedlg->selectedFiles();
		if(!files.size())
			return;

		QString filename = files[0];
		if(filename == "" || !QFile::exists(filename))
			return;
		m_sett->setValue("dir", QFileInfo(filename).path());

		New();


		std::shared_ptr<QProgressDialog> dlgProgress = std::make_shared<QProgressDialog>(
			"Loading \"" + QFileInfo(filename).fileName() + "\"...", "Cancel", 0, 1000, this);
		bool bCancelled = 0;
		auto progressHandler = [dlgProgress, &bCancelled](t_real percentage) -> bool
		{
			dlgProgress->setValue(int(percentage*10));
			bCancelled = dlgProgress->wasCanceled();
			return !bCancelled;
		};
		m_mol.SubscribeToLoadProgress(progressHandler);
		dlgProgress->setWindowModality(Qt::WindowModal);


		if(!m_mol.LoadFile(filename.toStdString(), filedlg->GetFrameSkip()))
		{
			// only show error if not explicitely cancelled
			if(!bCancelled)
				QMessageBox::critical(this, PROG_NAME, "Error loading file.");
			return;
		}


		m_mol.UnsubscribeFromLoadProgress(&progressHandler);
		m_sliderFrame->setMaximum(m_mol.GetFrameCount() - 1);


		// crystal A and B matrices
		const t_vec& _a = m_mol.GetBaseA();
		const t_vec& _b = m_mol.GetBaseB();
		const t_vec& _c = m_mol.GetBaseC();

		m_crystA = tl2::create<t_mat>({
			_a[0],	_b[0],	_c[0],
			_a[1],	_b[1],	_c[1],
			_a[2], 	_b[2],	_c[2] });

		bool ok = true;
		std::tie(m_crystB, ok) = tl2::inv(m_crystA);
		if(!ok)
		{
			m_crystB = tl2::unit<t_mat>();
			QMessageBox::critical(this, PROG_NAME, "Error: Cannot invert A matrix.");
		}

		m_crystB /= t_real_gl(2)*tl2::pi<t_real_gl>;
		t_mat_gl matA{m_crystA};
		m_plot->GetRenderer()->SetBTrafo(m_crystB, &matA);

		std::cout << "A matrix: " << m_crystA << ", \n"
			<< "B matrix: " << m_crystB << "." << std::endl;


		// atom colors
		std::vector<t_vec> cols =
		{
			tl2::create<t_vec>({1, 0, 0}),
			tl2::create<t_vec>({0, 0, 1}),
			tl2::create<t_vec>({0, 0.5, 0}),
			tl2::create<t_vec>({0, 0.5, 0.5}),
			tl2::create<t_vec>({0.5, 0.5, 0}),
			tl2::create<t_vec>({0, 0, 0}),
		};

		// add atoms to 3d view
		if(m_mol.GetFrameCount())
		{
			const auto& frame = m_mol.GetFrame(0);
			m_sphereHandles.reserve(frame.GetNumAtomTypes());

			for(std::size_t atomtypeidx=0; atomtypeidx<frame.GetNumAtomTypes(); ++atomtypeidx)
			{
				const auto& coords = frame.GetCoords(atomtypeidx);

				int atomidx = 0;
				for(const t_vec& vec : coords)
				{
					t_real atomscale = m_spinScale->value();

					std::size_t handle = Add3DAtom(vec, cols[atomtypeidx % cols.size()], atomscale,
						m_mol.GetAtomName(atomtypeidx), atomidx);
					m_sphereHandles.push_back(handle);

					++atomidx;
				}
			}
		}
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, PROG_NAME, ex.what());
	}


	UpdateAtomsStatusMsg();
	m_plot->update();
}


void MolDynDlg::SaveAs()
{
	try
	{
		QString dirLast = m_sett->value("dir", "").toString();
		QString filename = QFileDialog::getSaveFileName(this, "Save File", dirLast, "Molecular Dynamics File (*)");
		if(filename == "")
			return;
		m_sett->setValue("dir", QFileInfo(filename).path());


		std::shared_ptr<QProgressDialog> dlgProgress = std::make_shared<QProgressDialog>(
			"Saving \"" + QFileInfo(filename).fileName() + "\"...", "Cancel", 0, 1000, this);
		bool bCancelled = 0;
		auto progressHandler = [dlgProgress, &bCancelled](t_real percentage) -> bool
		{
			dlgProgress->setValue(int(percentage*10));
			bCancelled = dlgProgress->wasCanceled();
			return !bCancelled;
		};
		m_mol.SubscribeToSaveProgress(progressHandler);
		dlgProgress->setWindowModality(Qt::WindowModal);


		if(!m_mol.SaveFile(filename.toStdString()))
		{
			QMessageBox::critical(this, PROG_NAME, "Error saving file.");
		}


		m_mol.UnsubscribeFromSaveProgress(&progressHandler);
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, PROG_NAME, ex.what());
	}
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
/**
 * mouse hovers over 3d object
 */
void MolDynDlg::PickerIntersection(const t_vec3_gl* pos,
	std::size_t objIdx, std::size_t /*triagIdx*/,
	const t_vec3_gl* posSphere)
{
	if(!m_plot) return;

	if(pos)
		m_curPickedObj = long(objIdx);
	else
		m_curPickedObj = -1;


	if(m_curPickedObj > 0)
	{
		const auto [typelabel, atomidx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(m_curPickedObj));

		std::ostringstream ostrLabel;
		ostrLabel << "Current Atom: " << typelabel << " #" << (atomidx+1);

		m_statusCurAtom->setText(ostrLabel.str().c_str());
		//SetStatusMsg(label);
	}
	else
	{
		m_statusCurAtom->setText("");
		//SetStatusMsg("");
	}
}



/**
 * set status label text in 3d dialog
 */
void MolDynDlg::SetStatusMsg(const std::string& msg)
{
	if(!m_status) return;
	m_status->showMessage(msg.c_str(), 2000);
}



void  MolDynDlg::UpdateAtomsStatusMsg()
{
	if(!m_statusAtoms || !m_sliderFrame) return;

	// Atoms
	std::string numAtoms = std::to_string(m_mol.GetNumAtomsTotal()) + " atoms.";

	// Selected
	if(m_plot)
	{
		std::size_t numSelected = 0;

		for(auto handle : m_sphereHandles)
        {
			if(m_plot->GetRenderer()->GetObjectHighlight(handle))
				++numSelected;
        }

		numAtoms += " " + std::to_string(numSelected) + " selected.";
	}

	// Frames
	numAtoms += " Frame " + std::to_string(m_sliderFrame->value()+1) + " of " + std::to_string(m_mol.GetFrameCount()) + ".";

	m_statusAtoms->setText(numAtoms.c_str());
}



/**
 * mouse button pressed
 */
void MolDynDlg::PlotMouseDown(bool left, bool mid, bool right)
{
	if(!m_plot) return;

	if(left && m_curPickedObj > 0)
	{
		m_plot->GetRenderer()->SetObjectHighlight(m_curPickedObj, !m_plot->GetRenderer()->GetObjectHighlight(m_curPickedObj));
		UpdateAtomsStatusMsg();
		m_plot->update();
	}
}


/**
 * mouse button released
 */
void MolDynDlg::PlotMouseUp(bool left, bool mid, bool right)
{
}


/**
 * mouse button clicked
 */
void MolDynDlg::PlotMouseClick(bool left, bool mid, bool right)
{
	// show context menu
	if(right && m_curPickedObj > 0)
	{
		const auto [typelabel, atomidx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(m_curPickedObj));

		QString atomLabel = typelabel.c_str();
		m_atomContextMenu->actions()[0]->setText("Delete This \"" + atomLabel + "\" Atom");
		m_atomContextMenu->actions()[5]->setText("Delete All \"" + atomLabel + "\" Atoms");
		m_atomContextMenu->actions()[6]->setText("Delete All But \"" + atomLabel + "\" Atoms");
		m_atomContextMenu->actions()[8]->setText("Select All \"" + atomLabel + "\" Atoms");

		auto ptGlob = QCursor::pos();
		ptGlob.setY(ptGlob.y() + 8);
		m_atomContextMenu->popup(ptGlob);
	}
}


// ----------------------------------------------------------------------------


/**
 * select all atoms
 */
void MolDynDlg::SelectAll()
{
	if(!m_plot) return;

	for(auto handle : m_sphereHandles)
		m_plot->GetRenderer()->SetObjectHighlight(handle, 1);

	UpdateAtomsStatusMsg();
	m_plot->update();
}


/**
 * unselect all atoms
 */
void MolDynDlg::SelectNone()
{
	if(!m_plot) return;

	for(auto handle : m_sphereHandles)
		m_plot->GetRenderer()->SetObjectHighlight(handle, 0);

	UpdateAtomsStatusMsg();
	m_plot->update();
}


void MolDynDlg::SliderValueChanged(int val)
{
	if(!m_plot) return;

	if(val < 0 || val >= m_mol.GetFrameCount())
		return;


	// update atom position with selected frame
	const auto& frame = m_mol.GetFrame(val);
	t_real atomscale = m_spinScale->value();

	std::size_t counter = 0;
	for(std::size_t atomtypeidx=0; atomtypeidx<frame.GetNumAtomTypes(); ++atomtypeidx)
	{
		const auto& coords = frame.GetCoords(atomtypeidx);

		int atomidx = 0;
		for(const t_vec& vec : coords)
		{
			std::size_t obj = m_sphereHandles[counter];
			Change3DAtom(obj, &vec, nullptr, &atomscale, nullptr, atomidx);

			++counter;
			++atomidx;
		}
	}

#ifdef USE_QHULL
	CalculateConvexHulls();
#endif
	UpdateAtomsStatusMsg();
	m_plot->update();
}


// ----------------------------------------------------------------------------

/**
 * extract [typelabel, atomindex] from the atom data strings
 */
std::tuple<std::string, int> MolDynDlg::SplitDataString(const std::string& data) const
{
	std::size_t idx = data.rfind(DATA_SEP);

	// no separator found, just return atom type string
	if(idx == std::string::npos)
		return std::make_tuple(data, -1);

	// split string
	std::string atomtype = data.substr(0, idx);
	std::string atomidx = data.substr(idx + std::strlen(DATA_SEP));

	return std::make_tuple(atomtype, tl2::stoval<int>(atomidx));
}


/**
 * get the index of the atom in the m_mol data structure
 * from the handle of the displayed 3d object
 */
std::tuple<bool, std::size_t, std::size_t, std::size_t>
MolDynDlg::GetAtomIndexFromHandle(std::size_t handle) const
{
	// find handle in sphere handle vector
	auto iter = std::find(m_sphereHandles.begin(), m_sphereHandles.end(), handle);
	if(iter == m_sphereHandles.end())
		return std::make_tuple(0, 0, 0, 0);

	std::size_t sphereIdx = iter - m_sphereHandles.begin();

	std::size_t atomCountsSoFar = 0;
	std::size_t atomTypeIdx = 0;
	for(atomTypeIdx=0; atomTypeIdx<m_mol.GetNumAtomTypes(); ++atomTypeIdx)
	{
		std::size_t numAtoms = m_mol.GetAtomNum(atomTypeIdx);
		if(atomCountsSoFar + numAtoms > sphereIdx)
			break;

		atomCountsSoFar += numAtoms;
	}

	std::size_t atomSubTypeIdx = sphereIdx-atomCountsSoFar;

	return std::make_tuple(1, atomTypeIdx, atomSubTypeIdx, sphereIdx);
}


/**
 * select all atoms of the same type as the one under the cursor
 */
void MolDynDlg::SelectAtomsOfSameType()
{
	// nothing under cursor
	if(m_curPickedObj <= 0)
		return;

	// atom type to be selected
	const auto [atomLabel, atomidx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(m_curPickedObj));

	for(auto handle : m_sphereHandles)
    {
		const auto [typelabel, atomidx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(handle));
		if(typelabel == atomLabel)
			m_plot->GetRenderer()->SetObjectHighlight(handle, 1);
    }

	UpdateAtomsStatusMsg();
	m_plot->update();
}


/**
 * delete the selected atoms
 */
void MolDynDlg::DeleteSelectedAtoms()
{
	std::size_t totalRemoved = 0;

	for(auto iter=m_sphereHandles.begin(); iter!=m_sphereHandles.end();)
	{
		auto handle = *iter;
		if(m_plot->GetRenderer()->GetObjectHighlight(handle))
		{
			const auto [bOk, atomTypeIdx, atomSubTypeIdx, sphereIdx] = GetAtomIndexFromHandle(handle);
			if(!bOk)
			{
				QMessageBox::critical(this, PROG_NAME, "Atom handle not found, data is corrupted.");
				return;
			}

			// remove 3d objects
			m_plot->GetRenderer()->RemoveObject(handle);
			iter = m_sphereHandles.erase(iter);

			// remove atom
			m_mol.RemoveAtom(atomTypeIdx, atomSubTypeIdx);

			++totalRemoved;
			continue;
		}

		++iter;
	}

	SetStatusMsg(std::to_string(totalRemoved) + " atoms removed.");
	UpdateAtomsStatusMsg();
	m_plot->update();
}


/**
 * delete all but the selected atoms
 */
void MolDynDlg::OnlyKeepSelectedAtoms()
{
	std::size_t totalRemoved = 0;

	for(auto iter=m_sphereHandles.begin(); iter!=m_sphereHandles.end();)
	{
		auto handle = *iter;
		if(!m_plot->GetRenderer()->GetObjectHighlight(handle))
		{
			const auto [bOk, atomTypeIdx, atomSubTypeIdx, sphereIdx] = GetAtomIndexFromHandle(handle);
			if(!bOk)
			{
				QMessageBox::critical(this, PROG_NAME, "Atom handle not found, data is corrupted.");
				return;
			}

			// remove 3d objects
			m_plot->GetRenderer()->RemoveObject(handle);
			iter = m_sphereHandles.erase(iter);

			// remove atom
			m_mol.RemoveAtom(atomTypeIdx, atomSubTypeIdx);

			++totalRemoved;
			continue;
		}

		++iter;
	}

	SetStatusMsg(std::to_string(totalRemoved) + " atoms removed.");
	UpdateAtomsStatusMsg();
	m_plot->update();
}


/**
 * delete one atom
 */
void MolDynDlg::DeleteAtomUnderCursor()
{
	// nothing under cursor
	if(m_curPickedObj <= 0)
		return;

	// atom type to be deleted
	const auto [atomLabel, atomidx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(m_curPickedObj));
	const auto [bOk, atomTypeIdx, atomSubTypeIdx, sphereIdx] = GetAtomIndexFromHandle(m_curPickedObj);
	if(!bOk)
	{
		QMessageBox::critical(this, PROG_NAME, "Atom handle not found, data is corrupted.");
		return;
	}

	if(m_mol.GetAtomName(atomTypeIdx) != atomLabel)
	{
		QMessageBox::critical(this, PROG_NAME, "Mismatch in atom type, data is corrupted.");
		return;
	}

	// remove 3d objects
	m_plot->GetRenderer()->RemoveObject(m_sphereHandles[sphereIdx]);
	m_sphereHandles.erase(m_sphereHandles.begin()+sphereIdx);

	// remove atom
	m_mol.RemoveAtom(atomTypeIdx, atomSubTypeIdx);

	SetStatusMsg("1 atom removed.");
	UpdateAtomsStatusMsg();
	m_plot->update();
}


/**
 * delete all atoms of the type under the cursor
 */
void MolDynDlg::DeleteAllAtomsOfSameType()
{
	// nothing under cursor
	if(m_curPickedObj <= 0)
		return;

	// atom type to be deleted
	const auto [atomLabel, atomidx] = SplitDataString(m_plot->GetRenderer()->GetObjectDataString(m_curPickedObj));

	std::size_t startIdx = 0;
	std::size_t totalRemoved = 0;
	for(std::size_t atomIdx=0; atomIdx<m_mol.GetNumAtomTypes();)
	{
		std::size_t numAtoms = m_mol.GetAtomNum(atomIdx);

		if(m_mol.GetAtomName(atomIdx) == atomLabel)
		{
			// remove 3d objects
			for(std::size_t sphereIdx=startIdx; sphereIdx<startIdx+numAtoms; ++sphereIdx)
				m_plot->GetRenderer()->RemoveObject(m_sphereHandles[sphereIdx]);
			m_sphereHandles.erase(m_sphereHandles.begin()+startIdx, m_sphereHandles.begin()+startIdx+numAtoms);

			// remove atoms
			m_mol.RemoveAtoms(atomIdx);

			totalRemoved += numAtoms;
		}
		else
		{
			startIdx += numAtoms;
			++atomIdx;
		}
	}

	SetStatusMsg(std::to_string(totalRemoved) + " atoms removed.");
	UpdateAtomsStatusMsg();
	m_plot->update();
}


/**
 * delete all atoms NOT of the type under the cursor
 */
void MolDynDlg::KeepAtomsOfSameType()
{
	// nothing under cursor
	if(m_curPickedObj <= 0)
		return;

	// atom type to be deleted
	const std::string& atomLabel = m_plot->GetRenderer()->GetObjectDataString(m_curPickedObj);

	std::size_t startIdx = 0;
	std::size_t totalRemoved = 0;
	for(std::size_t atomIdx=0; atomIdx<m_mol.GetNumAtomTypes();)
	{
		std::size_t numAtoms = m_mol.GetAtomNum(atomIdx);

		if(m_mol.GetAtomName(atomIdx) != atomLabel)
		{
			// remove 3d objects
			for(std::size_t sphereIdx=startIdx; sphereIdx<startIdx+numAtoms; ++sphereIdx)
				m_plot->GetRenderer()->RemoveObject(m_sphereHandles[sphereIdx]);
			m_sphereHandles.erase(m_sphereHandles.begin()+startIdx, m_sphereHandles.begin()+startIdx+numAtoms);

			// remove atoms
			m_mol.RemoveAtoms(atomIdx);

			totalRemoved += numAtoms;
		}
		else
		{
			startIdx += numAtoms;
			++atomIdx;
		}
	}

	SetStatusMsg(std::to_string(totalRemoved) + " atoms removed.");
	UpdateAtomsStatusMsg();
	m_plot->update();
}



// ----------------------------------------------------------------------------
void MolDynDlg::AfterGLInitialisation()
{
	if(!m_plot) return;

	// reference sphere for linked objects
	m_sphere = m_plot->GetRenderer()->AddSphere(0.05, 0.,0.,0., 1.,1.,1.,1.);
	m_plot->GetRenderer()->SetObjectVisible(m_sphere, false);

	// B matrix
	m_plot->GetRenderer()->SetBTrafo(m_crystB);

	// GL device info
	// TODO: outputting something here to stdout interferes with pipeproc process creation in takin main gui...
	//auto [strGlVer, strGlShaderVer, strGlVendor, strGlRenderer]
	//	= m_plot->GetRenderer()->GetGlDescr();
	//std::cout << "GL Version: " << strGlVer << ", Shader Version: " << strGlShaderVer << "." << std::endl;
	//std::cout << "GL Device: " << strGlRenderer << ", " << strGlVendor << "." << std::endl;
}


void MolDynDlg::GLInitialisationFailed()
{
	std::string err = "GL initialisation failed.";
	err += " Need GL version " + std::to_string(_GL_MAJ_VER) + "." + std::to_string(_GL_MIN_VER);
	err += " and shader version " + std::to_string(_GLSL_MAJ_VER) + std::to_string(_GLSL_MIN_VER) + "0.";
	QMessageBox::critical(this, PROG_NAME, err.c_str());
}


void MolDynDlg::closeEvent(QCloseEvent *evt)
{
	if(m_sett)
	{
		m_sett->setValue("geo", saveGeometry());
	}
}


void MolDynDlg::keyPressEvent(QKeyEvent *evt)
{
	if(evt->key()==Qt::Key_Left || evt->key()==Qt::Key_Down)
		m_sliderFrame->setValue(m_sliderFrame->value() - m_sliderFrame->singleStep());
	else if(evt->key()==Qt::Key_Right || evt->key()==Qt::Key_Up)
		m_sliderFrame->setValue(m_sliderFrame->value() + m_sliderFrame->singleStep());
	else if(evt->key()==Qt::Key_PageUp)
		m_sliderFrame->setValue(m_sliderFrame->value() + m_sliderFrame->pageStep());
	else if(evt->key()==Qt::Key_PageDown)
		m_sliderFrame->setValue(m_sliderFrame->value() - m_sliderFrame->pageStep());
	else if(evt->key()==Qt::Key_Home)
		m_sliderFrame->setValue(m_sliderFrame->minimum());
	else if(evt->key()==Qt::Key_End)
		m_sliderFrame->setValue(m_sliderFrame->maximum());

	QMainWindow::keyPressEvent(evt);
}
// ----------------------------------------------------------------------------





// ----------------------------------------------------------------------------

int main(int argc, char** argv)
{
	tl2::set_gl_format(1, _GL_MAJ_VER, _GL_MIN_VER, 0);
	tl2::set_locales();

	QApplication::addLibraryPath(QString(".") + QDir::separator() + "qtplugins");
	auto app = std::make_unique<QApplication>(argc, argv);
	QApplication::addLibraryPath(QApplication::applicationDirPath() + QDir::separator() + "qtplugins");

	auto dlg = std::make_unique<MolDynDlg>(nullptr);
	dlg->show();

	return app->exec();
}
// ----------------------------------------------------------------------------
