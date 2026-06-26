/**
 * exports
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date dec-2015
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include "taz.h"
#include "tlibs/math/geo_prim.h"
#include "tlibs/phys/atoms.h"
#include "tlibs/file/x3d.h"
#include "tlibs/log/log.h"
#include "tlibs/time/chrono.h"
#include "libs/version.h"

#include <QMessageBox>
#include <QFileDialog>
#include <QtSvg/QSvgGenerator>


using t_real = t_real_glob;

//--------------------------------------------------------------------------------
// image exports

void TazDlg::ExportReal()
{
	TasLayout *pTas = m_sceneReal.GetTasLayout();
	if(!pTas) return;

	const t_real dZoom = pTas->GetZoom();
	pTas->SetZoom(1.);
	ExportSceneSVG(m_sceneReal);
	pTas->SetZoom(dZoom);
}


void TazDlg::ExportTof()
{
	TofLayout *pTof = m_sceneTof.GetTofLayout();
	if(!pTof) return;

	const t_real dZoom = pTof->GetZoom();
	pTof->SetZoom(1.);
	ExportSceneSVG(m_sceneTof);
	pTof->SetZoom(dZoom);
}


void TazDlg::ExportRealLattice()
{
	RealLattice *pLatt = m_sceneRealLattice.GetLattice();
	if(!pLatt) return;

	const t_real dZoom = pLatt->GetZoom();
	pLatt->SetZoom(1.);
	ExportSceneSVG(m_sceneRealLattice);
	pLatt->SetZoom(dZoom);
}


void TazDlg::ExportRecip()
{
	ScatteringTriangle *pTri = m_sceneRecip.GetTriangle();
	if(!pTri) return;

	const t_real dZoom = pTri->GetZoom();
	pTri->SetZoom(1.);
	ExportSceneSVG(m_sceneRecip);
	pTri->SetZoom(dZoom);
}


void TazDlg::ExportProj()
{
	ProjLattice *pLatt = m_sceneProjRecip.GetLattice();
	if(!pLatt) return;

	const t_real dZoom = pLatt->GetZoom();
	pLatt->SetZoom(1.);
	ExportSceneSVG(m_sceneProjRecip);
	pLatt->SetZoom(dZoom);
}


void TazDlg::ExportSceneSVG(QGraphicsScene& scene)
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = m_settings.value("main/last_dir_export", "~").toString();
	QString strFile = QFileDialog::getSaveFileName(this,
		"Export SVG", strDirLast, "SVG files (*.svg *.SVG)", nullptr, fileopt);
	if(strFile == "")
		return;
	if(!strFile.endsWith(".svg", Qt::CaseInsensitive))
		strFile += ".svg";

	QRectF rect = scene.sceneRect();

	QSvgGenerator svg;
	svg.setFileName(strFile);
	svg.setSize(QSize(rect.width(), rect.height()));
	//svg.setResolution(300);
	svg.setViewBox(QRectF(0, 0, rect.width(), rect.height()));
	svg.setDescription("Created with Takin (https://dx.doi.org/10.5281/zenodo.4117437)");

	QPainter painter;
	painter.begin(&svg);
	scene.render(&painter);
	painter.end();

	std::string strDir = tl::get_dir(strFile.toStdString());
	m_settings.setValue("main/last_dir_export", QString(strDir.c_str()));
}


/**
 * export the current cut of the brillouin zone with the scattering plane
 */
void TazDlg::ExportBZCut()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = m_settings.value("main/last_dir_export", "~").toString();
	QString strFile = QFileDialog::getSaveFileName(this,
		"Export Data", strDirLast, "Data files (*.dat *.DAT)", nullptr, fileopt);
	if(strFile == "")
		return;
	if(!strFile.endsWith(".dat", Qt::CaseInsensitive))
		strFile += ".dat";


	const ScatteringTriangle *pTri = m_sceneRecip.GetTriangle();
	if(!pTri) return;

	const auto& bz = pTri->GetBZ3D();
	if(!bz.IsValid())
	{
		QMessageBox::critical(this, "Error", "3D Brillouin zone calculation is disabled or results are invalid.");
		return;
	}

	std::ofstream ofstr(strFile.toStdString());
	ofstr.precision(g_iPrec);
	ofstr << "# Brillouin zone cut vertices.\n";
	ofstr << "# Created with Takin " + std::string(TAKIN_VER) + " (https://dx.doi.org/10.5281/zenodo.4117437).\n";
	ofstr << "\n";

	const std::vector<ublas::vector<t_real>>& verts = pTri->GetBZ3DPlaneVerts();
	for(ublas::vector<t_real> vert : verts)
	{
		if(vert.size() < 2)
		{
			tl::log_err("Invalid vertex in Brillouin zone cut.");
			continue;
		}

		tl::set_eps_0(vert, g_dEps);
		for(int coordidx=0; coordidx<2; ++coordidx)
			ofstr << std::setw(g_iPrec*2.5) << std::left << vert[coordidx] << " ";
		ofstr << "\n";
	}
	ofstr.flush();


	std::string strDir = tl::get_dir(strFile.toStdString());
	m_settings.setValue("main/last_dir_export", QString(strDir.c_str()));
}


/**
 * export the 3d BZ as 3D model
 */
void TazDlg::ExportBZ3DModel()
{
	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = m_settings.value("main/last_dir_export", "~").toString();
	QString strFile = QFileDialog::getSaveFileName(this,
		"Export X3D", strDirLast, "X3D files (*.x3d *.X3D)", nullptr, fileopt);
	if(strFile == "")
		return;
	if(!strFile.endsWith(".x3d", Qt::CaseInsensitive))
		strFile += ".x3d";


	const ScatteringTriangle *pTri = m_sceneRecip.GetTriangle();
	if(!pTri) return;

	const auto& bz = pTri->GetBZ3D();
	if(!bz.IsValid())
	{
		QMessageBox::critical(this, "Error", "3D Brillouin zone calculation is disabled or results are invalid.");
		return;
	}


	using t_vec = ublas::vector<t_real>;
	tl::X3d<t_real> x3d;

	/*// vertices
	for(const t_vec& vec : bz.GetVertices())
	{
		tl::X3dTrafo<t_real> *pTrafo = new tl::X3dTrafo<t_real>();
		pTrafo->SetTrans(vec - bz.GetCentralReflex());

		tl::X3dSphere<t_real> *pSphere = new tl::X3dSphere<t_real>(0.025);
		pSphere->SetColor(tl::make_vec({1., 0., 0.}));
		pTrafo->AddChild(pSphere);

		x3d.GetScene().AddChild(pTrafo);
	}*/

	// symmetry points
	for(const t_vec& vec : pTri->GetBZ3DSymmVerts())
	{
		tl::X3dTrafo<t_real> *pTrafo = new tl::X3dTrafo<t_real>();
		pTrafo->SetTrans(vec - bz.GetCentralReflex());

		tl::X3dSphere<t_real> *pSphere = new tl::X3dSphere<t_real>(0.025);
		pSphere->SetColor(tl::make_vec({1., 0., 0.}));
		pTrafo->AddChild(pSphere);

		x3d.GetScene().AddChild(pTrafo);
	}

	// polygons
	for(const std::vector<t_vec>& vecPoly : bz.GetPolys())
	{
		tl::X3dPolygon<t_real> *pPoly = new tl::X3dPolygon<t_real>();
		pPoly->SetColor(tl::make_vec({0., 0., 1.}));

		for(const t_vec& vec : vecPoly)
			pPoly->AddVertex(vec - bz.GetCentralReflex());

		x3d.GetScene().AddChild(pPoly);
	}


	// Comment
	std::ostringstream ostrComment;
	ostrComment << "Brillouin zone contains: ";
	ostrComment << bz.GetPolys().size() << " faces, ";
	ostrComment << bz.GetVertices().size() << " vertices, ";
	ostrComment << pTri->GetBZ3DSymmVerts().size() << " symmetry points.\n";

	std::string strTakin = "\nCreated with Takin " + std::string(TAKIN_VER) + " (https://dx.doi.org/10.5281/zenodo.4117437).\n";
	strTakin += "Timestamp: " + tl::epoch_to_str<t_real>(tl::epoch<t_real>()) + "\n\n";
	x3d.SetComment(strTakin + ostrComment.str());

	bool bOk = x3d.Save(strFile.toStdString().c_str());


	if(!bOk)
		QMessageBox::critical(this, "Error", "Could not export 3D model.");

	if(bOk)
	{
		std::string strDir = tl::get_dir(strFile.toStdString());
		m_settings.setValue("main/last_dir_export", QString(strDir.c_str()));
	}
}



//--------------------------------------------------------------------------------
// 3d model exports

void TazDlg::ExportUCModel()
{
	using t_mat = ublas::matrix<t_real>;
	using t_vec = ublas::vector<t_real>;


	if(m_vecAtoms.size() == 0)
	{
		QMessageBox::critical(this, "Error", "No atom positions defined for unit cell.");
		return;
	}

	RealLattice *pReal = m_sceneRealLattice.GetLattice();
	if(!pReal) return;

	const auto& ws = pReal->GetWS3D();
	if(!ws.IsValid())
	{
		QMessageBox::critical(this, "Error", "3D unit cell calculation is disabled or results are invalid.");
		return;
	}

	xtl::SpaceGroup<t_real>* pSpaceGroup = nullptr;
	int iSpaceGroupIdx = comboSpaceGroups->currentIndex();
	if(iSpaceGroupIdx != 0)
		pSpaceGroup = (xtl::SpaceGroup<t_real>*)comboSpaceGroups->itemData(iSpaceGroupIdx).value<void*>();
	if(!pSpaceGroup)
	{
		QMessageBox::critical(this, "Error", "Invalid space group.");
		return;
	}

	std::vector<t_vec> vecAtoms;
	std::vector<std::string> vecAtomNames;
	for(const xtl::AtomPos<t_real>& atom : m_vecAtoms)
	{
		t_vec vecAtomPos = atom.vecPos;
		vecAtomPos.resize(4, 1);
		vecAtomPos[3] = 1.;
		vecAtoms.push_back(std::move(vecAtomPos));
		vecAtomNames.push_back(atom.strAtomName);
	}


	const std::vector<t_vec> vecColors =
	{
		tl::make_vec({0., 0., 1.}), tl::make_vec({1., 0., 0.}), tl::make_vec({0., 1., 0.}),
		tl::make_vec({1., 1., 0.}), tl::make_vec({0., 1., 1.}), tl::make_vec({1., 0., 1.}),
		tl::make_vec({0.5, 0., 0.}), tl::make_vec({0., 0.5, 0.}), tl::make_vec({0., 0., 0.5}),
		tl::make_vec({0.5, 0.5, 0.}), tl::make_vec({0., 0.5, 0.5}), tl::make_vec({0.5, 0., 0.5}),
		tl::make_vec({1., 1., 1.}), tl::make_vec({0., 0., 0.}), tl::make_vec({0.5, 0.5, 0.5}),
	};

	const std::vector<t_real> vecRadii =
	{
		0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2
	};
	const t_real dRadScale = 1.;

	// to transform into program-specific coordinate systems
	const t_mat matGlobal = tl::make_mat(
	{	{-1., 0., 0., 0.},
		{ 0., 0., 1., 0.},
		{ 0., 1., 0., 0.},
		{ 0., 0., 0., 1.}	});


	QFileDialog::Option fileopt = QFileDialog::Option(0);
	if(!m_settings.value("main/native_dialogs", 1).toBool())
		fileopt = QFileDialog::DontUseNativeDialog;

	QString strDirLast = m_settings.value("main/last_dir_export", "~").toString();
	QString strFile = QFileDialog::getSaveFileName(this,
		"Export X3D", strDirLast, "X3D files (*.x3d *.X3D)", nullptr, fileopt);
	if(strFile == "")
		return;
	if(!strFile.endsWith(".x3d", Qt::CaseInsensitive))
		strFile += ".x3d";

	tl::X3d<t_real> x3d;

	for(std::size_t iAtom=0; iAtom<m_latticecommon.vecAllAtoms.size(); ++iAtom)
	//for(std::size_t iAtom=0; iAtom<m_latticecommon.vecIdxSC.size(); ++iAtom)
	{
		std::size_t iAtomType = m_latticecommon.vecAllAtomTypes[iAtom];
		t_vec vecCoord = m_latticecommon.vecAllAtoms[iAtom];

		vecCoord.resize(4,1); vecCoord[3] = 1.;

		tl::X3dTrafo<t_real> *pTrafo = new tl::X3dTrafo<t_real>();
		pTrafo->SetTrans(vecCoord);

		t_real dRadius = vecRadii[iAtomType % vecRadii.size()] * dRadScale;
		tl::X3dSphere<t_real> *pSphere = new tl::X3dSphere<t_real>(dRadius);
		pSphere->SetColor(vecColors[iAtomType % vecColors.size()]);
		pTrafo->AddChild(pSphere);


		// coordination polyhedra
		if(m_latticecommon.vecAllAtomPosAux.size() == m_latticecommon.vecAllAtoms.size())
		{
			for(const auto& vecPoly : m_latticecommon.vecAllAtomPosAux[iAtom].vecPolys)
			{
				tl::X3dPolygon<t_real> *pPoly = new tl::X3dPolygon<t_real>();
				pPoly->SetColor(vecColors[iAtomType % vecColors.size()]);

				for(const t_vec& vecVertex : vecPoly)
					pPoly->AddVertex(vecVertex);

				pTrafo->AddChild(pPoly);
			}
		}

		x3d.GetScene().AddChild(pTrafo);
	}


	// unit cell
	for(const std::vector<t_vec>& vecPoly : ws.GetPolys())
	{
		tl::X3dLines<t_real> *pPoly = new tl::X3dLines<t_real>();
		pPoly->SetColor(tl::make_vec({0., 0., 0., 1.}));

		for(const t_vec& vec : vecPoly)
			pPoly->AddVertex(vec - ws.GetCentralReflex());

		x3d.GetScene().AddChild(pPoly);
	}


	// comment
	std::ostringstream ostrComment;
	ostrComment << "Unit cell.\n";

	for(std::size_t iAtom=0; iAtom<vecAtoms.size(); ++iAtom)
	{
		ostrComment << "Atom " << (iAtom+1) << ": " << vecAtomNames[iAtom]
			<< ", colour: " << vecColors[iAtom % vecColors.size()] << "\n";
	}
	ostrComment << "Unit cell contains " << m_latticecommon.vecAllAtoms.size() << " atoms.\n";
	ostrComment << "Super cell contains " << m_latticecommon.vecIdxSC.size() << " atoms.\n";

	tl::log_info(ostrComment.str());
	std::string strTakin = "\nCreated with Takin " + std::string(TAKIN_VER) + " (https://dx.doi.org/10.5281/zenodo.4117437).\n";
	strTakin += "Timestamp: " + tl::epoch_to_str<t_real>(tl::epoch<t_real>()) + "\n\n";
	x3d.SetComment(strTakin + ostrComment.str());


	bool bOk = x3d.Save(strFile.toStdString().c_str());

	if(bOk)
	{
		std::string strDir = tl::get_dir(strFile.toStdString());
		m_settings.setValue("main/last_dir_export", QString(strDir.c_str()));
	}
	else
	{
		QMessageBox::critical(this, "Error", "Error exporting x3d file.");
	}
}
