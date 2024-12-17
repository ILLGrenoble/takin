/**
 * brillouin zone tool -- 2d plot
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

#include "plot_cut.h"
#include <QtWidgets/QApplication>


// --------------------------------------------------------------------------------
BZCutScene::BZCutScene(QWidget *parent) : QGraphicsScene(parent)
{
}


BZCutScene::~BZCutScene()
{
	ClearAll();
}


/**
 * adds the line segments of a brillouin zone cut
 * @arg [start, end, Q]
 */
void BZCutScene::AddCut(
	const std::vector<std::tuple<t_vec, t_vec, std::array<t_real, 3>>>& lines)
{
	// (000) brillouin zone
	std::vector<const std::tuple<t_vec, t_vec, std::array<t_real, 3>>*> lines000;

	QPen pen;
	pen.setCosmetic(true);
	pen.setColor(qApp->palette().color(QPalette::WindowText));
	pen.setWidthF(1.);

	m_bzcut.reserve(lines.size() * 2);

	// draw brillouin zones
	for(const auto& line : lines)
	{
		const auto& Q = std::get<2>(line);

		// (000) BZ?
		if(tl2::equals_0(Q[0], g_eps) &&
			tl2::equals_0(Q[1], g_eps) &&
			tl2::equals_0(Q[2], g_eps))
		{
			lines000.push_back(&line);
			continue;
		}

		QGraphicsLineItem *plot_line = addLine(QLineF(
			std::get<0>(line)[0]*m_scale, std::get<0>(line)[1]*m_scale,
			std::get<1>(line)[0]*m_scale, std::get<1>(line)[1]*m_scale),
			pen);

		m_bzcut.push_back(plot_line);
	}


	// draw (000) brillouin zone
	pen.setColor(QColor(0xff, 0x00, 0x00));
	pen.setWidthF(2.);

	for(const auto* line : lines000)
	{
		QGraphicsLineItem *plot_line = addLine(QLineF(
			std::get<0>(*line)[0]*m_scale, std::get<0>(*line)[1]*m_scale,
			std::get<1>(*line)[0]*m_scale, std::get<1>(*line)[1]*m_scale),
			pen);

		m_bzcut.push_back(plot_line);
	}
}


/**
 * adds a plot curve from a set of points
 */
void BZCutScene::AddCurve(const std::vector<t_vec>& points)
{
	if(points.size() < 2)
		return;

	QPen pen;
	pen.setCosmetic(true);
	pen.setColor(QColor(0x00, 0x00, 0xff));
	pen.setWidthF(2.);

	for(std::size_t i = 0; i < points.size() - 1; ++i)
	{
		std::size_t j = i + 1;

		QGraphicsLineItem *plot_line = addLine(QLineF(
			points[i][0]*m_scale, points[i][1]*m_scale,
			points[j][0]*m_scale, points[j][1]*m_scale),
			pen);

		m_curves.push_back(plot_line);
	}
}


void BZCutScene::ClearAll()
{
	ClearCut();
	ClearCurves();
	clear();
}


void BZCutScene::ClearCut()
{
	for(QGraphicsItem *item : m_bzcut)
		delete item;
	m_bzcut.clear();
}


void BZCutScene::ClearCurves()
{
	for(QGraphicsItem *item : m_curves)
		delete item;
	m_curves.clear();
}


/**
 * get the centre of the brillouin zone
 */
QPointF BZCutScene::GetCentre() const
{
	QPointF centre{0, 0};

	for(const QGraphicsItem* item : m_bzcut)
		centre += item->pos();
	centre /= qreal(m_bzcut.size());

	return centre;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
BZCutView::BZCutView(BZCutScene* scene)
	: QGraphicsView(scene, static_cast<QWidget*>(scene->parent())),
	m_scene(scene)
{
	setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	setDragMode(QGraphicsView::ScrollHandDrag);
	setInteractive(true);
	setMouseTracking(true);
	scale(1., -1.);
}


BZCutView::~BZCutView()
{
}


/**
 * centre the view
 */
void BZCutView::Centre()
{
	QPointF centre{0, 0};
	if(m_scene)
		centre = m_scene->GetCentre();

	centerOn(centre);
}


void BZCutView::mouseMoveEvent(QMouseEvent *evt)
{
	QPointF pos = mapToScene(evt->pos());
	t_real scale = m_scene->GetScale();
	emit SignalMouseCoordinates(pos.x()/scale, pos.y()/scale);

	QGraphicsView::mouseMoveEvent(evt);
}


void BZCutView::wheelEvent(QWheelEvent *evt)
{
	t_real sc = std::pow(2., evt->angleDelta().y()/8.*0.01);
	QGraphicsView::scale(sc, sc);
}
// --------------------------------------------------------------------------------
