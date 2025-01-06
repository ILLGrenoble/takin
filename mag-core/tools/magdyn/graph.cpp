/**
 * graph with weight factors
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "graph.h"

#include <algorithm>

#include "tlibs2/libs/maths.h"



/**
 * a graph with weight factors per data point
 */
GraphWithWeights::GraphWithWeights(QCPAxis *x, QCPAxis *y)
	: QCPGraph(x, y)
{
	// can't use adaptive sampling because we have to match weights to data points
	setAdaptiveSampling(false);
}



GraphWithWeights::~GraphWithWeights()
{}



/**
 * sets the symbol sizes
 * setData() needs to be called with already_sorted=true,
 * otherwise points and weights don't match
 */
void GraphWithWeights::SetWeights(const QVector<qreal>& weights)
{
	m_weights = weights;
}



void GraphWithWeights::SetWeightScale(qreal sc, qreal min, qreal max)
{
	m_weight_scale = sc;
	m_weight_min = min;
	m_weight_max = max;
}



void GraphWithWeights::SetWeightAsPointSize(bool b)
{
	m_weight_as_point_size = b;
}



void GraphWithWeights::SetWeightAsAlpha(bool b)
{
	m_weight_as_alpha = b;
}



void GraphWithWeights::AddColour(const QColor& col)
{
	m_colours.push_back(col);
}



void GraphWithWeights::SetColourIndices(const QVector<int>& cols)
{
	m_colour_indices = cols;
}



/**
 * scatter plot with variable symbol sizes
 */
void GraphWithWeights::drawScatterPlot(
	QCPPainter* paint,
	const QVector<QPointF>& points,
	const QCPScatterStyle& _style) const
{
	const int num_points = points.size();
	if(!num_points)
		return;

	const qreal eps = 1e-3;

	// find the absolute data start index for indexing the weight data
	// (more elegant would be if we could get the indices for the data to be drawn directly from qcp)
	int data_start_idx = -1;
	for(auto iter = mDataContainer->constBegin(); iter != mDataContainer->constEnd(); ++iter)
	{
		qreal pt_x = keyAxis()->coordToPixel(iter->mainKey());
		qreal pt_y = valueAxis()->coordToPixel(iter->mainValue());

		if(tl2::equals(pt_x, points[0].x(), eps) &&
			tl2::equals(pt_y, points[0].y(), eps))
		{
			data_start_idx = iter - mDataContainer->constBegin();
			break;
		}
	}

	// need to overwrite point size
	QCPScatterStyle& style = const_cast<QCPScatterStyle&>(_style);

	// see: QCPGraph::drawScatterPlot
	QCPGraph::applyScattersAntialiasingHint(paint);

	QPen pen = this->pen();
	style.applyTo(paint, pen);

	const bool has_weights = (data_start_idx >= 0 && m_weights.size()-data_start_idx >= num_points);
	const qreal size_saved = style.size();

	// iterate all data points
	for(int idx = 0; idx < num_points; ++idx)
	{
		// set colour from index if given
		if(idx < m_colour_indices.size())
		{
			int col_idx = m_colour_indices[idx];
			if(col_idx >= m_colours.size())
				col_idx = m_colours.size() - 1;

			if(col_idx >= 0)
				pen.setColor(m_colours[col_idx]);

			style.setPen(pen);
			style.applyTo(paint, pen);
		}

		// data point
		const QPointF& pt = points[idx];
		bool weight_range_valid = m_weight_max >= 0. && m_weight_min >= 0. && m_weight_min <= m_weight_max;

		if(m_weight_as_point_size)
		{
			// size-based weight factor
			qreal weight = has_weights ? m_weights[idx + data_start_idx] : size_saved;
			weight *= m_weight_scale;
			if(weight_range_valid)
				weight = std::clamp(weight, m_weight_min, m_weight_max);

			// set symbol sizes per point
			style.setSize(weight);
		}
		else
		{
			style.setSize(2.);
		}

		if(m_weight_as_alpha)
		{
			// colour-based weight factor
			qreal weight_c = has_weights
				? m_weights[idx + data_start_idx] * m_weight_scale
				: 1.;
			if(has_weights && weight_range_valid)
			{
				weight_c = std::clamp(weight_c, m_weight_min, m_weight_max);
				weight_c /= m_weight_max;

				//weight_c = (weight_c - m_weight_min) / (m_weight_max - m_weight_min);
				//weight_c = std::clamp(weight_c, 0., 1.);
			}

			QColor col = pen.color();
			col.setAlphaF(weight_c);
			pen.setColor(col);
			style.setPen(pen);
			style.applyTo(paint, pen);
		}

		// draw the symbol with the modified size
		style.drawShape(paint, pt);
	}

	// restore original symbol size
	style.setSize(size_saved);
}
