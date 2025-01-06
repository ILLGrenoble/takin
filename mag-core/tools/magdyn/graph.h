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

#ifndef __MAG_DYN_GRAPH_H__
#define __MAG_DYN_GRAPH_H__

#include <qcustomplot.h>
#include <QtCore/QVector>



/**
 * a graph with weight factors per data point
 */
class GraphWithWeights : public QCPGraph
{
public:
	GraphWithWeights(QCPAxis *x, QCPAxis *y);
	virtual ~GraphWithWeights();

	/**
	 * sets the symbol sizes
	 * setData() needs to be called with already_sorted=true,
	 * otherwise points and weights don't match
	 */
	void SetWeights(const QVector<qreal>& weights);
	void SetWeightScale(qreal sc, qreal min, qreal max);
	void SetWeightAsPointSize(bool b);
	void SetWeightAsAlpha(bool b);

	void AddColour(const QColor& col);
	void SetColourIndices(const QVector<int>& cols);


	/**
	 * scatter plot with variable symbol sizes
	 */
	virtual void drawScatterPlot(
		QCPPainter* paint,
		const QVector<QPointF>& points,
		const QCPScatterStyle& _style) const override;


private:
	// symbol sizes
	QVector<qreal> m_weights{};

	// colour palette
	QVector<QColor> m_colours{};
	QVector<int> m_colour_indices{};

	qreal m_weight_scale = 1;
	qreal m_weight_min = -1;
	qreal m_weight_max = -1;

	bool m_weight_as_point_size = true;
	bool m_weight_as_alpha = false;
};


#endif
