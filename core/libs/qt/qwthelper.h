/**
 * qwt helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date feb-2016
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

#ifndef __QWT_HELPER_H__
#define __QWT_HELPER_H__

#include <vector>
#include <tuple>
#include <string>
#include <type_traits>
#include <mutex>

#include <qwt_plot.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_raster_data.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_picker.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#include <qwt_legend.h>
#include <qwt_color_map.h>

#include "tlibs/string/string.h"


using t_real_qwt = double;   // qwt's intrinsic value type


class MyQwtRasterData : public QwtRasterData
{
protected:
	std::shared_ptr<t_real_qwt> m_pData;
	std::size_t m_iW = 0, m_iH = 0;
	t_real_qwt m_dXRange[2], m_dYRange[2], m_dZRange[2];

public:
	void Init(std::size_t iW, std::size_t iH)
	{
		if(m_iW == iW && m_iH == iH) return;

		m_iW = iW;
		m_iH = iH;

		m_pData.reset(new t_real_qwt[iW*iH], [](t_real_qwt *pArr) { delete[] pArr; } );
		std::fill(m_pData.get(), m_pData.get()+m_iW*m_iH, t_real_qwt(0));
	}

	MyQwtRasterData(std::size_t iW=0, std::size_t iH=0)
	{
		Init(iW, iH);
	}

	MyQwtRasterData(const MyQwtRasterData& dat);
	const MyQwtRasterData& operator=(const MyQwtRasterData& dat);

	virtual ~MyQwtRasterData();

	void SetXRange(t_real_qwt dMin, t_real_qwt dMax);
	void SetYRange(t_real_qwt dMin, t_real_qwt dMax);
	void SetZRange(t_real_qwt dMin, t_real_qwt dMax);
	void SetZRange();	// automatically determined range

	t_real_qwt GetXMin() const { return m_dXRange[0]; }
	t_real_qwt GetXMax() const { return m_dXRange[1]; }
	t_real_qwt GetYMin() const { return m_dYRange[0]; }
	t_real_qwt GetYMax() const { return m_dYRange[1]; }
	t_real_qwt GetZMin() const { return m_dZRange[0]; }
	t_real_qwt GetZMax() const { return m_dZRange[1]; }

	std::size_t GetWidth() const { return m_iW; }
	std::size_t GetHeight() const { return m_iH; }

	void SetPixel(std::size_t iX, std::size_t iY, t_real_qwt dVal)
	{
		if(iX<m_iW && iY<m_iH)
			m_pData.get()[iY*m_iW + iX] = dVal;
	}

	t_real_qwt GetPixel(std::size_t iX, std::size_t iY) const
	{
		if(!m_pData) return t_real_qwt(0);
		return m_pData.get()[iY*m_iW + iX];
	}

	virtual t_real_qwt value(t_real_qwt dx, t_real_qwt dy) const override;

	virtual QwtRasterData* clone() const;
	virtual QwtRasterData* copy() const /*override only for qwt5*/;
	virtual QwtInterval range() const /*override only for qwt5*/;

#if QWT_VERSION >= 0x060200
	virtual QwtInterval interval(Qt::Axis ax) const override;
#endif
};


// ----------------------------------------------------------------------------



class MyQwtCurve : public QwtPlotCurve
{
protected:
	bool m_bShowErrors = false;
	const std::vector<t_real_qwt>* m_pvecYErr = nullptr;

protected:
	virtual void drawDots(QPainter* pPainter,
		const QwtScaleMap& scX, const QwtScaleMap& scY,
		const QRectF& rect, int iPtFirst, int iPtLast) const /*override*/;

public:
	MyQwtCurve();
	virtual ~MyQwtCurve();

	void SetShowErrors(bool bErr) { m_bShowErrors = bErr; }
	bool GetShowErrors() const { return m_bShowErrors; }

	void SetYErr(const std::vector<t_real_qwt>* pvecYErr) { m_pvecYErr = pvecYErr; }
};


// ----------------------------------------------------------------------------


class QwtPlotWrapper : public QObject
{ Q_OBJECT
protected:
	mutable std::mutex m_mutex;

	QwtPlot *m_pPlot = nullptr;
	QwtPlotGrid *m_pGrid = nullptr;
	QwtPlotPicker *m_pPicker = nullptr;
	QwtPlotZoomer *m_pZoomer = nullptr;
	QwtPlotPanner *m_pPanner = nullptr;
	QwtLegend *m_pLegend = nullptr;

	int m_curYScaler = 0;

	std::vector<MyQwtCurve*> m_vecCurves;	// 1d plots
	QwtPlotSpectrogram *m_pSpec = nullptr;	// 2d plot

	bool m_bHasDataPtrs = true;
	std::vector<std::tuple<
		const std::vector<t_real_qwt>*,		// x data
		const std::vector<t_real_qwt>*,		// y data
		const std::vector<t_real_qwt>*		// y error
			>> m_vecDataPtrs;
	MyQwtRasterData *m_pRaster = nullptr;

public:
	QwtPlotWrapper(QwtPlot *pPlot, unsigned int iNumCurves = 1,
		bool bNoTrackerSignal = false, bool bUseSpline = false,
		bool bUseSpectrogram = false, bool bShowLegend = false);
	virtual ~QwtPlotWrapper();

	QwtPlot* GetPlot() { return m_pPlot; }
	MyQwtCurve* GetCurve(unsigned int iCurve = 0) { return m_vecCurves[iCurve]; }
	unsigned int GetNumCurves() const { return m_vecCurves.size(); }
	QwtPlotZoomer* GetZoomer() { return m_pZoomer; }
	QwtPlotPicker* GetPicker() { return m_pPicker; }
	MyQwtRasterData* GetRaster() { return m_pRaster; }
	bool HasTrackerSignal() const;

	void SetData(const std::vector<t_real_qwt>& vecX, const std::vector<t_real_qwt>& vecY,
		unsigned int iCurve = 0, bool bReplot = true, bool bCopy = false,
		const std::vector<t_real_qwt>* pvecYErr = nullptr);

	std::mutex& GetMutex() { return m_mutex; }


public slots:
	void setAxisTitle(int iAxis, const QString& str);
	void scaleColorBar();
	void setZoomBase(const QRectF&);
	void ToggleLogY();
	void doUpdate();

	void SavePlotGraphics() const;	// save plot graphics
	void SavePlot() const;		// save plot data
	void ExportGpl() const;
};


// ----------------------------------------------------------------------------


template<typename t_real, bool bSetDirectly = std::is_same<t_real, t_real_qwt>::value>
struct set_qwt_data
{
	void operator()(QwtPlotWrapper& plot,
		const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
		unsigned int iCurve = 0, bool bReplot = true,
		const std::vector<t_real>* pvecYErr = nullptr) {}
};

// same types -> set data directly
template<typename t_real>
struct set_qwt_data<t_real, 1>
{
	void operator()(QwtPlotWrapper& plot,
		const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
		unsigned int iCurve = 0, bool bReplot = true,
		const std::vector<t_real>* pvecYErr = nullptr)
	{
		// copy pointers
		plot.SetData(vecX, vecY, iCurve, bReplot, false, pvecYErr);
	}
};

// different types -> copy & convert data first
// TODO: errorbars
template<typename t_real>
struct set_qwt_data<t_real, 0>
{
	void operator()(QwtPlotWrapper& plot,
		const std::vector<t_real>& vecX, const std::vector<t_real>& vecY,
		unsigned int iCurve = 0, bool bReplot = true,
		const std::vector<t_real>* pvecYErr = nullptr)
	{
		std::vector<t_real_qwt> vecNewX, vecNewY;
		vecNewX.reserve(vecX.size());
		vecNewY.reserve(vecY.size());

		for(t_real d : vecX) vecNewX.push_back(d);
		for(t_real d : vecY) vecNewY.push_back(d);

		// copy data
		plot.SetData(vecNewX, vecNewY, iCurve, bReplot, true);
	}
};


// ----------------------------------------------------------------------------


extern void set_zoomer_base(QwtPlotZoomer *pZoomer,
	t_real_qwt dL, t_real_qwt dR, t_real_qwt dT, t_real_qwt dB,
	bool bMetaCall = false, QwtPlotWrapper* pPlotWrap = nullptr);

extern void set_zoomer_base(QwtPlotZoomer *pZoomer,
	const std::vector<t_real_qwt>& vecX, const std::vector<t_real_qwt>& vecY,
	bool bMetaCall = false, QwtPlotWrapper* pPlotWrap = nullptr,
	bool bUseYErrs = false, const std::vector<t_real_qwt> *vecYErr = nullptr);

extern void set_zoomer_base(QwtPlotZoomer *pZoomer,
	const std::vector<std::vector<t_real_qwt>>& vecvecX,
	const std::vector<std::vector<t_real_qwt>>& vecvecY,
	bool bMetaCall = false, QwtPlotWrapper* pPlotWrap = nullptr,
	bool bUseYErrs = false, const std::vector<std::vector<t_real_qwt>> *vecvecYErr = nullptr);

extern void set_zoomer_base(QwtPlotZoomer *pZoomer,
	const std::vector<t_real_qwt>& vecX,
	const std::vector<std::vector<t_real_qwt>>& vecvecY,
	bool bMetaCall = false, QwtPlotWrapper* pPlotWrap = nullptr,
	bool bUseYErrs = false);


#endif
