import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {
  statisticsProperties,
  FitConfidenceIntervals,
  IFitChartData,
  CONFIDENCE_INTERVAL_FILL_COLOR,
  CONFIDENCE_INTERVAL_STROKE_COLOR,
  CURVE_CONFIDENCE_INTERVAL_BOUNDS,
  FIT_CELL_TYPE,
  IFitSeries,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {BoxPlotStatistics, calculateBoxPlotStatistics} from '@datagrok-libraries/statistics/src/box-plot-statistics';
import {Viewport} from '@datagrok-libraries/utils/src/transform';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';

import {
  fitSeries,
  getChartData,
  getChartBounds,
  getSeriesFitFunction,
  getSeriesConfidenceInterval,
  getSeriesStatistics,
  getCurve,
} from '@datagrok-libraries/statistics/src/fit/fit-data';

import {convertXMLToIFitChartData} from './fit-parser';
import {MultiCurveViewer} from './multi-curve-viewer';


export const TAG_FIT_CHART_FORMAT = '.fitChartFormat';
export const TAG_FIT_CHART_FORMAT_3DX = '3dx';
export const CANDLESTICK_BORDER_PX_SIZE = 4;
export const OUTLIER_HITBOX_RADIUS = 2;


/** Performs a chart layout, returning [viewport, xAxis, yAxis] */
export function layoutChart(rect: DG.Rect): [DG.Rect, DG.Rect?, DG.Rect?] {
  if (rect.width < 100 || rect.height < 100)
    return [rect, undefined, undefined];
  return [
    rect.cutLeft(30).cutBottom(30),
    rect.getBottom(30).cutLeft(30),
    rect.getLeft(30).cutBottom(30)
  ];
}

/** Performs candlestick border drawing */
function drawCandlestickBorder(g: CanvasRenderingContext2D, x: number, adjacentValue: number, transform: Viewport): void {
  g.moveTo(transform.xToScreen(x) - (CANDLESTICK_BORDER_PX_SIZE / 2), transform.yToScreen(adjacentValue));
  g.lineTo(transform.xToScreen(x) + (CANDLESTICK_BORDER_PX_SIZE / 2), transform.yToScreen(adjacentValue));
}

/** Performs candlestick drawing */
function drawCandlestick(g: CanvasRenderingContext2D, x: number, boxPlotStats: BoxPlotStatistics,
  transform: Viewport, ratio: number, markerColor: number): void {
  drawCandlestickBorder(g, x, boxPlotStats.lowerAdjacentValue, transform);
  g.moveTo(transform.xToScreen(x), transform.yToScreen(boxPlotStats.lowerAdjacentValue));
  g.lineTo(transform.xToScreen(x), transform.yToScreen(boxPlotStats.upperAdjacentValue));
  drawCandlestickBorder(g, x, boxPlotStats.upperAdjacentValue, transform);
  DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, transform.xToScreen(x), transform.yToScreen(boxPlotStats.q2),
    markerColor, 3.5 * ratio);
}

/** Performs points drawing */
function drawPoints(g: CanvasRenderingContext2D, series: IFitSeries,
  transform: Viewport, ratio: number): void {
  for (let i = 0; i < series.points.length!; i++) {
    const p = series.points[i];
    DG.Paint.marker(g,
      p.outlier ? DG.MARKER_TYPE.OUTLIER : DG.MARKER_TYPE.CIRCLE,
      transform.xToScreen(p.x), transform.yToScreen(p.y),
      series.pointColor ? DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker,
      (p.outlier ? 6 : 4) * ratio);
  }
}

/** Performs candles drawing */
function drawCandles(g: CanvasRenderingContext2D, series: IFitSeries,
  transform: Viewport, ratio: number) : void {
  for (let i = 0, candleStart = null; i < series.points.length!; i++) {
    const p = series.points[i];
    if (p.outlier)
      continue;
    const nextSame = i + 1 < series.points.length && series.points[i + 1].x === p.x;
    if (!candleStart && nextSame)
      candleStart = i;
    else if (candleStart !== null && !nextSame) {
      const values: number[] = [];
      for (let j = candleStart, ind = 0; j <= i; j++, ind++) {
        values[ind] = series.points[j].y;
      }
      const boxPlotStats = calculateBoxPlotStatistics(values);

      g.beginPath();
      drawCandlestick(g, p.x, boxPlotStats, transform, ratio, series.pointColor ?
        DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker);
      g.stroke();

      if (series.showPoints === 'both') {
        for (let ind = 0; ind < values.length; ind++) {
          if (values[ind] < boxPlotStats.lowerAdjacentValue || values[ind] > boxPlotStats.upperAdjacentValue) {
            DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER,
              transform.xToScreen(p.x), transform.yToScreen(values[ind]),
              series.pointColor ? DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker,
              6 * ratio);
          }
        }
      }

      candleStart = null;
    }
  }
}

/** Performs a curve confidence interval drawing */
function drawConfidenceInterval(g: CanvasRenderingContext2D, confIntervals: FitConfidenceIntervals,
  screenBounds: DG.Rect, transform: Viewport, confidenceType: string): void {
  g.beginPath();
  for (let i = 0; i <= screenBounds.width; i++) {
    const x = screenBounds.x + i;
    const y = confidenceType === CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
      transform.yToScreen(confIntervals.confidenceTop(transform.xToWorld(x))) :
      transform.yToScreen(confIntervals.confidenceBottom(transform.xToWorld(x)));
    if (i === 0)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }
  g.stroke();
}

// TODO: Denis will look at axes
// TODO: automatic output of parameters
// TODO: add candlesticks to Curves documentation
// TODO: add how to control options - series opitons, chart options, etc.

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, confIntervals: FitConfidenceIntervals,
  screenBounds: DG.Rect, transform: Viewport): void {
  g.beginPath();
  for (let i = 0; i <= screenBounds.width; i++) {
    const x = screenBounds.x + i;
    const y = transform.yToScreen(confIntervals.confidenceTop(transform.xToWorld(x)));
    if (i === 0)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }

  // reverse traverse to make a shape of confidence interval to fill it
  for (let i = screenBounds.width; i >= 0; i--) {
    const x = screenBounds.x + i;
    const y = transform.yToScreen(confIntervals.confidenceBottom(transform.xToWorld(x)));
    g.lineTo(x, y);
  }
  g.closePath();
  g.fill();
}

export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return FIT_CELL_TYPE; }

  get cellType() { return FIT_CELL_TYPE; }

  getDefaultSize(gridColumn: DG.GridColumn): {width?: number | null, height?: number | null} {
    return {width: 160, height: 100};
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    grok.shell.o = gridCell;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    const screenBounds = gridCell.bounds.inflate(-6, -6);
    const dataBox = layoutChart(screenBounds)[0];
    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);

    for (let i = 0; i < data.series?.length!; i++) {
      if (!data.series![i].clickToToggle || data.series![i].showPoints === 'candlesticks')
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        const screenX = viewport.xToScreen(p.x);
        const screenY = viewport.yToScreen(p.y);
        const pxPerMarkerType = ((p.outlier ? 6 : 4) / 2) + OUTLIER_HITBOX_RADIUS;
        if (e.offsetX >= screenX - pxPerMarkerType && e.offsetX <= screenX + pxPerMarkerType &&
          e.offsetY >= screenY - pxPerMarkerType && e.offsetY <= screenY + pxPerMarkerType) {
          p.outlier = !p.outlier;
          // temporarily works only for JSON structure
          if (gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) !== TAG_FIT_CHART_FORMAT_3DX) {
            const gridCellValue = JSON.parse(gridCell.cell.value) as IFitChartData;
            gridCellValue.series![i].points[j].outlier = p.outlier;
            gridCell.cell.value = JSON.stringify(gridCellValue);
          }
          return;
        }
      }
    }
  }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    ui.dialog({title: 'Edit chart'})
      .add(MultiCurveViewer.fromChartData(getChartData(gridCell)).root)
      .show();
  }

  renderCurves(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, data: IFitChartData): void {
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();

    const screenBounds = new DG.Rect(x, y, w, h).inflate(-6, -6);
    const [dataBox, xAxisBox, yAxisBox] = layoutChart(screenBounds);

    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);
    const minSize = Math.min(dataBox.width, dataBox.height);
    const ratio = minSize > 100 ? 1 : 0.2 + (minSize / 100) * 0.8;

    viewport.drawCoordinateGrid(g, xAxisBox, yAxisBox);

    for (const series of data.series!) {
      if (w < 70 || h < 45) {
        series.showPoints = '';
        if (data.chartOptions)
          data.chartOptions.showStatistics = [];
      }
      series.points.sort((a, b) => a.x - b.x);
      let userParamsFlag = true;
      const fitFunc = getSeriesFitFunction(series);
      let curve: (x: number) => number;
      if (series.parameters)
        curve = getCurve(series, fitFunc);
      else {
        const fitResult = fitSeries(series, fitFunc);
        curve = fitResult.fittedCurve;
        series.parameters = fitResult.parameters;
        userParamsFlag = false;
      }

      if (series.showPoints ?? 'points') {
        g.strokeStyle = series.pointColor ?? '0xFF40699c';
        if (series.showPoints === 'points')
          drawPoints(g, series, viewport, ratio);
        else if (['candlesticks', 'both'].includes(series.showPoints!))
          drawCandles(g, series, viewport, ratio);
      }

      if (series.showFitLine ?? true) {
        g.strokeStyle = series.fitLineColor ?? 'black';
        g.lineWidth = 2 * ratio;

        g.beginPath();
        for (let i = 0; i <= screenBounds.width; i++) {
          const x = screenBounds.x + i;
          const y = viewport.yToScreen(curve(viewport.xToWorld(x)));
          if (i === 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();
      }

      if ((series.showFitLine ?? true) && (series.showCurveConfidenceInterval ?? true)) {
        g.strokeStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_STROKE_COLOR;
        g.fillStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_FILL_COLOR;

        const confidenceIntervals = getSeriesConfidenceInterval(series, fitFunc, userParamsFlag);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM);
        fillConfidenceInterval(g, confidenceIntervals, screenBounds, viewport);
      }


      if (data.chartOptions?.showStatistics) {
        const statistics = getSeriesStatistics(series, fitFunc);
        for (let i = 0; i < data.chartOptions.showStatistics.length; i++) {
          const statName = data.chartOptions.showStatistics[i];
          const prop = statisticsProperties.find(p => p.name === statName);
          if (prop) {
            const s = StringUtils.formatNumber(prop.get(statistics));
            g.fillStyle = series.fitLineColor ?? 'black';
            g.textAlign = 'left';
            g.fillText(prop.name + ': ' + s, dataBox.x + 5, dataBox.y + 20 + 20 * i);
          }
        }
      }
    }
    g.restore();
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    if (w < 20 || h < 10) return;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    this.renderCurves(g, x, y, w, h, data);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.bounds.width < 50) {
      const canvas = ui.canvas(300, 200);
      this.render(canvas.getContext('2d')!, 0, 0, 300, 200, gridCell, null as any);
      const content = ui.divV([canvas]);
      ui.tooltip.show(content, e.x, e.y);
    }

    // TODO: add caching
    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    const screenBounds = gridCell.bounds.inflate(-6, -6);
    const dataBox = layoutChart(screenBounds)[0];
    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);

    for (let i = 0; i < data.series?.length!; i++) {
      if (!data.series![i].clickToToggle || data.series![i].showPoints === 'candlesticks')
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        const screenX = viewport.xToScreen(p.x);
        const screenY = viewport.yToScreen(p.y);
        const pxPerMarkerType = ((p.outlier ? 6 : 4) / 2) + OUTLIER_HITBOX_RADIUS;
        if (e.offsetX >= screenX - pxPerMarkerType && e.offsetX <= screenX + pxPerMarkerType &&
          e.offsetY >= screenY - pxPerMarkerType && e.offsetY <= screenY + pxPerMarkerType) {
          document.body.style.cursor = 'pointer';
          return;
        }
      }
    }
    document.body.style.cursor = 'default';
  }
}

const sample: IFitChartData = {
  // chartOptions could be retrieved either from the column, or from the cell
  'chartOptions': {
    'minX': 0, 'minY': 0, 'maxX': 5, 'maxY': 10,
    'xAxisName': 'concentration',
    'yAxisName': 'activity',
    'logX': false,
    'logY': false,
  },
  // These options are used as default options for the series. They could be overridden in series.
  'seriesOptions': {
    'fitFunction': 'sigmoid',
    // parameters not specified -> auto-fitting by default
    'pointColor': 'blue',
    'fitLineColor': 'red',
    'clickToToggle': true,
    'showPoints': 'points',
    'showFitLine': true,
    'showCurveConfidenceInterval': true,
  },
  'series': [
    {
      'fitFunction': 'sigmoid',
      // parameters specified -> use them, no autofitting
      'parameters': [1.86011e-07, -0.900, 103.748, -0.001],
      'points': [
        {'x': 0, 'y': 0},
        {'x': 1, 'y': 0.5},
        {'x': 2, 'y': 1},
        {'x': 3, 'y': 10, 'outlier': true},
        {'x': 4, 'y': 0},
      ],
    },
  ],
};
