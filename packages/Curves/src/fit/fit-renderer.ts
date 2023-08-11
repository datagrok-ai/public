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
  FitStatistics,
  fitChartDataProperties,
  fitSeriesProperties,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {BoxPlotStatistics, calculateBoxPlotStatistics} from '@datagrok-libraries/statistics/src/box-plot-statistics';
import {Viewport} from '@datagrok-libraries/utils/src/transform';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';

import {
  fitSeries,
  createDefaultChartData,
  getChartBounds,
  getSeriesFitFunction,
  getSeriesConfidenceInterval,
  getSeriesStatistics,
  getCurve,
  getColumnChartOptions,
  LogOptions
} from '@datagrok-libraries/statistics/src/fit/fit-data';

import {convertXMLToIFitChartData} from './fit-parser';
import {MultiCurveViewer} from './multi-curve-viewer';


export const TAG_FIT_CHART_FORMAT = '.fitChartFormat';
export const TAG_FIT_CHART_FORMAT_3DX = '3dx';
const MIN_CELL_RENDERER_PX_WIDTH = 20;
const MIN_CELL_RENDERER_PX_HEIGHT = 10;
const MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH = 70;
const MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT = 45;
const OUTLIER_PX_SIZE = 12;
const POINT_PX_SIZE = 4;
const OUTLIER_HITBOX_RADIUS = 2;
const MIN_AXES_CELL_PX_WIDTH = 70;
const MIN_AXES_CELL_PX_HEIGHT = 55;
const MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH = 150;
const MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT = 100;
const AXES_LEFT_PX_MARGIN = 30;
const AXES_TOP_PX_MARGIN = 5;
const AXES_BOTTOM_PX_MARGIN = 15;
const CANDLESTICK_BORDER_PX_SIZE = 4;
const CANDLESTICK_MEDIAN_PX_SIZE = 3.5;
const CANDLESTICK_OUTLIER_PX_SIZE = 6;


/** Merges properties of the two objects by iterating over the specified {@link properties}
 * and assigning properties from {@link source} to {@link target} only when
 * the property is not defined in target and is defined in source. */
export function mergeProperties(properties: DG.Property[], source: any, target: any): void {
  if (!source || !target)
    return;

  for (const p of properties) {
    if (!(p.name in target) && p.name in source)
      target[p.name] = source[p.name];
  }
}

/** Constructs {@link IFitChartData} from the grid cell, taking into account
 * chart and fit settings potentially defined on the dataframe and column level. */
export function getChartData(gridCell: DG.GridCell): IFitChartData {
  const cellChartData: IFitChartData = gridCell.cell?.column?.type === DG.TYPE.STRING ?
    (gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX ?
    convertXMLToIFitChartData(gridCell.cell.value) :
    JSON.parse(gridCell.cell.value ?? '{}') ?? {}) : createDefaultChartData();

  const columnChartOptions = getColumnChartOptions(gridCell.gridColumn);

  cellChartData.series ??= [];
  cellChartData.chartOptions ??= columnChartOptions.chartOptions;

  // merge cell options with column options
  mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions, cellChartData.chartOptions);
  for (const series of cellChartData.series)
    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, series);

  return cellChartData;
}

/** Performs a chart layout, returning [viewport, xAxis, yAxis] */
export function layoutChart(rect: DG.Rect): [DG.Rect, DG.Rect?, DG.Rect?] {
  if (rect.width < MIN_AXES_CELL_PX_WIDTH || rect.height < MIN_AXES_CELL_PX_HEIGHT)
    return [rect, undefined, undefined];
  return [
    rect.cutLeft(AXES_LEFT_PX_MARGIN).cutBottom(AXES_BOTTOM_PX_MARGIN).cutTop(AXES_TOP_PX_MARGIN),
    rect.getBottom(AXES_BOTTOM_PX_MARGIN).getTop(AXES_TOP_PX_MARGIN).cutLeft(AXES_LEFT_PX_MARGIN),
    rect.getLeft(AXES_LEFT_PX_MARGIN).cutBottom(AXES_BOTTOM_PX_MARGIN).cutTop(AXES_TOP_PX_MARGIN)
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
    markerColor, CANDLESTICK_MEDIAN_PX_SIZE * ratio);
}

/** Performs points drawing */
function drawPoints(g: CanvasRenderingContext2D, series: IFitSeries,
  transform: Viewport, ratio: number, logOptions: LogOptions, pointColor: number): void {
  for (let i = 0; i < series.points.length!; i++) {
    const p = series.points[i];
    const color = p.outlier ? DG.Color.red : pointColor;
    DG.Paint.marker(g,
      p.outlier ? DG.MARKER_TYPE.OUTLIER : (series.markerType as DG.MARKER_TYPE),
      transform.xToScreen(p.x), transform.yToScreen(p.y), color,
      (p.outlier ? OUTLIER_PX_SIZE : POINT_PX_SIZE) * ratio);
  }
}

/** Performs candles drawing */
function drawCandles(g: CanvasRenderingContext2D, series: IFitSeries,
  transform: Viewport, ratio: number, markerColor: number) : void {
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
      drawCandlestick(g, p.x, boxPlotStats, transform, ratio, markerColor);
      g.stroke();

      if (series.showPoints === 'both') {
        for (let ind = 0; ind < values.length; ind++) {
          if (values[ind] < boxPlotStats.lowerAdjacentValue || values[ind] > boxPlotStats.upperAdjacentValue) {
            DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER,
              transform.xToScreen(p.x), transform.yToScreen(values[ind]),
              series.pointColor ? DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker,
              CANDLESTICK_OUTLIER_PX_SIZE * ratio);
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
  for (let i = AXES_LEFT_PX_MARGIN; i <= screenBounds.width; i++) {
    const x = screenBounds.x + i;
    const y = confidenceType === CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
      transform.yToScreen(confIntervals.confidenceTop(transform.xToWorld(x))) :
      transform.yToScreen(confIntervals.confidenceBottom(transform.xToWorld(x)));
    if (i === AXES_LEFT_PX_MARGIN)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }
  g.stroke();
}

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, confIntervals: FitConfidenceIntervals,
  screenBounds: DG.Rect, transform: Viewport): void {
  g.beginPath();
  for (let i = AXES_LEFT_PX_MARGIN; i <= screenBounds.width; i++) {
    const x = screenBounds.x + i;
    const y = transform.yToScreen(confIntervals.confidenceTop(transform.xToWorld(x)));
    if (i === AXES_LEFT_PX_MARGIN)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }

  // reverse traverse to make a shape of confidence interval to fill it
  for (let i = screenBounds.width; i >= AXES_LEFT_PX_MARGIN; i--) {
    const x = screenBounds.x + i;
    const y = transform.yToScreen(confIntervals.confidenceBottom(transform.xToWorld(x)));
    g.lineTo(x, y);
  }
  g.closePath();
  g.fill();
}

/** Performs a dropline drawing */
function drawDropline(g: CanvasRenderingContext2D, transform: Viewport, xValue: number, dataBounds: DG.Rect,
  curve: (x: number) => number): void {
  g.moveTo(transform.xToScreen(dataBounds.minX), transform.yToScreen(curve(xValue)));
  g.lineTo(transform.xToScreen(xValue), transform.yToScreen(curve(xValue)));
  g.lineTo(transform.xToScreen(xValue), transform.yToScreen(dataBounds.minY));
}

export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return FIT_CELL_TYPE; }

  get cellType() { return FIT_CELL_TYPE; }

  getDefaultSize(gridColumn: DG.GridColumn): {width?: number | null, height?: number | null} {
    return {width: 160, height: 120};
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell.value)
      return;

    grok.shell.o = gridCell;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    const screenBounds = gridCell.bounds.inflate(-6, -6);
    const dataBox = layoutChart(screenBounds)[0];
    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);

    for (let i = 0; i < data.series?.length!; i++) {
      if (!data.series![i].clickToToggle || data.series![i].showPoints !== 'points')
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        const screenX = viewport.xToScreen(p.x);
        const screenY = viewport.yToScreen(p.y);
        const pxPerMarkerType = ((p.outlier ? OUTLIER_PX_SIZE : POINT_PX_SIZE) / 2) + OUTLIER_HITBOX_RADIUS;
        if (e.offsetX >= screenX - pxPerMarkerType && e.offsetX <= screenX + pxPerMarkerType &&
          e.offsetY >= screenY - pxPerMarkerType && e.offsetY <= screenY + pxPerMarkerType) {
          p.outlier = !p.outlier;
          const columns = gridCell.grid.dataFrame.columns.byTags({'.sourceColumn':
            gridCell.cell.column.name, '.seriesNumber': i});
          if (columns) {
            const stats = getSeriesStatistics(data.series![i], getSeriesFitFunction(data.series![i]));
            for (const column of columns) {
              column.set(gridCell.cell.rowIndex, stats[column.name as keyof FitStatistics]);
            }
          }
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
    if (!gridCell.cell.value)
      return;

    ui.dialog({title: 'Edit chart'})
      .add(MultiCurveViewer.fromChartData(getChartData(gridCell)).root)
      .show({resizable: true});
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
    const chartLogOptions: LogOptions = {logX: data.chartOptions?.logX, logY: data.chartOptions?.logY};

    viewport.drawCoordinateGrid(g, xAxisBox, yAxisBox);

    for (let i = 0; i < data.series?.length!; i++) {
      const series = data.series![i];
      if (w < MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH || h < MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT) {
        series.showPoints = '';
        if (data.chartOptions)
          data.chartOptions.showStatistics = [];
      }
      series.points.sort((a, b) => a.x - b.x);
      let userParamsFlag = true;
      const fitFunc = getSeriesFitFunction(series);
      let curve: (x: number) => number;
      if (series.parameters) {
        if (data.chartOptions?.logX)
          series.parameters[2] = Math.log10(series.parameters[2]);
        curve = getCurve(series, fitFunc);
      }
      else {
        const fitResult = fitSeries(series, fitFunc, chartLogOptions);
        curve = fitResult.fittedCurve;
        series.parameters = fitResult.parameters;
        userParamsFlag = false;
      }
  
      if (series.showPoints ?? 'points') {
        const pointColor = series.pointColor ? DG.Color.fromHtml(series.pointColor) ?
          series.pointColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        g.strokeStyle = pointColor;
        if (series.showPoints === 'points')
          drawPoints(g, series, viewport, ratio, chartLogOptions, DG.Color.fromHtml(pointColor));
        else if (['candlesticks', 'both'].includes(series.showPoints!))
          drawCandles(g, series, viewport, ratio, DG.Color.fromHtml(pointColor));
      }
  
      if (series.showFitLine ?? true) {
        const lineColor = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ?
          series.fitLineColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        g.strokeStyle = lineColor;
        g.lineWidth = 2 * ratio;
  
        g.beginPath();
        for (let j = AXES_LEFT_PX_MARGIN; j <= screenBounds.width; j++) {
          const x = screenBounds.x + j;
          const xForY = data.chartOptions?.logX ? Math.log10(viewport.xToWorld(x)) : viewport.xToWorld(x);
          const y = data.chartOptions?.logY ? viewport.yToScreen(Math.pow(10, curve(xForY))) : viewport.yToScreen(curve(xForY));
            viewport.yToScreen(curve(viewport.xToWorld(x)));
          if (j === AXES_LEFT_PX_MARGIN)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();
      }
  
      if ((series.showFitLine ?? true) && (series.showCurveConfidenceInterval ?? false)) {
        g.strokeStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_STROKE_COLOR;
        g.fillStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_FILL_COLOR;
  
        const confidenceIntervals = getSeriesConfidenceInterval(series, fitFunc, userParamsFlag);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM);
        fillConfidenceInterval(g, confidenceIntervals, screenBounds, viewport);
      }
  
      if (series.droplines) {
        g.save();
        g.strokeStyle = 'blue';
        g.lineWidth = ratio;
        g.beginPath();
        g.setLineDash([5, 5]);
        for (let j = 0; j < series.droplines.length; j++) {
          const droplineName = series.droplines[j];
          if (droplineName === 'IC50')
            drawDropline(g, viewport, series.parameters[2], dataBounds, curve);
        }
        g.stroke();
        g.restore();
      }

      if (data.chartOptions?.showStatistics) {
        const statistics = getSeriesStatistics(series, fitFunc);
        for (let j = 0; j < data.chartOptions.showStatistics.length; j++) {
          const statName = data.chartOptions.showStatistics[j];
          const prop = statisticsProperties.find(p => p.name === statName);
          if (prop) {
            const s = StringUtils.formatNumber(prop.get(statistics));
            const statColor = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ?
              series.fitLineColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
            g.fillStyle = statColor;
            g.textAlign = 'left';
            g.fillText(prop.name + ': ' + s, dataBox.x + 5, dataBox.y + 20 + 20 * j);
          }
        }
      }
    }

    g.restore();
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    if (w < MIN_CELL_RENDERER_PX_WIDTH || h < MIN_CELL_RENDERER_PX_HEIGHT)
      return;

    if (!gridCell.cell.value)
      return;
    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    this.renderCurves(g, x, y, w, h, data);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell.value)
      return;

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
      if (!data.series![i].clickToToggle || data.series![i].showPoints !== 'points')
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        const screenX = viewport.xToScreen(p.x);
        const screenY = viewport.yToScreen(p.y);
        const pxPerMarkerType = ((p.outlier ? OUTLIER_PX_SIZE : POINT_PX_SIZE) / 2) + OUTLIER_HITBOX_RADIUS;
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
