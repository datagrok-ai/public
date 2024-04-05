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
  LogOptions,
  getDataFrameChartOptions
} from '@datagrok-libraries/statistics/src/fit/fit-data';

import {convertXMLToIFitChartData} from './fit-parser';
import {CellRenderViewer} from './cell-render-viewer';
import { calculateSeriesStats, getChartDataAggrStats } from './fit-grid-cell-handler';


export const TAG_FIT_CHART_FORMAT = '.fitChartFormat';
export const TAG_FIT_CHART_FORMAT_3DX = '3dx';
const MIN_CELL_RENDERER_PX_WIDTH = 20;
const MIN_CELL_RENDERER_PX_HEIGHT = 10;
const MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH = 70;
const MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT = 45;
const OUTLIER_PX_SIZE = 12;
const POINT_PX_SIZE = 4;
const OUTLIER_HITBOX_RADIUS = 2;
const MIN_AXES_CELL_PX_WIDTH = 100;
const MIN_AXES_CELL_PX_HEIGHT = 55;
const MIN_TITLE_PX_WIDTH = 275;
const MIN_TITLE_PX_HEIGHT = 225;
const MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH = 180;
const MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT = 140;
const MIN_LEGEND_PX_WIDTH = 325;
const MIN_LEGEND_PX_HEIGHT = 275;
const MIN_DROPLINES_VISIBILITY_PX_WIDTH = 120;
const MIN_DROPLINES_VISIBILITY_PX_HEIGHT = 110;
const AXES_LEFT_PX_MARGIN = 38;
const AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS = 45;
const AXES_TOP_PX_MARGIN = 5;
const AXES_TOP_PX_MARGIN_WITH_TITLE = 25;
const AXES_RIGHT_PX_MARGIN = 18;
const AXES_BOTTOM_PX_MARGIN = 15;
const AXES_BOTTOM_PX_MARGIN_WITH_AXES_LABELS = 30;
const CANDLESTICK_BORDER_PX_SIZE = 4;
const CANDLESTICK_MEDIAN_PX_SIZE = 3.5;
const CANDLESTICK_OUTLIER_PX_SIZE = 6;
export const INFLATE_SIZE = -12;
export const LINE_STYLES: {[key: string]: number[]} = {
  'solid': [],
  'dotted': [1, 1],
  'dashed': [5, 5],
  'dashdotted': [5, 5, 2, 5],
};


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

export function mergeSeries(series: IFitSeries[]): IFitSeries | null {
  if (series.length === 0)
    return null;
  const mergedSeries: IFitSeries = {
    points: [],
    name: series[0].name,
    fitFunction: series[0].fitFunction,
    markerType: series[0].markerType,
    lineStyle: series[0].lineStyle,
    pointColor: series[0].pointColor,
    fitLineColor: series[0].fitLineColor,
    confidenceIntervalColor: series[0].confidenceIntervalColor,
    outlierColor: series[0].outlierColor,
    connectDots: series[0].connectDots,
    showFitLine: series[0].showFitLine,
    showPoints: series[0].showPoints,
    showCurveConfidenceInterval: series[0].showCurveConfidenceInterval,
    errorModel: series[0].errorModel,
    clickToToggle: series[0].clickToToggle,
    labels: series[0].labels,
    droplines: series[0].droplines,
    columnName: series[0].columnName,
  };
  for (const s of series)
    mergedSeries.points = [...mergedSeries.points, ...s.points];
  return mergedSeries;
}

/** Constructs {@link IFitChartData} from the grid cell, taking into account
 * chart and fit settings potentially defined on the dataframe and column level. */
export function getChartData(gridCell: DG.GridCell): IFitChartData {
  // removing '|' from JSON (how did it get here?)
  let cellValue = gridCell.cell.value as string;
  if (cellValue.includes('|'))
    cellValue = cellValue.replaceAll('|', '');
  const cellChartData: IFitChartData = gridCell.cell?.column?.type === DG.TYPE.STRING ?
    (gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX ?
    convertXMLToIFitChartData(cellValue) :
    JSON.parse(cellValue ?? '{}') ?? {}) : createDefaultChartData();

  const columnChartOptions = gridCell.cell.column ? getColumnChartOptions(gridCell.cell.column) : {};
  const dfChartOptions = gridCell.cell.column ? getDataFrameChartOptions(gridCell.cell.dataFrame) : {};

  cellChartData.series ??= [];
  cellChartData.chartOptions ??= columnChartOptions.chartOptions;

  // merge cell options with column options
  mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions, cellChartData.chartOptions);
  mergeProperties(fitChartDataProperties, dfChartOptions.chartOptions, cellChartData.chartOptions);
  for (const series of cellChartData.series) {
    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, series);
    mergeProperties(fitSeriesProperties, dfChartOptions.seriesOptions, series);
  }

  return cellChartData;
}

/** Performs a chart layout, returning [viewport, xAxis, yAxis] */
export function layoutChart(rect: DG.Rect, showAxesLabels: boolean, showTitle: boolean): [DG.Rect, DG.Rect?, DG.Rect?] {
  if (rect.width < MIN_AXES_CELL_PX_WIDTH || rect.height < MIN_AXES_CELL_PX_HEIGHT)
    return [rect, undefined, undefined];
  const axesLeftPxMargin = showAxesLabels ? AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : AXES_LEFT_PX_MARGIN;
  const axesBottomPxMargin = showAxesLabels ? AXES_BOTTOM_PX_MARGIN_WITH_AXES_LABELS : AXES_BOTTOM_PX_MARGIN;
  const axesTopPxMargin = showTitle ? AXES_TOP_PX_MARGIN_WITH_TITLE : AXES_TOP_PX_MARGIN;
  return [
    rect.cutLeft(axesLeftPxMargin).cutBottom(axesBottomPxMargin).cutTop(axesTopPxMargin).cutRight(AXES_RIGHT_PX_MARGIN),
    rect.getBottom(axesBottomPxMargin).getTop(axesTopPxMargin).cutLeft(axesLeftPxMargin).cutRight(AXES_RIGHT_PX_MARGIN),
    rect.getLeft(axesLeftPxMargin).getRight(AXES_RIGHT_PX_MARGIN).cutBottom(axesBottomPxMargin).cutTop(axesTopPxMargin)
  ];
}

/** Checks if the color is valid */
export function isColorValid(color: string | null | undefined): boolean {
  if (color === undefined || color === null || color === '')
    return false;
  return DG.Color.fromHtml(color) !== undefined;
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
  transform: Viewport, ratio: number, logOptions: LogOptions, pointColor: number, outlierColor: number): void {
  for (let i = 0; i < series.points.length!; i++) {
    const p = series.points[i];
    const color = !series.connectDots ? p.outlier ? (p.outlierColor ? DG.Color.fromHtml(p.outlierColor) ?
      DG.Color.fromHtml(p.outlierColor) : outlierColor : outlierColor) : p.color ? DG.Color.fromHtml(p.color) ?
      DG.Color.fromHtml(p.color) : pointColor : pointColor : pointColor;
    const marker = p.marker ? p.marker as DG.MARKER_TYPE : series.markerType as DG.MARKER_TYPE;
    const size = !series.connectDots ? p.outlier ? OUTLIER_PX_SIZE * ratio : p.size ? p.size : POINT_PX_SIZE * ratio : POINT_PX_SIZE * ratio;
    DG.Paint.marker(g,
      !series.connectDots ? p.outlier ? DG.MARKER_TYPE.OUTLIER : marker : marker,
      transform.xToScreen(p.x), transform.yToScreen(p.y), color, size);
    if (p.stdev && !p.outlier) {
      g.strokeStyle = DG.Color.toHtml(color);
      g.beginPath();
      g.moveTo(transform.xToScreen(p.x), transform.yToScreen(p.y + p.stdev));
      g.lineTo(transform.xToScreen(p.x), transform.yToScreen(p.y - p.stdev));
      g.stroke();
    }
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
function drawConfidenceInterval(g: CanvasRenderingContext2D, confIntervals: FitConfidenceIntervals, screenBounds: DG.Rect,
  transform: Viewport, confidenceType: string, showAxes: boolean, showAxesLabels: boolean, logOptions: LogOptions): void {
  g.beginPath();
  const axesLeftPxMargin = showAxes ? showAxesLabels ? AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : AXES_LEFT_PX_MARGIN : 0;
  const axesRightPxMargin = showAxes ? AXES_RIGHT_PX_MARGIN : 0;
  const confIntervalFunc = confidenceType === CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
    confIntervals.confidenceTop : confIntervals.confidenceBottom;
  for (let i = axesLeftPxMargin; i <= screenBounds.width - axesRightPxMargin; i++) {
    const x = screenBounds.x + i;
    const xForY = logOptions.logX ? Math.log10(transform.xToWorld(x)) : transform.xToWorld(x);
    const y = logOptions.logY ? transform.yToScreen(Math.pow(10, confIntervalFunc(xForY))) : transform.yToScreen(confIntervalFunc(xForY));
    if (i === axesLeftPxMargin)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }
  g.stroke();
}

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, confIntervals: FitConfidenceIntervals, screenBounds: DG.Rect,
  transform: Viewport, showAxes: boolean, showAxesLabels: boolean, logOptions: LogOptions): void {
  g.beginPath();
  const axesLeftPxMargin = showAxes ? showAxesLabels ? AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : AXES_LEFT_PX_MARGIN : 0;
  const axesRightPxMargin = showAxes ? AXES_RIGHT_PX_MARGIN : 0;
  for (let i = axesLeftPxMargin; i <= screenBounds.width - axesRightPxMargin; i++) {
    const x = screenBounds.x + i;
    const xForY = logOptions.logX ? Math.log10(transform.xToWorld(x)) : transform.xToWorld(x);
    const y = logOptions.logY ? transform.yToScreen(Math.pow(10, confIntervals.confidenceTop(xForY))) :
      transform.yToScreen(confIntervals.confidenceTop(xForY));
    if (i === axesLeftPxMargin)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }

  // reverse traverse to make a shape of confidence interval to fill it
  for (let i = screenBounds.width - axesRightPxMargin; i >= axesLeftPxMargin; i--) {
    const x = screenBounds.x + i;
    const xForY = logOptions.logX ? Math.log10(transform.xToWorld(x)) : transform.xToWorld(x);
    const y = logOptions.logY ? transform.yToScreen(Math.pow(10, confIntervals.confidenceBottom(xForY))) :
      transform.yToScreen(confIntervals.confidenceBottom(xForY));
    g.lineTo(x, y);
  }
  g.closePath();
  g.fill();
}

/** Performs a dropline drawing */
function drawDropline(g: CanvasRenderingContext2D, transform: Viewport, xValue: number, dataBounds: DG.Rect,
  curve: (x: number) => number, logOptions: LogOptions): void {
  if (logOptions.logX)
    xValue = Math.pow(10, xValue);
  const xForY = logOptions.logX ? Math.log10(xValue) : xValue;
  const y = logOptions.logY ? Math.pow(10, curve(xForY)) : curve(xForY);
  g.moveTo(transform.xToScreen(dataBounds.minX), transform.yToScreen(y));
  g.lineTo(transform.xToScreen(xValue), transform.yToScreen(y));
  g.lineTo(transform.xToScreen(xValue), transform.yToScreen(dataBounds.minY));
}

/** Performs x zeroes substitution if log x */
export function substituteZeroes(data: IFitChartData): void {
  for (let i = 0; i < data.series?.length!; i++) {
    const series = data.series![i];
    if (series.points.every((p) => p.x !== 0))
      continue;
    let minNonZeroX = Number.MAX_VALUE;
    let maxNonZeroX = 0;
    let countOfDistNonZeroX = 0;
    const uniqueArr: number[] = [];
    for (let j = 0; j < series.points.length; j++) {
      if (series.points[j].x < minNonZeroX && series.points[j].x !== 0)
        minNonZeroX = series.points[j].x;
      if (series.points[j].x > maxNonZeroX && series.points[j].x !== 0)
        maxNonZeroX = series.points[j].x;
      if (!uniqueArr.includes(series.points[j].x)) {
        uniqueArr[uniqueArr.length] = series.points[j].x;
        countOfDistNonZeroX++;
      }
    }
    const zeroSubstitute = Math.pow(10, Math.log10(minNonZeroX) - (Math.log10(maxNonZeroX) - Math.log10(minNonZeroX) / (countOfDistNonZeroX - 1)));
    for (let j = 0; j < series.points.length; j++) {
      if (series.points[j].x === 0)
        series.points[j].x = zeroSubstitute;
    }
  }
}

export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return FIT_CELL_TYPE; }

  get cellType() { return FIT_CELL_TYPE; }

  getDefaultSize(gridColumn: DG.GridColumn): {width?: number | null, height?: number | null} {
    return {width: 220, height: 150};
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell.value)
      return;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    if (data.series?.some((series) => series.points.length === 0))
      return;

    grok.shell.o = gridCell;

    const screenBounds = gridCell.bounds.inflate(INFLATE_SIZE / 2, INFLATE_SIZE / 2);
    const showAxesLabels = gridCell.bounds.width >= MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH && gridCell.bounds.height >= MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT;
    const showTitle = gridCell.bounds.width >= MIN_TITLE_PX_WIDTH && gridCell.bounds.height >= MIN_TITLE_PX_HEIGHT;
    const dataBox = layoutChart(screenBounds, showAxesLabels, showTitle)[0];
    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);

    for (let i = 0; i < data.series?.length!; i++) {
      if (data.series![i].connectDots || !data.series![i].clickToToggle || data.series![i].showPoints !== 'points' ||
        screenBounds.width < MIN_AXES_CELL_PX_WIDTH || screenBounds.height < MIN_AXES_CELL_PX_HEIGHT)
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        const screenX = viewport.xToScreen(p.x);
        const screenY = viewport.yToScreen(p.y);
        const pxPerMarkerType = ((p.outlier ? OUTLIER_PX_SIZE : POINT_PX_SIZE) / 2) + OUTLIER_HITBOX_RADIUS;
        if (e.offsetX >= screenX - pxPerMarkerType && e.offsetX <= screenX + pxPerMarkerType &&
          e.offsetY >= screenY - pxPerMarkerType && e.offsetY <= screenY + pxPerMarkerType) {
          p.outlier = !p.outlier;
          const columns = gridCell.grid.dataFrame.columns.byTags({'.sourceColumn': gridCell.cell.column.name});
          if (columns) {
            for (const column of columns) {
              const chartLogOptions: LogOptions = {logX: data.chartOptions?.logX, logY: data.chartOptions?.logY};
              const stats = column.tags['.seriesAggregation'] !== null ?
                getChartDataAggrStats(data, column.tags['.seriesAggregation']) :
                column.tags['.seriesNumber'] === i ? calculateSeriesStats(data.series![i], chartLogOptions) : null;
              if (stats === null)
                continue;
              column.set(gridCell.cell.rowIndex, stats[column.tags['.statistics'] as keyof FitStatistics]);  
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

    const cellRenderViewer = CellRenderViewer.fromGridCell(gridCell);
    const dlg = ui.dialog({title: 'Edit chart'})
      .add(cellRenderViewer.root)
      .show({resizable: true});

    // canvas is created as (300, 150), so we change its size to the dialog contents box size
    const dlgContentsBox = dlg.root.getElementsByClassName('d4-dialog-contents dlg-edit-chart')[0].firstChild as HTMLElement;
    cellRenderViewer.canvas.width = dlgContentsBox.clientWidth;
    cellRenderViewer.canvas.height = dlgContentsBox.clientHeight;
    cellRenderViewer.render();

    // contents ui-box isn't resizable by default
    dlgContentsBox.style.width = '100%';
    dlgContentsBox.style.height = '100%';

    ui.tools.handleResize(dlgContentsBox, (w: number, h: number) => {
      cellRenderViewer.canvas.width = w;
      cellRenderViewer.canvas.height = h;
      cellRenderViewer.render();
    });
  }

  renderCurves(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, data: IFitChartData): void {
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();

    if (data.chartOptions?.allowXZeroes && data.chartOptions?.logX && data.series?.some((series) => series.points.some((p) => p.x === 0)))
      substituteZeroes(data);
    const screenBounds = new DG.Rect(x, y, w, h).inflate(INFLATE_SIZE / 2, INFLATE_SIZE / 2);
    const showAxes = screenBounds.width >= MIN_AXES_CELL_PX_WIDTH && screenBounds.height >= MIN_AXES_CELL_PX_HEIGHT;
    // TODO: make bigger sizes
    const showAxesLabels = w >= MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH && h >= MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT
      && !!data.chartOptions?.xAxisName && !!data.chartOptions.yAxisName;
    const showTitle = w >= MIN_TITLE_PX_WIDTH && h >= MIN_TITLE_PX_HEIGHT && !!data.chartOptions?.title;
    const showLegend = w >= MIN_LEGEND_PX_WIDTH && h >= MIN_LEGEND_PX_HEIGHT;
    const showDroplines =  w >= MIN_DROPLINES_VISIBILITY_PX_WIDTH && h >= MIN_DROPLINES_VISIBILITY_PX_HEIGHT;
    const [dataBox, xAxisBox, yAxisBox] = layoutChart(screenBounds, showAxesLabels, showTitle);

    const dataBounds = getChartBounds(data);
    if ((dataBounds.x < 0 && data.chartOptions) || (dataBounds.x === 0 && data.chartOptions && !data.chartOptions.allowXZeroes))
      data.chartOptions.logX = false;
    if (dataBounds.y <= 0 && data.chartOptions) data.chartOptions.logY = false;
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);
    const minSize = Math.min(dataBox.width, dataBox.height);
    // TODO: make thinner
    const ratio = minSize > 100 ? 1 : 0.2 + (minSize / 100) * 0.8;
    const chartLogOptions: LogOptions = {logX: data.chartOptions?.logX, logY: data.chartOptions?.logY};

    g.save();
    g.font = '11px Roboto, "Roboto Local"';
    viewport.drawCoordinateGrid(g, xAxisBox, yAxisBox);
    g.restore();

    const mergedSeries = data.chartOptions?.mergeSeries ? mergeSeries(data.series!) : null;
    if (data.chartOptions?.mergeSeries && mergedSeries === null)
      return;
    for (let i = 0; i < (data.chartOptions?.mergeSeries ? 1 : data.series?.length!); i++) {
      const series = mergedSeries ?? data.series![i];
      if (series.points.some((point) => point.x === undefined || point.y === undefined))
        continue;
      if (w < MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH || h < MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT) {
        series.showPoints = '';
        if (data.chartOptions)
          data.chartOptions.showStatistics = [];
      }
      series.points.sort((a, b) => a.x - b.x);

      if (series.connectDots ?? false) {
        g.strokeStyle = series.pointColor ? DG.Color.fromHtml(series.pointColor) ?
          series.pointColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        g.lineWidth = 2 * ratio;
        g.beginPath();
        for (let j = 0; j < series.points.length; j++) {
          const x = series.points[j].x;
          const y = series.points[j].y;
          const screenX = viewport.xToScreen(x);
          const screenY = viewport.yToScreen(y);
          if (j === 0)
            g.moveTo(screenX, screenY);
          else
            g.lineTo(screenX, screenY);
        }
        g.stroke();
      }

      let userParamsFlag = true;
      const fitFunc = getSeriesFitFunction(series);
      let curve: ((x: number) => number) | null = null;
      if (!(series.connectDots && !series.showFitLine)) {
        if (series.parameters) {
          if (data.chartOptions?.logX) {
            if (series.parameters[2] > 0)
              series.parameters[2] = Math.log10(series.parameters[2]);
          }
          curve = getCurve(series, fitFunc);
        }
        else {
          const fitResult = fitSeries(series, fitFunc, chartLogOptions);
          curve = fitResult.fittedCurve;
          series.parameters = fitResult.parameters;
          userParamsFlag = false;
        }
      }

      if (series.showPoints ?? 'points') {
        const pointColor = series.pointColor ? DG.Color.fromHtml(series.pointColor) ?
          series.pointColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        const outlierColor = series.outlierColor ? DG.Color.fromHtml(series.outlierColor) ?
          DG.Color.fromHtml(series.outlierColor) : DG.Color.red : DG.Color.red;
        g.strokeStyle = pointColor;
        if (series.connectDots && series.showPoints != '')
          drawPoints(g, series, viewport, ratio, chartLogOptions, DG.Color.fromHtml(pointColor), outlierColor);
        else if (series.showPoints === 'points')
          drawPoints(g, series, viewport, ratio, chartLogOptions, DG.Color.fromHtml(pointColor), outlierColor);
        else if (['candlesticks', 'both'].includes(series.showPoints!))
          drawCandles(g, series, viewport, ratio, DG.Color.fromHtml(pointColor));
      }
  
      if (series.showFitLine ?? true) {
        const lineColor = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ?
          series.fitLineColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        g.save();
        g.strokeStyle = lineColor;
        g.lineWidth = 2 * ratio;
  
        g.beginPath();
        if (series.lineStyle)
          g.setLineDash(LINE_STYLES[series.lineStyle]);
        const axesLeftPxMargin = showAxes ? showAxesLabels ? AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : AXES_LEFT_PX_MARGIN : 0;
        const axesRightPxMargin = showAxes ? AXES_RIGHT_PX_MARGIN : 0;
        for (let j = axesLeftPxMargin; j <= screenBounds.width - axesRightPxMargin; j++) {
          const x = screenBounds.x + j;
          const xForY = data.chartOptions?.logX ? Math.log10(viewport.xToWorld(x)) : viewport.xToWorld(x);
          const y = data.chartOptions?.logY ? viewport.yToScreen(Math.pow(10, curve!(xForY))) : viewport.yToScreen(curve!(xForY));
          if (j === axesLeftPxMargin)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();
        g.restore();
      }
  
      if ((series.showFitLine ?? true) && (series.showCurveConfidenceInterval ?? false)) {
        g.strokeStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_STROKE_COLOR;
        g.fillStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_FILL_COLOR;
  
        const confidenceIntervals = getSeriesConfidenceInterval(series, fitFunc, userParamsFlag, chartLogOptions);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport,
          CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP, showAxes, showAxesLabels, chartLogOptions);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport,
          CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM, showAxes, showAxesLabels, chartLogOptions);
        fillConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, showAxes, showAxesLabels, chartLogOptions);
      }
  
      if ((series.showFitLine ?? true) && series.droplines && showDroplines) {
        g.save();
        g.strokeStyle = 'blue';
        g.lineWidth = ratio;
        g.beginPath();
        g.setLineDash([5, 5]);
        for (let j = 0; j < series.droplines.length; j++) {
          const droplineName = series.droplines[j];
          if (droplineName === 'IC50')
            drawDropline(g, viewport, series.parameters![2], dataBounds, curve!, chartLogOptions);
        }
        g.stroke();
        g.restore();
      }

      if ((series.showFitLine ?? true) && data.chartOptions?.showStatistics) {
        const statistics = getSeriesStatistics(series, fitFunc, chartLogOptions);
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

    if (showTitle) {
      g.font = '12px Roboto, "Roboto Local"';
      g.textAlign = 'center';
      g.fillStyle = 'black';
      g.fillText(data.chartOptions?.title!, dataBox.midX - 5, y + 15);
    }    

    if (showAxesLabels) {
      g.font = '11px Roboto, "Roboto Local"';
      g.textAlign = 'center';
      g.fillStyle = 'black';
      g.fillText(data.chartOptions?.xAxisName!, dataBox.midX - 5, y + h - 10);
      g.translate(x, y);
      g.rotate(-Math.PI / 2);
      const axesTopPxMargin = showTitle ? AXES_TOP_PX_MARGIN_WITH_TITLE : AXES_TOP_PX_MARGIN;
      g.fillText(data.chartOptions?.yAxisName!, -(dataBox.height / 2 + axesTopPxMargin + 15), 15);
      g.restore();
    }

    if (showLegend) {
      g.font = '11px Roboto, "Roboto Local"';
      for (let i = 0; i < data.series?.length!; i++) {
        g.beginPath();
        const series = data.series![i];
        const lineColor = series.fitLineColor ? DG.Color.fromHtml(series.fitLineColor) ?
          series.fitLineColor : DG.Color.toHtml(DG.Color.getCategoricalColor(i)) : DG.Color.toHtml(DG.Color.getCategoricalColor(i));
        const markerColor = series.pointColor ? DG.Color.fromHtml(series.pointColor) ?
          DG.Color.fromHtml(series.pointColor) : DG.Color.getCategoricalColor(i) : DG.Color.getCategoricalColor(i);
        g.strokeStyle = lineColor;
        g.lineWidth = 2 * ratio;
        const text = `${data.chartOptions?.showColumnLabel ? series.columnName ?? '' : ''} ${series.name ?? ''}`;
        const textWidth = g.measureText(text).width;
        g.moveTo(dataBox.maxX - textWidth - 25, dataBox.y + 20 - 4 + 20 * i);
        g.lineTo(dataBox.maxX - textWidth - 5, dataBox.y + 20 - 4 + 20 * i);
        const marker = series.markerType ? series.markerType as DG.MARKER_TYPE : DG.MARKER_TYPE.CIRCLE;
        DG.Paint.marker(g, marker, dataBox.maxX - textWidth - 25 + (20 / 2), dataBox.y + 20 - 4 + 20 * i, markerColor, POINT_PX_SIZE * ratio);
        g.fillStyle = 'black';
        g.fillText(text, dataBox.maxX - textWidth, dataBox.y + 20 + 20 * i);
        g.stroke();
      }
      g.restore();
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
    const data = gridCell.cell.column?.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    if (data.series?.some((series) => series.points.length === 0))
      return;

    if (data.series?.some((series) => series.points.every((point) => point.x === 0)))
      return;

    if (data.series?.some((series) => series.points.every((point) => point.y === 0)))
      return;

    data.series?.forEach((series) => series.columnName = gridCell.cell.column.name);

    this.renderCurves(g, x, y, w, h, data);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell.value)
      return;

    if (gridCell.bounds.width < MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH || gridCell.bounds.height < MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT) {
      const canvas = ui.canvas(300, 200);
      this.render(canvas.getContext('2d')!, 0, 0, 300, 200, gridCell, null as any);
      const content = ui.divV([canvas]);
      ui.tooltip.show(content, e.x, e.y);
    }

    // TODO: add caching
    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    const screenBounds = gridCell.bounds.inflate(INFLATE_SIZE / 2, INFLATE_SIZE / 2);
    if (screenBounds.width >= MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH && screenBounds.height >= MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT) {
      const showAxesLabels = gridCell.bounds.width >= MIN_X_AXIS_NAME_VISIBILITY_PX_WIDTH && gridCell.bounds.height >= MIN_Y_AXIS_NAME_VISIBILITY_PX_HEIGHT;
      const showTitle = gridCell.bounds.width >= MIN_TITLE_PX_WIDTH && gridCell.bounds.height >= MIN_TITLE_PX_HEIGHT;
      const dataBox = layoutChart(screenBounds, showAxesLabels, showTitle)[0];
      const dataBounds = getChartBounds(data);
      const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);

      for (let i = 0; i < data.series?.length!; i++) {
        if (data.series![i].showPoints !== 'points')
          continue;
        for (let j = 0; j < data.series![i].points.length!; j++) {
          const p = data.series![i].points[j];
          const screenX = viewport.xToScreen(p.x);
          const screenY = viewport.yToScreen(p.y);
          const pxPerMarkerType = ((p.outlier ? OUTLIER_PX_SIZE : POINT_PX_SIZE) / 2) + OUTLIER_HITBOX_RADIUS;
          if (e.offsetX >= screenX - pxPerMarkerType && e.offsetX <= screenX + pxPerMarkerType &&
              e.offsetY >= screenY - pxPerMarkerType && e.offsetY <= screenY + pxPerMarkerType) {
            ui.tooltip.show(ui.divV([ui.divText(`x: ${DG.format(p.x, '#0.000')}`),
              ui.divText(`y: ${DG.format(p.y, '#0.000')}`)]), e.x + 16, e.y + 16);
            if (!data.series![i].connectDots && data.series![i].clickToToggle && screenBounds.width >= MIN_AXES_CELL_PX_WIDTH &&
              screenBounds.height >= MIN_AXES_CELL_PX_HEIGHT)
              document.body.style.cursor = 'pointer';
            return;
          }
        }
      }
      ui.tooltip.hide();
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
