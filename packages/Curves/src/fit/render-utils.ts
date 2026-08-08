/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {
  FitConfidenceIntervals,
  IFitChartData,
  IFitSeries,
  LogOptions,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {
  getSeriesConfidenceInterval,
  getSeriesFit,
  toDataSpace,
} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {Viewport} from '@datagrok-libraries/utils/src/transform';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {BoxPlotStatistics, calculateBoxPlotStatistics} from '@datagrok-libraries/statistics/src/box-plot-statistics';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';
import {FitFunction, getStatistic, getStatisticProperty} from '@datagrok-libraries/statistics/src/fit/fit-engine';


export enum ColorType {
  POINT = 'pointColor',
  OUTLIER = 'outlierColor',
  FIT_LINE = 'fitLineColor',
}
export type SeriesColorType = `${ColorType}`;

interface FitRenderOptions {
    viewport: Viewport;
    screenBounds?: DG.Rect;
    ratio?: number;
    seriesIdx?: number;
}

interface FitPointRenderOptions extends FitRenderOptions {
    x: number;
    color: string;
    boxPlotStats: BoxPlotStatistics;
}

interface FitLineRenderOptions extends FitRenderOptions {
    logOptions: LogOptions;
    showAxes: boolean;
    showAxesLabels: boolean;
    screenBounds: DG.Rect;
    curveFunc?: (x: number) => number;
}

interface FitConfidenceIntervalRenderOptions extends FitLineRenderOptions {
    fitFunc?: FitFunction;
    userParamsFlag?: boolean;
    confidenceIntervals?: FitConfidenceIntervals;
    confidenceType?: string;
    dataPoints?: {x: number[], y: number[]};
}

interface FitDroplineRenderOptions extends FitRenderOptions {
    showDroplines?: boolean;
    showDroplineLabels?: boolean;
    fitFunc: FitFunction;
    dataBounds: DG.Rect;
    curveFunc: (x: number) => number;
    logOptions: LogOptions;
}

interface FitLabelsRenderOptions {
    names?: string[];
    startLine?: number;
    color: string;
    dataBox: DG.Rect;
    screenBounds?: DG.Rect;
}

interface FitStatisticsRenderOptions {
    statistics?: string[];
    startLine?: number;
    fitFunc: FitFunction;
    logOptions: LogOptions;
    dataBox: DG.Rect;
    screenBounds?: DG.Rect;
    dataPoints?: {x: number[], y: number[]};
    seriesIdx?: number;
}

interface FitTitleRenderOptions {
    showTitle: boolean;
    title?: string;
    dataBox: DG.Rect;
    screenBounds: DG.Rect;
}

interface FitAxesLabelsRenderOptions extends FitTitleRenderOptions {
    showAxesLabels: boolean;
    xAxisName?: string;
    yAxisName?: string;
}

interface FitLegendRenderOptions {
    showLegend?: boolean;
    dataBox: DG.Rect;
    ratio?: number;
}

interface FitLegendColumnlabelSeriesRenderOptions extends FitLegendRenderOptions {
    columnIdx: number;
    drawnCurvesInLegend: number;
    showColumnLabel?: boolean;
    seriesIdx?: number;
    useAuxLegendNames?: boolean;
}


/** Checks if the color is valid */
export function isColorValid(color: string | null | undefined): boolean {
  if (color === undefined || color === null || color === '')
    return false;
  return DG.Color.fromHtml(color) !== undefined;
}

export function getSeriesColor(series: IFitSeries, seriesIdx: number, colorType: SeriesColorType): string {
  const color = DG.Color.toHtml(colorType === 'outlierColor' ? DG.Color.red : DG.Color.getCategoricalColor(seriesIdx));
  return series[colorType] ? DG.Color.fromHtml(series[colorType]) ? series[colorType] : color : color;
}

export function renderConnectDots(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitRenderOptions): void {
  if (series.connectDots ?? false) {
    const viewport = renderOptions.viewport;
    g.strokeStyle = getSeriesColor(series, renderOptions.seriesIdx!, 'pointColor');
    g.lineWidth = 2 * renderOptions.ratio!;
    g.beginPath();
    for (let i = 0; i < series.points.length; i++) {
      const x = series.points[i].x;
      const y = series.points[i].y;
      const screenX = viewport.xToScreen(x);
      const screenY = viewport.yToScreen(y);
      if (i === 0)
        g.moveTo(screenX, screenY);
      else
        g.lineTo(screenX, screenY);
    }
    g.stroke();
  }
}

export function renderPoints(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions) {
  const screenBounds = options.screenBounds!;
  const showPoints = screenBounds.width < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH ||
      screenBounds.height < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT ? '' : series.showPoints ?? 'points';
  if (showPoints) {
    g.strokeStyle = getSeriesColor(series, options.seriesIdx!, ColorType.POINT);
    if (series.connectDots || series.showPoints === 'points')
      drawPoints(g, series, options);
    else if (['candlesticks', 'both'].includes(series.showPoints!))
      drawCandles(g, series, options);
  }
}

/** Performs points drawing */
/** The two segments of an error bar in screen coordinates, running from the marker's edge outwards so
 * the bar stays readable at any marker size. Null when the deviation is smaller than the marker and
 * would be drawn inside it. */
export function stdevWhisker(yCenter: number, yTop: number, yBottom: number, radius: number):
  {top: [number, number], bottom: [number, number]} | null {
  if (yCenter - yTop <= radius)
    return null;
  return {top: [yTop, yCenter - radius], bottom: [yCenter + radius, yBottom]};
}

function drawPoints(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions): void {
  const ratio = options.ratio!;
  const defaultSize = FitConstants.POINT_PX_SIZE * ratio;
  const viewport = options.viewport;
  const connectDots = series.connectDots;
  const pointColor = getSeriesColor(series, options.seriesIdx!, ColorType.POINT);
  const outlierColor = getSeriesColor(series, options.seriesIdx!, ColorType.OUTLIER);

  for (let i = 0; i < series.points.length!; i++) {
    const p = series.points[i];
    if (p.outlier && !series.showOutliers)
      continue;
    const xScreen = viewport.xToScreen(p.x);
    const yScreen = viewport.yToScreen(p.y);
    const color = connectDots ? pointColor :
      p.outlier ? (p.outlierColor ? DG.Color.fromHtml(p.outlierColor) ? p.outlierColor : outlierColor : outlierColor) :
        (p.color ? DG.Color.fromHtml(p.color) ? p.color : pointColor : pointColor);
    const marker = p.marker ? p.marker as DG.MARKER_TYPE : series.markerType as DG.MARKER_TYPE;
    const outlierMarker = p.outlierMarker ? p.outlierMarker as DG.MARKER_TYPE : series.outlierMarkerType as DG.MARKER_TYPE;
    const size = !connectDots ? p.outlier ? FitConstants.OUTLIER_PX_SIZE * ratio : p.size ? p.size : defaultSize : defaultSize;
    const markerToDraw = !connectDots ? p.outlier ? outlierMarker : marker : marker;

    DG.Paint.marker(g, markerToDraw, xScreen, yScreen, color, size);
    if (p.stdev && !p.outlier) {
      const radius = size / 2;
      const whisker = stdevWhisker(yScreen, viewport.yToScreen(p.y + p.stdev),
        viewport.yToScreen(p.y - p.stdev), radius);
      if (whisker) {
        g.strokeStyle = color;
        // set rather than inherited: the fit line and the droplines both leave a width behind
        g.lineWidth = ratio;
        g.beginPath();
        g.moveTo(xScreen, whisker.top[0]);
        g.lineTo(xScreen, whisker.top[1]);
        g.moveTo(xScreen, whisker.bottom[0]);
        g.lineTo(xScreen, whisker.bottom[1]);
        // caps as wide as the marker, so the extent reads at a glance
        g.moveTo(xScreen - radius, whisker.top[0]);
        g.lineTo(xScreen + radius, whisker.top[0]);
        g.moveTo(xScreen - radius, whisker.bottom[1]);
        g.lineTo(xScreen + radius, whisker.bottom[1]);
        g.stroke();
      }
    }
  }
}

/** Performs candles drawing */
function drawCandles(g: CanvasRenderingContext2D, series: IFitSeries, options: FitRenderOptions) : void {
  const ratio = options.ratio!;
  const viewport = options.viewport;
  const pointColor = getSeriesColor(series, options.seriesIdx!, ColorType.POINT);
  for (let i = 0, candleStart = null; i < series.points.length!; i++) {
    const p = series.points[i];
    if (p.outlier)
      continue;
    const nextSame = i + 1 < series.points.length && series.points[i + 1].x === p.x;
    if (!candleStart && nextSame) { candleStart = i; } else if (candleStart !== null && !nextSame) {
      const values: number[] = [];
      for (let j = candleStart, ind = 0; j <= i; j++, ind++)
        values[ind] = series.points[j].y;
      const boxPlotStats = calculateBoxPlotStatistics(values);

      g.beginPath();
      drawCandlestick(g, {viewport: viewport, ratio: ratio, x: p.x, boxPlotStats: boxPlotStats, color: pointColor});
      g.stroke();

      if (series.showPoints === 'both') {
        for (let ind = 0; ind < values.length; ind++) {
          if (values[ind] < boxPlotStats.lowerAdjacentValue || values[ind] > boxPlotStats.upperAdjacentValue) {
            DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER, viewport.xToScreen(p.x),
              viewport.yToScreen(values[ind]), pointColor, FitConstants.CANDLESTICK_OUTLIER_PX_SIZE * ratio);
          }
        }
      }

      candleStart = null;
    }
  }
}

/** Performs candlestick drawing */
function drawCandlestick(g: CanvasRenderingContext2D, renderOptions: FitPointRenderOptions): void {
  const viewport = renderOptions.viewport;
  const x = renderOptions.x;
  const screenX = viewport.xToScreen(x);
  const boxPlotStats = renderOptions.boxPlotStats;

  drawCandlestickBorder(g, x, boxPlotStats.lowerAdjacentValue, viewport);
  g.moveTo(screenX, viewport.yToScreen(boxPlotStats.lowerAdjacentValue));
  g.lineTo(screenX, viewport.yToScreen(boxPlotStats.upperAdjacentValue));
  drawCandlestickBorder(g, x, boxPlotStats.upperAdjacentValue, viewport);
  DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE, screenX, viewport.yToScreen(boxPlotStats.q2), renderOptions.color,
    FitConstants.CANDLESTICK_MEDIAN_PX_SIZE * renderOptions.ratio!);
}

/** Performs candlestick border drawing */
function drawCandlestickBorder(g: CanvasRenderingContext2D, x: number, adjacentValue: number, transform: Viewport): void {
  const xScreen = transform.xToScreen(x);
  const yScreen = transform.yToScreen(adjacentValue);
  g.moveTo(xScreen - (FitConstants.CANDLESTICK_BORDER_PX_SIZE / 2), yScreen);
  g.lineTo(xScreen + (FitConstants.CANDLESTICK_BORDER_PX_SIZE / 2), yScreen);
}

/** Gets rendering variables */
function getRenderingVariables(renderOptions: FitConfidenceIntervalRenderOptions | FitLineRenderOptions) {
  const logX = renderOptions.logOptions.logX;
  const logY = renderOptions.logOptions.logY;
  const screenBounds = renderOptions.screenBounds;
  const axesLeftPxMargin = renderOptions.showAxes ? renderOptions.showAxesLabels ?
    FitConstants.AXES_LEFT_PX_MARGIN_WITH_AXES_LABELS : FitConstants.AXES_LEFT_PX_MARGIN : 0;
  const axesRightPxMargin = renderOptions.showAxes ? FitConstants.AXES_RIGHT_PX_MARGIN : 0;
  const xMin = screenBounds.x + axesLeftPxMargin;
  const xMax = screenBounds.x + screenBounds.width - axesRightPxMargin;
  const viewport = renderOptions.viewport;
  return {logX, logY, xMin, xMax, viewport};
}

export function renderFitLine(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitLineRenderOptions): void {
  if (series.showFitLine ?? true) {
    g.save();
    g.strokeStyle = getSeriesColor(series, renderOptions.seriesIdx!, ColorType.FIT_LINE);
    g.lineWidth = 2 * renderOptions.ratio!;
    g.beginPath();
    if (series.lineStyle)
      g.setLineDash(FitConstants.LINE_STYLES[series.lineStyle]);
    const {logX, logY, xMin, xMax, viewport} = getRenderingVariables(renderOptions);
    const curveFunc = renderOptions.curveFunc!;

    for (let i = xMin; i <= xMax; i++) {
      const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
      const y = logY ? viewport.yToScreen(Math.pow(10, curveFunc(xForY))) : viewport.yToScreen(curveFunc(xForY));
      if (i === xMin)
        g.moveTo(i, y);
      else
        g.lineTo(i, y);
    }
    g.stroke();
    g.restore();
  }
}

export function renderConfidenceIntervals(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitConfidenceIntervalRenderOptions): void {
  if ((series.showFitLine ?? true) && (series.showCurveConfidenceInterval ?? false)) {
    // TODO: improve confidence intervals colors
    g.strokeStyle = series.confidenceIntervalColor ?? FitConstants.CONFIDENCE_INTERVAL_STROKE_COLOR;
    g.fillStyle = series.confidenceIntervalColor ?? FitConstants.CONFIDENCE_INTERVAL_FILL_COLOR;
    const showAxes = renderOptions.showAxes;
    const showAxesLabels = renderOptions.showAxesLabels;
    const logOptions = renderOptions.logOptions;
    const viewport = renderOptions.viewport;
    const screenBounds = renderOptions.screenBounds;
    const dataPoints = renderOptions.dataPoints;

    const confidenceIntervals = getSeriesConfidenceInterval(series, renderOptions.fitFunc!, renderOptions.userParamsFlag!, dataPoints, logOptions);
    drawConfidenceInterval(g, {viewport: viewport, logOptions: logOptions, showAxes: showAxes, showAxesLabels: showAxesLabels,
      screenBounds: screenBounds, confidenceIntervals: confidenceIntervals, confidenceType: FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP});
    drawConfidenceInterval(g, {viewport: viewport, logOptions: logOptions, showAxes: showAxes, showAxesLabels: showAxesLabels,
      screenBounds: screenBounds, confidenceIntervals: confidenceIntervals, confidenceType: FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM});
    fillConfidenceInterval(g, {viewport: viewport, logOptions: logOptions, showAxes: showAxes, showAxesLabels: showAxesLabels,
      screenBounds: screenBounds, confidenceIntervals: confidenceIntervals});
  }
}

/** Performs a curve confidence interval drawing */
function drawConfidenceInterval(g: CanvasRenderingContext2D, renderOptions: FitConfidenceIntervalRenderOptions): void {
  g.beginPath();

  const {logX, logY, xMin, xMax, viewport} = getRenderingVariables(renderOptions);
  const confIntervalFunc = renderOptions.confidenceType === FitConstants.CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
      renderOptions.confidenceIntervals!.confidenceTop : renderOptions.confidenceIntervals!.confidenceBottom;

  for (let i = xMin; i <= xMax; i++) {
    const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
    const y = logY ? viewport.yToScreen(Math.pow(10, confIntervalFunc(xForY))) : viewport.yToScreen(confIntervalFunc(xForY));
    if (i === xMin)
      g.moveTo(i, y);
    else
      g.lineTo(i, y);
  }
  g.stroke();
}

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, renderOptions: FitConfidenceIntervalRenderOptions): void {
  g.beginPath();

  const {logX, logY, xMin, xMax, viewport} = getRenderingVariables(renderOptions);
  const confTop = renderOptions.confidenceIntervals!.confidenceTop;
  const confBottom = renderOptions.confidenceIntervals!.confidenceBottom;

  for (let i = xMin; i <= xMax; i++) {
    const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
    const y = logY ? viewport.yToScreen(Math.pow(10, confTop(xForY))) : viewport.yToScreen(confTop(xForY));
    if (i === xMin)
      g.moveTo(i, y);
    else
      g.lineTo(i, y);
  }

  // reverse traverse to make a shape of confidence interval to fill it
  for (let i = xMax; i >= xMin; i--) {
    const xForY = logX ? Math.log10(viewport.xToWorld(i)) : viewport.xToWorld(i);
    const y = logY ? viewport.yToScreen(Math.pow(10, confBottom(xForY))) : viewport.yToScreen(confBottom(xForY));
    g.lineTo(i, y);
  }
  g.closePath();
  g.fill();
}

/** The travelled fraction an ICxx / ECxx name asks for, or undefined when the name is not one. */
export function droplineFraction(name: string): number | undefined {
  const match = /^(?:IC|EC)(\d{1,2}(?:\.\d+)?)$/i.exec(name);
  const percent = match ? parseFloat(match[1]) : NaN;
  return percent > 0 && percent < 100 ? percent / 100 : undefined;
}

export function renderDroplines(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitDroplineRenderOptions): void {
  if (!(series.showFitLine ?? true) || !series.droplines?.length || !series.parameters || !renderOptions.showDroplines)
    return;
  g.save();
  g.strokeStyle = 'blue';
  g.lineWidth = renderOptions.ratio!;
  g.beginPath();
  g.setLineDash([5, 5]);
  const drawn: {name: string, screenY: number}[] = [];
  for (const name of series.droplines) {
    const fraction = droplineFraction(name);
    const x = fraction === undefined ? undefined :
      renderOptions.fitFunc.inverse(Float32Array.from(series.parameters), fraction);
    const screenY = x === undefined ? null : drawDropline(g, x, renderOptions);
    if (screenY !== null)
      drawn.push({name: name, screenY: screenY});
  }
  g.stroke();
  // a single line is unambiguous; several are told apart by where each one leaves the axis, which
  // is why the caption sits there rather than under the curve
  if (series.droplines.length > 1 && renderOptions.showDroplineLabels) {
    g.setLineDash([]);
    g.fillStyle = 'blue';
    g.font = '10px Roboto, "Roboto Local"';
    g.textAlign = 'left';
    const left = renderOptions.viewport.xToScreen(renderOptions.dataBounds.minX) + 2;
    for (const dropline of drawn)
      g.fillText(dropline.name, left, dropline.screenY - 3);
  }
  g.restore();
}

/** Performs a dropline drawing */
/** Draws one dropline, returning the screen y it leaves the axis at - or null when it falls off the plot. */
function drawDropline(g: CanvasRenderingContext2D, x: number, renderOptions: FitDroplineRenderOptions): number | null {
  const logX = renderOptions.logOptions.logX;
  const logY = renderOptions.logOptions.logY;
  const xValue = logX ? Math.pow(10, x) : x;
  const viewport = renderOptions.viewport;
  const dataBounds = renderOptions.dataBounds;
  if (xValue < dataBounds.minX || xValue > dataBounds.maxX)
    return null;

  const xForY = logX ? Math.log10(xValue) : xValue;
  const y = logY ? Math.pow(10, renderOptions.curveFunc(xForY)) : renderOptions.curveFunc(xForY);
  const screenX = viewport.xToScreen(xValue);
  const screenY = viewport.yToScreen(y);
  g.moveTo(viewport.xToScreen(dataBounds.minX), screenY);
  g.lineTo(screenX, screenY);
  g.lineTo(screenX, viewport.yToScreen(dataBounds.minY));
  return screenY;
}

// formatNumber is fixed at 2 decimals, which renders a sub-micromolar IC50 as "0.00"
function formatStatistic(value: number): string {
  return value !== 0 && (Math.abs(value) < 0.001 || Math.abs(value) >= 1e6) ?
    DG.format(value, 'scientific') : StringUtils.formatNumber(value);
}

/** Draws one series' statistics from `startLine`, returning the lines used so the next can continue
 * below. Lines falling outside the plot are dropped. */
/** Draws the named labels from `labels` starting at `startLine`, returning the lines used so the
 * caller can continue below. Same line budget as the statistics - a cell has room for about five. */
export function renderLabels(g: CanvasRenderingContext2D, labels: {[key: string]: string | number | boolean} | undefined,
  renderOptions: FitLabelsRenderOptions): number {
  const screenBounds = renderOptions.screenBounds!;
  const names = screenBounds.width < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH ||
      screenBounds.height < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT ? [] : renderOptions.names;
  if (!labels || !names || names.length === 0)
    return 0;
  const dataBox = renderOptions.dataBox;
  let line = renderOptions.startLine ?? 0;
  for (const name of names) {
    const value = labels[name];
    if (value === undefined || value === null)
      continue;
    const y = dataBox.y + 20 + 20 * line;
    if (y > dataBox.maxY)
      break;
    g.fillStyle = renderOptions.color;
    g.textAlign = 'left';
    g.fillText(`${name}: ${typeof value === 'number' ? formatStatistic(value) : value}`, dataBox.x + 5, y);
    line++;
  }
  return line - (renderOptions.startLine ?? 0);
}

export function renderStatistics(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitStatisticsRenderOptions): number {
  const screenBounds = renderOptions.screenBounds!;
  const statistics = screenBounds.width < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH ||
      screenBounds.height < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT ? [] : renderOptions.statistics;
  if (!((series.showFitLine ?? true) && statistics && statistics.length > 0))
    return 0;
  const dataBox = renderOptions.dataBox;
  const dataPoints = renderOptions.dataPoints;
  const fitFunc = renderOptions.fitFunc!;
  const fit = toDataSpace(getSeriesFit(series, fitFunc, dataPoints, renderOptions.logOptions),
    renderOptions.logOptions);
  const color = getSeriesColor(series, renderOptions.seriesIdx!, ColorType.FIT_LINE);
  let line = renderOptions.startLine ?? 0;
  for (const statName of statistics) {
    const value = getStatistic(fit, statName);
    const prop = getStatisticProperty(fitFunc, statName);
    // skip statistics this fit function does not produce instead of rendering NaN
    if (value === undefined || !prop)
      continue;
    const y = dataBox.y + 20 + 20 * line;
    if (y > dataBox.maxY)
      break;
    g.fillStyle = color;
    g.textAlign = 'left';
    g.fillText(`${prop.friendlyName}: ${formatStatistic(value)}`, dataBox.x + 5, y);
    line++;
  }
  return line - (renderOptions.startLine ?? 0);
}

export function renderTitle(g: CanvasRenderingContext2D, renderOptions: FitTitleRenderOptions): void {
  if (renderOptions.showTitle) {
    g.font = '12px Roboto, "Roboto Local"';
    g.textAlign = 'center';
    g.fillStyle = 'black';
    g.fillText(renderOptions.title!, renderOptions.dataBox.midX - 5, renderOptions.screenBounds.y + 15);
  }
}

export function renderAxesLabels(g: CanvasRenderingContext2D, renderOptions: FitAxesLabelsRenderOptions): void {
  if (renderOptions.showAxesLabels) {
    const screenBounds = renderOptions.screenBounds;
    const dataBox = renderOptions.dataBox;
    g.font = '11px Roboto, "Roboto Local"';
    g.textAlign = 'center';
    g.fillStyle = 'black';
    g.fillText(renderOptions.xAxisName!, dataBox.midX - 5, screenBounds.maxY - FitConstants.X_AXIS_LABEL_BOTTOM_PX_MARGIN);
    g.translate(screenBounds.x, screenBounds.y);
    g.rotate(-Math.PI / 2);
    const axesTopPxMargin = renderOptions.showTitle ? FitConstants.AXES_TOP_PX_MARGIN_WITH_TITLE : FitConstants.AXES_TOP_PX_MARGIN;
    g.fillText(renderOptions.yAxisName!, -(dataBox.height / 2 + axesTopPxMargin + 15), 15);
    g.restore();
  }
}

export function renderLegend(g: CanvasRenderingContext2D, data: IFitChartData, renderOptions: FitLegendRenderOptions): void {
  if (renderOptions.showLegend) {
    const dataBox = renderOptions.dataBox;
    const ratio = renderOptions.ratio;
    g.font = '11px Roboto, "Roboto Local"';
    // the legend positions text by measured width, so it must own its alignment - the axes labels
    // rendered just before leave textAlign centred
    g.textAlign = 'left';
    const columnNames = [...new Set(data.series?.map((series) => series.columnName))]
      .filter((colName) => colName !== null && colName !== undefined);
    let drawnCurvesInLegend = 0;
    for (let i = 0; i < columnNames.length; i++) {
      const colName = columnNames[i];
      if (data.chartOptions?.showColumnLabel)
        renderColumnLabel(g, colName!, {dataBox: dataBox, columnIdx: i, drawnCurvesInLegend});

      const series = data.series?.filter((series) => series.columnName === colName);
      for (let j = 0; j < series?.length!; j++) {
        const currentSeries = series![j];
        if (currentSeries.name === '' || currentSeries.name === null || currentSeries.name === undefined)
          continue;
        renderLegendSeries(g, currentSeries, {dataBox: dataBox, ratio: ratio, useAuxLegendNames: data.chartOptions?.useAuxLegendNames,
          columnIdx: i, drawnCurvesInLegend, showColumnLabel: data.chartOptions?.showColumnLabel, seriesIdx: j});
        drawnCurvesInLegend++;
      }
    }
    g.restore();
  }
}

function renderColumnLabel(g: CanvasRenderingContext2D, colName: string, renderOptions: FitLegendColumnlabelSeriesRenderOptions): void {
  const dataBox = renderOptions.dataBox;
  g.beginPath();
  g.fillStyle = 'black';
  const colNameWidth = g.measureText(colName!).width;
  g.fillText(colName!, dataBox.maxX - colNameWidth, dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN +
      (renderOptions.columnIdx * FitConstants.LEGEND_RECORD_PX_HEIGHT + FitConstants.LEGEND_RECORD_PX_HEIGHT * renderOptions.drawnCurvesInLegend));
}

function renderLegendSeries(g: CanvasRenderingContext2D, series: IFitSeries, renderOptions: FitLegendColumnlabelSeriesRenderOptions): void {
  const dataBox = renderOptions.dataBox;
  const columnIdx = renderOptions.columnIdx;
  const ratio = renderOptions.ratio!;
  const showColumnLabel = renderOptions.showColumnLabel;
  const drawnCurvesInLegend = renderOptions.drawnCurvesInLegend;
  const pointColor = getSeriesColor(series, renderOptions.seriesIdx!, ColorType.POINT);
  const fitLineColor = getSeriesColor(series, renderOptions.seriesIdx!, ColorType.FIT_LINE);
  g.beginPath();
  g.strokeStyle = fitLineColor;
  g.lineWidth = 2 * ratio;
  const nameToRender = renderOptions.useAuxLegendNames ? series.auxLegendName ?? series.name! : series.name!;
  const textWidth = g.measureText(nameToRender).width;
  g.moveTo(dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_PX_WIDTH - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN,
    dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN +
      (showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) +
      FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend);
  g.lineTo(dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN,
    dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN +
      (showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) +
      FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend);
  let pointsMarkerValue = series.points && series.points.length > 0 ? series.points[0].marker : null;
  pointsMarkerValue = series.points && series.points.every((p) => p.marker === pointsMarkerValue && p.marker != null) ? pointsMarkerValue : null;
  const marker = pointsMarkerValue != null ? pointsMarkerValue as DG.MARKER_TYPE : series.markerType ? series.markerType as DG.MARKER_TYPE : DG.MARKER_TYPE.CIRCLE;
  const legendMarkerSize = Math.floor(FitConstants.POINT_PX_SIZE * ratio * 1.5);
  DG.Paint.marker(g, marker,
    dataBox.maxX - textWidth - FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_PX_WIDTH - legendMarkerSize,
    dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN + (showColumnLabel ? FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) + FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend,
    pointColor, legendMarkerSize);
  g.fillStyle = fitLineColor;
  g.fillText(nameToRender, dataBox.maxX - textWidth, dataBox.y + FitConstants.LEGEND_TOP_PX_MARGIN + (showColumnLabel ?
    FitConstants.LEGEND_RECORD_PX_HEIGHT * (columnIdx + 1) : 0) + FitConstants.LEGEND_RECORD_PX_HEIGHT * drawnCurvesInLegend);
  g.stroke();
}
