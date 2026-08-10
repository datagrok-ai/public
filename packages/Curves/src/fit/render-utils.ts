/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {
  FitConfidenceIntervals,
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


/** How often what is drawn reports where it went, which is enough to tell a busy corner from a free one. */
export const CURVE_SAMPLE_PX_STEP = 8;

/** Where a piece of text ends up, so the legend does not take a corner the plot already writes in. */
function reportText(drawnAt: DG.Point[] | undefined, g: CanvasRenderingContext2D, text: string,
  x: number, y: number): void {
  if (!drawnAt)
    return;
  const width = g.measureText(text).width;
  for (let at = x; at <= x + width; at += CURVE_SAMPLE_PX_STEP)
    drawnAt.push(new DG.Point(at, y));
}

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
    /** Collects where the chart put ink, which is how the legend finds a free corner. */
    drawnAt?: DG.Point[];
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
    drawnAt?: DG.Point[];
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
    drawnAt?: DG.Point[];
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
    const outlierSize = Math.min(FitConstants.OUTLIER_PX_SIZE * ratio, viewport.screen.height / 10);
    const size = !connectDots ? p.outlier ? outlierSize : p.size ? p.size : defaultSize : defaultSize;
    const markerToDraw = !connectDots ? p.outlier ? outlierMarker : marker : marker;

    DG.Paint.marker(g, markerToDraw, xScreen, yScreen, color, size);
    if (p.stdev && !p.outlier) {
      const radius = size / 2;
      const whisker = stdevWhisker(yScreen, viewport.yToScreen(p.y + p.stdev),
        viewport.yToScreen(p.y - p.stdev), radius);
      if (whisker) {
        options.drawnAt?.push(new DG.Point(xScreen, whisker.top[0]), new DG.Point(xScreen, whisker.bottom[1]));
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
  const viewport = renderOptions.viewport;
  // the box the viewport maps onto, so a curve stops where the axes do instead of running under
  // whatever the layout keeps to the side of it
  const xMin = viewport.screen.left;
  const xMax = viewport.screen.maxX;
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
      if (renderOptions.drawnAt && i % CURVE_SAMPLE_PX_STEP === 0)
        renderOptions.drawnAt.push(new DG.Point(i, y));

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
    if (renderOptions.drawnAt && i % CURVE_SAMPLE_PX_STEP === 0)
      renderOptions.drawnAt.push(new DG.Point(i, y));
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
  const left = viewport.xToScreen(dataBounds.minX);
  const bottom = viewport.yToScreen(dataBounds.minY);
  g.moveTo(left, screenY);
  g.lineTo(screenX, screenY);
  g.lineTo(screenX, bottom);
  if (renderOptions.drawnAt) {
    for (let x = left; x <= screenX; x += CURVE_SAMPLE_PX_STEP)
      renderOptions.drawnAt.push(new DG.Point(x, screenY));
    for (let y2 = screenY; y2 <= bottom; y2 += CURVE_SAMPLE_PX_STEP)
      renderOptions.drawnAt.push(new DG.Point(screenX, y2));
  }
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
    const text = `${name}: ${typeof value === 'number' ? formatStatistic(value) : value}`;
    g.fillText(text, dataBox.x + 5, y);
    reportText(renderOptions.drawnAt, g, text, dataBox.x + 5, y);
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
    const text = `${prop.friendlyName}: ${formatStatistic(value)}`;
    g.fillText(text, dataBox.x + 5, y);
    reportText(renderOptions.drawnAt, g, text, dataBox.x + 5, y);
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
