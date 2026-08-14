/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {IFitChartData, IFitChartOptions, IFitSeries, LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {Viewport} from '@datagrok-libraries/utils/src/transform';

import {
  getChartBounds,
  getSeriesFitFunction,
  getCurve,
  seriesInFitSpace,
} from '@datagrok-libraries/statistics/src/fit/fit-data';

import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {FitFunction} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {
  renderAxesLabels,
  renderConfidenceIntervals, renderConnectDots,
  renderDroplines,
  renderFitLine,
  renderPoints, renderStatistics, renderTitle, renderLabels, getSeriesColor, ColorType, CURVE_SAMPLE_PX_STEP
} from './render-utils';
import {isLegendVisible, renderLegend} from './fit-legend';
import {getChartDataAggrStats} from './fit-statistics';
import {
  fittedCurves, parsedCurves, curvesDataPoints, getOrCreateParsedChartData, getOrCreateCachedFitCurve,
  getOrCreateCachedCurvesDataPoints, mergeSeries, substituteZeroes
} from './fit-chart-data';
import {
  areAxesLabelsShown, areAxesShown, areDroplineLabelsShown, areDroplinesShown, inflateScreenBounds, isTitleShown,
  layoutChart, shownAxesNames,
} from './fit-layout';
import {handleClick, handleMouseMove, inspectCurve} from './fit-interaction';

function clamp(value: number, min: number, max: number): number {
  return Math.min(Math.max(value, min), max);
}

/** How far past its data a curve may stretch the axis before it is left cut instead. */
const CURVE_HEADROOM = 0.15;

/** The plot has to hold the curves it draws, not only the points they were fitted to. */
function withCurvesInView(bounds: DG.Rect, curves: (((x: number) => number) | null)[],
  options: IFitChartOptions | undefined): DG.Rect {
  const logX = options?.logX ?? false;
  const logY = options?.logY ?? false;
  const SAMPLES = 100;
  let min = bounds.minY;
  let max = bounds.maxY;
  for (const curve of curves) {
    if (!curve)
      continue;
    for (let k = 0; k <= SAMPLES; k++) {
      const t = k / SAMPLES;
      const x = logX ? bounds.minX * Math.pow(bounds.maxX / bounds.minX, t) : bounds.minX + bounds.width * t;
      const raw = curve(logX ? Math.log10(x) : x);
      const y = logY ? Math.pow(10, raw) : raw;
      if (isFinite(y)) {
        min = Math.min(min, y);
        max = Math.max(max, y);
      }
    }
  }
  // a bound the user set is the area they asked for, so a curve outside it stays cut
  if (options?.minY != null)
    min = bounds.minY;
  if (options?.maxY != null)
    max = bounds.maxY;
  if (min >= bounds.minY && max <= bounds.maxY)
    return bounds;
  // worth showing whole while it stays near its data; one that runs away is left cut instead
  if (logY && bounds.minY > 0) {
    const room = (Math.log10(bounds.maxY) - Math.log10(bounds.minY)) * CURVE_HEADROOM;
    min = Math.max(min, Math.pow(10, Math.log10(bounds.minY) - room));
    max = Math.min(max, Math.pow(10, Math.log10(bounds.maxY) + room));
  } else {
    min = Math.max(min, bounds.minY - bounds.height * CURVE_HEADROOM);
    max = Math.min(max, bounds.maxY + bounds.height * CURVE_HEADROOM);
  }
  if (min >= bounds.minY && max <= bounds.maxY)
    return bounds;
  const span = max - min;
  const bottom = min < bounds.minY ? (logY && min > 0 ? min / 1.05 : min - span * 0.05) : bounds.minY;
  const top = max > bounds.maxY ? (logY ? max * 1.05 : max + span * 0.05) : bounds.maxY;
  return new DG.Rect(bounds.x, bottom, bounds.width, top - bottom);
}

function showIncorrectFitCell(g: CanvasRenderingContext2D, screenBounds: DG.Rect): void {
  DG.Paint.marker(g, DG.MARKER_TYPE.OUTLIER, screenBounds.midX, screenBounds.midY, DG.Color.red,
    clamp(Math.min(screenBounds.width, screenBounds.height) * 0.8, 0, 30));
}


@grok.decorators.cellRenderer({
  name: 'Fit',
  cellType: 'fit',
  virtual: true,
})
export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return FitConstants.FIT_CELL_TYPE; }

  get cellType() { return FitConstants.FIT_CELL_TYPE; }

  get defaultHeight() { return FitConstants.CELL_DEFAULT_HEIGHT; }

  get defaultWidth() { return FitConstants.CELL_DEFAULT_WIDTH; }

  getDefaultSize(_gridColumn: DG.GridColumn): {width?: number | null, height?: number | null} {
    return {width: FitConstants.CELL_DEFAULT_WIDTH, height: FitConstants.CELL_DEFAULT_HEIGHT};
  }

  // kept as aliases of the shared caches - these are public surface
  static fittedCurves = fittedCurves;
  static parsedCurves = parsedCurves;
  static curvesDataPoints = curvesDataPoints;

  static inflateScreenBounds(rect: DG.Rect): DG.Rect {
    return inflateScreenBounds(rect);
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    handleClick(gridCell, e);
  }

  onDoubleClick(gridCell: DG.GridCell, _e: MouseEvent): void {
    inspectCurve(gridCell, undefined, true);
    _e.preventDefault();
    _e.stopPropagation();
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    handleMouseMove(gridCell, e, this);
  }

  /** [fitIdentities] key the fits for a caller with no grid cell; only here are the log options final. */
  renderCurves(g: CanvasRenderingContext2D, screenBounds: DG.Rect, data: IFitChartData, gridCell?: DG.GridCell,
    fitIdentities?: (string | undefined)[]): void {
    g.save();
    g.beginPath();
    g.rect(screenBounds.x, screenBounds.y, screenBounds.width, screenBounds.height);
    g.clip();

    const isRenderedOnGrid = gridCell?.grid && gridCell?.grid.dart && g.canvas === gridCell?.grid?.canvas; // only use cache for grid, because there is no guarantee that rowIdx will be correct for other places
    // a merged series is fitted at index 0, the key real series 0 uses
    const useFitCache = isRenderedOnGrid && !data.chartOptions?.mergeSeries;
    const tableCell = gridCell?.cell;
    if (data.chartOptions?.allowXZeroes && data.chartOptions?.logX &&
      data.series?.some((series) => series.points.some((p) => p.x === 0)))
      substituteZeroes(data);
    const axesNames = shownAxesNames(screenBounds, data);
    const [dataBox, xAxisBox, yAxisBox] = layoutChart(screenBounds,
      axesNames.x, axesNames.y, isTitleShown(screenBounds, data));

    let dataBounds = getChartBounds(data);
    if ((dataBounds.x < 0 && data.chartOptions) || (dataBounds.x === 0 && data.chartOptions && !data.chartOptions.allowXZeroes))
      data.chartOptions.logX = false;
    if (dataBounds.y <= 0 && data.chartOptions)
      data.chartOptions.logY = false;
    const chartLogOptions: LogOptions = {logX: data.chartOptions?.logX, logY: data.chartOptions?.logY};
    // fitted before the plot is scaled, so that a curve sitting above every point is still in view
    const fits: ({fitFunc: FitFunction, curve: ((x: number) => number) | null, fitSpaceSeries: IFitSeries,
      userParamsFlag: boolean} | null)[] = [];
    for (let i = 0; i < data.series?.length!; i++) {
      const series = data.series![i];
      if (series.points.some((point) => point.x === undefined || point.y === undefined) || series.points.length <= 1) {
        fits.push(null);
        continue;
      }
      series.points.sort((a, b) => a.x - b.x);
      const fitFunc = getSeriesFitFunction(series);
      let curve: ((x: number) => number) | null = null;
      // everything downstream needs fit-space parameters; the cached series keeps its data-space ones
      let fitSpaceSeries = series;
      let userParamsFlag = true;
      if (!(series.connectDots && !(series.showFitLine ?? true))) {
        if (series.parameters) {
          fitSpaceSeries = seriesInFitSpace(series, chartLogOptions);
          curve = getCurve(fitSpaceSeries, fitFunc);
        } else {
          const fitResult = getOrCreateCachedFitCurve(series, i, fitFunc, chartLogOptions, tableCell, useFitCache,
            fitIdentities?.[i]);
          curve = fitResult.fittedCurve;
          fitSpaceSeries = {...series, parameters: [...fitResult.parameters]};
          userParamsFlag = false;
        }
      }
      fits.push({fitFunc, curve, fitSpaceSeries, userParamsFlag});
    }
    // only a curve that levels off at both ends: an unbounded one is an exponent under a log axis,
    // so fitting the plot around it would bury the data
    dataBounds = withCurvesInView(dataBounds,
      fits.map((fit) => fit?.fitFunc.hasAsymptotes ? fit.curve : null), data.chartOptions);

    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);
    const minSize = Math.min(dataBox.width, dataBox.height);
    // TODO: make thinner
    const ratio = minSize > 100 ? 1 : 0.2 + (minSize / 100) * 0.8;
    // where the curves and points land, so the legend can take a corner they leave free
    const drawnAt = isLegendVisible(data, screenBounds) ? [] as DG.Point[] : undefined;

    g.save();
    g.font = '11px Roboto, "Roboto Local"';
    viewport.drawCoordinateGrid(g, xAxisBox, yAxisBox);
    g.restore();

    // statistics and labels of every series share one column, so each continues below the previous.
    // Plot-level labels come first and once, since they describe the cell rather than any one curve
    let statisticsLine = renderLabels(g, data.chartOptions?.labels, {names: data.chartOptions?.showLabels,
      color: FitConstants.PLOT_LABEL_COLOR, dataBox, screenBounds, startLine: 0, drawnAt});

    // with a single series the aggregate is just that series' own value, so the mode only applies to
    // a cell holding several curves
    const mode = (data.series?.length ?? 0) > 1 ? data.chartOptions?.statisticsMode ?? 'series' : 'series';
    const statistics = data.chartOptions?.showStatistics;
    if (statistics?.length && mode !== 'series') {
      const aggrType = data.chartOptions?.aggrType ?? 'med';
      const stats = getChartDataAggrStats(data, aggrType, tableCell);
      const summarised = Object.fromEntries(statistics.map((name) => [`${aggrType} ${name}`, stats[name]]));
      statisticsLine += renderLabels(g, summarised as {[key: string]: number}, {names: Object.keys(summarised),
        color: FitConstants.PLOT_LABEL_COLOR, dataBox, screenBounds, startLine: statisticsLine, drawnAt});
    }
    for (let i = 0; i < data.series?.length!; i++) {
      const series = data.series![i];
      const fit = fits[i];
      if (!fit)
        continue;
      const {fitFunc, curve, fitSpaceSeries, userParamsFlag} = fit;

      // a bound narrows the plot, and what falls outside it is out of view rather than out of bounds
      g.save();
      g.beginPath();
      g.rect(dataBox.x, dataBox.y, dataBox.width, dataBox.height);
      g.clip();
      renderFitLine(g, series, {viewport, ratio, logOptions: chartLogOptions, showAxes: areAxesShown(screenBounds),
        showAxesLabels: areAxesLabelsShown(screenBounds, data), screenBounds, curveFunc: curve!, seriesIdx: i, drawnAt});
      renderConnectDots(g, series, {viewport, ratio, seriesIdx: i});
      renderPoints(g, series, {viewport, ratio, screenBounds, seriesIdx: i});
      // the points and the path they trace, however they are drawn - markers, candlesticks or
      // connected dots. Sampled like a fit line, so a series is not mistaken for a free corner
      if (drawnAt) {
        let previous: DG.Point | null = null;
        for (const point of series.points) {
          if (point.outlier && !(series.showOutliers ?? true))
            continue;
          const at = new DG.Point(viewport.xToScreen(point.x), viewport.yToScreen(point.y));
          drawnAt.push(at);
          if (previous) {
            const steps = Math.floor(Math.abs(at.x - previous.x) / CURVE_SAMPLE_PX_STEP);
            for (let step = 1; step < steps; step++) {
              drawnAt.push(new DG.Point(previous.x + (at.x - previous.x) * step / steps,
                previous.y + (at.y - previous.y) * step / steps));
            }
          }
          previous = at;
        }
      }
      renderConfidenceIntervals(g, fitSpaceSeries, {viewport, logOptions: chartLogOptions, showAxes: areAxesShown(screenBounds),
        showAxesLabels: areAxesLabelsShown(screenBounds, data), screenBounds, fitFunc, userParamsFlag, drawnAt,
        dataPoints: getOrCreateCachedCurvesDataPoints(series, i, chartLogOptions, userParamsFlag, tableCell, useFitCache)});
      renderDroplines(g, fitSpaceSeries, {viewport, ratio, showDroplines: areDroplinesShown(screenBounds),
        showDroplineLabels: areDroplineLabelsShown(screenBounds), fitFunc, dataBounds, curveFunc: curve!,
        logOptions: chartLogOptions, drawnAt});
      g.restore();
      statisticsLine += renderStatistics(g, fitSpaceSeries, {statistics: mode === 'aggregated' ? [] : statistics, fitFunc,
        logOptions: chartLogOptions, dataBox, screenBounds, seriesIdx: i, startLine: statisticsLine, drawnAt});
      statisticsLine += renderLabels(g, series.labels, {names: data.chartOptions?.showLabels,
        color: getSeriesColor(series, i, ColorType.FIT_LINE), dataBox, screenBounds, startLine: statisticsLine, drawnAt});
    }

    renderTitle(g, {showTitle: isTitleShown(screenBounds, data), title: data.chartOptions?.title, dataBox, screenBounds});
    renderAxesLabels(g, {showTitle: isTitleShown(screenBounds, data), dataBox, screenBounds,
      showAxesLabels: areAxesLabelsShown(screenBounds, data), xAxisName: data.chartOptions?.xAxisName,
      yAxisName: data.chartOptions?.yAxisName});
    if (drawnAt)
      renderLegend(g, data, dataBox, ratio, drawnAt);

    g.restore();
  }

  // TODO: Curves: make less margins in the right, also in the left
  // TODO: Curves: add a warning or error in top right if some mistakes in data
  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle): void {
    if (!gridCell.cell.value)
      return;
    if (w < FitConstants.MIN_CELL_RENDERER_PX_WIDTH || h < FitConstants.MIN_CELL_RENDERER_PX_HEIGHT)
      return;

    const isRenderedOnGrid = gridCell?.grid && gridCell?.grid.dart && g.canvas === gridCell?.grid?.canvas; // only use cache for grid, because there is no guarantee that rowIdx will be correct for other places
    let data = getOrCreateParsedChartData(gridCell.cell, isRenderedOnGrid);
    const screenBounds = FitChartCellRenderer.inflateScreenBounds(new DG.Rect(x, y, w, h));

    for (const [_message, condition] of Object.entries(FitConstants.CONDITION_MAP)) {
      if (condition(data.series)) {
        showIncorrectFitCell(g, screenBounds);
        return;
      }
    }

    if (gridCell.cell.column?.name)
      data.series?.forEach((series) => series.columnName = gridCell.cell.column.name);
    // on a copy: this chart data is cached and shared with the statistics
    if (data.chartOptions?.mergeSeries)
      data = {...data, series: [mergeSeries(data.series!)!]};

    g.clearRect(screenBounds.x, screenBounds.y, screenBounds.width, screenBounds.height);
    this.renderCurves(g, screenBounds, data, gridCell);
  }
}
