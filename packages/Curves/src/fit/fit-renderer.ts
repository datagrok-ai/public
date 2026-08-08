/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {IFitChartData, LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {Viewport} from '@datagrok-libraries/utils/src/transform';

import {
  getChartBounds,
  getSeriesFitFunction,
  getCurve,
  seriesInFitSpace,
} from '@datagrok-libraries/statistics/src/fit/fit-data';

import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {
  renderAxesLabels,
  renderConfidenceIntervals, renderConnectDots,
  renderDroplines,
  renderFitLine, renderLegend,
  renderPoints, renderStatistics, renderTitle, renderLabels, getSeriesColor, ColorType
} from './render-utils';
import {getChartDataAggrStats} from './fit-statistics';
import {
  fittedCurves, parsedCurves, curvesDataPoints, getOrCreateParsedChartData, getOrCreateCachedFitCurve,
  getOrCreateCachedCurvesDataPoints, mergeSeries, substituteZeroes
} from './fit-chart-data';
import {
  areAxesLabelsShown, areAxesShown, areDroplinesShown, inflateScreenBounds, isLegendShown, isTitleShown, layoutChart,
} from './fit-layout';
import {handleClick, handleMouseMove, inspectCurve} from './fit-interaction';

function clamp(value: number, min: number, max: number): number {
  return Math.min(Math.max(value, min), max);
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

  renderCurves(g: CanvasRenderingContext2D, screenBounds: DG.Rect, data: IFitChartData, gridCell?: DG.GridCell): void {
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
    const [dataBox, xAxisBox, yAxisBox] = layoutChart(screenBounds,
      areAxesLabelsShown(screenBounds, data), isTitleShown(screenBounds, data));

    const dataBounds = getChartBounds(data);
    if ((dataBounds.x < 0 && data.chartOptions) || (dataBounds.x === 0 && data.chartOptions && !data.chartOptions.allowXZeroes))
      data.chartOptions.logX = false;
    if (dataBounds.y <= 0 && data.chartOptions)
      data.chartOptions.logY = false;
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);
    const minSize = Math.min(dataBox.width, dataBox.height);
    // TODO: make thinner
    const ratio = minSize > 100 ? 1 : 0.2 + (minSize / 100) * 0.8;
    const chartLogOptions: LogOptions = {logX: data.chartOptions?.logX, logY: data.chartOptions?.logY};

    g.save();
    g.font = '11px Roboto, "Roboto Local"';
    viewport.drawCoordinateGrid(g, xAxisBox, yAxisBox);
    g.restore();

    // statistics and labels of every series share one column, so each continues below the previous.
    // Plot-level labels come first and once, since they describe the cell rather than any one curve
    let statisticsLine = renderLabels(g, data.chartOptions?.labels, {names: data.chartOptions?.showLabels,
      color: FitConstants.PLOT_LABEL_COLOR, dataBox, screenBounds, startLine: 0});

    // with a single series the aggregate is just that series' own value, so the mode only applies to
    // a cell holding several curves
    const mode = (data.series?.length ?? 0) > 1 ? data.chartOptions?.statisticsMode ?? 'series' : 'series';
    const statistics = data.chartOptions?.showStatistics;
    if (statistics?.length && mode !== 'series') {
      const aggrType = data.chartOptions?.aggrType ?? 'med';
      const stats = getChartDataAggrStats(data, aggrType, tableCell);
      const summarised = Object.fromEntries(statistics.map((name) => [`${aggrType} ${name}`, stats[name]]));
      statisticsLine += renderLabels(g, summarised as {[key: string]: number}, {names: Object.keys(summarised),
        color: FitConstants.PLOT_LABEL_COLOR, dataBox, screenBounds, startLine: statisticsLine});
    }
    for (let i = 0; i < data.series?.length!; i++) {
      const series = data.series![i];
      if (series.points.some((point) => point.x === undefined || point.y === undefined) || series.points.length <= 1)
        continue;
      series.points.sort((a, b) => a.x - b.x);

      let userParamsFlag = true;
      const fitFunc = getSeriesFitFunction(series);
      let curve: ((x: number) => number) | null = null;
      // everything downstream needs fit-space parameters; the cached series keeps its data-space ones
      let fitSpaceSeries = series;
      if (!(series.connectDots && !(series.showFitLine ?? true))) {
        if (series.parameters) {
          fitSpaceSeries = seriesInFitSpace(series, chartLogOptions);
          curve = getCurve(fitSpaceSeries, fitFunc);
        } else {
          const fitResult = getOrCreateCachedFitCurve(series, i, fitFunc, chartLogOptions, tableCell, useFitCache);
          curve = fitResult.fittedCurve;
          fitSpaceSeries = {...series, parameters: [...fitResult.parameters]};
          userParamsFlag = false;
        }
      }

      renderFitLine(g, series, {viewport, ratio, logOptions: chartLogOptions, showAxes: areAxesShown(screenBounds),
        showAxesLabels: areAxesLabelsShown(screenBounds, data), screenBounds, curveFunc: curve!, seriesIdx: i});
      renderConnectDots(g, series, {viewport, ratio, seriesIdx: i});
      renderPoints(g, series, {viewport, ratio, screenBounds, seriesIdx: i});
      renderConfidenceIntervals(g, fitSpaceSeries, {viewport, logOptions: chartLogOptions, showAxes: areAxesShown(screenBounds),
        showAxesLabels: areAxesLabelsShown(screenBounds, data), screenBounds, fitFunc, userParamsFlag,
        dataPoints: getOrCreateCachedCurvesDataPoints(series, i, chartLogOptions, userParamsFlag, tableCell, useFitCache)});
      // the inflection point, by name - it is slot 2 on the dose-response families and absent on the
      // two-parameter ones, where reading a slot gave NaN
      const inflectionSlot = fitFunc.statisticFields.findIndex((f) => f === 'ic50' || f === 'ec50');
      if (fitSpaceSeries.parameters && inflectionSlot >= 0) {
        renderDroplines(g, fitSpaceSeries, {viewport, ratio, showDroplines: areDroplinesShown(screenBounds),
          xValue: fitSpaceSeries.parameters![inflectionSlot], dataBounds, curveFunc: curve!, logOptions: chartLogOptions});
      }
      statisticsLine += renderStatistics(g, fitSpaceSeries, {statistics: mode === 'aggregated' ? [] : statistics, fitFunc,
        logOptions: chartLogOptions, dataBox, screenBounds, seriesIdx: i, startLine: statisticsLine});
      statisticsLine += renderLabels(g, series.labels, {names: data.chartOptions?.showLabels,
        color: getSeriesColor(series, i, ColorType.FIT_LINE), dataBox, screenBounds, startLine: statisticsLine});
    }

    renderTitle(g, {showTitle: isTitleShown(screenBounds, data), title: data.chartOptions?.title, dataBox, screenBounds});
    renderAxesLabels(g, {showTitle: isTitleShown(screenBounds, data), dataBox, screenBounds,
      showAxesLabels: areAxesLabelsShown(screenBounds, data), xAxisName: data.chartOptions?.xAxisName,
      yAxisName: data.chartOptions?.yAxisName});
    renderLegend(g, data, {showLegend: isLegendShown(screenBounds), dataBox, ratio});

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
