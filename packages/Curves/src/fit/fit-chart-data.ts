/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {
  IFitChartData,
  IFitSeries,
  IFitChartOptions,
  FitCurve,
  fitChartDataProperties,
  LogOptions,
} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {createDefaultChartData} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {Fit, FitFunction, fitSeries, fitSeriesProperties} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {getDataPoints} from '@datagrok-libraries/statistics/src/fit/fit-points';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {parseCellValue} from './curve-converter';

/** Chart data assembly and the fit caches. A leaf, so the renderer and the statistics do not import
 * each other through it. */

export const fittedCurves: DG.LruCache<string, FitCurve> = new DG.LruCache<string, FitCurve>(1000);
export const parsedCurves: DG.LruCache<string, IFitChartData> = new DG.LruCache<string, IFitChartData>(1000);
export const curvesDataPoints: DG.LruCache<string, {x: number[], y: number[]}> = new DG.LruCache<string, {x: number[], y: number[]}>(2000);

/** Copies the {@link properties} that {@link target} does not define yet from {@link source}. */
export function mergeProperties(properties: DG.Property[], source: any, target: any): void {
  if (!source || !target)
    return;

  for (const p of properties) {
    if (!(p.name in target) && p.name in source)
      target[p.name] = source[p.name];
  }
}

export const CHART_OPTIONS = 'chartOptions';
export const SERIES_OPTIONS = 'seriesOptions';

/** Assigns over {@link target}, unlike {@link mergeProperties} - this is how an explicitly set level
 * outranks the value a series declares for itself. */
function applyExplicit(source: any, names: string[] | undefined, target: any, claimedByCell?: string[]): void {
  if (!source || !target || !names)
    return;
  for (const name of names) {
    if (name in source && !claimedByCell?.includes(name))
      target[name] = source[name];
  }
}

/** Column and dataframe levels applied over a cell, in that order, skipping anything the cell claims. */
function applySourceExplicit(cellChartData: IFitChartData, columnChartData: IFitChartData,
  dfChartData: IFitChartData, section: 'chartOptions' | 'seriesOptions', target: any): void {
  const claimed = cellChartData.explicit?.[section];
  applyExplicit(dfChartData[section], dfChartData.explicit?.[section], target, claimed);
  applyExplicit(columnChartData[section], columnChartData.explicit?.[section], target, claimed);
}

/** Labels are data rather than a setting, so the levels combine per key instead of the nearest one
 * winning outright - an assay name set on the column and a plate statistic set on the cell both
 * belong on the plot. */
function mergeLabels(target: IFitChartOptions | undefined, ...sources: (IFitChartOptions | undefined)[]): void {
  if (!target)
    return;
  for (const source of sources) {
    if (source?.labels)
      target.labels = {...source.labels, ...target.labels};
  }
}

export function mergeChartOptions(chartOptions: IFitChartOptions[]): IFitChartOptions {
  if (chartOptions.length === 0)
    return {};

  let minX = Number.MAX_VALUE;
  let minY = Number.MAX_VALUE;
  let maxX = Number.MIN_VALUE;
  let maxY = Number.MIN_VALUE;
  let xAxisName: string | undefined;
  let yAxisName: string | undefined;
  let title: string | undefined;
  let logX: boolean = false;
  let logY: boolean = false;
  let allowXZeroes: boolean = false;

  for (const options of chartOptions) {
    if (options.minX !== null && options.minX !== undefined)
      minX = Math.min(minX, options.minX);
    if (options.minY !== null && options.minY !== undefined)
      minY = Math.min(minY, options.minY);
    if (options.maxX !== null && options.maxX !== undefined)
      maxX = Math.max(maxX, options.maxX);
    if (options.maxY !== null && options.maxY !== undefined)
      maxY = Math.max(maxY, options.maxY);
    if (options.title !== null && options.title !== undefined)
      title ??= options.title;
    if (options.xAxisName !== null && options.xAxisName !== undefined)
      xAxisName ??= options.xAxisName;
    if (options.yAxisName !== null && options.yAxisName !== undefined)
      yAxisName ??= options.yAxisName;
    if (options.logX !== null && options.logX !== undefined && options.logX)
      logX = true;
    if (options.logY !== null && options.logY !== undefined && options.logY)
      logY = true;
    if (options.allowXZeroes !== null && options.allowXZeroes !== undefined && options.allowXZeroes)
      allowXZeroes = true;
  }

  return {
    minX: minX === Number.MAX_VALUE ? undefined : minX,
    minY: minY === Number.MAX_VALUE ? undefined : minY,
    maxX: maxX === Number.MIN_VALUE ? undefined : maxX,
    maxY: maxY === Number.MIN_VALUE ? undefined : maxY,
    title: title,
    xAxisName: xAxisName,
    yAxisName: yAxisName,
    logX: logX,
    logY: logY,
    allowXZeroes: allowXZeroes,
  };
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
    showOutliers: series[0].showOutliers,
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

function parsedChartDataKey(column: DG.Column, rowIdx: number): string {
  return `tableId: ${column.dataFrame.id} || tableName: ${column.dataFrame.name} || colName: ${column.name} || colVersion: ${column.version} || rowIdx: ${rowIdx}`;
}

/** Returns either cached or constructed chart data for the specified table cell. */
export function getOrCreateParsedChartData(tableCell: DG.Cell, useCache = true): IFitChartData {
  const column = tableCell?.column;
  return (useCache && column && tableCell) ?
    parsedCurves.getOrCreate(parsedChartDataKey(column, tableCell.rowIndex), () => getChartData(tableCell)) :
    getChartData(tableCell);
}

/** The cascade applied without a cell, for the detached column the recalculation path hands us. */
export function mergeSourceChartOptions(cellChartData: IFitChartData, column: DG.Column,
  df?: DG.DataFrame): IFitChartData {
  const columnChartOptions = getColumnChartOptions(column);
  const dfChartOptions = df ? getDataFrameChartOptions(df) : {};
  cellChartData.series ??= [];
  cellChartData.chartOptions ??= columnChartOptions.chartOptions;
  mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions, cellChartData.chartOptions);
  mergeProperties(fitChartDataProperties, dfChartOptions.chartOptions, cellChartData.chartOptions);
  mergeLabels(cellChartData.chartOptions, columnChartOptions.chartOptions, dfChartOptions.chartOptions);
  applySourceExplicit(cellChartData, columnChartOptions, dfChartOptions, CHART_OPTIONS, cellChartData.chartOptions);
  for (const series of cellChartData.series) {
    mergeProperties(fitSeriesProperties, cellChartData.seriesOptions, series);
    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, series);
    mergeProperties(fitSeriesProperties, dfChartOptions.seriesOptions, series);
    applySourceExplicit(cellChartData, columnChartOptions, dfChartOptions, SERIES_OPTIONS, series);
  }
  return cellChartData;
}

/** Strips the stray '|' some cell values carry. Both parse paths must agree. */
export function sanitizeCellValue(value: string): string {
  return value.includes('|') ? value.replaceAll('|', '') : value;
}

/** Constructs {@link IFitChartData} for a cell, applying the column and dataframe levels. */
function getChartData(tableCell: DG.Cell): IFitChartData {
  const cellValue = tableCell.value as string;
  const column = tableCell.column;
  let cellChartData: IFitChartData = column ? (column.type === DG.TYPE.STRING ?
    parseCellValue(cellValue, column) :
    createDefaultChartData()) : JSON.parse(sanitizeCellValue(cellValue) ?? '{}') ?? {};
  if (column?.type === DG.TYPE.STRING && !cellChartData.series?.length)
    cellChartData = parseCellValue(sanitizeCellValue(cellValue), column);

  const columnChartOptions = tableCell.column ? getColumnChartOptions(tableCell.column) : {};
  const dfChartOptions = tableCell.column ? getDataFrameChartOptions(tableCell.dataFrame) : {};

  cellChartData.series ??= [];
  cellChartData.chartOptions ??= columnChartOptions.chartOptions;

  mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions, cellChartData.chartOptions);
  mergeProperties(fitChartDataProperties, dfChartOptions.chartOptions, cellChartData.chartOptions);
  mergeLabels(cellChartData.chartOptions, columnChartOptions.chartOptions, dfChartOptions.chartOptions);
  applySourceExplicit(cellChartData, columnChartOptions, dfChartOptions, CHART_OPTIONS, cellChartData.chartOptions);
  for (const series of cellChartData.series) {
    mergeProperties(fitSeriesProperties, cellChartData.seriesOptions, series);
    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, series);
    mergeProperties(fitSeriesProperties, dfChartOptions.seriesOptions, series);
    applySourceExplicit(cellChartData, columnChartOptions, dfChartOptions, SERIES_OPTIONS, series);
  }

  return cellChartData;
}

/** Kept apart from [fittedCurves]: a viewer merges a cell's series and applies its own log options. */
export const viewerFits: DG.LruCache<string, FitCurve> = new DG.LruCache<string, FitCurve>(2000);
const chartDataIds: WeakMap<IFitChartData, number> = new WeakMap<IFitChartData, number>();
let nextChartDataId = 0;

/** A parsed cell's identity: re-parsing yields a new object, and so a new id, retiring its fits. */
export function chartDataId(data: IFitChartData): number {
  let id = chartDataIds.get(data);
  if (id === undefined)
    chartDataIds.set(data, id = ++nextChartDataId);
  return id;
}

/** The log options belong here with the cell: the same curve fits differently on a log axis. */
function cellCurveKey(column: DG.Column, tableCell: DG.Cell, idx: number, logOptions?: LogOptions): string {
  return `tableId: ${column.dataFrame.id} || tableName: ${column.dataFrame.name} || colName: ${column.name} || ` +
    `colVersion: ${column.version} || rowIdx: ${tableCell.rowIndex} || idx: ${idx} || ` +
    `logX: ${logOptions?.logX} || logY: ${logOptions?.logY}`;
}

/** Returns existing, or fits curve for the specified grid cell and series. */
export function getOrCreateCachedFitCurve(series: IFitSeries, seriesIdx: number, fitFunc: FitFunction<Fit>,
  chartLogOptions: LogOptions, tableCell?: DG.Cell, useCache = true, identity?: string): FitCurve {
  const dataPoints = getOrCreateCachedCurvesDataPoints(series, seriesIdx, chartLogOptions, false, tableCell, useCache);
  const column = tableCell?.column;
  if (identity)
    return viewerFits.getOrCreate(`${identity} || logX: ${chartLogOptions?.logX} || logY: ${chartLogOptions?.logY}`,
      () => fitSeries(series, fitFunc, dataPoints, chartLogOptions));
  return (useCache && column && tableCell) ?
    fittedCurves.getOrCreate(cellCurveKey(column, tableCell, seriesIdx, chartLogOptions), () => {
      return fitSeries(series, fitFunc, dataPoints, chartLogOptions);
    }) : fitSeries(series, fitFunc, dataPoints, chartLogOptions);
}

/** Returns existing, or maps new data points for the specified series. */
export function getOrCreateCachedCurvesDataPoints(series: IFitSeries, idx: number, logOptions?: LogOptions,
  userParamsFlag?: boolean, tableCell?: DG.Cell, useCache = true): {x: number[], y: number[]} {
  const column = tableCell?.column;
  return (useCache && column && tableCell) ?
    curvesDataPoints.getOrCreate(`${cellCurveKey(column, tableCell, idx, logOptions)} || userParamsFlag: ${userParamsFlag}`, () => {
      return getDataPoints(series, logOptions, userParamsFlag);
    }) : getDataPoints(series, logOptions, userParamsFlag);
}

/** Reads a level's options, migrating a value stored under the pre-`.%` tag name so that only one of
 * the two ever holds the options. */
function readChartOptions(tags: any): IFitChartData {
  const stored = tags[FitConstants.TAG_FIT];
  if (stored !== null && stored !== undefined)
    return JSON.parse(stored);
  const legacy = tags[FitConstants.TAG_FIT_LEGACY];
  const migrated = legacy !== null && legacy !== undefined;
  tags[FitConstants.TAG_FIT] = migrated ? legacy : JSON.stringify(createDefaultChartData());
  if (migrated)
    delete tags[FitConstants.TAG_FIT_LEGACY];
  return JSON.parse(tags[FitConstants.TAG_FIT]);
}

/** Returns existing, or creates new dataframe default chart options. */
export function getDataFrameChartOptions(df: DG.DataFrame): IFitChartData {
  return readChartOptions(df.tags);
}

/** Returns existing, or creates new column default chart options. */
export function getColumnChartOptions(column: DG.Column): IFitChartData {
  return readChartOptions(column.tags);
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
