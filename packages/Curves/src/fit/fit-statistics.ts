/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {getSeriesFit, getSeriesFitFunction, toDataSpace}
  from '@datagrok-libraries/statistics/src/fit/fit-data';
import {statisticsProperties, IFitChartData, IFitSeries, FitStatistics, LEGACY_FIT_STATISTICS, LogOptions}
  from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {Fit, FitFunction, getStatistic, getStatisticProperty}
  from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {seriesInFitSpace} from '@datagrok-libraries/statistics/src/fit/fit-points';
import {getOrCreateCachedFitCurve, getOrCreateCachedCurvesDataPoints, getOrCreateParsedChartData,
  substituteZeroes, mergeSourceChartOptions} from './fit-chart-data';
import {parseCellValue} from './curve-converter';

const AGGREGATION_TYPES: {[key: string]: string} = {
  'count': 'totalCount', 'nulls': 'missingValueCount', 'unique': 'uniqueCount', 'values': 'valueCount',
  'min': 'min', 'max': 'max', 'sum': 'sum', 'avg': 'avg', 'stdev': 'stdev', 'variance': 'variance',
  'skew': 'skew', 'kurt': 'kurt', 'med': 'med', 'q1': 'q1', 'q2': 'q2', 'q3': 'q3',
};

/** Aggregations that stay on the statistic's own scale. Unlogging a count of 3 would report 1000. */
const DATA_SPACE_AGGREGATIONS = new Set(['min', 'max', 'avg', 'med', 'q1', 'q2', 'q3']);

/** Returns the typed fit for a series, with the inflection point reported in data space. */
export function calculateSeriesFit(series: IFitSeries, seriesIdx: number, chartLogOptions: LogOptions,
  tableCell?: DG.Cell, useCache: boolean = true, inDataSpace: boolean = true): Fit {
  const fitFunction = getSeriesFitFunction(series);
  const fitInput = series.parameters ? seriesInFitSpace(series, chartLogOptions) :
    {...series, parameters: [...getOrCreateCachedFitCurve(series, seriesIdx, fitFunction,
      chartLogOptions, tableCell, useCache).parameters]};

  const fit = getSeriesFit(fitInput, fitFunction,
    getOrCreateCachedCurvesDataPoints(series, seriesIdx, chartLogOptions, false, tableCell), chartLogOptions);
  return inDataSpace ? toDataSpace(fit, chartLogOptions) : fit;
}

export type AggregatedFitStatistics = FitStatistics & {[name: string]: number | undefined};

/** Statistics viable for every fit function in the cell - aggregating one only some series produce
 * would average over a subset while still labelling it "across series". */
export function aggregatedStatisticsProperties(chartData: IFitChartData): DG.Property[] {
  const seriesList = chartData.series ?? [];
  if (!seriesList.length)
    return statisticsProperties;
  let common = [...getSeriesFitFunction(seriesList[0]).statisticsProperties];
  for (let i = 1; i < seriesList.length; i++) {
    const names = new Set(getSeriesFitFunction(seriesList[i]).statisticsProperties.map((p) => p.name));
    common = common.filter((p) => names.has(p.name));
  }
  return common;
}

/** Aggregates in fit space and converts once at the end, so an averaged IC50 is a geometric mean. */
export function getChartDataAggrStats(chartData: IFitChartData, aggrType: string,
  tableCell?: DG.Cell): AggregatedFitStatistics {
  const chartLogOptions: LogOptions = {logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY};
  const values: Map<string, (number | undefined)[]> = new Map();
  const fitFunctionsUsed: FitFunction[] = [];
  const common = new Set(aggregatedStatisticsProperties(chartData).map((p) => p.name));

  for (let i = 0; i < chartData.series?.length!; i++) {
    const series = chartData.series![i];
    if (series.points.every((p) => p.outlier))
      continue;
    const fitFunction = getSeriesFitFunction(series);
    fitFunctionsUsed.push(fitFunction);
    const fit = calculateSeriesFit(series, i, chartLogOptions, tableCell, true, false);
    for (const prop of fitFunction.statisticsProperties) {
      if (!common.has(prop.name))
        continue;
      if (!values.has(prop.name))
        values.set(prop.name, []);
      values.get(prop.name)!.push(getStatistic(fit, prop.name));
    }
  }

  const aggregated: AggregatedFitStatistics = {};
  for (const [name, seriesValues] of values) {
    aggregated[name] = seriesValues.some((v) => v === undefined || v === null || isNaN(v!)) ? undefined :
      DG.Stats.fromValues(seriesValues as number[])[AGGREGATION_TYPES[aggrType] as keyof DG.Stats] as number;
  }

  if (DATA_SPACE_AGGREGATIONS.has(aggrType))
    toDataSpace(aggregated as unknown as Fit, chartLogOptions);

  for (const legacyName of LEGACY_FIT_STATISTICS) {
    if (aggregated[legacyName] !== undefined)
      continue;
    for (const fitFunction of fitFunctionsUsed) {
      const prop = getStatisticProperty(fitFunction, legacyName);
      if (prop && aggregated[prop.name] !== undefined) {
        aggregated[legacyName] = aggregated[prop.name];
        break;
      }
    }
  }
  return aggregated;
}

/** Parsed chart data for one cell of a curve column, with x zeroes substituted when needed.
 * Recalculation hands us a detached column holding only the changed rows, so there is no dataframe
 * cell to key the caches on - fall back to parsing the value directly in that case. */
function chartDataAt(curveColumn: DG.Column, rowIdx: number, table?: DG.DataFrame): {data: IFitChartData, cell?: DG.Cell} | null {
  const value = curveColumn.get(rowIdx);
  if (value === null || value === undefined || value === '')
    return null;
  const cell = curveColumn.dataFrame ? curveColumn.dataFrame.cell(rowIdx, curveColumn.name) : undefined;
  // the detached column keeps its tags, and the dataframe comes in as an argument, so both levels of
  // the cascade still apply on the recalculation path
  const data = cell ? getOrCreateParsedChartData(cell, true) :
    mergeSourceChartOptions(parseCellValue(value, curveColumn), curveColumn, table);
  if (data.chartOptions?.allowXZeroes && data.chartOptions?.logX &&
    data.series?.some((series) => series.points.some((p) => p.x === 0)))
    substituteZeroes(data);
  return {data, cell};
}

/** One statistic of one series, resolved by name so legacy names keep working. */
export function curveStatisticAt(curveColumn: DG.Column, rowIdx: number, propName: string, seriesNumber: number,
  table?: DG.DataFrame): number | null {
  // guard per row: one malformed cell must not abort the whole column
  try {
    const parsed = chartDataAt(curveColumn, rowIdx, table);
    if (!parsed)
      return null;
    const series = parsed.data.series?.[seriesNumber];
    if (!series || series.points.every((p) => p.outlier))
      return null;
    const logOptions: LogOptions = {logX: parsed.data.chartOptions?.logX, logY: parsed.data.chartOptions?.logY};
    return getStatistic(calculateSeriesFit(series, seriesNumber, logOptions, parsed.cell, true), propName) ?? null;
  } catch (_) {
    return null;
  }
}

/** One statistic aggregated across every series of a curve. */
export function curveAggrStatisticAt(curveColumn: DG.Column, rowIdx: number, propName: string, aggrType: string,
  table?: DG.DataFrame): number | null {
  try {
    const parsed = chartDataAt(curveColumn, rowIdx, table);
    if (!parsed)
      return null;
    if (parsed.data.series?.every((series) => series.points.every((p) => p.outlier)))
      return null;
    return getChartDataAggrStats(parsed.data, aggrType, parsed.cell)[propName as keyof FitStatistics] ?? null;
  } catch (_) {
    return null;
  }
}
