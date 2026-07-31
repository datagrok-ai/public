/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {getSeriesFit, getSeriesFitFunction, toFitStatistics, toDataSpace}
  from '@datagrok-libraries/statistics/src/fit/fit-data';
import {statisticsProperties, IFitChartData, IFitSeries, FitStatistics, LEGACY_FIT_STATISTICS, LogOptions}
  from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {Fit, FitFunction, getStatistic, getStatisticProperty}
  from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {getOrCreateCachedFitCurve, getOrCreateCachedCurvesDataPoints} from './fit-renderer';

const AGGREGATION_TYPES: {[key: string]: string} = {
  'count': 'totalCount', 'nulls': 'missingValueCount', 'unique': 'uniqueCount', 'values': 'valueCount',
  'min': 'min', 'max': 'max', 'sum': 'sum', 'avg': 'avg', 'stdev': 'stdev', 'variance': 'variance',
  'skew': 'skew', 'kurt': 'kurt', 'med': 'med', 'q1': 'q1', 'q2': 'q2', 'q3': 'q3',
};
/** Returns the typed fit for a series, with the inflection point reported in data space. */
export function calculateSeriesFit(series: IFitSeries, seriesIdx: number, chartLogOptions: LogOptions,
  tableCell?: DG.Cell, useCache: boolean = true, inDataSpace: boolean = true): Fit {
  const fitFunction = getSeriesFitFunction(series);
  if (series.parameters) {
    if (chartLogOptions.logX) {
      if (series.parameters[2] > 0)
        series.parameters[2] = Math.log10(series.parameters[2]);
    }
  } else {
    const params = getOrCreateCachedFitCurve(series, seriesIdx, fitFunction, chartLogOptions, tableCell, useCache).parameters;
    series.parameters = [...params];
  }

  const fit = getSeriesFit(series, fitFunction,
    getOrCreateCachedCurvesDataPoints(series, seriesIdx, chartLogOptions, false, tableCell), chartLogOptions);
  return inDataSpace ? toDataSpace(fit, chartLogOptions) : fit;
}

/** Returns series statistics in the legacy shape. Prefer {@link calculateSeriesFit}. */
export function calculateSeriesStats(series: IFitSeries, seriesIdx: number, chartLogOptions: LogOptions,
  tableCell?: DG.Cell, useCache: boolean = true): FitStatistics {
  return toFitStatistics(calculateSeriesFit(series, seriesIdx, chartLogOptions, tableCell, useCache));
}
export type AggregatedFitStatistics = FitStatistics & {[name: string]: number | undefined};

/** Statistics descriptors covering every fit function used in the cell, without duplicates. */
export function aggregatedStatisticsProperties(chartData: IFitChartData): DG.Property[] {
  const props: DG.Property[] = [];
  for (const series of chartData.series ?? []) {
    for (const prop of getSeriesFitFunction(series).statisticsProperties) {
      if (!props.some((p) => p.name === prop.name))
        props.push(prop);
    }
  }
  return props.length ? props : statisticsProperties;
}

/** Aggregates across all series every statistic that the cell's fit functions produce - not just the
 * seven legacy ones. Values are aggregated in fit space and converted once at the end, so an averaged
 * IC50 is a geometric mean. Legacy names stay available so recorded transforms keep resolving. */
export function getChartDataAggrStats(chartData: IFitChartData, aggrType: string,
  tableCell?: DG.Cell): AggregatedFitStatistics {
  const chartLogOptions: LogOptions = {logX: chartData.chartOptions?.logX, logY: chartData.chartOptions?.logY};
  const values: Map<string, (number | undefined)[]> = new Map();
  const fitFunctionsUsed: FitFunction[] = [];

  for (let i = 0; i < chartData.series?.length!; i++) {
    const series = chartData.series![i];
    if (series.points.every((p) => p.outlier))
      continue;
    const fitFunction = getSeriesFitFunction(series);
    fitFunctionsUsed.push(fitFunction);
    const fit = calculateSeriesFit(series, i, chartLogOptions, tableCell, true, false);
    for (const prop of fitFunction.statisticsProperties) {
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

  // convert once, after aggregating - the same boundary a single series goes through
  toDataSpace(aggregated as unknown as Fit, chartLogOptions);

  // legacy names resolve to whatever the series' fit functions call them today
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

