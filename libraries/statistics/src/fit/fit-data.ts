/* eslint-disable valid-jsdoc */
import * as DG from 'datagrok-api/dg';

import {
  FitErrorModel,
  fitData,
  getCurveConfidenceIntervals,
  getStatistics,
  FitFunction,
  getFittedCurve,
  FitStatistics,
  FitConfidenceIntervals,
  FitCurve,
  getOrCreateFitFunction,
  IFitPoint,
  IFitChartData,
  IFitSeries,
  fitSeriesProperties,
  fitChartDataProperties,
  TAG_FIT,
} from './fit-curve';


/** Creates new object with the default values specified in {@link properties} */
function createFromProperties(properties: DG.Property[]): any {
  const o: any = {};
  for (const p of properties) {
    if (p.defaultValue !== null)
      o[p.name] = p.defaultValue;
  }
  return o;
}

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

/** Creates default {@link IFitChartData} object */
function createDefaultChartData(): IFitChartData {
  return {
    chartOptions: createFromProperties(fitChartDataProperties),
    seriesOptions: createFromProperties(fitSeriesProperties),
  };
}

/** Returns existing, or creates new column default chart options. */
export function getColumnChartOptions(gridColumn: DG.GridColumn): IFitChartData {
  return gridColumn.temp[TAG_FIT] ??= createDefaultChartData();
}

/** Returns points arrays from {@link IFitPoint} array */
export function getPointsArrays(points: IFitPoint[]): {xs: number[], ys: number[]} {
  const xs: number[] = [];
  const ys: number[] = [];
  for (let i = 0; i < points.length; i++) {
    xs[i] = points[i].x;
    ys[i] = points[i].y;
  }
  return {xs: xs, ys: ys};
}

/** Returns the bounds of an {@link IFitChartData} object */
export function getChartBounds(chartData: IFitChartData): DG.Rect {
  const o = chartData.chartOptions;
  if (o?.minX && o.minY && o.maxX && o.maxY)
    return new DG.Rect(o.minX, o.minY, o.maxX - o.minX, o.maxY - o.minY);
  if (!chartData.series?.length || chartData.series.length === 0)
    return new DG.Rect(0, 0, 1, 1);
  else {
    const {xs, ys} = getPointsArrays(chartData.series[0].points);
    let bounds = DG.Rect.fromXYArrays(xs, ys);
    for (let i = 1; i < chartData.series!.length; i++) {
      const {xs, ys} = getPointsArrays(chartData.series[i].points);
      bounds = bounds.union(DG.Rect.fromXYArrays(xs, ys));
    }
    return bounds;
  }
}

/** Constructs {@link IFitChartData} from the grid cell, taking into account
 * chart and fit settings potentially defined on the dataframe and column level. */
export function getChartData(gridCell: DG.GridCell): IFitChartData {
  const cellChartData: IFitChartData = gridCell.cell?.column?.type === DG.TYPE.STRING ?
    (JSON.parse(gridCell.cell.value ?? '{}') ?? {}) : createDefaultChartData();

  const columnChartOptions = getColumnChartOptions(gridCell.gridColumn);

  cellChartData.series ??= [];
  cellChartData.chartOptions ??= columnChartOptions.chartOptions;

  // merge cell options with column options
  mergeProperties(fitChartDataProperties, columnChartOptions.chartOptions, cellChartData.chartOptions);
  for (const series of cellChartData.series)
    mergeProperties(fitSeriesProperties, columnChartOptions.seriesOptions, series);

  return cellChartData;
}

/** Returns series fit function */
export function getSeriesFitFunction(series: IFitSeries): FitFunction {
  return getOrCreateFitFunction(series.fitFunction!);
}

/** Returns a curve function, either using the pre-computed parameters or by fitting on-the-fly */
export function getCurve(series: IFitSeries, fitFunc: FitFunction): (x: number) => number {
  return getFittedCurve(fitFunc.y, series.parameters!);
}

/** Fits the series data according to the series fitting settings */
export function fitSeries(series: IFitSeries, fitFunc: FitFunction): FitCurve {
  const data = {x: series.points.filter((p) => !p.outlier).map((p) => p.x),
    y: series.points.filter((p) => !p.outlier).map((p) => p.y)};
  return fitData(data, fitFunc, FitErrorModel.Constant);
}

/** Returns series confidence interval functions */
export function getSeriesConfidenceInterval(series: IFitSeries, fitFunc: FitFunction,
  userParamsFlag: boolean): FitConfidenceIntervals {
  const data = userParamsFlag ? {x: series.points.map((p) => p.x), y: series.points.map((p) => p.y)} :
    {x: series.points.filter((p) => !p.outlier).map((p) => p.x),
      y: series.points.filter((p) => !p.outlier).map((p) => p.y)};
  if (!series.parameters)
    series.parameters = fitSeries(series, fitFunc).parameters;
  return getCurveConfidenceIntervals(data, series.parameters, fitFunc.y, 0.05, FitErrorModel.Constant);
}

/** Returns series statistics */
export function getSeriesStatistics(series: IFitSeries, fitFunc: FitFunction): FitStatistics {
  const data = {x: series.points.filter((p) => !p.outlier).map((p) => p.x),
    y: series.points.filter((p) => !p.outlier).map((p) => p.y)};
  if (!series.parameters)
    series.parameters = fitSeries(series, fitFunc).parameters;
  return getStatistics(data, series.parameters, fitFunc.y, true);
}
