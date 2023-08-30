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
  FitParamBounds,
  IFitChartOptions,
} from './fit-curve';

export type LogOptions = {
  logX: boolean | undefined,
  logY: boolean | undefined
};


/** Creates new object with the default values specified in {@link properties} */
function createFromProperties(properties: DG.Property[]): any {
  const o: any = {};
  for (const p of properties) {
    if (p.defaultValue !== null)
      o[p.name] = p.defaultValue;
  }
  return o;
}

// TODO: set column with fit readonly value (in detectors) - try to only show chart - remove editable or prevent it??
/** Creates default {@link IFitChartData} object */
export function createDefaultChartData(): IFitChartData {
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

/** Returns median points from within multiple points with the same x. */
function getMedianPoints(data: {x: number[], y: number[]}): {x: number[], y: number[]} {
  const medianPoints: {x: number[], y: number[]} = {x: [], y: []};
  const currentPoints: {x: number[], y: number[]} = {x: [data.x[0]], y: [data.y[0]]};
  for (let i = 1; i < data.x.length; i++) {
    if (data.x[i] === currentPoints.x[0]) {
      currentPoints.x[currentPoints.x.length] = data.x[i];
      currentPoints.y[currentPoints.y.length] = data.y[i];
      continue;
    }
    const mid = Math.floor(currentPoints.y.length / 2);
    const sortedPoints = currentPoints.y.sort((a, b) => a - b);
    const median = sortedPoints.length % 2 === 0 ? (sortedPoints[mid - 1] + sortedPoints[mid]) / 2 : sortedPoints[mid];

    medianPoints.x[medianPoints.x.length] = currentPoints.x[0];
    medianPoints.y[medianPoints.y.length] = median;
    currentPoints.x = [data.x[i]];
    currentPoints.y = [data.y[i]];
  }

  return medianPoints;
}

/** Returns logarithmic IC50 parameter bounds. */
function logIC50ParameterBounds(ic50Bounds: FitParamBounds): FitParamBounds {
  if (ic50Bounds) {
    if (ic50Bounds.max !== undefined)
      ic50Bounds.max = Math.log10(ic50Bounds.max);
    if (ic50Bounds.min !== undefined) {
      ic50Bounds.min = ic50Bounds.min === 0 ?
        -Number.MAX_VALUE : Math.log10(ic50Bounds.min);
    }
  }
  return ic50Bounds;
}

function changeBounds(bounds: DG.Rect, chartOptions: IFitChartOptions): DG.Rect {
  let x = bounds.x;
  let y = bounds.y;
  let width = bounds.width;
  let height = bounds.height;

  if (chartOptions.minX !== undefined && chartOptions.minX > 0) {
    width += x - chartOptions.minX;
    x = chartOptions.minX;
  }
  if (chartOptions.maxX !== undefined && chartOptions.maxX > 0)
    width += chartOptions.maxX - (x + width);
  if (chartOptions.minY !== undefined && chartOptions.minY > 0) {
    height += y - chartOptions.minY;
    y = chartOptions.minY;
  }
  if (chartOptions.maxY !== undefined && chartOptions.maxY > 0)
    height += chartOptions.maxY - (y + height);

  return new DG.Rect(x, y, width, height);
}

/** Returns the bounds of an {@link IFitChartData} object */
export function getChartBounds(chartData: IFitChartData): DG.Rect {
  const o = chartData.chartOptions;
  if (!chartData.series?.length || chartData.series.length === 0)
    return new DG.Rect(0, 0, 1, 1);
  else {
    const {xs, ys} = getPointsArrays(chartData.series[0].points);
    let bounds = DG.Rect.fromXYArrays(xs, ys);
    for (let i = 1; i < chartData.series!.length; i++) {
      const {xs, ys} = getPointsArrays(chartData.series[i].points);
      bounds = bounds.union(DG.Rect.fromXYArrays(xs, ys));
    }
    return o ? changeBounds(bounds, o!): bounds;
  }
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
export function fitSeries(series: IFitSeries, fitFunc: FitFunction, logOptions?: LogOptions): FitCurve {
  const data = {x: series.points.filter((p) => !p.outlier).map((p) => logOptions?.logX ? Math.log10(p.x) : p.x),
    y: series.points.filter((p) => !p.outlier).map((p) => logOptions?.logY ? Math.log10(p.y) : p.y)};
  if (series.parameterBounds && logOptions?.logX)
    series.parameterBounds[2] = logIC50ParameterBounds(series.parameterBounds[2]);
  return fitData(getMedianPoints(data), fitFunc, FitErrorModel.Constant, series.parameterBounds);
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
