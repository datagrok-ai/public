/* eslint-disable valid-jsdoc */
import * as DG from 'datagrok-api/dg';

import {
  FitErrorModel,
  getStatistics,
  getFittedCurve,
  FitStatistics,
  FitConfidenceIntervals,
  FitCurve,
  IFitPoint,
  IFitChartData,
  IFitSeries,
  fitChartDataProperties,
  FitParamBounds,
  IFitChartOptions,
  FitErrorModelType,
} from './fit-curve';
import {fitData, FitFunction, fitSeriesProperties,
  getCurveConfidenceIntervals, getOrCreateFitFunction} from './new-fit-API';

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

/** Returns median from within multiple points */
function getMedian(points: {x: number[], y: number[]}): number {
  const mid = Math.floor(points.y.length / 2);
  const sortedPoints = points.y.sort((a, b) => a - b);
  const median = sortedPoints.length % 2 === 0 ? (sortedPoints[mid - 1] + sortedPoints[mid]) / 2 : sortedPoints[mid];
  return median;
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
    const median = getMedian(currentPoints);
    medianPoints.x[medianPoints.x.length] = currentPoints.x[0];
    medianPoints.y[medianPoints.y.length] = median;
    currentPoints.x = [data.x[i]];
    currentPoints.y = [data.y[i]];
  }
  const median = getMedian(currentPoints);
  medianPoints.x[medianPoints.x.length] = currentPoints.x[0];
  medianPoints.y[medianPoints.y.length] = median;

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
  let boundsAdjusted = false;

  if (chartOptions.minX != null && chartOptions.minX <= bounds.maxX &&
    ((!chartOptions.logX) || (chartOptions.logX && chartOptions.minX > 0))) {
    width += x - chartOptions.minX;
    x = chartOptions.minX;
    boundsAdjusted = true;
  }
  if (chartOptions.maxX != null && chartOptions.maxX >= bounds.minX &&
    ((!chartOptions.logX) || (chartOptions.logX && chartOptions.maxX > 0))) {
    width += chartOptions.maxX - (x + width);
    boundsAdjusted = true;
  }
  if (chartOptions.minY != null && chartOptions.minY <= bounds.maxY &&
    ((!chartOptions.logY) || (chartOptions.logY && chartOptions.minY > 0))) {
    height += y - chartOptions.minY;
    y = chartOptions.minY;
    boundsAdjusted = true;
  }
  if (chartOptions.maxY != null && chartOptions.maxY >= bounds.minY &&
    ((!chartOptions.logY) || (chartOptions.logY && chartOptions.maxY > 0))) {
    height += chartOptions.maxY - (y + height);
    boundsAdjusted = true;
  }

  return boundsAdjusted ? new DG.Rect(x, y, width, height) : padBounds(bounds, chartOptions);
}

function padBounds(bounds: DG.Rect, chartOptions: IFitChartOptions, ratio = 0.05): DG.Rect {
  let left = bounds.left;
  let right = bounds.right;
  let top = bounds.top;
  let bottom = bounds.bottom;

  if (chartOptions.logX) {
    left = left > 0 ? left / (1 + ratio) : left;
    right = right * (1 + ratio);
  }
  else {
    const dx = bounds.width * ratio;
    left -= dx;
    right += dx;
  }
  if (chartOptions.logY) {
    top = top > 0 ? top / (1 + ratio) : top;
    bottom = bottom * (1 + ratio);
  }
  else {
    const dy = bounds.height * ratio;
    top -= dy;
    bottom += dy;
  }

  return new DG.Rect(left, top, right - left, bottom - top);
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
      if (xs.some((x) => x === undefined) || ys.some((y) => y === undefined))
        continue;
      bounds = bounds.union(DG.Rect.fromXYArrays(xs, ys));
    }
    return o ? changeBounds(bounds, o!): padBounds(bounds, o!);
  }
}

/** Returns series fit function */
export function getSeriesFitFunction(series: IFitSeries): FitFunction {
  return getOrCreateFitFunction(series.fitFunction!);
}

/** Returns a curve function, either using the pre-computed parameters or by fitting on-the-fly */
export function getCurve(series: IFitSeries, fitFunc: FitFunction): (x: number) => number {
  const params = new Float32Array(series.parameters?.length!);
  params.set(series.parameters!);
  return getFittedCurve(fitFunc.y, params);
}

/** Returns the data points of a series with filtered outliers and logarithmic data if needed */
export function getDataPoints(series: IFitSeries, logOptions?: LogOptions, userParamsFlag?: boolean):
  {x: number[], y: number[]} {
  const pointsToMap = userParamsFlag ? series.points : series.points.filter((p) => !p.outlier);
  return {x: pointsToMap.map((p) => logOptions?.logX ? Math.log10(p.x) : p.x),
    y: pointsToMap.map((p) => logOptions?.logY ? Math.log10(p.y) : p.y)};
}

/** Fits the series data according to the series fitting settings */
export function fitSeries(series: IFitSeries, fitFunc: FitFunction, dataPoints?: {x: number[], y: number[]},
  logOptions?: LogOptions): FitCurve {
  dataPoints ??= getDataPoints(series, logOptions, false);
  if (series.parameterBounds && logOptions?.logX)
    series.parameterBounds[2] = logIC50ParameterBounds(series.parameterBounds[2]);
  return fitData(getMedianPoints(dataPoints), fitFunc, series.errorModel ?? FitErrorModel.CONSTANT as FitErrorModelType,
    series.parameterBounds);
}

/** Returns series confidence interval functions */
export function getSeriesConfidenceInterval(series: IFitSeries, fitFunc: FitFunction, userParamsFlag: boolean,
  dataPoints?: {x: number[], y: number[]}, logOptions?: LogOptions): FitConfidenceIntervals {
  dataPoints ??= getDataPoints(series, logOptions, userParamsFlag);
  if (!series.parameters) {
    const params = fitSeries(series, fitFunc, dataPoints).parameters;
    series.parameters = [...params];
  }
  const params = new Float32Array(series.parameters?.length!);
  params.set(series.parameters!);
  return getCurveConfidenceIntervals(dataPoints, params, fitFunc.y, 0.05,
    series.errorModel ?? FitErrorModel.CONSTANT as FitErrorModelType);
}

/** Returns series statistics */
export function getSeriesStatistics(series: IFitSeries, fitFunc: FitFunction, dataPoints?: {x: number[], y: number[]},
  logOptions?: LogOptions): FitStatistics {
  dataPoints ??= getDataPoints(series, logOptions, false);
  if (!series.parameters) {
    const params = fitSeries(series, fitFunc, dataPoints).parameters;
    series.parameters = [...params];
  }
  const params = new Float32Array(series.parameters?.length!);
  params.set(series.parameters!);
  return getStatistics(dataPoints, params, fitFunc.y, true);
}
