/* eslint-disable valid-jsdoc */
import {FitParamBounds, IFitSeries, LogOptions} from './fit-curve';

/** Low-level point helpers. A leaf, so the engine and the series API can both use it without a cycle. */

/** Returns the data points of a series with filtered outliers and logarithmic data if needed */
export function getDataPoints(series: IFitSeries, logOptions?: LogOptions, userParamsFlag?: boolean):
  {x: number[], y: number[]} {
  const pointsToMap = userParamsFlag ? series.points : series.points.filter((p) => !p.outlier);
  return {x: pointsToMap.map((p) => logOptions?.logX ? Math.log10(p.x) : p.x),
    y: pointsToMap.map((p) => logOptions?.logY ? Math.log10(p.y) : p.y)};
}

/** Returns median from within multiple points */
export function getMedian(points: {x: number[], y: number[]}): number {
  const mid = Math.floor(points.y.length / 2);
  const sortedPoints = points.y.sort((a, b) => a - b);
  const median = sortedPoints.length % 2 === 0 ? (sortedPoints[mid - 1] + sortedPoints[mid]) / 2 : sortedPoints[mid];
  return median;
}

/** Returns median points from within multiple points with the same x. */
export function getMedianPoints(data: {x: number[], y: number[]}): {x: number[], y: number[]} {
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
export function logIC50ParameterBounds(ic50Bounds: FitParamBounds): FitParamBounds {
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
