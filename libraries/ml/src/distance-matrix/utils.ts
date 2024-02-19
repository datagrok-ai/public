import {DistanceAggregationMethod, DistanceAggregationMethods} from './types';

export const isNil = (x: any) => x === null || x === undefined;

export function insertSmaller(distancesAr: number[], indexes: number[], num: number, index: number) {
  if (num > distancesAr[distancesAr.length-1])
    return;

  const newPosition = distancesAr.findIndex((v) => num < v);
  distancesAr.pop();
  distancesAr.splice(newPosition, 0, num);
  indexes.pop();
  indexes.splice(newPosition, 0, index);
}

export function insertLarger(distancesAr: number[], indexes: number[], num: number, index: number) {
  if (num < distancesAr[distancesAr.length-1])
    return;

  const newPosition = distancesAr.findIndex((v) => num > v);
  distancesAr.pop();
  distancesAr.splice(newPosition, 0, num);
  indexes.pop();
  indexes.splice(newPosition, 0, index);
}

export function getAggregationFunction(
  aggregationMethod: DistanceAggregationMethod, weights: number[]
): (values: number[]) => number {
  switch (aggregationMethod) {
    case DistanceAggregationMethods.MANHATTAN:
      return (vs: number[]) => vs.reduce((acc, val, idx) => acc + val * weights[idx], 0);
    default:
      return (vs: number[]) => {
        // euclidean
        const sum = vs.reduce((acc, val, idx) => acc + (val * weights[idx]) ** 2, 0);
        return Math.sqrt(sum);
      };
  }
}
