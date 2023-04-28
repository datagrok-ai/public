import {distance} from 'fastest-levenshtein';

function levensteinDistance(a: string, b: string): number {
  return distance(a, b);
};

function numericDistance(a: number, b: number): number {
  return Math.abs(a - b);
}

function hammingDistance(a: string, b: string) {
  // hamming distance should be using with same size strings,
  // but still, lets add a check and if so add the difference to the result
  let diff = 0;
  if (a.length !== b.length)
    diff = Math.abs(a.length - b.length);

  let result = 0;
  for (let i = 0; i < Math.min(a.length, b.length); i++) {
    if (a[i] !== b[i])
      result++;
  }
  return result + diff;
}

type StringDistanceFunctions = 'levenstein' | 'hamming';
type NumberDistanceFunctions = 'numericDistance';

export type DistanceFunctionNames<T> = T extends string ? StringDistanceFunctions : NumberDistanceFunctions;

export const distanceFunctions: Record<DistanceFunctionNames<any>, (a: any, b: any) => number> = {
  levenstein: levensteinDistance,
  hamming: hammingDistance,
  numericDistance: numericDistance,
};
