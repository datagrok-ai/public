import {
  mmDistanceFunctions,
  MmDistanceFunctionsNames,
} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

function numericDistance(a: number, b: number): number {
  return Math.abs(a - b);
}

type NumberDistanceFunctions = 'numericDistance';

export type DistanceFunctionNames<T> = T extends string ? MmDistanceFunctionsNames : NumberDistanceFunctions;

export const distanceFunctions: Record<DistanceFunctionNames<any>, (a: any, b: any) => number> = {
  'numericDistance': numericDistance,
  [MmDistanceFunctionsNames.HAMMING]: mmDistanceFunctions[MmDistanceFunctionsNames.HAMMING](),
  [MmDistanceFunctionsNames.LEVENSHTEIN]: mmDistanceFunctions[MmDistanceFunctionsNames.LEVENSHTEIN](),
  [MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]: mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](),
};
