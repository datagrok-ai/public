import {
  mmDistanceFunctions,
  mmDistanceFunctionsNames
} from '@datagrok-libraries/bio/src/distance-functions/macromolecule-distance-functions';

function numericDistance(a: number, b: number): number {
  return Math.abs(a - b);
}

type NumberDistanceFunctions = 'numericDistance';

export type DistanceFunctionNames<T> = T extends string ? mmDistanceFunctionsNames : NumberDistanceFunctions;

export const distanceFunctions: Record<DistanceFunctionNames<any>, (a: any, b: any) => number> = {
  'numericDistance': numericDistance,
  [mmDistanceFunctionsNames.HAMMING]: mmDistanceFunctions[mmDistanceFunctionsNames.HAMMING](),
  [mmDistanceFunctionsNames.LEVENSHTEIN]: mmDistanceFunctions[mmDistanceFunctionsNames.LEVENSHTEIN](),
  [mmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]: mmDistanceFunctions[mmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](),
};
