import {hamming} from './hamming';
import {levenstein} from './levenstein';
import {needlemanWunsch} from './needleman-wunsch';
import {mmDistanceFunctionType} from './types';

/** Enum containing currently supported macromolecule distance functions
 * Hamming distance will be used if the sequences are already aligned
 * Needleman distance will be used for protein sequences with known BLOSUM62 matrix
 * Levenshtein distance will be used for nucleotide sequences as for them substitution matrix is same as identity matrix
 */
export enum MmDistanceFunctionsNames {
  HAMMING = 'Hamming',
  LEVENSHTEIN = 'Levenshtein',
  NEEDLEMANN_WUNSCH = 'Needlemann-Wunsch',
  MONOMER_CHEMICAL_DISTANCE = 'Monomer chemical distance'
};

export const mmDistanceFunctions: Record<MmDistanceFunctionsNames, (value?: any) => mmDistanceFunctionType> = {
  [MmDistanceFunctionsNames.HAMMING]: hamming,
  [MmDistanceFunctionsNames.LEVENSHTEIN]: levenstein,
  [MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]: needlemanWunsch,
  [MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE]: hamming
};
