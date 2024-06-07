import {Nucleotide} from './types';

export const NUCLEOTIDE_SEMTYPE = 'dna_nucleotide';

export const COMPLEMENT_TABLE: Record<Nucleotide, Nucleotide> = {
  'A': 'T',
  'T': 'A',
  'G': 'C',
  'C': 'G',
};

export const complement = (nucleotides: string): string => (
  Array.from(nucleotides).map((char: string): string => {
    const upperChar = char.toUpperCase();
    const replaceable = upperChar in COMPLEMENT_TABLE;
    if (!replaceable) return char;
    const isUpperCase = upperChar === char;
    const replaced = COMPLEMENT_TABLE[upperChar as Nucleotide];
    return isUpperCase ? replaced : replaced.toLowerCase();
  }).join('')
);

export const generateSubsequences = (sequence: string, n: number): string[] => {
  if (n < 1) return [];

  const subsequences = new Array<string>(sequence.length - n + 1);
  for (let i = 0; i < subsequences.length; ++i)
    subsequences[i] = sequence.substring(i, i + n);

  return subsequences;
};
