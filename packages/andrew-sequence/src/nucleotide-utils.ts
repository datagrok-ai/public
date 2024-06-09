import {AminoAcid, Nucleotide} from './types';

export const NUCLEOTIDE_SEMTYPE = 'dna_nucleotide';
export const ENA_ID_SEMTYPE = 'EnaID';

export const COMPLEMENT_TABLE: Partial<Record<Nucleotide, Nucleotide>> = {
  'A': 'T',
  'T': 'A',
  'G': 'C',
  'C': 'G',
};

export const NUCLEOTIDE_COLORS: Partial<Record<AminoAcid, string>> = {
  'A': 'green',
  'T': 'red',
  'G': 'black',
  'C': 'blue',
  'R': 'blue',
  'N': 'aqua',
  'D': 'red',
};

export const checkIfNucleotide = (value: string): boolean => value.toUpperCase() in COMPLEMENT_TABLE;

export const complement = (nucleotides: string): string => (
  Array.from(nucleotides).map((char: string): string => {
    const upperChar = char.toUpperCase();
    if (!checkIfNucleotide(upperChar)) return char;
    const isUpperCase = upperChar === char;
    const replaced = COMPLEMENT_TABLE[upperChar as Nucleotide] ?? char;
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

export const countSubsequences = (sequence: string, n: number): Map<string, number> => {
  const subsequences = new Map<string, number>();
  if (n < 1) return subsequences;

  const seqLen = sequence.length - n + 1;
  for (let i = 0; i < seqLen; ++i) {
    const subsequence = sequence.substring(i, i + n);
    const prevCount = subsequences.get(subsequence) ?? 0;
    subsequences.set(subsequence, prevCount + 1);
  }

  return subsequences;
};
