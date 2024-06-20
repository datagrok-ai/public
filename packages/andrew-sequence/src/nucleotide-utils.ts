export type Nucleotide = 'A'| 'T' | 'C' | 'G' | 'a' | 't' | 'c' | 'g';
export type AminoAcid = Nucleotide | 'R' | 'N' | 'D' | 'r' | 'n' | 'd';

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

export const conutMatchedBySeparatorArrayStrategy = (seqences: string[], n: number, separator: number): number[] => {
  const subsequenceLists = seqences.map((seq: string): string[] => (
    generateSubsequences(seq.replaceAll(/\s+/g, ''), n)
  ));

  const matchesList = new Array<number>(subsequenceLists.length).fill(0);
  for (let i = 0; i < subsequenceLists.length; ++i) {
    const isFirstfHalf = i < separator;
    const jEnd = isFirstfHalf ? subsequenceLists.length : separator;
    for (let j = isFirstfHalf ? separator : 0; j < jEnd; ++j) {
      for (const seqEl of subsequenceLists[i])
        matchesList[i] += subsequenceLists[j].filter((n: string): boolean => n === seqEl).length;
    }
  }

  return matchesList;
};

export const conutMatchedBySeparatorMapStrategy = (seqences: string[], n: number, separator: number): number[] => {
  const subsequenceTables = seqences.map((seq: string): Map<string, number> => (
    countSubsequences(seq.replaceAll(/\s+/g, ''), n)
  ));

  const matchesList = new Array<number>(subsequenceTables.length).fill(0);
  for (let i = 0; i < subsequenceTables.length; ++i) {
    const isFirstfHalf = i < separator;
    const jEnd = isFirstfHalf ? subsequenceTables.length : separator;
    for (let j = isFirstfHalf ? separator : 0; j < jEnd; ++j) {
      for (const [seqEl, seqCount] of subsequenceTables[i])
        matchesList[i] += seqCount * (subsequenceTables[j].get(seqEl) ?? 0);
    }
  }

  return matchesList;
};


