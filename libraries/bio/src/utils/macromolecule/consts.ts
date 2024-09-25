import {CandidateType} from './types';

/** enum type to simplify setting "user-friendly" notation if necessary */
export enum NOTATION {
  FASTA = 'fasta',
  SEPARATOR = 'separator',
  HELM = 'helm',
  /* Requires notation handler */ CUSTOM = 'custom',
}

export const enum ALIGNMENT {
  SEQ_MSA = 'SEQ.MSA',
  SEQ = 'SEQ',
}

export enum ALPHABET {
  DNA = 'DNA',
  RNA = 'RNA',
  PT = 'PT',
  /** Unknown */
  UN = 'UN',
}

export enum TAGS {
  aligned = 'aligned',
  alphabet = 'alphabet',
  alphabetSize = '.alphabetSize',
  alphabetIsMultichar = '.alphabetIsMultichar',
  separator = 'separator',
  isHelmCompatible = '.isHelmCompatible',
  positionNames = '.positionNames',
  positionLabels = '.positionLabels',
  regions = '.regions',
}

export const positionSeparator: string = ', ';

export const monomerRe: RegExp = /(?:\[([A-Za-z0-9_\-,()]+)\])|([A-Za-z\-])/g;

export const helmRe: RegExp = /(PEPTIDE1|DNA1|RNA1)\{([^}]+)}/g;
export const helmPp1Re: RegExp = /\[([^\[\]]+)]/g;

export const Alphabets = new class {
  fasta = {
    peptide: new Set<string>([
      'G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
      'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T',
    ]),
    dna: new Set<string>(['A', 'C', 'G', 'T']),
    rna: new Set<string>(['A', 'C', 'G', 'U']),
  };
}();

export const candidateAlphabets: CandidateType[] = [
  new CandidateType(ALPHABET.PT, Alphabets.fasta.peptide, 0.50),
  new CandidateType(ALPHABET.DNA, Alphabets.fasta.dna, 0.55),
  new CandidateType(ALPHABET.RNA, Alphabets.fasta.rna, 0.55),
];

/** Canonical gap symbol */
export const GAP_SYMBOL: string = '';

export const GapOriginals: {
  [units: string]: string
} = {
  [NOTATION.FASTA]: '-',
  [NOTATION.SEPARATOR]: '',
  [NOTATION.HELM]: '*',
};
