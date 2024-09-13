import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ALPHABET, NOTATION} from './consts';
import {SeqHandler} from '../seq-handler';

/** Canonical gap symbol */
export const GAP_SYMBOL: string = '';

export type SeqSplittedBase = ArrayLike<string> & Iterable<string>;

export interface ISeqSplitted {
  get canonicals(): SeqSplittedBase;
  get originals(): SeqSplittedBase;

  isGap(posIdx: number): boolean;

  /** */
  getCanonical(posIdx: number): string;

  /** For fasta and Helm must not be enclosed to square brackets [meA].*/
  getOriginal(posIdx: number): string;

  length: number;
}

export interface INotationProvider {
  get splitter(): SplitterFunc;

  getHelm(seqCol: DG.Column<string>, options?: any): Promise<DG.Column<string>>;
}

export type SeqColStats = { freq: MonomerFreqs, sameLength: boolean }
export type SplitterFunc = (seq: string) => ISeqSplitted;
export type MonomerFreqs = { [m: string]: number };

/** Alphabet candidate type */
export class CandidateType {
  name: string;
  alphabet: Set<string>;
  cutoff: number;

  constructor(name: string, alphabet: Set<string>, cutoff: number) {
    this.name = name;
    this.alphabet = alphabet;
    this.cutoff = cutoff;
  }
}

/** Alphabet candidate similarity type */
export class CandidateSimType extends CandidateType {
  freq: MonomerFreqs;
  /** Cos, max = 1, min = 0 */
  similarity: number;

  constructor(candidate: CandidateType, freq: MonomerFreqs, similarity: number) {
    super(candidate.name, candidate.alphabet, candidate.cutoff);
    this.freq = freq;
    this.similarity = similarity;
  }
}
