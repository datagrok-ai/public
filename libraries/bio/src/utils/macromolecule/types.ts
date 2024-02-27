import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ALPHABET, NOTATION} from './consts';
import {UnitsHandler} from '../units-handler';

/** Canonical gap symbol */
export const GAP_SYMBOL: string = '';

export interface ISeqMonomer {
  get canonical(): string;

  /** For fasta and Helm must not be enclosed to square brackets [meA].*/
  get original(): string;
}

export class GapSeqMonomer implements ISeqMonomer {
  get canonical(): string { return GAP_SYMBOL; }

  get original(): string { return this._original; }

  get isGap(): boolean { return true; };

  constructor(
    private readonly _original: string
  ) {}
}

export interface ISeqSplitted extends Iterable<ISeqMonomer> {
  [jPos: number]: ISeqMonomer;

  length: number;
}

// export class SeqSplitted extends Array<ISeqMonomer> implements ISeqSplitted {
//   [jPos: number]: ISeqMonomer;
//
//   constructor(
//     mList: ISeqMonomer[],
//     public readonly notation: NOTATION
//   ) {
//     super(...mList);
//   }
// }

export type SeqColStats = { freq: MonomerFreqs, sameLength: boolean }
export type MonomerFunc = (original: string, jPos: number) => ISeqMonomer;
export type SplitterFunc = (seq: string, getMonomer: MonomerFunc) => ISeqSplitted;
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
