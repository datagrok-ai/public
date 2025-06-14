import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {CellRendererBackBase} from '../cell-renderer-back-base';
import {ISeqHandler} from './seq-handler';

export type SeqSplittedBase = ArrayLike<string> & Iterable<string>;

export interface ISeqSplitted {
  isGap(posIdx: number): boolean;

  /** */
  getCanonical(posIdx: number): string;

  /** For fasta and Helm must not be enclosed to square brackets [meA].*/
  getOriginal(posIdx: number): string;

  // TODO: Get ISeqMonomer for seq position
  // get(posIdx: number): ISeqMonomer;

  length: number;

  /** Returns list of canonical monomers in the region specified */
  getCanonicalRegion(start: number, end: number): string[];

  /** Returns the list of original monomers in the region specified */
  getOriginalRegion(start: number, end: number): string[];

  get gapOriginal(): string;
}

export interface INotationProvider {
  get defaultGapOriginal(): string;

  /** Adjust {@link seqHandler} units, {@link seqHandler.column.tags} by {@link seqHandler} constructor */
  setUnits(seqHandler: ISeqHandler): void;

  get splitter(): SplitterFunc;

  /** Any Macromolecule can be presented as Helm notation */
  getHelm(seq: string, options: any): string;

  createCellRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>): CellRendererBackBase<string>;
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
