import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {CellRendererBackBase} from '../cell-renderer-back-base';
import {ISeqHandler} from './seq-handler';
import {PolymerType} from '../../helm/types';
import {NOTATION_PROVIDER_CONSTRUCTOR_ROLE} from './consts';

export type SeqSplittedBase = ArrayLike<string> & Iterable<string>;

export interface ISeqConnection {
  seq1Type: PolymerType;
  seq2Type: PolymerType;
  seqIndex1: number;
  seqIndex2: number;
  monomerIndex1: number;
  monomerIndex2: number;
  rGroup1: number;
  rGroup2: number;
}

export interface ISeqGraphInfo {
  /**  Connections between monomers with numbers - note that monomer numbers are 0-based*/
  connections: ISeqConnection[];
  /** Start positions of disjoint sequence parts in the original sequence */
  disjointSeqStarts: number[];
  /** Polimer types same as the indexes. only applicapble to helm */
  polymerTypes?: PolymerType[];
}

export interface ISeqSplitted {

  graphInfo?: ISeqGraphInfo;

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

export abstract class NotationProviderBase {
  /** Name of the custom notation */
  static get notationName(): string {
    return 'Custom';
  };

  /** flag to let bio know if this provider implements method for converting helm to it */
  static get implementsFromHelm(): boolean {
    return false;
  };

  /** Method for converting HELM to this notation */
  static convertFromHelm(helm: string, options: any): string {
    throw new Error(`Method convertFromHelm not implemented for this notation provider`);
  };

  static async getProviderConstructors(): Promise<typeof NotationProviderBase[]> {
    const constFuncs = DG.Func.find({meta: {role: NOTATION_PROVIDER_CONSTRUCTOR_ROLE}});
    return Promise.all(constFuncs.map((f) => f.apply({})));
  }
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
