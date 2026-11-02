import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {CellRendererBackBase} from '../cell-renderer-back-base';
import {ISeqHandler} from './seq-handler';
import {HelmType, PolymerType} from '../../helm/types';

/** meta.role of a custom notation provider's constructor function (re-exported by consts.ts).
 * Defined here so that types.ts has no run-time import of consts.ts: consts.ts builds `candidateAlphabets` from
 * CandidateType when it loads, so the cycle threw "Cannot access 'CandidateType' before initialization" in any bundle
 * that reached it through types.ts or utils.ts first (Helm's package-test.js did). */
export const NOTATION_PROVIDER_CONSTRUCTOR_ROLE = 'notationProviderConstructor';

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

/** A HELM section-1 connection with its kind: hydrogen bonds (`pair`) are 'hbond' and carry rGroup 0 */
export interface ISeqConnectionExt extends ISeqConnection {
  kind?: 'covalent' | 'hbond';
}

/** A fixed repeat `A'3'`, `(A.G)'2'`, `(R(U)P)'3'` over the flat monomer range [start, end), listed once in the
 * splitted; the assembly input copies the range `min` (=== `max`) times in place */
export interface ISeqRepeat {
  start: number;
  end: number;
  min: number;
  max: number;
}

/** Graph info of a HELM sequence split by the hwe parser; extra fields are optional, old consumers ignore them.
 * Mirrors IAsmGraphInfoExt of monomer-works/assembly-input. */
export interface ISeqGraphInfoHelm extends ISeqGraphInfo {
  /** Covalent connections only (kind 'covalent'): hydrogen bonds are in {@link hbonds}. Ambiguous connections
   * (multi-site `(1,3):R3`, ranges `2-5`, `?`, a symbol `C:R3` found more than once, group endpoints) are left out
   * and listed in {@link errors}: hwe draws them at the first candidate, but an old consumer (joinToHelm, joinToBiln,
   * the atomic-level converter) would take that for a definite bond. */
  connections: ISeqConnectionExt[];
  /** Hydrogen bonds (`pair`), kind 'hbond'; kept out of {@link connections} so old consumers never bond them */
  hbonds?: ISeqConnectionExt[];
  /** Per flat monomer: HELM_AA, HELM_SUGAR, HELM_BASE, HELM_LINKER, HELM_CHEM, HELM_BLOB; unset for RNA rows with a
   * lone sugar-slot monomer (e.g. a 5' phosphate), whose role hwe cannot tell */
  biotypes?: HelmType[];
  /** Polymer ids per chain in source order ('PEPTIDE1', 'CHEM1', ...) */
  chainIds?: string[];
  /** Fixed repeats in source order, non-overlapping (hwe's expandRepeats decides which are fixed); ranges `'2-5'`,
   * CHEM and nested repeats are in {@link errors} instead */
  repeats?: ISeqRepeat[];
  /** Problems that make the row unbuildable at atomic level (ambiguous monomers or connections, repeat ranges,
   * dropped connections, ...) */
  errors?: string[];
  /** hwe recovered from a syntax error (listed in {@link errors}): its parse may miss monomers after it, or have
   * none */
  recovered?: boolean;
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

/**
 * Extracts displayable monomer symbol(s) from a raw position value.
 * Used by the Monomer cell renderer when a notation provider needs custom
 * interpretation of split-to-monomer position values.
 *
 * For example, a Parabilis position `[VHL1379]*dLys` contains two monomers
 * (`VHL1379` and `dLys`). The canonicalizer separates them with
 * {@link MONOMER_MOTIF_SPLITTER} so the renderer can display/tooltip each one.
 */
export interface IMonomerCanonicalizer {
  /** Convert a raw position value to a displayable monomer string.
   *  If the position contains multiple monomers, join them with
   *  {@link MONOMER_MOTIF_SPLITTER} (` , `). */
  canonicalize(original: string): string;

  /** Returns true if the raw position value represents a gap in this notation. */
  isGap(original: string): boolean;
}

export interface INotationProvider {
  get defaultGapOriginal(): string;

  /** Adjust {@link seqHandler} units, {@link seqHandler.column.tags} by {@link seqHandler} constructor */
  setUnits(seqHandler: ISeqHandler): void;

  get splitter(): SplitterFunc;

  /** Any Macromolecule can be presented as Helm notation */
  getHelm(seq: string, options: any): string;

  createCellRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>): CellRendererBackBase<string>;

  /** nqName of a registered function that returns an {@link IMonomerCanonicalizer} instance.
   *  When set, split-to-monomers will tag Monomer columns with this function name so the
   *  monomer cell renderer can resolve display symbols from raw position values.
   *  Return null/undefined if no custom canonicalization is needed. */
  readonly monomerCanonicalizerFuncName?: string | null;
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
    // this is terrible, I know, but otherwise this gets put in webworkers and fails due to DG resolution)))
    // @ts-ignore
    if (window?.DG) {
      // @ts-ignore
      const constFuncs: any[] = window.DG.Func.find({meta: {role: NOTATION_PROVIDER_CONSTRUCTOR_ROLE}});
      return Promise.all(constFuncs.map((f) => f.apply({})));
    }
    return [];
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
