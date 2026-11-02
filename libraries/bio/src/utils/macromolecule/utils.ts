/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {Vector} from '@datagrok-libraries/utils/src/type-declarations';
import {vectorDotProduct, vectorLength} from '@datagrok-libraries/utils/src/vector-operations';

import {
  Atom, Bond, expandRepeats, formatConnectionPosition, Group, isAmbiguousBond, Mol, parseHelm, ParseResult, RGroupRef,
  ValidationDiagnostic,
} from '@datagrok-libraries/hwe';

import {
  CandidateSimType,
  CandidateType, ISeqConnectionExt, ISeqGraphInfo, ISeqGraphInfoHelm, ISeqRepeat, ISeqSplitted,
  MonomerFreqs,
  SeqColStats, SeqSplittedBase,
  SplitterFunc
} from './types';
import {ALPHABET, Alphabets, candidateAlphabets, GAP_SYMBOL, GapOriginals, monomerRe, NOTATION} from './consts';
import {SeqPalette} from '../../seq-palettes';
import {AminoacidsPalettes} from '../../aminoacids';
import {NucleotidesPalettes} from '../../nucleotides';
import {UnknownSeqPalettes} from '../../unknown';
import {ISeqHelper} from '../seq-helper';

import {ISeqHandler} from './seq-handler';
import {HelmTypes, PolymerTypes} from '../../helm/consts';
import {HelmType, PolymerType} from '../../helm/types';


export class StringListSeqSplitted implements ISeqSplitted {
  get length(): number { return this.mList.length; }

  get graphInfo(): ISeqGraphInfo | undefined {
    return undefined;
  };

  isGap(posIdx: number): boolean {
    return this.getOriginal(posIdx) === this.gapOriginalMonomer;
  }

  /** @param {number} posIdx monomer position 0-based index */
  getCanonical(posIdx: number): string {
    if (this.length <= posIdx)
      throw new Error('Index out of bounds');
    const om = this.mList[posIdx];
    return om !== this.gapOriginalMonomer ? om : GAP_SYMBOL;
  }

  getOriginal(posIdx: number): string {
    if (this.length <= posIdx)
      throw new Error('Index out of bounds');
    return this.mList[posIdx];
  }

  getCanonicalRegion(start: number, end: number): string[] {
    const actStart = Math.min(Math.max(0, start), this.length);
    const actEnd = Math.max(Math.min(this.length, end), 0);
    const len = actEnd - actStart;
    const gap = this.gapOriginalMonomer;
    return new Array(len).fill(null).map((_, i) => { const om = this.mList[i + actStart]; return om === gap ? GAP_SYMBOL : om; });
  }

  getOriginalRegion(start: number, end: number): string[] {
    const actStart = Math.min(Math.max(0, start), this.length);
    const actEnd = Math.max(Math.min(this.length, end), 0);
    const len = actEnd - actStart;
    return new Array(len).fill(null).map((_, i) => this.mList[i + actStart]);
  }
  get gapOriginal(): string {
    return this.gapOriginalMonomer;
  }

  constructor(
    private readonly mList: SeqSplittedBase,
    private readonly gapOriginalMonomer: string
  ) {}
}


/** Same as string list splitted, but also encodes the graph info (chains, connections), built by
 * {@link parseConnections} on first access */
export abstract class HelmSplitted extends StringListSeqSplitted {
  protected _graphInfo: ISeqGraphInfo | undefined = undefined;
  get graphInfo(): ISeqGraphInfo {
    this._graphInfo ??= this.parseConnections();
    return this._graphInfo;
  }

  protected abstract parseConnections(): ISeqGraphInfo;
}

/** Gets canonical monomers for original ones with cyclization marks */
export class BilnSeqSplitted extends HelmSplitted {
  override getCanonical(posIdx: number): string {
    if (this.isGap(posIdx)) return GAP_SYMBOL;
    const om = this.getOriginal(posIdx);
    let cmRes = om;
    const re = /\(\d{1,2},\d{1,2}\)$/; // there can be multiple cyclization marks
    while (re.test(cmRes))
      cmRes = cmRes.replace(re, '');
    return cmRes;
  }

  override getOriginal(posIdx: number): string {
    const om = super.getOriginal(posIdx);
    // replace [] with '' for biln
    return om.replace(/^\[(.*)\]/, '$1');
  }

  /** @param {string[][]} mListSeparated monomers of each disjoint part ('.'-separated in BILN) */
  constructor(protected readonly mListSeparated: string[][], gapOriginalMonomer: string) {
    super(mListSeparated.flat(), gapOriginalMonomer);
  }

  protected override parseConnections(): ISeqGraphInfo {
    // create graph info for cyclization marks
    // EG: A(1,2).G(1,1)(2,2).C(2,1) -> first number is the bond identifier, second is the R group (1-based)
    // . is used as separator between different sequences
    const graphInfo: ISeqGraphInfo = {connections: [], disjointSeqStarts: []};
    const connectionStorage: {[bondId: string]: {seqIndex: number, monomerIndex: number, rGroup: number}} = {};
    let seqLengthCounter = 0;
    for (let i = 0; i < this.mListSeparated.length; i++) {
      const ss = this.mListSeparated[i];
      graphInfo.disjointSeqStarts.push(seqLengthCounter);
      seqLengthCounter += ss.length;
      for (let j = 0; j < ss.length; j++) {
        const om = ss[j];
        const re = /\(\d{1,2},\d{1,2}\)$/; // there can be multiple cyclization marks
        let cmRes = om;
        while (re.test(cmRes)) {
          const con = /\(\d{1,2},\d{1,2}\)$/.exec(cmRes)?.[0];
          cmRes = cmRes.replace(re, '');
          if (con != null) {
            const nums = con.substring(1, con.length - 1).split(',').map((s) => parseInt(s, 10));
            if (nums.length !== 2 || nums.some((n) => n == null || Number.isNaN(n)))
              continue; // this should never happen, but just in case
            const bondId = nums[0].toString();
            const rGroup = nums[1];
            if (bondId in connectionStorage) {
              // there is already a monomer with this bondId
              const prev = connectionStorage[bondId];
              graphInfo.connections.push({
                seq1Type: PolymerTypes.PEPTIDE, // for biln, its always peptides
                seq2Type: PolymerTypes.PEPTIDE,
                seqIndex1: prev.seqIndex,
                seqIndex2: i,
                monomerIndex1: prev.monomerIndex,
                monomerIndex2: j,
                rGroup1: prev.rGroup,
                rGroup2: rGroup
              });
              delete connectionStorage[bondId];
            } else
              connectionStorage[bondId] = {seqIndex: i, monomerIndex: j, rGroup: rGroup}; // store for later
          }
        }
      }
    }
    return graphInfo;
  }
}


export class FastaSimpleSeqSplitted implements ISeqSplitted {
  get length(): number { return this.seqS.length; }

  isGap(posIdx: number): boolean {
    return this.getOriginal(posIdx) === GapOriginals[NOTATION.FASTA];
  }

  getCanonical(posIdx: number): string {
    if (this.length <= posIdx)
      throw new Error('Index out of bounds');
    if (this.isGap(posIdx))
      return GAP_SYMBOL;
    return this.seqS[posIdx];
  }

  getOriginal(posIdx: number): string {
    if (this.length <= posIdx)
      throw new Error('Index out of bounds');
    return this.seqS[posIdx];
  }

  getCanonicalRegion(start: number, end: number): string[] {
    const actStart = Math.max(0, start);
    const actEnd = Math.min(this.length, end);
    const len = actEnd - actStart;
    const gap = GapOriginals[NOTATION.FASTA];
    return new Array(len).fill(null).map((_, i) => {
      const om = this.seqS[i + actStart];
      return om === gap ? GAP_SYMBOL : om;
    });
  }

  getOriginalRegion(start: number, end: number): string[] {
    const actStart = Math.max(0, start);
    const actEnd = Math.min(this.length, end);
    return this.seqS.slice(actStart, actEnd) as unknown as string[]; // hacky way but quite effective
  }

  get gapOriginal(): string {
    return GapOriginals[NOTATION.FASTA];
  }

  constructor(
    private readonly seqS: string
  ) {}
}

/** Stats of sequences with specified splitter func, returns { freq, sameLength }.
 * @param {DG.Column} seqCol
 * @param {number} minLength
 * @param {SplitterFunc} splitter
 * @return { SeqColStats }, sameLength: boolean } stats of column sequences
 */
export function getStatsForCol(seqCol: DG.Column<string>, minLength: number, splitter: SplitterFunc): SeqColStats {
  const cats = seqCol.categories;
  const splitted: Iterable<ISeqSplitted> = wu.enumerate(seqCol.getRawData())
    .map(([catI, rowIdx]) => splitter(cats[catI]));
  return getStats(splitted, minLength);
}

function getStats(splitted: Iterable<ISeqSplitted>, minLength: number): SeqColStats {
  const freq: { [m: string]: number } = {};
  let sameLength = true;
  let firstLength = null;

  for (const seqSS of splitted) {
    if (firstLength == null)
      firstLength = seqSS.length;
    else if (seqSS.length !== firstLength)
      sameLength = false;

    if (seqSS.length >= minLength) {
      for (let posIdx = 0; posIdx < seqSS.length; ++posIdx) {
        const cm = seqSS.getCanonical(posIdx);
        if (!(cm in freq))
          freq[cm] = 0;
        freq[cm] += 1;
      }
    }
  }
  return {freq: freq, sameLength: sameLength};
}

/** Split sequence for single character monomers, square brackets multichar monomer names or gap symbol.
 * @param seq object with sequence
 * @param getMonomer Source of the {@link seq} string
 * @return {string[]} array of monomers
 */
export const splitterAsFasta: SplitterFunc = (seq: string): ISeqSplitted => {
  const mmList = wu<RegExpMatchArray>(seq.toString().matchAll(monomerRe))
    .map((ma: RegExpMatchArray) => {
      return ma[2] ?? ma[1]; // preserve '-' as gap symbol for compatibility with simpleAsFastaSimple
    }).toArray();

  return new StringListSeqSplitted(mmList, GapOriginals[NOTATION.FASTA]);


  // return new Proxy(splittedList as object, {
  //   get(target: string[], p: string | symbol, receiver: any): any {
  //     const k = 11;
  //   }
  // }) as ISeqSplitted;
};

export const splitterAsFastaSimple: SplitterFunc = (seq: string): ISeqSplitted => {
  return !!seq ? new FastaSimpleSeqSplitted(seq) : new StringListSeqSplitted([], GapOriginals[NOTATION.FASTA]);
};

/** Gets method to split sequence by separator
 * @param separator Monomer separator
 * @param limit
 * @return {SplitterFunc}
 */
export function getSplitterWithSeparator(separator: string, limit: number | undefined = undefined): SplitterFunc {
  return (seq: string) => {
    if (!seq)
      return new StringListSeqSplitted([], GapOriginals[NOTATION.SEPARATOR]);
    else {
      let mmList: string[];
      const mRe = new RegExp(`(?<=^|\\${separator})("-"|'-'|[^\\${separator}]*)(?=\\${separator}|$)`, 'g'); // depends on separator args
      if (limit != undefined) {
        mRe.lastIndex = 0;
        mmList = wu(seq.matchAll(mRe)).take(limit).map((ea) => ea[0]).toArray();
      } else
        mmList = seq.replaceAll('\"-\"', '').replaceAll('\'-\'', '').split(separator, limit);

      return new StringListSeqSplitted(mmList, GapOriginals[NOTATION.SEPARATOR]);
    }
  };
}

/** Non-implicit hwe bond in flat positions; graphInfo is built from these on first access. `ambiguous` is the
 * connection as written when its ends are not one monomer each (hwe's `isAmbiguousBond`). */
type HweLink = { a1: number, a2: number, r1: RGroupRef, r2: RGroupRef, kind: string, ambiguous?: string };

/** What {@link HweHelmSplitted} keeps of an hwe parse to build its graph info on first access */
type HweRow = {
  monomerPolymerTypes: PolymerType[],
  chainStarts: number[],
  chainIds: string[],
  links: HweLink[],
  biotypes: HelmType[],
  biotypesKnown: boolean,
  /** The parsed molecule, kept only when it has repeat brackets: they expand on first access */
  repeatMol: Mol | null,
  diagnostics: readonly ValidationDiagnostic[],
  ambiguities: string[] | null,
  /** An hwe parse the adapter cannot map to positions (not expected from hwe): the row is empty, with this error */
  unmapped?: string,
};

/** HELM split by the hwe parser (PLAN §5.2, D4): one position per hwe atom in source order (RNA sugar, base and
 * linker are separate positions), position = `continuousId - 1`. The extended graph info is built on first access. */
class HweHelmSplitted extends HelmSplitted {
  constructor(mList: string[], private readonly row: HweRow) {
    super(mList, GapOriginals[NOTATION.HELM]);
  }

  protected override parseConnections(): ISeqGraphInfo {
    const row = this.row;
    const starts = row.chainStarts;
    const gi: ISeqGraphInfoHelm = {connections: [], disjointSeqStarts: starts, polymerTypes: row.monomerPolymerTypes};
    const chainOf = (pos: number): number => {
      let c = starts.length - 1;
      while (c > 0 && starts[c] > pos) --c;
      return c;
    };
    const errors: string[] = row.unmapped !== undefined ? [row.unmapped] : [];
    if (row.ambiguities) errors.push(...row.ambiguities);
    const hbonds: ISeqConnectionExt[] = [];
    for (const l of row.links) {
      // Multi-site positions, ranges, '?', symbol positions found more than once, group endpoints: hwe draws such a
      // bond at the first candidate, which is not a bond of the molecule; it is in neither list (see
      // ISeqGraphInfoHelm.connections).
      if (l.ambiguous !== undefined) {
        errors.push(`ambiguous connection '${l.ambiguous}' is not supported`);
        continue;
      }
      const c1 = chainOf(l.a1);
      const c2 = chainOf(l.a2);
      if (l.kind === 'hbond') {
        hbonds.push({
          seq1Type: gi.polymerTypes![l.a1], seq2Type: gi.polymerTypes![l.a2], seqIndex1: c1, seqIndex2: c2,
          monomerIndex1: l.a1 - starts[c1], monomerIndex2: l.a2 - starts[c2], rGroup1: 0, rGroup2: 0, kind: 'hbond',
        });
        continue;
      }
      if (l.kind !== 'cross-polymer' || l.r1.kind !== 'numeric' || l.r2.kind !== 'numeric') {
        const desc = `${row.chainIds[c1]},${row.chainIds[c2]},${l.a1 - starts[c1] + 1}-${l.a2 - starts[c2] + 1}`;
        errors.push(`connection ${desc} with R-groups '${rGroupRefToString(l.r1)}', '${rGroupRefToString(l.r2)}' ` +
          `is not supported`);
        continue;
      }
      gi.connections.push({
        seq1Type: gi.polymerTypes![l.a1], seq2Type: gi.polymerTypes![l.a2], seqIndex1: c1, seqIndex2: c2,
        monomerIndex1: l.a1 - starts[c1], monomerIndex2: l.a2 - starts[c2], rGroup1: l.r1.r, rGroup2: l.r2.r,
        kind: 'covalent',
      });
    }
    // HELM has no DNA polymer; hwe reads DNA1{..} as RNA, but its monomers are not RNA library keys
    for (const id of row.chainIds)
      if (id.startsWith('DNA')) errors.push(`polymer type 'DNA' is not supported (${id})`);
    // Lenient-parse errors that change the graph (dropped connections, unterminated braces, CHEM repeats, ...);
    // missing '$' sections and the trailing tokens of V1.0 strings are harmless. Whitespace, multi-site positions
    // and group endpoints are warnings and infos in hwe 1.1.0.
    const recovered = row.unmapped !== undefined || row.diagnostics.some(isHweRecoveryError);
    const repeatTooLarge = row.diagnostics.some((d) => d.code === 'bracket-repeat-too-large');
    errors.push(...hweDiagnosticErrors(row.diagnostics));

    gi.hbonds = hbonds;
    gi.chainIds = row.chainIds;
    if (row.biotypesKnown) gi.biotypes = row.biotypes;
    gi.repeats = row.repeatMol && !repeatTooLarge ? hweRepeats(row.repeatMol, this.length, errors) : [];
    gi.errors = errors;
    if (recovered) gi.recovered = true;
    return gi;
  }
}

/** An empty row with the error {@link unmapped} and hwe's diagnostics: an hwe parse whose atoms cannot be mapped to
 * positions (not expected from hwe) */
function hweUnmappedRow(unmapped: string, diagnostics: readonly ValidationDiagnostic[]): HweHelmSplitted {
  return new HweHelmSplitted([], {monomerPolymerTypes: [], chainStarts: [], chainIds: [], links: [], biotypes: [],
    biotypesKnown: false, repeatMol: null, diagnostics, ambiguities: null, unmapped});
}

/** Lenient-parse error codes of hwe that leave the parsed graph right */
const HWE_HARMLESS_ERRORS = new Set<string>(['missing-section-separator']);
/** hwe warnings whose graph is not what the HELM says: a zero repeat count (the monomers stay once), a repeat over
 * 10000 (not expanded), two polymers with one id (connections go to the last) */
const HWE_ERROR_WARNINGS = new Set<string>(['repeat-zero-count', 'bracket-repeat-too-large', 'duplicate-polymer-id']);

/** An hwe error after which the lenient parser skipped part of the string (an unterminated brace, an unexpected
 * token, ...): monomers after it may be missing from the parse */
function isHweRecoveryError(d: ValidationDiagnostic): boolean {
  return d.severity === 'error' && !HWE_HARMLESS_ERRORS.has(d.code);
}

/** Row errors from hwe's diagnostics: the recovery errors and the warnings of {@link HWE_ERROR_WARNINGS} */
function hweDiagnosticErrors(diagnostics: readonly ValidationDiagnostic[]): string[] {
  const errors: string[] = [];
  for (const d of diagnostics) {
    if (isHweRecoveryError(d) || HWE_ERROR_WARNINGS.has(d.code))
      errors.push(`HELM ${d.code}: ${d.message}`);
  }
  return errors;
}

function rGroupRefToString(r: RGroupRef): string {
  switch (r.kind) {
  case 'numeric': return `R${r.r}`;
  case 'positional': return `${r.position}:R${r.r}`;
  case 'unspecified': return '?';
  case 'wildcard': return '*';
  case 'pair': return 'pair';
  }
}

/** A connection with {@link Bond.endpoints} as the HELM wrote it, e.g. `PEPTIDE1,CHEM1,(1,3):R3-1:R1` */
function hweConnectionToString(b: Bond): string {
  const [e1, e2] = b.endpoints!;
  return `${e1.polymerId},${e2.polymerId},${formatConnectionPosition(e1.position)}:${rGroupRefToString(b.r1)}-` +
    `${formatConnectionPosition(e2.position)}:${rGroupRefToString(b.r2)}`;
}

/** True when every end of an ambiguous-looking bond is one monomer after all: an index, or a symbol position
 * (`C:R3`) whose symbol occurs once in its polymer, where hwe draws the bond. Group ids, `?`, lists and ranges are
 * never definite. */
function isDefiniteBySymbol(b: Bond, atoms: readonly Atom[]): boolean {
  for (const e of b.endpoints!) {
    if (e.isGroup) return false;
    const p = e.position;
    if (p.kind === 'index') continue;
    if (p.kind !== 'symbol') return false;
    let count = 0;
    for (const a of atoms) if (a.symbol === p.value && a.bio.chainId === e.polymerId) ++count;
    if (count !== 1) return false;
  }
  return true;
}

const repeatCountToString = (r: { min: number, max: number }): string =>
  r.min === r.max ? `'${r.min}'` : `'${r.min}-${r.max}'`;

/** The fixed repeats of {@link mol} as flat source ranges [start, end) for assembly-input's `expandAsmRepeats`, which
 * copies each range in place, the copies keeping the source positions (D9). Which repeats are fixed comes from hwe's
 * `expandRepeats`; the ranges are checked against its expansion (`bio.sourceContinuousId` of every expanded atom).
 * Ranges `'2-5'`, CHEM repeats, nested repeats and zero counts add to {@link errors}.
 * @param {Mol} mol hwe parse with repeat brackets
 * @param {number} n flat monomer count (the atoms of {@link mol})
 * @param {string[]} errors the row's errors, appended to
 * @return {ISeqRepeat[]} the fixed repeats, in source order */
function hweRepeats(mol: Mol, n: number, errors: string[]): ISeqRepeat[] {
  const atomOfId = new Map<number, Atom>();
  for (const a of mol.atoms) atomOfId.set(a.id, a);
  const spanOf = (g: Group): { start: number, end: number } => {
    let start = n;
    let end = 0;
    for (const m of g.members) {
      const a = atomOfId.get(parseInt(m.target, 10));
      if (a === undefined) continue;
      start = Math.min(start, a.bio.continuousId - 1);
      // An RNA sugar member brings its base (the next position)
      end = Math.max(end, a.bio.continuousId, a.branchAtom?.bio.continuousId ?? 0);
    }
    return {start, end};
  };
  const {mol: expanded, unexpanded} = expandRepeats(mol);
  const repeats: ISeqRepeat[] = [];
  for (const g of mol.groups) {
    if (g.kind !== 'bracket' || g.repeat === undefined) continue;
    const {start, end} = spanOf(g);
    const where = end > start ? ` at position ${start + 1}` : '';
    if (unexpanded.includes(g)) {
      errors.push(g.repeat.min !== g.repeat.max ? `repeat range ${repeatCountToString(g.repeat)}${where} is not supported` :
        `repeat ${repeatCountToString(g.repeat)}${where} is not supported`);
    } else if (!(g.repeat.min >= 1))
      errors.push(`repeat ${repeatCountToString(g.repeat)}${where} is not supported`);
    else if (end > start)
      repeats.push({start, end, min: g.repeat.min, max: g.repeat.max});
  }
  repeats.sort((x, y) => x.start - y.start);
  for (let i = 1; i < repeats.length; ++i) {
    if (repeats[i].start < repeats[i - 1].end) {
      errors.push(`nested repeats at position ${repeats[i - 1].start + 1} are not supported`);
      return repeats;
    }
  }
  // The expansion the ranges give, in source positions, must be hwe's
  const expected: number[] = [];
  for (let i = 0, ri = 0; i < n;) {
    const r = ri < repeats.length && repeats[ri].start === i ? repeats[ri++] : null;
    if (r === null) { expected.push(i++); continue; }
    for (let copy = 0; copy < r.min; ++copy) for (let k = r.start; k < r.end; ++k) expected.push(k);
    i = r.end;
  }
  const got: number[] = new Array(expanded.atoms.length);
  for (const a of expanded.atoms) got[a.bio.continuousId - 1] = (a.bio.sourceContinuousId ?? a.bio.continuousId) - 1;
  if (got.length !== expected.length || got.some((p, i) => p !== expected[i]))
    errors.push(`repeats at positions ${repeats.map((r) => r.start + 1).join(', ')} expand differently in hwe`);
  return repeats;
}

/** Parses HELM with hwe, leniently and without a monomer library (no resolution, no library side effects) */
function hweParseHelm(seq: string): ParseResult {
  return parseHelm(seq, null, {mode: 'lenient'});
}

/** Splits HELM to monomers with the hwe parser (PLAN §5.2, D4), keeping the {@link ISeqSplitted} contract:
 * one position per hwe atom in source order, so RNA `R(A)P` gives the positions R, A, P (position = `continuousId - 1`
 * of hwe); originals without the enclosing square brackets; `*` is the gap. Repeats are listed once, as written
 * (`A.(G.C)'2'` is A, G, C): renderers and WebLogo see the source tokens, and only the assembly input expands them,
 * from `graphInfo.repeats`. `graphInfo` (built on first access) holds the section-1 connections with the chain index
 * in source order and any position or R-group number, chain starts, per-monomer polymer types, and the
 * {@link ISeqGraphInfoHelm} extras: hydrogen bonds, chain ids, biotypes, repeats and the errors that make the row
 * unbuildable at atomic level (ambiguity, ambiguous connections, repeat ranges, ...).
 * Invalid HELM keeps what hwe read: a row hwe cannot parse at all has no monomers, and one hwe recovered from a syntax
 * error in may miss the monomers after it (hwe skips the rest of the polymer, so `PEPTIDE1{Ac(1).A.A}$$$$` reads as
 * `Ac` alone). Both carry hwe's errors in `graphInfo.errors` and `graphInfo.recovered`, which keep the row
 * unbuildable at atomic level.
 * EG byciclic helm: PEPTIDE1{F}|PEPTIDE2{[dI].[Trp_Ome].[Asp_OMe].[Cys_Bn].[meG].[Phe_3Cl].[dD].T.[dI].T.[dK].[aG].[3Pal].[xiIle].[meD].[Ala_tBu]}|PEPTIDE3{L.[Pro_4Me3OH].S.[NMe2Abz].Q.[3Pal].[xiIle].[D-Hyp].[Ala_tBu].[dI].[Trp_Ome].[Asp_OMe].N.[meG].[Phe_34diCl].[Phe_34diCl]}$PEPTIDE2,PEPTIDE2,16:R2-1:R1|PEPTIDE3,PEPTIDE3,16:R2-1:R1|PEPTIDE3,PEPTIDE2,10:R3-1:R3|PEPTIDE1,PEPTIDE2,1:R2-9:R3$$$V2.0
 * Also supports smiles in monomer names, like: PEPTIDE1{[ac].L.N.[*N[C@H](C(=O)*)Cc1ccc(cc1)OP(=O)(O)O |$_R1;;;;;_R2;;;;;;;;;;;;$|].I.D.L.D.L.V.[am]}$$$$
 * @param {string} seq Source string of HELM notation
 * @return {ISeqSplitted}
 */
export const splitterAsHelm: SplitterFunc = (seq: string): ISeqSplitted => {
  const {mol, diagnostics} = hweParseHelm(seq);
  const atoms = mol.atoms;
  const chains = mol.chains;
  const chainCount = chains.length;
  const n = atoms.length;
  const mList: string[] = new Array(n);
  const monomerPolymerTypes: PolymerType[] = new Array(n);
  const biotypes: HelmType[] = new Array(n);
  const chainIds: string[] = new Array(chainCount);
  const chainStarts: number[] = new Array(chainCount);
  // Atoms come chain after chain, in the order of the chains: a chain covers its backbone atoms and their branches
  // (RNA bases). Sizes, not ids, cut the chains, so two chains with one id (invalid, an error below) stay two.
  let total = 0;
  for (let c = 0; c < chainCount; ++c) {
    chainIds[c] = chains[c].id;
    chainStarts[c] = total;
    for (const a of chains[c].backboneAtoms) total += a.branchAtom === undefined ? 1 : 2;
  }
  if (total !== n) return hweUnmappedRow('HELM parse: hwe chain sizes do not add up to its monomers', diagnostics);
  let biotypesKnown = true;
  let ambiguities: string[] | null = null;
  let inOrder = true;
  let c = -1;
  let chainId = '';
  let chainPt: PolymerType = PolymerTypes.PEPTIDE;
  for (let i = 0; i < n; ++i) {
    const a = atoms[i];
    const bio = a.bio;
    if (bio.continuousId !== i + 1) inOrder = false;
    // a chain without atoms (G1{}) starts where the next one does
    while (c + 1 < chainCount && chainStarts[c + 1] <= i) {
      ++c;
      chainId = chainIds[c];
      chainPt = chains[c].polymerType;
    }
    if (bio.chainId !== chainId) inOrder = false;
    mList[i] = a.symbol;
    monomerPolymerTypes[i] = chainPt;
    biotypes[i] = a.biotype;
    // hwe (1.1.0 too) reads a lone RNA token (a 5' phosphate `P.R(A)P`, a spacer) as a sugar slot without a base
    if (a.biotype === HelmTypes.SUGAR && a.branchAtom === undefined) biotypesKnown = false;
    if (bio.ambiguity !== undefined)
      (ambiguities ??= []).push(`ambiguous monomer '${bio.ambiguity}' at position ${i + 1} is not supported`);
  }
  // the position invariant or the chain cuts would not hold; not expected from hwe
  if (!inOrder) return hweUnmappedRow('HELM parse: hwe chains out of order', diagnostics);

  const links: HweLink[] = [];
  for (const b of mol.bonds) {
    if (b.kind === 'backbone' || b.kind === 'branch') continue;
    const l: HweLink = {a1: b.a1.bio.continuousId - 1, a2: b.a2.bio.continuousId - 1, r1: b.r1, r2: b.r2, kind: b.kind};
    if (b.endpoints !== undefined && isAmbiguousBond(b) && !isDefiniteBySymbol(b, atoms))
      l.ambiguous = hweConnectionToString(b);
    links.push(l);
  }
  const repeatMol = mol.groups.some((g) => g.kind === 'bracket' && g.repeat !== undefined) ? mol : null;
  return new HweHelmSplitted(mList, {monomerPolymerTypes, chainStarts, chainIds, links, biotypes, biotypesKnown,
    repeatMol, diagnostics, ambiguities});
};

export function splitterAsBiln(seq: string): ISeqSplitted {
  const disjointPartsRaw = seq.split('.')
    .map((ds) => ds.replaceAll('\"-\"', '').replaceAll('\'-\'', ''));
    // some monomers might contain '-', and need to be enclosed in [], e.g. A-[D-Arg]-C
    // replace all separators inside [] with some temporary string, then split by separator, then restore
  const tempSep = '___TEMP___';
  for (let i = 0; i < disjointPartsRaw.length; i++) {
    const re = /\[([^\]]+)\]/g;
    let match;
    let dp = disjointPartsRaw[i];
    while ((match = re.exec(disjointPartsRaw[i])) !== null) {
      const inside = match[1];
      const replaced = inside.replaceAll('-', tempSep);
      dp = dp.replace(inside, replaced);
    }
    disjointPartsRaw[i] = dp;
  }
  const disjointParts: string[][] = disjointPartsRaw
    .map((ds) => ds.split('-')
      .map((m) => m.replaceAll(tempSep, '-')));

  //const totalLength = disjointParts.reduce((acc, ss) => acc + ss.length, 0);
  const res = new BilnSeqSplitted(
    disjointParts,
    GapOriginals[NOTATION.SEPARATOR]);

  return res;
}


/** Func type to shorten a {@link monomerLabel} with length {@link limit} */
export type MonomerToShortFunc = (monomerLabel: string, limit: number) => string;

/** Get splitter method to split sequences to monomers (required for MacromoleculeDifferenceCellRenderer)
 * @param {string} units
 * @param {string} separator
 * @param limit
 * @return {SplitterFunc}
 */
export function getSplitter(units: string, separator: string, limit: number | undefined = undefined): SplitterFunc {
  if (units.toLowerCase().startsWith(NOTATION.FASTA))
    return splitterAsFasta;
  else if (units.toLowerCase().startsWith(NOTATION.SEPARATOR))
    return getSplitterWithSeparator(separator, limit);
  else if (units.toLowerCase().startsWith(NOTATION.HELM))
    return splitterAsHelm;
  else if (units.toLowerCase().startsWith(NOTATION.BILN))
    return splitterAsBiln;
  else
    throw new Error(`Unexpected units ${units} .`);

  // TODO: Splitter for HELM
}

const longMonomerPartRe: RegExp = /([^\W_]+)/g;

/** Convert long monomer names to short ones */
export function monomerToShort(amino: string, maxLengthOfMonomer: number): string {
  if (amino.length <= maxLengthOfMonomer)
    return amino;
  //const kebabAmino = amino.replace(/[A-Z]+(?![a-z])|[A-Z]/g, ($, ofs) => (ofs ? '-' : '') + $);
  const shortAminoMatch: RegExpMatchArray | null = amino.match(longMonomerPartRe);
  const needAddDots: boolean = amino.length > maxLengthOfMonomer || (shortAminoMatch?.length ?? 0) > 1;
  const shortAmino = shortAminoMatch?.[0] ?? ' ';
  return !needAddDots ? shortAmino : shortAmino.substring(0, maxLengthOfMonomer - 1) + '…';
}

/** */
export function getAlphabet(alphabet: ALPHABET): Set<string> {
  switch (alphabet) {
  case ALPHABET.DNA:
    return Alphabets.fasta.dna;
  case ALPHABET.RNA:
    return Alphabets.fasta.rna;
  case ALPHABET.PT:
    return Alphabets.fasta.peptide;
  default:
    throw new Error(`Unsupported alphabet '${alphabet}'.`);
  }
}

/** Calculate similarity in current sequence and alphabet.
 * @param {MonomerFreqs} freq
 * @param {Set<string>} alphabet
 * @param {string} gapSymbol
 * @return {number} Cosine similarity
 */
export function getAlphabetSimilarity(freq: MonomerFreqs, alphabet: Set<string>, gapSymbol: string = '-'): number {
  const keys = new Set<string>([...new Set(Object.keys(freq)), ...alphabet]);
  keys.delete(gapSymbol);

  const freqA: number[] = [];
  const alphabetA: number[] = [];
  for (const m of keys) {
    freqA.push(m in freq ? freq[m] : 0);
    alphabetA.push(alphabet.has(m) ? 1 : 0);
  }
  /* There were a few ideas: chi-squared, pearson correlation (variance?), scalar product */
  const freqV: Vector = new Vector(freqA);
  const alphabetV: Vector = new Vector(alphabetA);
  return vectorDotProduct(freqV, alphabetV) / (vectorLength(freqV) * vectorLength(alphabetV));
}

/** From detectMacromolecule */
export function detectAlphabet(freq: MonomerFreqs, candidates: CandidateType[], gapSymbol: string = '-') {
  const candidatesSims: CandidateSimType[] = candidates.map((c) => {
    const sim = getAlphabetSimilarity(freq, c.alphabet, gapSymbol);
    return new CandidateSimType(c, freq, sim);
  });

  let alphabetName: string;
  const maxSim = Math.max(...candidatesSims.map(
    (cs) => cs.similarity > cs.cutoff ? cs.similarity : -1));
  if (maxSim > 0) {
    const sim = candidatesSims.find((cs) => cs.similarity === maxSim)!;
    alphabetName = sim.name;
  } else
    alphabetName = ALPHABET.UN;
  return alphabetName;
}

export function detectHelmAlphabet(freq: MonomerFreqs, candidates: CandidateType[], gapSymbol: string = '-') {
  // helm can contain DNA/RNA and other monomers. need to check that
  const mons = Object.keys(freq);
  const helmNucleotideFullRe = /\(|\)/;
  // heuristic for detecting dna/rna
  const hasBranching = mons.filter((m) => m.split(helmNucleotideFullRe).filter((p) => !!p).length === 3).length > mons.length * 0.8;
  const correctedFreqs = hasBranching ? Object.entries(freq)
    .reduce((acc, [m, f]) => {
      const split = m.split(helmNucleotideFullRe);
      const actualMonomer = split[1]; // the inside part is the actual monomer (for example r(A)p)
      if (actualMonomer)
        acc[actualMonomer] = f; return acc;
    },
      {} as Record<string, number>) : freq;
  return detectAlphabet(correctedFreqs, candidates, gapSymbol);
}

/** Selects a suitable palette based on column data
 * @param {DG.Column} seqCol Column to look for a palette
 * @param {number}  minLength minimum length of sequence to detect palette (empty strings are allowed)
 * @return {SeqPalette} Palette corresponding to the alphabet of the sequences in the column
 */
export function pickUpPalette(seqCol: DG.Column, seqHelper: ISeqHelper, minLength: number = 5): SeqPalette {
  let alphabet: string;
  if (seqCol.semType == DG.SEMTYPE.MACROMOLECULE) {
    const sh: ISeqHandler = seqHelper.getSeqHandler(seqCol);
    alphabet = sh.alphabet;
  } else {
    const stats: SeqColStats = getStatsForCol(seqCol, minLength, splitterAsFasta);
    alphabet = detectAlphabet(stats.freq, candidateAlphabets);
  }

  const res = getPaletteByType(alphabet);
  return res;
}

export function getPaletteByType(paletteType: string): SeqPalette {
  switch (paletteType) {
  case ALPHABET.PT:
    return AminoacidsPalettes.GrokGroups;
  case ALPHABET.DNA:
  case ALPHABET.RNA:
    return NucleotidesPalettes.Chromatogram;
    // other
  default:
    return UnknownSeqPalettes.Color;
  }
}

export function pickUpSeqCol(df: DG.DataFrame): DG.Column<string> | null {
  const semTypeColList = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  let resCol: DG.Column | null = semTypeColList.find((col) => {
    const units = col.meta.units;
    return units ? units.indexOf('MSA') !== -1 : false;
  }) ?? null;
  if (!resCol && semTypeColList.length > 0)
    resCol = semTypeColList[0];
  return resCol;
}

export function polymerTypeToHelmType(pt: PolymerType): HelmType {
  switch (pt) {
  case PolymerTypes.PEPTIDE:
    return HelmTypes.AA;
  case PolymerTypes.RNA:
    return HelmTypes.NUCLEOTIDE;
  case PolymerTypes.CHEM:
    return HelmTypes.CHEM;
  case PolymerTypes.BLOB:
    return HelmTypes.BLOB;
  default:
    return HelmTypes.AA;
  }
}
