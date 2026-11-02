import * as DG from 'datagrok-api/dg';

import {
  AsmBatch, AsmConnection, AsmMonomer, AsmPolymerType, AsmRole, AsmSequence, packBatch, validateAsmSequence,
} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';

import {HelmType, PolymerType} from '../helm/types';
import {HelmTypes, MonomerTypes, PolymerTypes} from '../helm/consts';
import {IMonomerLibBase, Monomer} from '../types/monomer-library';
import {ISeqConnection, ISeqGraphInfo, ISeqSplitted} from '../utils/macromolecule/types';
import {ISeqHandler} from '../utils/macromolecule/seq-handler';
import {ALPHABET, NOTATION} from '../utils/macromolecule/consts';
import {splitterAsHelm} from '../utils/macromolecule/utils';
import {DEOXYRIBOSE_SYMBOL, PHOSPHATE_SYMBOL, RIBOSE_SYMBOL} from '../utils/const';
import {isLibraryMonomer} from './library-json';

/* Input of the macromolecule assembler: any notation to the parsed sequence of chem-meta's macromolecule-assembler
 * (Bio on Crux, PLAN §3, §5.4). Datagrok resolves every symbol to its exact library key (polymer type + symbol, with
 * the library's quirks: RNA R/P to r/p, the cross-type fallback, override libraries); the assembler only reads keys.
 * Rows the assembler cannot build (unknown monomers, ambiguity, ranges, '?' endpoints, BLOB) get a per-row error. */

/** A fixed or ranged repeat of a flat monomer range [start, end) (HELM `A'3'`, `(A.G)'2'`, `A'2-5'`). */
export interface AsmRepeat {
  /** First flat monomer index of the repeated range */
  start: number;
  /** Flat index after the range */
  end: number;
  /** Repeat counts; min === max is a fixed repeat, which expands */
  min: number;
  max: number;
}

/** Optional graphInfo fields a splitter may provide beyond {@link ISeqGraphInfo}: the HELM splitter (hwe,
 * `splitterAsHelm`, PLAN §5.2) fills them, the BILN one sets none. */
export interface IAsmGraphInfoExt extends ISeqGraphInfo {
  connections: (ISeqConnection & { kind?: 'covalent' | 'hbond' })[];
  /** Per flat monomer: HELM_AA, HELM_SUGAR, HELM_BASE, HELM_LINKER, HELM_CHEM, HELM_BLOB */
  biotypes?: HelmType[];
  /** Polymer ids per chain in source order ('PEPTIDE1', 'CHEM1', ...) */
  chainIds?: string[];
  /** Fixed repeats, expanded in place ({@link expandAsmRepeats}); the splitter lists ranges in {@link errors} */
  repeats?: AsmRepeat[];
  /** Problems the splitter found that make the row unbuildable (ambiguity, '?', ranges, ...) */
  errors?: string[];
  /** Hydrogen bonds (kind 'hbond'): the hwe splitter keeps them out of {@link connections} for old consumers */
  hbonds?: (ISeqConnection & { kind?: 'covalent' | 'hbond' })[];
}

export type AsmSequenceRes = { seq: AsmSequence, error: null } | { seq: null, error: string };

export interface AsmInputOptions {
  /** The row's HELM. Required for a custom notation (its notation provider's HELM); for HELM it is the raw cell
   * (a non-empty HELM without monomers is an error). */
  helm?: string;
  /** Shared symbol resolution (and its cache) across rows */
  resolver?: MonomerResolver;
}

class AsmInputError extends Error {}

const ASM_POLYMER_TYPES = new Set<string>([PolymerTypes.PEPTIDE, PolymerTypes.RNA, PolymerTypes.CHEM]);

/** One resolved symbol: an exact library key, or an inline SMILES. */
export type AsmResolvedSymbol = {
  polymerType: AsmPolymerType, symbol: string, inlineSmiles?: string, monomer: Monomer | null,
};

/** Attachment points as inline SMILES writes them: [*:1], [1*], [1*:1], [R1], CXSMILES _R1 labels. */
const INLINE_SMILES_ATTACHMENT_RE = /\[\*:\d+\]|\[\d+\*(?::\d+)?\]|\[R\d+\]|_R\d+/g;
/** A symbol with a repeat suffix, as a splitter that does not read repeats leaves it: A'3', C)'2', A'2-5' */
const REPEAT_TOKEN_RE = /'(\d+)(?:-(\d+))?'$/;
/** An ambiguity token or its fragments: (A,G), (A+G), (A:1.1,G:69.5) split by '.' */
const AMBIGUITY_TOKEN_RE = /^\(|\)$|[,+]/;

/** Resolves symbols through the library's lookup ({@link IMonomerLibBase.getMonomer}) to exact library keys, with
 * a cache. Lookup never registers anything in the library. */
export class MonomerResolver {
  private readonly cache = new Map<string, AsmResolvedSymbol | string>();

  constructor(public readonly monomerLib: IMonomerLibBase) {}

  /** @return {AsmResolvedSymbol | string} the resolved symbol or an error message */
  resolve(polymerType: PolymerType, symbol: string): AsmResolvedSymbol | string {
    const cacheKey = `${polymerType}\u0001${symbol}`;
    let res = this.cache.get(cacheKey);
    if (res === undefined) {
      res = this.resolveInt(polymerType, symbol);
      this.cache.set(cacheKey, res);
    }
    return res;
  }

  private resolveInt(polymerType: PolymerType, symbol: string): AsmResolvedSymbol | string {
    if (!symbol) return `empty monomer symbol in ${polymerType}`;
    const m = this.monomerLib.getMonomer(polymerType, symbol);
    if (isLibraryMonomer(m)) {
      const mPt = m.polymerType as string;
      return {polymerType: (ASM_POLYMER_TYPES.has(mPt) ? mPt : polymerType) as AsmPolymerType, symbol: m.symbol,
        monomer: m};
    }
    if (countInlineSmilesAttachments(symbol) > 0)
      return {polymerType: polymerType as AsmPolymerType, symbol: '', inlineSmiles: symbol, monomer: null};
    const repeatM = REPEAT_TOKEN_RE.exec(symbol);
    if (repeatM) {
      return repeatM[2] !== undefined ?
        `repeat range '${repeatM[1]}-${repeatM[2]}' at '${symbol}' is not supported` :
        `repeat '${symbol}' is not supported by this splitter`;
    }
    if (symbol === '?') return `ambiguous monomer '?' in ${polymerType}`;
    if (AMBIGUITY_TOKEN_RE.test(symbol)) return `ambiguous monomer '${symbol}' is not supported`;
    return `unknown monomer ${polymerType}:${symbol}`;
  }
}

export function countInlineSmilesAttachments(smiles: string): number {
  return smiles.match(INLINE_SMILES_ATTACHMENT_RE)?.length ?? 0;
}

function defaultRole(polymerType: string): AsmRole {
  switch (polymerType) {
  case PolymerTypes.RNA: return AsmRole.LINKER;
  case PolymerTypes.CHEM: return AsmRole.CHEM;
  default: return AsmRole.AA;
  }
}

function roleOfBiotype(biotype: HelmType | undefined): AsmRole | null {
  switch (biotype) {
  case HelmTypes.AA: return AsmRole.AA;
  case HelmTypes.SUGAR: return AsmRole.SUGAR;
  case HelmTypes.BASE: return AsmRole.BASE;
  case HelmTypes.LINKER: return AsmRole.LINKER;
  case HelmTypes.CHEM: return AsmRole.CHEM;
  default: return null;
  }
}

/** The role of an RNA monomer from its chemistry, as the linear path's triples mode classifies it: a Branch monomer
 * is a base, one with R3 (or three R-groups) a sugar, anything else (phosphates, terminal modifiers) a linker. */
function classifyRnaMonomer(res: AsmResolvedSymbol, prevRole: AsmRole | null): AsmRole {
  if (res.monomer) {
    if (res.monomer.monomerType === MonomerTypes.BRANCH) return AsmRole.BASE;
    const rgroups = res.monomer.rgroups ?? [];
    const hasBranchR = rgroups.some((rg) => rg.label === 'R3') || rgroups.length >= 3;
    return hasBranchR ? AsmRole.SUGAR : AsmRole.LINKER;
  }
  const rCount = countInlineSmilesAttachments(res.inlineSmiles ?? '');
  if (rCount >= 3) return AsmRole.SUGAR;
  return rCount === 1 && prevRole === AsmRole.SUGAR ? AsmRole.BASE : AsmRole.LINKER;
}

function toAsmMonomer(res: AsmResolvedSymbol, role: AsmRole, position: number): AsmMonomer {
  const m: AsmMonomer = {symbol: res.symbol, polymerType: res.polymerType, role, position};
  if (res.inlineSmiles) m.inlineSmiles = res.inlineSmiles;
  return m;
}

function resolveOrThrow(resolver: MonomerResolver, polymerType: PolymerType, symbol: string): AsmResolvedSymbol {
  const res = resolver.resolve(polymerType, symbol);
  if (typeof res === 'string') throw new AsmInputError(res);
  return res;
}

/** Expands fixed repeats: the copies keep the source monomers' positions, chain starts and connection indices move
 * to the expanded indices (a connection to a repeated monomer goes to its first copy). */
export function expandAsmRepeats(seq: AsmSequence, repeats: AsmRepeat[]): AsmSequence {
  if (repeats.length === 0) return seq;
  const src = seq.monomers;
  const sorted = [...repeats].sort((a, b) => a.start - b.start);
  for (let i = 0; i < sorted.length; ++i) {
    const r = sorted[i];
    if (r.min !== r.max) throw new AsmInputError(`repeat range '${r.min}-${r.max}' is not supported`);
    if (!(r.min >= 1) || r.start < 0 || r.end <= r.start || r.end > src.length ||
      (i > 0 && r.start < sorted[i - 1].end))
      throw new AsmInputError(`invalid repeat [${r.start}, ${r.end}) x${r.min}`);
  }
  const newIdx = new Int32Array(src.length);
  const monomers: AsmMonomer[] = [];
  let ri = 0;
  for (let i = 0; i < src.length;) {
    const r = ri < sorted.length && sorted[ri].start === i ? sorted[ri] : null;
    if (r) {
      for (let k = r.start; k < r.end; ++k) newIdx[k] = monomers.length + (k - r.start);
      for (let copy = 0; copy < r.min; ++copy)
        for (let k = r.start; k < r.end; ++k) monomers.push({...src[k]});

      i = r.end;
      ++ri;
    } else {
      newIdx[i] = monomers.length;
      monomers.push(src[i]);
      ++i;
    }
  }
  const srcStarts = seq.graphInfo.disjointSeqStarts;
  const disjointSeqStarts = srcStarts.map((s) => newIdx[s]);
  const connections: AsmConnection[] = seq.graphInfo.connections.map((c) => ({
    ...c,
    monomerIndex1: newIdx[srcStarts[c.seqIndex1] + c.monomerIndex1] - disjointSeqStarts[c.seqIndex1],
    monomerIndex2: newIdx[srcStarts[c.seqIndex2] + c.monomerIndex2] - disjointSeqStarts[c.seqIndex2],
  }));
  return {monomers, graphInfo: {...seq.graphInfo, connections, disjointSeqStarts}};
}

/** HELM (and BILN, and custom notations through their HELM): chains, connections and polymer types from graphInfo. */
function helmToAsm(splitted: ISeqSplitted, resolver: MonomerResolver, helm: string | undefined): AsmSequence {
  const n = splitted.length;
  const gi = splitted.graphInfo as IAsmGraphInfoExt | undefined;
  if (n === 0) {
    // HELM hwe cannot read gives no monomers and hwe's error
    if (helm && helm.trim().length > 0) throw new AsmInputError(gi?.errors?.[0] ?? 'no monomers found in the HELM');
    return {monomers: [], graphInfo: {connections: [], disjointSeqStarts: []}};
  }
  if (gi?.errors && gi.errors.length > 0) throw new AsmInputError(gi.errors[0]);
  const starts = (gi?.disjointSeqStarts && gi.disjointSeqStarts.length > 0) ? gi.disjointSeqStarts : [0];
  const chainIds: string[] | undefined = gi?.chainIds;

  const monomers: AsmMonomer[] = new Array(n);
  for (let c = 0; c < starts.length; ++c) {
    const start = starts[c];
    const end = c + 1 < starts.length ? starts[c + 1] : n;
    let prevRole: AsmRole | null = null;
    const chainPt: string = gi?.polymerTypes?.[start] ?? PolymerTypes.PEPTIDE;
    if (chainPt === PolymerTypes.CHEM && end - start > 1)
      throw new AsmInputError(`CHEM chain ${chainIds?.[c] ?? c + 1} has ${end - start} monomers`);
    for (let i = start; i < end; ++i) {
      const pt: string = gi?.polymerTypes?.[i] ?? chainPt;
      if (pt === PolymerTypes.BLOB) throw new AsmInputError('BLOB polymers are not supported');
      if (!ASM_POLYMER_TYPES.has(pt)) throw new AsmInputError(`polymer type '${pt}' is not supported`);
      if (gi?.biotypes?.[i] === HelmTypes.BLOB) throw new AsmInputError('BLOB monomers are not supported');
      if (splitted.isGap(i)) {
        monomers[i] = {symbol: '', polymerType: pt as AsmPolymerType, role: defaultRole(pt), gap: true, position: i};
        continue;
      }
      const res = resolveOrThrow(resolver, pt as PolymerType, splitted.getCanonical(i));
      let role: AsmRole = roleOfBiotype(gi?.biotypes?.[i]) ?? defaultRole(pt);
      if (pt === PolymerTypes.RNA && !gi?.biotypes) role = classifyRnaMonomer(res, prevRole);
      if (role === AsmRole.BASE && prevRole !== AsmRole.SUGAR)
        throw new AsmInputError(`base '${splitted.getOriginal(i)}' at position ${i + 1} follows no sugar`);
      monomers[i] = toAsmMonomer(res, role, i);
      prevRole = role;
    }
  }

  const connections: AsmConnection[] = [...(gi?.connections ?? []), ...(gi?.hbonds ?? [])].map((sc) => {
    const ac: AsmConnection = {
      seqIndex1: sc.seqIndex1, monomerIndex1: sc.monomerIndex1, rGroup1: sc.rGroup1,
      seqIndex2: sc.seqIndex2, monomerIndex2: sc.monomerIndex2, rGroup2: sc.rGroup2,
    };
    if (sc.kind === 'hbond') ac.kind = 'hbond';
    return ac;
  });
  let seq: AsmSequence = {monomers, graphInfo: {connections, disjointSeqStarts: starts.slice()}};
  if (chainIds) seq.graphInfo.chainIds = chainIds.slice();
  if (gi?.repeats && gi.repeats.length > 0) seq = expandAsmRepeats(seq, gi.repeats);
  return seq;
}

/** FASTA and separator: one chain. Peptides map one to one. DNA/RNA expand as the linear path does: every base
 * gets the default sugar (d for DNA, r for RNA) and, after the first nucleotide, a phosphate p before its sugar; all
 * three carry the base's index as position. Gaps stay gaps (no sugar, no phosphate); the backbone continues across
 * them. Unlike the old linear path, a token that is not a base (r, p, ...) is an error, not a base on a sugar. */
function linearToAsm(splitted: ISeqSplitted, seqHandler: ISeqHandler, resolver: MonomerResolver): AsmSequence {
  const n = splitted.length;
  const monomers: AsmMonomer[] = [];
  const alphabet = seqHandler.alphabet;
  if (alphabet === ALPHABET.DNA || alphabet === ALPHABET.RNA) {
    const pt = PolymerTypes.RNA;
    const sugarSymbol = alphabet === ALPHABET.DNA ? DEOXYRIBOSE_SYMBOL : RIBOSE_SYMBOL;
    let sugar: AsmResolvedSymbol | null = null;
    let phosphate: AsmResolvedSymbol | null = null;
    let hasPrev = false;
    for (let i = 0; i < n; ++i) {
      if (splitted.isGap(i)) {
        monomers.push({symbol: '', polymerType: pt, role: AsmRole.LINKER, gap: true, position: i});
        continue;
      }
      const base = resolveOrThrow(resolver, pt, splitted.getCanonical(i));
      const bm = base.monomer;
      if (bm && bm.monomerType && bm.monomerType !== MonomerTypes.BRANCH)
        throw new AsmInputError(`'${splitted.getOriginal(i)}' at position ${i + 1} is not a nucleobase`);
      sugar ??= resolveOrThrow(resolver, pt, sugarSymbol);
      if (hasPrev) {
        phosphate ??= resolveOrThrow(resolver, pt, PHOSPHATE_SYMBOL);
        monomers.push(toAsmMonomer(phosphate, AsmRole.LINKER, i));
      }
      monomers.push(toAsmMonomer(sugar, AsmRole.SUGAR, i));
      monomers.push(toAsmMonomer(base, AsmRole.BASE, i));
      hasPrev = true;
    }
  } else {
    const pt = PolymerTypes.PEPTIDE;
    for (let i = 0; i < n; ++i) {
      if (splitted.isGap(i))
        monomers.push({symbol: '', polymerType: pt, role: AsmRole.AA, gap: true, position: i});
      else
        monomers.push(toAsmMonomer(resolveOrThrow(resolver, pt, splitted.getCanonical(i)), AsmRole.AA, i));
    }
  }
  return {monomers, graphInfo: {connections: [], disjointSeqStarts: monomers.length > 0 ? [0] : []}};
}

/** Turns one split sequence of any notation into the assembler's parsed sequence.
 * - HELM: monomers, chains, connections and polymer types from {@link ISeqSplitted.graphInfo}; RNA roles from the
 *   splitter's biotypes or, without them (a lone sugar-slot monomer in the row), from the library (base, sugar,
 *   linker); repeats via {@link IAsmGraphInfoExt.repeats}.
 * - BILN: as HELM, all PEPTIDE, connections from the cyclisation marks.
 * - custom: the notation provider's HELM ({@link AsmInputOptions.helm}, required), split as HELM.
 * - FASTA / separator: one chain; DNA/RNA expanded to sugar/base/phosphate.
 * @param {ISeqSplitted} splitted the row's split sequence (ignored for a custom notation)
 * @param {ISeqHandler} seqHandler the column's handler (notation, alphabet)
 * @param {IMonomerLibBase} monomerLib the library to resolve symbols with (an override library included)
 * @param {AsmInputOptions} options the row's HELM, a shared resolver
 * @return {AsmSequenceRes} the sequence, or the row's error */
export function toAsmSequence(
  splitted: ISeqSplitted, seqHandler: ISeqHandler, monomerLib: IMonomerLibBase, options?: AsmInputOptions,
): AsmSequenceRes {
  const resolver = options?.resolver ?? new MonomerResolver(monomerLib);
  try {
    let seq: AsmSequence;
    switch (seqHandler.notation) {
    case NOTATION.HELM:
      seq = helmToAsm(splitted, resolver, options?.helm);
      break;
    case NOTATION.CUSTOM: {
      const helm = options?.helm;
      if (helm === undefined) throw new AsmInputError('a custom notation needs the row\'s HELM');
      seq = helmToAsm(splitterAsHelm(helm), resolver, helm);
      break;
    }
    case NOTATION.BILN:
      seq = helmToAsm(splitted, resolver, undefined);
      break;
    default:
      seq = linearToAsm(splitted, seqHandler, resolver);
    }
    const invalid = validateAsmSequence(seq);
    if (invalid) return {seq: null, error: invalid};
    return {seq, error: null};
  } catch (err: any) {
    if (err instanceof AsmInputError) return {seq: null, error: err.message};
    return {seq: null, error: err instanceof Error ? err.message : String(err)};
  }
}

/** Builds the packed batch of the given rows of a macromolecule column. Empty cells give empty sequences (no
 * monomers, no error); rows that cannot be built get {@link AsmBatch.rowError}.
 * @param {ISeqHandler} seqHandler the column's handler
 * @param {ArrayLike<number> | null} rowIdxs rows to take (batch row i is rowIdxs[i]); null for all rows
 * @param {IMonomerLibBase} monomerLib the library to resolve symbols with
 * @return {AsmBatch} */
export function buildSeqBatch(
  seqHandler: ISeqHandler, rowIdxs: ArrayLike<number> | null, monomerLib: IMonomerLibBase,
): AsmBatch {
  const col: DG.Column<string> = seqHandler.column;
  const rowCount = rowIdxs ? rowIdxs.length : col.length;
  const resolver = new MonomerResolver(monomerLib);
  const seqs: (AsmSequence | null)[] = new Array(rowCount);
  const errors: (string | null)[] = new Array(rowCount);
  const notation = seqHandler.notation;
  for (let i = 0; i < rowCount; ++i) {
    const rowIdx = rowIdxs ? rowIdxs[i] : i;
    seqs[i] = null;
    errors[i] = null;
    try {
      const raw = col.get(rowIdx);
      if (!raw) {
        seqs[i] = {monomers: [], graphInfo: {connections: [], disjointSeqStarts: []}};
        continue;
      }
      let res: AsmSequenceRes;
      if (notation === NOTATION.CUSTOM) {
        const helm = seqHandler.getHelm(rowIdx);
        res = toAsmSequence(splitterAsHelm(helm), seqHandler, monomerLib, {helm, resolver});
      } else {
        const helm = notation === NOTATION.HELM ? raw : undefined;
        res = toAsmSequence(seqHandler.getSplitted(rowIdx), seqHandler, monomerLib, {helm, resolver});
      }
      seqs[i] = res.seq;
      errors[i] = res.error;
    } catch (err: any) {
      errors[i] = err instanceof Error ? err.message : String(err);
    }
  }
  return packBatch(seqs, errors);
}
