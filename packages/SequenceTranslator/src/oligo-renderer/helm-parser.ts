/**
 * Lightweight HELM parser dedicated to RNA / DNA duplex visualization.
 *
 * Parses strings of the shape `RNA1{...}|RNA2{...}$$$$` (single or two
 * chain), with each monomer in `sugar(base)phosphate` form (with brackets
 * around multi-char codes). Conjugates / linkers without a base are
 * recognized as standalone units (e.g. `[GalNAc]`, `[L3]`, `[Chol]`).
 *
 * Why a custom parser rather than reusing bio's `HelmHelper.parse()`:
 *  - we want plain `{sugar, base, phosphate}` triples, not a JSDraw2 mol graph
 *  - parsing is on the cell-render hot path; this is ~zero-allocation
 *  - bio's HELM splitter has had several normalization changes recently
 *    (see Apr 2026 commits around RNA triplet splitting)
 */

import {cleanupHelmSymbol} from '@datagrok-libraries/bio/src/helm/utils';
import {
  canonicalPhosphateSymbol, canonicalSugarSymbol, isLinkerSymbol,
  DuplexPair, ParsedConjugate, ParsedDuplex, ParsedLinker, ParsedMonomer, ParsedNucleotide, ParsedStrand,
} from './types';

/** Parse the whole HELM string. Tolerant of whitespace and missing `$$$$`.
 *
 * Beyond the polymer bodies, this also reads the HELM `$connections$` and
 * extended-annotations sections when present, and uses them as the *source of
 * truth* for strand roles and base-pair register:
 *   - `strandtype` annotations (`ss` / `as`) decide which chain is the sense
 *     strand (the parser swaps chains so `sense` is always the sense strand);
 *   - `pair` connections fix the alignment shift between sense and antisense.
 * When that info is absent, `sense` defaults to the first chain and the
 * renderer auto-aligns by complementarity (see `alignment.ts`). */
export function parseHelmDuplex(helm: string): ParsedDuplex {
  const polymerSection = (helm ?? '').split('$')[0] ?? '';
  const chains = polymerSection.split('|').map((c) => c.trim()).filter((c) => c.length > 0);

  const parsed = chains.map(parseChainWithId);
  let senseEntry: ParsedChain = parsed[0] ?? {strand: {type: 'RNA', monomers: []}, id: undefined, body: ''};
  let antiEntry: ParsedChain | null = parsed[1] ?? null;

  const strandTypes = parseHelmAnnotations(helm);

  // Strand directions: if annotations name the sense / antisense strand, make
  // `sense` actually be that strand (swap chains when the HELM disagrees with
  // the default first-chain-is-sense assumption).
  if (antiEntry && Object.keys(strandTypes).length) {
    const senseId = findStrandByType(strandTypes, SENSE_TYPES);
    const antiId = findStrandByType(strandTypes, ANTISENSE_TYPES);
    if ((senseId && senseId === antiEntry.id) || (antiId && antiId === senseEntry.id))
      [senseEntry, antiEntry] = [antiEntry, senseEntry];
  }

  const sense: ParsedStrand = {...senseEntry.strand, id: senseEntry.id};
  const antisense: ParsedStrand | null = antiEntry ? {...antiEntry.strand, id: antiEntry.id} : null;

  const result: ParsedDuplex = {sense, antisense, raw: helm};
  if (Object.keys(strandTypes).length) result.strandTypes = strandTypes;

  // Explicit base pairs → authoritative alignment shift.
  if (antiEntry) {
    const conns = parseHelmConnections(helm).filter((c) => isPairAnnotation(c.srcAnn) || isPairAnnotation(c.tgtAnn));
    if (conns.length) {
      const senseMap = nucleotideOrdinalByHelmPos(senseEntry.body);
      const antiMap = nucleotideOrdinalByHelmPos(antiEntry.body);
      const pairs = connectionsToPairs(conns, senseEntry.id, antiEntry.id, senseMap, antiMap);
      if (pairs.length) {
        result.pairs = pairs;
        result.alignment = {shift: shiftFromPairs(pairs, countNucleotides(antisense!)), source: 'explicit'};
      }
    }
  }
  return result;
}

interface ParsedChain { strand: ParsedStrand; id: string | undefined; body: string; }

function parseChainWithId(chain: string): ParsedChain {
  // RNA1{...}, DNA2{...}, CHEM3{...}
  const match = chain.match(/^((?:RNA|DNA|CHEM|PEPTIDE|BLOB)\d+)\{(.*)\}$/);
  if (!match)
    return {strand: {type: 'CHEM', monomers: []}, id: undefined, body: ''};
  const id = match[1];
  const body = match[2];
  const type = id.replace(/\d+$/, '');
  return {strand: {type, monomers: parseMonomers(body), id}, id, body};
}

/** Split chain contents on top-level dots (depth-aware to preserve `[a.b]`/`(c.d)`). */
function splitSegments(content: string): string[] {
  const out: string[] = [];
  let depth = 0;
  let cur = '';
  for (const ch of content) {
    if (ch === '[' || ch === '(') {
      depth++;
    } else if (ch === ']' || ch === ')') {
      depth--;
    } else if (ch === '.' && depth === 0) {
      if (cur) out.push(cur);
      cur = '';
      continue;
    }
    cur += ch;
  }
  if (cur) out.push(cur);
  return out;
}

function parseMonomers(content: string): ParsedMonomer[] {
  return splitSegments(content).map((seg, i) => parseMonomer(seg, i));
}

/** Parse one monomer text into a structured triple OR conjugate. */
function parseMonomer(s: string, position: number): ParsedMonomer {
  let i = 0;
  let sugar = '';
  let base: string | null = null;
  let phosphate = '';
  const raw = s;

  // sugar (optional brackets)
  if (s[0] === '[') {
    const end = s.indexOf(']');
    sugar = s.substring(1, end);
    i = end + 1;
  } else {
    sugar = s[0] ?? '';
    i = 1;
  }

  // base in parens (optional)
  if (s[i] === '(') {
    const end = s.indexOf(')', i);
    base = cleanupHelmSymbol(s.substring(i + 1, end));
    i = end + 1;
  }

  // phosphate / linkage (optional, possibly bracketed)
  if (i < s.length) {
    if (s[i] === '[') {
      const end = s.indexOf(']', i);
      phosphate = s.substring(i + 1, end);
      i = end + 1;
    } else {
      phosphate = s.substring(i);
    }
  }

  // No base → either a standalone backbone linker (`p`, `[sp]`, …) or a
  // conjugate / cap (`[GalNAc]`, `[Chol]`, `[L3]`). Known phosphate symbols
  // become linkers (drawn as arcs); everything else is a conjugate.
  if (base === null) {
    if (isLinkerSymbol(sugar)) {
      const lnk: ParsedLinker = {kind: 'linker', position, symbol: sugar, raw};
      return lnk;
    }
    const conj: ParsedConjugate = {kind: 'conjugate', position, symbol: sugar, raw};
    return conj;
  }

  const nt: ParsedNucleotide = {kind: 'nucleotide', position, sugar, base, phosphate, raw};
  return nt;
}

/** Quick check that a string is HELM-like enough to attempt parsing. */
export function looksLikeHelm(s: string): boolean {
  return /^\s*(?:RNA|DNA|PEPTIDE|CHEM|BLOB)\d+\s*\{/.test(s ?? '');
}

/* ---------------------------------------------------------------- *
 * HELM connection + annotation parsing
 *
 * A full HELM string is `polymers$connections$groups$annotations$version`.
 * We read two of those extra sections:
 *   - connections (section 1): `RNA1,RNA2,2:pair-57:pair|...` — base pairs;
 *   - extended annotations (section 3): `{"RNA1":{"strandtype":"ss"}}`.
 * Monomer positions in connections are 1-based HELM monomer indices that count
 * each sub-monomer (sugar, base, linker) separately — so we map them back to
 * nucleotide ordinals to express pairs in nucleotide coordinates.
 * ---------------------------------------------------------------- */

const SENSE_TYPES = ['ss', 'sense'];
const ANTISENSE_TYPES = ['as', 'antisense'];

interface RawConnection {
  srcId: string; srcPos: number; srcAnn: string;
  tgtId: string; tgtPos: number; tgtAnn: string;
}

/** Parse the `$connections$` section into raw (id, position, annotation) tuples. */
export function parseHelmConnections(helm: string): RawConnection[] {
  const section = (helm ?? '').split('$')[1] ?? '';
  if (!section.trim()) return [];
  const out: RawConnection[] = [];
  for (const entry of section.split('|')) {
    const e = entry.trim();
    if (!e) continue;
    // `SrcID,TgtID,SrcPos:SrcAnn-TgtPos:TgtAnn`
    const parts = e.split(',');
    if (parts.length < 3) continue;
    const srcId = parts[0].trim();
    const tgtId = parts[1].trim();
    const dash = parts.slice(2).join(',').split('-');
    if (dash.length < 2) continue;
    const [srcPos, srcAnn] = splitPosAnnotation(dash[0]);
    const [tgtPos, tgtAnn] = splitPosAnnotation(dash[1]);
    if (srcPos == null || tgtPos == null) continue;
    out.push({srcId, srcPos, srcAnn, tgtId, tgtPos, tgtAnn});
  }
  return out;
}

function splitPosAnnotation(s: string): [number | null, string] {
  const [pos, ann] = s.split(':');
  const n = parseInt(pos, 10);
  return [Number.isFinite(n) ? n : null, (ann ?? '').trim()];
}

function isPairAnnotation(a: string): boolean { return /pair/i.test(a); }

/** Parse the extended-annotations section into `{polymerId: strandtype}`. */
export function parseHelmAnnotations(helm: string): Record<string, string> {
  const section = (helm ?? '').split('$')[3] ?? '';
  if (!section.trim()) return {};
  let json: unknown;
  try { json = JSON.parse(section); } catch { return {}; }
  if (!json || typeof json !== 'object') return {};
  const out: Record<string, string> = {};
  for (const [id, v] of Object.entries(json as Record<string, unknown>)) {
    const st = (v as {strandtype?: unknown})?.strandtype;
    if (typeof st === 'string') out[id] = st;
  }
  return out;
}

function findStrandByType(types: Record<string, string>, want: string[]): string | undefined {
  for (const [id, t] of Object.entries(types))
    if (want.includes(t.toLowerCase())) return id;
  return undefined;
}

function countNucleotides(strand: ParsedStrand): number {
  let n = 0;
  for (const m of strand.monomers) if (m.kind === 'nucleotide') n++;
  return n;
}

/** Map every 1-based HELM monomer position of a chain body to the 0-based
 * nucleotide ordinal that owns it. Non-nucleotide monomers (bare terminal
 * phosphates, conjugates) consume a position but map to nothing. */
function nucleotideOrdinalByHelmPos(content: string): Map<number, number> {
  const map = new Map<number, number>();
  let pos = 0; // running 1-based HELM monomer counter (incremented per component)
  let ord = 0; // nucleotide ordinal
  for (const seg of splitSegments(content)) {
    const m = parseMonomer(seg, 0);
    if (m.kind === 'nucleotide') {
      const nt = m as ParsedNucleotide;
      const count = 2 + (nt.phosphate ? 1 : 0); // sugar + base (+ linker)
      for (let k = 1; k <= count; k++) map.set(pos + k, ord);
      ord++;
      pos += count;
    } else {
      pos += 1; // bare phosphate / conjugate — single HELM monomer
    }
  }
  return map;
}

/** Convert raw `pair` connections to nucleotide-index pairs on sense/antisense. */
function connectionsToPairs(
  conns: RawConnection[], senseId: string | undefined, antiId: string | undefined,
  senseMap: Map<number, number>, antiMap: Map<number, number>,
): DuplexPair[] {
  const pairs: DuplexPair[] = [];
  for (const c of conns) {
    let senseOrd: number | undefined;
    let antiOrd: number | undefined;
    if (c.srcId === senseId && c.tgtId === antiId) {
      senseOrd = senseMap.get(c.srcPos); antiOrd = antiMap.get(c.tgtPos);
    } else if (c.srcId === antiId && c.tgtId === senseId) {
      senseOrd = senseMap.get(c.tgtPos); antiOrd = antiMap.get(c.srcPos);
    } else {
      continue;
    }
    if (senseOrd === undefined || antiOrd === undefined) continue;
    pairs.push({senseIdx: senseOrd, antiIdx: antiOrd});
  }
  return pairs;
}

/** Most-common shift implied by the pairs (robust to a stray mis-pair).
 * shift = senseIdx − (antiLen − 1 − antiIdx) for each pair. */
function shiftFromPairs(pairs: DuplexPair[], antiLen: number): number {
  const counts = new Map<number, number>();
  for (const p of pairs) {
    const s = p.senseIdx - (antiLen - 1 - p.antiIdx);
    counts.set(s, (counts.get(s) ?? 0) + 1);
  }
  let bestShift = 0;
  let bestN = -1;
  for (const [s, n] of counts.entries()) {
    if (n > bestN || (n === bestN && Math.abs(s) < Math.abs(bestShift))) {
      bestShift = s; bestN = n;
    }
  }
  return bestShift;
}

/* ---------------------------------------------------------------- *
 * Canonicalization — rewrite aliased sugar / phosphate symbols
 * (e.g. `mR` → `m`, `fR` → `fl2r`, `LR` → `lna`, `sP` → `sp`) to the
 * HELMCore canonical forms before sending the HELM to Bio's
 * monomer-library-driven pipelines (toAtomicLevelPanel, helmToMolfileV3K, …).
 * Without this, the central Bio library can't find aliased monomers.
 * ---------------------------------------------------------------- */

/** Serialize a parsed monomer back to HELM, bracketing multi-char symbols
 * and using canonical (HELMCore) sugar / phosphate codes. */
function serializeCanonicalMonomer(m: ParsedMonomer): string {
  if (m.kind === 'conjugate') return `[${m.symbol}]`;
  if (m.kind === 'linker') {
    const phos = canonicalPhosphateSymbol(m.symbol);
    return phos.length === 1 ? phos : `[${phos}]`;
  }
  const nt = m as ParsedNucleotide;
  const sugar = canonicalSugarSymbol(nt.sugar);
  const phos = canonicalPhosphateSymbol(nt.phosphate);
  const sugarPart = sugar.length === 1 ? sugar : `[${sugar}]`;
  const phosPart = !phos ? '' : (phos.length === 1 ? phos : `[${phos}]`);
  const baseStr = !nt.base ? '' : (nt.base.length === 1 ? `(${nt.base})` : `([${nt.base}])`);
  return `${sugarPart}${baseStr}${phosPart}`;
}

/** Canonicalize the body (between `{...}`) of a HELM chain. */
export function canonicalizeChainBody(body: string): string {
  return parseMonomers(body).map(serializeCanonicalMonomer).join('.');
}

/** Canonicalize a full HELM string (any number of `|`-separated chains). */
export function canonicalizeHelm(helm: string): string {
  const polymerSection = (helm ?? '').split('$')[0];
  const chains = polymerSection.split('|').map((c) => c.trim()).filter(Boolean);
  const out = chains.map((c) => {
    const m = c.match(/^((?:RNA|DNA|PEPTIDE|CHEM|BLOB)\d+)\{(.*)\}$/);
    if (!m) return c;
    const [, prefix, body] = m;
    return `${prefix}{${canonicalizeChainBody(body)}}`;
  });
  return `${out.join('|')}$$$$`;
}
