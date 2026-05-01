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

import {
  canonicalPhosphateSymbol, canonicalSugarSymbol,
  ParsedConjugate, ParsedDuplex, ParsedMonomer, ParsedNucleotide, ParsedStrand,
} from './types';

/** Parse the whole HELM string. Tolerant of whitespace and missing `$$$$`. */
export function parseHelmDuplex(helm: string): ParsedDuplex {
  const polymerSection = (helm ?? '').split('$')[0] ?? '';
  const chains = polymerSection.split('|').map((c) => c.trim()).filter((c) => c.length > 0);

  const parsedChains: ParsedStrand[] = chains.map(parseChain);
  return {
    sense: parsedChains[0] ?? {type: 'RNA', monomers: []},
    antisense: parsedChains[1] ?? null,
    raw: helm,
  };
}

function parseChain(chain: string): ParsedStrand {
  // RNA1{...}, DNA2{...}, CHEM3{...}
  const match = chain.match(/^(RNA|DNA|CHEM|PEPTIDE|BLOB)\d+\{(.*)\}$/);
  if (!match)
    return {type: 'CHEM', monomers: []};
  const type = match[1];
  const content = match[2];
  return {type, monomers: parseMonomers(content)};
}

/** Split chain contents on top-level dots (depth-aware to preserve `[a.b]`/`(c.d)`). */
function parseMonomers(content: string): ParsedMonomer[] {
  const out: ParsedMonomer[] = [];
  let depth = 0;
  let cur = '';
  let pos = 0;
  for (const ch of content) {
    if (ch === '[' || ch === '(') {
      depth++;
    } else if (ch === ']' || ch === ')') {
      depth--;
    } else if (ch === '.' && depth === 0) {
      if (cur) out.push(parseMonomer(cur, pos++));
      cur = '';
      continue;
    }
    cur += ch;
  }
  if (cur) out.push(parseMonomer(cur, pos));
  return out;
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
    base = s.substring(i + 1, end);
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

  // No base → it's a conjugate / linker / cap (e.g. [GalNAc], [Chol], [L3])
  if (base === null) {
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
  const nt = m as ParsedNucleotide;
  const sugar = canonicalSugarSymbol(nt.sugar);
  const phos = canonicalPhosphateSymbol(nt.phosphate);
  const sugarPart = sugar.length === 1 ? sugar : `[${sugar}]`;
  const phosPart = !phos ? '' : (phos.length === 1 ? phos : `[${phos}]`);
  const baseStr = nt.base ? `(${nt.base})` : '';
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
