/**
 * Types and modification dictionary for the OligoNucleotide cell renderer.
 *
 * The renderer is intentionally self-contained — it does NOT depend on
 * `_package.initLibData()` having completed, so the grid can render oligo
 * cells before any app has opened. Only the tooltip's RDKit structure
 * depends on the monomer library (which is loaded async there).
 */

export const OLIGO_SEM_TYPE = 'OligoNucleotide';
export const OLIGO_UNITS = 'helm';

/** Per-nucleotide parsed model. One per chip in the rendered output. */
export interface ParsedNucleotide {
  kind: 'nucleotide';
  /** Index inside the strand (0-based). */
  position: number;
  /** Sugar HELM symbol, normalized (no brackets). HELMCore canonical: 'r', 'd', 'm', 'fl2r', 'lna'. */
  sugar: string;
  /** Single-letter base — A / C / G / U / T or null when missing. */
  base: string | null;
  /** Phosphate HELM symbol normalized (no brackets). HELMCore canonical: 'p', 'sp'. Empty for terminal. */
  phosphate: string;
  /** Original HELM monomer text for tooltip + structure lookup. */
  raw: string;
}

/** Terminal conjugate (GalNAc, cholesterol, biotin, L3 linker, …). */
export interface ParsedConjugate {
  kind: 'conjugate';
  position: number;
  symbol: string;
  raw: string;
}

/** A standalone backbone linker monomer (`p`, `[sp]`, …) — a phosphate-type
 * unit with no base of its own. Rendered as an arc (no chip), like the trailing
 * phosphate of a nucleotide, and counted as a linkage (not a conjugate). */
export interface ParsedLinker {
  kind: 'linker';
  position: number;
  /** Phosphate HELM symbol, normalized (no brackets) — e.g. `p`, `sp`. */
  symbol: string;
  raw: string;
}

export type ParsedMonomer = ParsedNucleotide | ParsedConjugate | ParsedLinker;

export interface ParsedStrand {
  /** 'RNA' | 'DNA' | 'CHEM' */
  type: string;
  monomers: ParsedMonomer[];
  /** HELM polymer id this strand came from (e.g. `RNA1`). Undefined for
   * synthetic strands. Used to map HELM `pair` connections back to a strand. */
  id?: string;
}

/** A base-pair connection parsed from the HELM connection section, expressed
 * in *nucleotide* indices (0-based, 5'→3') on the sense and antisense strands.
 * Both indices count nucleotides only (conjugates / bare phosphates skipped). */
export interface DuplexPair {
  /** Nucleotide index on the sense strand (5'→3'). */
  senseIdx: number;
  /** Nucleotide index on the antisense strand (5'→3'). */
  antiIdx: number;
}

/** How the sense / antisense strands line up for display.
 *
 * `shift` is the column offset of the antisense *display* (rendered 3'→5')
 * relative to the sense strand. With antisense reversed for display, sense
 * nucleotide `j` pairs antisense nucleotide `k` when
 *   `j === (antiLen - 1 - k) + shift`.
 * `shift === 0` is a blunt duplex with both 5' ends flush. Positive shift
 * means a 5' sense overhang on the left; negative means a 3' antisense
 * overhang on the left (and the mirror overhangs on the right). */
export interface DuplexAlignment {
  shift: number;
  /** Where the shift came from: explicit HELM `pair` connections, the
   * complementary-base auto-aligner, or trivially none (single strand). */
  source: 'explicit' | 'auto' | 'none';
}

export interface ParsedDuplex {
  sense: ParsedStrand;
  antisense: ParsedStrand | null;
  /** Original HELM string (kept verbatim for tooltip / debug). */
  raw: string;
  /** Base-pair connections parsed from the HELM `$connections$` section, in
   * nucleotide indices on sense/antisense. Empty / absent when the HELM has
   * no explicit pairing info. */
  pairs?: DuplexPair[];
  /** Per-polymer strand-type annotations from the HELM extended-annotations
   * section, keyed by polymer id (e.g. `{RNA1: 'ss', RNA2: 'as'}`). */
  strandTypes?: Record<string, string>;
  /** Strand alignment derived from explicit pairs when present. Absent when
   * no explicit info exists — the renderer then auto-aligns by complementarity. */
  alignment?: DuplexAlignment;
}

/** Visual classification used for chip coloring. */
export type ModCategory = 'sugar' | 'phosphate' | 'conjugate' | 'base';

export interface ModMeta {
  /** Full user-facing name, used in tooltip + legend. */
  name: string;
  /** 3–8 char display label for chips that don't show a base letter (conjugates). */
  short: string;
  /** Canonical CSS color. */
  color: string;
  category: ModCategory;
}

/* ---------------------------------------------------------------- *
 * Modification dictionary — flat lookups by HELM symbol.
 *
 * Canonical symbols are HELMCore (Pistoia / Bio package):
 *   sugars      — `r`, `d`, `m`, `25r`, `fl2r`, `lna`, `moe`, ...
 *   phosphates  — `p`, `sp`, ...
 *   bases       — `A`, `C`, `G`, `U`, `T`
 *
 * Some vendor / Axolabs-style strings (`R`, `mR`, `fR`, `LR`, `dR`, `sP`)
 * appear in legacy data and in older internal tools. The SYMBOL_ALIASES
 * map normalizes them to canonical at lookup time, so the renderer
 * accepts both. Tooltips / structure look-up use the canonical symbol
 * (which is what the central Bio monomer library actually contains).
 *
 * Anything not in the dictionary falls back to a deterministic hash color.
 * ---------------------------------------------------------------- */

export const SUGAR_MODS: Readonly<Record<string, ModMeta>> = Object.freeze({
  // Canonical sugars — chip uses base-canonical color (see BASE_COLORS), not a sugar color
  'r': {name: 'Ribose', short: 'RNA', color: '#D0D0D0', category: 'sugar'},
  'd': {name: 'Deoxyribose', short: 'DNA', color: '#BDBDBD', category: 'sugar'},

  // 2'-O-Methyl
  'm': {name: '2\'-O-Methyl', short: '2\'-OMe', color: '#4A90E2', category: 'sugar'},
  '25r': {name: '2\'-O-Methyl (2,5)', short: '2\'-OMe', color: '#4A90E2', category: 'sugar'},

  // 2'-Fluoro
  'fl2r': {name: '2\'-Fluoro', short: '2\'-F', color: '#50C878', category: 'sugar'},

  // LNA / cEt / ENA — bridged
  'lna': {name: 'LNA', short: 'LNA', color: '#E89D3C', category: 'sugar'},
  'cet': {name: 'cEt', short: 'cEt', color: '#E07B30', category: 'sugar'},
  'ena': {name: 'ENA', short: 'ENA', color: '#D26C20', category: 'sugar'},

  // 2'-MOE
  'moe': {name: '2\'-MOE', short: '2\'-MOE', color: '#66CDAA', category: 'sugar'},

  // GNA / UNA / FANA / hexitol-NA
  'Rgna': {name: '(R)-GNA', short: 'GNA', color: '#C7A981', category: 'sugar'},
  'Sgna': {name: '(S)-GNA', short: 'GNA', color: '#C7A981', category: 'sugar'},
  'una': {name: 'UNA', short: 'UNA', color: '#ED68B8', category: 'sugar'},
  'fana': {name: 'FANA', short: 'FANA', color: '#76B947', category: 'sugar'},
  'hna': {name: 'HNA', short: 'HNA', color: '#A4D080', category: 'sugar'},
});

/** Vendor / legacy / Pistoia-style sugar symbols mapped to canonical HELMCore. */
export const SUGAR_ALIASES: Readonly<Record<string, string>> = Object.freeze({
  'R': 'r',
  'dR': 'd',
  'mR': 'm',
  'fR': 'fl2r',
  'LR': 'lna',
  'L': 'lna',
  'MOE': 'moe',
  'GNA': 'Sgna',
  'UNA': 'una',
});

export const PHOSPHATE_MODS: Readonly<Record<string, ModMeta>> = Object.freeze({
  'p': {name: 'Phosphate', short: 'PO', color: '#888888', category: 'phosphate'},
  'sp': {name: 'Phosphorothioate', short: 'PS', color: '#6B46C1', category: 'phosphate'},
  's2p': {name: 'Phosphorodithioate', short: 'PS₂', color: '#7E2BC4', category: 'phosphate'},
  'mp': {name: 'Methylphosphonate', short: 'MeP', color: '#9B5DE5', category: 'phosphate'},
});

export const PHOSPHATE_ALIASES: Readonly<Record<string, string>> = Object.freeze({
  'P': 'p',
  'sP': 'sp',
});

export const CONJUGATE_MODS: Readonly<Record<string, ModMeta>> = Object.freeze({
  'GalNAc': {name: 'GalNAc', short: 'GalNAc', color: '#E91E63', category: 'conjugate'},
  'L3': {name: 'GalNAc-L3 linker', short: 'L3', color: '#E91E63', category: 'conjugate'},
  'Chol': {name: 'Cholesterol', short: 'Chol', color: '#B57F50', category: 'conjugate'},
  'Bio': {name: 'Biotin', short: 'Bio', color: '#FFC107', category: 'conjugate'},
});

/** Pale per-base colors used when sugar is unmodified (R / dR).
 * Modifications override these. */
export const BASE_COLORS: Readonly<Record<string, string>> = Object.freeze({
  A: '#C1E7C5', // pale green
  C: '#C5DCF0', // pale blue
  G: '#E8DCA9', // pale tan
  U: '#F0C5C5', // pale pink
  T: '#F0C5C5',
});

/** The canonical single-letter base symbols. Anything else is a custom /
 * modified base (e.g. `cpm6A`, `5BrU`, `psiU`) and gets its color via natural
 * analog lookup against the central Bio monomer library. */
export const CANONICAL_BASES: Readonly<Set<string>> = new Set(['A', 'C', 'G', 'U', 'T']);

export const FALLBACK_COLOR = '#BCBCBC';

/** True if `base` is one of the canonical single-letter symbols. */
export function isCanonicalBase(base: string | null | undefined): boolean {
  return !!base && CANONICAL_BASES.has(base);
}

/** Shorten a multi-character base symbol to a chip-friendly label.
 *  - 1-2 char symbols: returned as-is
 *  - 3+ char symbols: first letter + ellipsis (e.g. `cpm6A` → `c…`)
 * The renderer also has a width-aware path that shows the full symbol when
 * the chip can fit it; this helper is the safe fallback. */
export function displayBase(base: string | null): string {
  if (!base) return '';
  if (base.length <= 2) return base;
  return base[0] + '…';
}

/** Deterministic color for unknown HELM monomer symbols. Stable across cells. */
export function hashColor(symbol: string): string {
  let h = 0;
  for (let i = 0; i < symbol.length; i++)
    h = ((h * 31) + symbol.charCodeAt(i)) | 0;
  const hue = ((h >>> 0) % 360);
  return `hsl(${hue}, 42%, 60%)`;
}

/** Canonicalize a sugar symbol via the alias map, then return the canonical form. */
export function canonicalSugarSymbol(sugar: string): string {
  return SUGAR_ALIASES[sugar] ?? sugar;
}

/** Canonicalize a phosphate symbol via the alias map. */
export function canonicalPhosphateSymbol(phosphate: string): string {
  return PHOSPHATE_ALIASES[phosphate] ?? phosphate;
}

/** True if `symbol` is a known backbone linker / phosphate (`p`, `sp`, `s2p`,
 * `mp`, or a legacy alias like `P`, `sP`). Static — no monomer-library lookup.
 * Library-driven detection (custom symbols whose natural analog is `p`) is
 * layered on top of this in `alignment.isLinkerMonomer`. */
export function isLinkerSymbol(symbol: string): boolean {
  if (!symbol) return false;
  return canonicalPhosphateSymbol(symbol) in PHOSPHATE_MODS;
}

/** Resolve a sugar HELM symbol to (color + meta). Accepts canonical or aliased input. */
export function resolveSugar(sugar: string, base: string | null): { color: string; meta: ModMeta } {
  const canonical = canonicalSugarSymbol(sugar);
  const known = SUGAR_MODS[canonical];
  if (known) {
    // For ribose / deoxyribose use the per-base canonical color rather than a flat gray
    if ((canonical === 'r' || canonical === 'd') && base && BASE_COLORS[base])
      return {color: BASE_COLORS[base], meta: known};
    return {color: known.color, meta: known};
  }
  const c = hashColor(canonical);
  return {color: c, meta: {name: canonical, short: canonical, color: c, category: 'sugar'}};
}

export function resolvePhosphate(phosphate: string): { color: string; meta: ModMeta } {
  if (!phosphate)
    return {color: '#888', meta: {name: '(none)', short: '', color: '#888', category: 'phosphate'}};
  const canonical = canonicalPhosphateSymbol(phosphate);
  const known = PHOSPHATE_MODS[canonical];
  if (known)
    return {color: known.color, meta: known};
  const c = hashColor(canonical);
  return {color: c, meta: {name: canonical, short: canonical, color: c, category: 'phosphate'}};
}

export function resolveConjugate(symbol: string): { color: string; meta: ModMeta } {
  const known = CONJUGATE_MODS[symbol];
  if (known)
    return {color: known.color, meta: known};
  return {
    color: hashColor(symbol),
    meta: {name: symbol, short: symbol.length > 6 ? symbol.slice(0, 6) : symbol,
      color: hashColor(symbol), category: 'conjugate'},
  };
}

/** Pick black or white text for readability on a given background hex/rgb/hsl color. */
export function contrastTextColor(bg: string): string {
  // Resolve any CSS color string to RGB through a temp DOM element
  const probe = document.createElement('span');
  probe.style.color = bg;
  document.body.appendChild(probe);
  const computed = getComputedStyle(probe).color;
  document.body.removeChild(probe);
  const m = computed.match(/\d+/g);
  if (!m) return '#000';
  const [r, g, b] = m.slice(0, 3).map(Number);
  // perceptual luminance — black on light, white on dark
  const lum = 0.299 * r + 0.587 * g + 0.114 * b;
  return lum > 160 ? '#1a1a1a' : '#ffffff';
}
