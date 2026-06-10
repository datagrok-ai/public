/* eslint-disable max-params */
/**
 * Canvas drawing for the OligoNucleotide cell renderer.
 *
 * Pure functions — no DG / DOM dependencies — so the renderer is easy to
 * unit-test and to drive from the HTML prototype as well.
 *
 * Visual model:
 *   - Chip body — colored from the central Bio monomer library
 *     (`getMonomerColors(HELM_BASE, base)` → background + text). Falls back to
 *     local BASE_COLORS / analog colors if the lib has no entry.
 *   - Sugar modifications — narrow colored stripe glued to the *outside* edge
 *     of the chip row (top for sense / single-strand, bottom for antisense).
 *     Color from `getMonomerColors(HELM_SUGAR, sugar).backgroundcolor`.
 *   - Phosphate linkage modifications (PS / s2p / mp / …) — drawn as 45°
 *     apex triangles in the inter-chip gap. Apex points *outward* (up for
 *     sense / single-strand, down for antisense). Color from
 *     `getMonomerColors(HELM_LINKER, phos).backgroundcolor`.
 *   - Conjugates — rounded pills at chain ends; their actual width is
 *     propagated through layout so adjacent chips don't overlap. Color from
 *     `getMonomerColors(HELM_CHEM, symbol)`.
 *   - Antisense is rendered 3'→5' (reversed) so position N of sense visually
 *     pairs with position N of antisense (anti-parallel base-pair register).
 *
 * Sizing is fully adaptive: chip dimensions scale to fit the cell, preserving
 * aspect ratio, down to a minimum below which the cell falls back to a text
 * summary. The layout reserves vertical room for the apex zones above sense
 * (and below antisense, when present).
 */

import * as DG from 'datagrok-api/dg';

import {
  BASE_COLORS, DuplexAlignment, FALLBACK_COLOR,
  contrastTextColor,
  displayBase, isCanonicalBase,
  ParsedDuplex, ParsedMonomer, ParsedNucleotide, ParsedStrand,
  resolveConjugate, resolvePhosphate, resolveSugar,
} from './types';
import {getNaturalAnalog} from './analog-cache';
import {getMonomerColors} from './monomer-colors';
import {columnOffsets, reduceToCanonicalBase, resolveDuplexAlignment} from './alignment';

export interface RenderOpts {
  /** Show base letter inside chip. False at very small sizes. */
  showLetters: boolean;
  /** Draw antisense 3'→5' for base-pair alignment with sense. Default true. */
  pairAlign: boolean;
  /** Optional column-level scheme name (currently ignored — palette is fixed). */
  scheme?: string;
}

export const DEFAULT_OPTS: RenderOpts = {showLetters: true, pairAlign: true};

/* Visual tuning. */
const ASPECT_H_OVER_W = 1.25;
const STRAND_GAP_RATIO = 0.5;
const CHIP_GAP_RATIO = 0.40; // wide gaps give apex triangles room to breathe
const MIN_CHIP_W = 5;
// No hard upper cap on chip width — chipW grows until either (a) the duplex
// no longer fits in the cell's available width with all monomers shown at
// their natural widened size, or (b) the cell's height budget runs out.
// Show-everything is prioritised over making chips big.
const PAD = 4;
/** Vertical breathing room above the top apex zone and below the bottom one,
 * so the linkage marks never butt right up against the cell border. */
const V_PAD = 10;
const LABEL_W = 30;
/** Chip fill opacity — softens the often-saturated library backgrounds so the
 * sugar stripe and base label stay readable on top. */
const CHIP_FILL_ALPHA = 0.72;
const SUGAR_STRIPE_ALPHA = 1;
const CHIP_BORDER_W = 0.5;
const SUGAR_STRIPE_RATIO = 0.22; // sugar-mod stripe height as fraction of chipH
const APEX_RATIO = 0.45; // apex height as fraction of chipH (45° → also half-width)
const APEX_LINE_W = 2.5;

/** Cached layout for one rendered cell. */
export interface DuplexLayout {
  chipW: number;
  chipH: number;
  chipGap: number;
  strandGap: number;
  fontSize: number;
  labelW: number;
  padding: number;
  senseY: number;
  antiY: number; // -1 if no antisense
  seqX: number;
  /** Height of the apex zone above sense (and below antisense, if present). */
  apexH: number;
  textOnlyFallback: boolean;
  senseChips: ChipPos[];
  antiChips: ChipPos[];
  senseLinks: LinkagePos[];
  antiLinks: LinkagePos[];
  /** Whether antisense was rendered in reversed (3'→5') order. */
  antiReversed: boolean;
  /** Column offset of the antisense display relative to sense (see
   * `DuplexAlignment.shift`). 0 for a blunt / single-strand layout. */
  shift: number;
  /** Where the shift came from (explicit HELM pairs, auto-aligner, or none). */
  shiftSource: DuplexAlignment['source'];
}

/** The size / position fields of `DuplexLayout` available before chips and
 * linkages are placed (passed down to `placeStrand`). */
export type LayoutBase = Omit<DuplexLayout,
  'senseChips' | 'antiChips' | 'senseLinks' | 'antiLinks' | 'antiReversed' | 'shift' | 'shiftSource'>;

export interface ChipPos {
  x: number;
  w: number;
  monomer: ParsedMonomer;
  /** Original 0-based index in the data strand (regardless of display order). */
  origIdx: number;
  strand: StrandSide;
  /** Pairing column index (shared coordinate across strands). `-1` for
   * conjugates / non-nucleotides. Used to align base pairs across overhangs. */
  col: number;
}

export interface LinkagePos {
  x: number; w: number; y: number; h: number;
  phosphateSymbol: string;
  /** Original index of the nucleotide that owns the 3' phosphate. */
  ownerOrigIdx: number;
  strand: StrandSide;
}

export type StrandSide = 'sense' | 'antisense';

/* ---------------------------------------------------------------- *
 * Layout
 * ---------------------------------------------------------------- */

export function computeLayout(
  cellW: number, cellH: number, model: ParsedDuplex, opts: Partial<RenderOpts> = {},
): DuplexLayout {
  const o: RenderOpts = {...DEFAULT_OPTS, ...opts};
  const senseLen = model.sense.monomers.length;
  const antiLen = model.antisense ? model.antisense.monomers.length : 0;
  const hasAnti = antiLen > 0;
  const strandsCount = hasAnti ? 2 : 1;
  const maxLen = Math.max(1, senseLen, antiLen);

  // Strand alignment: explicit HELM pairs win; otherwise auto-align by
  // complementarity. The shift is the antisense display's column offset
  // relative to sense (overhangs ride on the column-offset difference). It is
  // only meaningful when antisense is rendered anti-parallel (pair-aligned).
  const alignment = resolveDuplexAlignment(model);
  const shift = (hasAnti && o.pairAlign) ? alignment.shift : 0;

  // Height budget caps chipW from above (chip aspect ratio is fixed). The
  // vertical pad reserves breathing room above/below the apex zones so the
  // linkage marks never touch the cell border.
  const apexCount = hasAnti ? 2 : 1;
  const heightFactor =
    strandsCount + (strandsCount - 1) * STRAND_GAP_RATIO + apexCount * APEX_RATIO;
  const hChipH = (cellH - 2 * V_PAD) / heightFactor;
  const heightCap = hChipH / ASPECT_H_OVER_W;

  // Pick chipW as the LARGEST value in [MIN_CHIP_W, heightCap] for which the
  // *uniform* layout (including any alignment shift / overhang columns) fits.
  // We size based on uniform so canonical single-char chips stay as big as the
  // cell allows; the body below then decides per cell whether the widened
  // multi-char labels also fit at this chipW. If even uniform doesn't fit at
  // MIN_CHIP_W, we fall back to a text-only summary.
  let chipW: number;
  let textOnlyFallback = false;
  // `maxLen` is read so this remains a no-op binding when callers reference it later.
  void maxLen;
  if (heightCap < MIN_CHIP_W) {
    chipW = MIN_CHIP_W;
    textOnlyFallback = true;
  } else if (computePlacement(MIN_CHIP_W, cellW, model, o, shift, 'uniform').fits) {
    chipW = findMaxChipW(cellW, model, o, shift, 'uniform', MIN_CHIP_W, heightCap);
  } else {
    chipW = MIN_CHIP_W;
    textOnlyFallback = true;
  }

  const chipH = chipW * ASPECT_H_OVER_W;
  const chipGap = Math.max(3, chipW * CHIP_GAP_RATIO);
  const strandGap = chipH * STRAND_GAP_RATIO;
  const fontSize = Math.max(7, Math.min(13, chipW * 0.62));
  const apexH = chipH * APEX_RATIO;

  const blockH = strandsCount * chipH + (strandsCount - 1) * strandGap + apexCount * apexH;
  const blockTop = Math.max(V_PAD, (cellH - blockH) / 2);
  // Sense always sits below a top apex zone; antisense (when present) has
  // its apex zone below it. Single-strand cases get a top apex zone only.
  const senseY = blockTop + apexH;
  const antiY = hasAnti ? senseY + chipH + strandGap : -1;
  const seqX = PAD + LABEL_W;

  const layoutBase: LayoutBase = {
    chipW, chipH, chipGap, strandGap, fontSize,
    labelW: LABEL_W, padding: PAD,
    senseY, antiY, seqX, apexH, textOnlyFallback,
  };

  if (textOnlyFallback) {
    return {
      ...layoutBase,
      senseChips: [], antiChips: [], senseLinks: [], antiLinks: [],
      antiReversed: false, shift: 0, shiftSource: alignment.source,
    };
  }

  const antiReversed = hasAnti && o.pairAlign;

  // Final placement. A non-zero shift forces uniform widths so the integer
  // column offsets land exactly on the chip grid (overhangs stay column-true);
  // a blunt duplex tries widened widths first and falls back to uniform when
  // the multi-char labels would overflow the cell.
  let mode: 'widened' | 'uniform' = shift !== 0 ? 'uniform' : 'widened';
  let place = computePlacement(chipW, cellW, model, o, shift, mode);
  if (mode === 'widened' && !place.fits) {
    mode = 'uniform';
    place = computePlacement(chipW, cellW, model, o, shift, 'uniform');
  }

  const senseRes = placeStrand(
    model.sense, false, senseY, place.senseStartX, place.seqEndX, layoutBase, 'sense',
    place.senseWidths, place.senseColOffset);
  const antiRes = hasAnti ? placeStrand(
    model.antisense!, antiReversed, antiY, place.antiStartX, place.seqEndX, layoutBase, 'antisense',
    place.antiWidths, place.antiColOffset,
  ) : {chips: [], links: []};

  return {
    ...layoutBase,
    senseChips: senseRes.chips,
    antiChips: antiRes.chips,
    senseLinks: senseRes.links,
    antiLinks: antiRes.links,
    antiReversed,
    shift,
    shiftSource: alignment.source,
  };
}

interface Placement {
  fits: boolean;
  senseStartX: number;
  antiStartX: number;
  seqEndX: number;
  senseWidths: number[];
  antiWidths: number[];
  senseColOffset: number;
  antiColOffset: number;
}

/** Compute strand start-X coordinates, per-chip widths, and whether everything
 * fits horizontally, for a candidate `chipW`, alignment `shift`, and width
 * `mode`. Shared by the chipW binary search and the final layout so sizing and
 * placement never diverge.
 *
 * Geometry: both strands' nucleotides sit on a shared column grid of pitch
 * `colW = chipW + chipGap`. Sense nucleotide `k` is at column `k + senseColOff`;
 * antisense display nucleotide `d` at `d + antiColOff`. The leftmost element
 * (conjugate or overhang nucleotide) is pinned to `seqX`. With `shift = 0` and
 * no offsets this reduces exactly to "align both strands' first nucleotides". */
function computePlacement(
  chipW: number, cellW: number, model: ParsedDuplex, opts: RenderOpts,
  shift: number, mode: 'widened' | 'uniform',
): Placement {
  const fontSize = Math.max(7, Math.min(13, chipW * 0.62));
  const chipGap = Math.max(3, chipW * CHIP_GAP_RATIO);
  const colW = chipW + chipGap;
  const hasAnti = !!model.antisense && model.antisense.monomers.length > 0;
  const antiReversed = hasAnti && opts.pairAlign;
  const {sense: senseColOffset, anti: antiColOffset} = columnOffsets(hasAnti ? shift : 0);

  const senseDisplay = model.sense.monomers;
  const antiDisplay = hasAnti ?
    (antiReversed ? model.antisense!.monomers.slice().reverse() : model.antisense!.monomers) :
    [];

  const senseLeadW = leadingConjugateWidth(model.sense.monomers, false, chipW, fontSize, chipGap);
  const antiLeadW = hasAnti ?
    leadingConjugateWidth(model.antisense!.monomers, antiReversed, chipW, fontSize, chipGap) : 0;

  const widths = mode === 'widened' ?
    computePairSyncedWidths(senseDisplay, antiDisplay, chipW, fontSize) :
    uniformWidths(senseDisplay, antiDisplay, chipW, fontSize);

  const seqX = PAD + LABEL_W;
  const seqEndX = cellW - PAD;
  // Pin the strand sticking out furthest left to seqX.
  const N0 = seqX + Math.max(senseLeadW - senseColOffset * colW, antiLeadW - antiColOffset * colW);
  const senseStartX = N0 + senseColOffset * colW - senseLeadW;
  const antiStartX = N0 + antiColOffset * colW - antiLeadW;

  const senseRightX = senseStartX + spanWidth(widths.sense, chipGap);
  const antiRightX = hasAnti ? antiStartX + spanWidth(widths.anti, chipGap) : senseStartX;
  const fits = seqEndX - seqX >= 0 && Math.max(senseRightX, antiRightX) <= seqEndX;

  return {
    fits, senseStartX, antiStartX, seqEndX,
    senseWidths: widths.sense, antiWidths: widths.anti,
    senseColOffset, antiColOffset,
  };
}

/** Total horizontal span of a per-chip widths array (chips + inter-chip gaps). */
function spanWidth(widths: number[], chipGap: number): number {
  if (widths.length === 0) return 0;
  let total = 0;
  for (const w of widths) total += w;
  return total + (widths.length - 1) * chipGap;
}

/** Place chips for one strand, optionally reversed. Returns chip and linkage positions.
 * `chipWidths` is per-display-monomer (i.e. matches the order of `strand.monomers` after
 * reversal if `reverse=true`). `colOffset` is the nucleotide column index of this strand's
 * first nucleotide on the shared duplex grid. */
function placeStrand(
  strand: ParsedStrand, reverse: boolean, y: number, startX: number, endX: number,
  layout: LayoutBase,
  side: StrandSide,
  chipWidths: number[],
  colOffset: number,
): { chips: ChipPos[]; links: LinkagePos[] } {
  const monomers = reverse ? strand.monomers.slice().reverse() : strand.monomers;
  const chips: ChipPos[] = [];
  const links: LinkagePos[] = [];
  let x = startX;
  let nucOrdinal = 0;

  for (let i = 0; i < monomers.length; i++) {
    const m = monomers[i];
    const w = chipWidths[i] ?? layout.chipW;

    if (x + w > endX + 0.5) break; // truncate at cell edge (0.5px float slack)

    const col = m.kind === 'nucleotide' ? colOffset + nucOrdinal : -1;
    if (m.kind === 'nucleotide') nucOrdinal++;
    chips.push({x, w, monomer: m, origIdx: m.position, strand: side, col});

    // Linkage in the gap to the right of display index i (between chips i
    // and i+1). The `phosphate` field on a nucleotide always means "what
    // comes AFTER this monomer in 5'→3' data order" — i.e. it lives on the
    // lower-indexed end of the bond. So:
    //   - not reversed: gap-to-right of display i pairs data (p, p+1), and
    //     the phosphate lives on `m` (data p);
    //   - reversed: gap-to-right of display i pairs data (p, p-1), and the
    //     phosphate lives on the lower-indexed end, which is monomers[i+1].
    if (i < monomers.length - 1) {
      const linkOwner = reverse ? monomers[i + 1] : m;
      if (linkOwner.kind === 'nucleotide') {
        const nt = linkOwner as ParsedNucleotide;
        if (nt.phosphate) {
          links.push({
            x: x + w, w: layout.chipGap, y, h: layout.chipH,
            phosphateSymbol: nt.phosphate,
            ownerOrigIdx: nt.position,
            strand: side,
          });
        }
      }
    }
    x += w + layout.chipGap;
  }
  return {chips, links};
}

/** Width estimate for a conjugate pill, using only chip metrics (no canvas measure). */
function estimateConjugateWidth(symbol: string, chipW: number, fontSize: number): number {
  const meta = resolveConjugate(symbol).meta;
  const charW = fontSize * 0.55;
  const textW = meta.short.length * charW;
  const padding = chipW * 0.6;
  return Math.max(chipW, Math.min(chipW * 4, textW + padding));
}

/** Sum of widths of leading conjugates in display order (after reversal if any),
 * including their trailing chipGap, so antisense and sense can be shifted
 * independently to align their first nucleotides at the same X coordinate. */
function leadingConjugateWidth(
  monomers: ParsedMonomer[], reversed: boolean,
  chipW: number, fontSize: number, chipGap: number,
): number {
  const seq = reversed ? monomers.slice().reverse() : monomers;
  let w = 0;
  for (const m of seq) {
    if (m.kind !== 'conjugate') break;
    w += estimateConjugateWidth(m.symbol, chipW, fontSize) + chipGap;
  }
  return w;
}

/** Desired width for one monomer if we render its base label in full (no
 * ellipsis). Conjugates get their pill width. Single/two-char bases just use
 * `chipW`; multi-char bases widen to fit the full label, capped at 3× chipW. */
function desiredChipWidth(m: ParsedMonomer, chipW: number, fontSize: number): number {
  if (m.kind === 'conjugate') return estimateConjugateWidth(m.symbol, chipW, fontSize);
  const base = (m as ParsedNucleotide).base ?? '';
  if (base.length <= 2) return chipW;
  // Empirical: glyph width ≈ fontSize * 0.55 for system-ui at our weights.
  const charW = fontSize * 0.55;
  const textW = base.length * charW;
  const padding = chipW * 0.4;
  return Math.max(chipW, Math.min(chipW * 3, textW + padding));
}

interface SyncedWidths { sense: number[]; anti: number[]; }

/** For each strand returns a per-chip width array; when both strands have a
 * chip at the same pair-index (counted past leading conjugates), the column
 * width is `max(senseDesired, antiDesired)` so pair-aligned positions stay
 * column-locked even when one side has a long base name. */
function computePairSyncedWidths(
  senseDisplay: ParsedMonomer[], antiDisplay: ParsedMonomer[], chipW: number, fontSize: number,
): SyncedWidths {
  const senseW = senseDisplay.map((m) => desiredChipWidth(m, chipW, fontSize));
  const antiW = antiDisplay.map((m) => desiredChipWidth(m, chipW, fontSize));

  // Strip leading conjugates: the first non-conjugate in display order.
  const senseStart = senseDisplay.findIndex((m) => m.kind === 'nucleotide');
  const antiStart = antiDisplay.findIndex((m) => m.kind === 'nucleotide');
  if (senseStart >= 0 && antiStart >= 0) {
    const pairLen = Math.min(senseDisplay.length - senseStart, antiDisplay.length - antiStart);
    for (let i = 0; i < pairLen; i++) {
      const si = senseStart + i;
      const ai = antiStart + i;
      const w = Math.max(senseW[si], antiW[ai]);
      senseW[si] = w; antiW[ai] = w;
    }
  }
  return {sense: senseW, anti: antiW};
}

/** Uniform-width fallback (every chip = chipW; conjugates still pill-wide). */
function uniformWidths(
  senseDisplay: ParsedMonomer[], antiDisplay: ParsedMonomer[], chipW: number, fontSize: number,
): SyncedWidths {
  const map = (m: ParsedMonomer) => m.kind === 'conjugate' ?
    estimateConjugateWidth(m.symbol, chipW, fontSize) : chipW;
  return {sense: senseDisplay.map(map), anti: antiDisplay.map(map)};
}

/** True iff a layout at this `chipW` (using either widened or uniform widths,
 * including any alignment shift) places every chip inside the cell. Thin
 * wrapper over `computePlacement` so the binary search and the final layout
 * agree exactly. */
function fitsAtChipW(
  chipW: number, cellW: number, model: ParsedDuplex, opts: RenderOpts,
  shift: number, mode: 'widened' | 'uniform',
): boolean {
  return computePlacement(chipW, cellW, model, opts, shift, mode).fits;
}

/** Binary-search the largest chipW ∈ [lo, hi] where `fitsAtChipW(..., mode)`
 * returns true. Assumes `fitsAtChipW(lo, ..., mode) === true` (caller checks).
 * Layout fit is monotone-decreasing in chipW, so a 0.25-px termination gives
 * sub-pixel resolution in a handful of iterations. */
function findMaxChipW(
  cellW: number, model: ParsedDuplex, opts: RenderOpts,
  shift: number, mode: 'widened' | 'uniform', lo: number, hi: number,
): number {
  if (fitsAtChipW(hi, cellW, model, opts, shift, mode)) return hi;
  while (hi - lo > 0.25) {
    const mid = (lo + hi) / 2;
    if (fitsAtChipW(mid, cellW, model, opts, shift, mode)) lo = mid;
    else hi = mid;
  }
  return lo;
}

/* ---------------------------------------------------------------- *
 * Drawing
 * ---------------------------------------------------------------- */

export function drawDuplex(
  g: CanvasRenderingContext2D, cellX: number, cellY: number,
  cellW: number, cellH: number, model: ParsedDuplex,
  opts: Partial<RenderOpts> = {}, skipDrawing = false,
): DuplexLayout {
  const o: RenderOpts = {...DEFAULT_OPTS, ...opts};
  const layout = computeLayout(cellW, cellH, model, o);
  if (skipDrawing) return layout;
  g.save();
  g.beginPath();
  g.rect(cellX, cellY, cellW, cellH);
  g.clip();
  g.translate(cellX, cellY);

  if (layout.textOnlyFallback) {
    drawFallbackText(g, cellW, cellH, model, layout);
    g.restore();
    return layout;
  }

  // Strand label "S 5'" left of sense
  drawStrandLabel(g, 'S', '5\'', layout.padding, layout.senseY + layout.chipH / 2, layout);
  // first draw links so they are behind
  for (const link of layout.senseLinks) drawLinkageApex(g, link, layout);
  // chip body's anti-aliased corners.
  drawChips(g, layout.senseChips, layout, o);
  drawTruncationMarker(g, layout.senseChips, model.sense.monomers.length, layout);

  if (layout.antiY >= 0 && model.antisense) {
    // When reversed, the leftmost chip in display is the 3' end of antisense.
    const leftLabel = layout.antiReversed ? '3\'' : '5\'';
    drawStrandLabel(g, 'AS', leftLabel, layout.padding, layout.antiY + layout.chipH / 2, layout);
    for (const link of layout.antiLinks) drawLinkageApex(g, link, layout);
    drawChips(g, layout.antiChips, layout, o);
    drawTruncationMarker(g, layout.antiChips, model.antisense.monomers.length, layout);

    // Watson-Crick base-pair indicators in the strand gap. Only meaningful
    // when antisense is rendered anti-parallel (i.e. reversed for display) —
    // that's what makes the vertical column alignment actually represent the
    // duplex partner pairs.
    if (layout.antiReversed) drawBasePairings(g, layout);
  }

  g.restore();
  return layout;
}

function drawStrandLabel(
  g: CanvasRenderingContext2D, strand: string, terminus: string, x: number, y: number, layout: DuplexLayout,
): void {
  g.fillStyle = '#8b949e';
  g.font = `${Math.max(8, layout.fontSize - 1)}px ui-monospace, Menlo, monospace`;
  g.textBaseline = 'middle';
  g.textAlign = 'left';
  g.fillText(`${strand} ${terminus}`, x, y);
}

function drawTruncationMarker(
  g: CanvasRenderingContext2D, chips: ChipPos[], totalCount: number, layout: DuplexLayout,
): void {
  if (chips.length >= totalCount) return;
  const last = chips[chips.length - 1];
  if (!last) return;
  g.fillStyle = '#6e7681';
  g.font = `${Math.max(7, layout.fontSize - 1)}px ui-monospace, Menlo, monospace`;
  g.textBaseline = 'middle';
  g.textAlign = 'left';
  g.fillText('…', last.x + last.w + 1, layout.senseY + layout.chipH / 2);
}

function drawChips(g: CanvasRenderingContext2D, chips: ChipPos[], layout: DuplexLayout, opts: RenderOpts): void {
  const side = chips[0]?.strand ?? 'sense';
  const y = side === 'sense' ? layout.senseY : layout.antiY;
  const decoSide = decorationSide(side);
  for (const cp of chips) {
    if (cp.monomer.kind === 'conjugate')
      drawConjugate(g, cp.monomer.symbol, cp.x, y, cp.w, layout.chipH, layout.fontSize);
    else
      drawChip(g, cp.monomer as ParsedNucleotide, cp.x, y, cp.w, layout.chipH, layout.fontSize, opts, decoSide);
  }
}

/** Which outside edge of a strand's chip row gets the sugar stripe and apex.
 * Sense (or single-strand) → 'top'; antisense → 'bottom'. */
function decorationSide(strand: StrandSide): 'top' | 'bottom' {
  return strand === 'antisense' ? 'bottom' : 'top';
}

function drawChip(
  g: CanvasRenderingContext2D, m: ParsedNucleotide, x: number, y: number,
  w: number, h: number, fontSize: number, opts: RenderOpts, decoSide: 'top' | 'bottom',
): void {
  const baseColors = m.base ? getMonomerColors('base', m.base) : null;
  const bg = baseColors?.backgroundcolor ?? resolveBaseColor(m.base);
  const textC = baseColors?.textcolor ?? contrastTextColor(bg);
  // Corner radius scales with the chip's smaller dimension so chips stay
  // pleasantly rounded at every size — the previous hardcoded 2.5px cap made
  // larger chips read as squares.
  const r = Math.min(w, h) / 4;
  const stripeH = Math.max(2, h * SUGAR_STRIPE_RATIO);

  // Chip body — background from monomer library (HELM_BASE), softened with
  // a light alpha so the sugar stripe and base label stay readable on top.
  drawRoundRect(g, x, y, w, h, r);
  g.fillStyle = withAlpha(bg, CHIP_FILL_ALPHA);
  g.fill();
  g.lineWidth = CHIP_BORDER_W;
  g.strokeStyle = 'rgba(0,0,0,0.22)';
  g.stroke();

  // Sugar stripe — drawn for every sugar (including canonical `r` / `d`) so
  // the sugar identity is always visually readable. Flips side per strand
  // (top for sense, bottom for antisense). Clipped to the rounded shape so
  // the stripe follows the chip's corners.
  const sugarColors = getMonomerColors('sugar', m.sugar);
  const stripeColor =
    sugarColors.backgroundcolor ?? resolveSugar(m.sugar, m.base).color;
  g.save();
  drawRoundRect(g, x, y, w, h, r);
  g.clip();
  g.fillStyle = withAlpha(stripeColor, SUGAR_STRIPE_ALPHA);
  if (decoSide === 'top')
    g.fillRect(x, y, w, stripeH);
  else
    g.fillRect(x, y + h - stripeH, w, stripeH);
  g.restore();

  // Base label — biased AWAY from the stripe edge so it stays centered in
  // the visible body. Full HELM symbol when the chip is wide enough; first-
  // letter + ellipsis otherwise.
  if (opts.showLetters && m.base && fontSize >= 8) {
    const label = pickBaseLabel(m.base, w, fontSize);
    g.fillStyle = textC;
    g.font = `600 ${fontSize}px system-ui, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif`;
    g.textBaseline = 'middle';
    g.textAlign = 'center';
    // Shift the label towards the unstriped half of the chip
    const yShift = decoSide === 'top' ? stripeH / 2 : -stripeH / 2;
    g.fillText(label, x + w / 2, y + h / 2 + yShift + 0.5);
  }
}

/** Fallback base-color resolution when the central monomer library has no
 * `backgroundcolor` for the base. Tries the canonical palette, then the
 * natural analog's palette, then a neutral fallback. */
function resolveBaseColor(base: string | null): string {
  if (!base) return FALLBACK_COLOR;
  if (BASE_COLORS[base]) return BASE_COLORS[base];
  if (isCanonicalBase(base)) return BASE_COLORS[base] ?? FALLBACK_COLOR;
  const analog = getNaturalAnalog(base);
  if (analog && BASE_COLORS[analog]) return BASE_COLORS[analog];
  return FALLBACK_COLOR;
}

/** Decide the on-chip label given the chip's actual width: full base symbol
 * if it fits, else `firstLetter + …`. Single/two-char bases always show fully. */
function pickBaseLabel(base: string, chipW: number, fontSize: number): string {
  if (base.length <= 2) return base;
  const charW = fontSize * 0.55;
  const fullW = base.length * charW + 4; // small horizontal padding
  if (fullW <= chipW) return base;
  return displayBase(base);
}

/** Draw the linkage marker for `link` — a soft arch anchored to the strand's
 * outside edge (top for sense, bottom for antisense). A single quadratic
 * curve sweeps from base-left to base-right with a rounded summit, matching
 * the rounded chip style. Color comes from the central Bio monomer library's
 * `backgroundcolor` (linecolor tends to be flat black across the lib, which
 * would wash differentiation out). Canonical `p` is drawn too — every
 * linkage is visually accounted for. */
function drawLinkageApex(g: CanvasRenderingContext2D, link: LinkagePos, layout: DuplexLayout): void {
  const linkerColors = getMonomerColors('linker', link.phosphateSymbol);
  const color = linkerColors.backgroundcolor ?? resolvePhosphate(link.phosphateSymbol).color;
  const apexH = layout.apexH;

  const apexWidthMult = Math.max(apexH / 10, 1);
  const halfW = apexH; // 45° → height equals half-base
  const centerX = link.x + link.w / 2;
  const decoSide = decorationSide(link.strand);
  const baseY = decoSide === 'top' ? link.y : link.y + link.h;
  // Quadratic-curve midpoint Y sits halfway between baseY and the control Y,
  // so overshooting the control by 2× lands the visual peak at apexH from baseY.
  const ctrlY = decoSide === 'top' ? baseY - apexH * 2 : baseY + apexH * 2;

  g.save();
  g.beginPath();
  g.moveTo(centerX - halfW, baseY);
  g.quadraticCurveTo(centerX, ctrlY, centerX + halfW, baseY);
  g.lineWidth = APEX_LINE_W * apexWidthMult;
  g.lineCap = 'round';
  g.strokeStyle = color;
  g.stroke();
  g.restore();
}

/* ---------------------------------------------------------------- *
 * Watson-Crick base pairing
 *
 * For each display column where a sense chip sits above an antisense chip
 * (counted past leading conjugates), determine the pair kind and draw the
 * canonical biology shorthand in the inter-strand gap:
 *   - G ↔ C   → 3 vertical lines (3 hydrogen bonds — stronger color)
 *   - A ↔ U / A ↔ T → 2 vertical lines (2 hydrogen bonds — lighter color)
 *   - anything else → dashed line (mismatch / bulge)
 *
 * Non-canonical bases (modified analogs like `5meC`, `psiU`, `cpm6A`) are
 * resolved via `getNaturalAnalog` against the central Bio monomer library;
 * pairing is determined on the natural-analog letter, so a `2'-OMe-5meC`
 * still pairs as `C` with a `G` partner.
 * ---------------------------------------------------------------- */
type PairKind = 'GC' | 'AU' | 'mismatch';
const PAIR_COLOR_AU = '#4a5da8'; // 3 H-bonds — bolder indigo
const PAIR_COLOR_GC = '#8c9dc8'; // 2 H-bonds — lighter blue-violet
const PAIR_COLOR_MISMATCH = '#a5a5a5'; // dashed neutral gray
const PAIR_LINE_W = 1.1;
const MISMATCH_LINE_W = 1;

function drawBasePairings(g: CanvasRenderingContext2D, layout: DuplexLayout): void {
  const yTop = layout.senseY + layout.chipH;
  const yBot = layout.antiY;
  if (yBot <= yTop) return;

  // Pair by shared column index, so overhang columns (present on only one
  // strand) simply draw no mark, and a shift lines partners up correctly.
  const antiByCol = new Map<number, ChipPos>();
  for (const c of layout.antiChips)
    if (c.monomer.kind === 'nucleotide' && c.col >= 0) antiByCol.set(c.col, c);

  for (const sc of layout.senseChips) {
    if (sc.monomer.kind !== 'nucleotide' || sc.col < 0) continue;
    const ac = antiByCol.get(sc.col);
    if (!ac) continue;
    const kind = basePairKind(
      (sc.monomer as ParsedNucleotide).base,
      (ac.monomer as ParsedNucleotide).base,
    );
    // Pair markers are anchored on the SENSE chip's center because both
    // chips share the same X column when pair-aligned, but multi-char bases
    // can give them slightly different widths.
    const x = sc.x + sc.w / 2;
    drawPairingMark(g, x, yTop, yBot, kind, sc.w);
  }
}

function drawPairingMark(
  g: CanvasRenderingContext2D, x: number, yTop: number, yBot: number,
  kind: PairKind, chipW: number,
): void {
  g.save();
  g.lineCap = 'round';
  const lineWidthMult = Math.max(chipW / 20, 1);
  if (kind === 'mismatch') {
    g.strokeStyle = PAIR_COLOR_MISMATCH;
    g.lineWidth = MISMATCH_LINE_W * lineWidthMult;
    g.setLineDash([1, 1.6]);
    g.beginPath();
    g.moveTo(x, yTop);
    g.lineTo(x, yBot);
    g.stroke();
  } else if (kind === 'GC') {
    // Three lines — spacing scales with chip width to stay readable on
    // both narrow and wide chips, capped so they don't bleed past the chip.
    const off = chipW * 0.22;
    g.strokeStyle = PAIR_COLOR_GC;
    g.lineWidth = PAIR_LINE_W * lineWidthMult;
    for (const dx of [-off, 0, off]) {
      g.beginPath();
      g.moveTo(x + dx, yTop + 1);
      g.lineTo(x + dx, yBot - 1);
      g.stroke();
    }
  } else { // AU / AT
    const off = chipW * 0.14;
    g.strokeStyle = PAIR_COLOR_AU;
    g.lineWidth = PAIR_LINE_W * lineWidthMult;
    for (const dx of [-off, off]) {
      g.beginPath();
      g.moveTo(x + dx, yTop + 1);
      g.lineTo(x + dx, yBot - 1);
      g.stroke();
    }
  }
  g.restore();
}

function basePairKind(senseBase: string | null, antiBase: string | null): PairKind {
  const s = reduceToCanonicalBase(senseBase);
  const a = reduceToCanonicalBase(antiBase);
  if (!s || !a) return 'mismatch';
  if ((s === 'G' && a === 'C') || (s === 'C' && a === 'G')) return 'GC';
  if ((s === 'A' && (a === 'U' || a === 'T')) ||
      ((s === 'U' || s === 'T') && a === 'A')) return 'AU';
  return 'mismatch';
}

function drawConjugate(
  g: CanvasRenderingContext2D, symbol: string, x: number, y: number,
  w: number, chipH: number, fontSize: number,
): void {
  const conjColors = getMonomerColors('chem', symbol);
  const conj = resolveConjugate(symbol);
  const fill = conjColors.backgroundcolor ?? conj.color;
  const textC = conjColors.textcolor ?? '#ffffff';
  const r = chipH / 2;
  drawRoundRect(g, x, y, w, chipH, r);
  g.fillStyle = fill;
  g.fill();
  g.lineWidth = 0.5;
  g.strokeStyle = 'rgba(0,0,0,0.2)';
  g.stroke();

  if (fontSize >= 8) {
    g.fillStyle = textC;
    g.font = `600 ${Math.max(8, fontSize - 1)}px system-ui, sans-serif`;
    g.textBaseline = 'middle';
    g.textAlign = 'center';
    g.fillText(conj.meta.short, x + w / 2, y + chipH / 2 + 0.5);
  }
}

function drawRoundRect(
  g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, r: number,
): void {
  g.beginPath();
  g.moveTo(x + r, y);
  g.lineTo(x + w - r, y);
  g.quadraticCurveTo(x + w, y, x + w, y + r);
  g.lineTo(x + w, y + h - r);
  g.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
  g.lineTo(x + r, y + h);
  g.quadraticCurveTo(x, y + h, x, y + h - r);
  g.lineTo(x, y + r);
  g.quadraticCurveTo(x, y, x + r, y);
  g.closePath();
}

function drawFallbackText(
  g: CanvasRenderingContext2D, _w: number, h: number, model: ParsedDuplex, layout: DuplexLayout,
): void {
  const sLen = model.sense.monomers.length;
  const aLen = model.antisense ? model.antisense.monomers.length : 0;
  const summary = model.antisense ? `${sLen}+${aLen} nt duplex` : `${sLen} nt`;
  g.fillStyle = '#8b949e';
  g.font = `${Math.max(8, layout.fontSize)}px ui-monospace, Menlo, monospace`;
  g.textBaseline = 'middle';
  g.textAlign = 'left';
  g.fillText(summary, 4, h / 2);
}

/** Apply alpha to any CSS color string, returning rgba(...). Memoized. */
const _alphaCache = new DG.LruCache<string, string>(256);
function withAlpha(color: string, alpha: number): string {
  const key = `${color}|${alpha}`;
  const cached = _alphaCache.get(key);
  if (cached) return cached;
  const probe = document.createElement('span');
  probe.style.color = color;
  document.body.appendChild(probe);
  const rgb = getComputedStyle(probe).color;
  document.body.removeChild(probe);
  const m = rgb.match(/\d+/g);
  const out = m ? `rgba(${m[0]},${m[1]},${m[2]},${alpha})` : color;
  _alphaCache.set(key, out);
  return out;
}

/* ---------------------------------------------------------------- *
 * Hit testing — uses the cached chip / linkage positions so it works
 * correctly with variable widths (conjugate pills) and reversed antisense.
 * ---------------------------------------------------------------- */

export interface HitResult {
  strand: StrandSide;
  /** Original 0-based index in the data strand. */
  position: number;
  monomer: ParsedMonomer;
  /** Set if the hit is on an inter-chip linkage marker (PS, etc.). */
  linkage?: { phosphateSymbol: string };
}

export function hitTest(
  localX: number, localY: number, model: ParsedDuplex, layout: DuplexLayout,
): HitResult | null {
  if (layout.textOnlyFallback) return null;

  const apexH = layout.apexH;
  const chipH = layout.chipH;

  // Sense band: chip row + top apex zone (apex sits above the chip row).
  if (localY >= layout.senseY - apexH && localY <= layout.senseY + chipH) {
    if (localY >= layout.senseY) {
      const cp = findChip(localX, layout.senseChips);
      if (cp) return {strand: 'sense', position: cp.origIdx, monomer: cp.monomer};
    }
    // Apex zone above sense
    if (localY < layout.senseY) {
      const link = findApex(localX, localY, layout.senseLinks, 'top', apexH);
      if (link) return resolveLinkHit(link, model.sense, 'sense');
    }
  }
  // Antisense band: chip row + bottom apex zone (apex sits below the chip row).
  if (layout.antiY >= 0 && localY >= layout.antiY && localY <= layout.antiY + chipH + apexH) {
    if (localY <= layout.antiY + chipH) {
      const cp = findChip(localX, layout.antiChips);
      if (cp) return {strand: 'antisense', position: cp.origIdx, monomer: cp.monomer};
    }
    if (localY > layout.antiY + chipH && model.antisense) {
      const link = findApex(localX, localY, layout.antiLinks, 'bottom', apexH);
      if (link) return resolveLinkHit(link, model.antisense, 'antisense');
    }
  }
  return null;
}

function findChip(x: number, chips: ChipPos[]): ChipPos | null {
  for (const cp of chips) if (x >= cp.x && x < cp.x + cp.w) return cp;
  return null;
}

/** Find a linkage whose apex triangle covers (x, y). For each link, the apex
 * is centered on `link.x + link.w/2`, with 45° slopes — total base width is
 * `2 * apexH`. The bounding box is used (slight over-inclusion at corners is
 * fine for hover targets and matches what users expect). */
function findApex(
  x: number, y: number, links: LinkagePos[], side: 'top' | 'bottom', apexH: number,
): LinkagePos | null {
  for (const l of links) {
    const centerX = l.x + l.w / 2;
    if (x < centerX - apexH || x > centerX + apexH) continue;
    if (side === 'top') {
      if (y >= l.y - apexH && y <= l.y) return l;
    } else {
      if (y >= l.y + l.h && y <= l.y + l.h + apexH) return l;
    }
  }
  return null;
}

function resolveLinkHit(link: LinkagePos, strand: ParsedStrand, side: StrandSide): HitResult {
  const owner = strand.monomers.find((m) => m.position === link.ownerOrigIdx)!;
  return {
    strand: side,
    position: link.ownerOrigIdx,
    monomer: owner,
    linkage: {phosphateSymbol: link.phosphateSymbol},
  };
}
