/**
 * Canvas drawing for the OligoNucleotide cell renderer.
 *
 * Pure functions — no DG / DOM dependencies — so the renderer is easy to
 * unit-test and to drive from the HTML prototype as well.
 *
 * Visual model:
 *   - Chip body uses base-canonical color (A green / C blue / G tan / U pink).
 *   - Sugar modifications are shown as a colored stripe at the bottom of the
 *     chip (so the chip itself stays readable for the base sequence, and the
 *     modification track scans easily horizontally).
 *   - Phosphate linkage modifications (PS) are shown as a saturated bar in
 *     the gap *between* chips — this is the actual chemistry: the linkage
 *     belongs to the bond, not to either nucleotide.
 *   - Wide gaps between chips give the modification markers room to breathe.
 *   - Antisense is rendered 3'→5' (reversed) so position N of sense visually
 *     pairs with position N of antisense (anti-parallel base-pair register).
 *   - Conjugates are rendered as wider pills at the chain ends; their actual
 *     width is propagated through layout so adjacent chips don't overlap.
 *
 * Sizing is fully adaptive: chip dimensions scale to fit the cell, preserving
 * aspect ratio, down to a minimum below which the cell falls back to a text
 * summary.
 */

import {
  BASE_COLORS, FALLBACK_COLOR,
  canonicalPhosphateSymbol, canonicalSugarSymbol,
  displayBase, isCanonicalBase,
  ParsedDuplex, ParsedMonomer, ParsedNucleotide, ParsedStrand,
  resolveConjugate, resolvePhosphate, resolveSugar,
} from './types';
import {getNaturalAnalog} from './analog-cache';

export interface RenderOpts {
  /** Show base letter inside chip. False at very small sizes. */
  showLetters: boolean;
  /** Draw antisense 3'→5' for base-pair alignment with sense. Default true. */
  pairAlign: boolean;
  /** Optional column-level scheme name (currently ignored — palette is fixed). */
  scheme?: string;
}

const DEFAULT_OPTS: RenderOpts = {showLetters: true, pairAlign: true};

/* Visual tuning. */
const ASPECT_H_OVER_W = 1.25;
const STRAND_GAP_RATIO = 0.5;
const CHIP_GAP_RATIO = 0.32; // wider — gives PS bar room to breathe
const MIN_CHIP_W = 5;
const MAX_CHIP_W = 17;
const PAD = 4;
const LABEL_W = 30;
const CHIP_FILL_ALPHA = 0.85; // chip body, base-canonical pale color
const CHIP_BORDER_W = 0.5;
const SUGAR_STRIPE_RATIO = 0.22; // height of bottom sugar-mod stripe
const PS_BAR_RATIO = 0.55; // PS bar width as fraction of chip gap

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
  textOnlyFallback: boolean;
  senseChips: ChipPos[];
  antiChips: ChipPos[];
  senseLinks: LinkagePos[];
  antiLinks: LinkagePos[];
  /** Whether antisense was rendered in reversed (3'→5') order. */
  antiReversed: boolean;
}

export interface ChipPos {
  x: number;
  w: number;
  monomer: ParsedMonomer;
  /** Original 0-based index in the data strand (regardless of display order). */
  origIdx: number;
  strand: StrandSide;
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

  const availW = Math.max(0, cellW - LABEL_W - 2 * PAD);
  const wChipW = availW / (maxLen + Math.max(0, maxLen - 1) * CHIP_GAP_RATIO);

  const heightFactor = strandsCount + Math.max(0, strandsCount - 1) * STRAND_GAP_RATIO;
  const hChipH = (cellH - 2 * PAD) / heightFactor;
  const hChipW = hChipH / ASPECT_H_OVER_W;

  let chipW = Math.min(MAX_CHIP_W, wChipW, hChipW);
  const textOnlyFallback = chipW < MIN_CHIP_W;
  if (chipW < MIN_CHIP_W) chipW = MIN_CHIP_W;

  const chipH = chipW * ASPECT_H_OVER_W;
  const chipGap = Math.max(2, chipW * CHIP_GAP_RATIO);
  const strandGap = chipH * STRAND_GAP_RATIO;
  const fontSize = Math.max(7, Math.min(13, chipW * 0.62));

  const blockH = strandsCount * chipH + (strandsCount - 1) * strandGap;
  const blockTop = Math.max(PAD, (cellH - blockH) / 2);
  const senseY = blockTop;
  const antiY = hasAnti ? blockTop + chipH + strandGap : -1;
  const seqX = PAD + LABEL_W;
  const seqEndX = cellW - PAD;

  const layoutBase: Omit<DuplexLayout, 'senseChips' | 'antiChips' | 'senseLinks' | 'antiLinks' | 'antiReversed'> = {
    chipW, chipH, chipGap, strandGap, fontSize,
    labelW: LABEL_W, padding: PAD,
    senseY, antiY, seqX, textOnlyFallback,
  };

  if (textOnlyFallback) {
    return {
      ...layoutBase,
      senseChips: [], antiChips: [], senseLinks: [], antiLinks: [],
      antiReversed: false,
    };
  }

  // Conjugates on either strand's leftmost end push back the start of
  // *that* strand's nucleotide region. To keep base pairs aligned across
  // strands, we measure the leading-conjugate width of each strand (in
  // its display order — antisense is reversed when pair-aligned) and
  // shift each strand right by the difference, so the first NUCLEOTIDE
  // of each strand sits at the same X coordinate.
  const antiReversed = hasAnti && o.pairAlign;
  const senseLeadW = leadingConjugateWidth(model.sense.monomers, false, chipW, fontSize, chipGap);
  const antiLeadW = hasAnti ?
    leadingConjugateWidth(model.antisense!.monomers, antiReversed, chipW, fontSize, chipGap) : 0;
  const alignAt = Math.max(senseLeadW, antiLeadW);
  const senseStartX = seqX + (alignAt - senseLeadW);
  const antiStartX = seqX + (alignAt - antiLeadW);

  // Per-chip widths: nucleotides with multi-char bases (e.g. `cpm6A`) want a
  // wider chip so the full symbol fits. We try the wide layout first; if the
  // total doesn't fit in the cell, we fall back to uniform chipW + ellipsis.
  // For pair alignment, columns sync across strands (pair-aligned column
  // takes the max width of the two strands' chips at that pair-index).
  const senseDisplay = model.sense.monomers;
  const antiDisplay = hasAnti ?
    (antiReversed ? model.antisense!.monomers.slice().reverse() : model.antisense!.monomers) :
    [];
  let widths = computePairSyncedWidths(senseDisplay, antiDisplay, chipW, fontSize);
  // Falls back to uniform chipW if either strand's chips don't fit.
  const widthBudget = seqEndX - (seqX + alignAt);
  if (!fitsInBudget(widths.sense, senseLeadW, chipGap, widthBudget) ||
      !fitsInBudget(widths.anti, antiLeadW, chipGap, widthBudget)) {
    widths = uniformWidths(senseDisplay, antiDisplay, chipW, fontSize);
  }

  const senseRes = placeStrand(
    model.sense, false, senseY, senseStartX, seqEndX, layoutBase, 'sense', widths.sense);
  const antiRes = hasAnti ? placeStrand(
    model.antisense!, antiReversed, antiY, antiStartX, seqEndX, layoutBase, 'antisense', widths.anti,
  ) : {chips: [], links: []};

  return {
    ...layoutBase,
    senseChips: senseRes.chips,
    antiChips: antiRes.chips,
    senseLinks: senseRes.links,
    antiLinks: antiRes.links,
    antiReversed,
  };
}

/** Place chips for one strand, optionally reversed. Returns chip and linkage positions.
 * `chipWidths` is per-display-monomer (i.e. matches the order of `strand.monomers` after
 * reversal if `reverse=true`). */
function placeStrand(
  strand: ParsedStrand, reverse: boolean, y: number, startX: number, endX: number,
  layout: Omit<DuplexLayout, 'senseChips' | 'antiChips' | 'senseLinks' | 'antiLinks' | 'antiReversed'>,
  side: StrandSide,
  chipWidths: number[],
): { chips: ChipPos[]; links: LinkagePos[] } {
  const monomers = reverse ? strand.monomers.slice().reverse() : strand.monomers;
  const chips: ChipPos[] = [];
  const links: LinkagePos[] = [];
  let x = startX;

  for (let i = 0; i < monomers.length; i++) {
    const m = monomers[i];
    const w = chipWidths[i] ?? layout.chipW;

    if (x + w > endX) break; // truncate at cell edge

    chips.push({x, w, monomer: m, origIdx: m.position, strand: side});

    // Linkage from this monomer's 3'-phosphate goes in the gap to the right.
    // When the strand is reversed (anti-parallel), the linkage that was
    // "owned by monomer N's 3' end" connects monomer N to monomer N+1 in the
    // data; in display these are still adjacent, so the gap-to-right is still
    // the right place to draw it. We just need to look up the correct owner.
    if (m.kind === 'nucleotide' && i < monomers.length - 1) {
      const nt = m as ParsedNucleotide;
      if (nt.phosphate) {
        links.push({
          x: x + w, w: layout.chipGap, y, h: layout.chipH,
          phosphateSymbol: nt.phosphate,
          ownerOrigIdx: nt.position,
          strand: side,
        });
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

/** Whether the per-chip widths array fits in `budget` (after the leading shift). */
function fitsInBudget(widths: number[], leadW: number, chipGap: number, budget: number): boolean {
  let total = leadW;
  for (let i = 0; i < widths.length; i++) {
    total += widths[i];
    if (i < widths.length - 1) total += chipGap;
  }
  return total <= budget;
}

/* ---------------------------------------------------------------- *
 * Drawing
 * ---------------------------------------------------------------- */

export function drawDuplex(
  g: CanvasRenderingContext2D, cellX: number, cellY: number,
  cellW: number, cellH: number, model: ParsedDuplex,
  opts: Partial<RenderOpts> = {},
): DuplexLayout {
  const o: RenderOpts = {...DEFAULT_OPTS, ...opts};
  const layout = computeLayout(cellW, cellH, model, o);

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

  // Linkages first (so chips paint over their rounded edges cleanly)
  for (const link of layout.senseLinks) drawLinkage(g, link);
  drawChips(g, layout.senseChips, layout, o);
  drawTruncationMarker(g, layout.senseChips, model.sense.monomers.length, layout);

  if (layout.antiY >= 0 && model.antisense) {
    // When reversed, the leftmost chip in display is the 3' end of antisense.
    const leftLabel = layout.antiReversed ? '3\'' : '5\'';
    drawStrandLabel(g, 'AS', leftLabel, layout.padding, layout.antiY + layout.chipH / 2, layout);
    for (const link of layout.antiLinks) drawLinkage(g, link);
    drawChips(g, layout.antiChips, layout, o);
    drawTruncationMarker(g, layout.antiChips, model.antisense.monomers.length, layout);
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
  const y = chips[0]?.strand === 'sense' ? layout.senseY : layout.antiY;
  for (const cp of chips) {
    if (cp.monomer.kind === 'conjugate')
      drawConjugate(g, cp.monomer.symbol, cp.x, y, cp.w, layout.chipH, layout.fontSize);
    else
      drawChip(g, cp.monomer as ParsedNucleotide, cp.x, y, cp.w, layout.chipH, layout.fontSize, opts);
  }
}

function drawChip(
  g: CanvasRenderingContext2D, m: ParsedNucleotide, x: number, y: number,
  w: number, h: number, fontSize: number, opts: RenderOpts,
): void {
  const sugarRes = resolveSugar(m.sugar, m.base);
  const baseColor = resolveBaseColor(m.base);
  const isModSugar = isModifiedSugar(m.sugar);
  const r = Math.min(2.5, w / 4);

  // Chip body — base-canonical pale color (or analog's color for custom bases)
  drawRoundRect(g, x, y, w, h, r);
  g.fillStyle = withAlpha(baseColor, CHIP_FILL_ALPHA);
  g.fill();
  g.lineWidth = CHIP_BORDER_W;
  g.strokeStyle = 'rgba(0,0,0,0.22)';
  g.stroke();

  // Sugar modification stripe at chip bottom — clipped to rounded shape
  if (isModSugar) {
    const stripeH = Math.max(2, h * SUGAR_STRIPE_RATIO);
    g.save();
    drawRoundRect(g, x, y, w, h, r);
    g.clip();
    g.fillStyle = sugarRes.color;
    g.fillRect(x, y + h - stripeH, w, stripeH);
    g.restore();
  }

  // Base label — biased upward to leave room for stripe.
  // Pick full label or shortened (`X…`) based on whether the chip is wide
  // enough at the current fontSize. This way, when the layout granted us a
  // wider chip, the user sees the full HELM symbol; otherwise it's clipped
  // to first-letter + ellipsis for legibility.
  if (opts.showLetters && m.base && fontSize >= 8) {
    const stripeH = isModSugar ? Math.max(2, h * SUGAR_STRIPE_RATIO) : 0;
    const label = pickBaseLabel(m.base, w, fontSize);
    g.fillStyle = '#1a1a1a';
    g.font = `600 ${fontSize}px system-ui, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif`;
    g.textBaseline = 'middle';
    g.textAlign = 'center';
    g.fillText(label, x + w / 2, y + (h - stripeH) / 2 + 0.5);
  }
}

/** Resolve a chip background color for any base, including custom (multi-char)
 * symbols whose color follows their natural analog (`A` / `C` / `G` / `U` / `T`).
 * Sync — backed by the central Bio monomer library that's wired up at
 * package init. */
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

function isModifiedSugar(sugar: string): boolean {
  const c = canonicalSugarSymbol(sugar);
  return c !== 'r' && c !== 'd';
}

function drawLinkage(g: CanvasRenderingContext2D, link: LinkagePos): void {
  // Only the canonical phosphate (`p` / aliased `P`) is treated as "no marker".
  // Every other linkage — known PS / PS₂ / MeP, or unknown custom symbol that
  // got a hash-derived color — gets a bar in the inter-chip gap so the user
  // can see and hover it. Color comes from resolvePhosphate which is
  // deterministic per symbol, so two distinct unknown symbols get distinct bars.
  const canonical = canonicalPhosphateSymbol(link.phosphateSymbol);
  if (canonical === 'p') return;
  const ps = resolvePhosphate(link.phosphateSymbol);
  const barW = Math.max(2.5, link.w * PS_BAR_RATIO);
  const barX = link.x + (link.w - barW) / 2;
  g.fillStyle = ps.color;
  g.fillRect(barX, link.y + 1, barW, link.h - 2);
}

function drawConjugate(
  g: CanvasRenderingContext2D, symbol: string, x: number, y: number,
  w: number, chipH: number, fontSize: number,
): void {
  const conj = resolveConjugate(symbol);
  const r = chipH / 2;
  drawRoundRect(g, x, y, w, chipH, r);
  g.fillStyle = conj.color;
  g.fill();
  g.lineWidth = 0.5;
  g.strokeStyle = 'rgba(0,0,0,0.2)';
  g.stroke();

  if (fontSize >= 8) {
    g.fillStyle = '#ffffff';
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
const _alphaCache = new Map<string, string>();
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
  if (_alphaCache.size > 256) _alphaCache.clear();
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

  // Sense chip?
  if (localY >= layout.senseY && localY <= layout.senseY + layout.chipH) {
    const cp = findChip(localX, layout.senseChips);
    if (cp) return {strand: 'sense', position: cp.origIdx, monomer: cp.monomer};
    const link = findLink(localX, localY, layout.senseLinks);
    if (link) return resolveLinkHit(link, model.sense, 'sense');
  }
  // Antisense chip?
  if (layout.antiY >= 0 && localY >= layout.antiY && localY <= layout.antiY + layout.chipH) {
    const cp = findChip(localX, layout.antiChips);
    if (cp) return {strand: 'antisense', position: cp.origIdx, monomer: cp.monomer};
    const link = findLink(localX, localY, layout.antiLinks);
    if (link && model.antisense) return resolveLinkHit(link, model.antisense, 'antisense');
  }
  return null;
}

function findChip(x: number, chips: ChipPos[]): ChipPos | null {
  for (const cp of chips) if (x >= cp.x && x < cp.x + cp.w) return cp;
  return null;
}

function findLink(x: number, y: number, links: LinkagePos[]): LinkagePos | null {
  for (const l of links)
    if (x >= l.x && x < l.x + l.w && y >= l.y && y < l.y + l.h) return l;
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
