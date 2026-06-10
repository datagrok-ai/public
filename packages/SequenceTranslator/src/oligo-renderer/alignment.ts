/**
 * Sense / antisense strand alignment for the OligoNucleotide renderer.
 *
 * Two sources of truth, in priority order:
 *
 *  1. **Explicit HELM pairing** — when the HELM carries a `$connections$`
 *     section with `pair` annotations (and/or strand-type annotations), the
 *     parser turns those into `model.pairs` + `model.alignment`. That is the
 *     authoritative register; we render exactly what the HELM says.
 *
 *  2. **Auto-alignment** — when there's no explicit info, we slide the
 *     (reversed) antisense across the sense and pick the column shift that
 *     maximizes Watson-Crick complementary base pairs. This is *not* a gapped
 *     alignment: it's a single rigid offset in `[-(antiLen-1), senseLen-1]`,
 *     which is exactly what an siRNA / ASO duplex with terminal overhangs needs.
 *     Custom / modified bases are reduced to their natural analog first
 *     (`cpm6A` → `A`, `5meC` → `C`, …) so a modified strand still pairs.
 *
 * The shift convention matches `DuplexAlignment.shift` in `types.ts`: it is the
 * column offset of the antisense *display* (3'→5') relative to sense.
 *
 * Pure module — no DG / DOM. The natural-analog lookup goes through
 * `analog-cache`, which gracefully returns `null` before the monomer library
 * is loaded (so canonical bases still align in that window, and in unit tests).
 */

import {
  DuplexAlignment, ParsedDuplex, ParsedNucleotide, ParsedStrand,
  isCanonicalBase,
} from './types';
import {getNaturalAnalog} from './analog-cache';

/** Reduce any base symbol to its canonical A/C/G/U/T letter (uppercase) for
 * pairing, via the natural-analog lookup. Returns null when missing / no
 * resolvable analog. Shared with the canvas renderer's pair-kind logic. */
export function reduceToCanonicalBase(base: string | null | undefined): string | null {
  if (!base) return null;
  if (isCanonicalBase(base)) return base;
  const analog = getNaturalAnalog(base);
  return analog && isCanonicalBase(analog) ? analog : null;
}

/** Watson-Crick complementarity on canonical-reduced bases (A↔U/T, G↔C).
 * Wobble (G·U) is intentionally NOT counted — overhang detection wants the
 * cleanest complementary register, and counting wobble blurs the argmax. */
export function basesComplementary(a: string | null | undefined, b: string | null | undefined): boolean {
  const ca = reduceToCanonicalBase(a);
  const cb = reduceToCanonicalBase(b);
  if (!ca || !cb) return false;
  return (
    (ca === 'A' && (cb === 'U' || cb === 'T')) ||
    ((ca === 'U' || ca === 'T') && cb === 'A') ||
    (ca === 'G' && cb === 'C') ||
    (ca === 'C' && cb === 'G')
  );
}

/** Bases (in 5'→3' data order) of a strand's nucleotides, conjugates skipped. */
function strandBases(strand: ParsedStrand): (string | null)[] {
  const out: (string | null)[] = [];
  for (const m of strand.monomers)
    if (m.kind === 'nucleotide') out.push((m as ParsedNucleotide).base);
  return out;
}

/** Number of overlapping columns (positions where both strands have a
 * nucleotide) for a given shift. */
export function overlapLength(shift: number, senseLen: number, antiLen: number): number {
  // sense column j is in [0, senseLen-1] after offset; antisense data k maps to
  // sense column (antiLen-1-k)+shift. Overlap = count of j with a partner k in range.
  let n = 0;
  for (let j = 0; j < senseLen; j++) {
    const k = antiLen - 1 - j + shift;
    if (k >= 0 && k < antiLen) n++;
  }
  return n;
}

export interface ShiftScan {
  shift: number;
  /** Complementary base pairs at `shift`. */
  score: number;
  /** Complementary base pairs at the blunt (shift = 0) register. */
  score0: number;
}

/** Slide the (reversed) antisense across sense and find the shift with the most
 * Watson-Crick complementary pairs. Tie-break: fewer-overhang (smaller |shift|),
 * then lower shift, for determinism. Operates on canonical-reduced bases.
 *
 * `senseBases` / `antiBases` are in 5'→3' data order. For each candidate shift,
 * sense data index `j` pairs antisense data index `k = antiLen-1-j+shift`. */
export function bestComplementaryShift(
  senseBases: (string | null)[], antiBases: (string | null)[],
): ShiftScan {
  const sLen = senseBases.length;
  const aLen = antiBases.length;
  const scoreAt = (shift: number): number => {
    let score = 0;
    for (let j = 0; j < sLen; j++) {
      const k = aLen - 1 - j + shift;
      if (k < 0 || k >= aLen) continue;
      if (basesComplementary(senseBases[j], antiBases[k])) score++;
    }
    return score;
  };

  const score0 = scoreAt(0);
  let best = {shift: 0, score: score0};
  // Range covers every register from full left-overhang to full right-overhang.
  for (let shift = -(aLen - 1); shift <= sLen - 1; shift++) {
    const score = scoreAt(shift);
    if (
      score > best.score ||
      (score === best.score && Math.abs(shift) < Math.abs(best.shift)) ||
      (score === best.score && Math.abs(shift) === Math.abs(best.shift) && shift < best.shift)
    )
      best = {shift, score};
  }
  return {shift: best.shift, score: best.score, score0};
}

/** Minimum complementary pairs before we trust a *non-zero* auto shift. Below
 * this, or when the shifted register is only marginally better than blunt, we
 * fall back to a blunt (shift 0) duplex — this keeps poorly-complementary or
 * synthetic inputs rendering predictably instead of snapping to a spurious
 * offset. A genuine siRNA / ASO duplex clears this comfortably. */
const MIN_ACCEPT_PAIRS = 3;
const MIN_ACCEPT_FRACTION = 0.6;
const MIN_GAIN_OVER_BLUNT = 2;

/** Auto-align two strands by complementarity. Returns the column shift to apply
 * to the antisense display. Returns 0 (blunt) when there's no antisense, no
 * clear complementary register, or the best register barely beats blunt. */
export function computeAutoShift(sense: ParsedStrand, antisense: ParsedStrand | null): number {
  if (!antisense) return 0;
  const senseBases = strandBases(sense);
  const antiBases = strandBases(antisense);
  if (senseBases.length === 0 || antiBases.length === 0) return 0;

  const {shift, score, score0} = bestComplementaryShift(senseBases, antiBases);
  if (shift === 0) return 0;

  const overlap = overlapLength(shift, senseBases.length, antiBases.length);
  const strongEnough =
    score >= MIN_ACCEPT_PAIRS &&
    score >= MIN_ACCEPT_FRACTION * overlap &&
    score >= score0 + MIN_GAIN_OVER_BLUNT;
  return strongEnough ? shift : 0;
}

/** Resolve the alignment shift to render with. Explicit HELM pairing wins;
 * otherwise auto-align by complementarity. Single-strand → blunt/none. */
export function resolveDuplexAlignment(model: ParsedDuplex): DuplexAlignment {
  if (!model.antisense) return {shift: 0, source: 'none'};
  if (model.alignment && model.alignment.source === 'explicit')
    return model.alignment;
  return {shift: computeAutoShift(model.sense, model.antisense), source: 'auto'};
}

export interface OverhangInfo {
  shift: number;
  source: DuplexAlignment['source'];
  /** Unpaired nucleotides hanging off the sense 5' end (left, sense side). */
  sense5: number;
  /** Unpaired nucleotides hanging off the sense 3' end (right, sense side). */
  sense3: number;
  /** Unpaired nucleotides hanging off the antisense 5' end (right in display). */
  anti5: number;
  /** Unpaired nucleotides hanging off the antisense 3' end (left in display). */
  anti3: number;
  /** Number of base-paired columns. */
  paired: number;
}

/** Column offsets for a given shift (see `DuplexAlignment.shift`). The strand
 * that sticks out furthest left gets offset 0; the other is pushed right. */
export function columnOffsets(shift: number): { sense: number; anti: number } {
  return shift >= 0 ? {sense: 0, anti: shift} : {sense: -shift, anti: 0};
}

/** Describe overhangs implied by an alignment, for panels / tooltips. */
export function describeOverhangs(model: ParsedDuplex): OverhangInfo | null {
  if (!model.antisense) return null;
  const senseLen = strandBases(model.sense).length;
  const antiLen = strandBases(model.antisense).length;
  const {shift, source} = resolveDuplexAlignment(model);
  const {sense: sOff, anti: aOff} = columnOffsets(shift);

  const senseStart = sOff; const senseEnd = sOff + senseLen; // [start, end)
  const antiStart = aOff; const antiEnd = aOff + antiLen;
  const overlapStart = Math.max(senseStart, antiStart);
  const overlapEnd = Math.min(senseEnd, antiEnd);
  const paired = Math.max(0, overlapEnd - overlapStart);

  // Left side: whichever strand starts further left overhangs there.
  const left = Math.abs(senseStart - antiStart);
  // Right side: whichever strand ends further right overhangs there.
  const right = Math.abs(senseEnd - antiEnd);

  // Sense reads 5'→3' left→right; antisense display reads 3'→5' left→right,
  // so the antisense *left* end is its 3' end and its *right* end is its 5' end.
  const senseLeftOwner = senseStart <= antiStart;
  const senseRightOwner = senseEnd >= antiEnd;
  return {
    shift, source,
    sense5: senseLeftOwner ? left : 0,
    anti3: senseLeftOwner ? 0 : left,
    sense3: senseRightOwner ? right : 0,
    anti5: senseRightOwner ? 0 : right,
    paired,
  };
}
