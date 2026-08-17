import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

import {Fingerprint, rdKitFingerprintToBitArray} from '../../utils/chem-common';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {SarMatrix, SarMatrixCell, logSarTime} from './sar-matrix-types';

/**
 * SAR-transfer detection: a transfer holds between two analog series (each a core varied at one
 * R-position) whose potencies track each other across the compounds the two series have in common.
 *
 * Two conditions, and both are load-bearing.
 *
 * The SUBSTITUENT must be the same R-group, matched by identity. That is what makes the comparison
 * controlled: the claim is that one particular change did one particular thing on two scaffolds, and it
 * only holds if both series actually made that change. A merely similar R-group does not support it.
 *
 * The COMPOUNDS must clear a Morgan-similarity floor. Two molecules can share a substituent and have
 * nothing else in common, and carrying a result between unrelated chemotypes on that basis is not a
 * transfer — over a handful of points their agreement is coincidence. Since the substituent is held
 * identical, what this threshold really measures is how related the two cores are, which is exactly the
 * question of whether a transfer between them is plausible.
 *
 * Only pairs from DIFFERENT matrices are reported. Two cores of one matrix are held together by that
 * matrix's own additive row-plus-column model, so their trends run parallel by construction — reporting
 * it states the assumption back, and since it scores near-perfectly it would crowd the genuinely
 * informative cross-chemotype pairs out of the list.
 *
 * Agreement is Spearman rank correlation: a transfer claims the two cores rank their shared substituents
 * the same way, which rank correlation tests directly and a single outlier pair cannot fake.
 */

/** Minimum shared substituents between two series to compare their trends. This is an overlap between
 *  different-scaffold series, not a series length — held low so genuine cross-chemotype transfers are
 *  not rejected. At the 0.7 floor, three points already demand a perfectly monotone trend. */
const MIN_COMMON = 3;
/** Minimum Spearman rank correlation for a transfer to be reported. */
const MIN_CORRELATION = 0.7;
/** Cap on reported transfers, strongest first. */
const MAX_TRANSFERS = 16;
/** Display cap on predicted columns. A series can carry hundreds of predictions and they are read
 *  alongside the measured ones, so the trend view takes only the most promising — ordered first, so
 *  the cut can never hide the analog most worth making. */
const MAX_PREDICTED = 6;
/** Morgan Tanimoto a matched pair's compounds must reach. Held well below the value that would mean
 *  "nearly the same molecule": the two carry the same R-group on deliberately different scaffolds, so
 *  demanding near-identity would reject every cross-chemotype pair the scan exists to find. */
export const DEFAULT_TRANSFER_SIMILARITY = 0.5;

/** One side of a transfer: a core (a matrix row) varied at one R-position. */
export interface TransferSide {
  matrixIndex: number;
  rowIndex: number;
  position: string;
}

/** A predicted analog one side of a transfer still has open — what the transfer argues for making. */
export interface VirtualAnalog {
  col: number;
  substSmiles: string;
  value: number;
}

/** A detected transfer between two cores, possibly in different matrices. */
export interface Transfer {
  a: TransferSide;
  b: TransferSide;
  /** The R-group of each matched pair. Both sides carry it — that is what being matched means — so the
   *  one label is true of the columns under it on both rows. */
  substituents: string[];
  /** Column indices in matrix a / b for each matched pair (aligned to `substituents`). */
  aCols: number[];
  bCols: number[];
  /** a's / b's cell for each matched pair (aligned to `substituents`). */
  aCells: SarMatrixCell[];
  bCells: SarMatrixCell[];
  /** Mean compound Tanimoto over the matched pairs. The R-groups are identical by construction, so
   *  what this measures is how alike the two scaffolds carrying them are. */
  similarity: number;
  correlation: number;
  /** Predicted analogs open on each side, in column order. */
  aVirtual: VirtualAnalog[];
  bVirtual: VirtualAnalog[];
  /** R-groups one side measured and the other has only predicted — what the transfer argues for
   *  making. Kept apart from the matched pairs: these carry no second observation, so they are shown
   *  but take no part in the correlation or the fold-change match. */
  predictedSubstituents: string[];
  predictedACols: number[];
  predictedBCols: number[];
}

/** One observed compound of a series: where it sits, what it measured, and what it varies. */
interface SeriesPoint {
  col: number;
  cell: SarMatrixCell;
  value: number;
  substSmiles: string;
  /** Index into the shared compound fingerprint table. */
  fp: number;
  /** The compound itself, compared for identity — one compound cannot corroborate itself. */
  smiles: string;
}

/** One core varied at one position: its observed compounds and the analogs it still has open. */
interface RowSeries {
  matrixIndex: number;
  rowIndex: number;
  position: string;
  points: SeriesPoint[];
  virtual: VirtualAnalog[];
}

function pearson(xs: number[], ys: number[]): number | null {
  const n = xs.length;
  if (n < MIN_COMMON)
    return null;
  const mx = xs.reduce((s, v) => s + v, 0) / n;
  const my = ys.reduce((s, v) => s + v, 0) / n;
  let cov = 0;
  let vx = 0;
  let vy = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx;
    const dy = ys[i] - my;
    cov += dx * dy;
    vx += dx * dx;
    vy += dy * dy;
  }
  if (vx === 0 || vy === 0)
    return null;
  return cov / Math.sqrt(vx * vy);
}

/** Average ranks (1-based; ties share their mean rank) — the transform that turns Pearson into
 *  Spearman. */
function averageRanks(values: number[]): number[] {
  const order = values.map((v, i) => ({v, i})).sort((a, b) => a.v - b.v);
  const ranks = new Array<number>(values.length);
  let i = 0;
  while (i < order.length) {
    let j = i;
    while (j + 1 < order.length && order[j + 1].v === order[i].v)
      j++;
    const rank = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++)
      ranks[order[k].i] = rank;
    i = j + 1;
  }
  return ranks;
}

/** Spearman rank correlation: Pearson over the two value orderings. Null on too few points or a
 *  constant side. */
export function spearman(xs: number[], ys: number[]): number | null {
  return pearson(averageRanks(xs), averageRanks(ys));
}

/** For one matrix: each R-position's substituents (first column wins per substituent), in column
 *  order, with the column index of each. Independent of the row, so it is computed once per matrix. */
function positionColumns(matrix: SarMatrix): Map<string, {order: string[], cols: number[]}> {
  const byPosition = new Map<string, {order: string[], cols: number[], seen: Set<string>}>();
  matrix.columns.forEach((col, ci) => {
    let group = byPosition.get(col.position);
    if (!group) {
      group = {order: [], cols: [], seen: new Set()};
      byPosition.set(col.position, group);
    }
    if (group.seen.has(col.substSmiles))
      return;
    group.seen.add(col.substSmiles);
    group.order.push(col.substSmiles);
    group.cols.push(ci);
  });
  return byPosition;
}

/**
 * Morgan fingerprints for every observed compound in the matrices, keyed by SMILES.
 *
 * Built once for the whole scan: the same compound sits in as many series as it has cut sites, and
 * fingerprinting is the one part of this that pays a round trip to the RDKit workers.
 */
async function fingerprintTable(matrices: SarMatrix[]):
  Promise<{index: Map<string, number>, fps: (BitArray | null)[]}> {
  const index = new Map<string, number>();
  const smiles: string[] = [];
  for (const matrix of matrices) {
    for (const rowCells of matrix.cells) {
      for (const cell of rowCells) {
        if (cell.kind === 'real' && cell.smiles !== null && !index.has(cell.smiles)) {
          index.set(cell.smiles, smiles.length);
          smiles.push(cell.smiles);
        }
      }
    }
  }
  if (!smiles.length)
    return {index, fps: []};
  const service = await getRdKitService();
  const raw = (await service.getFingerprints(Fingerprint.Morgan, smiles, false)).fps;
  return {index, fps: raw.map((fp) => fp ? rdKitFingerprintToBitArray(fp) : null)};
}

/**
 * Every (core, position) as an analog series, carrying its observed compounds and its open analogs.
 *
 * A series holding fewer observations than a transfer needs is dropped here rather than compared:
 * matching is one-to-one, so it could never reach `MIN_COMMON` pairs however good its partner is.
 */
function gatherSeries(matrices: SarMatrix[], fpIndex: Map<string, number>): RowSeries[] {
  const series: RowSeries[] = [];
  matrices.forEach((matrix, matrixIndex) => {
    for (const [position, {order, cols}] of positionColumns(matrix)) {
      matrix.rows.forEach((_row, rowIndex) => {
        const points: SeriesPoint[] = [];
        const virtual: VirtualAnalog[] = [];
        order.forEach((substSmiles, k) => {
          const col = cols[k];
          const cell = matrix.cells[rowIndex][col];
          if (cell.value === null)
            return;
          if (cell.kind === 'real' && cell.smiles !== null) {
            const fp = fpIndex.get(cell.smiles);
            if (fp !== undefined)
              points.push({col, cell, value: cell.value, substSmiles, fp, smiles: cell.smiles});
          } else if (cell.kind === 'virtual')
            virtual.push({col, substSmiles, value: cell.value});
        });
        if (points.length >= MIN_COMMON)
          series.push({matrixIndex, rowIndex, position, points, virtual});
      });
    }
  });
  return series;
}

/**
 * Pair the two series up: one pairing per R-group both have explored, kept only when the two compounds
 * carrying it are alike enough to count as related chemotypes.
 *
 * Substituents are matched by identity, which makes the correspondence one-to-one for free — a series
 * lists each R-group once — so nothing has to arbitrate between competing partners. The result comes
 * back in a's column order, which the trend view reads left to right and the fold-change match walks in
 * consecutive steps.
 */
function matchPoints(a: RowSeries, b: RowSeries, threshold: number,
  sim: (i: number, j: number) => number): {ai: number, bi: number, s: number}[] {
  const bBySubst = new Map<string, number>();
  b.points.forEach((point, j) => bBySubst.set(point.substSmiles, j));
  const matched: {ai: number, bi: number, s: number}[] = [];
  for (let i = 0; i < a.points.length; i++) {
    const j = bBySubst.get(a.points[i].substSmiles);
    if (j === undefined)
      continue;
    // One compound is not evidence that two cores track each other. Matrices are an overlapping cover,
    // so the same compound reaches both series routinely; matched against itself it puts its own
    // potency on both axes, which forces r to 1. Compared by structure rather than by row, since one
    // compound measured on two rows is still one compound.
    if (a.points[i].smiles === b.points[j].smiles)
      continue;
    const s = sim(a.points[i].fp, b.points[j].fp);
    if (s >= threshold)
      matched.push({ai: i, bi: j, s});
  }
  return matched;
}

/**
 * The R-groups one series has measured and the other has only predicted — the analogs this transfer
 * argues for making. Most promising first, so the display cap cannot hide the best one.
 */
function predictedPairs(a: RowSeries, b: RowSeries, higherIsBetter: boolean):
  {subst: string, aCol: number, bCol: number, value: number}[] {
  const out: {subst: string, aCol: number, bCol: number, value: number}[] = [];
  const collect = (measured: RowSeries, open: RowSeries, measuredIsA: boolean): void => {
    const bySubst = new Map<string, VirtualAnalog>();
    for (const analog of open.virtual)
      bySubst.set(analog.substSmiles, analog);
    for (const point of measured.points) {
      const analog = bySubst.get(point.substSmiles);
      if (analog === undefined)
        continue;
      out.push({
        subst: point.substSmiles,
        aCol: measuredIsA ? point.col : analog.col,
        bCol: measuredIsA ? analog.col : point.col,
        value: analog.value,
      });
    }
  };
  collect(a, b, true);
  collect(b, a, false);
  out.sort((p, q) => higherIsBetter ? q.value - p.value : p.value - q.value);
  return out.slice(0, MAX_PREDICTED);
}

/** Unordered key for a core pair, so R1 and R2 transfers between the same two cores dedupe to one. */
function pairKey(a: TransferSide, b: TransferSide): string {
  const x = `${a.matrixIndex}:${a.rowIndex}`;
  const y = `${b.matrixIndex}:${b.rowIndex}`;
  return x < y ? `${x}|${y}` : `${y}|${x}`;
}

/**
 * Every SAR transfer across the whole set of matrices, strongest first: for each pair of series at the
 * same R-position, every R-group both explored is paired up, the pairs whose compounds are too unalike
 * are dropped, and the Pearson correlation of the rest is kept when it clears `MIN_CORRELATION` over at
 * least `MIN_COMMON` pairs. The best-correlating R-position is kept per core pair.
 *
 * @param matrices Matrices to scan.
 * @param similarity Compound Tanimoto a pair must reach for the two scaffolds to count as related.
 * @param higherIsBetter Which end of the activity is more potent, for ordering the predicted analogs.
 */
export async function computeAllTransfers(matrices: SarMatrix[],
  similarity: number = DEFAULT_TRANSFER_SIMILARITY, higherIsBetter = true): Promise<Transfer[]> {
  const t0 = performance.now();
  const {index, fps} = await fingerprintTable(matrices);
  const sim = (i: number, j: number): number => {
    const x = fps[i];
    const y = fps[j];
    return x && y ? tanimotoSimilarity(x, y) : 0;
  };
  const series = gatherSeries(matrices, index);
  const bestByPair = new Map<string, Transfer>();
  for (let i = 0; i < series.length; i++) {
    for (let j = i + 1; j < series.length; j++) {
      const s1 = series[i];
      const s2 = series[j];
      // Two cores of one matrix track each other by construction (one additive model across its rows).
      if (s1.matrixIndex === s2.matrixIndex || s1.position !== s2.position)
        continue;
      const matched = matchPoints(s1, s2, similarity, sim);
      if (matched.length < MIN_COMMON)
        continue;
      const corr = spearman(matched.map((m) => s1.points[m.ai].value), matched.map((m) => s2.points[m.bi].value));
      if (corr === null || corr < MIN_CORRELATION)
        continue;
      const a: TransferSide = {matrixIndex: s1.matrixIndex, rowIndex: s1.rowIndex, position: s1.position};
      const b: TransferSide = {matrixIndex: s2.matrixIndex, rowIndex: s2.rowIndex, position: s2.position};
      const key = pairKey(a, b);
      const existing = bestByPair.get(key);
      if (existing && existing.correlation >= corr)
        continue; // keep the best-correlated R-position for this core pair
      const predicted = predictedPairs(s1, s2, higherIsBetter);
      bestByPair.set(key, {
        a, b,
        substituents: matched.map((m) => s1.points[m.ai].substSmiles),
        aCols: matched.map((m) => s1.points[m.ai].col),
        bCols: matched.map((m) => s2.points[m.bi].col),
        aCells: matched.map((m) => s1.points[m.ai].cell),
        bCells: matched.map((m) => s2.points[m.bi].cell),
        similarity: matched.reduce((sum, m) => sum + m.s, 0) / matched.length,
        correlation: corr,
        aVirtual: s1.virtual,
        bVirtual: s2.virtual,
        predictedSubstituents: predicted.map((p) => p.subst),
        predictedACols: predicted.map((p) => p.aCol),
        predictedBCols: predicted.map((p) => p.bCol),
      });
    }
  }
  logSarTime(`transfers (${series.length} series, ${fps.length} compounds)`, t0);
  return [...bestByPair.values()].sort((p, q) => q.correlation - p.correlation).slice(0, MAX_TRANSFERS);
}

export interface TransferStats {
  correlation: number;
  /** Per-step effect-size agreement (0-1) over consecutive matched pairs, or null. */
  foldMatch: number | null;
  /** The follower core's predicted analog the transfer argues for making, or null. */
  benefiting: {side: 'a' | 'b', substSmiles: string, value: number} | null;
}

function observed(cell: SarMatrixCell): number | null {
  return cell.kind === 'real' && cell.value !== null ? cell.value : null;
}

/**
 * Detailed statistics for a transfer: the correlation, the fold-change match (per-step effect-size
 * agreement between the two cores), and the "benefiting" analog — the untested compound one core is
 * still missing that the transferred rule predicts, preferring the follower core `b`.
 */
export function transferStats(transfer: Transfer, higherIsBetter: boolean): TransferStats {
  const deltas: {a: number, b: number}[] = [];
  for (let k = 1; k < transfer.substituents.length; k++) {
    const a0 = observed(transfer.aCells[k - 1]);
    const a1 = observed(transfer.aCells[k]);
    const b0 = observed(transfer.bCells[k - 1]);
    const b1 = observed(transfer.bCells[k]);
    if (a0 !== null && a1 !== null && b0 !== null && b1 !== null)
      deltas.push({a: a1 - a0, b: b1 - b0});
  }
  let foldMatch: number | null = null;
  if (deltas.length) {
    let sum = 0;
    for (const {a, b} of deltas) {
      const sameDirection = (a >= 0) === (b >= 0);
      const mag = Math.max(Math.abs(a), Math.abs(b));
      sum += sameDirection ? (mag === 0 ? 1 : Math.min(Math.abs(a), Math.abs(b)) / mag) : 0;
    }
    foldMatch = sum / deltas.length;
  }

  let benefiting: {side: 'a' | 'b', substSmiles: string, value: number} | null = null;
  for (const [side, analogs] of [['b', transfer.bVirtual], ['a', transfer.aVirtual]] as const) {
    let best = higherIsBetter ? -Infinity : Infinity;
    let pick: VirtualAnalog | null = null;
    for (const analog of analogs) {
      if (higherIsBetter ? analog.value > best : analog.value < best) {
        best = analog.value;
        pick = analog;
      }
    }
    if (pick) {
      benefiting = {side, substSmiles: pick.substSmiles, value: pick.value};
      break; // prefer the follower core's analogs
    }
  }
  return {correlation: transfer.correlation, foldMatch, benefiting};
}
