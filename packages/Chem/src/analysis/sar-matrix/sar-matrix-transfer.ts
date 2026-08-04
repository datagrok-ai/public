import {SarMatrix, SarMatrixCell} from './sar-matrix-types';

/**
 * SAR-transfer detection, literature-style: a transfer holds between two analog series (each a core
 * varied at one R-position) whose potencies are correlated across the substituents they SHARE —
 * regardless of whether the two cores are structurally similar or landed in the same cluster. So this
 * scans EVERY core pair across all matrices, catching transfers between differently-scaffolded series
 * (the valuable "carry the SAR to another chemotype" case), not only cores grouped in one matrix.
 */

/** Minimum shared, both-observed substituents for a pair to be considered. */
const MIN_COMMON = 3;
/** Minimum Pearson correlation for a transfer to be reported. */
const MIN_CORRELATION = 0.7;
/** Cap on reported transfers, strongest first. */
const MAX_TRANSFERS = 16;

/** One side of a transfer: a core (a matrix row) varied at one R-position. */
export interface TransferSide {
  matrixIndex: number;
  rowIndex: number;
  position: string;
}

/** A detected transfer between two cores, possibly in different matrices. */
export interface Transfer {
  a: TransferSide;
  b: TransferSide;
  /** Shared substituent SMILES, in a's column order. */
  substituents: string[];
  /** Column indices in matrix a / b for each shared substituent (aligned to `substituents`). */
  aCols: number[];
  bCols: number[];
  /** a's / b's cell for each shared substituent (aligned to `substituents`). */
  aCells: SarMatrixCell[];
  bCells: SarMatrixCell[];
  correlation: number;
  /** True when the two cores come from different matrices (a cross-series transfer). */
  crossMatrix: boolean;
}

/** One core varied at one position, indexed by substituent. */
interface RowSeries {
  matrixIndex: number;
  rowIndex: number;
  position: string;
  order: string[];
  cellBySubst: Map<string, SarMatrixCell>;
  colBySubst: Map<string, number>;
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

/** Every (core, position) as an analog series indexed by substituent. */
function gatherSeries(matrices: SarMatrix[]): RowSeries[] {
  const series: RowSeries[] = [];
  matrices.forEach((matrix, matrixIndex) => {
    for (const [position, {order, cols}] of positionColumns(matrix)) {
      matrix.rows.forEach((_row, rowIndex) => {
        const cellBySubst = new Map<string, SarMatrixCell>();
        const colBySubst = new Map<string, number>();
        order.forEach((sub, k) => {
          cellBySubst.set(sub, matrix.cells[rowIndex][cols[k]]);
          colBySubst.set(sub, cols[k]);
        });
        series.push({matrixIndex, rowIndex, position, order, cellBySubst, colBySubst});
      });
    }
  });
  return series;
}

/** Unordered key for a core pair, so R1 and R2 transfers between the same two cores dedupe to one. */
function pairKey(a: TransferSide, b: TransferSide): string {
  const x = `${a.matrixIndex}:${a.rowIndex}`;
  const y = `${b.matrixIndex}:${b.rowIndex}`;
  return x < y ? `${x}|${y}` : `${y}|${x}`;
}

/**
 * Every SAR transfer across the whole set of matrices, strongest first: for each pair of series at
 * the same R-position (within OR across matrices) that share at least `MIN_COMMON` observed
 * substituents, the Pearson correlation of their potencies; kept when it clears `MIN_CORRELATION`.
 * The best-correlating R-position is kept per core pair.
 */
export function computeAllTransfers(matrices: SarMatrix[]): Transfer[] {
  const series = gatherSeries(matrices);
  const bestByPair = new Map<string, Transfer>();
  for (let i = 0; i < series.length; i++) {
    for (let j = i + 1; j < series.length; j++) {
      const s1 = series[i];
      const s2 = series[j];
      if (s1.matrixIndex === s2.matrixIndex && s1.rowIndex === s2.rowIndex)
        continue; // same core — not a transfer
      if (s1.position !== s2.position)
        continue; // compare like positions so the substituent trends are comparable
      const shared = s1.order.filter((sub) => s2.cellBySubst.has(sub));
      if (shared.length < MIN_COMMON)
        continue;
      const xs: number[] = [];
      const ys: number[] = [];
      for (const sub of shared) {
        const c1 = s1.cellBySubst.get(sub)!;
        const c2 = s2.cellBySubst.get(sub)!;
        if (c1.kind === 'real' && c1.value !== null && c2.kind === 'real' && c2.value !== null) {
          xs.push(c1.value);
          ys.push(c2.value);
        }
      }
      const corr = pearson(xs, ys);
      if (corr === null || corr < MIN_CORRELATION)
        continue;
      const a: TransferSide = {matrixIndex: s1.matrixIndex, rowIndex: s1.rowIndex, position: s1.position};
      const b: TransferSide = {matrixIndex: s2.matrixIndex, rowIndex: s2.rowIndex, position: s2.position};
      const key = pairKey(a, b);
      const existing = bestByPair.get(key);
      if (existing && existing.correlation >= corr)
        continue;
      bestByPair.set(key, {
        a, b, substituents: shared,
        aCols: shared.map((sub) => s1.colBySubst.get(sub)!),
        bCols: shared.map((sub) => s2.colBySubst.get(sub)!),
        aCells: shared.map((sub) => s1.cellBySubst.get(sub)!),
        bCells: shared.map((sub) => s2.cellBySubst.get(sub)!),
        correlation: corr,
        crossMatrix: s1.matrixIndex !== s2.matrixIndex,
      });
    }
  }
  return [...bestByPair.values()].sort((p, q) => q.correlation - p.correlation).slice(0, MAX_TRANSFERS);
}

export interface TransferStats {
  correlation: number;
  /** Per-step effect-size agreement (0-1) over consecutive shared substituents, or null. */
  foldMatch: number | null;
  /** The follower core's predicted (virtual) analog the transfer lets you fill in, or null. */
  benefiting: {side: 'a' | 'b', substIndex: number} | null;
}

function observed(cell: SarMatrixCell): number | null {
  return cell.kind === 'real' && cell.value !== null ? cell.value : null;
}

/**
 * Detailed statistics for a transfer: the correlation, the fold-change match (per-step effect-size
 * agreement between the two cores), and the "benefiting" virtual cell — the untested analog in one
 * core that the transferred rule predicts (preferring the follower core `b`).
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

  let benefiting: {side: 'a' | 'b', substIndex: number} | null = null;
  let best = higherIsBetter ? -Infinity : Infinity;
  for (const [side, cells] of [['b', transfer.bCells], ['a', transfer.aCells]] as const) {
    for (let k = 0; k < cells.length; k++) {
      const cell = cells[k];
      if (cell.kind !== 'virtual' || cell.value === null)
        continue;
      if (higherIsBetter ? cell.value > best : cell.value < best) {
        best = cell.value;
        benefiting = {side, substIndex: k};
      }
    }
    if (benefiting)
      break; // prefer the follower core's analogs
  }
  return {correlation: transfer.correlation, foldMatch, benefiting};
}
