import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

import {Fingerprint, rdKitFingerprintToBitArray} from '../../utils/chem-common';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {SarMatrix, SarMatrixCell, logSarTime} from './sar-matrix-types';

/**
 * SAR-transfer detection: a transfer rank-correlates the potencies of shared substituents across
 * related-core rows. Same-matrix pairs (two rows = related-core analog series) are the canonical case
 * (Gupta-Ostermann & Bajorath, F1000Research 2014); cross-matrix pairs are compared too. The
 * `similarity` floor is optional (0 = off) — relatedness comes from the shared substituents, not a
 * whole-molecule Tanimoto cutoff.
 */

const MIN_COMMON = 3;
const MIN_CORRELATION = 0.7;
const MAX_TRANSFERS = 30;
const MAX_PREDICTED = 6;
export const DEFAULT_TRANSFER_SIMILARITY = 0;

export interface TransferSide {
  matrixIndex: number;
  rowIndex: number;
  position: string;
}

export interface VirtualAnalog {
  col: number;
  substSmiles: string;
  value: number;
}

export interface Transfer {
  a: TransferSide;
  b: TransferSide;
  substituents: string[];
  aCols: number[];
  bCols: number[];
  aCells: SarMatrixCell[];
  bCells: SarMatrixCell[];
  similarity: number;
  correlation: number;
  aVirtual: VirtualAnalog[];
  bVirtual: VirtualAnalog[];
  // Kept apart from the matched pairs: predicted analogs have no second observation, so they take no
  // part in the correlation or the fold-change match.
  predictedSubstituents: string[];
  predictedACols: number[];
  predictedBCols: number[];
}

interface SeriesPoint {
  col: number;
  cell: SarMatrixCell;
  value: number;
  substSmiles: string;
  fp: number;
  smiles: string;
}

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

export function spearman(xs: number[], ys: number[]): number | null {
  return pearson(averageRanks(xs), averageRanks(ys));
}

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

function pairKey(a: TransferSide, b: TransferSide): string {
  const x = `${a.matrixIndex}:${a.rowIndex}`;
  const y = `${b.matrixIndex}:${b.rowIndex}`;
  return x < y ? `${x}|${y}` : `${y}|${x}`;
}

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
      // Same-matrix pairs (rows = related-core analog series) are the canonical transfer; only like
      // positions are comparable.
      if (s1.position !== s2.position)
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
        continue;
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
  // Strongest first; among equal correlations the longer overlap is stronger evidence (a 3-point ρ=1
  // can be coincidence, a 7-point one cannot), then the more alike scaffolds.
  return [...bestByPair.values()]
    .sort((p, q) => q.correlation - p.correlation ||
      q.substituents.length - p.substituents.length || q.similarity - p.similarity)
    .slice(0, MAX_TRANSFERS);
}

export interface TransferStats {
  correlation: number;
  foldMatch: number | null;
  benefiting: {side: 'a' | 'b', substSmiles: string, value: number} | null;
}

function observed(cell: SarMatrixCell): number | null {
  return cell.kind === 'real' && cell.value !== null ? cell.value : null;
}

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
      break;
    }
  }
  return {correlation: transfer.correlation, foldMatch, benefiting};
}
