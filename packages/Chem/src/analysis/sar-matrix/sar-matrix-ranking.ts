import {SarMatrix} from './sar-matrix-types';

export enum SarRankScheme {
  Potency = 'Potent compounds',
  Discontinuity = 'SAR discontinuity',
  Preferred = 'Preferred substituent',
}

/** Observed activity values of one matrix row (real cells only). */
function rowValues(matrix: SarMatrix, rowIdx: number): number[] {
  const values: number[] = [];
  for (const cell of matrix.cells[rowIdx]) {
    if (cell.kind === 'real' && cell.value !== null)
      values.push(cell.value);
  }
  return values;
}

/** Best observed potency, as a higher-is-better score: the max scaled activity when higher means
 *  more potent (`-lg`), else the min (raw `none`/`lg`, where lower nM is more potent). */
function potencyScore(matrix: SarMatrix, higherIsBetter: boolean): number {
  if (!matrix.realCount)
    return 0;
  return higherIsBetter ? matrix.maxActivity : -matrix.minActivity;
}

/** Largest activity spread within any single row — a proxy for SAR discontinuity. */
function discontinuityScore(matrix: SarMatrix): number {
  let best = 0;
  for (let ri = 0; ri < matrix.rows.length; ri++) {
    const values = rowValues(matrix, ri);
    if (values.length < 2)
      continue;
    best = Math.max(best, Math.max(...values) - Math.min(...values));
  }
  return best;
}

/** Best mean-potency of any single column — a substituent good across cores. Scored higher-is-
 *  better, so for a lower-is-better activity the most negative column mean wins. */
function preferredScore(matrix: SarMatrix, higherIsBetter: boolean): number {
  let best = -Infinity;
  for (let ci = 0; ci < matrix.columns.length; ci++) {
    let sum = 0;
    let n = 0;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      const cell = matrix.cells[ri][ci];
      if (cell.kind === 'real' && cell.value !== null) {
        sum += cell.value;
        n++;
      }
    }
    if (n >= 2)
      best = Math.max(best, higherIsBetter ? sum / n : -(sum / n));
  }
  return best === -Infinity ? 0 : best;
}

/**
 * Step 8 — score every matrix under all schemes and sort by the chosen one. (SAR transfer is its own
 * navigator section, computed globally across matrices in `sar-matrix-transfer.ts`, not a rank mode.)
 *
 * @param higherIsBetter Whether a higher scaled activity means a more potent compound. Potency and
 *   preferred-substituent scores flip with it (the discontinuity spread is direction-agnostic). Every
 *   score is higher-is-better so the top of the list is always best.
 */
export function rankMatrices(matrices: SarMatrix[], scheme: SarRankScheme,
  higherIsBetter: boolean = true): SarMatrix[] {
  for (const matrix of matrices) {
    matrix.scores[SarRankScheme.Potency] = potencyScore(matrix, higherIsBetter);
    matrix.scores[SarRankScheme.Discontinuity] = discontinuityScore(matrix);
    matrix.scores[SarRankScheme.Preferred] = preferredScore(matrix, higherIsBetter);
  }
  return [...matrices].sort((a, b) => (b.scores[scheme] ?? 0) - (a.scores[scheme] ?? 0));
}
