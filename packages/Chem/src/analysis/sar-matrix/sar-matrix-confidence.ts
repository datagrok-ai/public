import {fitAdditiveModel} from './sar-matrix-assemble';
import {SarMatrix, SarMatrixCell, SarMatrixCellKind} from './sar-matrix-types';

/** Below this many cross-validatable observed cells a leave-one-out estimate is too noisy to report. */
const MIN_CV_POINTS = 4;

/** Cross-validated fit quality: `n` cross-validatable observed cells out of `total`. */
type Confidence = {r2: number, rmse: number, n: number, total: number};

/** Leave-one-out prediction of one observed cell within a single R-position's column slice — the same
 *  per-position additive (Free-Wilson) model assembly uses to fill cells. The held-out cell is blanked
 *  so the refit can't see it; returns null when its row or column has no other observation left in the
 *  slice. Fitting per position (not one global model over every column) means the reference-substituent
 *  "twin" cells of an all-reference compound live in different slices and can't leak into each other. */
function looPredictSlice(cells: SarMatrixCell[][], rowCount: number, colIdxs: number[],
  targetRow: number, targetK: number): number | null {
  const slice = cells.map((row, ri) => colIdxs.map((ci, k) =>
    (ri === targetRow && k === targetK) ? {...row[ci], kind: 'empty' as SarMatrixCellKind} : row[ci]));
  const predicted = fitAdditiveModel(slice, rowCount, colIdxs.length)(targetRow, targetK);
  return predicted ? predicted.value : null;
}

/**
 * Leave-one-out cross-validated quality of the two-way additive (Free-Wilson) fit that actually fills
 * the matrix — validated per R-position, exactly as the cells were predicted. For each observed cell:
 * refit its position's slice with that cell held out, predict it, and score the held-out residuals.
 * R² compares those residuals to just predicting the mean — near 1 means the additive assumption holds
 * and the virtual predictions are trustworthy; near 0 (or negative) means substituent effects are
 * non-additive here and the predictions should be doubted.
 *
 * As a side effect every observed real cell gets its LOO fitted value stored on `cell.fit` (or cleared
 * when it isn't cross-validatable), so the per-cell "non-additive" flag reflects an honest out-of-sample
 * residual instead of an in-sample one that would let a cliff hide its own deviation.
 *
 * Returns null when too few observed cells are cross-validatable. `n` cross-validatable of `total`
 * observed cells lets the caller show how much of the matrix the R² actually covers.
 */
export function computeMatrixConfidence(matrix: SarMatrix): Confidence | null {
  const rowCount = matrix.rows.length;
  const pairs: {observed: number, predicted: number}[] = [];
  let total = 0;
  for (const position of matrix.positions) {
    const colIdxs = matrix.columns
      .map((c, ci) => (c.position === position ? ci : -1)).filter((ci) => ci >= 0);
    for (let ri = 0; ri < rowCount; ri++) {
      for (let k = 0; k < colIdxs.length; k++) {
        const cell = matrix.cells[ri][colIdxs[k]];
        if (cell.kind !== 'real' || cell.value === null)
          continue;
        total++;
        const predicted = looPredictSlice(matrix.cells, rowCount, colIdxs, ri, k);
        if (predicted !== null) {
          pairs.push({observed: cell.value, predicted});
          cell.fit = predicted; // honest out-of-sample fitted value — drives the non-additive flag
        } else
          cell.fit = undefined; // can't cross-validate this cell — show no fit verdict for it
      }
    }
  }
  if (pairs.length < MIN_CV_POINTS)
    return null;

  const meanObserved = pairs.reduce((s, p) => s + p.observed, 0) / pairs.length;
  let ssRes = 0;
  let ssTot = 0;
  for (const {observed, predicted} of pairs) {
    ssRes += (observed - predicted) ** 2;
    ssTot += (observed - meanObserved) ** 2;
  }
  if (ssTot === 0)
    return null; // all observed values identical — R² is undefined, nothing to validate

  return {r2: 1 - ssRes / ssTot, rmse: Math.sqrt(ssRes / pairs.length), n: pairs.length, total};
}
