import {fitAdditiveModel} from './sar-matrix-assemble';
import {SarMatrix, SarMatrixCell, SarMatrixCellKind} from './sar-matrix-types';

/** Below this many observed cells a leave-one-out estimate is too noisy to report. */
const MIN_CV_POINTS = 4;

/** Predict one observed cell from a Free-Wilson fit that excludes it, or null when its row or column
 *  has no other observation left to fit against. */
function looPredict(cells: SarMatrixCell[][], rowCount: number, columnCount: number,
  targetRow: number, targetCol: number): number | null {
  // Shallow-clone the grid with the held-out cell turned 'empty' so the refit can't see it.
  const heldOut = cells.map((row, ri) => row.map((cell, ci) =>
    (ri === targetRow && ci === targetCol) ? {...cell, kind: 'empty' as SarMatrixCellKind} : cell));
  const predicted = fitAdditiveModel(heldOut, rowCount, columnCount)(targetRow, targetCol);
  return predicted ? predicted.value : null;
}

/**
 * Leave-one-out cross-validated quality of the two-way additive (Free-Wilson) fit over a matrix's
 * observed cells: refit the model with each observed cell held out, predict it, and score the
 * held-out residuals. R² compares those residuals to just predicting the mean — a value near 1 means
 * the additive assumption holds and the virtual predictions are trustworthy; a value near 0 (or
 * negative) means substituent effects are non-additive here and the predictions should be doubted.
 * Returns null when there are too few observations, or too few of them are cross-validatable.
 */
export function computeMatrixConfidence(matrix: SarMatrix): {r2: number, rmse: number, n: number} | null {
  const rowCount = matrix.rows.length;
  const columnCount = matrix.columns.length;
  const pairs: {observed: number, predicted: number}[] = [];
  for (let ri = 0; ri < rowCount; ri++) {
    for (let ci = 0; ci < columnCount; ci++) {
      const cell = matrix.cells[ri][ci];
      if (cell.kind !== 'real' || cell.value === null)
        continue;
      const predicted = looPredict(matrix.cells, rowCount, columnCount, ri, ci);
      if (predicted !== null)
        pairs.push({observed: cell.value, predicted});
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

  return {r2: 1 - ssRes / ssTot, rmse: Math.sqrt(ssRes / pairs.length), n: pairs.length};
}
