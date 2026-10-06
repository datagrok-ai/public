import {fitAdditiveModel} from './sar-matrix-assemble';
import {SarMatrix, SarMatrixCell, SarMatrixCellKind} from './sar-matrix-types';

/** Below this many cross-validatable observed cells a leave-one-out estimate is too noisy to report. */
const MIN_CV_POINTS = 4;

/** One cross-validated cell and its raw signed `observed − predicted`. Signed rather than
 *  direction-adjusted: the caller owns the activity direction, and folding it in here would put the
 *  same convention in two modules. */
type Residual = {ri: number, ci: number, residual: number};

type Confidence = {r2: number, rmse: number, n: number, total: number,
  hi: Residual | null, lo: Residual | null};

/** Leave-one-out prediction of one observed cell. The held-out cell is blanked so the refit can't see
 *  it; returns null when its row or column has no other observation left. Every column of a matrix
 *  carries the one position the matrix varies, so the slice is the whole grid and this validates
 *  exactly the model that fills it. */
function looPredictSlice(cells: SarMatrixCell[][], rowCount: number, colIdxs: number[],
  targetRow: number, targetK: number): number | null {
  const slice = cells.map((row, ri) => colIdxs.map((ci, k) =>
    (ri === targetRow && k === targetK) ? {...row[ci], kind: 'empty' as SarMatrixCellKind} : row[ci]));
  const predicted = fitAdditiveModel(slice, rowCount, colIdxs.length)(targetRow, targetK);
  return predicted ? predicted.value : null;
}

/**
 * Leave-one-out cross-validated quality of the Free-Wilson fit that fills the matrix. R² near 1 means
 * the additive assumption holds and the virtual predictions are trustworthy; near 0 or negative means
 * substituent effects are non-additive here.
 *
 * Out-of-sample deliberately: an in-sample fit would let a cliff pull the model toward itself and hide
 * its own deviation. Returns null when too few cells are cross-validatable.
 */
export function computeMatrixConfidence(matrix: SarMatrix): Confidence | null {
  const rowCount = matrix.rows.length;
  const pairs: {observed: number, predicted: number}[] = [];
  let total = 0;
  let hi: Residual | null = null;
  let lo: Residual | null = null;
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
        if (predicted === null)
          continue;
        pairs.push({observed: cell.value, predicted});
        const residual = cell.value - predicted;
        if (hi === null || residual > hi.residual)
          hi = {ri, ci: colIdxs[k], residual};
        if (lo === null || residual < lo.residual)
          lo = {ri, ci: colIdxs[k], residual};
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
    return null; // all observed values identical — R² is undefined

  return {r2: 1 - ssRes / ssTot, rmse: Math.sqrt(ssRes / pairs.length), n: pairs.length, total, hi, lo};
}
