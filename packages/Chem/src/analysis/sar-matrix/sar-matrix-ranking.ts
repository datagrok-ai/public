import {SarMatrix} from './sar-matrix-types';

export enum SarRankScheme {
  Potency = 'Potent compounds',
  Discontinuity = 'SAR discontinuity',
  Preferred = 'Preferred substituent',
}

function rowValues(matrix: SarMatrix, rowIdx: number): number[] {
  const values: number[] = [];
  for (const cell of matrix.cells[rowIdx]) {
    if (cell.kind === 'real' && cell.value !== null)
      values.push(cell.value);
  }
  return values;
}

/** Best observed potency as a higher-is-better score: max scaled activity when higher is more
 *  potent, else min. */
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

/** Best mean-potency of any single column — a substituent good across cores. */
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

/** Score every matrix under all schemes and sort by the chosen one. Every score is higher-is-better,
 *  so the top of the list is always best. */
export function rankMatrices(matrices: SarMatrix[], scheme: SarRankScheme,
  higherIsBetter: boolean = true): SarMatrix[] {
  for (const matrix of matrices) {
    matrix.scores[SarRankScheme.Potency] = potencyScore(matrix, higherIsBetter);
    matrix.scores[SarRankScheme.Discontinuity] = discontinuityScore(matrix);
    matrix.scores[SarRankScheme.Preferred] = preferredScore(matrix, higherIsBetter);
  }
  return [...matrices].sort((a, b) => (b.scores[scheme] ?? 0) - (a.scores[scheme] ?? 0));
}

function observedMolecules(matrix: SarMatrix): Set<number> {
  const mols = new Set<number>();
  for (const row of matrix.cells) {
    for (const cell of row) {
      if (cell.kind === 'real' && cell.molIdx !== null)
        mols.add(cell.molIdx);
    }
  }
  return mols;
}

/**
 * Nest matrices by compound containment: a matrix whose compounds are all present in a larger one
 * becomes its child. The parent chosen is the smallest strict superset, so a chain nests one level at
 * a time instead of collapsing onto its root. For identical compound sets the earlier index keeps the
 * parent role, making the result depend on input order alone and not on scoring.
 * @returns Parent index per matrix, `-1` for a root.
 */
export function nestByContainment(matrices: SarMatrix[]): number[] {
  const mols = matrices.map(observedMolecules);
  const parent = new Array<number>(matrices.length).fill(-1);
  for (let i = 0; i < matrices.length; i++) {
    for (let j = 0; j < matrices.length; j++) {
      if (i === j || mols[j].size < mols[i].size)
        continue;
      // Equal sets: only the earlier index may parent the later one, so the two can't adopt each other.
      if (mols[j].size === mols[i].size && j > i)
        continue;
      let contained = true;
      for (const m of mols[i]) {
        if (!mols[j].has(m)) {
          contained = false;
          break;
        }
      }
      if (contained && (parent[i] < 0 || mols[parent[i]].size > mols[j].size))
        parent[i] = j;
    }
  }
  return parent;
}
