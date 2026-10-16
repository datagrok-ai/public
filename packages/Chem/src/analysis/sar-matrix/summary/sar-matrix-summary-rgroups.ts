import {AdditiveFit} from '../build/sar-matrix-assemble';
import {median} from '../build/sar-matrix-decompose';
import {SarMatrix, SarMatrixCell} from '../sar-matrix-types';
import {LOSER_ROWS, RGROUP_MIN_SERIES, RGroupAcc, RGroupRow, RGroupWin, SUM_ROWS,
  SWAP_MIN_SERIES} from './sar-matrix-summary-types';

function bestMeasured(cells: SarMatrixCell[], dir: number): number {
  let best = -1;
  cells.forEach((cell, i) => {
    if (cell.kind === 'real' && cell.value !== null && (best < 0 || dir * cell.value > dir * cells[best].value!))
      best = i;
  });
  return best;
}

/** Ranked on the fit, not the column mean, which confounds a group with its cores. */
export function recordRGroupExtremes(acc: RGroupAcc, matrix: SarMatrix, fit: AdditiveFit, dir: number,
  root: string, fitHolds: boolean): void {
  // Only a converged fit may be compared with another matrix's.
  if (!fit.converged)
    return;
  const reference = matrix.refValues[matrix.positions[0] ?? ''];
  const found = reference ? matrix.columns.findIndex((c) => c.substSmiles === reference) : -1;
  // A column measured once has a fitted effect made of one residual, too noisy to compare.
  const refCi = found >= 0 && fit.colN[found] >= 2 ? found : -1;
  let bestCi = -1;
  let worstCi = -1;
  for (let c = 0; c < matrix.columns.length; c++) {
    const subst = matrix.columns[c].substSmiles;
    if (fit.colN[c] < 2)
      continue;
    let attempts = acc.tried.get(subst);
    if (attempts === undefined)
      acc.tried.set(subst, attempts = new Set());
    attempts.add(root);
    if (bestCi < 0 || dir * fit.colEffect[c] > dir * fit.colEffect[bestCi] ||
      (fit.colEffect[c] === fit.colEffect[bestCi] && subst < matrix.columns[bestCi].substSmiles))
      bestCi = c;
    if (worstCi < 0 || dir * fit.colEffect[c] < dir * fit.colEffect[worstCi] ||
      (fit.colEffect[c] === fit.colEffect[worstCi] && subst < matrix.columns[worstCi].substSmiles))
      worstCi = c;
  }
  if (bestCi >= 0)
    pushExtreme(acc.wins, matrix, fit, dir, root, fitHolds, bestCi, refCi, 1);
  // One comparable column is first and last at once, which is not a loss.
  if (worstCi >= 0 && worstCi !== bestCi)
    pushExtreme(acc.losses, matrix, fit, dir, root, fitHolds, worstCi, refCi, -1);
}

/** `keep` is +1 for winners, −1 for losers. */
function pushExtreme(pool: Map<string, RGroupWin[]>, matrix: SarMatrix, fit: AdditiveFit, dir: number,
  root: string, fitHolds: boolean, ci: number, refCi: number, keep: number): void {
  const ri = bestMeasured(matrix.cells.map((cells) => cells[ci]), dir);
  if (ri < 0)
    return;
  const win: RGroupWin = {
    matrix, ri, ci, root, n: fit.colN[ci],
    refDelta: refCi >= 0 && refCi !== ci ? dir * (fit.colEffect[ci] - fit.colEffect[refCi]) : null,
    fitHolds,
  };
  const subst = matrix.columns[ci].substSmiles;
  const held = pool.get(subst);
  if (held === undefined) {
    pool.set(subst, [win]);
    return;
  }
  // One entry per lineage: a folded tier repeats the same compounds.
  const at = held.findIndex((other) => other.root === root);
  if (at < 0)
    held.push(win);
  // A trusted occurrence beats an untrusted one, or the lineage drops out of the trusted count.
  else if (win.fitHolds !== held[at].fitHolds ? win.fitHolds : strongerWin(win, held[at], keep))
    held[at] = win;
}

/** Never on the fitted level, whose zero differs between matrices. */
function strongerWin(a: RGroupWin, b: RGroupWin, keep: number): boolean {
  if (a.refDelta !== null && b.refDelta !== null && a.refDelta !== b.refDelta)
    return keep * a.refDelta > keep * b.refDelta;
  if ((a.refDelta === null) !== (b.refDelta === null))
    return b.refDelta === null;
  return a.n !== b.n ? a.n > b.n : a.matrix.id < b.matrix.id;
}

export function rankRGroups(acc: RGroupAcc):
  {rgroups: RGroupRow[], rgroupsThin: RGroupRow[], rgroupLosers: RGroupRow[]} {
  const winners = rankPool(acc.wins, acc.tried, 1);
  const losers = rankPool(acc.losses, acc.tried, -1);
  return {
    rgroups: winners.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, SUM_ROWS + 2),
    rgroupsThin: winners.filter((row) => row.k < RGROUP_MIN_SERIES).slice(0, SUM_ROWS),
    // No thin tier for losers: "stop making" over two lineages is not a rule.
    rgroupLosers: losers.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, LOSER_ROWS),
  };
}

function rankPool(pool: Map<string, RGroupWin[]>, tried: Map<string, Set<string>>,
  keep: number): RGroupRow[] {
  const rows: RGroupRow[] = [];
  for (const [subst, entries] of pool) {
    const trusted = entries.filter((win) => win.fitHolds);
    if (trusted.length < SWAP_MIN_SERIES)
      continue;
    const deltas = trusted.map((win) => win.refDelta).filter((d): d is number => d !== null);
    trusted.sort((a, b) => strongerWin(a, b, keep) ? -1 : 1);
    rows.push({
      subst, win: trusted[0], k: trusted.length, m: entries.length,
      tried: tried.get(subst)!.size,
      magnitude: deltas.length ? median(deltas) : null,
      lo: deltas.length ? Math.min(...deltas) : 0,
      hi: deltas.length ? Math.max(...deltas) : 0,
      magnitudeSeries: deltas.length,
    });
  }
  const strength = (row: RGroupRow): number =>
    row.magnitude === null ? Number.NEGATIVE_INFINITY : keep * row.magnitude;
  rows.sort((a, b) => b.k - a.k || strength(b) - strength(a) || (a.subst < b.subst ? -1 : 1));
  return rows;
}
