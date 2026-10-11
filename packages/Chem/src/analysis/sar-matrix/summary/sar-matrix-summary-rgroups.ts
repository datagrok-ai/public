/* Which substituent comes first, and which reliably last, at each position across the series whose
   additive fit holds. Pure functions of their arguments. */
import {AdditiveFit} from '../build/sar-matrix-assemble';
import {median} from '../build/sar-matrix-decompose';
import {SarMatrix} from '../sar-matrix-types';
import {bestMeasured, LOSER_ROWS, RGROUP_MIN_SERIES, RGroupAcc, RGroupRow, RGroupWin, SUM_ROWS,
  SWAP_MIN_SERIES} from './sar-matrix-summary-types';

/** The substituents the fit puts first and last at this matrix's position, once per lineage, plus every
 *  substituent tried. The fit, not the column mean, which confounds a group with its cores; the
 *  magnitude is against the series' own reference, so it pools across matrices. */
export function recordRGroupExtremes(acc: RGroupAcc, matrix: SarMatrix, fit: AdditiveFit, dir: number,
  root: string, fitHolds: boolean): void {
  // Only a converged fit may be compared with another matrix's, so a series whose fit stopped short
  // contributes to no pool.
  if (!fit.converged)
    return;
  const reference = matrix.refValues[matrix.positions[0] ?? ''];
  const found = reference ? matrix.columns.findIndex((c) => c.substSmiles === reference) : -1;
  // A reference measured once is one residual, too noisy to compare against.
  const refCi = found >= 0 && fit.colN[found] >= 2 ? found : -1;
  let bestCi = -1;
  let worstCi = -1;
  for (let c = 0; c < matrix.columns.length; c++) {
    const subst = matrix.columns[c].substSmiles;
    // A column measured once has a fitted effect made of one residual.
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

/** Record one column as this lineage's extreme. `keep` is +1 for winners, −1 for losers. */
function pushExtreme(pool: Map<string, RGroupWin[]>, matrix: SarMatrix, fit: AdditiveFit, dir: number,
  root: string, fitHolds: boolean, ci: number, refCi: number, keep: number): void {
  // The column's most potent measured cell, so the row lands on a compound.
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
  else if (strongerWin(win, held[at], keep))
    held[at] = win;
}

/** Orders two occurrences of one substituent. On the difference to each matrix's reference, never the
 *  fitted level, whose zero differs between matrices; without a reference, on cell count. */
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

/** One pool ordered strongest first; `keep` is +1 for wins, −1 for losses. */
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
  // Rows with no comparator sort last in either pool: they still rank on how often the group placed.
  const strength = (row: RGroupRow): number =>
    row.magnitude === null ? Number.NEGATIVE_INFINITY : keep * row.magnitude;
  rows.sort((a, b) => b.k - a.k || strength(b) - strength(a) || (a.subst < b.subst ? -1 : 1));
  return rows;
}
