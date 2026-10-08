/* Which substituent comes first, and which reliably last, at each position across the series whose
   additive fit holds. Pure functions of their arguments. */
import {AdditiveFit} from '../build/sar-matrix-assemble';
import {median} from '../build/sar-matrix-decompose';
import {SarMatrix} from '../sar-matrix-types';
import {bestMeasured, LOSER_ROWS, RGROUP_MIN_SERIES, RGroupAcc, RGroupRow, RGroupWin, StripCols, SUM_ROWS,
  SWAP_MIN_SERIES} from './sar-matrix-summary-types';

/**
 * The substituents the fitted model puts first and last at this matrix's position, recorded once per
 * lineage, plus every substituent this series tried and — where a strip will draw — its columns.
 *
 * The fit, not the raw column mean: a column mean confounds a substituent with whichever cores it
 * happened to be made on, which is the normal case in a modular library. The magnitude is a
 * within-series difference against that series' own reference substituent, so it survives the fact
 * that two matrices centre their effects over different substituent menus. An ordinal survives it
 * too: adding a constant to every `colEffect` of a matrix does not reorder them, which is what lets
 * the last place pool the same way the first does.
 */
export function recordRGroupExtremes(acc: RGroupAcc, matrix: SarMatrix, fit: AdditiveFit, dir: number,
  root: string, fitHolds: boolean, strip: boolean): StripCols | null {
  // Only a converged fit may be compared with another matrix's, so a series whose fit stopped short
  // contributes to no pool and takes no strip slot either: an empty strip entry would draw the "never
  // tried" mark over a series that tried the group and was dropped.
  if (!fit.converged)
    return null;
  const reference = matrix.refValues[matrix.positions[0] ?? ''];
  const found = reference ? matrix.columns.findIndex((c) => c.substSmiles === reference) : -1;
  // A reference measured once has a fitted effect made of one residual, so a difference against it
  // is noise wearing a comparator's name.
  const refCi = found >= 0 && fit.colN[found] >= 2 ? found : -1;
  const cols: StripCols | null = strip ? new Map() : null;
  let bestCi = -1;
  let worstCi = -1;
  for (let c = 0; c < matrix.columns.length; c++) {
    const subst = matrix.columns[c].substSmiles;
    // A column measured once has a fitted effect made of one residual. The strip is filled under the
    // same test, so a square and the coverage sentence under it can never describe different sets of
    // columns.
    if (fit.colN[c] < 2)
      continue;
    if (cols !== null) {
      // 0 where the column IS the reference — the honest difference against itself, which is not the
      // "no comparator at all" the null stands for.
      cols.set(subst, {ci: c, n: fit.colN[c],
        refDelta: refCi < 0 ? null : refCi === c ? 0 :
          dir * (fit.colEffect[c] - fit.colEffect[refCi])});
    }
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
  return cols;
}

/** Record one column as this lineage's extreme, keeping the occurrence with the stronger
 *  within-series margin. `keep` is +1 for the winner pool and −1 for the loser pool, so one dedup
 *  serves both. */
function pushExtreme(pool: Map<string, RGroupWin[]>, matrix: SarMatrix, fit: AdditiveFit, dir: number,
  root: string, fitHolds: boolean, ci: number, refCi: number, keep: number): void {
  // The most potent measured cell of the column, so the row lands on a compound rather than on a
  // hole the reader has to hunt through.
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
  // One entry per lineage: a substituent extreme in a matrix and again in the tier that folded it has
  // been confirmed once, over one set of compounds, and must not read as two series.
  const at = held.findIndex((other) => other.root === root);
  if (at < 0)
    held.push(win);
  else if (strongerWin(win, held[at], keep))
    held[at] = win;
}

/**
 * Orders two occurrences of one substituent, best of `keep`'s pool first.
 *
 * Never on the fitted level itself: `colEffect` is centred over whichever substituents its own
 * matrix happened to make, so two matrices' levels have different zeros and picking a maximum over
 * them is arithmetic on two different gauges. The difference against each matrix's own reference has
 * no gauge in it; where a matrix records no reference, the measured cells behind the column decide.
 */
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
    // No thin tier for losers: a loss over two lineages is not a rule, and striking a building block
    // from the next enumeration is acted on without being checked.
    rgroupLosers: losers.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, LOSER_ROWS),
  };
}

/** One pool ordered best-first, `keep` being +1 for wins and −1 for losses so the strongest of each
 *  pool comes out on top of its own list. */
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
