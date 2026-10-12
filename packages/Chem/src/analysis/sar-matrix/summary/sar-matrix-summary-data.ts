import {_package} from '../../../package';
import {fitAdditiveFromTriples} from '../build/sar-matrix-assemble';
import {median} from '../build/sar-matrix-decompose';
import {fitRoleEffects} from './sar-matrix-role-fit';
import {logSarTime, SarMatrix} from '../sar-matrix-types';
import {MatrixCellRef} from '../sar-matrix-ui-common';
import {poolRoleSwaps, poolSwaps, rankSwaps} from './sar-matrix-summary-swaps';
import {rankRGroups, recordRGroupExtremes} from './sar-matrix-summary-rgroups';
import {ANALOG_LIST_MAX, BEST_FIT_MIN_N, MIN_LISTED_SUPPORT, MIN_SUPPORT, RGroupAcc, RowCell, SeriesStat,
  SERIES_FIT_SWEEPS, StartRow, SummaryData, SummaryHost, SummaryRow, supportOf, SWAP_ROW_CAP, SwapPool, TopList,
  TRUST_R2} from './sar-matrix-summary-types';

export * from './sar-matrix-summary-types';

/** Ties end at `matrix.id`: matrices arrive in worker order, which varies run to run. */
function finerSeries(a: SummaryRow, b: SummaryRow): boolean {
  return a.matrix.level !== b.matrix.level ? a.matrix.level < b.matrix.level : a.matrix.id < b.matrix.id;
}

function betterSupported(a: SummaryRow, b: SummaryRow): boolean {
  const sa = supportOf(a);
  const sb = supportOf(b);
  return sa !== sb ? sa > sb : finerSeries(a, b);
}

function betterEvidenced(a: SummaryRow, b: SummaryRow): boolean {
  const sa = supportOf(a);
  const sb = supportOf(b);
  if (sa !== sb)
    return sa > sb;
  const ra = a.matrix.confidence?.r2 ?? Number.NEGATIVE_INFINITY;
  const rb = b.matrix.confidence?.r2 ?? Number.NEGATIVE_INFINITY;
  return ra !== rb ? ra > rb : a.matrix.id < b.matrix.id;
}

export class SummaryCollector {
  constructor(private readonly host: SummaryHost, private readonly tierFilter: number | null) {}

  collect(): SummaryData {
    const t0 = performance.now();
    const host = this.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const log = host.activityIsLog;
    const roots = host.matrixRoots;
    const tiers = host.matrixTiers;
    let ceiling = -Infinity;
    for (const matrix of host.matrices) {
      for (const row of matrix.cells) {
        for (const cell of row) {
          if (cell.kind === 'real' && cell.value !== null)
            ceiling = Math.max(ceiling, dir * cell.value);
        }
      }
    }

    const compounds = new Set<number>();
    const untested = new Set<number>();
    // Measured and virtual SMILES are canonicalised differently, so this catches matches but proves no novelty.
    const measuredStructures = new Set<string>();
    const test = new TopList(betterSupported);
    const analogs = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const analogsThin = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const analogsAny = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const swaps = new Map<string, SwapPool>();
    const acc: RGroupAcc = {wins: new Map(), losses: new Map(), tried: new Map()};
    const series: SeriesStat[] = [];
    let lowR2 = 0;

    let measuredCells = 0;
    let trustedNoStructure = 0;
    let withheldThinSupport = 0;
    let withheldFitFails = 0;
    let withheldUnchecked = 0;
    let unchecked = 0;
    let fitHolds = 0;
    let nonConverged = 0;
    let sampledRows = 0;
    let swapCandidates = 0;

    const obsRow: number[] = [];
    const obsCol: number[] = [];
    const obsVal: number[] = [];
    const mols = new Set<number>();
    const rowCells: RowCell[] = [];
    // An analog's gain needs the matrix's best measured cell, final only after its rows are walked.
    const riBuf: number[] = [];
    const ciBuf: number[] = [];
    const riAlt: number[] = [];
    const ciAlt: number[] = [];

    const roleNames = host.axisRole === null ? [] : [host.coreRole!, ...host.roleColumns];
    const roleKind = roleNames.map((name) => name === host.coreRole ? 0 : name === host.axisRole ? 1 : 2);
    const roleValues: string[][] = roleNames.map(() => []);
    const roleActivity: number[] = [];
    const roleMol: number[] = [];
    const roleBest = roleNames.map(() => new Map<string, MatrixCellRef & {value: number}>());

    // Coverage is of the whole analysis: a tier filter must not move compounds out of "measured".
    if (this.tierFilter !== null) {
      for (const matrix of host.matrices) {
        for (const row of matrix.cells) {
          for (const cell of row) {
            if (cell.molIdx !== null && cell.kind === 'real' && cell.value !== null)
              compounds.add(cell.molIdx);
            else if (cell.molIdx !== null && cell.kind === 'unmeasured')
              untested.add(cell.molIdx);
          }
        }
      }
    }

    for (let mi = 0; mi < host.matrices.length; mi++) {
      if (this.tierFilter !== null && tiers[mi] !== this.tierFilter)
        continue;
      const matrix = host.matrices[mi];
      const nRows = matrix.rows.length;
      const nCols = matrix.columns.length;
      const r2 = matrix.confidence?.r2 ?? null;
      const holds = r2 !== null && r2 >= TRUST_R2;

      obsRow.length = 0;
      obsCol.length = 0;
      obsVal.length = 0;
      riBuf.length = 0;
      ciBuf.length = 0;
      riAlt.length = 0;
      ciAlt.length = 0;
      mols.clear();
      let realCells = 0;
      let impossibleCells = 0;
      let trusted = 0;
      let lo = Infinity;
      let hi = -Infinity;
      let best: {ri: number, ci: number, value: number} | null = null;
      let bestVirtual: {ri: number, ci: number, value: number} | null = null;
      let virtualCells = 0;

      for (let ri = 0; ri < nRows; ri++) {
        rowCells.length = 0;
        for (let ci = 0; ci < nCols; ci++) {
          const cell = matrix.cells[ri][ci];
          const value = cell.value;
          if (cell.kind === 'real' && value !== null && cell.molIdx !== null) {
            compounds.add(cell.molIdx);
            mols.add(cell.molIdx);
            measuredCells++;
            realCells++;
            obsRow.push(ri);
            obsCol.push(ci);
            obsVal.push(value);
            rowCells.push({ci, value, molIdx: cell.molIdx});
            if (roleNames.length > 0) {
              for (let r = 0; r < roleNames.length; r++) {
                const level = roleKind[r] === 0 ? matrix.rows[ri].coreSmiles :
                  roleKind[r] === 1 ? matrix.columns[ci].substSmiles :
                    matrix.rows[ri].foldedValues[roleNames[r]] ?? '';
                roleValues[r].push(level);
                const held = roleBest[r].get(level);
                if (held === undefined || dir * value > dir * held.value)
                  roleBest[r].set(level, {matrix, ri, ci, value});
              }
              roleActivity.push(value);
              roleMol.push(cell.molIdx);
            }
            if (value < lo)
              lo = value;
            if (value > hi)
              hi = value;
            if (best === null || dir * value > dir * best.value)
              best = {ri, ci, value};
            if (cell.smiles !== null)
              measuredStructures.add(cell.smiles);
          } else if (cell.kind === 'unmeasured' && cell.molIdx !== null) {
            untested.add(cell.molIdx);
            if (cell.smiles !== null)
              measuredStructures.add(cell.smiles);
            if (value !== null && (cell.support ?? 0) >= MIN_LISTED_SUPPORT)
              test.offer(String(cell.molIdx), dir * value, matrix, ri, ci);
          } else if (cell.kind === 'virtual' && value !== null) {
            const support = cell.support ?? 0;
            virtualCells++;
            const passes = support >= MIN_SUPPORT && holds;
            if (passes && (bestVirtual === null || dir * value > dir * bestVirtual.value))
              bestVirtual = {ri, ci, value};
            // A structureless prediction is a settings fix (an unfilled attachment point), not withheld.
            if (cell.smiles === null) {
              if (passes)
                trustedNoStructure++;
            } else if (passes) {
              trusted++;
              if (!measuredStructures.has(cell.smiles)) {
                riBuf.push(ri);
                ciBuf.push(ci);
              }
            } else {
              if (support >= MIN_LISTED_SUPPORT && !measuredStructures.has(cell.smiles)) {
                riAlt.push(ri);
                ciAlt.push(ci);
              }
              if (support < MIN_SUPPORT)
                withheldThinSupport++;
              else if (r2 === null)
                withheldUnchecked++;
              else
                withheldFitFails++;
            }
          } else if (cell.kind === 'impossible')
            impossibleCells++;
        }
        if (rowCells.length >= 2) {
          swapCandidates++;
          if (rowCells.length > SWAP_ROW_CAP)
            sampledRows++;
          poolSwaps(swaps, matrix, ri, rowCells, dir, log, roots[mi]);
        }
      }

      const fit = fitAdditiveFromTriples(obsRow, obsCol, obsVal, nRows, nCols, SERIES_FIT_SWEEPS);
      if (!fit.converged)
        nonConverged++;
      // A fit that stopped short is shrunk by an amount of its own, so its gains do not compare.
      const conf = matrix.confidence ?? null;
      if (fit.converged && best !== null && conf !== null) {
        // An error from a handful of cells is no scale to divide by: it is tiny and would top every ranking.
        const thin = conf.n < BEST_FIT_MIN_N || conf.rmse <= 0;
        for (let k = 0; k < riBuf.length; k++) {
          const cell = matrix.cells[riBuf[k]][ciBuf[k]];
          const gain = dir * (cell.value! - best.value);
          if (thin) {
            if (gain > 0)
              analogsThin.offer(cell.smiles!, gain, matrix, riBuf[k], ciBuf[k]);
          } else if (gain >= conf.rmse)
            analogs.offer(cell.smiles!, gain / conf.rmse, matrix, riBuf[k], ciBuf[k]);
        }
      }
      if (best !== null) {
        for (let k = 0; k < riAlt.length; k++) {
          const cell = matrix.cells[riAlt[k]][ciAlt[k]];
          const gain = dir * (cell.value! - best.value);
          // An untrusted fit predicting past every measured value is extrapolating, not proposing.
          if (gain > 0 && dir * cell.value! <= ceiling)
            analogsAny.offer(cell.smiles!, gain, matrix, riAlt[k], ciAlt[k]);
        }
      }
      const stat: SeriesStat = {
        matrix, root: roots[mi], tier: tiers[mi], cpd: mols.size, realCells, impossibleCells,
        totalCells: nRows * nCols, lo: realCells ? lo : 0, hi: realCells ? hi : 0, best, bestVirtual,
        typical: fit.grandMean, trusted, virtualCells,
        bestRow: this.bestRow(matrix, fit, dir),
      };
      series.push(stat);
      recordRGroupExtremes(acc, matrix, fit, dir, roots[mi], holds);

      if (r2 === null)
        unchecked++;
      else if (r2 < TRUST_R2)
        lowR2++;
      else
        fitHolds++;
    }

    if (nonConverged > 0) {
      _package.logger.warning(`SAR Matrix | the additive fit of ${nonConverged} series did not reach ` +
        'its tolerance; they are left out of the pooled R-group comparison');
    }
    if (sampledRows > 0) {
      _package.logger.warning(`SAR Matrix | ${sampledRows} rows carry more than ${SWAP_ROW_CAP} measured ` +
        'cells; the swap pool keeps their most and least potent halves, so mid-range pairs are dropped');
    }

    // A molecule first measured in a later matrix is only caught after the walk.
    const dropHeld = (key: string): boolean => measuredStructures.has(key);
    analogs.prune(dropHeld);
    analogsThin.prune(dropHeld);
    analogsAny.prune(dropHeld);

    const measured = series.filter((stat) => stat.realCells > 0);
    const confs = series.flatMap((stat) => stat.matrix.confidence ?? []);
    const roleFit = fitRoleEffects({names: roleNames, values: roleValues, activity: roleActivity,
      molIdx: roleMol, minSupport: MIN_SUPPORT, higherIsBetter: host.higherIsBetter});
    const data: SummaryData = {
      compounds: compounds.size,
      untested: untested.size,
      measuredCells,
      minObserved: measured.length ? Math.min(...measured.map((stat) => stat.lo)) : null,
      maxObserved: measured.length ? Math.max(...measured.map((stat) => stat.hi)) : null,
      trustedNoStructure,
      modelError: confs.length ? median(confs.map((conf) => conf.rmse)) : null,
      fitHolds, unchecked, lowR2,
      axisRole: host.axisRole,
      coreRole: host.coreRole,
      coresAreSeries: host.coresAreSeries,
      roleFit,
      roleBest: new Map(roleNames.map((name, r) => [name, roleBest[r]])),
      nonConverged,
      series,
      startHere: this.mergeStartHere(series, confs.map((conf) => conf.n), dir),
      swaps: rankSwaps(swaps),
      swapsByRole: poolRoleSwaps(roleNames, roleValues, roleActivity, roleMol, dir, log),
      swapCandidates,
      ...rankRGroups(acc),
      test, analogs, analogsThin, analogsAny,
      withheldThinSupport, withheldFitFails, withheldUnchecked,
      analogOverflow: analogs.overflow + analogsThin.overflow,
    };
    logSarTime('summary', t0);
    return data;
  }

  private bestRow(matrix: SarMatrix, fit: ReturnType<typeof fitAdditiveFromTriples>, dir: number):
    {ri: number, effect: number, n: number} | null {
    let best: {ri: number, effect: number, n: number} | null = null;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      if (fit.rowN[ri] >= MIN_SUPPORT && (best === null || dir * fit.rowEffect[ri] > dir * best.effect))
        best = {ri, effect: fit.rowEffect[ri], n: fit.rowN[ri]};
    }
    return best;
  }

  /** Merged by lineage root: one series can win several reasons, and a lineage is one place to look. */
  private mergeStartHere(series: SeriesStat[], confN: number[], dir: number): StartRow[] {
    if (series.length === 0)
      return [];
    // Floored well above the CV minimum: the best R² of tiny fits is just the noisiest one.
    const nFloor = Math.max(BEST_FIT_MIN_N, median(confN));
    const pick = (score: (s: SeriesStat) => number | null): SeriesStat | null => {
      let winner: SeriesStat | null = null;
      let bestScore = 0;
      for (const stat of series) {
        const value = score(stat);
        if (value === null)
          continue;
        if (winner === null || value > bestScore ||
          (value === bestScore && stat.matrix.id < winner.matrix.id)) {
          winner = stat;
          bestScore = value;
        }
      }
      return winner;
    };
    // Index-aligned with REASON_WORDS.
    const slots: (SeriesStat | null)[] = [
      pick((s) => s.cpd || null),
      pick((s) => s.realCells ? s.hi - s.lo : null),
      pick((s) => s.best === null ? null : dir * s.best.value),
      pick((s) => s.trusted || null),
      pick((s) => s.matrix.confidence && s.matrix.confidence.n >= nFloor ? s.matrix.confidence.r2 : null),
    ];
    const byRoot = new Map<string, StartRow>();
    for (let slot = 0; slot < slots.length; slot++) {
      const stat = slots[slot];
      if (stat === null)
        continue;
      const held = byRoot.get(stat.root);
      if (held === undefined)
        byRoot.set(stat.root, {primary: stat, reasons: [slot]});
      else
        held.reasons.push(slot);
    }
    return [...byRoot.values()];
  }
}
