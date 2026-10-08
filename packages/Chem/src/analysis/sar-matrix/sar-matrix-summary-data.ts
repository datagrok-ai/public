/* The Summary tab's computation: one walk over every matrix that produces everything the tab
   shows. Kept apart from the rendering because it touches no DOM: a SummaryCollector is handed the
   viewer and the fold tier, and returns SummaryData. The algorithms it calls and the types it
   produces live in their own modules; the types are re-exported so callers import from here. */
import {_package} from '../../package';
import {fitAdditiveFromTriples} from './sar-matrix-assemble';
import {median} from './sar-matrix-decompose';
import {fitRoleEffects} from './sar-matrix-role-fit';
import {logSarTime, SarMatrix} from './sar-matrix-types';
import {MatrixCellRef} from './sar-matrix-ui-common';
import {largestSwap, poolRoleSwaps, poolSwaps, rankSwaps} from './sar-matrix-summary-swaps';
import {rankRGroups, recordRGroupExtremes} from './sar-matrix-summary-rgroups';
import {ANALOG_LIST_MAX, BEST_FIT_MIN_N, MIN_SUPPORT, RGroupAcc, RowCell, SeriesStat, StartRow, STRIP_SLOTS,
  SummaryData, SummaryHost, SummaryRow, supportOf, SWAP_ROW_CAP, SwapPool, TopList,
  TRUST_R2} from './sar-matrix-summary-types';

export * from './sar-matrix-summary-types';

/** Ties end at `matrix.id`: fragments arrive in worker-completion order, so anything resolved by
 *  index or Map order would make the cards depend on scheduling. */
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

  /**
   * One pass over every cell of every matrix, feeding the totals, the two pools and the per-series
   * statistics.
   *
   * The observed cells are collected as triples on the way past and the additive fit is run from
   * those, so the fit does not walk the grid a second time. The buffers are reused across matrices:
   * one fresh set per matrix would trade the scan cost for GC cost, and nothing caps the SUM of cells across
   * matrices — only each matrix.
   */
  collect(): SummaryData {
    const t0 = performance.now();
    const host = this.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const log = host.activityIsLog;
    const roots = host.matrixRoots;
    const tiers = host.matrixTiers;

    const compounds = new Set<number>();
    const untested = new Set<number>();
    const trustedStructures = new Set<string>();
    // Assembled from the source table's molecule column while a virtual SMILES comes from RDKit's
    // linker, so the two canonicalisations need not agree: this catches matches, and is not a proof of
    // novelty — which is why a listed analog says "not in this dataset" and never "never made".
    const measuredStructures = new Set<string>();
    const test = new TopList(betterSupported);
    const analogs = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const analogsThin = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    // Everything the gate turned down, ranked on raw gain. A screen that answers "nothing qualifies"
    // and stops leaves the reader with no next structure at all, when what the gate actually withheld
    // is the confidence, not the ranking.
    const analogsAny = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const swaps = new Map<string, SwapPool>();
    const acc: RGroupAcc = {wins: new Map(), losses: new Map(), tried: new Map()};
    const series: SeriesStat[] = [];
    const lowR2: SarMatrix[] = [];
    // A strip needs a fixed partner axis and few enough slots to read as a mark; fragmented matrices
    // overlap through tiers and have neither.
    const wantStrip = host.axisRole !== null && host.matrices.length <= STRIP_SLOTS;

    let measuredCells = 0;
    let trustedNoStructure = 0;
    let trustedCells = 0;
    let withheldBelowError = 0;
    let withheldThinSupport = 0;
    let withheldFitFails = 0;
    let withheldUnchecked = 0;
    let withheldNotConverged = 0;
    const alreadyHeld = new Set<string>();
    let lowR2Virtual = 0;
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
    // Admitting an analog on its gain needs this matrix's best measured cell, which is final only when
    // its row loop ends — so the candidates wait here rather than being offered inside it.
    const riBuf: number[] = [];
    const ciBuf: number[] = [];
    // The same, for the candidates the trust gate turned down.
    const riAlt: number[] = [];
    const ciAlt: number[] = [];

    // `roleColumns` excludes the core column, so the core role is prepended: left out, the report omits
    // the one role the request names first. Empty outside fragment-columns mode, which is what keeps the
    // fit off a decomposition whose labels are local to one series.
    const roleNames = host.axisRole === null ? [] : [host.coreRole!, ...host.roleColumns];
    // Where each role's value lives, decided once: the core on the row, the axis on the column, and
    // every other role folded into the row identity under its own column name.
    const roleKind = roleNames.map((name) => name === host.coreRole ? 0 : name === host.axisRole ? 1 : 2);
    const roleValues: string[][] = roleNames.map(() => []);
    const roleActivity: number[] = [];
    const roleMol: number[] = [];
    const roleBest = roleNames.map(() => new Map<string, MatrixCellRef & {value: number}>());

    for (let mi = 0; mi < host.matrices.length; mi++) {
      // Reading at one tier: skipping the others here is what keeps every ranking below — components,
      // swaps, cores, Start here, Worth making — on the same set of series, since all of them are built
      // from this one walk.
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
            if (value !== null)
              test.offer(String(cell.molIdx), dir * value, matrix, ri, ci);
          } else if (cell.kind === 'virtual' && value !== null) {
            const support = cell.support ?? 0;
            virtualCells++;
            const passes = support >= MIN_SUPPORT && holds;
            if (passes && (bestVirtual === null || dir * value > dir * bestVirtual.value))
              bestVirtual = {ri, ci, value};
            // A structureless prediction is neither actionable nor withheld for want of evidence: the
            // core carries an attachment point no picked column fills, which is a settings fix.
            if (cell.smiles === null) {
              if (passes)
                trustedNoStructure++;
            } else if (passes) {
              trustedCells++;
              trusted++;
              trustedStructures.add(cell.smiles);
              if (measuredStructures.has(cell.smiles))
                alreadyHeld.add(cell.smiles);
              else {
                riBuf.push(ri);
                ciBuf.push(ci);
              }
            } else {
              if (!measuredStructures.has(cell.smiles)) {
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

      const fit = fitAdditiveFromTriples(obsRow, obsCol, obsVal, nRows, nCols);
      if (!fit.converged)
        nonConverged++;
      // Shrinkage toward the zero initialisation differs per design, so a gain read off a fit that
      // stopped short is not comparable with one read off a converged fit.
      const conf = matrix.confidence ?? null;
      if (!fit.converged || best === null || conf === null)
        withheldNotConverged += riBuf.length;
      else {
        // A fit whose error rests on a handful of cross-validatable cells is not a scale to divide by:
        // a plane through four points has a tiny RMSE and would top every ranking.
        const thin = conf.n < BEST_FIT_MIN_N || conf.rmse <= 0;
        for (let k = 0; k < riBuf.length; k++) {
          const cell = matrix.cells[riBuf[k]][ciBuf[k]];
          const gain = dir * (cell.value! - best.value);
          if (thin) {
            if (gain > 0)
              analogsThin.offer(cell.smiles!, gain, matrix, riBuf[k], ciBuf[k]);
          } else if (gain < conf.rmse)
            withheldBelowError++;
          else
            analogs.offer(cell.smiles!, gain / conf.rmse, matrix, riBuf[k], ciBuf[k]);
        }
      }
      // Ranked on raw gain, never on gain over error: the fits behind these are the ones that were not
      // trusted, so their error is not a scale to divide by.
      if (best !== null) {
        for (let k = 0; k < riAlt.length; k++) {
          const cell = matrix.cells[riAlt[k]][ciAlt[k]];
          const gain = dir * (cell.value! - best.value);
          if (gain > 0)
            analogsAny.offer(cell.smiles!, gain, matrix, riAlt[k], ciAlt[k]);
        }
      }
      const stat: SeriesStat = {
        matrix, root: roots[mi], tier: tiers[mi], cpd: mols.size, realCells, impossibleCells,
        totalCells: nRows * nCols, lo: realCells ? lo : 0, hi: realCells ? hi : 0, best, bestVirtual,
        typical: fit.grandMean, trusted, virtualCells,
        stripCols: recordRGroupExtremes(acc, matrix, fit, dir, roots[mi], holds, wantStrip),
        ...this.extractRanges(matrix, fit, dir),
      };
      series.push(stat);

      if (r2 === null)
        unchecked++;
      else if (r2 < TRUST_R2) {
        lowR2.push(matrix);
        lowR2Virtual += matrix.virtualCount;
      } else
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

    // The walk offers a candidate against the structures seen so far, so a molecule first measured in a
    // later matrix is only caught here.
    const dropHeld = (key: string): boolean => {
      if (!measuredStructures.has(key))
        return false;
      alreadyHeld.add(key);
      return true;
    };
    analogs.prune(dropHeld);
    analogsThin.prune(dropHeld);
    analogsAny.prune(dropHeld);

    const measured = series.filter((stat) => stat.realCells > 0);
    // Over the matrices this walk read, not over every matrix: at a tier filter a model error pooled
    // across tiers is the error of no ranking shown — it sets the band every effect is read against.
    const confs = series.flatMap((stat) => stat.matrix.confidence ?? []);
    const roleFit = fitRoleEffects({names: roleNames, values: roleValues, activity: roleActivity,
      molIdx: roleMol, minSupport: MIN_SUPPORT, higherIsBetter: host.higherIsBetter});
    const data: SummaryData = {
      compounds: compounds.size,
      untested: untested.size,
      measuredCells,
      minObserved: measured.length ? Math.min(...measured.map((stat) => stat.lo)) : null,
      maxObserved: measured.length ? Math.max(...measured.map((stat) => stat.hi)) : null,
      trustedCells,
      trustedStructures: trustedStructures.size,
      trustedNoStructure,
      modelError: confs.length ? median(confs.map((conf) => conf.rmse)) : null,
      fitHolds, unchecked, lowR2, lowR2Virtual,
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
      swapBest: largestSwap(swaps),
      ...rankRGroups(acc),
      test, analogs, analogsThin, analogsAny,
      withheldBelowError, withheldThinSupport, withheldFitFails, withheldUnchecked, withheldNotConverged,
      alreadyHeld: alreadyHeld.size,
      analogOverflow: analogs.overflow + analogsThin.overflow,
    };
    logSarTime('summary', t0);
    return data;
  }

  /** The scalars the screen needs out of one fit, so the Float64Arrays can be dropped with it. Both
   *  ranges are of centred effects, so they are comparable to each other inside this matrix — and to
   *  nothing outside it. One support floor on both axes: an effect estimated from two observations is
   *  noisier than one estimated from three, so admitting columns at a lower floor than rows would
   *  widen the column range before any chemistry entered. */
  private extractRanges(matrix: SarMatrix, fit: ReturnType<typeof fitAdditiveFromTriples>, dir: number):
    {colRange: number | null, rowRange: number | null, bestRow: {ri: number, effect: number, n: number} | null} {
    const range = (effect: Float64Array, n: Int32Array, count: number): number | null => {
      const kept = effect.subarray(0, count).filter((_e, i) => n[i] >= MIN_SUPPORT);
      return kept.length >= 2 ? Math.max(...kept) - Math.min(...kept) : null;
    };
    let bestRow: {ri: number, effect: number, n: number} | null = null;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      if (fit.rowN[ri] >= MIN_SUPPORT && (bestRow === null || dir * fit.rowEffect[ri] > dir * bestRow.effect))
        bestRow = {ri, effect: fit.rowEffect[ri], n: fit.rowN[ri]};
    }
    return {
      colRange: range(fit.colEffect, fit.colN, matrix.columns.length),
      rowRange: range(fit.rowEffect, fit.rowN, matrix.rows.length),
      bestRow,
    };
  }

  /**
   * Five independent argmaxes, then merged by lineage root.
   *
   * Independent because a series holding the most compounds, the widest range AND the best compound is
   * the answer, and pushing it out of a taken slot would put a weaker series under a label it does not
   * earn. Merged because those are one place to look, not three.
   */
  private mergeStartHere(series: SeriesStat[], confN: number[], dir: number): StartRow[] {
    if (series.length === 0)
      return [];
    // Floored well above MIN_CV_POINTS: an R² ranked at the cross-validation floor is an extreme-value
    // draw from the noisiest end, on a row captioned "best validated".
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
    // Index-aligned with REASON_GLYPHS / REASON_WORDS: the lane is a fixed set of slots, so the order
    // here is the order on screen.
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
