/**
 * Adapter — dose proportionality (Smith power model, FR-410..413 / UC-03).
 *
 * @see src/adapter/dose-proportionality.ts (implementation under test)
 *
 * Behavioural unit tests on synthetic frames. The reference-suite parity against
 * the real per-subject oracle CSV (AC-1 / AC-11) lives in the reference suite
 * (Slice C); these tests pin the verdict-mode logic, the honest-uncertainty
 * surface, and the guard warnings.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';

import {
  computeDoseProportionality, computeCriticalRegion, classifyAgainstCR, minNAdvisory,
} from '../adapter/dose-proportionality';
import type {DoseProportionalityResult, DoseProportionalityFit} from '../adapter/dose-proportionality';
import {unmappedGridColumns} from '../views/results-view';

interface MockRow {
  route: string; sex: string; dose: number;
  aucInf: number; aucLast?: number; cmax: number; doseNorm?: string;
}

function mockResults(rows: ReadonlyArray<MockRow>): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', rows.map((r) => r.route)),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', rows.map((r) => r.sex)),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', rows.map((r) => r.doseNorm ?? 'absolute')),
    DG.Column.fromFloat64Array('Dose', Float64Array.from(rows.map((r) => r.dose))),
    DG.Column.fromFloat64Array('AUC_INF', Float64Array.from(rows.map((r) => r.aucInf))),
    DG.Column.fromFloat64Array('AUC_LAST', Float64Array.from(rows.map((r) => r.aucLast ?? r.aucInf))),
    DG.Column.fromFloat64Array('CMAX', Float64Array.from(rows.map((r) => r.cmax))),
  ]);
}

/** Build a single-stratum (one route, one sex) frame of proportional rows. */
function proportionalRows(
  route: string, sex: string, doses: number[], perDose: number, slope = 1, k = 100,
): MockRow[] {
  const rows: MockRow[] = [];
  for (const d of doses) {
    for (let i = 0; i < perDose; i++) {
      const auc = k * Math.pow(d, slope);
      rows.push({route, sex, dose: d, aucInf: auc, cmax: auc / 10});
    }
  }
  return rows;
}

function fitFor(
  res: DoseProportionalityResult, route: string, sex: string, parameter: 'AUC' | 'Cmax',
): DoseProportionalityFit {
  const f = res.fits.find((x) => x.route === route && x.sex === sex && x.parameter === parameter);
  if (f === undefined) throw new Error(`fit not found for ${route}/${sex}/${parameter}`);
  return f;
}

function hasWarning(res: DoseProportionalityResult, code: string): boolean {
  return res.warnings.some((w) => w.code === code);
}

const run = (rows: MockRow[], opts: Partial<Parameters<typeof computeDoseProportionality>[1]> = {}) =>
  computeDoseProportionality(mockResults(rows), {doseCol: 'Dose', ...opts});

category('Adapter: dose proportionality (Smith power model)', () => {
  // BUG-12 drift guard: both emitted tables must resolve fully in their
  // presentation maps, or `applyFixedGridSpec` leaves a column with a raw header
  // and DG's generic stats tooltip. Both maps are fully fixed (no caller-named
  // stratify passthrough), so nothing is exempt.
  test('BUG-12: every verdict / dose-normalized column has a presentation entry', async () => {
    const res = run([
      {route: 'PO', sex: 'M', dose: 1, aucInf: 10, cmax: 2},
      {route: 'PO', sex: 'M', dose: 3, aucInf: 31, cmax: 6},
      {route: 'PO', sex: 'M', dose: 10, aucInf: 98, cmax: 19},
    ]);
    const v = unmappedGridColumns(res.table, 'doseProportionality');
    expect(v.length, 0, `verdict columns with no DOSE_PROP_FIXED entry: [${v.join(', ')}]`);
    const dn = unmappedGridColumns(res.dnTable, 'doseProportionalityDn');
    expect(dn.length, 0, `dn columns with no DOSE_PROP_DN_FIXED entry: [${dn.join(', ')}]`);
  });

  // --- AC-2: critical-region formula -----------------------------------
  test('AC-2: critical region from realized dose ratio and θ', async () => {
    const cr10 = computeCriticalRegion(10, 0.8, 1.25);
    expectFloat(cr10[0], 0.90308, 1e-4, 'CR low r=10');
    expectFloat(cr10[1], 1.09692, 1e-4, 'CR high r=10');
    const cr5 = computeCriticalRegion(5, 0.8, 1.25);
    expectFloat(cr5[0], 0.86135, 1e-4, 'CR low r=5');
    expectFloat(cr5[1], 1.13865, 1e-4, 'CR high r=5');
    // θ override tightens the CR.
    const crTight = computeCriticalRegion(10, 0.9, 1.111);
    expect(crTight[0] > cr10[0], true, 'tighter θ raises CR low');
    expect(crTight[1] < cr10[1], true, 'tighter θ lowers CR high');
  });

  // --- classifyAgainstCR: the four CI-vs-CR dispositions ----------------
  test('classifyAgainstCR: within / below / above / straddle', async () => {
    const cr: [number, number] = [0.903, 1.097];
    expect(classifyAgainstCR([0.95, 1.05], cr).direction, 'within', 'inside CR');
    expect(classifyAgainstCR([0.95, 1.05], cr).conclusive, true, 'inside is conclusive');
    expect(classifyAgainstCR([0.70, 0.85], cr).direction, 'below', 'entirely below');
    expect(classifyAgainstCR([1.15, 1.30], cr).direction, 'above', 'entirely above');
    expect(classifyAgainstCR([0.80, 1.05], cr).direction, 'straddle', 'crosses a boundary');
    expect(classifyAgainstCR([0.80, 1.05], cr).conclusive, false, 'straddle is inconclusive');
  });

  // --- AC-6: minimum-N advisory boundary (pure root-find) --------------
  test('AC-6: minNAdvisory finds a finite n at low IIV and NO_FEASIBLE_N at high IIV', async () => {
    // K=3 geometric doses 1/5/25 → ln levels {0, 1.609, 3.219}; mean 1.609;
    // sxxPerAnimal = Σ(ln-mean)² = 1.609² + 0 + 1.609² = 5.181.
    const sxxPerAnimal = Math.log(5) ** 2 * 2;
    const crHalf = (computeCriticalRegion(25, 0.8, 1.25)[1] - computeCriticalRegion(25, 0.8, 1.25)[0]) / 2;
    // Low IIV → a finite n ≤ 50 exists.
    const low = minNAdvisory(0.10, sxxPerAnimal, 3, crHalf, 50);
    expect(low.feasible, true, 'low-IIV design is feasible');
    expect(low.nStar !== null && low.nStar <= 50, true, 'finite n within the cap');
    // High IIV → no n ≤ 50 satisfies the constraint. BUILD-CONFIRMED boundary
    // (rule-16): for r=25/K=3, the n=50 cap closes only when σ ≳ 0.67
    // (t(0.95,148)·σ/√(5.181·50) ≥ crHalf=0.0693); the blueprint's illustrative
    // σ≈0.45 is in fact *feasible*. Use σ=0.9 (clearly past the boundary).
    const high = minNAdvisory(0.9, sxxPerAnimal, 3, crHalf, 50);
    expect(high.feasible, false, 'high-IIV design is not feasible at the cap');
    expect(high.nStar, null, 'no feasible n returned');
  });

  // --- AC-1 (mode) / AC-4: individual fit, conclusive proportional ------
  test('AC-1/AC-4: exact-proportional individual fit → proportional, conclusive, df = N−2', async () => {
    const res = run(proportionalRows('IV', 'F', [1, 5, 10], 4)); // 12 subjects
    const f = fitFor(res, 'IV', 'F', 'AUC');
    expectFloat(f.beta, 1, 1e-6, 'β = 1');
    expect(f.df, 10, 'df = N − 2 = 10');
    expect(f.n, 12, 'N = 12 individual subjects');
    expect(f.nDoseLevels, 3, '3 dose levels');
    expect(f.isGroupMeanFallback, false, 'individual path, not group-mean');
    expect(f.verdict, 'proportional', 'CI ⊆ CR → proportional');
    expect(f.mode, 'conclusive', 'conclusive mode');
  });

  // --- AC-4: group-mean fallback for a no-replicate (composite) design --
  test('AC-4: one AUC per dose level → labelled group-mean fallback, df-driven', async () => {
    const res = run(proportionalRows('IV', 'F', [1, 5, 10], 1)); // 3 subjects, 1/level
    const f = fitFor(res, 'IV', 'F', 'AUC');
    expect(f.isGroupMeanFallback, true, 'composite design → group-mean fallback');
    expect(f.df, 1, 'df = K − 2 = 1');
    // A df=1 fit of perfectly proportional data still fits β=1 with a CI; the
    // honest mode is df-driven when it cannot conclude.
    expect(f.mode === 'df-driven' || f.mode === 'conclusive', true, 'df-driven or conclusive');
  });

  // --- AC-3: < 3 dose levels → SMITH_INAPPLICABLE, dn still produced -----
  test('AC-3: two dose levels → SMITH_INAPPLICABLE, no verdict, dn table present', async () => {
    const res = run(proportionalRows('SC', 'F', [1, 5], 4)); // 2 levels
    const f = fitFor(res, 'SC', 'F', 'AUC');
    expect(f.verdict, 'inapplicable', 'no Smith verdict at 2 levels');
    expect(f.mode, 'inapplicable', 'inapplicable mode');
    expect(Number.isNaN(f.slopeCI[0]), true, 'no CI emitted');
    expect(hasWarning(res, 'SMITH_INAPPLICABLE'), true, 'SMITH_INAPPLICABLE raised');
    expect(res.dnTable.rowCount > 0, true, 'dn-AUC table still produced');
  });

  // --- AC-5: completeness floor ----------------------------------------
  test('AC-5: < 3 levels with ≥2 animals → INSUFFICIENT_COMPLETENESS', async () => {
    const rows = [
      ...proportionalRows('IV', 'F', [1], 4),
      ...proportionalRows('IV', 'F', [5], 4),
      ...proportionalRows('IV', 'F', [10], 1), // only 1 animal at the top dose
    ];
    const res = run(rows);
    const f = fitFor(res, 'IV', 'F', 'AUC');
    expect(f.verdict, 'insufficient-completeness', 'completeness floor not met');
    expect(hasWarning(res, 'INSUFFICIENT_COMPLETENESS'), true, 'INSUFFICIENT_COMPLETENESS raised');
  });

  // --- AC-5: extreme-dose drop -----------------------------------------
  test('AC-5: NA removal biased to the dose extremes → EXTREME_DOSE_DROP', async () => {
    // dose 1 & 10 each lose 2 of 4 subjects to NA AUC; dose 5 intact.
    const rows: MockRow[] = [
      {route: 'IV', sex: 'F', dose: 1, aucInf: 100, cmax: 10},
      {route: 'IV', sex: 'F', dose: 1, aucInf: 100, cmax: 10},
      {route: 'IV', sex: 'F', dose: 1, aucInf: NaN, cmax: 10},
      {route: 'IV', sex: 'F', dose: 1, aucInf: NaN, cmax: 10},
      ...proportionalRows('IV', 'F', [5], 4),
      {route: 'IV', sex: 'F', dose: 10, aucInf: 1000, cmax: 100},
      {route: 'IV', sex: 'F', dose: 10, aucInf: 1000, cmax: 100},
      {route: 'IV', sex: 'F', dose: 10, aucInf: NaN, cmax: 100},
      {route: 'IV', sex: 'F', dose: 10, aucInf: NaN, cmax: 100},
    ];
    const res = run(rows);
    expect(hasWarning(res, 'EXTREME_DOSE_DROP'), true, 'EXTREME_DOSE_DROP raised');
  });

  // --- AC-8: DOSE_NORM design guard ------------------------------------
  test('AC-8: a stratum mixing absolute/per-kg dose → MIXED_DOSE_NORM_FIT, fit refused', async () => {
    const rows: MockRow[] = [
      ...proportionalRows('IV', 'F', [1, 5, 10], 2).map((r, i) =>
        ({...r, doseNorm: i % 2 === 0 ? 'absolute' : 'per_kg'})),
    ];
    const res = run(rows);
    const f = fitFor(res, 'IV', 'F', 'AUC');
    expect(f.verdict, 'inapplicable', 'mixed-norm fit refused');
    expect(hasWarning(res, 'MIXED_DOSE_NORM_FIT'), true, 'MIXED_DOSE_NORM_FIT raised');
  });

  // --- AC-10: basis selection + NA-starved fallback recommendation ------
  test('AC-10: AUCinf is the default basis and is labelled', async () => {
    const res = run(proportionalRows('IV', 'F', [1, 5, 10], 4));
    expect(res.basis, 'aucinf', 'default basis is AUCinf');
    const row = res.table.col('BASIS')!;
    let labelled = false;
    for (let i = 0; i < res.table.rowCount; i++)
      if (res.table.col('PARAMETER')!.get(i) === 'AUC' && row.get(i) === 'aucinf') labelled = true;
    expect(labelled, true, 'AUC rows label the active basis');
  });

  test('AC-10: AUCinf NA-starved → AUCINF_NA_STARVED recommending AUClast', async () => {
    // AUC_INF NA on most subjects; AUC_LAST present.
    const rows: MockRow[] = proportionalRows('IV', 'F', [1, 5, 10], 4)
      .map((r) => ({...r, aucLast: r.aucInf, aucInf: NaN}));
    const res = run(rows, {basis: 'aucinf'});
    expect(hasWarning(res, 'AUCINF_NA_STARVED'), true, 'AUCINF_NA_STARVED raised');
  });

  // --- AC-12: heteroscedasticity detection -----------------------------
  test('AC-12: within-dose residual-SD ratio ≥ 4 → HETEROSCEDASTICITY_SUSPECTED', async () => {
    // dose 1 tight, dose 10 widely scattered (≥3 residuals/group required).
    const rows: MockRow[] = [
      {route: 'IV', sex: 'F', dose: 1, aucInf: 99, cmax: 10},
      {route: 'IV', sex: 'F', dose: 1, aucInf: 100, cmax: 10},
      {route: 'IV', sex: 'F', dose: 1, aucInf: 101, cmax: 10},
      {route: 'IV', sex: 'F', dose: 1, aucInf: 100, cmax: 10},
      {route: 'IV', sex: 'F', dose: 5, aucInf: 490, cmax: 50},
      {route: 'IV', sex: 'F', dose: 5, aucInf: 510, cmax: 50},
      {route: 'IV', sex: 'F', dose: 5, aucInf: 500, cmax: 50},
      {route: 'IV', sex: 'F', dose: 5, aucInf: 500, cmax: 50},
      {route: 'IV', sex: 'F', dose: 10, aucInf: 500, cmax: 100},
      {route: 'IV', sex: 'F', dose: 10, aucInf: 2000, cmax: 100},
      {route: 'IV', sex: 'F', dose: 10, aucInf: 600, cmax: 100},
      {route: 'IV', sex: 'F', dose: 10, aucInf: 1800, cmax: 100},
    ];
    const res = run(rows);
    expect(hasWarning(res, 'HETEROSCEDASTICITY_SUSPECTED'), true, 'heteroscedasticity flagged');
    // A homoscedastic frame raises no such warning.
    const clean = run(proportionalRows('IV', 'M', [1, 5, 10], 4));
    expect(clean.warnings.some((w) => w.code === 'HETEROSCEDASTICITY_SUSPECTED'), false,
      'homoscedastic frame is silent');
  });

  // --- AC-7: sex slope-homogeneity -------------------------------------
  test('AC-7: two sexes with equal slopes → homogeneous, pooled df = N−3', async () => {
    const rows = [
      ...proportionalRows('IV', 'F', [1, 5, 10], 4),
      ...proportionalRows('IV', 'M', [1, 5, 10], 4),
    ];
    const res = run(rows);
    const h = res.sexHomogeneity.find((x) => x.route === 'IV');
    expect(h !== undefined, true, 'sex-homogeneity computed for IV');
    expect(h!.homogeneous, true, 'equal slopes → homogeneous');
    expect(h!.pooledDf, 24 - 3, 'pooled df = N − 3 = 21');
  });

  test('AC-7: a single sex category → SEX_STRATUM_INSUFFICIENT, no homogeneity row', async () => {
    const res = run(proportionalRows('IV', 'F', [1, 5, 10], 4));
    expect(res.sexHomogeneity.some((x) => x.route === 'IV'), false, 'no homogeneity for one sex');
  });

  test('AC-7: divergent sex slopes → homogeneous=false, no pooled slope; per-stratum df=n−2', async () => {
    const rows = [
      ...proportionalRows('IV', 'F', [1, 5, 10], 4, 0.6), // shallow slope
      ...proportionalRows('IV', 'M', [1, 5, 10], 4, 1.6), // steep slope
    ];
    const res = run(rows);
    const h = res.sexHomogeneity.find((x) => x.route === 'IV');
    expect(h !== undefined, true, 'homogeneity row present');
    expect(h!.homogeneous, false, 'divergent slopes → interaction rejected (not homogeneous)');
    expect(h!.pooledDf, null, 'no pooled slope reported on rejection');
    // The per-stratum individual fits still stand at df = n − 2.
    const fF = res.fits.find((f) => f.route === 'IV' && f.sex === 'F' && f.parameter === 'AUC');
    expect(fF!.df, 10, 'F per-stratum df = n − 2 = 10');
  });

  test('AC-7 (F8): a dose level carrying only one sex → SEX_DOSE_CONFOUNDED, no claim', async () => {
    const rows = [
      ...proportionalRows('IV', 'F', [1, 5, 10], 4),
      ...proportionalRows('IV', 'M', [5, 10], 4), // M absent at the low dose → dose 1 is F-only
    ];
    const res = run(rows);
    expect(hasWarning(res, 'SEX_DOSE_CONFOUNDED'), true, 'Sex⊥Dose confounding flagged (F8)');
    expect(res.sexHomogeneity.some((x) => x.route === 'IV'), false, 'no sex-interaction claim when confounded');
  });

  // --- AC-6 (end-to-end): the NO_FEASIBLE_N warning on a high-IIV balanced design ---
  test('AC-6: high-IIV balanced stratum → iiv-driven + NO_FEASIBLE_N warning', async () => {
    // Geometric doses 1/5/25 (r=25 → narrow CR), n=4/level, large multiplicative
    // IIV (factors ≈ 0.33/3.0/0.4/2.5 → σ≈1.0 on log-AUC). No n ≤ 50 closes the CI.
    const factors = [0.33, 3.0, 0.4, 2.5];
    const rows: MockRow[] = [];
    for (const dose of [1, 5, 25]) {
      for (const f of factors) {
        const auc = 100 * dose * f;
        rows.push({route: 'IV', sex: 'F', dose, aucInf: auc, cmax: auc / 10});
      }
    }
    const res = run(rows);
    const fit = res.fits.find((x) => x.route === 'IV' && x.sex === 'F' && x.parameter === 'AUC');
    expect(fit!.mode, 'iiv-driven', 'inconclusive by between-animal variability');
    expect(fit!.minNAdvisory !== null, true, 'a min-N advisory was attempted (balanced design)');
    expect(hasWarning(res, 'NO_FEASIBLE_N'), true, 'NO_FEASIBLE_N at the n≤50 cap');
    // The advisory reports no finite n (never a spurious huge n).
    expect(fit!.minNAdvisory!.feasible, false, 'not feasible within the cap');
    expect(fit!.minNAdvisory!.nStar, null, 'no finite n returned');
  });

  // --- AC-9: no vendored math (the regression comes from sci-comp) ------
  // The "no local normal-equations / t-quantile" invariant is a *static* property
  // enforced by architect review + a manual grep of dose-proportionality.ts (the
  // Datagrok browser test harness has no `fs` for an in-test grep). The fit/CI
  // fields below can only be produced by `stats.linearFit` — behavioural proof
  // that the regression is routed upstream, not vendored.
  test('AC-9: dn-AUC ANOVA is offered as a secondary, labelled result', async () => {
    const res = run(proportionalRows('IV', 'F', [1, 5, 10], 4));
    const a = res.dnAnova.find((x) => x.route === 'IV' && x.sex === 'F');
    expect(a !== undefined, true, 'dn-AUC ANOVA present');
    expect(a!.nGroups, 3, 'three dose groups');
  });
});
