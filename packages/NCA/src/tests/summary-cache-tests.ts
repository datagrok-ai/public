/**
 * NCA Studies — buildSummaryCache per-(compound × route × dose-norm) category
 * (study-unified-grid AD-1 / AD-2).
 *
 * @see src/views/nca-analysis-open.ts (buildSummaryCache under test)
 *
 * The compound tier is now derived by the VAL-01-validated `aggregateNcaResults`
 * grouped by `[groupByCol, ROUTE, DOSE_NORM]` — so the grid computes NOTHING new
 * (thesis #1 / rule 18). These tests assert against REAL generated pipeline
 * output (CLAUDE.md rule 16), not a mirror:
 *  - AC-1a: 05_multicompound_batch_screen (20 compounds) → 20 per-compound rows.
 *  - AC-1b: a mixed-route analysis splits into route-homogeneous rows (IV basis
 *    never pooled with EV).
 *  - AC-1c: the cached mean/n equals a direct `aggregateNcaResults` for the same
 *    group (single source of truth — no re-derivation).
 *  - Per-param n (na.rm): a partial group's `MeanCv.n` = the aggregator's `N_P`.
 *
 * Runtime: `grok test --host local` (fixtures load via `grok.dapi.files`).
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {buildSummaryCache, buildCurvesByCompound, RESULTS_BASIS} from '../adapter';
import type {ObservationsCtx} from '../adapter';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {f64Column} from '../adapter/adapter-utils';
import {
  computeBioavailability, buildStudyGrid, CURVE_COL, buildMeanSdProfile, ROUTE_CODE_TO_STRING,
  normalizeGroupKey, NO_GROUP,
} from '../adapter';
import {DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';
import type {NcaStudy, SummaryCache, KeyValueCount} from '../adapter';
import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, multicompoundContract, multicompoundScalars,
  ratMultirouteContract, setupMultirouteRoute, ratScalars, runPipelineOrThrow,
  theophContract, theophScalars,
} from './test-utils';
import type {AnalysisRecipe} from '../adapter';
import type {PerProfileResult} from '../controller';

const ISO = '2026-07-03T00:00:00.000Z';

/** Minimal recipe carrying only the fields `buildSummaryCache` reads (the group-by
 *  set + the dose descriptor). The scalar param is a test convenience widened to the
 *  v8 `groupByCols` set. Cast through unknown (mirrors the coordinator tests). */
function recipeWith(groupByCol: string | null): AnalysisRecipe {
  return {
    columns: {doseSource: 'column', doseScalar: 0},
    units: {dose: 'mg/kg'},
    modes: {}, ncaRules: {},
    groupByCols: groupByCol === null ? [] : [groupByCol],
  } as unknown as AnalysisRecipe;
}

/** Recipe variant that also designates the per-row dose + sex columns %F needs
 *  (`computeFpctByCompound` reads `columns.doseCol`/`columns.sexCol`). */
function recipeWithFpct(
  groupByCol: string | null, doseCol: string, sexCol: string,
): AnalysisRecipe {
  return {
    columns: {doseSource: 'column', doseScalar: 0, doseCol, sexCol},
    units: {dose: 'mg/kg'},
    modes: {}, ncaRules: {},
    groupByCols: groupByCol === null ? [] : [groupByCol],
  } as unknown as AnalysisRecipe;
}


category('NCA Studies: buildSummaryCache per-(compound × route × dose-norm)', () => {
  test('AC-1a: 05 multicompound batch → 20 per-compound rows (real generated output)', async () => {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const contract = multicompoundContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, multicompoundScalars(), DEFAULT_RULES);
    const added = appendGroupingColumns(
      results, df, outcome.perProfile.map((p) => p.index), contract);
    expect(added.includes('CompoundID'), true, 'CompoundID carried into results');

    const cache = buildSummaryCache(
      results, outcome.perProfile, recipeWith('CompoundID'), 'Batch screen');
    expect(cache.compoundTier.length, 20, '20 compounds → 20 compound-tier rows');
    expect(cache.analysisTier.nCompounds, 20, 'nCompounds = 20 distinct compounds');
    // Single route (IV) + uniformly per-kg (no body weight) → each row is one
    // route-/unit-homogeneous compound of 3 subjects.
    for (const r of cache.compoundTier) {
      expect(r.route, 'IV-bolus', `compound ${r.compoundId} route = IV-bolus`);
      expect(r.doseNorm, 'per_kg', `compound ${r.compoundId} dose-norm = per_kg`);
      expect(r.n, 3, `compound ${r.compoundId} has 3 subjects`);
      // Grouped by CompoundID → each row is a real group, not the whole analysis
      // (v7, group-col-stale-on-rename).
      expect(r.isWholeAnalysis, false, `compound ${r.compoundId} isWholeAnalysis = false (grouped)`);
      // v3 (chemist-sar-grid) — MRT/Vss must be carried into the tier from real
      // pipeline output (rule 16: verified against generated fixture-05, not a mirror).
      // All-IV screen → Vss IS estimable (sci-comp gates Vss to NaN only for non-IV),
      // so both are finite with a per-param n of 3.
      expect(Number.isFinite(r.mrt.mean), true, `compound ${r.compoundId} MRT mean is finite`);
      expect(r.mrt.n, 3, `compound ${r.compoundId} MRT contributor n = 3`);
      expect(Number.isFinite(r.vss.mean), true, `compound ${r.compoundId} Vss mean is finite (all-IV)`);
      expect(r.vss.n, 3, `compound ${r.compoundId} Vss contributor n = 3 (IV → estimable)`);
    }
  });

  test('AC-1c: cached mean/n equals aggregateNcaResults for the same group (no re-derivation)', async () => {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const contract = multicompoundContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, multicompoundScalars(), DEFAULT_RULES);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);

    const cache = buildSummaryCache(
      results, outcome.perProfile, recipeWith('CompoundID'), 'Batch');
    // Independent aggregate on the same frame + keys — the single source of truth.
    const agg = aggregateNcaResults(results, ['CompoundID', 'ROUTE', 'DOSE_NORM']);
    const cidCol = agg.col('CompoundID')!;
    const meanClCol = agg.col('MEAN_CL')!;
    const nClCol = agg.col('N_CL')!;

    const row = cache.compoundTier[0];
    let m = -1;
    for (let i = 0; i < agg.rowCount; i++)
      if (cidCol.getString(i) === row.compoundId) {m = i; break;}
    expect(m >= 0, true, `aggregate row for ${row.compoundId} exists`);
    const expectedCl = meanClCol.isNone(m) ? NaN : (meanClCol.get(m) as number);
    if (Number.isNaN(expectedCl))
      expect(Number.isNaN(row.cl.mean), true, 'null aggregate CL ⇒ NaN cached mean');
    else
      expectFloat(row.cl.mean, expectedCl, 1e-9, 'cl.mean equals aggregator MEAN_CL (grid computes nothing)');
    expect(row.cl.n, nClCol.get(m), 'cl.n equals aggregator N_CL (per-param contributor count)');
  });

  test('AC-1b: mixed-route analysis → route-homogeneous rows (IV basis never pools EV)', async () => {
    const df = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(df);
    const contract = ratMultirouteContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);

    // No group-by → the tier splits by route × dose-norm alone (body weight is
    // present, so DOSE_NORM is uniformly absolute → one row per route).
    const cache = buildSummaryCache(results, outcome.perProfile, recipeWith(null), 'Multiroute');
    const routes = new Set(cache.compoundTier.map((r) => r.route));
    expect(routes.has('IV-bolus'), true, 'IV-bolus case present');
    expect(routes.has('PO'), true, 'PO case present');
    expect(routes.has('SC'), true, 'SC case present');
    // Every row is route-homogeneous: no row spans two routes (route is a key).
    expect(cache.compoundTier.length >= 3, true, 'at least one row per route (no pooling)');
    // No group-by → every row IS the whole analysis (v7); the study grid shows "—"
    // and this stays true across a later rename (label-independent sentinel).
    for (const r of cache.compoundTier)
      expect(r.isWholeAnalysis, true, `${r.route} row isWholeAnalysis = true (ungrouped)`);
  });

  test('per-param n: a partial group carries MeanCv.n = aggregator N_P (na.rm)', async () => {
    // Results-shaped frame: 1 compound, 3 subjects, one unfit on AUC_INF (NaN).
    // buildSummaryCache must report AUC∞ mean 150 over the 2 valid values with
    // n=2, while CMAX (all fit) keeps n=3 — the honest Mean|%CV|n triple.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['X', 'X', 'X']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus', 'IV-bolus']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['absolute', 'absolute', 'absolute']),
      // Built via `f64Column`, exactly as the production writer builds a results
      // frame (nca-dg-native-missing-values F1) — an unfit AUC_INF reaches the
      // aggregator as the PLATFORM null, not as a raw NaN that DG would count as a
      // present value and let poison the group mean.
      f64Column('AUC_INF', [100, NaN, 200]),
      f64Column('CMAX', [10, 20, 30]),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const perProfile = [stub(), stub(), stub()];

    const cache = buildSummaryCache(results, perProfile, recipeWith('CompoundID'), 'X');
    expect(cache.compoundTier.length, 1, 'one compound group');
    const r = cache.compoundTier[0];
    expectFloat(r.aucInf.mean, 150, 1e-6, 'AUC∞ mean = mean of [100, 200] (NaN dropped)');
    expect(r.aucInf.n, 2, 'AUC∞ n = 2 (only 2 of 3 subjects fit)');
    expect(r.cmax.n, 3, 'CMAX n = 3 (all subjects fit)');
    expect(r.n, 3, 'group subject count = 3');
  });

  test('AC-1d: one compound, mixed body-weight availability → separate rows per DOSE_NORM', async () => {
    // Same CompoundID + route, but subjects mix DOSE_NORM (absolute vs per_kg —
    // the mixed-body-weight case). buildSummaryCache MUST split into 2 rows so no
    // group averages L/h with L/h/kg CL (the unit-mixing R1-F1 guards against).
    // Synthetic results frame — the same technique the per-param-n test uses; a
    // full pipeline-driven mixed-BW CSV fixture is the separate DG-4 test-gap.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['X', 'X', 'X', 'X']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE',
        ['IV-bolus', 'IV-bolus', 'IV-bolus', 'IV-bolus']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM',
        ['absolute', 'absolute', 'per_kg', 'per_kg']),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 1.2, 5.0, 5.4])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 12, 10, 12])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 120, 100, 120])),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const perProfile = [stub(), stub(), stub(), stub()];

    const cache = buildSummaryCache(results, perProfile, recipeWith('CompoundID'), 'X');
    expect(cache.compoundTier.length, 2, 'mixed dose-norm → 2 rows (split per DOSE_NORM)');
    const abs = cache.compoundTier.find((r) => r.doseNorm === 'absolute');
    const perKg = cache.compoundTier.find((r) => r.doseNorm === 'per_kg');
    expect(abs !== undefined && perKg !== undefined, true, 'both dose-norm rows present');
    // 1e-6 (not 1e-9): DG column storage is Float32-influenced, so a mean of a
    // non-exactly-representable value (1.1) carries ~1e-8 slack — the split
    // (2 unpooled rows), not bit-exactness, is what AC-1d proves.
    expectFloat(abs!.cl.mean, 1.1, 1e-6, 'absolute CL mean = mean of [1.0,1.2] (NOT pooled with per-kg)');
    expectFloat(perKg!.cl.mean, 5.2, 1e-6, 'per-kg CL mean = mean of [5.0,5.4] (NOT pooled with absolute)');
  });

  test('S2 %F: fixture 05 (all-IV) → every compound-tier row fPct === null (no EV arm)', async () => {
    // Rule 16 against real generated output: fixture 05 is 20 IV compounds, no EV
    // arm anywhere. Every tier row is IV, so %F is null on all (F ≡ 1 by definition
    // on IV) — and the study grid therefore drops the whole %F column (Decision 2).
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const contract = multicompoundContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, multicompoundScalars(), DEFAULT_RULES);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);
    const cache = buildSummaryCache(
      results, outcome.perProfile, recipeWith('CompoundID'), 'Batch');
    for (const r of cache.compoundTier)
      expect(r.fPct, null, `compound ${r.compoundId} (IV) → fPct null (no EV arm)`);
    // End-to-end: the study grid drops the %F column entirely on this all-IV screen.
    const grid = buildStudyGrid({
      schemaVersion: 4, id: 's', label: 'S', createdAt: ISO, updatedAt: ISO,
      permissions: {view: [], edit: []},
      analyses: [{
        id: 'a', label: 'Batch', source: {kind: 'tableId', ref: 't'},
        columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols: ['CompoundID'],
        sourceSignature: {kind: 'content-hash', rowCount: 600, mtime: null, hash: 'h', computedAt: ISO},
        summaryCache: cache,
      }],
    } as unknown as NcaStudy);
    expect(grid.col('%F Mean'), null, 'all-IV screen → no %F column (Decision 2)');
  });

  test('S2 %F: paired IV+EV → EV compound row fPct equals computeBioavailability (no re-derivation)', async () => {
    // Synthetic paired frame (same technique as AC-1d): compound A with IV + PO
    // subjects across M/F. The EV (PO) tier row's fPct must equal the SAME
    // computeBioavailability → aggregate F% (rule 5: the ratio is not re-derived).
    // IV reference is sex-matched WITHIN the compound (mandatory: F% pairs like with like).
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A', 'A', 'A']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['per_kg', 'per_kg', 'per_kg', 'per_kg']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'F', 'M', 'F']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 10, 10])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 200, 50, 40])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 20, 5, 4])),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const perProfile = [stub(), stub(), stub(), stub()];

    const cache = buildSummaryCache(
      results, perProfile, recipeWithFpct('CompoundID', 'Dose', 'Sex'), 'Paired');

    // Independent computeBioavailability + aggregate — the single source of truth.
    const {table} = computeBioavailability(results, {
      doseCol: 'Dose', routeCol: 'ROUTE', stratifyCols: ['CompoundID', 'Sex'],
    });
    const fAgg = aggregateNcaResults(table, ['CompoundID', 'ROUTE'], {params: ['F_PCT']});
    // A/PO: F%(M)=100·(50/10)/(100/10)=50; F%(F)=100·(40/10)/(200/10)=20; mean=35.
    let expMean = NaN;
    const cid = fAgg.col('CompoundID')!; const rte = fAgg.col('ROUTE')!;
    const meanF = fAgg.col('MEAN_F_PCT')!;
    for (let i = 0; i < fAgg.rowCount; i++)
      if (cid.getString(i) === 'A' && rte.getString(i) === 'PO') expMean = meanF.get(i) as number;
    expectFloat(expMean, 35, 1e-9, 'sanity: direct aggregate A/PO %F mean = 35');

    const poRow = cache.compoundTier.find((r) => r.route === 'PO')!;
    const ivRow = cache.compoundTier.find((r) => r.route === 'IV-bolus')!;
    expect(poRow.fPct !== null, true, 'EV (PO) row carries a %F');
    expectFloat(poRow.fPct!.mean, expMean, 1e-9,
      'EV row fPct.mean equals computeBioavailability aggregate (no re-derivation)');
    expect(poRow.fPct!.n, 2, 'EV row %F n = 2 EV strata (M, F)');
    expect(ivRow.fPct, null, 'IV row fPct is null (F ≡ 1 by definition)');
  });

  test('SF-1 AC-F3.5 composite CompoundID-Sex confines fPct to the cell no cross-sex pooling', async () => {
    // Same paired frame as the %F test, but grouped by the COMPOSITE (CompoundID, Sex).
    // %F(A/M/PO) = 100·(50/10)/(100/10) = 50; %F(A/F/PO) = 100·(40/10)/(200/10) = 20.
    // If SF-1 were broken (naive groupByCols[0]='CompoundID' passed to computeFpctByCompound),
    // BOTH composite rows would carry the POOLED mean 35 while their Group labels show M and
    // F separately — the exact BUG-08-class silent cross-dimension pooling this cycle prevents.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A', 'A', 'A']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['per_kg', 'per_kg', 'per_kg', 'per_kg']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'F', 'M', 'F']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 10, 10])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 200, 50, 40])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 20, 5, 4])),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const recipe = {
      columns: {doseSource: 'column', doseScalar: 0, doseCol: 'Dose', sexCol: 'Sex'},
      units: {dose: 'mg/kg'}, modes: {}, ncaRules: {},
      groupByCols: ['CompoundID', 'Sex'],
    } as unknown as AnalysisRecipe;
    const cache = buildSummaryCache(
      results, [stub(), stub(), stub(), stub()], recipe, 'Composite');
    // The PO tier is now split by Sex — each composite cell's %F confined to its own subjects.
    const poRows = cache.compoundTier.filter((r) => r.route === 'PO');
    expect(poRows.length, 2, 'PO tier split into (A,M) and (A,F) composite cells');
    const means = poRows.map((r) => r.fPct!.mean).sort((a, b) => a - b);
    expectFloat(means[0], 20, 1e-9, 'A/F/PO %F = 20 (confined to F, not pooled)');
    expectFloat(means[1], 50, 1e-9, 'A/M/PO %F = 50 (confined to M, not pooled)');
    expect(means.some((m) => Math.abs(m - 35) < 1e-9), false,
      'the cross-Sex POOLED mean (35) does NOT appear — SF-1 confinement holds (no BUG-08 pooling)');
    // Each composite row carries per-key values for BOTH keys (Variant-2 projection).
    expect(poRows[0].groupValues.length, 2, 'composite row carries per-key values for CompoundID + Sex');
  });

  test('SF-1 AC-F3.5 composite CompoundID-Sex confines the curve thumbnail to the cell', async () => {
    // The curve half of SF-1 (the %F half is the test above). Group compound A by the
    // COMPOSITE (CompoundID, Sex). M and F subjects have DELIBERATELY different curves.
    // buildSummaryCache must synthesize the composite key onto the OBSERVATIONS frame and
    // key buildCurvesByCompound by it, so each tier row's thumbnail is confined to its own
    // Sex. If SF-1 were broken (a single 'CompoundID' column keyed the curve), M and F
    // would share ONE pooled series (the BUG-08 cross-dimension pooling this guards) —
    // and the composite lookup would miss, yielding null curves. Either way this fails.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['per_kg', 'per_kg']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'F']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 50])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([50, 5])),
    ]);
    // Per-(subject × time) observations: M curve is ~10× the F curve at every timepoint.
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A', 'A', 'A']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'M', 'F', 'F']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([50, 25, 5, 2.5])),
    ]);
    const obs: ObservationsCtx = {
      df: obsDf, timeName: 'Time', concName: 'Conc', doseName: null, routeCodeName: null};
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const recipe = {
      columns: {doseSource: 'column', doseScalar: 0, doseCol: 'Dose', sexCol: 'Sex'},
      units: {dose: 'mg/kg'}, modes: {}, ncaRules: {},
      groupByCols: ['CompoundID', 'Sex'],
    } as unknown as AnalysisRecipe;
    const cache = buildSummaryCache(results, [stub(), stub()], recipe, 'Composite', obs);

    const rowM = cache.compoundTier.find((r) => r.groupValues.some((g) => g.value === 'M'))!;
    const rowF = cache.compoundTier.find((r) => r.groupValues.some((g) => g.value === 'F'))!;
    expect(rowM.curve !== null && rowM.curve.series.length > 0, true, 'M cell carries a curve (not null)');
    expect(rowF.curve !== null && rowF.curve.series.length > 0, true, 'F cell carries a curve (not null)');
    // Confinement: the M cell's t=1 geomean is 50 (M subject only), the F cell's is 5 —
    // NOT the pooled geomean of {50,5} (≈15.8). Distinct curves ⇒ no cross-Sex pooling.
    const t1 = (r: typeof rowM): number => {
      const s = r.curve!.series[0];
      return s.geoMean[s.time.indexOf(1)];
    };
    expectFloat(t1(rowM), 50, 1e-9, 'M cell curve confined to M subjects (geomean 50 at t=1)');
    expectFloat(t1(rowF), 5, 1e-9, 'F cell curve confined to F subjects (geomean 5 at t=1)');
    expect(Math.abs(t1(rowM) - t1(rowF)) > 1e-6, true, 'M and F curves are distinct (not one pooled series)');
    // E3 fresh-curve relabel through buildSummaryCache (R2 finding 2): the pooled series'
    // `group` is the row's group values joined by ' · ' (composite key rendered readably),
    // NOT the internal ␁-join — the SAME label the v8→v9 migration writes on a preserved curve.
    expect(rowM.curve!.series[0].group, 'A · M', 'M cell series group relabeled to the row group (A · M)');
    expect(rowF.curve!.series[0].group, 'A · F', 'F cell series group relabeled to the row group (A · F)');
  });

  test('a results-only key (CL_BASIS) keys the curve thumbnail per profile', async () => {
    // CL_BASIS lives on the results frame only — the observations carry no such column. Each
    // tier row's curve must still be its own cell's (projected per profile through
    // perProfile[i].index.observationRows), not null (the pre-fix lookup missed every row).
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['per_kg', 'per_kg', 'per_kg', 'per_kg']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL', 'CL', 'CL/F', 'CL/F']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'F', 'M', 'F']),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 90, 50, 40])),
    ]);
    // Two rows per profile, in profile order; the t=1 value identifies the profile.
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'M', 'F', 'F', 'M', 'M', 'F', 'F']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 1, 2, 1, 2, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([80, 40, 60, 30, 20, 10, 8, 4])),
    ]);
    const obs: ObservationsCtx = {
      df: obsDf, timeName: 'Time', concName: 'Conc', doseName: null, routeCodeName: null};
    const stub = (i: number): PerProfileResult => ({
      index: {subjectKey: `S${i}`, observationRows: Int32Array.from([2 * i, 2 * i + 1])},
      route: i < 2 ? nca.ROUTE_IV_BOLUS : nca.ROUTE_PO,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const t1 = (c: SummaryCache, basis: string, sex: string): number | null => {
      const r = c.compoundTier.find((x) => x.groupValues.some((g) => g.value === basis) &&
        (x.groupValues.length === 1 || x.groupValues.some((g) => g.value === sex)));
      if (r === undefined || r.curve === null) return null;
      const sr = r.curve.series[0];
      return sr.geoMean[sr.time.indexOf(1)];
    };
    const recipe = (groupByCols: string[]): AnalysisRecipe => ({
      columns: {doseSource: 'value', doseScalar: 1}, units: {dose: 'mg/kg'}, modes: {}, ncaRules: {},
      groupByCols,
    } as unknown as AnalysisRecipe);
    const profiles = [0, 1, 2, 3].map(stub);

    const composite = buildSummaryCache(results, profiles, recipe(['CL_BASIS', 'Sex']), 'X', obs);
    expect(composite.compoundTier.length, 4, 'one tier row per CL basis × Sex');
    expectFloat(t1(composite, 'CL', 'M') ?? NaN, 80, 1e-9, 'CL · M: its own profile');
    expectFloat(t1(composite, 'CL', 'F') ?? NaN, 60, 1e-9, 'CL · F: its own profile');
    expectFloat(t1(composite, 'CL/F', 'M') ?? NaN, 20, 1e-9, 'CL/F · M: its own profile');
    expectFloat(t1(composite, 'CL/F', 'F') ?? NaN, 8, 1e-9, 'CL/F · F: its own profile');

    const single = buildSummaryCache(results, profiles, recipe(['CL_BASIS']), 'X', obs);
    expectFloat(t1(single, 'CL', '') ?? NaN, Math.sqrt(80 * 60), 1e-9, 'CL: the two IV profiles');
    expectFloat(t1(single, 'CL/F', '') ?? NaN, Math.sqrt(20 * 8), 1e-9, 'CL/F: the two PO profiles');
  });

  test('S2 %F: EV compound with NO IV reference → NE (NaN) fPct + surfaced no-IV-reference warning', async () => {
    // Two compounds: A is paired (IV+PO → finite %F), B is EV-only (PO, no IV arm
    // at all → its %F is not estimable). B's PO tier row must be the NE sentinel
    // (NaN mean) AND carry NO_IV_REFERENCE_FOR_STRATUM in its warnings so a chemist
    // sees the blank %F explained on the Flags tooltip — never a silent blank (rule
    // 14). A's presence keeps the %F column visible (Decision 2). This exercises the
    // no-IV-reference path end-to-end (the coverage gap the review flagged).
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A', 'B', 'B']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'PO', 'PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['per_kg', 'per_kg', 'per_kg', 'per_kg']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'M', 'M', 'M']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 10, 10])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 60, 80, 90])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 6, 8, 9])),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const cache = buildSummaryCache(
      results, [stub(), stub(), stub(), stub()],
      recipeWithFpct('CompoundID', 'Dose', 'Sex'), 'Mixed');

    const aPo = cache.compoundTier.find((r) => r.compoundId === 'A' && r.route === 'PO')!;
    const bPo = cache.compoundTier.find((r) => r.compoundId === 'B' && r.route === 'PO')!;
    // A/PO paired → finite %F, no no-IV warning.
    expect(aPo.fPct !== null && Number.isFinite(aPo.fPct.mean), true, 'A/PO has a finite %F');
    expect(aPo.warnings.some((w) => w.code === 'NO_IV_REFERENCE_FOR_STRATUM'), false,
      'A/PO carries no no-IV warning');
    // B/PO no IV arm → NE sentinel (non-null MeanCv with NaN mean) + surfaced warning.
    expect(bPo.fPct !== null, true, 'B/PO fPct is the NE sentinel (non-null), not IV-style null');
    expect(Number.isNaN(bPo.fPct!.mean), true, 'B/PO %F mean is NaN (not estimable — no IV reference)');
    const bPoNoIv = bPo.warnings.find((w) => w.code === 'NO_IV_REFERENCE_FOR_STRATUM');
    expect(bPoNoIv !== undefined, true,
      'B/PO surfaces NO_IV_REFERENCE_FOR_STRATUM (blank %F is explained, not silently dropped)');
    // v6: the warning also carries a human-readable message for the Flags tooltip.
    expect((bPoNoIv?.message ?? '').length > 0, true, 'the no-IV warning carries a human message');
    // The column is still present overall (A carries a finite %F).
    expect(cache.compoundTier.some((r) => r.fPct !== null && Number.isFinite(r.fPct.mean)), true,
      '%F column stays visible — at least one compound is paired');
  });

  test('group-by an implicit sub-key (ROUTE) → compoundId is the analysis label, not the route value', async () => {
    // Grouping by ROUTE (already an implicit sub-key) must NOT put the route
    // string into compoundId — that would duplicate the Route column (design B2:
    // route is a sub-key, not a compound). The rows still split by route.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM',
        ['absolute', 'absolute', 'absolute', 'absolute']),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 12, 8, 9])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 120, 80, 90])),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1, 1.2, 4, 4.4])),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const cache = buildSummaryCache(
      results, [stub(), stub(), stub(), stub()], recipeWith('ROUTE'), 'Multiroute');
    expect(cache.compoundTier.length, 2, 'two route rows (IV / PO)');
    for (const r of cache.compoundTier)
      expect(r.compoundId, 'Multiroute', 'compoundId = analysis label, not the route value');
    const routes = new Set(cache.compoundTier.map((r) => r.route));
    expect(routes.has('IV-bolus') && routes.has('PO'), true, 'rows still split by route');
  });

  test('E4 (AC-E4f/c): per-tier-row Sex inventory — route-split rows disclose DIFFERENT compositions', async () => {
    // Grouped ONLY by ROUTE → whole-analysis rows, one per route. Sex is present but NOT
    // grouped. IV arm = {F×1, M×1}; PO arm = {F×2, M×1}. Each route row's keyInventory must
    // carry ITS OWN Sex composition (per-subpopulation scope, R2 finding 3) — never one
    // analysis-wide count. NCA result params (AUC_INF/CMAX) must NOT enter the inventory.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO', 'PO']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM',
        ['per_kg', 'per_kg', 'per_kg', 'per_kg', 'per_kg']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['F', 'M', 'F', 'F', 'M']),
      // A low-cardinality Dose column (a grouping factor, NOT grouped here) → must ENTER the
      // inventory (AC-E4c "Dose IS included"); IV arm dosed 10, PO arm dosed 30.
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 30, 30, 30])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 110, 80, 85, 90])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 11, 8, 8.5, 9])),
    ]);
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const cache = buildSummaryCache(
      results, [stub(), stub(), stub(), stub(), stub()], recipeWith('ROUTE'), 'RouteSplit');
    const iv = cache.compoundTier.find((r) => r.route === 'IV-bolus')!;
    const po = cache.compoundTier.find((r) => r.route === 'PO')!;
    const invMap = (r: typeof iv, key: string): Record<string, number> => {
      const vs: readonly KeyValueCount[] = r.keyInventory.find((e) => e.key === key)?.values ?? [];
      return Object.fromEntries(vs.map((v) => [v.value, v.count]));
    };
    expect(invMap(iv, 'Sex').F, 1, 'IV row: F×1 (its OWN subpopulation)');
    expect(invMap(iv, 'Sex').M, 1, 'IV row: M×1');
    expect(invMap(po, 'Sex').F, 2, 'PO row: F×2 — DIFFERENT from IV (per-subpopulation scope, AC-E4f)');
    expect(invMap(po, 'Sex').M, 1, 'PO row: M×1');
    // AC-E4c: a low-cardinality Dose column IS in the inventory, per-subpopulation. Assert
    // presence + summed subject count (format-robust — a numeric Dose stringifies to "10"/"10.0").
    const ivDose = iv.keyInventory.find((e) => e.key === 'Dose');
    const poDose = po.keyInventory.find((e) => e.key === 'Dose');
    expect(ivDose !== undefined && poDose !== undefined, true, 'Dose IS a candidate in the inventory (AC-E4c)');
    expect(ivDose!.values.length, 1, 'IV row: one distinct Dose level');
    expect(ivDose!.values.reduce((a, v) => a + v.count, 0), 2, 'IV row: Dose covers its 2 subjects');
    expect(poDose!.values.reduce((a, v) => a + v.count, 0), 3, 'PO row: Dose covers its 3 subjects');
    // AC-E4c: NCA result parameters / continuous columns are NOT in the inventory.
    expect(iv.keyInventory.some((e) => e.key === 'AUC_INF' || e.key === 'CMAX'), false,
      'NCA result params excluded from the inventory (AC-E4c: role/param gate)');
  });
});

// ── S4a: inline per-compound curve thumbnail + Tmax carry-through (v5) ────────
category('NCA Studies: buildSummaryCache inline curve + Tmax (chemist-sar-grid S4a)', () => {
  test('fixture 05 → each compound row carries a semi-log geomean curve (1 series/dose) + Tmax', async () => {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const contract = multicompoundContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, multicompoundScalars(), DEFAULT_RULES);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);

    // Observations substrate = the raw fixture frame (SubjectID/Time_hr/Conc_ng_mL/
    // Dose_mg_kg/CompoundID) — the same shape the mapping pane hands over. Verified
    // against real generated output (rule 16), not a mirror.
    const obs: ObservationsCtx = {
      df, timeName: 'Time_hr', concName: 'Conc_ng_mL', doseName: 'Dose_mg_kg',
      routeCodeName: null, // fixture 05 is single-route (all IV)
    };
    const cache = buildSummaryCache(
      results, outcome.perProfile, recipeWith('CompoundID'), 'Batch screen', obs);
    expect(cache.compoundTier.length, 20, '20 compounds → 20 rows (each a thumbnail row)');
    for (const r of cache.compoundTier) {
      expect(r.curve !== null, true, `compound ${r.compoundId} has a cached curve`);
      // Fixture 05 is single-dose per compound → exactly ONE series (BUG-08 no pooling).
      expect(r.curve!.series.length, 1, `compound ${r.compoundId}: one series (single dose)`);
      const s = r.curve!.series[0];
      expect(s.time.length === s.geoMean.length && s.time.length >= 2, true,
        `compound ${r.compoundId}: ≥ 2 paired polyline points`);
      // Ascending time (polyline order) + all measurable (> 0, semi-log safe).
      for (let i = 1; i < s.time.length; i++)
        expect(s.time[i] > s.time[i - 1], true, `compound ${r.compoundId}: time strictly ascending`);
      for (const g of s.geoMean)
        expect(g > 0, true, `compound ${r.compoundId}: geoMean point positive (semi-log)`);
      // v5 Tmax carry-through (aggregate already present) for the drill-in card.
      expect(Number.isFinite(r.tmax.mean), true, `compound ${r.compoundId}: Tmax carried into tier`);
      expect(r.tmax.n, 3, `compound ${r.compoundId}: Tmax contributor n = 3`);
    }

    // No observations plumbed → curve null (honest "not-yet-computed", never faked).
    const cacheNoObs = buildSummaryCache(
      results, outcome.perProfile, recipeWith('CompoundID'), 'Batch screen');
    for (const r of cacheNoObs.compoundTier)
      expect(r.curve, null, `compound ${r.compoundId}: no observations → curve null`);
  });

  test('BUG-08 (E3): a multi-dose compound pools into ONE series + discloses CROSS_DOSE_POOL', async () => {
    // E3 (study-grid-grouping-columns) drops the vestigial inner dose-split: a compound
    // sampled across ≥2 dose levels pools into ONE geomean series (matching the analysis
    // view), and the ≥2-dose-cluster pool is DISCLOSED via `crossDoseOuters` — BUG-08 honesty
    // by DISCLOSURE, not by splitting (the rule-14-ratified pool+disclose decision). This
    // REPLACES the pre-E3 "one series per dose" contract. Two jittered dose arms (~10 and
    // ~100) so the log-space clusterer resolves ≥2 clusters. Real function output (rule 16).
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound',
        ['C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1']),
      // ~10 arm (9.5/10/10.5) + ~100 arm (95/100/105), jittered per-kg → 2 clusters.
      DG.Column.fromFloat64Array('Dose',
        Float64Array.from([9.5, 10, 10.5, 95, 100, 105, 10, 100])),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 1, 1, 1, 1, 1, 2, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([100, 110, 120, 900, 1000, 1100, 50, 500])),
    ]);
    const curves = buildCurvesByCompound(
      {df: obsDf, timeName: 'Time', concName: 'Conc', doseName: 'Dose',
        routeCodeName: null}, 'Compound', 'X');
    expect(curves.routeKeyed, false, 'no route column → compound-level keys');
    const c1 = curves.byKey.get('C1');
    expect(c1 !== undefined, true, 'C1 curve built');
    expect(c1!.series.length, 1, 'ONE pooled series (E3 drops the dose-split — never ~N dose lines)');
    expect(c1!.series[0].time.length, 2, 'the pooled series spans both timepoints');
    // The ≥2-dose-cluster pool is DISCLOSED (BUG-08 honesty via disclosure, not splitting).
    expect(curves.crossDoseOuters.has('C1'), true,
      'multi-dose-cluster pool disclosed via crossDoseOuters (→ CROSS_DOSE_POOL warning)');
  });

  test('no dose column → a single pooled series (group = compound identity), no cross-dose disclosure', async () => {
    // When no dose is mapped there is nothing to split on — the honest result is one pooled
    // series per compound. E3 labels the pooled series' `group` by the compound identity
    // (buildSummaryCache later relabels it to the tier's group value); no dose → no disclosure.
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound', ['C1', 'C1', 'C1', 'C1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([100, 50, 120, 60])),
    ]);
    const curves = buildCurvesByCompound(
      {df: obsDf, timeName: 'Time', concName: 'Conc', doseName: null,
        routeCodeName: null}, 'Compound', 'X');
    const c1 = curves.byKey.get('C1');
    expect(c1 !== undefined && c1!.series.length === 1, true, 'one pooled series when no dose mapped');
    expect(c1!.series[0].group, 'C1', 'pooled series group = the compound identity (E3)');
    expect(c1!.series[0].time.length, 2, 'both timepoints present in the pooled series');
    expect(curves.crossDoseOuters.has('C1'), false, 'no dose column → no cross-dose disclosure');
  });

  test('AC-E3e: theoph per-kg doses → ONE pooled series, NO false CROSS_DOSE_POOL (geometric spacing)', async () => {
    // The Rule-16 empirical pin: theoph's 10 distinct per-kg dose values sit in one narrow
    // band (1 nominal level, a covariate) → the log-space clusterer resolves 1 cluster, so a
    // whole-analysis theoph thumbnail is ONE pooled curve with NO spurious cross-dose warning.
    const df = await loadFixture('theoph.csv');
    const obs: ObservationsCtx = {
      df, timeName: 'Time', concName: 'conc', doseName: 'Dose', routeCodeName: null};
    const curves = buildCurvesByCompound(obs, null, 'theoph'); // null compound col → whole-analysis
    const c = curves.byKey.get('theoph');
    expect(c !== undefined && c!.series.length === 1, true, 'theoph → ONE pooled series (whole-analysis)');
    expect(curves.crossDoseOuters.has('theoph'), false,
      'theoph per-kg dose band is ONE cluster → NO false CROSS_DOSE_POOL (AC-E3e)');
  });

  test('CROSS_DOSE_POOL degenerate-case fix: a CLEAN 2-level pool fires; a close 2-value pair does not', async () => {
    // The decision-audit fix (2026-07-13): `countDoseClusters` is degenerate at EXACTLY 2
    // distinct dose values (always returns 1), so a clean 2-arm design (10/100, no per-kg
    // jitter) previously pooled SILENTLY. The absolute max/min ≥ 2× ratio fallback (applied
    // only at exactly 2 distinct doses) now discloses it. A close 2-value pair (10/12, 1.2×,
    // one arm with two body weights) must NOT fire.
    const twoArm = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound', ['C1', 'C1', 'C1', 'C1']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 100, 100])), // 2 clean levels, 10×
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([100, 50, 900, 450])),
    ]);
    const cTwoArm = buildCurvesByCompound(
      {df: twoArm, timeName: 'Time', concName: 'Conc', doseName: 'Dose', routeCodeName: null}, 'Compound', 'X');
    expect(cTwoArm.byKey.get('C1')!.series.length, 1, 'clean 2-level pool → ONE pooled series');
    expect(cTwoArm.crossDoseOuters.has('C1'), true,
      'clean 2-level pool (10/100, 10×) → CROSS_DOSE_POOL fires (degenerate-case ratio fallback)');
    const closePair = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound', ['C1', 'C1']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 12])), // 2 values, 1.2× — one arm
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 1])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([100, 110])),
    ]);
    const cClose = buildCurvesByCompound(
      {df: closePair, timeName: 'Time', concName: 'Conc', doseName: 'Dose', routeCodeName: null}, 'Compound', 'X');
    expect(cClose.crossDoseOuters.has('C1'), false,
      'a close 2-value pair (10/12, 1.2× < 2×) does NOT fire (within-arm body-weight jitter)');
  });

  test('F2 (v11): buildCurvesByCompound persists the ×/÷ gSD band == buildMeanSdProfile output', async () => {
    // context-panel-ct-profile AC-F2.1/F2.4: a FRESH curve carries geoLower/geoUpper (PRESENT,
    // not undefined) equal to buildMeanSdProfile's band for the same pooled points; a genuine
    // n=1 timepoint carries FLOAT_NULL (honest no-whisker). Real function output cross-checked
    // against the validated aggregator (rule 16) — t=1 has 3 obs (band estimable), t=2 has 1
    // (band FLOAT_NULL).
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound', ['C1', 'C1', 'C1', 'C1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 1, 1, 2])), // t=1 ×3, t=2 ×1
      DG.Column.fromFloat64Array('Conc', Float64Array.from([100, 110, 120, 50])),
    ]);
    const curves = buildCurvesByCompound(
      {df: obsDf, timeName: 'Time', concName: 'Conc', doseName: null, routeCodeName: null}, 'Compound', 'X');
    const s = curves.byKey.get('C1')!.series[0];
    expect(s.geoLower !== undefined && s.geoUpper !== undefined, true,
      'fresh curve carries the band arrays (present — distinguishes it from a legacy band-absent curve)');
    expect(s.geoLower!.length, s.time.length, 'band parallel to time/geoMean');
    // Cross-check element-for-element against buildMeanSdProfile over the SAME pooled frame.
    const grp = obsDf.clone();
    grp.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'G', new Array<string>(4).fill('C1')));
    const prof = buildMeanSdProfile(grp, {timeCol: 'Time', concCol: 'Conc', groupCol: 'G'});
    const pT = prof.col('Time')!; const pLo = prof.col('geoLower')!; const pUp = prof.col('geoUpper')!;
    const bandByTime = new Map<number, {lo: number; up: number; loNull: boolean}>();
    for (let i = 0; i < prof.rowCount; i++) {
      bandByTime.set(pT.get(i) as number,
        {lo: pLo.get(i) as number, up: pUp.get(i) as number, loNull: pLo.isNone(i)});
    }
    for (let i = 0; i < s.time.length; i++) {
      const ref = bandByTime.get(s.time[i])!;
      if (ref.loNull)
        expect(s.geoLower![i], DG.FLOAT_NULL, `t=${s.time[i]} n=1 → band FLOAT_NULL (honest no-whisker)`);
      else {
        expectFloat(s.geoLower![i], ref.lo, 1e-9, `geoLower[${i}] == buildMeanSdProfile`);
        expectFloat(s.geoUpper![i], ref.up, 1e-9, `geoUpper[${i}] == buildMeanSdProfile`);
      }
    }
  });

  test('F2: a fresh present band survives the JSON round-trip of the cache column (AC-F2.5)', async () => {
    // AC-F2.5: a fresh curve's PRESENT band survives save→load. The summary cache is stored as an
    // EMS json column — a JSON round-trip — so serialize a study whose curve carries
    // geoLower/geoUpper and assert the band is intact.
    const v11 = {
      id: 'v11', label: 'Fresh', kind: 'study', createdAt: ISO, updatedAt: ISO,
      access: {canEdit: true, canDelete: true, canShare: true},
      analyses: [{
        id: 'a', label: 'A', source: {kind: 'tableId', ref: 't'},
        columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols: [], passthroughCols: [],
        sourceSignature: {kind: 'content-hash', rowCount: 6, mtime: null, hash: 'h', computedAt: ISO},
        summaryCache: {
          analysisTier: {nCompounds: 1, nSubjects: 3, route: 'IV-bolus', dose: '1 mg/kg',
            status: 'ok', flagSummary: 'No flags'},
          compoundTier: [{
            compoundId: 'C1', isWholeAnalysis: true, groupValues: [], keyInventory: [], passthrough: [],
            n: 3, route: 'IV-bolus', doseNorm: 'per_kg',
            cmax: {mean: 10, cvPct: 5, n: 3}, tmax: {mean: 1, cvPct: 5, n: 3},
            aucInf: {mean: 100, cvPct: 5, n: 3}, thalf: {mean: 2, cvPct: 5, n: 3},
            cl: {mean: 1, cvPct: 5, n: 3}, vz: {mean: 3, cvPct: 5, n: 3},
            mrt: {mean: 4, cvPct: 5, n: 3}, vss: {mean: 6, cvPct: 5, n: 3},
            fPct: null, pctExtrap: {mean: 5, cvPct: 5, n: 3},
            // A fresh v11 curve WITH the ×/÷ gSD band present.
            curve: {series: [{group: '', time: [0, 1, 2], geoMean: [100, 50, 25],
              geoLower: [90, 45, 22], geoUpper: [110, 55, 28]}]},
            warnings: [],
          }],
        },
      }],
    };
    const loaded = JSON.parse(JSON.stringify(v11)) as NcaStudy;
    const s = loaded.analyses[0].summaryCache!.compoundTier[0].curve!.series[0];
    expect(s.geoLower !== undefined && s.geoUpper !== undefined, true, 'the fresh band survives the round-trip');
    expect(s.geoLower![0], 90, 'lower band value intact');
    expect(s.geoUpper![2], 28, 'upper band value intact');
  });

  test('AC-E3c end-to-end: a multi-dose pooled tier row carries the CROSS_DOSE_POOL warning', async () => {
    // The end-to-end honesty pin (requirement-trace + decision-audit coverage gap): a
    // whole-analysis (ungrouped) analysis whose observations span ≥2 dose levels →
    // buildSummaryCache emits the CROSS_DOSE_POOL {code,message} onto the tier row's warnings
    // (which the study grid surfaces in the Warnings cell). Clean 2-arm doses (10/100).
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV-bolus', 'IV-bolus']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM', ['per_kg', 'per_kg']),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([100, 200])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 20])),
    ]);
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 100, 100])),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([100, 50, 900, 450])),
    ]);
    const obs: ObservationsCtx = {
      df: obsDf, timeName: 'Time', concName: 'Conc', doseName: 'Dose', routeCodeName: null};
    const stub = (): PerProfileResult => ({
      route: nca.ROUTE_IV_BOLUS,
      result: {status: 'ok', provenance: {warnings: []}},
    } as unknown as PerProfileResult);
    const cache = buildSummaryCache(results, [stub(), stub()], recipeWith(null), 'Pooled', obs);
    const row = cache.compoundTier[0];
    expect(row.warnings.some((w) => w.code === 'CROSS_DOSE_POOL'), true,
      'multi-dose pooled tier row carries the CROSS_DOSE_POOL warning (AC-E3c, end-to-end)');
    const msg = row.warnings.find((w) => w.code === 'CROSS_DOSE_POOL')?.message ?? '';
    expect(msg.length > 0 && msg.includes('dose'), true, 'the warning carries a human-readable message');
  });
});

category('NCA Studies: study grid end-to-end from real cache (AD-3)', () => {
  test('AD-3 e2e: real 05 pipeline → buildSummaryCache → buildStudyGrid = 20 rows, IV→CL routing', async () => {
    // Closes the Rule-16 loop end-to-end: real generated pipeline output flows
    // through buildSummaryCache into buildStudyGrid (not a hand-built cache).
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const contract = multicompoundContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, multicompoundScalars(), DEFAULT_RULES);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);
    // Plumb the observations so the v5 inline thumbnail column is exercised
    // end-to-end on the REAL 20-compound fixture (rule 16 — the "20 thumbnail rows"
    // claim, not a synthetic single-row check).
    const obs: ObservationsCtx = {
      df, timeName: 'Time_hr', concName: 'Conc_ng_mL', doseName: 'Dose_mg_kg',
      routeCodeName: null,
    };
    const cache = buildSummaryCache(
      results, outcome.perProfile, recipeWith('CompoundID'), 'Batch', obs);

    const study = {
      schemaVersion: 2, id: 's', label: 'Study',
      createdAt: ISO, updatedAt: ISO, permissions: {view: [], edit: []},
      analyses: [{
        id: 'a1', label: 'Batch', source: {kind: 'tableId', ref: 't'},
        columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols: ['CompoundID'],
        sourceSignature: {kind: 'content-hash', rowCount: 600, mtime: null, hash: 'h', computedAt: ISO},
        summaryCache: cache satisfies SummaryCache,
        resultsBasis: RESULTS_BASIS, // a freshly computed cache (AD-8)
      }],
    } as unknown as NcaStudy;

    const grid = buildStudyGrid(study);
    expect(grid.rowCount, 20, '20 real compound rows projected into the study grid');
    // Rule 16: the cached CL is the per-kg mean — 05 is /kg with no BW, so per kg = PKNCA cl.obs
    // (L/h/kg) × 1000/60; the first compound's mean is in the multicompound per-kg band (0.1–3.3 L/h/kg).
    const cl0 = grid.col('CL Mean')?.get(0) as number;
    expect(cl0 > 0.1 * 1000 / 60 && cl0 < 3.3 * 1000 / 60, true, `CL Mean in mL/min/kg (got ${cl0})`);
    expect(grid.col('CL Mean')?.getTag('units'), 'mL/min/kg', 'fixed per-kg unit');
    // All 20 compounds are IV-bolus → CL group present, CL/F group all-empty → dropped.
    expect(grid.col('CL Mean') !== null, true, 'CL column present (IV cases)');
    expect(grid.col('CL/F Mean'), null, 'CL/F group dropped (no extravascular cases, hide-empty)');
    // S4b: the inline C–t thumbnail DATA_FRAME column is present, with a real curve
    // frame in every one of the 20 rows (the "20 thumbnail rows" DoD, end-to-end).
    const curveCol = grid.col(CURVE_COL);
    expect(curveCol !== null, true, 'C–t thumbnail column present on real fixture-05 output');
    expect(curveCol!.type, DG.COLUMN_TYPE.DATA_FRAME, 'C–t is a DATA_FRAME column');
    let thumbRows = 0;
    for (let i = 0; i < grid.rowCount; i++) {
      const cell = curveCol!.get(i) as DG.DataFrame | null;
      if (cell !== null && cell !== undefined && cell.rowCount > 0) thumbRows++;
    }
    expect(thumbRows, 20, 'all 20 compound rows carry a curve thumbnail frame');
    // Provenance + no internal-id leak on real output.
    expect(grid.col('Analysis')!.get(0), 'Batch', 'provenance column carries the analysis label');
  });

  test('Variant A: Route is offerable and splits BOTH consumers on a multi-route fixture', async () => {
    // The F1 gate re-offers ROUTE iff BOTH the results frame AND the observations frame
    // can split by it. On 03_rat_multiroute (IV/PO/SC): results always carries the
    // canonical ROUTE, and denseDf.col('ROUTE') resolves the source `Route` column
    // (DG col() is case-insensitive) — so the aggregate grid and the Mean±SD chart split
    // by route consistently. Verified against real generated output (rule 16).
    const df = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(df);
    const contract = ratMultirouteContract(df);
    const {results} = runPipelineOrThrow(df, contract, ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    // The exact condition the picker gates ROUTE on:
    expect(results.col('ROUTE') !== null, true, 'results carries the canonical ROUTE (grid can split)');
    expect(df.col('ROUTE') !== null, true,
      'observations resolve a route column case-insensitively (chart can split) — no asymmetry');
    // Aggregate grid: grouping by ROUTE yields ≥2 route groups.
    const agg = aggregateNcaResults(results, ['ROUTE']);
    expect(agg.rowCount >= 2, true, 'aggregate grid splits by route (IV / PO / SC)');
    // Mean±SD chart: grouping the OBSERVATIONS by ROUTE also yields ≥2 groups (consistency).
    const profile = buildMeanSdProfile(df, {timeCol: 'Time_hr', concCol: 'Conc_ng_mL', groupCol: 'ROUTE'});
    const routeGroups = new Set<string>();
    const gCol = profile.col('Group')!;
    for (let i = 0; i < profile.rowCount; i++) routeGroups.add(gCol.getString(i));
    expect(routeGroups.size >= 2, true, 'mean±SD chart splits by route too (both consumers honour ROUTE)');

    // Grid↔chart LABEL consistency: the chart relabels the observations' route CODE via
    // ROUTE_CODE_TO_STRING (the SAME map results-writer uses for the results ROUTE column),
    // so its route labels equal the grid's canonical ROUTE values (e.g. "IV-bolus"), NOT the
    // raw source "IV" that denseDf.col('ROUTE') would otherwise resolve to. Verified against
    // real generated output.
    const codeCol = df.col(DERIVED_ROUTE_CODE_COL)!;
    const chartLabels = new Set<string>();
    for (let i = 0; i < df.rowCount; i++)
      if (!codeCol.isNone(i)) chartLabels.add(ROUTE_CODE_TO_STRING[codeCol.get(i) as number]);
    const gridLabels = new Set<string>();
    const rc = results.col('ROUTE')!;
    for (let i = 0; i < results.rowCount; i++) gridLabels.add(rc.getString(i));
    expect([...chartLabels].sort().join(','), [...gridLabels].sort().join(','),
      'chart route labels == grid canonical ROUTE labels (IV-bolus, not raw IV)');
    // Concretely: the canonical set includes the normalized "IV-bolus", never the raw "IV".
    expect(gridLabels.has('IV-bolus') && !gridLabels.has('IV'), true,
      'canonical labels are normalized (IV-bolus), matching what the chart now shows');
  });

  test('route-split: a multi-route analysis gives each route its OWN curve (no IV+PO pooling)', async () => {
    // 03_rat_multiroute (one compound, IV + PO + SC arms). The tier splits by route;
    // the curve MUST too — an IV-bolus decline curve is never pooled with a PO
    // absorption curve. Verified against real generated output (rule 16).
    const df = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(df);
    const contract = ratMultirouteContract(df);
    const {results, outcome} = runPipelineOrThrow(
      df, contract, ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);

    // Route-split key: the pipeline's own per-row derived RouteCode column
    // (DERIVED_ROUTE_CODE_COL, added by setupMultirouteRoute) — the SAME code that
    // set each tier row's ROUTE, read directly (no re-normalization). This is the
    // identical mechanism NONMEM uses (its transform encodes the same column).
    const obs: ObservationsCtx = {
      df, timeName: 'Time_hr', concName: 'Conc_ng_mL', doseName: 'Dose_mg_kg',
      routeCodeName: DERIVED_ROUTE_CODE_COL,
    };
    const cache = buildSummaryCache(results, outcome.perProfile, recipeWith(null), 'Multiroute', obs);

    const iv = cache.compoundTier.find((r) => r.route === 'IV-bolus');
    const po = cache.compoundTier.find((r) => r.route === 'PO');
    expect(iv !== undefined && po !== undefined, true, 'IV-bolus + PO route rows present');
    expect(iv!.curve !== null, true, 'IV-bolus row has its own curve');
    expect(po!.curve !== null, true, 'PO row has its own curve');
    // The two curves are route-distinct: an IV-bolus profile starts at its peak
    // (t≈0 is Cmax) while a PO profile absorbs up from a lower initial value, so the
    // earliest geomean of the pooled-together case cannot equal both. Assert the
    // first points differ — proving the curves were NOT pooled into one shared object.
    const ivFirst = iv!.curve!.series[0].geoMean[0];
    const poFirst = po!.curve!.series[0].geoMean[0];
    expect(ivFirst !== poFirst, true, 'IV and PO curves are distinct (not a shared pooled curve)');
  });

  test('route-split at an OVERLAPPING dose (the NONMEM-shape defect): routes never pool', async () => {
    // The residual gap the architect flagged: a compound dosed via >1 route at the
    // SAME dose value (real case: 03_rat_multiroute_nonmem shares dose 5.0 across
    // IV/PO/SC). With a compound␁dose key that would pool all routes at dose 5.0
    // into one series. The derived RouteCode column splits them. Synthetic obs frame
    // carrying the per-row RouteCode (the shape both NONMEM + standard produce) →
    // deterministic proof the route dimension separates same-dose routes.
    const obsDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound', ['C1', 'C1', 'C1', 'C1']),
      DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL,
        Int32Array.from([nca.ROUTE_IV_BOLUS, nca.ROUTE_IV_BOLUS, nca.ROUTE_PO, nca.ROUTE_PO])),
      // SAME dose (5) for both routes — the overlap that pooled them before the fix.
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Dose', ['5', '5', '5', '5']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([300, 150, 80, 120])), // IV declines; PO absorbs
    ]);
    const curves = buildCurvesByCompound(
      {df: obsDf, timeName: 'Time', concName: 'Conc', doseName: 'Dose',
        routeCodeName: DERIVED_ROUTE_CODE_COL}, 'Compound', 'X');
    expect(curves.routeKeyed, true, 'derived RouteCode column → route-keyed curves');
    const ivKey = `C1${String.fromCharCode(1)}IV-bolus`;
    const poKey = `C1${String.fromCharCode(1)}PO`;
    const ivCurve = curves.byKey.get(ivKey);
    const poCurve = curves.byKey.get(poKey);
    expect(ivCurve !== undefined && poCurve !== undefined, true,
      'IV-bolus and PO get SEPARATE curve keys at the SAME dose (not one pooled series)');
    // IV t=1 geomean = 300; PO t=1 geomean = 80 — proving no cross-route pooling.
    // (Geomean of a single value = exp(ln x) carries fp noise → tolerance compare.)
    expectFloat(ivCurve!.series[0].geoMean[0], 300, 1e-6, 'IV curve is the IV-only profile');
    expectFloat(poCurve!.series[0].geoMean[0], 80, 1e-6, 'PO curve is the PO-only profile');
  });
});

// ── DATA-GAP-32: multi-analysis cross-cased saved-study fixture (E2/E3/E4 real round-trip) ──
// Rule-16 review-gate fixture (study-grid-grouping-columns DG-1, user decision Variant-2
// 2026-07-13): FOUR analyses built from REAL pipeline output (not hand-typed summaryCache),
// assembled into one NcaStudy, JSON-serialized and reloaded (the summary cache is an EMS json
// column), then projected by the REAL buildStudyGrid. Exercises the cross-analysis
// merge + counted join + pooled theoph thumbnail against genuine generated output + the
// storage round-trip — the surface the hand-built projection tests could not reach.
category('NCA Studies: DATA-GAP-32 multi-analysis cross-cased saved study (E2/E3/E4)', () => {
  test('case-merge (Sex/sex) + counted join + absent + pooled theoph, through the JSON round-trip', async () => {
    // Run the proven multiroute pipeline on a df → real results + observations for the curve.
    const runRat = (df: DG.DataFrame): {
      results: DG.DataFrame; perProfile: ReadonlyArray<PerProfileResult>; obs: ObservationsCtx;
    } => {
      const opts = setupMultirouteRoute(df);
      const contract = ratMultirouteContract(df);
      const {results, outcome} = runPipelineOrThrow(df, contract, ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
      appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);
      const obs: ObservationsCtx = {
        df, timeName: 'Time_hr', concName: 'Conc_ng_mL', doseName: 'Dose_mg_kg',
        routeCodeName: DERIVED_ROUTE_CODE_COL};
      return {results, perProfile: outcome.perProfile, obs};
    };

    // A: 03_rat_multiroute grouped by Sex. C: SAME real results grouped by Strain (Sex present-
    // but-ungrouped → counted join). B: a FRESH load with Sex renamed to lowercase 'sex',
    // grouped by 'sex' → must MERGE with A's Sex column on the typographic canonical id.
    const ratA = await loadFixture('03_rat_multiroute.csv');
    const runA = runRat(ratA);
    const cacheA = buildSummaryCache(runA.results, runA.perProfile, recipeWith('Sex'), 'A_Sex', runA.obs);
    const cacheC = buildSummaryCache(runA.results, runA.perProfile, recipeWith('Strain'), 'C_Strain', runA.obs);

    const ratB = await loadFixture('03_rat_multiroute.csv');
    ratB.col('Sex')!.name = 'sex'; // a differently-CASED grouping column (real generated output)
    const runB = runRat(ratB);
    const cacheB = buildSummaryCache(runB.results, runB.perProfile, recipeWith('sex'), 'B_sex', runB.obs);

    // D: theoph whole-analysis (no Sex column) → Sex cell '—' + ONE pooled thumbnail series.
    const theophDf = await loadFixture('theoph.csv');
    const {results: tRes, outcome: tOut} =
      runPipelineOrThrow(theophDf, theophContract(theophDf), theophScalars(), DEFAULT_RULES);
    appendGroupingColumns(tRes, theophDf, tOut.perProfile.map((p) => p.index), theophContract(theophDf));
    const theophObs: ObservationsCtx = {
      df: theophDf, timeName: 'Time', concName: 'conc', doseName: 'Dose', routeCodeName: null};
    const cacheD = buildSummaryCache(tRes, tOut.perProfile, recipeWith(null), 'D_theoph', theophObs);

    const mkEntry = (
      id: string, label: string, groupByCols: string[], cache: SummaryCache,
    ): Record<string, unknown> => ({
      id, label, source: {kind: 'tableId', ref: 't'},
      columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols, passthroughCols: [],
      sourceSignature: {kind: 'content-hash', rowCount: 6, mtime: null, hash: 'h', computedAt: ISO},
      summaryCache: cache,
    });
    const study = {
      id: 'dg32', label: 'DATA-GAP-32', kind: 'study',
      createdAt: ISO, updatedAt: ISO, access: {canEdit: true, canDelete: true, canShare: true},
      analyses: [
        mkEntry('a', 'A_Sex', ['Sex'], cacheA),
        mkEntry('b', 'B_sex', ['sex'], cacheB),
        mkEntry('c', 'C_Strain', ['Strain'], cacheC),
        mkEntry('d', 'D_theoph', [], cacheD),
      ],
    };
    // The storage round-trip (the summary cache is an EMS json column). NaN means serialize to null
    // and read back as blank (documented).
    const loaded = JSON.parse(JSON.stringify(study)) as NcaStudy;
    const grid = buildStudyGrid(loaded);

    // AC-E2a: Sex (analysis A) and sex (analysis B) MERGE into ONE column (canonical 'sex').
    const sexFamily = grid.columns.names().filter((nm) => normalizeGroupKey(nm) === 'sex');
    expect(sexFamily.length, 1, 'Sex/sex from two analyses merge into ONE column (AC-E2a, real round-trip)');
    const sexCol = grid.col(sexFamily[0])!;
    const analysisCol = grid.col('Analysis')!;
    const cellsFor = (label: string): string[] => {
      const out: string[] = [];
      for (let i = 0; i < grid.rowCount; i++)
        if (analysisCol.get(i) === label) out.push(String(sexCol.get(i)));
      return out;
    };
    // A + B rows carry plain grouped M/F under the SAME merged column.
    const aCells = cellsFor('A_Sex'); const bCells = cellsFor('B_sex');
    expect(aCells.length > 0 && aCells.every((v) => v === 'M' || v === 'F'), true,
      'A grouped-by-Sex cells are plain M/F');
    expect(bCells.length > 0 && bCells.every((v) => v === 'M' || v === 'F'), true,
      'B grouped-by-sex cells merge under the SAME Sex column (case-merge)');
    // AC-E4a: C (grouped by Strain, Sex present-but-ungrouped) → counted composition 'X (n); …'.
    const cCells = cellsFor('C_Strain');
    expect(cCells.length > 0 && cCells.every((v) => /\(\d+\)/.test(v)), true,
      'C Sex cell = counted composition (present-but-ungrouped, AC-E4a)');
    // AC-E4b: D (theoph, no Sex column at all) → '—'.
    const dCells = cellsFor('D_theoph');
    expect(dCells.length > 0 && dCells.every((v) => v === NO_GROUP), true,
      'theoph has no Sex → — (genuinely absent, AC-E4b)');

    // AC-E3a: theoph's real thumbnail is ONE pooled geomean series (whole-analysis).
    const dRow = loaded.analyses.find((a) => a.label === 'D_theoph')!.summaryCache!.compoundTier[0];
    expect(dRow.curve !== null && dRow.curve.series.length === 1, true,
      'theoph → ONE pooled thumbnail series after real round-trip (AC-E3a)');
  });
});
