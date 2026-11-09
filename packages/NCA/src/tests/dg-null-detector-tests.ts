/**
 * The permanent raw-NaN drift guard over REAL pipeline output
 * (nca-dg-native-missing-values F4, AC-12 / AC-12b / AC-13 / AC-14) plus the
 * F3 aggregator-deletion evidence (AC-9) — all anchored on the DG-1 mixed-fit
 * fixture.
 *
 * **Why a MIXED fixture is a precondition, not a nicety (DG-1).** The migration's
 * whole payoff is how DG-native consumers behave on a frame where SOME subjects
 * fit and SOME fail. On a uniformly-fit frame every assertion here passes
 * vacuously and the change looks like a no-op. `indometh` at a raised λz
 * `minRSquared` produces exactly that shape, and {@link assertMixed} FAILS if it
 * ever stops doing so — a fixture that silently drifts to all-fit would turn this
 * whole file green while proving nothing.
 *
 * **Rule 16.** Every claim below is asserted against the frame the real pipeline
 * actually produced, never against a mirror of a writer's source.
 */

import {category, test, expect, before} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import {nca} from '@datagrok-libraries/sci-comp';

import {DEFAULT_RULES, withLambdaZ} from '../adapter/defaults';
import {
  scanNumericSentinels, describeSentinelEntries, toPlatformNull, log10OrNull,
} from '../adapter/adapter-utils';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {buildMeanSdProfile} from '../adapter/mean-sd-profile';
import {computeBioavailability} from '../adapter/bioavailability';
import {computeDoseProportionality} from '../adapter/dose-proportionality';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {
  loadFixture, indomethContract, indomethScalars, ratScalars,
  ratMultirouteContract, setupMultirouteRoute, runPipelineOrThrow,
} from './test-utils';

const TIMEOUT = 60000;

/**
 * **DG-1 — the mixed fit/unfit fixture.** `indometh` with the λz acceptance
 * raised to adj-R² ≥ 0.87 (default is 0.85). Live QA on 2026-07-30 observed rows
 * 3/4/5 failing and 1/2/6 fitting at exactly this threshold; the assertion below
 * re-derives that from the produced frame rather than trusting the number.
 */
const MIXED_FIT_MIN_R2 = 0.87;

let mixed: DG.DataFrame;
let uniform: DG.DataFrame;

/** Rows whose λz did not fit (LAMBDA_Z absent) / did fit, read through `isNone`. */
function fitSplit(results: DG.DataFrame): {fit: number[]; unfit: number[]} {
  const lz = results.col('LAMBDA_Z')!;
  const fit: number[] = []; const unfit: number[] = [];
  for (let i = 0; i < results.rowCount; i++) (lz.isNone(i) ? unfit : fit).push(i);
  return {fit, unfit};
}

/** DG-1's own precondition: the fixture must genuinely be MIXED. */
function assertMixed(results: DG.DataFrame, tag: string): {fit: number[]; unfit: number[]} {
  const s = fitSplit(results);
  expect(s.unfit.length > 0, true,
    `${tag}: DG-1 precondition violated — no UNFIT profile, every assertion downstream is vacuous`);
  expect(s.fit.length > 0, true,
    `${tag}: DG-1 precondition violated — no FIT profile, the frame is uniformly empty`);
  return s;
}

/** Fail on any raw NaN; report (do not fail on) any ±Infinity. AC-12 / AC-12b. */
function assertNoRawNan(df: DG.DataFrame, frameName: string): void {
  const entries = scanNumericSentinels(df);
  const nans = entries.filter((e) => e.kind === 'raw-nan');
  expect(nans.length, 0,
    `${frameName}: raw IEEE NaN in a DG column — invisible to isNone, counted as PRESENT by ` +
    `group-by/filter/export, and it poisons the column's statistics. ` +
    `Offenders: ${describeSentinelEntries(frameName, nans)}`);
  // AC-12b — report-only. `f64Column` deliberately passes ±Infinity through
  // (AC-3); the `ivRef === 0 → +Infinity` defect is filed as its own gap. But the
  // detector must not CERTIFY a frame as clean while it holds one, and the
  // migration makes that worse rather than neutral: a co-resident NaN used to
  // dominate every statistic and mask the Inf; now the NaN is excluded and the
  // Inf is the sole surviving poisoner.
  const infs = entries.filter((e) => e.kind === 'non-finite');
  if (infs.length > 0)
    console.log(`[AC-12b report-only] ${describeSentinelEntries(frameName, infs)}`);
}

category('DG null detector: Arm 1 over real pipeline output (AC-12/12b/13)', () => {
  before(async () => {
    const df = await loadFixture('indometh.csv');
    ({results: mixed} = runPipelineOrThrow(
      df, indomethContract(df), indomethScalars(),
      withLambdaZ(DEFAULT_RULES, {minRSquared: MIXED_FIT_MIN_R2})));
    const df2 = await loadFixture('indometh.csv');
    ({results: uniform} = runPipelineOrThrow(
      df2, indomethContract(df2), indomethScalars(), DEFAULT_RULES));
  });

  test('DG-1: the fixture is genuinely MIXED, and unfit rows read as MISSING not present', async () => {
    const {fit, unfit} = assertMixed(mixed, 'indometh@0.87');
    // The whole point of the migration, asserted on the produced frame: an unfit
    // λz parameter must read `isNone` — NOT merely "be NaN". Before the boundary
    // existed, `isNone` was FALSE here and every DG-native consumer read the cell
    // as a real value.
    for (const p of ['LAMBDA_Z', 'HALF_LIFE', 'AUC_INF', 'PCT_EXTRAP', 'LAMBDAZ_ADJ_R2']) {
      const col = mixed.col(p)!;
      for (const i of unfit)
        expect(col.isNone(i), true, `${p}[${i}]: an unfit profile's parameter must read isNone`);
    }
    for (const i of fit)
      expect(mixed.col('LAMBDA_Z')!.isNone(i), false, `LAMBDA_Z[${i}]: a fit profile keeps its value`);
  }, {timeout: TIMEOUT});

  test('DG-1: column statistics are no longer poisoned by the unfit rows', async () => {
    // Measured platform behaviour (probe PR-3): ONE raw NaN sends min/max/avg/
    // sum/stdev of the whole column to NaN, while FLOAT_NULL is excluded cleanly.
    // This is the user-visible payoff and it is asserted on real output.
    const {unfit} = assertMixed(mixed, 'indometh@0.87');
    const stats = mixed.col('HALF_LIFE')!.stats;
    expect(Number.isFinite(stats.avg), true,
      'HALF_LIFE mean must be a real number despite unfit profiles in the column');
    expect(stats.missingValueCount, unfit.length,
      'every unfit profile is counted as MISSING, not as a value');
    expect(stats.valueCount, mixed.rowCount - unfit.length, 'value count excludes the unfit rows');
  }, {timeout: TIMEOUT});

  test('AC-8: is the MRT/VSS fragility-tint break REACHABLE? (measure, then decide)', async () => {
    // The plan requires establishing reachability BEFORE shipping a test, so that a
    // criterion that cannot fail is never mistaken for one that passed.
    //
    // The three unguarded sibling reads at `results-view.ts:716-718` were masked by
    // the `own` guard at `:713-714`: the tint only evaluates when MRT/VSS is itself
    // present. So the break is observable ONLY on a row where MRT/VSS is PRESENT
    // while LAMBDAZ_ADJ_R2 is NULL. This test measures whether the pipeline can
    // produce such a row, on a frame that genuinely contains failed λz fits.
    const {unfit} = assertMixed(mixed, 'indometh@0.87');
    const adjR2 = mixed.col('LAMBDAZ_ADJ_R2')!;
    let unmasked = 0;
    for (const p of ['MRT', 'VSS']) {
      const col = mixed.col(p)!;
      for (let i = 0; i < mixed.rowCount; i++)
        if (!col.isNone(i) && adjR2.isNone(i)) unmasked++;
    }
    // OBSERVED 2026-08-03 on indometh @ minRSquared 0.87: 0 such rows. MRT and Vss
    // are both derived from AUCinf/AUMCinf, which require the λz extrapolation — so
    // when λz fails they are absent too, and the `own` guard returns early. The
    // dependency is structural, not incidental.
    expect(unmasked, 0,
      `MRT/VSS present while LAMBDAZ_ADJ_R2 is null on ${unfit.length} unfit rows — if this ` +
      `is ever non-zero the fragility tint becomes OBSERVABLE and needs a behavioural ` +
      `assertion, not just the isNone-first read`);
    // VERDICT: the path is UNREACHABLE through the current derivation chain, so the
    // `numAt` fix at those three reads is DEFENSIVE — it removes a correctness
    // dependency on a sibling column's nullity rather than fixing a visible symptom.
    // Recorded here and in the read-site verdict list instead of shipping a
    // behavioural test that could only ever pass vacuously.
  }, {timeout: TIMEOUT});

  test('AC-12/AC-13: no raw NaN in any result frame — per-subject writer + aggregator', async () => {
    assertNoRawNan(mixed, 'results-writer(indometh@0.87)');
    assertNoRawNan(uniform, 'results-writer(indometh@default)');
    assertNoRawNan(aggregateNcaResults(mixed, 'ROUTE'), 'results-aggregator(indometh@0.87)');
  }, {timeout: TIMEOUT});

  test('AC-13: no raw NaN — mean-sd-profile, bioavailability, dose-proportionality', async () => {
    // mean-sd-profile over a real observations frame.
    const obs = await loadFixture('indometh.csv');
    assertNoRawNan(
      buildMeanSdProfile(obs, {timeCol: 'time', concCol: 'conc', groupCol: 'Subject'}),
      'mean-sd-profile(indometh)');

    // bioavailability + dose-proportionality over the rat multi-route frame — the
    // only fixture with both an IV reference and extravascular arms.
    const rat = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(rat);
    const {results, outcome} = runPipelineOrThrow(
      rat, ratMultirouteContract(rat), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    appendGroupingColumns(
      results, rat, outcome.perProfile.map((p) => p.index), ratMultirouteContract(rat));
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['Sex']});
    assertNoRawNan(table, 'bioavailability(rat-multiroute)');
    const dp = computeDoseProportionality(results, {doseCol: 'Dose_mg_kg'});
    assertNoRawNan(dp.table, 'dose-proportionality(rat-multiroute)');
    assertNoRawNan(dp.dnTable, 'dose-proportionality-dn(rat-multiroute)');
  }, {timeout: TIMEOUT});
});

category('DG null detector: Arm 2 — Float32 narrowing regression pin (AC-14)', () => {
  // `/ops:bug-stress` closed the site inventory at exactly 3 narrowing sites and
  // all 3 are now routed through these two shared guards, so pinning the guards
  // pins the sites. Previously each site hand-spelled the predicate — which is how
  // one of them (`buildFacetDataFrames`) shipped unguarded until `39b693e`.
  test('AC-14: toPlatformNull maps BOTH sentinels to the platform null', async () => {
    expect(toPlatformNull(NaN) === DG.FLOAT_NULL, true, 'raw NaN → FLOAT_NULL (load-bearing half)');
    expect(toPlatformNull(DG.FLOAT_NULL) === DG.FLOAT_NULL, true, 'FLOAT_NULL → FLOAT_NULL');
    expect(toPlatformNull(7.5), 7.5, 'a real value passes through');
    expect(toPlatformNull(0), 0, 'a real zero is a value, not a null');
    // Narrowed to Float32 the sentinel is preserved exactly (probe PR-10) and the
    // resulting cell still reads isNone — the property the guard relies on.
    const narrowed = new Float32Array(1);
    narrowed[0] = toPlatformNull(NaN);
    expect(narrowed[0] === DG.FLOAT_NULL, true, 'FLOAT_NULL survives Float32 narrowing exactly');
    expect(DG.Column.fromFloat32Array('c', narrowed).isNone(0), true,
      'and the narrowed cell reads isNone');
  });

  test('AC-14: log10OrNull refuses the sentinel, which PASSES a naive `> 0`', async () => {
    // The separately load-bearing half: FLOAT_NULL is a tiny POSITIVE number, so
    // an unguarded `v > 0 ? Math.log10(v) : null` stores a REAL −33.6 and silently
    // rescales the axis 33 decades.
    expect(DG.FLOAT_NULL > 0, true, 'the sentinel is positive — this is why the check is needed');
    expect(log10OrNull(DG.FLOAT_NULL) === DG.FLOAT_NULL, true, 'sentinel → null, not log10(sentinel)');
    expect(log10OrNull(NaN) === DG.FLOAT_NULL, true, 'raw NaN → null');
    expect(log10OrNull(0) === DG.FLOAT_NULL, true, 'zero (BLQ) → null, never −Infinity');
    expect(log10OrNull(-1) === DG.FLOAT_NULL, true, 'negative → null');
    expect(log10OrNull(100), 2, 'a real positive value is log10-ed');
  });
});

category('DG null detector: F3 aggregator deletion evidence (AC-9/AC-10)', () => {
  test('AC-9: every aggregated param reaches groupBy already converted (coverage/ordering)', async () => {
    // **What this can and cannot prove.** The clone that F3 deleted applied the
    // predicate `Number.isNaN(v) ? FLOAT_NULL : v` — BYTE-IDENTICAL to f64Column's.
    // So before/after equality holds BY CONSTRUCTION and proves nothing semantic
    // (plan review R1, F-4). Its real discriminating power is COVERAGE and
    // ORDERING: it fails if a parameter column reaches the aggregator without
    // having passed through f64Column. That is what is asserted here, directly.
    const {unfit} = assertMixed(mixed, 'indometh@0.87');
    const params = ['CMAX', 'TMAX', 'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP',
      'LAMBDA_Z', 'HALF_LIFE', 'CL', 'VZ', 'AUMC_INF', 'MRT', 'VSS'];
    for (const p of params) {
      const col = mixed.col(p);
      if (col === null) continue;
      const raw = col.getRawData() as Float64Array;
      for (let i = 0; i < raw.length; i++) {
        expect(Number.isNaN(raw[i]), false,
          `${p}[${i}] reaches the aggregator as a RAW NaN — a writer bypassed f64Column, ` +
          `which the deleted shadow clone used to hide`);
      }
    }
    // And the na.rm actually excluded something — otherwise the whole criterion is
    // vacuous on this frame.
    const agg = aggregateNcaResults(mixed, 'ROUTE');
    const nHalfLife = agg.col('N_HALF_LIFE')!;
    const nCohort = agg.col('N')!;
    expect(agg.rowCount > 0, true, 'aggregate produced at least one group');
    expect((nHalfLife.get(0) as number) < (nCohort.get(0) as number), true,
      `N_HALF_LIFE must be below the cohort N on a mixed frame (${unfit.length} unfit rows) — ` +
      `if equal, the unfit profiles were counted as contributors`);
  }, {timeout: TIMEOUT});

  test('AC-9: MEAN/SD/N are identical with and without the deleted clone', async () => {
    // The clone, re-applied by hand to a copy of the source frame. Post-migration
    // it is a no-op, so the two aggregates must agree cell-for-cell. This is the
    // "true by construction" half — kept because it fails loudly if the writer's
    // conversion and the aggregator ever run in the wrong order.
    const cols: DG.Column[] = [];
    for (const c of mixed.columns as Iterable<DG.Column>) {
      if (c.type !== DG.COLUMN_TYPE.FLOAT) {cols.push(c); continue;}
      const raw = c.getRawData() as Float64Array;
      const cleaned = new Float64Array(raw.length);
      for (let i = 0; i < raw.length; i++)
        cleaned[i] = Number.isNaN(raw[i]) ? DG.FLOAT_NULL : raw[i];
      cols.push(DG.Column.fromFloat64Array(c.name, cleaned));
    }
    const legacy = aggregateNcaResults(DG.DataFrame.fromColumns(cols), 'ROUTE');
    const current = aggregateNcaResults(mixed, 'ROUTE');
    expect(current.rowCount, legacy.rowCount, 'same group count');
    for (const p of ['CMAX', 'AUC_INF', 'HALF_LIFE', 'CL', 'VZ']) {
      for (const stat of [`MEAN_${p}`, `SD_${p}`, `N_${p}`]) {
        const a = current.col(stat); const b = legacy.col(stat);
        expect(a !== null && b !== null, true, `${stat} present in both aggregates`);
        for (let i = 0; i < current.rowCount; i++) {
          expect(a!.isNone(i), b!.isNone(i), `${stat}[${i}] missingness identical`);
          if (!a!.isNone(i))
            expect(a!.get(i) as number, b!.get(i) as number, `${stat}[${i}] value identical`);
        }
      }
    }
    expect((current.col('N')!.get(0) as number), (legacy.col('N')!.get(0) as number),
      'cohort N identical');
  }, {timeout: TIMEOUT});

  test('AC-10: the %CV guard keeps isNone — Number.isFinite would NOT catch the sentinel', async () => {
    const agg = aggregateNcaResults(mixed, 'ROUTE');
    const cv = agg.col('CV_PCT_HALF_LIFE');
    expect(cv !== null, true, 'CV_PCT_HALF_LIFE emitted');
    // F1 scope: the CV column is built by f64Column, so an unestimable %CV is the
    // platform null and NOT the literal text "NaN" in the grid.
    const raw = cv!.getRawData() as Float64Array;
    for (let i = 0; i < raw.length; i++)
      expect(Number.isNaN(raw[i]), false, `CV_PCT_HALF_LIFE[${i}] must not be a raw NaN`);
    // And the reason the guard cannot be "simplified" to isFinite:
    expect(Number.isFinite(DG.FLOAT_NULL), true,
      'FLOAT_NULL is FINITE — an isFinite-based %CV guard would divide by the sentinel');
  }, {timeout: TIMEOUT});
});
