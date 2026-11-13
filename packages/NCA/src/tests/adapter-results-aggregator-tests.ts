/**
 * Adapter — results aggregator (SPIKE: per-compound-aggregation).
 *
 * @see src/adapter/results-aggregator.ts (implementation under test)
 *
 * Proves two spike unknowns:
 *  1. native DG `groupBy().aggregate()` produces a correct per-group
 *     ARITHMETIC summary (mean / SD / CV% / N) over NCA parameter columns;
 *  2. how DG `avg` treats a NaN parameter in a group (partial subjects) —
 *     the result informs the na.rm policy decision for the formalized feature.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';

import {
  aggregateNcaResults, resolveAggregateGroupKeys, DEFAULT_AGG_PARAMS,
  blankMixedBasisFamily, IV_EV_INCOMPARABLE_FAMILY,
} from '../adapter/results-aggregator';
import {f64Column} from '../adapter/adapter-utils';

interface MockRow {readonly group: string; readonly cmax: number; readonly auc: number}

/** A stand-in results frame built the way the PRODUCTION writers build one — via
 *  `f64Column`, so an unfit parameter arrives as the platform null rather than as a
 *  raw IEEE NaN (nca-dg-native-missing-values F1/F5).
 *
 *  Load-bearing, not cosmetic. The aggregator used to clone every param column and
 *  convert NaN to FLOAT_NULL itself; F3 deleted that clone because the writers now
 *  emit the sentinel at construction. A fixture still handing the aggregator a RAW
 *  NaN would be testing a frame shape production can no longer produce — and would
 *  fail, correctly, because DG `avg` does not skip a raw NaN. Tests below still
 *  pass `NaN` at the call site: that is the in-package spelling of "missing", and
 *  `f64Column` is exactly the boundary that translates it. */
function mockResults(rows: ReadonlyArray<MockRow>): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'COMPOUND', rows.map((r) => r.group)),
    f64Column('CMAX', rows.map((r) => r.cmax)),
    f64Column('AUC_INF', rows.map((r) => r.auc)),
  ]);
}

function rowByGroup(agg: DG.DataFrame, value: string, groupCol = 'COMPOUND'): number {
  const col = agg.col(groupCol)!;
  for (let i = 0; i < agg.rowCount; i++)
    if (col.get(i) === value) return i;
  throw new Error(`group '${value}' not found in aggregate`);
}

category('Adapter: results aggregator (spike)', () => {
  test('arithmetic mean / SD / CV% / N per group on a clean 2-compound set', async () => {
    const df = mockResults([
      {group: 'A', cmax: 10, auc: 100},
      {group: 'A', cmax: 20, auc: 200},
      {group: 'A', cmax: 30, auc: 300},
      {group: 'B', cmax: 5, auc: 50},
      {group: 'B', cmax: 7, auc: 70},
    ]);

    const agg = aggregateNcaResults(df, 'COMPOUND', {params: ['CMAX', 'AUC_INF']});
    expect(agg.rowCount, 2, 'two groups');

    const a = rowByGroup(agg, 'A');
    const b = rowByGroup(agg, 'B');

    // Means + counts — the core mechanism proof.
    expectFloat(agg.col('MEAN_CMAX')!.get(a) as number, 20, 1e-9, 'A mean Cmax');
    expectFloat(agg.col('MEAN_AUC_INF')!.get(a) as number, 200, 1e-9, 'A mean AUC');
    expect(agg.col('N')!.get(a), 3, 'A N=3');
    expectFloat(agg.col('MEAN_CMAX')!.get(b) as number, 6, 1e-9, 'B mean Cmax');
    expect(agg.col('N')!.get(b), 2, 'B N=2');

    // SD + CV%. A Cmax = [10,20,30]: sample SD (n-1) = 10 → CV% = 50.
    // SPIKE FINDING: if DG uses population SD this asserts ~8.165 instead —
    // adjust + document the convention when observed.
    expectFloat(agg.col('SD_CMAX')!.get(a) as number, 10, 1e-6, 'A sample SD Cmax');
    expectFloat(agg.col('CV_PCT_CMAX')!.get(a) as number, 50, 1e-6, 'A CV% Cmax');
  });

  test('na.rm semantics: an unfit parameter yields mean-of-valid + per-param N', async () => {
    // The compound has 3 subjects; one is partial on AUC_INF. That value is the
    // platform null by the time it reaches `groupBy` — written that way at
    // CONSTRUCTION by `f64Column` (F1), where it used to be converted by a clone
    // inside this module (F3 deleted it). DG `avg` skips the null (na.rm), giving
    // the mean of the 2 valid AUC values. CMAX is clean → unaffected. `N` is the
    // group row count (3); `N_AUC_INF` is the valid-fit count (2).
    const df = mockResults([
      {group: 'C', cmax: 10, auc: 100},
      {group: 'C', cmax: 20, auc: NaN}, // partial subject: AUC_INF unfit
      {group: 'C', cmax: 30, auc: 200},
    ]);

    const agg = aggregateNcaResults(df, 'COMPOUND', {params: ['CMAX', 'AUC_INF']});
    expect(agg.rowCount, 1, 'one group');

    // Clean param unaffected: mean of [10, 20, 30] = 20.
    expectFloat(agg.col('MEAN_CMAX')!.get(0) as number, 20, 1e-9, 'Cmax mean unaffected by AUC NaN');

    // Partial-aware aggregation: mean of [100, _, 200] = 150 (not null).
    expect(agg.col('MEAN_AUC_INF')!.isNone(0), false,
      'MEAN_AUC_INF must be a real value (na.rm), not null');
    expectFloat(agg.col('MEAN_AUC_INF')!.get(0) as number, 150, 1e-6,
      'AUC_INF mean = mean of the 2 valid values [100, 200] (NaN dropped per na.rm)');

    // N family — DMPK convention: N is cohort size, N_<P> is contributors.
    expect(agg.col('N')!.get(0), 3, 'N = group row count (3 subjects)');
    expect(agg.col('N_CMAX')!.get(0), 3, 'N_CMAX = 3 (all subjects fit Cmax)');
    expect(agg.col('N_AUC_INF')!.get(0), 2, 'N_AUC_INF = 2 (only 2 of 3 subjects fit AUC_INF)');

    // Source-frame immutability. The aggregator now groups the SOURCE frame
    // directly (no clone), so this matters MORE than it did, not less: nothing here
    // may write back into the input.
    expect(df.col('AUC_INF')!.isNone(1), true,
      'aggregator must NOT mutate the source frame — the null cell is still null');
    expect(df.col('AUC_INF')!.get(0) as number, 100,
      'aggregator must NOT mutate the source frame — real values untouched');
  });

  test('AC-C1 (R2-NF-2): DEFAULT_AGG_PARAMS contains the FR-200 moment headlines + the companions', async () => {
    // Guards the 9→12 expansion (AUMC_INF/MRT/VSS) and the 12→15 one (the absolute companions,
    // nca-cl-bw-normalised-units AD-4) so a future drop is caught.
    expect(DEFAULT_AGG_PARAMS.length, 15, 'DEFAULT_AGG_PARAMS must have 15 entries');
    for (const p of ['AUMC_INF', 'MRT', 'VSS', 'CL_ABS', 'VZ_ABS', 'VSS_ABS']) {
      expect(DEFAULT_AGG_PARAMS.includes(p), true,
        `DEFAULT_AGG_PARAMS must include ${p}`);
    }
    // Per AD-13 these are intentionally NOT aggregated by default.
    for (const p of ['AUMC_LAST', 'TLAG', 'PCT_AUMC_EXTRAP']) {
      expect(DEFAULT_AGG_PARAMS.includes(p), false,
        `DEFAULT_AGG_PARAMS must NOT include ${p} (AD-13)`);
    }
  });

  test('AC-C2: Mean Vss pools IV rows only (non-IV Vss = NaN → na.rm skipped)', async () => {
    // Mixed IV + EV frame: the route gate writes NaN Vss on the EV row, so
    // DG `avg` drops it (na.rm) — Mean Vss is the IV-only mean for free.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'COMPOUND', ['X', 'X', 'X']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM',
        ['absolute', 'absolute', 'absolute']),
      f64Column('VSS', [10, 20, NaN]),
    ]);
    const agg = aggregateNcaResults(df, 'COMPOUND', {params: ['VSS']});
    const r = rowByGroup(agg, 'X');
    expectFloat(agg.col('MEAN_VSS')!.get(r) as number, 15, 1e-9,
      'Mean Vss = mean of the 2 IV rows [10, 20]; the empty EV row is dropped');
    expect(agg.col('N_VSS')!.get(r), 2, 'N_VSS = 2 (IV rows only; EV Vss is empty)');
  });

  test('throws when groupColName is not a column of results', async () => {
    const df = mockResults([{group: 'A', cmax: 10, auc: 100}]);
    let threw = false;
    let message = '';
    try {
      aggregateNcaResults(df, 'NOT_A_COLUMN', {params: ['CMAX']});
    } catch (e) {
      threw = true;
      message = e instanceof Error ? e.message : String(e);
    }
    expect(threw, true, 'aggregateNcaResults must throw on a missing group column');
    expect(message.includes('NOT_A_COLUMN'), true,
      `error message must name the missing column; got: ${message}`);
  });

  // ---------------------------------------------------------------------
  // PLATFORM RECORD (AC-11). These two tests are the reason the sentinel migration
  // exists. RE-POINTED, never deleted: they document DG's MEASURED behaviour, which
  // is what makes `f64Column` load-bearing rather than a stylistic preference. If
  // the platform ever changes here, these fail first.
  // ---------------------------------------------------------------------
  test('PLATFORM RECORD: DG `avg` PROPAGATES a raw IEEE NaN — one cell nulls the group', async () => {
    // The original spike finding, now asserted directly instead of implied. This is
    // the defect: a raw NaN is NOT `none` to DG, so it is counted as a contributor
    // AND destroys the statistic.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'GRP', ['C', 'C', 'C']),
      DG.Column.fromFloat64Array('AUC', Float64Array.from([100, NaN, 200])),
    ]);
    expect(df.col('AUC')!.isNone(1), false,
      'a raw NaN is NOT null to DG — the whole premise of the migration');
    const agg = df.groupBy(['GRP']).avg('AUC', 'MEAN_AUC').valueCount('AUC', 'N_AUC').aggregate();
    // The group mean is DESTROYED — DG nulls it rather than returning mean-of-valid.
    // Assert it through `isNone`, NOT through `Number.isFinite(get(...))`: `get()` on
    // a nulled cell returns the raw `FLOAT_NULL`, which IS finite and positive
    // (probe PR-2). An `isFinite` spelling of this assertion reads "true" and tests
    // nothing — the exact trap this whole migration is about, reproduced inside its
    // own regression test.
    expect(agg.col('MEAN_AUC')!.isNone(0), true,
      'one raw NaN destroys the group mean (this is what f64Column prevents)');
    // MEASURED 2026-08-03, not inferred: the raw NaN is ALSO counted as a
    // contributor. So the pre-migration frame was the worst of both worlds — the
    // group mean unusable AND `N` overstating how many subjects supported it.
    expect(agg.col('N_AUC')!.get(0), 3, 'and it IS counted as a contributor — 3, not 2');
  });

  test('PLATFORM RECORD: DG `avg` SKIPS an explicit FLOAT_NULL (na.rm) and does not count it', async () => {
    // The counterpart. Together the pair is the evidence for writing the platform
    // sentinel at construction instead of compensating downstream.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'GRP', ['C', 'C', 'C']),
      DG.Column.fromFloat64Array('AUC', Float64Array.from([100, DG.FLOAT_NULL, 200])),
    ]);
    expect(df.col('AUC')!.isNone(1), true,
      'pre-check: middle position is null after FLOAT_NULL sentinel write');

    const agg = df.groupBy(['GRP']).avg('AUC', 'MEAN_AUC').valueCount('AUC', 'N_AUC').aggregate();
    expectFloat(agg.col('MEAN_AUC')!.get(0) as number, 150, 1e-6,
      'DG avg skips FLOAT_NULL → mean of [100, _, 200] = 150');
    expect(agg.col('N_AUC')!.get(0), 2, 'and it is NOT counted as a contributor — 2, not 3');
  });

  test('units / label carry-through: MEAN+SD inherit source units, CV_PCT="%", CL/F label prefixed', async () => {
    // Mirror writeResultsToDataFrame output for an extravascular run:
    // CMAX has units only; CL has both units and a CL/F label.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'COMPOUND', ['A', 'A', 'B']),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 20, 5])),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 1.5, 0.8])),
    ]);
    df.col('CMAX')!.setTag('units', 'mg/L');
    df.col('CL')!.setTag('units', 'L/h');
    df.col('CL')!.setTag('label', 'CL/F');

    const agg = aggregateNcaResults(df, 'COMPOUND', {params: ['CMAX', 'CL']});

    expect(agg.col('MEAN_CMAX')!.tags['units'], 'mg/L', 'MEAN_CMAX inherits source units');
    expect(agg.col('SD_CMAX')!.tags['units'], 'mg/L', 'SD_CMAX inherits source units');
    expect(agg.col('CV_PCT_CMAX')!.tags['units'], '%', 'CV_PCT_CMAX units = "%"');
    expect(agg.col('MEAN_CL')!.tags['units'], 'L/h', 'MEAN_CL inherits source units');
    expect(agg.col('SD_CL')!.tags['units'], 'L/h', 'SD_CL inherits source units');
    expect(agg.col('CV_PCT_CL')!.tags['units'], '%', 'CV_PCT_CL units = "%"');

    // CL/F-vs-CL distinction preserved via prefixed labels — the grid header
    // reads "Mean CL/F", not an opaque "MEAN_CL" that loses the route convention.
    expect(agg.col('MEAN_CL')!.tags['label'], 'Mean CL/F', 'MEAN_CL label prefixed from source');
    expect(agg.col('SD_CL')!.tags['label'], 'SD CL/F', 'SD_CL label prefixed from source');
    expect(agg.col('CV_PCT_CL')!.tags['label'], '%CV CL/F', 'CV_PCT_CL label prefixed from source');

    // Regression (surfaced by spike/groupby-grid-presentation): a param with NO
    // source label (CMAX) must NOT get a literal "Mean null" label. An unset DG
    // tag reads back as `null`, so the label guard must exclude null, not just
    // undefined — otherwise `Mean ${null}` = "Mean null" leaks into the header.
    expect(agg.col('MEAN_CMAX')!.tags['label'] ?? null, null, 'MEAN_CMAX has no label (not "Mean null")');
    expect(agg.col('SD_CMAX')!.tags['label'] ?? null, null, 'SD_CMAX has no label (not "SD null")');
    expect(agg.col('CV_PCT_CMAX')!.tags['label'] ?? null, null, 'CV_PCT_CMAX has no label (not "%CV null")');
  });
});

// =====================================================================
// group-by-cl-basis-explicit-control F1/F3 — mixed-basis is now WITHHELD, not
// auto-stratified (reversal of per-row-route AC-C4, user-ratified 2026-07-10).
// resolveAggregateGroupKeys honours the user's keys EXACTLY (no CL_BASIS/DOSE_NORM
// append) and returns the SET of mixed cell keys; the caller aggregates by the user's
// keys, then blankMixedBasisFamily withholds (FLOAT_NULL → NE) the corrupted family in
// exactly those pooled cells. These ARE the rule-18 science-preservation tests.
// =====================================================================

/** Build a results frame with a CL_BASIS column + a chosen group key. */
function mixedBasisFrame(): DG.DataFrame {
  // 2 sexes × (IV true-CL + PO apparent-CL/F). Grouping by SEX pools both bases in each
  // Sex cell → the family is withheld (NE) in both output rows.
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', ['M', 'M', 'F', 'F']),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL', 'CL/F', 'CL', 'CL/F']),
    DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 4.0, 1.2, 4.4])),
  ]);
}

/** A wide mixed-basis frame carrying the full IV↔EV-incomparable family + two non-family
 *  params (CMAX, HALF_LIFE) so a blank pass can be asserted param-by-param. */
function wideMixedBasisFrame(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', ['M', 'M', 'F', 'F']),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL', 'CL/F', 'CL', 'CL/F']),
    DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 4.0, 1.2, 4.4])),
    DG.Column.fromFloat64Array('VZ', Float64Array.from([2.0, 8.0, 2.2, 8.2])),
    DG.Column.fromFloat64Array('VSS', Float64Array.from([3.0, 6.0, 3.2, 6.2])),
    DG.Column.fromFloat64Array('MRT', Float64Array.from([4.0, 5.0, 4.2, 5.2])),
    DG.Column.fromFloat64Array('AUMC_INF', Float64Array.from([100, 120, 102, 122])),
    DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 12, 11, 13])),
    DG.Column.fromFloat64Array('HALF_LIFE', Float64Array.from([2.0, 2.5, 2.1, 2.6])),
  ]);
}

category('Adapter: results aggregator — mixed CL_BASIS withhold (F1/F3)', () => {
  test('F1: resolveAggregateGroupKeys does NOT append CL_BASIS; returns the mixed cell keys', async () => {
    const df = mixedBasisFrame();
    const r = resolveAggregateGroupKeys(df, 'SEX');
    expect(r.mixedBasis, true, 'mixed basis still detected');
    expect(r.keys.length, 1, 'keys = the user selection only (no auto-append)');
    expect(r.keys[0], 'SEX', 'the user key is preserved verbatim');
    // The mixed cell keys are the two Sex buckets (each pools CL + CL/F).
    expect(r.mixedBasisKeys.size, 2, 'both Sex cells are mixed-basis');
    expect(r.mixedBasisKeys.has('M'), true, 'M cell is mixed');
    expect(r.mixedBasisKeys.has('F'), true, 'F cell is mixed');
  });

  test('F1: single-basis frame → no detection, empty key set', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', ['M', 'F']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL/F', 'CL/F']),
      DG.Column.fromFloat64Array('CL', Float64Array.from([4.0, 4.4])),
    ]);
    const r = resolveAggregateGroupKeys(df, 'SEX');
    expect(r.mixedBasis, false, 'single basis → no mix');
    expect(r.mixedBasisKeys.size, 0, 'no mixed cells');
    expect(r.keys.length, 1, 'only the selected key');
  });

  test('F1: grouping by CL_BASIS itself self-suppresses (basis constant within each cell)', async () => {
    const r = resolveAggregateGroupKeys(mixedBasisFrame(), 'CL_BASIS');
    expect(r.mixedBasis, false, 'no mix when grouping by CL_BASIS directly');
    expect(r.mixedBasisKeys.size, 0, 'no mixed cells');
    expect(r.keys.length, 1, 'single key');
  });

  test('F3: blank pass withholds the whole family in mixed cells; N_ and non-family survive', async () => {
    const df = wideMixedBasisFrame();
    const {keys, mixedBasisKeys} = resolveAggregateGroupKeys(df, 'SEX');
    const params = ['CL', 'VZ', 'VSS', 'MRT', 'AUMC_INF', 'CMAX', 'HALF_LIFE'];
    const agg = aggregateNcaResults(df, keys, {params});
    expect(agg.rowCount, 2, 'user keys honoured → 2 Sex rows (not 4)');
    blankMixedBasisFamily(agg, mixedBasisKeys, keys, IV_EV_INCOMPARABLE_FAMILY);
    for (let i = 0; i < agg.rowCount; i++) {
      // Family MEAN/SD/CV withheld (NE) in both pooled rows (the family members this frame carries).
      for (const p of IV_EV_INCOMPARABLE_FAMILY.filter((x) => params.includes(x))) {
        for (const stat of [`MEAN_${p}`, `SD_${p}`, `CV_PCT_${p}`])
          expect(agg.col(stat)!.isNone(i), true, `${stat} withheld (NE) in mixed row ${i}`);

        // ...but the contributor count survives the withhold (R2-5).
        expect(agg.col(`N_${p}`)!.get(i), 2, `N_${p} = 2 kept in mixed row ${i}`);
      }
      // Non-family exposure params stay populated (identical units across route).
      expect(agg.col('MEAN_CMAX')!.isNone(i), false, `MEAN_CMAX finite in row ${i}`);
      expect(agg.col('MEAN_HALF_LIFE')!.isNone(i), false, `MEAN_HALF_LIFE finite (A′) in row ${i}`);
      expect(agg.col('N')!.get(i), 2, `group N = 2 in row ${i}`);
    }
  });

  test('F3.2: grouping by Sex × CL basis blanks NOTHING (every cell single-basis)', async () => {
    const df = wideMixedBasisFrame();
    const {keys, mixedBasisKeys} = resolveAggregateGroupKeys(df, ['SEX', 'CL_BASIS']);
    expect(mixedBasisKeys.size, 0, 'no cell mixes basis → nothing to withhold');
    const params = ['CL', 'VZ', 'VSS', 'MRT', 'AUMC_INF'];
    const agg = aggregateNcaResults(df, keys, {params});
    blankMixedBasisFamily(agg, mixedBasisKeys, keys, IV_EV_INCOMPARABLE_FAMILY);
    expect(agg.rowCount, 4, 'four single-basis cells');
    for (let i = 0; i < agg.rowCount; i++) {
      for (const p of params)
        expect(agg.col(`MEAN_${p}`)!.isNone(i), false, `MEAN_${p} finite in single-basis row ${i}`);
    }
  });

  test('F3.4 parity pin: a ROUTE-grouped (1:1 basis) frame is byte-identical — nothing withheld', async () => {
    // Each ROUTE group is single-basis (IV ⇒ CL, PO/SC ⇒ CL/F), so grouping by ROUTE never
    // pools basis → no detection, no blank → aggregate identical to today (PKNCA 60/60 pin).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', ['IV', 'IV', 'PO', 'SC']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL', 'CL', 'CL/F', 'CL/F']),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 1.2, 4.0, 4.4])),
    ]);
    const r = resolveAggregateGroupKeys(df, 'ROUTE');
    expect(r.mixedBasis, false, 'ROUTE 1:1-determines CL_BASIS → no mix');
    expect(r.mixedBasisKeys.size, 0, 'no mixed cells');
    expect(r.keys.length, 1, 'only ROUTE; CL_BASIS not appended');
    const agg = aggregateNcaResults(df, r.keys, {params: ['CL']});
    // Locate rows by the ROUTE group column (this aggregate is keyed on ROUTE, not COMPOUND).
    const routeCol = agg.col('ROUTE')!;
    const rowByRoute = (v: string): number => {
      for (let i = 0; i < agg.rowCount; i++) if (routeCol.get(i) === v) return i;
      throw new Error(`route '${v}' not found in aggregate`);
    };
    const before = agg.col('MEAN_CL')!.get(rowByRoute('IV')) as number;
    blankMixedBasisFamily(agg, r.mixedBasisKeys, r.keys, IV_EV_INCOMPARABLE_FAMILY);
    expect(agg.rowCount, 3, 'three single-basis route groups (IV/PO/SC)');
    // Blank pass is a no-op on the empty key set → the IV mean is untouched.
    expectFloat(agg.col('MEAN_CL')!.get(rowByRoute('IV')) as number, before, 1e-12,
      'MEAN_CL byte-unchanged by the (empty) blank pass');
    expect(agg.col('MEAN_CL')!.isNone(rowByRoute('PO')), false, 'PO mean finite (single basis)');
  });

  test('F2.3 numeric-key coupling (R1-1): blank lands on the mixed cell under a numeric composite key', async () => {
    // Sex × Dose with Dose as a NUMERIC column — the flagship composite. The (M,10) cell
    // pools IV + PO (mixed); (F,20) is single-basis. The blank must land on exactly (M,10),
    // which only holds if compositeGroupKey's getString yields the SAME "10" string on the
    // pre-aggregation frame (set construction) and the post-aggregation frame (per-row match)
    // — the type-preservation the architecture spine relies on.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', ['M', 'M', 'F', 'F']),
      DG.Column.fromFloat64Array('DOSE', Float64Array.from([10, 10, 20, 20])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL', 'CL/F', 'CL', 'CL']),
      DG.Column.fromFloat64Array('CL', Float64Array.from([1.0, 4.0, 1.2, 1.3])),
    ]);
    const r = resolveAggregateGroupKeys(df, ['SEX', 'DOSE']);
    expect(r.mixedBasis, true, '(M,10) mixes basis');
    expect(r.mixedBasisKeys.size, 1, 'exactly one mixed cell');
    const agg = aggregateNcaResults(df, r.keys, {params: ['CL']});
    expect(agg.rowCount, 2, 'two composite cells: (M,10) and (F,20)');
    blankMixedBasisFamily(agg, r.mixedBasisKeys, r.keys, IV_EV_INCOMPARABLE_FAMILY);
    // Locate rows by their (SEX,DOSE) identity on the aggregate frame.
    const sexCol = agg.col('SEX')!; const doseCol = agg.col('DOSE')!;
    let mRow = -1; let fRow = -1;
    for (let i = 0; i < agg.rowCount; i++) {
      if (sexCol.get(i) === 'M' && Number(doseCol.get(i)) === 10) mRow = i;
      if (sexCol.get(i) === 'F' && Number(doseCol.get(i)) === 20) fRow = i;
    }
    expect(mRow >= 0 && fRow >= 0, true, 'both composite cells present');
    expect(agg.col('MEAN_CL')!.isNone(mRow), true, 'mixed (M,10) cell → MEAN_CL withheld (NE)');
    expect(agg.col('MEAN_CL')!.isNone(fRow), false, 'single-basis (F,20) cell → MEAN_CL finite');
    expectFloat(agg.col('MEAN_CL')!.get(fRow) as number, 1.25, 1e-9, '(F,20) mean = mean of [1.2,1.3]');
  });
});

// AD-4 (nca-cl-bw-normalised-units) — CL / VZ / VSS are per kg on every row, so a group mixing
// `absolute` and `per_kg` ENGINE doses pools one basis; the former DOSE_NORM withholding is gone.
// The fixture is a mixed-DOSE_NORM frame as the writer now emits it — per-kg CL on both kinds of
// row, the uniform mL/min/kg tag, the absolute companion null on the per-kg rows.
function mixedDoseNormFrame(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', ['M', 'M', 'F', 'F']),
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'DOSE_NORM',
      ['absolute', 'per_kg', 'absolute', 'per_kg']),
    DG.Column.fromFloat64Array('CL', Float64Array.from([5.0, 6.0, 5.2, 6.4])),
    f64Column('CL_ABS', [1.25, NaN, 1.3, NaN]),
  ]);
  df.col('CL')!.setTag('units', 'mL/min/kg');
  df.col('CL_ABS')!.setTag('units', 'L/h');
  return df;
}

category('Adapter: results aggregator — per-kg identity across DOSE_NORM (AD-4)', () => {
  test('AC-3.1: a group mixing absolute and per_kg engine doses yields a per-kg CL mean over both', async () => {
    const df = mixedDoseNormFrame();
    const r = resolveAggregateGroupKeys(df, 'SEX');
    expect(r.mixedBasis, false, 'no CL_BASIS column → nothing withheld');
    const agg = aggregateNcaResults(df, r.keys, {params: ['CL', 'CL_ABS']});
    blankMixedBasisFamily(agg, r.mixedBasisKeys, r.keys, IV_EV_INCOMPARABLE_FAMILY);
    const m = rowByGroup(agg, 'M', 'SEX');
    expectFloat(agg.col('MEAN_CL')!.get(m) as number, 5.5, 1e-12, 'M per-kg mean = (5.0 + 6.0) / 2');
    expect(agg.col('N_CL')!.get(m), 2, 'both rows contribute to the per-kg mean');
    // The absolute companion averages only the rows that have it, with its own N.
    expectFloat(agg.col('MEAN_CL_ABS')!.get(m) as number, 1.25, 1e-12, 'CL_ABS mean over the absolute row');
    expect(agg.col('N_CL_ABS')!.get(m), 1, 'N_CL_ABS = 1 (null on the per-kg row)');
  });

  test('AC-3.2: the CL-basis guard covers the companions as it covers CL / VZ / VSS', async () => {
    for (const p of ['CL_ABS', 'VZ_ABS', 'VSS_ABS'])
      expect(IV_EV_INCOMPARABLE_FAMILY.includes(p), true, `${p} withheld with its per-kg sibling`);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SEX', ['M', 'M']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CL_BASIS', ['CL', 'CL/F']),
      DG.Column.fromFloat64Array('CL_ABS', Float64Array.from([1.0, 4.0])),
    ]);
    const r = resolveAggregateGroupKeys(df, 'SEX');
    const agg = aggregateNcaResults(df, r.keys, {params: ['CL_ABS']});
    blankMixedBasisFamily(agg, r.mixedBasisKeys, r.keys, IV_EV_INCOMPARABLE_FAMILY);
    expect(agg.col('MEAN_CL_ABS')!.isNone(0), true, 'CL + CL/F pooled → MEAN_CL_ABS withheld (NE)');
  });

  test('AC-3.3: MEAN_ / SD_ carry the uniform source units', async () => {
    const df = mixedDoseNormFrame();
    const agg = aggregateNcaResults(df, 'SEX', {params: ['CL', 'CL_ABS']});
    expect(agg.col('MEAN_CL')!.tags['units'], 'mL/min/kg', 'MEAN_CL carries mL/min/kg');
    expect(agg.col('SD_CL')!.tags['units'], 'mL/min/kg', 'SD_CL carries mL/min/kg');
    expect(agg.col('MEAN_CL_ABS')!.tags['units'], 'L/h', 'MEAN_CL_ABS carries L/h');
    expect(agg.col('CV_PCT_CL')!.tags['units'], '%', 'CV is always dimensionless %');
  });
});

category('Adapter: results aggregator — family membership (F3, rule-18 audit)', () => {
  test('IV_EV_INCOMPARABLE_FAMILY = CL, VZ, VSS, MRT, AUMC_INF + the absolute companions', async () => {
    expect(IV_EV_INCOMPARABLE_FAMILY.length, 8, 'five per-kg/moment params + three companions');
    for (const p of ['CL', 'VZ', 'VSS', 'MRT', 'AUMC_INF', 'CL_ABS', 'VZ_ABS', 'VSS_ABS'])
      expect(IV_EV_INCOMPARABLE_FAMILY.includes(p), true, `${p} in the CL-basis family`);
    // Exposure params must NOT be withheld (identical units across route; A′ populate).
    for (const p of ['CMAX', 'AUC_LAST', 'AUC_INF', 'TMAX', 'HALF_LIFE', 'LAMBDA_Z', 'PCT_EXTRAP'])
      expect(IV_EV_INCOMPARABLE_FAMILY.includes(p), false, `${p} NOT withheld`);
  });
});
