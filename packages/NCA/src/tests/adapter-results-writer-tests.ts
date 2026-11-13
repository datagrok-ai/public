/**
 * Adapter — results-writer category.
 *
 * @see src/adapter/results-writer.ts (implementation under test)
 *
 * Covers: 36-column schema (+ CL_BASIS + DOSE_NORM + 6 FR-200 moment/lag cols + the per-kg
 * companions), Float64 invariant (`fromFloat64Array`
 * storage verified via `getRawData() instanceof Float64Array`), NaN sentinel for
 * unfit parameters, per-row CL/Vz basis (neutral `'CL'`/`'Vz'` label + per-row
 * `CL_BASIS`), `LAMBDAZ_NPOINTS` `INT_NULL` for null λz (NOTE-02), and the `.d42`
 * NaN round-trip empirical check that satisfies the Phase 3 deferral.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {writeResultsToDataFrame, patchResultRow} from '../adapter/results-writer';
import {scanNumericSentinels} from '../adapter/adapter-utils';
import type {ProfileIndex} from '../adapter/profile-index';

const RESULT_COLUMN_NAMES: ReadonlyArray<string> = [
  'SUBJECT', 'ROUTE', 'CL_BASIS', 'DOSE_NORM', 'STATUS',
  'CMAX', 'C0_EXTRAP', 'TMAX', 'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP', 'LAMBDA_Z', 'HALF_LIFE',
  'CL', 'VZ',
  'AUMC_LAST', 'AUMC_INF', 'MRT', 'VSS', 'TLAG', 'PCT_AUMC_EXTRAP',
  'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG',
  'LAMBDAZ_R2', 'LAMBDAZ_ADJ_R2', 'LAMBDAZ_NPOINTS', 'SPAN_RATIO', 'LAMBDAZ_TFIRST',
  'LAMBDAZ_TLAST', 'LAMBDAZ_INTERCEPT',
  'AUC_METHOD', 'COMPENSATED',
  'WARNINGS_JSON', 'PROVENANCE_JSON',
];

function key(subject: string): ProfileIndex {
  return {subjectKey: subject, observationRows: Int32Array.from([0])};
}

function okValues(): nca.ParameterValues {
  return {
    cmax: 10, tmax: 1, aucLast: 50, aucInf: 60, pctExtrap: 16.67,
    lambdaZ: 0.1, halfLife: 6.93, cl: 1.67, vz: 16.7,
    aumcLast: 300, aumcInf: 420, mrt: 7, vss: 11.67, tlag: 0,
    pctExtrapAumc: 28.57,
  };
}

function okLambdaZ(): nca.LambdaZResult {
  return {
    lambdaZ: 0.1,
    intercept: 2.3,
    rSquared: 0.999,
    adjRSquared: 0.998,
    pointsUsed: Int32Array.from([2, 3, 4]),
    tStart: 2,
    tEnd: 8,
    // sci-comp 0.10.0. Consistent with the rest of this double rather than
    // arbitrary: t½ = ln2/0.1 = 6.93, window 2→8 = 6 h, so span ≈ 0.866.
    spanRatio: (8 - 2) / (Math.LN2 / 0.1),
  };
}

function okProvenance(): nca.ProfileProvenance {
  return {
    lambdaZ: okLambdaZ(),
    blqApplied: {conc: new Float64Array(0), excluded: new Int32Array(0)},
    aucMethod: 'linear-up-log-down',
    compensated: false,
    warnings: [],
  };
}

function okResult(): nca.ComputeResult {
  return {values: okValues(), provenance: okProvenance(), status: 'ok'};
}

function failedResult(): nca.ComputeResult {
  return {
    values: {
      cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
      lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN,
      aumcLast: NaN, aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN,
      pctExtrapAumc: NaN,
    },
    provenance: {
      lambdaZ: null,
      blqApplied: {conc: new Float64Array(0), excluded: new Int32Array(0)},
      aucMethod: 'linear-up-log-down',
      compensated: false,
      warnings: [{code: 'ENGINE_ERROR', severity: 'error', message: 'no measurable points'}],
    },
    status: 'failed',
  };
}

category('Adapter: results writer', () => {
  test('AC-2.1: writes the 36-column UPPERCASE schema in order', async () => {
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()],
    );
    expect(df.columns.length, 36,
      'expected 36 columns (22 + C0_EXTRAP + CL_BASIS + DOSE_NORM + 6 FR-200 ' +
      'moment/lag + SPAN_RATIO + CL_ABS/VZ_ABS/VSS_ABS/BW_KG)');
    const names = df.columns.names();
    expectArray(names, RESULT_COLUMN_NAMES as unknown as ArrayLike<unknown>);
  });

  test('numeric columns are Float64Array-backed (fromFloat64Array invariant)', async () => {
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()],
    );
    const numericNames = [
      'CMAX', 'C0_EXTRAP', 'TMAX', 'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP', 'LAMBDA_Z',
      'HALF_LIFE', 'CL', 'VZ', 'LAMBDAZ_R2', 'LAMBDAZ_ADJ_R2', 'SPAN_RATIO',
      'LAMBDAZ_TFIRST', 'LAMBDAZ_TLAST', 'LAMBDAZ_INTERCEPT',
      'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG',
    ];
    for (const name of numericNames) {
      const raw = df.col(name)!.getRawData();
      expect(raw instanceof Float64Array, true,
        `${name} must be Float64Array-backed; got ${raw.constructor.name}`);
    }
    const npoints = df.col('LAMBDAZ_NPOINTS')!.getRawData();
    expect(npoints instanceof Int32Array, true,
      `LAMBDAZ_NPOINTS must be Int32Array-backed; got ${npoints.constructor.name}`);
  });

  test('CL/VZ labels are neutral; CL_BASIS "CL" for all intravenous subjects', async () => {
    const df = writeResultsToDataFrame(
      [key('S1'), key('S2')],
      [nca.ROUTE_IV_BOLUS, nca.ROUTE_IV_INFUSION],
      [okResult(), okResult()],
    );
    // AD-4: neutral parameter-name labels — basis is now per-row.
    expect(df.col('CL')!.tags['label'], 'CL');
    expect(df.col('VZ')!.tags['label'], 'Vz');
    expect(df.col('CL_BASIS')!.get(0), 'CL', 'IV-bolus → CL');
    expect(df.col('CL_BASIS')!.get(1), 'CL', 'IV-infusion → CL');
  });

  test('mixed route: per-row CL_BASIS holds "CL" and "CL/F" simultaneously (AC-C1)', async () => {
    const df = writeResultsToDataFrame(
      [key('S1'), key('S2'), key('S3')],
      [nca.ROUTE_IV_BOLUS, nca.ROUTE_PO, nca.ROUTE_SC],
      [okResult(), okResult(), okResult()],
    );
    // Neutral label even on a mixed frame (no frame-level CL/F lie).
    expect(df.col('CL')!.tags['label'], 'CL');
    expect(df.col('VZ')!.tags['label'], 'Vz');
    // Per-row basis follows each subject's route, all in one frame.
    expect(df.col('CL_BASIS')!.get(0), 'CL', 'IV-bolus row → CL');
    expect(df.col('CL_BASIS')!.get(1), 'CL/F', 'PO row → CL/F');
    expect(df.col('CL_BASIS')!.get(2), 'CL/F', 'SC row → CL/F');
  });

  // ----- nca-cl-bw-normalised-units AD-1/AD-3: the column identity (CL / VZ / VSS per kg on every
  // row, absolute companions + BW_KG), uniform unit tags from the one map (AC-2.2). DOSE_NORM keeps
  // its engine meaning (nca-cl-vz-convention AD-3). -----

  test('AC-2.2: every frame carries the uniform tags — mixed DOSE_NORM included', async () => {
    for (const norms of [undefined, ['per_kg', 'per_kg'], ['absolute', 'per_kg']] as const) {
      const df = writeResultsToDataFrame(
        [key('S1'), key('S2')], [nca.ROUTE_IV_BOLUS, nca.ROUTE_IV_BOLUS], [okResult(), okResult()],
        undefined, norms === undefined ? undefined : [...norms]);
      const want: Record<string, string> = {
        CL: 'mL/min/kg', VZ: 'L/kg', VSS: 'L/kg', CL_ABS: 'L/h', VZ_ABS: 'L', VSS_ABS: 'L', BW_KG: 'kg'};
      for (const [c, u] of Object.entries(want))
        expect(df.col(c)!.tags['units'], u, `${c} tag on DOSE_NORM ${JSON.stringify(norms)}`);
    }
  });

  test('AC-2.2: per_kg row → per-kg values in mL/min/kg; absolute companions blank', async () => {
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()], undefined, ['per_kg']);
    expect(df.col('DOSE_NORM')!.get(0), 'per_kg', 'DOSE_NORM keeps the engine basis');
    expect(Math.abs((df.col('CL')!.get(0) as number) - 1.67 * 1000 / 60) < 1e-9, true, 'CL per kg × 1000/60');
    expect(df.col('VZ')!.get(0), 16.7, 'Vz per kg unscaled');
    expect(df.col('CL_ABS')!.isNone(0), true, 'no weight to form CL abs.');
    expect(df.col('BW_KG')!.isNone(0), true, 'no weight');
  });

  test('AC-2.2: absolute row with a body weight → per kg = abs / BW; companions = engine', async () => {
    const bw = {bodyWeightKg: 0.25, bodyWeightRangeKg: {min: 0.25, max: 0.25}, perKgDose: true};
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()], undefined, ['absolute'], [bw]);
    expect(Math.abs((df.col('CL')!.get(0) as number) - 1.67 / 0.25 * 1000 / 60) < 1e-9, true, 'CL per kg');
    expect(Math.abs((df.col('VZ')!.get(0) as number) - 16.7 / 0.25) < 1e-9, true, 'Vz per kg');
    expect(df.col('CL_ABS')!.get(0), 1.67, 'CL abs. = engine');
    expect(df.col('BW_KG')!.get(0), 0.25, 'BW_KG = the weight used');
  });

  test('AC-2.2: absolute row without a body weight → per kg blank, companions kept', async () => {
    const df = writeResultsToDataFrame([key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()]);
    expect(df.col('DOSE_NORM')!.get(0), 'absolute', 'default basis');
    expect(df.col('CL')!.isNone(0), true, 'no per kg without a weight');
    expect(df.col('CL_ABS')!.get(0), 1.67, 'CL abs. reported');
  });

  test('bodyWeights length mismatch throws', async () => {
    let threw = false;
    try {
      writeResultsToDataFrame([key('S1'), key('S2')], [nca.ROUTE_PO, nca.ROUTE_PO], [okResult(), okResult()],
        undefined, undefined, [{bodyWeightKg: null, bodyWeightRangeKg: null, perKgDose: false}]);
    } catch (e) {
      threw = e instanceof Error && e.message.includes('bodyWeights length mismatch');
    }
    expect(threw, true, 'mismatched bodyWeights must throw');
  });

  test('AC-4: doseNorms length mismatch throws', async () => {
    let threw = false;
    try {
      writeResultsToDataFrame(
        [key('S1'), key('S2')],
        [nca.ROUTE_IV_BOLUS, nca.ROUTE_IV_BOLUS],
        [okResult(), okResult()],
        undefined,
        ['per_kg'],
      );
    } catch (e) {
      threw = e instanceof Error && e.message.includes('doseNorms length mismatch');
    }
    expect(threw, true, 'mismatched doseNorms must throw');
  });

  test('failed result: the PLATFORM null is written, not a raw NaN', async () => {
    // **This assertion is INVERTED from what it used to be**
    // (nca-dg-native-missing-values, user decision D1, 2026-08-03). It previously
    // demanded `Number.isNaN(cmaxRaw[0])` and explicitly forbade `FLOAT_NULL`.
    // That convention was wrong at the platform level: `isNone` returns FALSE for
    // a raw NaN, so every DG-native consumer — grid sort, column filter,
    // `col.stats`, group-by, CSV/.d42 export — read an unfit parameter as a
    // PRESENT value, and a single such cell poisoned the whole column's
    // statistics. The old test passed while that defect shipped.
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [failedResult()],
    );
    const cmax = df.col('CMAX')!;
    expect(cmax.isNone(0), true, 'CMAX[0] must read as MISSING to the platform');
    const cmaxRaw = cmax.getRawData() as Float64Array;
    expect(cmaxRaw[0] === DG.FLOAT_NULL, true,
      `CMAX[0] must be DG.FLOAT_NULL in the Float64 buffer; got ${cmaxRaw[0]}`);
    expect(Number.isNaN(cmaxRaw[0]), false,
      'CMAX[0] must NOT be a raw IEEE NaN — DG does not treat that as missing');
    // And the payoff, asserted rather than assumed: the column's statistics are
    // computable. With a raw NaN they were all NaN (probe PR-3).
    expect(cmax.stats.missingValueCount, 1, 'the unfit cell is counted as missing');
    expect(cmax.stats.valueCount, 0, 'and not as a value');
  });

  test('failed result: LAMBDAZ_NPOINTS uses INT_NULL (not NaN) when λz=null (NOTE-02)', async () => {
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [failedResult()],
    );
    const np = df.col('LAMBDAZ_NPOINTS')!.getRawData() as Int32Array;
    expect(np[0], DG.INT_NULL, 'LAMBDAZ_NPOINTS unfit sentinel is INT_NULL');
  });

  test('SPAN_RATIO: carried from lz.spanRatio when fitted, isNone when λz=null', async () => {
    const ok = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()],
    );
    const okCol = ok.col('SPAN_RATIO')!;
    expect(okCol.isNone(0), false, 'a fitted profile reports a span ratio');
    // Read straight off the double, so this pins the CARRY (writer copies
    // lz.spanRatio verbatim, rule 5 — never recomputed here from tStart/tEnd/λz)
    // rather than re-deriving the formula and testing arithmetic against itself.
    expect(okCol.get(0), okLambdaZ().spanRatio, 'SPAN_RATIO is lz.spanRatio verbatim');
    expect(okCol.getRawData() instanceof Float64Array, true,
      'SPAN_RATIO is Float64Array-backed (NFR-04 Float64 invariant)');

    const failed = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [failedResult()],
    );
    const failedCol = failed.col('SPAN_RATIO')!;
    // The trap this asserts: `f64Column` must have translated the no-fit NaN to
    // DG.FLOAT_NULL at construction. A raw NaN would leave `isNone` FALSE, so the
    // grid would read a failed fit as a present value, render it as the literal
    // text "NaN", and collapse min/max/avg on the whole column.
    expect(failedCol.isNone(0), true,
      'an unfit profile reports SPAN_RATIO as platform-null, not a raw NaN');
    expect(failedCol.stats.valueCount, 0, 'and not as a value');
  });

  test('PROVENANCE_JSON carries pointsUsed and excluded as arrays', async () => {
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [okResult()],
    );
    const json = df.col('PROVENANCE_JSON')!.get(0) as string;
    const parsed = JSON.parse(json);
    expectArray(parsed.lambdaZ.pointsUsed, [2, 3, 4]);
    expectArray(parsed.blqApplied.excluded, []);
  });

  test('parallel-array length mismatch throws', async () => {
    let threw = false;
    try {
      writeResultsToDataFrame(
        [key('S1'), key('S2')], [nca.ROUTE_IV_BOLUS], [okResult()],
      );
    } catch (e) {
      threw = e instanceof Error && e.message.includes('parallel-array length mismatch');
    }
    expect(threw, true, 'mismatched parallel arrays must throw');
  });

  test('.d42 round-trip: the null sentinel survives toByteArray/fromByteArray bit-exact', async () => {
    // The question this test was originally written to answer ("does the unfit
    // sentinel survive serialization, or do we need parallel *_VALID bool-mask
    // columns?") is CLOSED and the answer is favourable: `DG.FLOAT_NULL` survives
    // bit-exact AND still reads `isNone`, so the migration required no schema
    // change (probe PR-6, 2026-07-30). Re-pointed from the NaN convention.
    const df = writeResultsToDataFrame(
      [key('S1')], [nca.ROUTE_IV_BOLUS], [failedResult()],
    );
    const bytes = df.toByteArray();
    const round = DG.DataFrame.fromByteArray(bytes);
    const cmax = round.col('CMAX')!;
    expect(cmax.isNone(0), true, 'CMAX[0] must still read as MISSING after the round-trip');
    const raw = cmax.getRawData() as Float64Array;
    expect(raw[0] === DG.FLOAT_NULL, true,
      `CMAX[0] must remain DG.FLOAT_NULL after .d42 round-trip; got ${raw[0]}`);
    // Also assert LAMBDAZ_NPOINTS INT_NULL survives.
    const np = round.col('LAMBDAZ_NPOINTS')!.getRawData() as Int32Array;
    expect(np[0], DG.INT_NULL,
      'LAMBDAZ_NPOINTS must remain INT_NULL after .d42 round-trip');
  });
});

category('Adapter: patchResultRow (in-place single-row update)', () => {
  test('patches only the target row; frame identity + other rows unchanged; FLOAT_NULL per-cell', async () => {
    const df = writeResultsToDataFrame(
      [key('A'), key('B'), key('C')],
      [nca.ROUTE_IV_BOLUS, nca.ROUTE_PO, nca.ROUTE_IV_BOLUS],
      [okResult(), okResult(), okResult()],
    );
    const before = df.rowCount;
    const cmaxA = df.col('CMAX')!.get(0);
    const cmaxC = df.col('CMAX')!.get(2);

    patchResultRow(df, key('B'), nca.ROUTE_PO, failedResult(), NaN, 'absolute');

    expect(df.rowCount, before, 'row count unchanged (in-place, no rebuild)');
    expect(df.col('STATUS')!.getString(1), 'failed', 'B status patched');
    expect(df.col('CMAX')!.isNone(1), true,
      'B CMAX now MISSING (per-cell NaN→FLOAT_NULL, not a raw NaN) — N-2');
    // No raw NaN leaked anywhere by the patch (the exact defect N-2 guards).
    const rawNaN = scanNumericSentinels(df).filter((e) => e.kind === 'raw-nan');
    expect(rawNaN.length, 0, 'patch must not write a raw NaN into any numeric column');
    // A and C untouched.
    expect(df.col('CMAX')!.get(0), cmaxA, 'A CMAX unchanged');
    expect(df.col('CMAX')!.get(2), cmaxC, 'C CMAX unchanged');
    expect(df.col('STATUS')!.getString(0), 'ok', 'A status unchanged');
    expect(df.col('STATUS')!.getString(2), 'ok', 'C status unchanged');
  });

  test('AC-2.3: patched row (non-null values) equals a full rebuild of that subject', async () => {
    // A body weight on B so the per-kg + companion columns carry values, not just nulls.
    const bw = {bodyWeightKg: 0.3, bodyWeightRangeKg: {min: 0.3, max: 0.3}, perKgDose: true};
    const none = {bodyWeightKg: null, bodyWeightRangeKg: null, perKgDose: false};
    const patched = writeResultsToDataFrame(
      [key('A'), key('B')], [nca.ROUTE_IV_BOLUS, nca.ROUTE_PO], [okResult(), failedResult()]);
    patchResultRow(patched, key('B'), nca.ROUTE_PO, okResult(), 5.0, 'absolute', bw);

    const rebuilt = writeResultsToDataFrame(
      [key('A'), key('B')], [nca.ROUTE_IV_BOLUS, nca.ROUTE_PO], [okResult(), okResult()],
      Float64Array.from([NaN, 5.0]), ['absolute', 'absolute'], [none, bw]);
    expect(rebuilt.col('CL')!.isNone(1), false, 'the rebuild formed a per-kg CL for B');

    for (const name of RESULT_COLUMN_NAMES) {
      const pc = patched.col(name)!;
      const rc = rebuilt.col(name)!;
      if (pc.isNone(1) || rc.isNone(1))
        expect(pc.isNone(1), rc.isNone(1), `${name}[B]: missing-ness must match rebuild`);
      else
        expect(pc.getString(1), rc.getString(1), `${name}[B]: value must match rebuild`);
    }
  });

  test('unknown subject throws', async () => {
    const df = writeResultsToDataFrame([key('A')], [nca.ROUTE_PO], [okResult()]);
    let threw = false;
    try {
      patchResultRow(df, key('ZZ'), nca.ROUTE_PO, okResult(), NaN, 'absolute');
    } catch (e) {
      threw = e instanceof Error && e.message.includes('not found');
    }
    expect(threw, true, 'patching a subject not in the frame must throw');
  });
});
