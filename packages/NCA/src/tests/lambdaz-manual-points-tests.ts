/**
 * Adapter — λz manual point-selection index translation (P-B1c) category.
 *
 * @see src/adapter/lambdaz-manual-points.ts (implementation under test)
 *
 * Pure unit tests — the module is DG-free, so ProfileInputs literals suffice.
 * Boundary cases: extravascular / IV-bolus / native-t=0 augmentation, embedded
 * BLQ (order preserved), DUPLICATE-TIME rows (distinct row ids resolve), and the
 * force-in/force-out/auto membership algebra.
 */

import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {DEFAULT_BLQ_STRATEGY} from '../adapter/defaults';
import {
  buildAugmentedProfile, buildManualPoints, pointsUsedToRows, subjectHasManualOverride,
  applyMinPointsFloor, LAMBDAZ_BELOW_MIN_POINTS, isLambdaZNonContiguous,
  formatOverrideTimestamp,
} from '../adapter/lambdaz-manual-points';
import type {AugmentedProfile} from '../adapter/lambdaz-manual-points';
import type {LzAction} from '../adapter/lambdaz-manual-points';

function resultWithLz(nPoints: number): nca.ComputeResult {
  return {
    status: 'ok',
    values: {
      cmax: 10, tmax: 1, aucLast: 50, aucInf: 60, pctExtrap: 16, lambdaZ: 0.1, halfLife: 6.93,
      cl: 1.6, vz: 16, aumcLast: 300, aumcInf: 420, mrt: 7, vss: 11, tlag: 0, pctExtrapAumc: 28,
    },
    provenance: {
      lambdaZ: {
        lambdaZ: 0.1, intercept: 2.3, rSquared: 0.99, adjRSquared: 0.98,
        pointsUsed: Int32Array.from(Array.from({length: nPoints}, (_, i) => i + 2)),
        tStart: 2, tEnd: 8, spanRatio: 0.86,
      },
      blqApplied: {conc: new Float64Array(0), excluded: new Int32Array(0)},
      aucMethod: 'linear-up-log-down', compensated: false, warnings: [],
    },
  };
}

function inputs(
  time: number[], conc: number[], route: nca.RouteCode, blq?: number[],
): nca.ProfileInputs {
  return {
    time: Float64Array.from(time),
    conc: Float64Array.from(conc),
    blqMask: Uint8Array.from(blq ?? time.map(() => 0)),
    lloq: 0,
    dose: 100,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    route,
    infusionDuration: null,
    bodyWeight: null,
  };
}

function rows(...r: number[]): Int32Array {
  return Int32Array.from(r);
}

category('Adapter: λz manual points (P-B1c)', () => {
  test('extravascular: (0,0) prepend, augRowId[0] = -1 (AC-2.1 b)', async () => {
    const aug = buildAugmentedProfile(
      inputs([1, 2, 4, 8], [10, 20, 15, 5], nca.ROUTE_PO), rows(3, 7, 11, 15), DEFAULT_BLQ_STRATEGY);
    expect(aug.rowId.length, 5, 'prepend adds one element');
    expectArray(Array.from(aug.rowId), [-1, 3, 7, 11, 15]);
    expectArray(Array.from(aug.time), [0, 1, 2, 4, 8]);
    expectArray(Array.from(aug.conc), [0, 10, 20, 15, 5]);
  });

  test('native measurable t=0: no prepend (AC-2.1 c)', async () => {
    const aug = buildAugmentedProfile(
      inputs([0, 1, 2, 4], [100, 50, 25, 12], nca.ROUTE_PO), rows(0, 1, 2, 3), DEFAULT_BLQ_STRATEGY);
    expect(aug.rowId.length, 4, 'no augmentation when a measurable t=0 exists');
    expectArray(Array.from(aug.rowId), [0, 1, 2, 3]);
  });

  test('IV-bolus: C0 prepend, augRowId[0] = -1 (AC-2.1 a)', async () => {
    const aug = buildAugmentedProfile(
      inputs([0.5, 1, 2, 4], [100, 60, 30, 10], nca.ROUTE_IV_BOLUS), rows(5, 6, 7, 8),
      DEFAULT_BLQ_STRATEGY);
    expect(aug.rowId.length, 5, 'IV-bolus with no t=0 prepends back-extrapolated C0');
    expect(aug.rowId[0], -1, 'synthetic prepend has no source row');
    expectArray(Array.from(aug.rowId.slice(1)), [5, 6, 7, 8]);
  });

  test('IV-bolus (0, BLQ) pre-dose row is REPLACED: rowId[0] = -1, source row absent (AC-N.2.1)', async () => {
    // Source rows 40..44: row 40 is a flagged pre-dose sample at t=0. The engine's
    // route-aware gate removes it and puts (0, c0) in its slot, so the augmented
    // profile has the INPUT's length and row 40 appears nowhere in rowId — a
    // positional mirror would have mapped augmented element 0 onto row 40.
    const aug = buildAugmentedProfile(
      inputs([0, 0.5, 1, 2, 4], [0, 100, 60, 30, 10], nca.ROUTE_IV_BOLUS, [1, 0, 0, 0, 0]),
      rows(40, 41, 42, 43, 44), DEFAULT_BLQ_STRATEGY);
    expect(aug.rowId.length, 5, 'REPLACE keeps the input length');
    expect(aug.rowId[0], -1, 'slot 0 is the synthetic (0, c0)');
    expectArray(Array.from(aug.rowId), [-1, 41, 42, 43, 44]);
    expect(Array.from(aug.rowId).includes(40), false, 'the pre-dose source row is absent from the fit space');
    expect(aug.time[0], 0, 'dose time');
    expect(aug.conc[0] > 100, true, 'slot 0 carries the back-extrapolated c0, not the BLQ value');
    // A force-in on the replaced row is ignored (it is not selectable)…
    const mp = buildManualPoints(aug.rowId, new Set<number>([43, 44]), new Map<number, LzAction>([[40, 'force-in']]));
    expectArray(Array.from(mp), [3, 4]);
    // …and pointsUsedToRows never returns it, even if the engine used slot 0.
    expectArray(pointsUsedToRows(Int32Array.from([0, 3, 4]), aug.rowId), [43, 44]);
    // The unflagged-zero pre-dose row takes the same path.
    const zero = buildAugmentedProfile(
      inputs([0, 0.5, 1, 2, 4], [0, 100, 60, 30, 10], nca.ROUTE_IV_BOLUS), rows(40, 41, 42, 43, 44),
      DEFAULT_BLQ_STRATEGY);
    expectArray(Array.from(zero.rowId), [-1, 41, 42, 43, 44]);
  });

  test('no measurable point (all BLQ): empty augmented profile, nothing to translate', async () => {
    const aug = buildAugmentedProfile(
      inputs([0, 1, 2], [0, 0, 0], nca.ROUTE_PO, [1, 1, 1]), rows(1, 2, 3), DEFAULT_BLQ_STRATEGY);
    expect(aug.rowId.length, 0, 'kernel null → empty');
    expectArray(pointsUsedToRows(Int32Array.from([0, 1]), aug.rowId), []);
  });

  test('duplicate-time rows get DISTINCT ids — the F-2 fix (AC-2.1 f)', async () => {
    // two source rows (20, 21) share time = 1
    const aug = buildAugmentedProfile(
      inputs([1, 1, 2, 4], [10, 12, 8, 4], nca.ROUTE_PO), rows(20, 21, 22, 23), DEFAULT_BLQ_STRATEGY);
    expectArray(Array.from(aug.rowId), [-1, 20, 21, 22, 23]);
    // force-OUT only row 20 (not its same-time sibling 21)
    const actions = new Map<number, LzAction>([[20, 'force-out']]);
    const autoIn = new Set<number>([20, 21, 22, 23]);
    const mp = buildManualPoints(aug.rowId, autoIn, actions);
    // augIdx 1 (row 20) dropped; 2,3,4 (rows 21,22,23) kept — a time-value key
    // could not have separated 20 from 21.
    expectArray(Array.from(mp), [2, 3, 4]);
  });

  test('embedded BLQ preserves order + positional ids (AC-2.1 e)', async () => {
    const aug = buildAugmentedProfile(
      inputs([1, 2, 3, 4], [10, 0, 8, 4], nca.ROUTE_PO, [0, 1, 0, 0]), rows(30, 31, 32, 33),
      DEFAULT_BLQ_STRATEGY);
    expectArray(Array.from(aug.rowId), [-1, 30, 31, 32, 33]);
  });

  test('pointsUsedToRows drops the synthetic point + maps the rest (AC-2.2)', async () => {
    const augRowId = Int32Array.from([-1, 5, 6, 7, 8]);
    // engine used augmented indices 0 (synthetic), 3, 4
    const usedRows = pointsUsedToRows(Int32Array.from([0, 3, 4]), augRowId);
    expectArray(usedRows, [7, 8]);
  });

  test('buildManualPoints: force-in adds, force-out removes, auto follows window (AC-2.1)', async () => {
    const augRowId = Int32Array.from([-1, 10, 11, 12, 13]);
    const autoIn = new Set<number>([12, 13]); // auto window = rows 12,13
    const actions = new Map<number, LzAction>([[11, 'force-in'], [12, 'force-out']]);
    // 11 forced in; 12 forced out; 13 auto-in; 10 auto-out
    const mp = buildManualPoints(augRowId, autoIn, actions);
    expectArray(Array.from(mp), [2, 4]); // augIdx of rows 11 and 13
  });

  test('all-auto subject: buildManualPoints returns exactly the auto window', async () => {
    const augRowId = Int32Array.from([-1, 40, 41, 42]);
    const autoIn = new Set<number>([41, 42]);
    const mp = buildManualPoints(augRowId, autoIn, new Map());
    expectArray(Array.from(mp), [2, 3]);
  });

  test('subjectHasManualOverride: only explicit force actions count', async () => {
    expect(subjectHasManualOverride(new Map()), false, 'empty → auto');
    expect(subjectHasManualOverride(new Map<number, LzAction>([[1, 'auto']])), false, 'auto → false');
    expect(subjectHasManualOverride(new Map<number, LzAction>([[1, 'force-out']])), true, 'force-out → true');
    expect(subjectHasManualOverride(new Map<number, LzAction>([[1, 'force-in']])), true, 'force-in → true');
  });

  test('applyMinPointsFloor: 2-point fit < minPoints(3) → suppressed to Not estimated (AC-4.3)', async () => {
    const out = applyMinPointsFloor(resultWithLz(2), 3);
    expect(out.status, 'partial', 'downgraded to partial (renders Not estimated)');
    expect(out.provenance.lambdaZ, null, 'λz provenance nulled');
    expect(Number.isNaN(out.values.halfLife), true, 'λz-derived t½ suppressed');
    expect(Number.isNaN(out.values.cl), true, 'CL suppressed');
    expect(Number.isNaN(out.values.aucInf), true, 'AUCinf suppressed');
    expect(out.values.cmax, 10, 'Cmax preserved (λz-independent)');
    expect(out.values.aucLast, 50, 'AUClast preserved');
    expect(out.values.tlag, 0, 'Tlag preserved');
    expect(out.provenance.warnings.some((w) => w.code === LAMBDAZ_BELOW_MIN_POINTS), true,
      'floor warning attached');
  });

  test('applyMinPointsFloor: 3-point fit meets floor → unchanged (same object)', async () => {
    const r = resultWithLz(3);
    expect(applyMinPointsFloor(r, 3), r, 'floor met → input returned unchanged');
  });

  test('applyMinPointsFloor: already-null λz → unchanged', async () => {
    const base = resultWithLz(3);
    const r: nca.ComputeResult = {
      status: 'partial', values: base.values,
      provenance: {...base.provenance, lambdaZ: null},
    };
    expect(applyMinPointsFloor(r, 3), r, 'no λz → input returned unchanged');
  });

  test('isLambdaZNonContiguous: contiguous run → false, interior gap → true (decision 2)', async () => {
    const aug: AugmentedProfile = {
      time: Float64Array.from([0, 1, 2, 4, 8]),
      conc: Float64Array.from([10, 5, 2.5, 0.625, 0.039]),
      blqMask: Uint8Array.from([0, 0, 0, 0, 0]),
      rowId: Int32Array.from([-1, 1, 2, 3, 4]),
    };
    // Indices 2,3,4 (t=2,4,8) — nothing measurable skipped inside → contiguous.
    expect(isLambdaZNonContiguous(Int32Array.from([2, 3, 4]), aug), false, 'contiguous terminal run');
    // Indices 1,3,4 (t=1,4,8) — the measurable t=2 (index 2) sits inside and is skipped.
    expect(isLambdaZNonContiguous(Int32Array.from([1, 3, 4]), aug), true, 'interior measurable point skipped');
    // A 1-point (or empty) fit cannot have an interior gap.
    expect(isLambdaZNonContiguous(Int32Array.from([3]), aug), false, 'single-point fit');
  });

  test('isLambdaZNonContiguous: a BLQ / ≤0 interior point is NOT counted as a gap', async () => {
    const aug: AugmentedProfile = {
      time: Float64Array.from([1, 2, 3, 4]),
      conc: Float64Array.from([10, 0, 8, 4]), // t=2 is BLQ (conc 0) → not measurable
      blqMask: Uint8Array.from([0, 1, 0, 0]),
      rowId: Int32Array.from([1, 2, 3, 4]),
    };
    // Fit uses t=1,3,4 (indices 0,2,3); the skipped interior t=2 is non-measurable.
    expect(isLambdaZNonContiguous(Int32Array.from([0, 2, 3]), aug), false,
      'a non-measurable interior point is a legitimate skip, not a contiguity gap');
  });

  test('applyMinPointsFloor: ≥3-point fit with NON-POSITIVE λz → Not estimated (AC-4.3 λz≤0)', async () => {
    const base = resultWithLz(5); // 5 points → meets the point floor
    // force a non-decaying (rising / absorption-phase) manual selection.
    const rising: nca.ComputeResult = {
      ...base,
      values: {...base.values, lambdaZ: -0.05, halfLife: -13.86},
      provenance: {...base.provenance, lambdaZ: {...base.provenance.lambdaZ!, lambdaZ: -0.05}},
    };
    const out = applyMinPointsFloor(rising, 3);
    expect(out.status, 'partial', 'non-decaying manual fit → Not estimated (not a trusted result)');
    expect(out.provenance.lambdaZ, null, 'λz provenance nulled');
    expect(Number.isNaN(out.values.lambdaZ), true, 'λz suppressed');
    expect(Number.isNaN(out.values.halfLife), true, 'no negative t½ leaks');
    expect(Number.isNaN(out.values.cl), true, 'CL suppressed');
    expect(out.values.cmax, 10, 'Cmax preserved (λz-independent)');
    expect(out.provenance.warnings.some((w) => w.code === LAMBDAZ_BELOW_MIN_POINTS), true,
      'warning attached');
  });
});

// ── formatOverrideTimestamp (λz-override who/when surfacing) ──────────────────
category('Adapter: formatOverrideTimestamp', () => {
  test('blank input → blank (no override)', async () => {
    expect(formatOverrideTimestamp(''), '', 'empty string stays empty');
  });

  test('valid ISO → local YYYY-MM-DD HH:MM (timezone-independent round-trip)', async () => {
    // Build a LOCAL date, serialize to UTC ISO (as writeLzActions does), and confirm the
    // formatter recovers the ORIGINAL local components — independent of the test machine's zone.
    const d = new Date(2026, 8, 8, 14, 30); // local 2026-09-08 14:30 (month is 0-based)
    expect(formatOverrideTimestamp(d.toISOString()), '2026-09-08 14:30', 'local components recovered');
  });

  test('zero-pads single-digit month/day/hour/minute', async () => {
    const d = new Date(2026, 0, 3, 4, 5); // local 2026-01-03 04:05
    expect(formatOverrideTimestamp(d.toISOString()), '2026-01-03 04:05', 'all fields padded to 2 digits');
  });

  test('present-but-unparseable stamp → raw string (visible, not silently dropped)', async () => {
    expect(formatOverrideTimestamp('not-a-date'), 'not-a-date', 'garbage is shown, not blanked');
  });
});
