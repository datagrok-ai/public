/**
 * Unit tests for the NONMEM event-driven normalization transform
 * (`adapter/event-driven.ts`). Synthetic, runtime-free — they exercise the
 * classification / ladder / TAD / route / EVID-SS logic in isolation, without a
 * Datagrok UI or the PKNCA oracle (the end-to-end oracle parity lives in
 * `reference-suite-tests.ts`, Slice E).
 */

import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import {nca} from '@datagrok-libraries/sci-comp';

import {
  normalizeEventDrivenObservations, expandDoseLadder, assertDoseLadder,
  TAD_COL, CONC_COL, DOSE_COL, SUBJECT_COL, OCCASION_COL,
  NONMEM_ROUTE_INFERRED, EVID_RESET_NEW_EPISODE, EVID_OTHER_IGNORED,
  NONMEM_STEADY_STATE_ASSUMED, NONMEM_INFUSION_RATE_INFERRED,
  type EventContract,
} from '../adapter/event-driven';
import {DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';

// ── Synthetic NONMEM DataFrame builder ──────────────────────────────────────
interface NmSpec {
  readonly id: ReadonlyArray<string>;
  readonly time: ReadonlyArray<number>;
  readonly dv: ReadonlyArray<number>; // NaN = missing (dose rows)
  readonly amt: ReadonlyArray<number>; // NaN = missing (obs rows)
  readonly evid: ReadonlyArray<number>;
  readonly mdv?: ReadonlyArray<number>;
  readonly cmt?: ReadonlyArray<number>; // INT_NULL = missing
  readonly rate?: ReadonlyArray<number>;
  readonly addl?: ReadonlyArray<number>; // INT_NULL = missing
  readonly ii?: ReadonlyArray<number>;
  readonly ss?: ReadonlyArray<number>; // INT_NULL = missing
  readonly route?: ReadonlyArray<string>;
  readonly day?: ReadonlyArray<number>; // INT_NULL = missing
  readonly bw?: ReadonlyArray<number>;
  readonly sex?: ReadonlyArray<string>;
}

function nmDf(s: NmSpec): DG.DataFrame {
  const cols: DG.Column[] = [
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', Array.from(s.id)),
    DG.Column.fromFloat64Array('TIME', Float64Array.from(s.time)),
    DG.Column.fromFloat64Array('DV', Float64Array.from(s.dv)),
    DG.Column.fromFloat64Array('AMT', Float64Array.from(s.amt)),
    DG.Column.fromInt32Array('EVID', Int32Array.from(s.evid)),
  ];
  if (s.mdv !== undefined) cols.push(DG.Column.fromInt32Array('MDV', Int32Array.from(s.mdv)));
  if (s.cmt !== undefined) cols.push(DG.Column.fromInt32Array('CMT', Int32Array.from(s.cmt)));
  if (s.rate !== undefined) cols.push(DG.Column.fromFloat64Array('RATE', Float64Array.from(s.rate)));
  if (s.addl !== undefined) cols.push(DG.Column.fromInt32Array('ADDL', Int32Array.from(s.addl)));
  if (s.ii !== undefined) cols.push(DG.Column.fromFloat64Array('II', Float64Array.from(s.ii)));
  if (s.ss !== undefined) cols.push(DG.Column.fromInt32Array('SS', Int32Array.from(s.ss)));
  if (s.route !== undefined) cols.push(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ROUTE', Array.from(s.route)));
  if (s.day !== undefined) cols.push(DG.Column.fromInt32Array('DAY', Int32Array.from(s.day)));
  if (s.bw !== undefined) cols.push(DG.Column.fromFloat64Array('BW_kg', Float64Array.from(s.bw)));
  if (s.sex !== undefined) cols.push(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', Array.from(s.sex)));
  return DG.DataFrame.fromColumns(cols);
}

const BASE_EC: EventContract = {
  id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID',
};

const NUL = DG.INT_NULL;

/** Collect (TAD, conc) pairs for a given pseudo-subject key. */
function pointsFor(df: DG.DataFrame, key: string): Array<{tad: number; conc: number; occ: number}> {
  const subj = df.col(SUBJECT_COL)!; const tad = df.col(TAD_COL)!;
  const conc = df.col(CONC_COL)!; const occ = df.col(OCCASION_COL)!;
  const out: Array<{tad: number; conc: number; occ: number}> = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (String(subj.get(i)) === key)
      out.push({tad: tad.get(i) as number, conc: conc.get(i) as number, occ: occ.get(i) as number});
  }
  return out;
}

category('Event-driven: ADDL/II ladder (AC-B4/B5)', () => {
  test('expandDoseLadder(0, 27, 24) → 28 doses, [0]=0, [27]=648', async () => {
    const l = expandDoseLadder(0, 27, 24);
    expect(l.length, 28, 'ladder length = 1 + ADDL');
    expect(l[0], 0, 'ladder starts at the dose time');
    expect(l[27], 648, 'last dose at 27·II = 648');
  });

  test('expandDoseLadder non-zero t_dose starts at t_dose', async () => {
    const l = expandDoseLadder(12, 2, 24);
    expect(l.length, 3, '1 + ADDL = 3');
    expect(l[0], 12, 'starts at t_dose=12');
    expect(l[2], 60, '12 + 2·24 = 60');
  });

  test('assertDoseLadder rejects a missing-initial-dose ladder (off-by-one, R2 N-1)', async () => {
    const missingInitial: number[] = [];
    for (let k = 1; k <= 27; k++) missingInitial.push(k * 24); // 24..648, length 27
    let threw = false;
    try {assertDoseLadder(missingInitial, 0, 27);} catch {threw = true;}
    expect(threw, true, 'length 27 ≠ 1 + 27 must throw');
  });

  test('assertDoseLadder rejects a phantom-dose ladder (length 1+ADDL+1)', async () => {
    const phantom: number[] = [];
    for (let k = 0; k <= 28; k++) phantom.push(k * 24); // 0..672, length 29
    let threw = false;
    try {assertDoseLadder(phantom, 0, 27);} catch {threw = true;}
    expect(threw, true, 'length 29 ≠ 1 + 27 must throw');
  });

  test('AC-B5: explicit second dose row is union-merged, not re-expanded', async () => {
    // Dose1 ADDL=3/II=24 at t=0 → {0,24,48,72}; explicit dose2 at t=12 (no ADDL).
    // Union = 5 distinct times, NOT a 4×N blow-up.
    const df = nmDf({
      id: ['s', 's', 's'],
      time: [0, 12, 100], dv: [NaN, NaN, 5], amt: [1, 1, NaN], evid: [1, 1, 0],
      addl: [3, NUL, NUL], ii: [24, NaN, NaN],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, addl: 'ADDL', ii: 'II'});
    const ladder = res.doseLadderBySubject.get('s')!;
    expect(ladder.length, 5, 'union of {0,24,48,72} and {12} = 5 doses (no re-expansion)');
  });
});

category('Event-driven: TAD assignment (AC-B3/B11)', () => {
  test('AC-B3: obs at TIME=672 with II=24 → TAD=24=II retained (closed interval)', async () => {
    const df = nmDf({
      id: ['a', 'a'],
      time: [0, 672], dv: [NaN, 8.7], amt: [2.87, NaN], evid: [1, 0],
      addl: [27, NUL], ii: [24, NaN], day: [28, 28],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, addl: 'ADDL', ii: 'II', day: 'DAY'});
    // Multi-dose ladder (ADDL=27) ⇒ occasion-scoped key (AD-8 / FM-NEW-2).
    const pts = pointsFor(res.observations, 'a_Day28');
    expect(pts.length, 1, 'the TIME=672 boundary obs is retained');
    expectFloat(pts[0].tad, 24, 1e-9, 'TAD = 672 − 648 = 24 = II (closed right)');
    expect(pts[0].occ, 28, 'DAY occasion carried');
  });

  test('Day-28 sample maps abs time → TAD on the 28th-dose anchor', async () => {
    const df = nmDf({
      id: ['a', 'a'],
      time: [0, 648.5], dv: [NaN, 100], amt: [2.87, NaN], evid: [1, 0],
      addl: [27, NUL], ii: [24, NaN], day: [28, 28],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, addl: 'ADDL', ii: 'II', day: 'DAY'});
    // Multi-dose ladder (ADDL=27) ⇒ occasion-scoped key (AD-8 / FM-NEW-2).
    const pts = pointsFor(res.observations, 'a_Day28');
    expectFloat(pts[0].tad, 0.5, 1e-9, 'TAD = 648.5 − 648 = 0.5');
    expect(pts[0].occ, 28, 'DAY occasion carried');
  });

  test('AC-B11: post-EVID=4 (reset+dose) obs anchors on the EVID=4 dose', async () => {
    // Dose at t=0, obs at t=2 (TAD 2); EVID=4 reset+dose at t=10, obs at t=12.
    const df = nmDf({
      id: ['x', 'x', 'x', 'x'],
      time: [0, 2, 10, 12], dv: [NaN, 5, NaN, 7], amt: [1, NaN, 1, NaN], evid: [1, 0, 4, 0],
    });
    const res = normalizeEventDrivenObservations(df, BASE_EC);
    // Two episodes (the reset splits them): _ep0 and _ep1.
    const ep0 = pointsFor(res.observations, 'x_ep0');
    const ep1 = pointsFor(res.observations, 'x_ep1');
    expect(ep0.length, 1, 'pre-reset obs in episode 0');
    expectFloat(ep0[0].tad, 2, 1e-9, 'pre-reset TAD = 2 − 0');
    expect(ep1.length, 1, 'post-reset obs in episode 1');
    expectFloat(ep1[0].tad, 2, 1e-9, 'post-EVID=4 TAD = 12 − 10 (anchored on the EVID=4 dose)');
  });
});

category('Event-driven: dose-of-record + vehicle control (AC-B7/B12)', () => {
  test('AC-B7: per-row dose = AMT exactly (not re-derived)', async () => {
    const df = nmDf({
      id: ['r', 'r'],
      time: [0, 0.25], dv: [NaN, 295.82], amt: [0.2441, NaN], evid: [1, 0],
    });
    const res = normalizeEventDrivenObservations(df, BASE_EC);
    const dose = res.observations.col(DOSE_COL)!.get(0) as number;
    expectFloat(dose, 0.2441, 1e-12, 'dose = AMT of record');
  });

  test('AC-B12: vehicle control AMT=0 classified as a normal dose, not dropped', async () => {
    const df = nmDf({
      id: ['v', 'v', 'v'],
      time: [0, 0.5, 2], dv: [NaN, 0, 0], amt: [0, NaN, NaN], evid: [1, 0, 0],
      addl: [27, NUL, NUL], ii: [24, NaN, NaN], day: [1, 1, 1],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, addl: 'ADDL', ii: 'II', day: 'DAY'});
    // Multi-dose ladder (ADDL=27) ⇒ occasion-scoped key (AD-8 / FM-NEW-2).
    const pts = pointsFor(res.observations, 'v_Day1');
    expect(pts.length, 2, 'AMT=0 control obs are retained');
    expect(res.observations.col(DOSE_COL)!.get(0) as number, 0, 'dose carried as 0 (not NaN/dropped)');
    expect(res.classificationCounts.dose, 1, 'the AMT=0 row counts as a dose event');
  });
});

category('Event-driven: EVID 2/3 + SS handling (AC-D1/D2/D3)', () => {
  test('AC-D1: EVID=3 reset → two episodes, pre/post AUCs not pooled', async () => {
    const df = nmDf({
      id: ['p', 'p', 'p', 'p', 'p'],
      time: [0, 1, 5, 6, 7], dv: [NaN, 9, NaN, NaN, 4], amt: [1, NaN, NaN, 1, NaN],
      evid: [1, 0, 3, 1, 0],
    });
    const res = normalizeEventDrivenObservations(df, BASE_EC);
    const ep0 = pointsFor(res.observations, 'p_ep0');
    const ep1 = pointsFor(res.observations, 'p_ep1');
    expect(ep0.length, 1, 'episode 0 obs');
    expect(ep1.length, 1, 'episode 1 obs (post-reset)');
    const warns = res.warningsBySubject.get('p') ?? [];
    expect(warns.some((w) => w.code === EVID_RESET_NEW_EPISODE), true, 'reset warning present');
  });

  test('AC-D2: EVID=2 dropped from conc stream + info warning, no obs row', async () => {
    const df = nmDf({
      id: ['q', 'q', 'q'],
      time: [0, 1, 2], dv: [NaN, 9, 5], amt: [1, NaN, NaN], evid: [1, 2, 0],
    });
    const res = normalizeEventDrivenObservations(df, BASE_EC);
    const pts = pointsFor(res.observations, 'q');
    expect(pts.length, 1, 'only the EVID=0 row becomes an obs');
    expect(res.classificationCounts.other, 1, 'EVID=2 counted as other');
    const warns = res.warningsBySubject.get('q') ?? [];
    expect(warns.some((w) => w.code === EVID_OTHER_IGNORED), true, 'EVID_OTHER_IGNORED info warning');
  });

  test('AC-D3: SS=1 + II honored; SS=1 no II rejected; SS=2 rejected', async () => {
    // SS=1 with II present → honor warning, subject computed.
    const honored = normalizeEventDrivenObservations(
      nmDf({id: ['h', 'h'], time: [0, 1], dv: [NaN, 9], amt: [1, NaN], evid: [1, 0],
        ss: [1, NUL], ii: [24, NaN]}),
      {...BASE_EC, ss: 'SS', ii: 'II'});
    expect(honored.errors.length, 0, 'SS=1 + II is honored (no blocking error)');
    expect((honored.warningsBySubject.get('h') ?? []).some((w) => w.code === NONMEM_STEADY_STATE_ASSUMED),
      true, 'NONMEM_STEADY_STATE_ASSUMED warning');

    // SS=1 without II → blocking reject.
    const noII = normalizeEventDrivenObservations(
      nmDf({id: ['n', 'n'], time: [0, 1], dv: [NaN, 9], amt: [1, NaN], evid: [1, 0], ss: [1, NUL]}),
      {...BASE_EC, ss: 'SS'});
    expect(noII.errors.some((e) => e.code === 'NONMEM_SS_UNSUPPORTED'), true, 'SS=1 no II rejected');

    // SS=2 (superposition) → always rejected, distinct code.
    const ss2 = normalizeEventDrivenObservations(
      nmDf({id: ['z', 'z'], time: [0, 1], dv: [NaN, 9], amt: [1, NaN], evid: [1, 0],
        ss: [2, NUL], ii: [24, NaN]}),
      {...BASE_EC, ss: 'SS', ii: 'II'});
    expect(ss2.errors.some((e) => e.code === 'NONMEM_SS_SUPERPOSITION_UNSUPPORTED'), true,
      'SS=2 rejected with the superposition code regardless of II');
  });
});

category('Event-driven: route resolution (AC-C1/C2/C4)', () => {
  test('AC-C1: explicit ROUTE column resolves without a confirm gate', async () => {
    const df = nmDf({
      id: ['R001', 'R001'], time: [0, 0.25], dv: [NaN, 295.82], amt: [0.2441, NaN], evid: [1, 0],
      cmt: [1, 1], route: ['IV', 'IV'],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, cmt: 'CMT', route: 'ROUTE'});
    expect(res.routeConfirmRequired, false, 'explicit ROUTE needs no confirm');
    expect(res.observations.col(DERIVED_ROUTE_CODE_COL)!.get(0) as number, nca.ROUTE_IV_BOLUS,
      'bare IV → IV-bolus code');
    const inf = res.inferredRoutes.find((r) => r.subjectKey === 'R001')!;
    expect(inf.inferred, false, 'route not flagged inferred');
  });

  test('AC-C2: RATE/CMT inference (no ROUTE col) → IV-bolus from central CMT + confirm required', async () => {
    const df = nmDf({
      id: ['s', 's'], time: [0, 0.25], dv: [NaN, 295], amt: [0.2441, NaN], evid: [1, 0],
      cmt: [1, 1], rate: [0, NaN],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, cmt: 'CMT', rate: 'RATE'});
    expect(res.observations.col(DERIVED_ROUTE_CODE_COL)!.get(0) as number, nca.ROUTE_IV_BOLUS,
      'central CMT + RATE=0 → IV-bolus');
    expect(res.routeConfirmRequired, true, 'inferred route requires confirm');
    expect((res.warningsBySubject.get('s') ?? []).some((w) => w.code === NONMEM_ROUTE_INFERRED),
      true, 'NONMEM_ROUTE_INFERRED warning');
  });

  test('AC-C4: depot CMT does NOT silently invert to IV — maps to extravascular (other)', async () => {
    const df = nmDf({
      id: ['d', 'd'], time: [0, 1], dv: [NaN, 50], amt: [1, NaN], evid: [1, 0],
      cmt: [2, 2], rate: [0, NaN],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, cmt: 'CMT', rate: 'RATE'});
    expect(res.observations.col(DERIVED_ROUTE_CODE_COL)!.get(0) as number, nca.ROUTE_OTHER,
      'depot CMT → extravascular (other), never a fixed CMT→route inversion');
  });

  test('AC-C3: RATE>0 → IV-infusion with duration + inferred (info) warning', async () => {
    const df = nmDf({
      id: ['i', 'i'], time: [0, 1], dv: [NaN, 50], amt: [10, NaN], evid: [1, 0],
      cmt: [1, 1], rate: [5, NaN],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, cmt: 'CMT', rate: 'RATE'});
    expect(res.observations.col(DERIVED_ROUTE_CODE_COL)!.get(0) as number, nca.ROUTE_IV_INFUSION,
      'RATE>0 → IV-infusion');
    const inf = res.inferredRoutes.find((r) => r.subjectKey === 'i')!;
    expectFloat(inf.infusionDuration!, 2, 1e-9, 'duration = AMT/RATE = 10/5 = 2');
    // BD-1 renamed NONMEM_INFUSION_RATE_UNVALIDATED → _INFERRED (severity info).
    const infWarn = (res.warningsBySubject.get('i') ?? []).find((w) => w.code === NONMEM_INFUSION_RATE_INFERRED);
    expect(infWarn !== undefined, true, 'NONMEM_INFUSION_RATE_INFERRED warning present');
    expect(infWarn!.severity, 'info', 'infusion-rate-inferred warning is informational (BD-1)');
  });

  test('confirmedRoute overrides inference and clears the confirm gate', async () => {
    const df = nmDf({
      id: ['s', 's'], time: [0, 1], dv: [NaN, 50], amt: [1, NaN], evid: [1, 0], cmt: [2, 2],
    });
    const res = normalizeEventDrivenObservations(
      df, {...BASE_EC, cmt: 'CMT'}, {confirmedRoute: nca.ROUTE_PO});
    expect(res.observations.col(DERIVED_ROUTE_CODE_COL)!.get(0) as number, nca.ROUTE_PO,
      'confirmed PO applied');
    expect(res.routeConfirmRequired, false, 'confirm satisfied');
  });
});

category('Event-driven: occasion keying + covariate hoist (AC-B8/B9)', () => {
  test('AC-B9: same animal Day1 + Day28 → two profiles; single-occasion keeps bare id', async () => {
    // Multi-occasion: DAY ∈ {1,28}.
    const multi = normalizeEventDrivenObservations(
      nmDf({
        id: ['A', 'A', 'A', 'A'],
        time: [0, 1, 648, 649], dv: [NaN, 5, NaN, 6], amt: [1, NaN, 1, NaN], evid: [1, 0, 1, 0],
        addl: [27, NUL, NUL, NUL], ii: [24, NaN, NaN, NaN], day: [1, 1, 28, 28], bw: [0.25, 0.25, 0.30, 0.30],
      }),
      {...BASE_EC, addl: 'ADDL', ii: 'II', day: 'DAY', bodyWeight: 'BW_kg'});
    const day1 = pointsFor(multi.observations, 'A_Day1');
    const day28 = pointsFor(multi.observations, 'A_Day28');
    expect(day1.length, 1, 'Day1 profile');
    expect(day28.length, 1, 'Day28 profile (two profiles, not collapsed)');

    // Single-occasion (single dose, no DAY variation) → bare id.
    const single = normalizeEventDrivenObservations(
      nmDf({id: ['R001', 'R001'], time: [0, 0.25], dv: [NaN, 295], amt: [0.2441, NaN], evid: [1, 0]}),
      BASE_EC);
    expect(pointsFor(single.observations, 'R001').length, 1, 'single-occasion keeps bare subject id');
  });

  test('FM-NEW-2: repeat-dose WITHOUT a DAY column still splits Day-1/Day-28 into 2 profiles', async () => {
    // No DAY column → occasion derived from the ladder index (AC-E3 / F-1). The
    // same animal sampled on Day-1 (TAD 1 @ TIME 1) and Day-28 (TAD 1 @ TIME 649)
    // must become TWO profiles, else buildProfileIndices merges them (wrong NCA).
    const df = nmDf({
      id: ['A', 'A', 'A'],
      time: [0, 1, 649], dv: [NaN, 5, 6], amt: [1, NaN, NaN], evid: [1, 0, 0],
      addl: [27, NUL, NUL], ii: [24, NaN, NaN],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, addl: 'ADDL', ii: 'II'});
    const day1 = pointsFor(res.observations, 'A_Day1');
    const day28 = pointsFor(res.observations, 'A_Day28');
    expect(day1.length, 1, 'Day-1 obs in its own profile (no DAY column)');
    expect(day28.length, 1, 'Day-28 obs split into a separate profile via the ladder index');
    expectFloat(day1[0].tad, 1, 1e-9, 'Day-1 TAD = 1 − 0');
    expectFloat(day28[0].tad, 1, 1e-9, 'Day-28 TAD = 649 − 648');
  });

  test('AC-B8: hoisted Sex/BW present on every obs row', async () => {
    const df = nmDf({
      id: ['A', 'A', 'A'],
      time: [0, 1, 2], dv: [NaN, 5, 3], amt: [1, NaN, NaN], evid: [1, 0, 0],
      bw: [0.25, 0.25, 0.25], sex: ['M', 'M', 'M'],
    });
    const res = normalizeEventDrivenObservations(df, {...BASE_EC, bodyWeight: 'BW_kg', covariates: ['Sex']});
    expect(res.observations.col('Sex') !== null, true, 'Sex hoisted');
    expect(res.observations.col('BW_kg') !== null, true, 'BW hoisted');
    for (let i = 0; i < res.observations.rowCount; i++) {
      expect(String(res.observations.col('Sex')!.get(i)), 'M', `Sex on obs row ${i}`);
      expectFloat(res.observations.col('BW_kg')!.get(i) as number, 0.25, 1e-9, `BW on obs row ${i}`);
    }
  });
});
