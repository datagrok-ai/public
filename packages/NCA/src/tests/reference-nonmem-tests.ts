/**
 * Reference suite — NONMEM event-driven ingestion parity (Slice E, Rule-18 gate).
 *
 * The transform `normalizeEventDrivenObservations` emits the canonical
 * observations frame; the EXISTING pipeline (`runNcaPipeline`) then runs on it
 * unchanged (AD-1). Two PKNCA NONMEM cards are the oracle:
 *  - `03_rat_multiroute_nonmem.csv` — serial / dense, 64 subjects, single dose
 *    at t=0, explicit ROUTE (IV/PO/SC) → per-subject parity vs the multiroute card.
 *  - `06_rat_tk28d_nonmem.csv` — destructive-composite 28-day repeat-dose
 *    (ADDL=27/II=24), 360 animals sampled once each → pooled to composite-mean
 *    profiles per dose×sex×day cell, the discriminating ADDL/II-expander test.
 */

import {expectNoRawNan} from './test-utils';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';
import * as DG from 'datagrok-api/dg';

import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, loadReference, loadReferenceCsv, indexFixtureProfiles, runPipelineOrThrow,
  type FixtureProfile,
} from './test-utils';
import type {ColumnContract} from '../adapter/column-contract';
import type {ProfileScalars} from '../adapter/profile-inputs';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {computeAccumulationRatios, computeExposureMargin} from '../adapter/toxicokinetics';
import {DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';
import {
  normalizeEventDrivenObservations, SUBJECT_COL, TAD_COL, CONC_COL, DOSE_COL, OCCASION_COL,
  DERIVED_INFUSION_DURATION_COL, NONMEM_INFUSION_RATE_INFERRED, EVID_RESET_NEW_EPISODE,
  NONMEM_STEADY_STATE_ASSUMED,
  type EventContract,
} from '../adapter/event-driven';
import {poolCompositeByCell, STUDY_DAY_COL, CELL_BW_COL} from '../adapter/composite-pool';

const TIMEOUT = 60000;

/** λz window-selection divergence (sci-comp best-fit vs PKNCA): R019 formerly
 *  selected 8 terminal points where the reference used 4 → λz/t½/Vz ~2.3% off.
 *  RESOLVED in sci-comp 0.9.1 (best-fit now matches PKNCA's flat adj-R² tolerance),
 *  so R019 asserts at strict tolerance. Set kept empty as the parameterized band
 *  for any future regression (VAL-01-LZ-R019). */
const RAT_LZ_DIVERGENT: ReadonlySet<string> = new Set([]);

/** Canonical observations contract over the normalized event-driven frame. AMT
 *  is ABSOLUTE mg (AC-A6) → doseUnits 'mg'; conc carries the source ng/mL. */
function nonmemNormalizedContract(
  withBw: boolean, withDose: boolean, concUnit: 'ng/mL' | 'mg/L' = 'ng/mL',
): ColumnContract {
  const mut: {-readonly [K in keyof ColumnContract['observations']]: ColumnContract['observations'][K]} = {
    subjectId: {name: SUBJECT_COL, required: true, rawType: 'string'},
    time: {name: TAD_COL, required: true, rawType: 'float64', units: 'h'},
    concentration: {name: CONC_COL, required: true, rawType: 'float64', units: concUnit},
    route: {name: DERIVED_ROUTE_CODE_COL, required: false, rawType: 'int32'},
    // Slice A2 / BD-2: the transform always emits `~InfusionDuration` (NaN for
    // RATE=0 subjects → extractor falls back to the scalar null, so 03/06 stay
    // byte-unchanged). Referencing it threads each RATE>0 subject's own duration.
    infusionDuration: {name: DERIVED_INFUSION_DURATION_COL, required: false, rawType: 'float64', units: 'h'},
  };
  if (withDose) mut.dose = {name: DOSE_COL, required: false, rawType: 'float64', units: 'mg'};
  if (withBw) mut.bodyWeight = {name: 'BW_kg', required: false, rawType: 'float64', units: 'kg'};
  return {observations: mut};
}

function nonmemScalars(concUnit: 'ng/mL' | 'mg/L' = 'ng/mL'): ProfileScalars {
  // `concentrationUnits` (NOT the contract's `units` field) drives the conc
  // rescale in scaleProfileInputs — the synthetic analytic fixtures (1/2/3/4/6/7)
  // carry mg/L directly (dose mg / AUC → CL), so they must pass 'mg/L' here;
  // 03/06 carry raw ng/mL (default).
  return {dose: 0, route: nca.ROUTE_OTHER, doseUnits: 'mg', concentrationUnits: concUnit,
    timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0};
}

/** Per-subject parity (cmax 0.1%, tmax exact, auclast/aucinf 0.1%, cl 1%,
 *  λz/t½/vz strict 0.5%/1%). The `RAT_LZ_DIVERGENT` set widens to 3% per listed
 *  subject; currently empty (VAL-01-LZ-R019 resolved in sci-comp 0.9.1).
 *  The NONMEM cards dose AMT in mg, so their cl / vz are ABSOLUTE and are asserted on the
 *  companions (nca-cl-bw-normalised-units probe BREAK 2). `bwKg` — the subject's weight as recorded
 *  in the fixture — adds the per-kg IDENTITY (PKNCA ÷ BW × 1000/60); without it the per-kg value is
 *  not formed (an absolute dose and no weight). */
function ratRowAssert(results: DG.DataFrame, i: number, subj: string, e: FixtureProfile, tag: string,
  bwKg?: number): void {
  const get = (n: string) => results.col(n)!.get(i) as number;
  const p = e.parameters;
  if (p.cmax !== null) expectFloat(get('CMAX'), p.cmax, Math.abs(p.cmax) * 0.001, `${tag}: cmax`);
  if (p.tmax !== null) expect(get('TMAX'), p.tmax, `${tag}: tmax`);
  if (p.auclast !== null) expectFloat(get('AUC_LAST'), p.auclast, Math.abs(p.auclast) * 0.001, `${tag}: auclast`);
  if (p.aucinf !== null) expectFloat(get('AUC_INF'), p.aucinf, Math.abs(p.aucinf) * 0.001, `${tag}: aucinf`);
  if (p.cl !== null) {
    expectFloat(get('CL_ABS'), p.cl, Math.abs(p.cl) * 0.01, `${tag}: cl abs.`);
    if (bwKg !== undefined) {
      const want = p.cl / bwKg * 1000 / 60;
      expectFloat(get('CL'), want, want * 0.01, `${tag}: cl per kg (IDENTITY)`);
    } else
      expect(results.col('CL')!.isNone(i), true, `${tag}: absolute dose, no BW → no per-kg cl`);
  }
  if (p.pct_aucextrap !== null) expectFloat(get('PCT_EXTRAP'), p.pct_aucextrap, 0.5, `${tag}: %extrap`);
  const lzTol = RAT_LZ_DIVERGENT.has(subj) ? 0.03 : 0.005;
  const derivedTol = RAT_LZ_DIVERGENT.has(subj) ? 0.03 : 0.01;
  if (p.lambda_z !== null) expectFloat(get('LAMBDA_Z'), p.lambda_z, Math.abs(p.lambda_z) * lzTol, `${tag}: λz`);
  if (p.half_life !== null)
    expectFloat(get('HALF_LIFE'), p.half_life, Math.abs(p.half_life) * derivedTol, `${tag}: t½`);
  if (p.vz !== null) {
    expectFloat(get('VZ_ABS'), p.vz, Math.abs(p.vz) * derivedTol, `${tag}: Vz abs.`);
    if (bwKg !== undefined)
      expectFloat(get('VZ'), p.vz / bwKg, Math.abs(p.vz / bwKg) * derivedTol, `${tag}: Vz per kg (IDENTITY)`);
  }
}

category('Reference suite: NONMEM event-driven (rat multiroute 03)', () => {
  const ec03: EventContract = {
    id: 'SubjectID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID',
    cmt: 'CMT', rate: 'RATE', mdv: 'MDV', route: 'ROUTE', bodyWeight: 'BW_kg',
    covariates: ['Sex', 'DOSE_MGKG'],
  };

  test('AC-B1/AC-C1/AC-E1: 03 NONMEM → 64 subjects within tolerances vs PKNCA multiroute card', async () => {
    const df = await loadFixture('03_rat_multiroute_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec03);
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);
    expect(norm.routeConfirmRequired, false, 'explicit ROUTE column → no confirm gate (AC-C1)');

    const contract = nonmemNormalizedContract(true, true);
    const {results, outcome} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars(), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    expect(results.rowCount, 64, `expected 64 subjects; got ${results.rowCount}`);

    // AC-B8 / AC-B10 (R1 S-2, blocking): the covariate hoist survives into the
    // results frame via `appendGroupingColumns`, fed the NORMALIZED frame. This
    // is the concrete four-consumer-threading gate: the profiles index the
    // 664-obs-row normalized frame, so passing the raw 729-row NONMEM df here
    // would throw a rowCount mismatch (grouping-columns.ts). It does not — and
    // the carried Sex equals the normalized obs-row value (R001 = 'M').
    const profiles = outcome.perProfile.map((p) => p.index);
    const carried = appendGroupingColumns(results, norm.observations, profiles, contract);
    expect(carried.includes('Sex'), true, `Sex carried into results; got ${carried.join(',')}`);
    let r001row = -1;
    for (let i = 0; i < results.rowCount; i++)
      if (String(results.col('SUBJECT')!.get(i)) === 'R001') {r001row = i; break;}
    expect(String(results.col('Sex')!.get(r001row)), 'M', 'R001 carried Sex = M (normalized obs-row value)');

    const fx = await loadReferenceCsv('rat_multiroute_nca_datagrok_match.csv');
    const byKey = indexFixtureProfiles(fx);
    // Each animal's weight as recorded in the NONMEM file (`BW_kg`), the per-kg IDENTITY input.
    const bwBySubject = new Map<string, number>();
    for (let r = 0; r < df.rowCount; r++)
      bwBySubject.set(String(df.col('SubjectID')!.get(r)), df.col('BW_kg')!.get(r) as number);
    let r001checked = false;
    for (let i = 0; i < results.rowCount; i++) {
      const subj = String(results.col('SUBJECT')!.get(i));
      const expected = byKey.get(subj);
      expect(expected !== undefined, true, `subject '${subj}' missing from PKNCA card`);
      expect(String(results.col('STATUS')!.get(i)), 'ok', `${subj}: STATUS must be ok`);
      ratRowAssert(results, i, subj, expected!, `03nm ${subj}`, bwBySubject.get(subj));
      if (subj === 'R001') {
        // AC-B1 (the discriminating co-timed t=0 dose+obs split): the dose row is
        // dropped, the first point is the t=0 obs DV=235.16; cmax=0.29582 @ TAD=0.25.
        expectFloat(results.col('CMAX')!.get(i) as number, 0.29582, 0.29582 * 0.001, 'R001 cmax (AC-B1)');
        expect(results.col('TMAX')!.get(i) as number, 0.25, 'R001 tmax (AC-B1)');
        expectFloat(results.col('AUC_LAST')!.get(i) as number, 2.74253917, 2.7425 * 0.001, 'R001 auclast (AC-B1)');
        expect(String(results.col('CL_BASIS')!.get(i)), 'CL', 'R001 IV → CL basis (AC-C1)');
        r001checked = true;
      }
    }
    expect(r001checked, true, 'R001 asserted (AC-B1)');
  }, {timeout: TIMEOUT});
});

category('Reference suite: NONMEM event-driven (rat TK28d 06 composite)', () => {
  const ec06 = (withDay: boolean): EventContract => {
    const base: EventContract = {
      id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
      addl: 'ADDL', ii: 'II', route: 'ROUTE', bodyWeight: 'BW_kg', covariates: ['Sex', 'DOSE_MGKG'],
    };
    return withDay ? {...base, day: 'DAY'} : base;
  };

  /** Pool the normalized per-animal frame into composite-mean profiles, one per
   *  (dose × sex × occasion) cell — the destructive-composite design's pooling.
   *  Emits `Subject` = the cell key (`D{dose}_{sex}_Day{occasion}`). */
  function poolByCell(norm: DG.DataFrame): DG.DataFrame {
    const tad = norm.col(TAD_COL)!; const conc = norm.col(CONC_COL)!;
    const dose = norm.col('DOSE_MGKG')!; const sex = norm.col('Sex')!;
    const occ = norm.col(OCCASION_COL)!; const route = norm.col(DERIVED_ROUTE_CODE_COL)!;
    const groups = new Map<string, {sum: number; nReal: number; tad: number; cell: string; route: number}>();
    for (let i = 0; i < norm.rowCount; i++) {
      const d = Math.round(dose.get(i) as number);
      const cell = `D${d}_${String(sex.get(i))}_Day${occ.get(i) as number}`;
      const t = tad.get(i) as number;
      const key = `${cell}|${t}`;
      let g = groups.get(key);
      if (g === undefined) {g = {sum: 0, nReal: 0, tad: t, cell, route: route.get(i) as number}; groups.set(key, g);}
      const c = conc.get(i) as number;
      if (Number.isFinite(c)) {g.sum += c; g.nReal++;}
    }
    const subj: string[] = []; const ot: number[] = []; const oc: number[] = []; const orr: number[] = [];
    for (const g of groups.values()) {
      subj.push(g.cell); ot.push(g.tad); oc.push(g.nReal > 0 ? g.sum / g.nReal : NaN); orr.push(g.route);
    }
    return DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, SUBJECT_COL, subj),
      DG.Column.fromFloat64Array(TAD_COL, Float64Array.from(ot)),
      DG.Column.fromFloat64Array(CONC_COL, Float64Array.from(oc)),
      DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, Int32Array.from(orr)),
    ]);
  }

  function cellRow(results: DG.DataFrame, cell: string): number {
    for (let i = 0; i < results.rowCount; i++)
      if (String(results.col('SUBJECT')!.get(i)) === cell) return i;
    return -1;
  }

  function day28Points(pooled: DG.DataFrame, cell: string): number {
    let n = 0;
    for (let i = 0; i < pooled.rowCount; i++)
      if (String(pooled.col(SUBJECT_COL)!.get(i)) === cell) n++;
    return n;
  }

  // The 16 reference cells (4 dose × 2 sex × 2 day); the vehicle control (0 mg/kg)
  // has no card row and is excluded (AC-B12 / AC-E1, verified vs the card).
  const CARD_CELLS: ReadonlyArray<{cell: string; auclast: number; cmax: number; tmax: number}> = [
    {cell: 'D10_F_Day1', auclast: 40.4501897648454, cmax: 2.2118, tmax: 2},
    {cell: 'D10_F_Day28', auclast: 219.451089543992, cmax: 10.21428, tmax: 8},
    {cell: 'D10_M_Day1', auclast: 42.7965912153502, cmax: 2.12221, tmax: 4},
    {cell: 'D10_M_Day28', auclast: 249.269352380811, cmax: 12.27822, tmax: 8},
    {cell: 'D30_F_Day1', auclast: 155.138097571997, cmax: 7.23590666666667, tmax: 2},
    {cell: 'D30_F_Day28', auclast: 592.534970440798, cmax: 29.36122, tmax: 8},
    {cell: 'D30_M_Day1', auclast: 144.54290336827, cmax: 6.92223, tmax: 8},
    {cell: 'D30_M_Day28', auclast: 565.277537613813, cmax: 30.7319666666667, tmax: 4},
    {cell: 'D100_F_Day1', auclast: 491.682832050664, cmax: 24.0324466666667, tmax: 1},
    {cell: 'D100_F_Day28', auclast: 2047.95429786658, cmax: 103.010643333333, tmax: 4},
    {cell: 'D100_M_Day1', auclast: 425.4836471084, cmax: 21.08948, tmax: 8},
    {cell: 'D100_M_Day28', auclast: 1991.81109477773, cmax: 98.2592066666667, tmax: 4},
    {cell: 'D300_F_Day1', auclast: 1426.05643143957, cmax: 74.6976, tmax: 4},
    {cell: 'D300_F_Day28', auclast: 6791.63439978575, cmax: 304.529166666667, tmax: 2},
    {cell: 'D300_M_Day1', auclast: 1616.59970414946, cmax: 86.9902166666667, tmax: 8},
    {cell: 'D300_M_Day28', auclast: 6848.64123311154, cmax: 325.604573333333, tmax: 1},
  ];

  test('AC-E1/AC-E2: 06 NONMEM → 16 per-cell composites within tolerances + Day-28 point-count = 6', async () => {
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec06(true));
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);

    const pooled = poolByCell(norm.observations);
    const {results} = runPipelineOrThrow(
      pooled, nonmemNormalizedContract(false, false), nonmemScalars(), DEFAULT_RULES);

    for (const c of CARD_CELLS) {
      const row = cellRow(results, c.cell);
      expect(row >= 0, true, `cell '${c.cell}' present in results`);
      // cmax/auclast carry ng/mL→mg/L rescale → 0.1% (NFR-04); tmax exact.
      expectFloat(results.col('CMAX')!.get(row) as number, c.cmax, Math.abs(c.cmax) * 0.001, `${c.cell} cmax`);
      expectFloat(results.col('AUC_LAST')!.get(row) as number, c.auclast, Math.abs(c.auclast) * 0.001,
        `${c.cell} auclast`);
      expect(results.col('TMAX')!.get(row) as number, c.tmax, `${c.cell} tmax`);
    }

    // AC-E2: the discriminating composite-point count for the Day-28 cell — a
    // missing-initial-dose 5-point ladder would slip an auclast-only check.
    expect(day28Points(pooled, 'D10_F_Day28'), 6,
      'D10_F_Day28 composite has exactly 6 TAD points {0.5,1,2,4,8,24}');
  }, {timeout: TIMEOUT});

  test('AC-E3: no-DAY variant reproduces Day-28 (auclast=219.451 + point-count=6) via the ADDL/II ladder', async () => {
    // Drop the DAY covariate → occasion MUST be derived from the dose-ladder
    // index (28th dose at 648 h), forcing the expander code path (a DAY shortcut
    // would silently bypass the expander while it sits dead).
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec06(false));
    const pooled = poolByCell(norm.observations);
    const {results} = runPipelineOrThrow(
      pooled, nonmemNormalizedContract(false, false), nonmemScalars(), DEFAULT_RULES);

    const row = cellRow(results, 'D10_F_Day28');
    expect(row >= 0, true, 'D10_F_Day28 reconstructed from the ladder index (no DAY column)');
    expectFloat(results.col('AUC_LAST')!.get(row) as number, 219.451089543992, 219.451 * 0.001,
      'no-DAY Day-28 auclast matches the oracle (expander ran)');
    expect(day28Points(pooled, 'D10_F_Day28'), 6, 'no-DAY Day-28 composite still has 6 TAD points');
  }, {timeout: TIMEOUT});

  test('AC-E4: NONMEM composite means equal the serial/composite encoding (06_tk_serial_demo)', async () => {
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec06(true));
    const pooled = poolByCell(norm.observations);

    // Reference composite-mean frame (mg/L), keyed by (PseudoSubject, NominalTime).
    const demo = await loadFixture('06_tk_serial_demo.csv');
    const dSubj = demo.col('Subject')!; const dTime = demo.col('Time_hr')!; const dConc = demo.col('Conc_mg_L')!;
    const demoMean = new Map<string, number>();
    for (let i = 0; i < demo.rowCount; i++)
      demoMean.set(`${String(dSubj.get(i))}|${dTime.get(i) as number}`, dConc.get(i) as number);

    const pSubj = pooled.col(SUBJECT_COL)!; const pTad = pooled.col(TAD_COL)!; const pConc = pooled.col(CONC_COL)!;
    let checked = 0;
    for (let i = 0; i < pooled.rowCount; i++) {
      const cell = String(pSubj.get(i));
      if (cell.startsWith('D0_')) continue; // vehicle control absent from the demo
      const key = `${cell}|${pTad.get(i) as number}`;
      const expected = demoMean.get(key);
      expect(expected !== undefined, true, `demo missing composite point ${key}`);
      // ng/mL → mg/L; identical per-TAD composite means (AC-E4).
      expectFloat((pConc.get(i) as number) / 1000, expected!, Math.max(Math.abs(expected!) * 1e-4, 1e-6),
        `composite mean mismatch at ${key}`);
      checked++;
    }
    expect(checked > 90, true, `enough composite points checked; got ${checked}`);
  }, {timeout: TIMEOUT});

  test('AC-E5 (Rule 18 / IMPL-GAP-2): PRODUCTION poolCompositeByCell reproduces the 16-cell card', async () => {
    // The live mapping pane routes destructive-composite NONMEM through this
    // adapter util — NOT the test-local poolByCell above — so the PRODUCTION
    // pooling code is the thing validated against the per-cell oracle (the
    // cycle-state Rule-18 gate: "run the FULL adapter pool, not the test-local
    // poolByCell"). poolByCell is retained above only as the math reference.
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec06(true));
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);

    const {pooled, cellCount} = poolCompositeByCell(norm.observations, {
      doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex',
    });
    // AC-13 (F4 Arm 1): the pooled composite frame must carry no raw IEEE NaN.
    expectNoRawNan(pooled, 'composite-pool.poolCompositeByCell');
    // 16 active cells (4 dose × 2 sex × 2 day) + 4 vehicle D0 cells = 20.
    expect(cellCount, 20, `expected 20 composite cells (16 active + 4 vehicle); got ${cellCount}`);

    const {results} = runPipelineOrThrow(
      pooled, nonmemNormalizedContract(false, false), nonmemScalars(), DEFAULT_RULES);

    for (const c of CARD_CELLS) {
      const row = cellRow(results, c.cell);
      expect(row >= 0, true, `cell '${c.cell}' present in results`);
      expectFloat(results.col('CMAX')!.get(row) as number, c.cmax, Math.abs(c.cmax) * 0.001, `${c.cell} cmax`);
      expectFloat(results.col('AUC_LAST')!.get(row) as number, c.auclast, Math.abs(c.auclast) * 0.001,
        `${c.cell} auclast`);
      expect(results.col('TMAX')!.get(row) as number, c.tmax, `${c.cell} tmax`);
    }
    expect(day28Points(pooled, 'D10_F_Day28'), 6,
      'D10_F_Day28 composite has exactly 6 TAD points {0.5,1,2,4,8,24}');

    // The pane's TK / dose-prop / F% decorators read these off the pooled frame:
    // StudyDay (the repeat-dose accumulation gate), the NOMINAL DOSE_MGKG group
    // key (stable Day-1↔Day-28), and Sex (stratifier). Their presence is the
    // contract the live composite path depends on.
    expect(pooled.col(STUDY_DAY_COL) !== null, true, 'pooled frame carries StudyDay');
    expect(pooled.col('DOSE_MGKG') !== null, true, 'pooled frame carries the nominal DOSE_MGKG group key');
    expect(pooled.col('Sex') !== null, true, 'pooled frame carries Sex');
  }, {timeout: TIMEOUT});

  // nca-cl-bw-normalised-units F8b / F6 (AD-5b, D-m): the pooled cell carries its mean body weight,
  // so the per-kg CL of a weight-based cell equals the card's per-kg cl.obs — and the exposure
  // margin, which reads per-kg CL, is the card's. At HEAD the pooler nulled BW and the margin read
  // the absolute CL (mean AMT dose) as per kg: 2.4270 against the card's 10.7838 (E6).
  // The card reports PKNCA's best λz fit whatever its R² (D10_F_Day1 adj. R² 0.66), so the cells are
  // run here at the card's rules (no R² gate); at the app's default 0.85 gate that cell has no λz,
  // hence no CL and no projected margin — at HEAD and now alike (measured at build, 2026-10-23).
  const CARD_RULES = {...DEFAULT_RULES, lambdaZ: {...DEFAULT_RULES.lambdaZ, minRSquared: 0}};
  test('F8b / AC-6.1-6.2: 06 NONMEM cell-mean BW → per-kg CL = card cl.obs; MARGIN_PROJECTED asserted', async () => {
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec06(true));
    const {pooled, warningsByCell} = poolCompositeByCell(norm.observations,
      {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex', bodyWeightCol: 'BW_kg'});
    // AC-8b.1: the cell-mean BW on every cell; the per-animal BW stays null on non-constant cells.
    const cellBw = pooled.col(CELL_BW_COL)!;
    for (let i = 0; i < pooled.rowCount; i++) expect(cellBw.isNone(i), false, `cell BW on row ${i}`);
    expect(pooled.col('BW_kg')!.isNone(0), true, 'per-animal BW keeps its honest null');
    // Phase-0 measurement: every cell of the gold fixture is weight-based at recorded precision.
    expect([...warningsByCell.values()].flat().some((w) => w.code === 'COMPOSITE_NOT_WEIGHT_BASED'), false,
      '06 is weight-based in every cell');

    const contract = nonmemNormalizedContract(false, true); // (withBw, withDose): the dose; BW = the cell's
    const withBw: ColumnContract = {observations: {...contract.observations,
      bodyWeight: {name: CELL_BW_COL, required: false, rawType: 'float64', units: 'kg'}}};
    const {results, outcome} = runPipelineOrThrow(pooled, withBw, nonmemScalars(), CARD_RULES);
    appendGroupingColumns(results, pooled, outcome.perProfile.map((p) => p.index), withBw);

    // AC-8b.2: per-kg CL vs the card's per-kg cl.obs within 1 % (NFR-05, fixed up front — R2 N-8).
    const card = await loadFixture('rat_tk28d_nca_per_cell.csv');
    // Only card cells with a usable fit: PKNCA also reports a cl.obs on a NEGATIVE adj-R² window
    // (D100_F_Day1: adj.R² −0.43, R² 0.05, %extrap 90 %), which the engine refuses even with no R²
    // gate — the "present ≠ valid" convention of METHOD-15, not a parity gap (measured at build).
    let checked = 0;
    for (let r = 0; r < card.rowCount; r++) {
      const clObs = card.col('cl.obs')!;
      if (clObs.isNone(r) || !((card.col('adj.r.squared')!.get(r) as number) >= 0)) continue;
      const cell = String(card.col('PseudoSubject')!.get(r));
      const want = (clObs.get(r) as number) * 1000 / 60;
      const row = cellRow(results, cell);
      expect(row >= 0, true, `card cell '${cell}' found in the results`);
      expectFloat(results.col('CL')!.get(row) as number, want, want * 0.01,
        `${cell}: per-kg CL vs card (1%) — CL abs. ${results.col('CL_ABS')!.get(row)}, ` +
        `BW_KG ${results.col('BW_KG')!.get(row)}, STATUS ${results.col('STATUS')!.get(row)}`);
      checked++;
    }
    expect(checked >= 1, true, `card cells with a usable cl.obs checked: ${checked} (D10_F_Day1)`);

    // The margin on the NONMEM path, against the value the card's per-kg CL gives under the same
    // anchor — the cohort's own mean body weight, read here from the raw file (F/10 mg/kg/Day 1).
    const ref = df.clone(DG.BitSet.create(df.rowCount, (i) => String(df.col('EVID')!.get(i)) === '0' &&
      Number(df.col('DOSE_MGKG')!.get(i)) === 10 && String(df.col('Sex')!.get(i)) === 'F' &&
      Number(df.col('DAY')!.get(i)) === 1));
    let bwSum = 0;
    for (let i = 0; i < ref.rowCount; i++) bwSum += ref.col('BW_kg')!.get(i) as number;
    const cohortBw = bwSum / ref.rowCount;
    const clKg = card.col('cl.obs')!.get(0) as number; // D10_F_Day1
    const noaelAuc = results.col('AUC_LAST')!.get(cellRow(results, 'D100_F_Day28')) as number;
    const expected = noaelAuc / ((200 / 70) / (clKg * Math.pow(70 / cohortBw, -0.25)));
    const margin = computeExposureMargin(results, {
      noaelMgKg: 100, humanDoseMg: 200, humanBodyWeightKg: 70,
      doseCol: 'DOSE_MGKG', studyDayCol: STUDY_DAY_COL, stratifyCols: ['Sex'],
    });
    const t = margin.table!;
    const f = [...Array(t.rowCount).keys()].find((i) => t.col('Sex')!.get(i) === 'F')!;
    const got = t.col('MARGIN_PROJECTED')!.get(f) as number;
    expectFloat(got, expected, expected * 0.01, `F MARGIN_PROJECTED ${got} vs ${expected} (1%)`);
    expect(got > 9, true, `the per-kg margin, not the HEAD 2.43 (got ${got})`);
    expect(String(t.col('RAT_REF_BW_SOURCE')!.get(f)).startsWith('cohort mean'), true, 'cohort anchor');
    return `NONMEM 06 F MARGIN_PROJECTED = ${got.toFixed(4)} (expected ${expected.toFixed(4)}; ` +
      `cohort BW ${cohortBw.toFixed(5)} kg; HEAD measured 2.4270)`;
  }, {timeout: TIMEOUT});

  test('AC-E6 (Rule 18): composite TK accumulation Rac matches oracle Day-28/Day-1 ratios', async () => {
    // The live pane surfaces an Accumulation tab on the composite frame. Rule 18:
    // a surfaced computation needs a reference assertion. The per-cell card has no
    // accumulation column, but Rac is derivable from it — Rac(AUC) =
    // auclast_Day28/auclast_Day1, Rac(Cmax) = cmax_Day28/cmax_Day1 per dose×sex.
    // This validates the grouping (StudyDay × nominal DOSE_MGKG × Sex) end-to-end
    // through the PRODUCTION computeAccumulationRatios, not just the division.
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec06(true));
    const {pooled} = poolCompositeByCell(norm.observations, {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex'});
    const contract = nonmemNormalizedContract(false, false);
    const {results, outcome} = runPipelineOrThrow(pooled, contract, nonmemScalars(), DEFAULT_RULES);
    appendGroupingColumns(results, pooled, outcome.perProfile.map((p) => p.index), contract);

    // tlastCol defaults to the absent 'TLAST' → guard unverified (not a hard NaN);
    // λz is NaN on the sparse composite → Formula A (AUClast/AUClast), the live path.
    const acc = computeAccumulationRatios(results, {
      tau: 24, doseCol: 'DOSE_MGKG', studyDayCol: STUDY_DAY_COL, stratifyCols: ['Sex'],
    });
    expect(acc.day1, 1, 'baseline day = 1');
    expect(acc.dayN, 28, 'steady-state day = 28');

    const card = (cell: string) => CARD_CELLS.find((c) => c.cell === cell)!;
    const accRow = (sex: string, dose: number): number => {
      const sx = acc.table.col('Sex')!; const ds = acc.table.col('DOSE')!;
      for (let i = 0; i < acc.table.rowCount; i++)
        if (String(sx.get(i)) === sex && Math.round(Number(ds.get(i))) === dose) return i;
      return -1;
    };
    // Both endpoints are pipeline auclast/cmax (each within 0.1% of the card), so
    // the ratio is within ~0.2% of the oracle ratio — assert at 0.5%.
    for (const [sex, dose] of [['F', 10], ['M', 300]] as ReadonlyArray<[string, number]>) {
      const row = accRow(sex, dose);
      expect(row >= 0, true, `accumulation row present for ${sex}/${dose}`);
      const d1 = card(`D${dose}_${sex}_Day1`); const dN = card(`D${dose}_${sex}_Day28`);
      const expRacAuc = dN.auclast / d1.auclast;
      const expRacCmax = dN.cmax / d1.cmax;
      expectFloat(acc.table.col('RAC_AUC')!.get(row) as number, expRacAuc, expRacAuc * 0.005,
        `${sex}/${dose} Rac(AUC)`);
      expectFloat(acc.table.col('RAC_CMAX')!.get(row) as number, expRacCmax, expRacCmax * 0.005,
        `${sex}/${dose} Rac(Cmax)`);
    }
  }, {timeout: TIMEOUT});
});

/** Locate a results row by its SUBJECT key (−1 if absent). */
function rowForSubject(results: DG.DataFrame, subj: string): number {
  for (let i = 0; i < results.rowCount; i++)
    if (String(results.col('SUBJECT')!.get(i)) === subj) return i;
  return -1;
}

category('Reference suite: NONMEM IV-infusion (RATE>0, fixtures 1 + 4b)', () => {
  // No ROUTE / ADDL / II columns → RATE>0 inference path. Conc is mg/L (the 1-cmt
  // analytic model's natural unit; dose 100 mg / AUC 50 → CL = 2), so the
  // contract pins concUnit = 'mg/L'. Reference JSONs are PKNCA 0.12.1, generated
  // by iv_infusion_nonmem_pknca_reference.R from the SAME analytic profile.
  const ecInf: EventContract = {
    id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID',
    cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
  };

  test('AC-1.1/1.2/1.3/1.4: fixture 1 RATE>0 → PKNCA infusion parity (MRT/Vss via per-subject duration)', async () => {
    const df = await loadFixture('04_iv_infusion_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ecInf);
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);

    // AC-1.1: route inferred from RATE>0 → IV-infusion; exactly one
    // NONMEM_INFUSION_RATE_INFERRED warning (severity info, BD-1) + the
    // route-confirm gate fires (routeConfirmRequired).
    expect(norm.routeConfirmRequired, true, 'RATE>0 inferred route → confirm gate fires');
    const w1 = norm.warningsBySubject.get('1') ?? [];
    const inf = w1.filter((w) => w.code === NONMEM_INFUSION_RATE_INFERRED);
    expect(inf.length, 1, 'exactly one NONMEM_INFUSION_RATE_INFERRED warning');
    expect(inf[0].severity, 'info', 'infusion-rate-inferred is informational (BD-1)');
    expect(norm.observations.col(DERIVED_ROUTE_CODE_COL)!.get(0) as number, nca.ROUTE_IV_INFUSION,
      'route resolves to IV-infusion');
    // AC-1.4: the production per-subject duration is threaded via the emitted
    // ~InfusionDuration column (= AMT/RATE = 100/100 = 1 h), NOT hand-set.
    const durCol = norm.observations.col(DERIVED_INFUSION_DURATION_COL)!;
    expectFloat(durCol.get(0) as number, 1, 1e-9, 'fixture 1 threaded duration = AMT/RATE = 1 h');

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    expect(results.rowCount, 1, 'one subject');

    const ref = await loadReference('04_iv_infusion_nonmem-reference.json');
    const prof = ref.profiles[0];
    // AC-1.2: cmax 0.1% / tmax exact / auclast,aucinf 0.1% / cl 1% / λz 0.5% vs PKNCA.
    ratRowAssert(results, 0, '1', prof, 'fixture1');
    // AC-1.3: MRT within 0.5%, Vss within 1% vs PKNCA mrt.iv.obs / vss.iv.obs.
    // A dropped duration would give the BOLUS MRT (≈ 5 + T_inf/2 = 5.5, ~10% high)
    // — excluded by the 0.5% band around 5.0028.
    expectFloat(results.col('MRT')!.get(0) as number, prof.parameters.mrt!,
      Math.abs(prof.parameters.mrt!) * 0.005, 'fixture1: MRT (AC-1.3)');
    expectFloat(results.col('VSS_ABS')!.get(0) as number, prof.parameters.vss!,
      Math.abs(prof.parameters.vss!) * 0.01, 'fixture1: Vss abs. (AC-1.3)');
    expect(Number.isFinite(results.col('MRT')!.get(0) as number), true,
      'fixture1: MRT finite (AC-1.4 — duration threaded through production path)');
  }, {timeout: TIMEOUT});

  test('AC-A2.1/2.2/2.3: fixture 4b two-subject infusion → per-subject duration (single pipeline)', async () => {
    // R2-C-1 (anti-vacuity): BOTH subjects processed in ONE pipeline invocation
    // on the multi-subject frame, no hand-set ProfileInputs.infusionDuration.
    // Subject 1 T_inf=1 h, subject 2 T_inf=2 h. A global-scalar impl applies one
    // subject's duration to both → the mismatched subject's −T_inf/2 correction
    // is wrong (off ≈0.5 h ≈10% in MRT) → fails the 0.5% MRT band. Both recover
    // MRT≈5 / Vss≈10 ONLY with per-subject threading (AC-A2.3 discriminator).
    const df = await loadFixture('04b_iv_infusion_multi_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ecInf);
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);

    // AC-A2.1: each RATE>0 subject carries its OWN AMT/RATE-derived duration.
    const durCol = norm.observations.col(DERIVED_INFUSION_DURATION_COL)!;
    const subjCol = norm.observations.col(SUBJECT_COL)!;
    const durBySubj = new Map<string, number>();
    for (let i = 0; i < norm.observations.rowCount; i++)
      durBySubj.set(String(subjCol.get(i)), durCol.get(i) as number);
    expectFloat(durBySubj.get('1')!, 1, 1e-9, 'subject 1 duration = 100/100 = 1 h');
    expectFloat(durBySubj.get('2')!, 2, 1e-9, 'subject 2 duration = 100/50 = 2 h');

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    expect(results.rowCount, 2, 'two subjects in one pipeline run');

    const ref = await loadReference('04b_iv_infusion_multi_nonmem-reference.json');
    const refBySubj = indexFixtureProfiles(ref);
    for (const subj of ['1', '2']) {
      const row = rowForSubject(results, subj);
      expect(row >= 0, true, `subject ${subj} present`);
      const prof = refBySubj.get(subj)!;
      ratRowAssert(results, row, subj, prof, `fixture4b subj${subj}`);
      // The load-bearing discriminator: each subject recovers MRT≈5 / Vss≈10
      // (within 0.5% / 1%) — only possible if its OWN duration was applied.
      expectFloat(results.col('MRT')!.get(row) as number, prof.parameters.mrt!,
        Math.abs(prof.parameters.mrt!) * 0.005, `fixture4b subj${subj}: MRT (AC-A2.3)`);
      expectFloat(results.col('VSS_ABS')!.get(row) as number, prof.parameters.vss!,
        Math.abs(prof.parameters.vss!) * 0.01, `fixture4b subj${subj}: Vss abs. (AC-A2.3)`);
    }
  }, {timeout: TIMEOUT});

  test('AC-A2.2: RATE=0 NONMEM (03) emits all-NaN ~InfusionDuration → scalar fallback (byte-unchanged)', async () => {
    // Science-preservation gate: the per-subject duration column is the ONLY new
    // input the BD-2 fix threads. For RATE=0 subjects resolveRoute yields
    // infusionDuration=null → NaN in ~InfusionDuration → extractor falls back to
    // scalars.infusionDuration (null). So every existing RATE=0 fixture (03/06)
    // feeds computeNca the identical inputs as before the fix — AC-E1..E6 stay
    // green. This asserts the mechanism directly: the 03 frame's duration column
    // is entirely NaN, so no RATE=0 profile can receive a spurious duration.
    const df = await loadFixture('03_rat_multiroute_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, {
      id: 'SubjectID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID',
      cmt: 'CMT', rate: 'RATE', mdv: 'MDV', route: 'ROUTE', bodyWeight: 'BW_kg',
      covariates: ['Sex', 'DOSE_MGKG'],
    });
    const durCol = norm.observations.col(DERIVED_INFUSION_DURATION_COL)!;
    let allNaN = true;
    for (let i = 0; i < durCol.length; i++)
      if (!Number.isNaN(durCol.get(i) as number)) {allNaN = false; break;}
    expect(allNaN, true, '03 (RATE=0) ~InfusionDuration is all-NaN → scalar fallback → byte-unchanged');
  }, {timeout: TIMEOUT});
});

/** Minimum TAD among a profile's rows in the normalized frame (−1 if absent). */
function minTadForSubject(norm: DG.DataFrame, subj: string): number {
  const s = norm.col(SUBJECT_COL)!; const t = norm.col(TAD_COL)!;
  let min = Number.POSITIVE_INFINITY;
  for (let i = 0; i < norm.rowCount; i++)
    if (String(s.get(i)) === subj) min = Math.min(min, t.get(i) as number);
  return Number.isFinite(min) ? min : -1;
}

category('Reference suite: NONMEM EVID reset → new episode (fixtures 3a + 3b)', () => {
  // CMT=2 depot (no RATE) → extravascular; no ROUTE column → route inferred.
  // EVID=3 (3a) and EVID=4 (3b) take distinct code paths; each yields ONE reset
  // and splits the record into two independent episodes (`S1_ep0`, `S1_ep1`).
  const ecReset: EventContract = {
    id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', mdv: 'MDV',
  };

  async function runReset(csv: string, jsonRef: string, tag: string): Promise<void> {
    const df = await loadFixture(csv);
    const norm = normalizeEventDrivenObservations(df, ecReset);
    expect(norm.errors.length, 0, `${tag}: no blocking errors; got ${JSON.stringify(norm.errors)}`);
    // AC-3.1: exactly one reset, EVID_RESET_NEW_EPISODE warning.
    expect(norm.classificationCounts.reset, 1, `${tag}: classificationCounts.reset == 1`);
    const w = norm.warningsBySubject.get('S1') ?? [];
    expect(w.some((x) => x.code === EVID_RESET_NEW_EPISODE), true,
      `${tag}: EVID_RESET_NEW_EPISODE warning`);

    // AC-3.3 (discriminator): episode-1 first observation TAD == 0 (anchor reset);
    // a no-anchor-reset impl keeps global-origin TAD and would only pass AUC by
    // coincidence. (Never expect(x, undefined) — exact 0 check.)
    expect(minTadForSubject(norm.observations, 'S1_ep1'), 0, `${tag}: episode-1 first obs TAD == 0`);
    expect(minTadForSubject(norm.observations, 'S1_ep0'), 0, `${tag}: episode-0 first obs TAD == 0`);

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    // AC-3.2: two pseudo-subjects, parity each vs its per-episode PKNCA reference.
    expect(results.rowCount, 2, `${tag}: two episode pseudo-subjects`);
    const ref = await loadReference(jsonRef);
    const byKey = indexFixtureProfiles(ref);
    for (const ep of ['S1_ep0', 'S1_ep1']) {
      const row = rowForSubject(results, ep);
      expect(row >= 0, true, `${tag}: ${ep} present`);
      ratRowAssert(results, row, ep, byKey.get(ep)!, `${tag} ${ep}`);
    }
    // AC-3.4: episode-1 AUC is NOT polluted by pre-reset exposure — the two
    // episodes are the same profile, so auclast must match within tolerance.
    const r0 = rowForSubject(results, 'S1_ep0'); const r1 = rowForSubject(results, 'S1_ep1');
    const a0 = results.col('AUC_LAST')!.get(r0) as number;
    const a1 = results.col('AUC_LAST')!.get(r1) as number;
    expectFloat(a1, a0, Math.abs(a0) * 0.001, `${tag}: episode-1 AUC un-polluted (== episode-0)`);
  }

  test('AC-3.1..3.4 (3a): EVID=3 pure reset → two episodes, reset count 1, TAD reset', async () => {
    await runReset('07a_evid3_reset_nonmem.csv', '07a_evid3_reset_nonmem-reference.json', '3a');
  }, {timeout: TIMEOUT});

  test('AC-3.1..3.4 (3b): EVID=4 reset+dose → two episodes, reset count 1, TAD reset', async () => {
    await runReset('07b_evid4_reset_nonmem.csv', '07b_evid4_reset_nonmem-reference.json', '3b');
  }, {timeout: TIMEOUT});
});

category('Reference suite: NONMEM non-zero / multi-dose ladders (fixtures 6 + 7)', () => {
  test('AC-6.1..6.3: fixture 6 washout-delayed first dose (TAD-rebased, first interval only)', async () => {
    // Dose @ TIME=72, ADDL=2, II=24 → ladder {72,96,120}; obs only in [72,96].
    const ec: EventContract = {
      id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', mdv: 'MDV',
      addl: 'ADDL', ii: 'II',
    };
    const df = await loadFixture('09_washout_dose_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec);
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);
    // AC-6.1: union ladder built off the non-zero first dose time.
    const ladder = norm.doseLadderBySubject.get('S1') ?? [];
    expect(JSON.stringify(Array.from(ladder)), JSON.stringify([72, 96, 120]),
      'doseLadderBySubject(S1) = {72,96,120}');
    // AC-6.3: a ladder hard-wired to start at 0 would push TAD=72.5 outside [0,II]
    // → dropped → empty profile. The correct ladder rebases the first obs to
    // TAD=0.5 (non-empty profile, the assertDoseLadder ladder[0]===tDose guard).
    expect(minTadForSubject(norm.observations, 'S1_Day1'), 0.5,
      'first interval first obs rebased to TAD=0.5 (non-empty profile)');

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    // R2-N-1: oracle scoped to the first interval only (intervals 2/3 have no obs
    // → exactly one profile; a hard-wired-0 ladder dropping TAD=72.5 would give 0).
    expect(results.rowCount, 1, 'only the first interval has observations → one profile');
    const row = rowForSubject(results, 'S1_Day1');
    expect(row >= 0, true, 'S1_Day1 present');
    const ref = await loadReference('09_washout_dose_nonmem-reference.json');
    // AC-6.2: parity vs the TAD-rebased PKNCA reference (first interval).
    ratRowAssert(results, row, 'S1_Day1', ref.profiles[0], 'fixture6 S1_Day1');
  }, {timeout: TIMEOUT});

  test('AC-7.1..7.3: fixture 7 loading+maintenance union ladder, per-interval dose-of-record', async () => {
    // Loading 200 @ t=0 (ADDL=0) + maintenance 100 @ t=24 (ADDL=2, II=24).
    const ec: EventContract = {
      id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', mdv: 'MDV',
      addl: 'ADDL', ii: 'II',
    };
    const df = await loadFixture('10_loading_maint_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ec);
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);
    // AC-7.1 / AC-7.3: union ladder {0,24,48,72}, no phantom blow-up (FM-2 dedup).
    const ladder = norm.doseLadderBySubject.get('S1') ?? [];
    expect(JSON.stringify(Array.from(ladder)), JSON.stringify([0, 24, 48, 72]),
      'doseLadderBySubject(S1) = union {0,24,48,72} (no phantom doses)');

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    // AC-7.3: exactly two interval profiles — a phantom-dose blow-up would add a
    // spurious S1_Day3+; assert the count, not just that the expected two exist.
    expect(results.rowCount, 2, 'exactly two interval profiles (no phantom doses)');
    const ref = await loadReference('10_loading_maint_nonmem-reference.json');
    const byKey = indexFixtureProfiles(ref);
    // AC-7.2: per-interval AUCτ + CL against the per-interval dose-of-record
    // (loading=200 for S1_Day1, maintenance=100 for S1_Day2). A single-AMT
    // oracle would mis-state the loading CL; the two-PKNCAdose reference does not.
    for (const ep of ['S1_Day1', 'S1_Day2']) {
      const row = rowForSubject(results, ep);
      expect(row >= 0, true, `${ep} present`);
      ratRowAssert(results, row, ep, byKey.get(ep)!, `fixture7 ${ep}`);
    }
  }, {timeout: TIMEOUT});
});

/** Max TAD among a profile's rows in the normalized frame (−1 if absent). */
function maxTadForSubject(norm: DG.DataFrame, subj: string): number {
  const s = norm.col(SUBJECT_COL)!; const t = norm.col(TAD_COL)!;
  let max = Number.NEGATIVE_INFINITY;
  for (let i = 0; i < norm.rowCount; i++)
    if (String(s.get(i)) === subj) max = Math.max(max, t.get(i) as number);
  return Number.isFinite(max) ? max : -1;
}

category('Reference suite: NONMEM SS=1 steady state (fixture 2)', () => {
  // SS=1 + II=24, one τ interval [0,24] of a 1-cmt IV-bolus at steady state
  // (k=0.2). Ground truth AUCτ,ss = Dose/CL = 50 (derived in the R script);
  // PKNCA integrates the SAME analytic C_ss(t) grid and lands on 50 too
  // (consistency 0%). The engine does NOT clip [0,II] under SS — the one-interval
  // construction is enforced structurally (AC-2.3). Reject branches (SS=1 no II;
  // SS=2) are already covered by AC-D3 in adapter-event-driven-tests.ts (rule 5,
  // no duplication) — fixture 2 adds the oracle-grounded AUCτ parity.
  const TAU = 24; const AUC_HAND = 50; // Dose/CL = 100/2
  const ecSs: EventContract = {
    id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', mdv: 'MDV',
    ss: 'SS', ii: 'II',
  };

  test('AC-2.1/2.2/2.3: SS=1 honored → AUCτ within 0.5% of Dose/CL (+ PKNCA co-primary)', async () => {
    const df = await loadFixture('05_ss_singledose_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ecSs);
    // AC-2.1: exactly the SS-assumed warning, no blocking error.
    expect(norm.errors.length, 0, `SS=1+II honored, no error; got ${JSON.stringify(norm.errors)}`);
    const w = norm.warningsBySubject.get('S1') ?? [];
    expect(w.filter((x) => x.code === NONMEM_STEADY_STATE_ASSUMED).length, 1,
      'exactly one NONMEM_STEADY_STATE_ASSUMED warning');

    // AC-2.3 (structural self-check): the engine does not enforce the [0,II]
    // window under SS, so the one-interval construction must be machine-verified.
    expect(maxTadForSubject(norm.observations, 'S1') <= TAU, true,
      'SS fixture obs constrained to one interval (lastObsTAD ≤ II)');

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});
    expect(results.rowCount, 1, 'one SS subject');

    const auc = results.col('AUC_LAST')!.get(0) as number;
    // AC-2.2: AUCτ within 0.5% of the hand-derived Dose/CL ground truth.
    expectFloat(auc, AUC_HAND, AUC_HAND * 0.005, 'AUCτ within 0.5% of Dose/CL = 50');
    // Co-primary: within 0.1% of the PKNCA reference (same analytic grid + method).
    const ref = await loadReference('05_ss_singledose_nonmem-reference.json');
    const aucRef = ref.profiles[0].parameters.auclast!;
    expectFloat(auc, aucRef, Math.abs(aucRef) * 0.001, 'AUCτ within 0.1% of PKNCA AUClast');
  }, {timeout: TIMEOUT});
});

category('Reference suite: NONMEM per-day-clock reset (fixture 4)', () => {
  // 08_perday_reset_nonmem.csv PROVENANCE: hand-sliced from the D10_F composite
  // means of dataset 06 (the Day-1 and Day-28 rows of `06_tk_serial_demo.csv`,
  // mg/L), re-encoded with the per-DAY clock and two explicit per-day dose rows
  // at TIME=0 (no ADDL). No generator R script — the oracle is dataset 06's
  // already-rule-16-verified per-cell card (CARD_CELLS / rat_tk28d_nca_per_cell),
  // so this fixture asserts cross-encoding equivalence, not a fresh PKNCA run.
  // per-DAY clock: TIME is the in-day clock (0…24 for BOTH Day 1 and Day 28),
  // DAY distinguishes occasions, and each day carries an EXPLICIT dose row at
  // TIME=0 (NO ADDL) — contrast 06's absolute clock + ADDL=27. The conc values
  // are 06's D10_F composite means (mg/L), so the per-cell composites must match
  // 06's card (cross-encoding equivalence, AC-E4 spirit; rule-16-verified
  // Day-28 auclast = 219.451089543992). BD-3: the per-day path resolves occasion
  // from the DAY column with TAD=TIME (dose at 0) — this test IS the rule-16
  // check that the per-day branch produces the right composites; it must be run
  // green on `grok test`.
  const ecDay: EventContract = {
    id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', mdv: 'MDV',
    day: 'DAY',
  };
  // D10_F per-cell card values (identical to CARD_CELLS / dataset 06).
  const D10F = {
    Day1: {auclast: 40.4501897648454, cmax: 2.2118, tmax: 2},
    Day28: {auclast: 219.451089543992, cmax: 10.21428, tmax: 8},
  } as const;

  test('AC-4.1/4.2: per-day clock → same D10_F composites as 06; no ADDL expansion', async () => {
    const df = await loadFixture('08_perday_reset_nonmem.csv');
    const norm = normalizeEventDrivenObservations(df, ecDay);
    expect(norm.errors.length, 0, `no blocking errors; got ${JSON.stringify(norm.errors)}`);

    // AC-4.2 (coverage discriminator): two explicit dose rows at the in-day TIME=0
    // dedup to a single anchor {0} — NO 28-element ADDL ladder (06's path). Each
    // dose's expandDoseLadder(0, NaN, NaN) = [0].
    const ladder = norm.doseLadderBySubject.get('D10_F') ?? [];
    expect(JSON.stringify(Array.from(ladder)), JSON.stringify([0]),
      'per-day explicit doses → ladder {0} (no ADDL expansion, unlike 06)');

    const contract = nonmemNormalizedContract(false, true, 'mg/L');
    const {results} = runPipelineOrThrow(
      norm.observations, contract, nonmemScalars('mg/L'), DEFAULT_RULES,
      {routeWarningsBySubject: norm.warningsBySubject});

    // AC-4.1: Day-1 and Day-28 occasions reproduce the 06 per-cell card.
    expect(results.rowCount, 2, 'two DAY occasions → exactly two profiles');
    for (const [day, exp] of [['Day1', D10F.Day1], ['Day28', D10F.Day28]] as const) {
      const row = rowForSubject(results, `D10_F_${day}`);
      expect(row >= 0, true, `D10_F_${day} occasion present`);
      expectFloat(results.col('AUC_LAST')!.get(row) as number, exp.auclast,
        Math.abs(exp.auclast) * 0.001, `D10_F_${day} auclast == 06 card`);
      expectFloat(results.col('CMAX')!.get(row) as number, exp.cmax,
        Math.abs(exp.cmax) * 0.001, `D10_F_${day} cmax == 06 card`);
      expect(results.col('TMAX')!.get(row) as number, exp.tmax, `D10_F_${day} tmax == 06 card`);
    }
  }, {timeout: TIMEOUT});
});
