/**
 * Conc-mis-scale guard — integration fire/no-fire vs real `files/*.csv`
 * (nca-conc-dose-misscale-guard AD-1 / AD-2; rule 16 — assert observable
 * fire/no-fire against actual generated magnitudes, never a mirror pattern).
 *
 * **The mis-scale exemplar.** `03_rat_iv.csv` carries native ng/mL in
 * `Conc_ng_mL`. Mapping it as `mg/L` (the canonical default — a user who left
 * the picker alone) reproduces the silent ng/mL→mg/L 1000× error exactly: the
 * raw values pass through `scaleProfileInputs` un-rescaled. This is DATA-GAP-15's
 * stand-in (a labelled synthetic fixture is a future nice-to-have).
 *
 * **P1 (raw-conc ceiling) reach vs P2 (implied-CL) reach — the headline.**
 *  - P1 fires on the egregious mid/high-dose mis-scales (5/10 mg/kg, median
 *    895–3126 mg/L > 500 ceiling) but STRUCTURALLY MISSES the 1 mg/kg low-dose
 *    subjects (median ≤388.6 < 500) — AC-3, an *expected* non-fire.
 *  - P2 fires on ALL 24 mis-scaled subjects INCLUDING the 1 mg/kg subjects P1
 *    misses (implied CL/kg ≈ 3.1e-4 < 1e-3 floor) — AC-7, the proof P2 reaches
 *    further. Neither replaces the unit picker.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';
import * as DG from 'datagrok-api/dg';

import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, withDoseMg,
  theophContract, theophScalars, indomethContract, indomethScalars,
  ratContract, ratScalars, runPipelineOrThrow,
} from './test-utils';
import type {ColumnContract} from '../adapter/column-contract';
import type {ProfileScalars} from '../adapter/profile-inputs';
import type {PerProfileResult} from '../controller/nca-controller';

const TIMEOUT = 30000;
const P1 = 'CONC_MAGNITUDE_HEURISTIC_WARNING';
const P2 = 'CL_IMPLAUSIBLE_HEURISTIC_WARNING';

/** Rat contract that MIS-LABELS the native ng/mL `Conc_ng_mL` column as the
 *  canonical `mg/L` — the silent 1000× error this guard exists to catch. */
function ratMisscaleContract(df: DG.DataFrame): ColumnContract {
  const base = ratContract(df);
  return {
    observations: {
      ...base.observations,
      concentration: {...base.observations.concentration!, units: 'mg/L'},
    },
  };
}

/** Rat scalars with the concentration unit left at the canonical-default mg/L
 *  (i.e. the user did NOT pick ng/mL) — drives the un-rescaled mis-scale. */
function ratMisscaleScalars(route: nca.RouteCode): ProfileScalars {
  return {...ratScalars(route), concentrationUnits: 'mg/L'};
}

/** Map SubjectID → Dose_mg_kg from the raw rat fixture. */
function subjectDoseMap(df: DG.DataFrame): Map<string, number> {
  const out = new Map<string, number>();
  const s = df.col('SubjectID')!;
  const d = df.col('Dose_mg_kg')!;
  for (let i = 0; i < df.rowCount; i++) out.set(String(s.get(i)), Number(d.get(i)));
  return out;
}

function hasCode(pp: PerProfileResult, code: string): boolean {
  return pp.result.provenance.warnings.some((w) => w.code === code);
}

category('Conc-mis-scale guard: P1 raw-conc ceiling (integration)', () => {
  test('AC-2/AC-3 — mis-scaled rat IV: 5/10 mg/kg fire P1, 1 mg/kg is the documented miss', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    const dose = subjectDoseMap(df);
    const {outcome} = runPipelineOrThrow(
      df, ratMisscaleContract(df), ratMisscaleScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);

    let lowChecked = 0;
    let midHighChecked = 0;
    for (const pp of outcome.perProfile) {
      const d = dose.get(pp.index.subjectKey)!;
      if (d === 1) {
        // AC-3: documented expected NON-fire (median ≤388.6 < 500 ceiling).
        expect(hasCode(pp, P1), false,
          `1 mg/kg subject ${pp.index.subjectKey} is P1's documented low-dose miss (must NOT fire)`);
        lowChecked++;
      } else {
        // AC-2: 5 / 10 mg/kg medians 895–3126 mg/L > 500 → fire.
        expect(hasCode(pp, P1), true,
          `${d} mg/kg subject ${pp.index.subjectKey} (median > 500) must fire P1`);
        midHighChecked++;
      }
    }
    expect(lowChecked, 8, 'eight 1 mg/kg subjects asserted as documented misses');
    expect(midHighChecked, 16, 'sixteen 5/10 mg/kg subjects asserted as fires');
  }, {timeout: TIMEOUT});

  test('AC-1 — legit theoph / indometh do NOT fire P1 (genuinely mg/L, below ceiling)', async () => {
    const theoph = await loadFixture('theoph.csv');
    const {outcome: tOut} = runPipelineOrThrow(
      withDoseMg(theoph), theophContract(theoph), theophScalars(), DEFAULT_RULES);
    for (const pp of tOut.perProfile)
      expect(hasCode(pp, P1), false, `theoph ${pp.index.subjectKey} (median ~6.9 mg/L) must NOT fire P1`);

    const indo = await loadFixture('indometh.csv');
    const {outcome: iOut} = runPipelineOrThrow(
      indo, indomethContract(indo), indomethScalars(), DEFAULT_RULES);
    for (const pp of iOut.perProfile)
      expect(hasCode(pp, P1), false, `indometh ${pp.index.subjectKey} (median ~0.4 mg/L) must NOT fire P1`);
  }, {timeout: TIMEOUT});

  test('AC-1 (REGRESSION_PIN) — legit high-dose 06_tk D300 (298 mg/L) does NOT fire P1', async () => {
    // The closest legitimate profile to the ceiling: 06_tk_serial 300 mg/kg
    // Day28 median 298.3 mg/L < 500. A false positive here would mean the
    // ceiling is mis-calibrated against real high-dose TK exposure (rule 16).
    const df = await loadFixture('06_tk_serial_demo.csv');
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time_hr', required: true, rawType: 'float64', units: 'h'},
        concentration: {name: 'Conc_mg_L', required: true, rawType: 'float64', units: 'mg/L'},
        dose: {name: 'Dose_mg_kg', required: false, rawType: 'float64', units: 'mg/kg'},
      },
    };
    const scalars: ProfileScalars = {
      dose: 0, route: nca.ROUTE_PO, doseUnits: 'mg/kg', concentrationUnits: 'mg/L',
      timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0,
    };
    const {outcome} = runPipelineOrThrow(df, contract, scalars, DEFAULT_RULES);
    for (const pp of outcome.perProfile) {
      expect(hasCode(pp, P1), false,
        `legit 06_tk ${pp.index.subjectKey} (max median 298 mg/L) must NOT fire P1`);
    }
  }, {timeout: TIMEOUT});

  test('AC-6 — P1 warning rides WARNINGS_JSON end-to-end (serialized in the results column)', async () => {
    // The full pipeline → WARNINGS_JSON path for P1 specifically: a mis-scaled
    // 5/10 mg/kg subject (median > 500) must serialize CONC_MAGNITUDE_HEURISTIC_WARNING
    // into the results column, exactly as AC-10 proves for P2.
    const df = await loadFixture('03_rat_iv.csv');
    const dose = subjectDoseMap(df);
    const {results, outcome} = runPipelineOrThrow(
      df, ratMisscaleContract(df), ratMisscaleScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);
    const warnCol = results.col('WARNINGS_JSON')!;
    const subjCol = results.col('SUBJECT')!;
    let serialized = 0;
    for (let i = 0; i < results.rowCount; i++) {
      if (dose.get(String(subjCol.get(i))) === 1) continue; // 1 mg/kg is P1's documented miss
      const warns = JSON.parse(String(warnCol.get(i))) as Array<{code: string}>;
      expect(warns.some((w) => w.code === P1), true,
        `mis-scaled row ${i} (5/10 mg/kg) must serialize P1 into WARNINGS_JSON`);
      serialized++;
    }
    expect(serialized, 16, 'all 16 mid/high-dose mis-scaled rows carry P1 in WARNINGS_JSON');
    expect(outcome.perProfile.length, 24);
  }, {timeout: TIMEOUT});
});

category('Conc-mis-scale guard: P2 implied-CL cross-check (integration)', () => {
  test('AC-7 — mis-scaled rat IV: P2 fires on ALL 24 subjects, incl. the 1 mg/kg P1 misses', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    const dose = subjectDoseMap(df);
    const {outcome} = runPipelineOrThrow(
      df, ratMisscaleContract(df), ratMisscaleScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);

    expect(outcome.perProfile.length, 24, 'all 24 rat IV subjects present');
    let lowDoseProof = 0;
    for (const pp of outcome.perProfile) {
      expect(hasCode(pp, P2), true,
        `mis-scaled subject ${pp.index.subjectKey} must fire P2 (implied CL/kg < 1e-3)`);
      // The reaches-further proof: the 1 mg/kg subjects P1 does NOT fire, P2 DOES.
      if (dose.get(pp.index.subjectKey) === 1) {
        expect(hasCode(pp, P1), false, `${pp.index.subjectKey} (1 mg/kg) P1 miss precondition`);
        expect(hasCode(pp, P2), true, `${pp.index.subjectKey} (1 mg/kg) P2 reaches it`);
        lowDoseProof++;
      }
    }
    expect(lowDoseProof, 8, 'eight 1 mg/kg subjects: P1 misses, P2 reaches (the headline)');
  }, {timeout: TIMEOUT});

  test('AC-8 — legit profiles do NOT fire P2 (CL/kg in band): rat IV (BW) + theoph (Wt)', async () => {
    // True band-pass no-fire: BW present, correctly-scaled ng/mL → CL/kg in band.
    const rat = await loadFixture('03_rat_iv.csv');
    const {outcome: rOut} = runPipelineOrThrow(
      rat, ratContract(rat), ratScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);
    for (const pp of rOut.perProfile) {
      expect(hasCode(pp, P2), false,
        `legit rat IV ${pp.index.subjectKey} (CL/kg ~0.31–0.99) must NOT fire P2`);
    }

    const theoph = await loadFixture('theoph.csv');
    const {outcome: tOut} = runPipelineOrThrow(
      withDoseMg(theoph), theophContract(theoph), theophScalars(), DEFAULT_RULES);
    for (const pp of tOut.perProfile) {
      expect(hasCode(pp, P2), false,
        `legit theoph ${pp.index.subjectKey} (CL/kg ~0.03–0.06) must NOT fire P2`);
    }
  }, {timeout: TIMEOUT});

  test('AC-9 — BW absent (indometh) → P2 returns null (no CL_IMPLAUSIBLE warning)', async () => {
    // indometh is Subject,time,conc only — no BW column. P2 cannot per-kg
    // normalize → null. Its no-fire is the BW-absent path, NOT a band-pass.
    const df = await loadFixture('indometh.csv');
    const {outcome} = runPipelineOrThrow(df, indomethContract(df), indomethScalars(), DEFAULT_RULES);
    for (const pp of outcome.perProfile)
      expect(hasCode(pp, P2), false, `indometh ${pp.index.subjectKey} (BW absent) → P2 null`);
  }, {timeout: TIMEOUT});

  test('AC-9 (SPEC CORRECTION) — legit 06_tk is per_kg/BW-absent → P2 does not fire', async () => {
    // RULE-16 CORRECTION: the synthesis AC-8 assumed 06_tk_serial "has BW" and
    // would be a band-pass no-fire. The actual CSV header carries only
    // `HumanBodyWeight_kg` (the human margin reference) — there is NO rat
    // subject body-weight column, so its dose (mg/kg, no BW) takes the per_kg path.
    // Since the S1 extension (nca-cl-bw-normalised-units §0b, ratified 2026-10-23) P2 DOES read a
    // per_kg engine CL. The no-fire stays, for a third reason the build-entry corpus check measured
    // (2026-10-23): none of the 16 τ-profiles has a finite CL (λz does not fit inside a 24 h dosing
    // interval at the default rules), and P2 never fabricates a warning from a missing CL. A per_kg
    // profile WITH a CL is checked against the band — asserted in cl-plausibility-tests (S1).
    const df = await loadFixture('06_tk_serial_demo.csv');
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time_hr', required: true, rawType: 'float64', units: 'h'},
        concentration: {name: 'Conc_mg_L', required: true, rawType: 'float64', units: 'mg/L'},
        dose: {name: 'Dose_mg_kg', required: false, rawType: 'float64', units: 'mg/kg'},
      },
    };
    const scalars: ProfileScalars = {
      dose: 0, route: nca.ROUTE_PO, doseUnits: 'mg/kg', concentrationUnits: 'mg/L',
      timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0,
    };
    const {outcome} = runPipelineOrThrow(df, contract, scalars, DEFAULT_RULES);
    for (const pp of outcome.perProfile) {
      expect(pp.doseNormalization, 'per_kg',
        `06_tk ${pp.index.subjectKey} has no rat BW → per_kg basis (the spec-correction premise)`);
      expect(hasCode(pp, P2), false, `06_tk ${pp.index.subjectKey} per_kg → no P2 warning`);
      // The measured reason (rule 16): no finite CL to judge — or, if a future engine fits one, it is
      // inside the band (the legit per-kg corpus spans 0.0187–3.25 L/h/kg).
      const cl = pp.result.values.cl;
      expect(!Number.isFinite(cl) || (cl >= 5e-3 && cl <= 15), true, `06_tk ${pp.index.subjectKey} CL ${cl}`);
    }
  }, {timeout: TIMEOUT});

  test('AC-10 — P2 warning rides WARNINGS_JSON; failed-profile all-NaN invariant intact', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    const {results} = runPipelineOrThrow(
      df, ratMisscaleContract(df), ratMisscaleScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);
    const warnCol = results.col('WARNINGS_JSON')!;
    let serialized = 0;
    for (let i = 0; i < results.rowCount; i++) {
      const warns = JSON.parse(String(warnCol.get(i))) as Array<{code: string}>;
      if (warns.some((w) => w.code === P2)) {
        serialized++;
        // Any P2-flagged row must still carry a real (non-NaN) CL — P2 only
        // INSPECTS the engine CL, never blanks values (science-preservation).
        const cl = results.col('CL')!.get(i) as number;
        const status = String(results.col('STATUS')!.get(i));
        if (status === 'ok')
          expect(Number.isFinite(cl), true, `row ${i}: ok status → CL finite (P2 never writes values)`);
      }
    }
    expect(serialized, 24, 'all 24 mis-scaled rows serialize the P2 warning into WARNINGS_JSON');
  }, {timeout: TIMEOUT});
});
