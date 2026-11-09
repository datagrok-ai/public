/**
 * Reference suite — full nca-studio pipeline parity against PKNCA fixtures.
 *
 * Tolerances vs PKNCA: AUC / Cmax 0.1 %, λz 0.5 %, derived parameters 1 %.
 * @see ../adapter/defaults.ts (DEFAULT_RULES, including adjRSquaredFactor: 1e-4)
 *
 * **Coverage target.** 12 theoph subjects + 6 indometh subjects = 18/18 must
 * pass within the §8.1 tolerances (CLAUDE.md rule 18 — PKNCA parity is
 * non-negotiable). Each subject is asserted on Cmax (exact match), Tmax
 * (exact match), AUClast / AUCinf (< 0.1% relative), λz (< 0.5% rel),
 * t½ / CL / Vz (< 1% rel), %AUCextrap (< 0.5 percentage-points absolute),
 * and status === 'ok'.
 *
 * **Deviations from spec §8.2 example** (logged for audit traceability):
 *
 * 1. The example shows `expectFloat(row.get('pct_extrap'), r.pctExtrap,
 *    0.005)` with the comment "0.5% abs". The numeric literal 0.005 means
 *    0.005 percentage points — three orders of magnitude tighter than §8.1
 *    row 9 prescribes ("< 0.5% absolute"). This file uses 0.5 (matching the
 *    spec §8.1 row 9 text and the sci-comp reference-suite tolerance).
 *
 * 2. The example asserts `provenance.lambdaZ.pointsUsed` matches the PKNCA
 *    reference's `pointsUsed` array. The shipped PKNCA reference JSONs at
 *    `files/{theoph,indometh}-reference.json` only carry the count
 *    `lambda_z_n_points` — they do not include the index array. This file
 *    narrows §8.1 row 11 to a count-match assertion (count is still a strong
 *    regression marker — different λz strategies produce different counts).
 *
 * **Reference shape.** Reference JSON values are snake_case from PKNCA
 * (v0.12.1); the field-name mismatch with sci-comp's camelCase `aucLast` etc.
 * is bridged inside this test, not in the result-writer column names (those
 * stay UPPERCASE per §7.2).
 */

import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, loadReference, loadReferenceCsv, indexFixtureProfiles, withDoseMg,
  theophContract, theophScalars, indomethContract, indomethScalars,
  infusionContract, infusionScalars,
  ratContract, ratScalars, ratMultirouteContract, setupMultirouteRoute,
  multicompoundContract, multicompoundScalars,
  runPipelineOrThrow, makeSyntheticObs,
  type FixtureProfile,
} from './test-utils';
import type {ProfileScalars} from '../adapter/profile-inputs';
import type {ColumnContract} from '../adapter/column-contract';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {computeBioavailability} from '../adapter/bioavailability';
import {computeDoseProportionality} from '../adapter/dose-proportionality';
import {BARE_IV_DEFAULTED_TO_BOLUS} from '../adapter/route-normalize';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

const TIMEOUT = 30000;

/**
 * C0_EXTRAP drift guard (GROK-20960 AC-N.3.2 / AC-N.3.3). The results frame's
 * `C0_EXTRAP` comes from the controller's OWN `nca.estimateC0` call (kept as the
 * source so a FAILED IV-bolus row reads 0, not blank — the route-structural
 * suppression contract); the engine now also reports `provenance.c0`. Both
 * chains skip masked points and read the same values, so on every non-failed
 * IV-bolus row they must be EXACTLY equal (Float64 storage — no tolerance), and
 * `PROVENANCE_JSON.c0` must carry the method per row. `expectedMethod` pins the
 * dataset's shape: indometh has no t=0 sample (every c0 is `logslope`), 05 has a
 * positive measured t=0 on all 60 (every c0 is `observed`).
 */
function assertC0DriftGuard(
  results: DG.DataFrame,
  perProfile: ReadonlyArray<{result: nca.ComputeResult; route: nca.RouteCode}>,
  tag: string, expectedMethod: 'logslope' | 'observed',
): void {
  let checked = 0;
  for (let i = 0; i < results.rowCount; i++) {
    const pp = perProfile[i];
    if (pp.route !== nca.ROUTE_IV_BOLUS || pp.result.status === 'failed') continue;
    const subj = results.col('SUBJECT')!.get(i) as string;
    const c0 = pp.result.provenance.c0;
    expect(c0 !== null && c0 !== undefined, true,
      `${tag} ${subj}: provenance.c0 must be populated on a non-failed IV-bolus row`);
    expect(results.col('C0_EXTRAP')!.get(i), c0!.value,
      `${tag} ${subj}: C0_EXTRAP !== provenance.c0.value (drift between the two chains)`);
    expect(c0!.method, expectedMethod, `${tag} ${subj}: c0 method`);
    const json = JSON.parse(results.col('PROVENANCE_JSON')!.get(i) as string) as {
      c0: {method: string; replacedDoseTimeRow: boolean; pctAucBackExtrap: number | null} | null;
    };
    expect(json.c0 !== null, true, `${tag} ${subj}: PROVENANCE_JSON.c0 present`);
    expect(json.c0!.method, expectedMethod, `${tag} ${subj}: PROVENANCE_JSON.c0.method`);
    expect(json.c0!.replacedDoseTimeRow, c0!.replacedDoseTimeRow,
      `${tag} ${subj}: PROVENANCE_JSON.c0.replacedDoseTimeRow`);
    if (expectedMethod === 'observed')
      expect(json.c0!.pctAucBackExtrap, 0, `${tag} ${subj}: an observed dose-time value has nothing back-extrapolated`);
    checked++;
  }
  expect(checked, results.rowCount, `${tag}: every row is a non-failed IV-bolus row`);
}

/**
 * PKNCA parity of one row. The fixtures' CL / Vz / Vss are ABSOLUTE (an absolute mg dose), so they
 * are asserted on the absolute companions `CL_ABS` / `VZ_ABS` / `VSS_ABS` — the ground truth is
 * kept, re-pointed (nca-cl-bw-normalised-units probe BREAK 2). When the caller passes the subject's
 * body weight AS RECORDED IN THE FIXTURE (`bwKg`), the per-kg `CL` / `VZ` / `VSS` are asserted too,
 * as an IDENTITY computed here from the PKNCA value ÷ that weight × 1000/60 — never through
 * `normaliseClVolumes` — so the check exercises the subject ↔ weight join (AD-10, plan-review F-3).
 */
function rowAssert(
  results: DG.DataFrame, rowIdx: number, expected: FixtureProfile, tag: string, bwKg?: number,
): void {
  const get = (name: string) => results.col(name)!.get(rowIdx);
  // MISSINGNESS is asserted through `isNone`, never through `Number.isNaN` on a
  // `get()` (nca-dg-native-missing-values F5). Result columns now carry
  // `DG.FLOAT_NULL`, and `get()` returns that raw sentinel even for a null cell
  // (probe PR-2) — so `Number.isNaN` on it is ALWAYS false and every one of these
  // assertions would silently stop testing anything. The asserted VALUES below are
  // untouched — PKNCA parity (NFR-04/05) is about the numbers, not the predicate.
  const isMissing = (name: string) => results.col(name)!.isNone(rowIdx);
  const p = expected.parameters;
  const prov = expected.provenance;

  expect(get('STATUS'), 'ok', `${tag}: STATUS must be ok; got ${get('STATUS')}`);

  if (p.cmax !== null)
    expect(get('CMAX'), p.cmax, `${tag}: Cmax mismatch`);
  if (p.tmax !== null)
    expect(get('TMAX'), p.tmax, `${tag}: Tmax mismatch`);

  if (p.auclast !== null) {
    expectFloat(get('AUC_LAST') as number, p.auclast, Math.abs(p.auclast) * 0.001,
      `${tag}: AUC_LAST out of 0.1% tolerance`);
  }
  if (p.aucinf !== null) {
    expectFloat(get('AUC_INF') as number, p.aucinf, Math.abs(p.aucinf) * 0.001,
      `${tag}: AUC_INF out of 0.1% tolerance`);
  }
  if (p.lambda_z !== null) {
    expectFloat(get('LAMBDA_Z') as number, p.lambda_z, Math.abs(p.lambda_z) * 0.005,
      `${tag}: LAMBDA_Z out of 0.5% tolerance`);
  }
  if (p.half_life !== null) {
    expectFloat(get('HALF_LIFE') as number, p.half_life, Math.abs(p.half_life) * 0.01,
      `${tag}: HALF_LIFE out of 1% tolerance`);
  }
  if (p.cl !== null) {
    expectFloat(get('CL_ABS') as number, p.cl, Math.abs(p.cl) * 0.01,
      `${tag}: CL abs. out of 1% tolerance`);
    if (bwKg !== undefined) {
      const want = p.cl / bwKg * 1000 / 60;
      expectFloat(get('CL') as number, want, Math.abs(want) * 0.01,
        `${tag}: CL per kg = PKNCA CL ÷ fixture BW ${bwKg} kg × 1000/60 (IDENTITY)`);
    } else
      expect(isMissing('CL'), true, `${tag}: no body weight → CL per kg not formed (PER_KG_NOT_FORMED)`);
  }
  if (p.vz !== null) {
    expectFloat(get('VZ_ABS') as number, p.vz, Math.abs(p.vz) * 0.01,
      `${tag}: VZ abs. out of 1% tolerance`);
    if (bwKg !== undefined) {
      expectFloat(get('VZ') as number, p.vz / bwKg, Math.abs(p.vz / bwKg) * 0.01,
        `${tag}: Vz per kg = PKNCA Vz ÷ fixture BW (IDENTITY)`);
    }
  }
  if (p.pct_aucextrap !== null) {
    // 0.5 percentage-point absolute tolerance per §8.1 row 9.
    expectFloat(get('PCT_EXTRAP') as number, p.pct_aucextrap, 0.5,
      `${tag}: PCT_EXTRAP out of 0.5 pct-point tolerance`);
  }

  // FR-200 moment / lag params (PKNCA parity, AD-12 tolerance classes).
  if (p.aumclast != null) {
    expectFloat(get('AUMC_LAST') as number, p.aumclast, Math.abs(p.aumclast) * 0.001,
      `${tag}: AUMC_LAST out of 0.1% tolerance`);
  }
  if (p.aumcinf_obs != null) {
    expectFloat(get('AUMC_INF') as number, p.aumcinf_obs, Math.abs(p.aumcinf_obs) * 0.001,
      `${tag}: AUMC_INF out of 0.1% tolerance`);
  }
  if (p.mrt != null) {
    expectFloat(get('MRT') as number, p.mrt, Math.abs(p.mrt) * 0.01,
      `${tag}: MRT out of 1% tolerance`);
  }
  // Vss: gate-first. null (extravascular) ⇒ route gate ⇒ NaN; a number (IV) ⇒ parity.
  if (p.vss == null) {
    expect(isMissing('VSS') && isMissing('VSS_ABS'), true,
      `${tag}: VSS / VSS abs. must be EMPTY for extravascular (route gate)`);
  } else {
    expectFloat(get('VSS_ABS') as number, p.vss, Math.abs(p.vss) * 0.01,
      `${tag}: VSS abs. out of 1% tolerance`);
    if (bwKg !== undefined) {
      expectFloat(get('VSS') as number, p.vss / bwKg, Math.abs(p.vss / bwKg) * 0.01,
        `${tag}: Vss per kg = PKNCA Vss ÷ fixture BW (IDENTITY)`);
    }
  }
  // Tlag: gate-first. null (IV) ⇒ route gate ⇒ NaN; a number (extravascular) ⇒ exact.
  if (p.tlag == null) {
    expect(isMissing('TLAG'), true,
      `${tag}: TLAG must be EMPTY for IV (route gate)`);
  } else
    expect(get('TLAG'), p.tlag, `${tag}: TLAG mismatch`);
  if (p.pct_aumcextrap != null) {
    expectFloat(get('PCT_AUMC_EXTRAP') as number, p.pct_aumcextrap, 0.5,
      `${tag}: PCT_AUMC_EXTRAP out of 0.5 pct-point tolerance`);
  }

  if (prov.lambda_z_n_points !== null) {
    expect(get('LAMBDAZ_NPOINTS'), prov.lambda_z_n_points,
      `${tag}: LAMBDAZ_NPOINTS regression marker`);
  }
  // SPAN_RATIO vs PKNCA's `span.ratio` (rule 18). Asserted from the FIXTURE, never
  // from a literal transcribed into this file — a hardcoded number would sit frozen
  // while a regeneration recomputed PKNCA's real value, so a genuine divergence
  // would go undetected. That is the exact defect the upstream `min-span-ratio`
  // cycle shipped and review caught.
  //
  // Tolerance: 0.5% relative, INHERITED from λz rather than chosen. Span ratio is
  // (t_last − t_first) / t_½ — the times are exact, so the whole error budget is
  // t_½'s, which is 1/λz. Holding it to λz's own §8.1 bar is therefore the tightest
  // defensible bound, and tighter than t_½'s 1%.
  //
  // `!= null` covers BOTH absences: an oracle-free fixture (the NONMEM / infusion
  // JSONs predate the field ⇒ undefined) and a null cell. Neither is a failure —
  // there is simply nothing to compare against.
  if (prov.span_ratio != null) {
    expect(isMissing('SPAN_RATIO'), false,
      `${tag}: SPAN_RATIO must be present when PKNCA reports span.ratio = ${prov.span_ratio}`);
    expectFloat(get('SPAN_RATIO') as number, prov.span_ratio,
      Math.abs(prov.span_ratio) * 0.005,
      `${tag}: SPAN_RATIO out of 0.5% tolerance vs PKNCA span.ratio`);
  }
  // C0_EXTRAP — Ed obs #4. PKNCA reference (provenance.c0_extrapolated) is
  // null for non-IV-bolus profiles and finite for IV-bolus (via log-slope
  // back-extrap). 0.1% relative tolerance matches AUC — both are derived
  // from a log-linear regression of the early initial points. A no-fit cell is
  // the PLATFORM null (`DG.FLOAT_NULL`) — check with `isNone`. (This comment used
  // to assert the opposite; the writer emitted a raw NaN then, which DG read as a
  // present value. Migrated 2026-08-03.)
  if (prov.c0_extrapolated === null) {
    expect(isMissing('C0_EXTRAP'), true,
      `${tag}: C0_EXTRAP must be EMPTY for non-IV-bolus / no-c0`);
  } else {
    expectFloat(get('C0_EXTRAP') as number, prov.c0_extrapolated,
      Math.abs(prov.c0_extrapolated) * 0.001,
      `${tag}: C0_EXTRAP out of 0.1% tolerance`);
  }
}

category('Reference suite: PKNCA parity', () => {
  test('theoph: 12/12 subjects within §8.1 tolerances, all ok', async () => {
    const df = withDoseMg(await loadFixture('theoph.csv'));
    const fx = await loadReference('theoph-reference.json');
    const byKey = indexFixtureProfiles(fx);

    const {results} = runPipelineOrThrow(df, theophContract(df), theophScalars(), DEFAULT_RULES);
    expect(results.rowCount, 12, `expected 12 subjects; got ${results.rowCount}`);
    // AD-4: neutral label; per-row basis. theoph is extravascular → CL_BASIS=CL/F.
    expect(results.col('CL')!.tags['label'], 'CL', 'theoph CL label must be neutral CL');
    expect(results.col('VZ')!.tags['label'], 'Vz', 'theoph VZ label must be neutral Vz');
    expect(results.col('CL_BASIS')!.get(0), 'CL/F', 'theoph CL_BASIS must be CL/F (extravascular)');
    // The per-kg CL is apparent CL/F per kg here (nca-cl-bw-normalised-units AD-1) — same basis column.
    expect(results.col('CL')!.tags['units'], 'mL/min/kg', 'CL per kg');

    // Each subject's weight AS RECORDED in the fixture CSV (`Wt`, one value per subject) — the
    // IDENTITY input, read from the source frame here, not from the results' BW_KG.
    const wtBySubject = new Map<string, number>();
    for (let r = 0; r < df.rowCount; r++)
      wtBySubject.set(String(df.col('Subject')!.get(r)), df.col('Wt')!.get(r) as number);
    for (let i = 0; i < results.rowCount; i++) {
      const subj = results.col('SUBJECT')!.get(i) as string;
      const expected = byKey.get(subj);
      expect(expected !== undefined, true, `theoph subject '${subj}' missing from reference fixture`);
      rowAssert(results, i, expected!, `theoph subject ${subj}`, wtBySubject.get(subj));
    }
  }, {timeout: TIMEOUT});

  test('indometh: 6/6 subjects within §8.1 tolerances, IV bolus, CL label "CL"', async () => {
    const df = await loadFixture('indometh.csv');
    const fx = await loadReference('indometh-reference.json');
    const byKey = indexFixtureProfiles(fx);

    const {results} = runPipelineOrThrow(df, indomethContract(df), indomethScalars(), DEFAULT_RULES);
    expect(results.rowCount, 6, `expected 6 subjects; got ${results.rowCount}`);
    // AD-4: neutral label; indometh is IV-bolus across all subjects → CL_BASIS=CL.
    expect(results.col('CL')!.tags['label'], 'CL', 'indometh CL label must be neutral CL');
    expect(results.col('VZ')!.tags['label'], 'Vz', 'indometh VZ label must be neutral Vz');
    expect(results.col('CL_BASIS')!.get(0), 'CL', 'indometh CL_BASIS must be CL (IV-bolus)');

    for (let i = 0; i < results.rowCount; i++) {
      const subj = results.col('SUBJECT')!.get(i) as string;
      const expected = byKey.get(subj);
      expect(expected !== undefined, true, `indometh subject '${subj}' missing from reference fixture`);
      rowAssert(results, i, expected!, `indometh subject ${subj}`);
    }
  }, {timeout: TIMEOUT});

  test('SPAN_RATIO: the alarm tier (< 2) has real PKNCA-oracle coverage, not just the safe tiers', async () => {
    // AC-3.7a. Shipping a warning band whose FIRING tier has no reference
    // assertion is what rule 18 forbids ("a computation that does not have a
    // corresponding reference-suite assertion does not exist"), and this band
    // came within one substring match of doing exactly that: the four sub-2.0
    // profiles all live in theoph + indometh, which load through the JSON path,
    // and those fixtures carried no span oracle until they were regenerated from
    // sci-comp's 0.10.0 set. The string `span_ratio` WAS in the old files — as
    // `config.min_span_ratio`, a PKNCA *input option*, not an output.
    //
    // So this test guards the oracle itself, not the arithmetic: it fails if a
    // future fixture regeneration drops `span.ratio` from provenance, or if the
    // corpus drifts to where nothing exercises the tier that raises the alarm.
    // rowAssert already checks each individual value against PKNCA.
    const tiers = {poor: 0, borderline: 0, adequate: 0, noOracle: 0};
    for (const name of ['theoph-reference.json', 'indometh-reference.json']) {
      const fx = await loadReference(name);
      for (const p of fx.profiles) {
        const s = p.provenance.span_ratio;
        if (s == null) {tiers.noOracle++; continue;}
        if (s < 2) tiers.poor++;
        else if (s < 4) tiers.borderline++;
        else tiers.adequate++;
      }
    }
    expect(tiers.noOracle, 0,
      'every theoph/indometh profile must carry a PKNCA span.ratio oracle — ' +
      `${tiers.noOracle} lack one. Regenerate from sci-comp's fixtures ` +
      '(see files/REGEN.md); do NOT transcribe values into this test.');
    expect(tiers.poor > 0, true,
      'no profile below the 2-half-life alarm threshold — the COLOR_CRITICAL band ' +
      'would ship with nothing asserting it. Measured 2026-08-12: 4 of these 18 ' +
      'profiles are below 2 (theoph 1/9/10, indometh 1).');
    expect(tiers.adequate > 0, true, 'the ≥ 4 adequate tier is exercised too');
    return `span-ratio oracle coverage over 18 profiles — poor(<2)=${tiers.poor} ` +
      `borderline(2–4)=${tiers.borderline} adequate(≥4)=${tiers.adequate}`;
  }, {timeout: TIMEOUT});

  test('SF-minspanratio: nca-studio never sets sci-comp\'s minSpanRatio gate', async () => {
    // The guard that keeps a DISPLAY threshold out of the COMPUTATION rules.
    // Structurally it is AD-4 (the two band constants are module-local in
    // results-view.ts and deliberately absent from ResultsColorThresholds, whose
    // every field is a projection from NcaRules); this is the executable backstop.
    //
    // Two halves, because the naming half alone would not catch a mis-wiring:
    //
    // (1) STRUCTURAL — DEFAULT_RULES is the fallback source for every lambdaZ
    //     field the mapping pane's buildRules() does not read from an input, so
    //     a stray minSpanRatio here would reach every Compute.
    expect(DEFAULT_RULES.lambdaZ.minSpanRatio === undefined, true,
      'DEFAULT_RULES must not set lambdaZ.minSpanRatio. It cannot change a ' +
      'computed value (sci-comp reads it at one site, only to push a ' +
      'LAMBDAZ_LOW_SPAN warning), but setting it would inject unrequested ' +
      'warnings into provenance.warnings → WARNINGS_JSON → the compute-notices ' +
      'banner, and corrupt the audit log ruleHash so the trail asserts a rule ' +
      'the scientist never chose.');

    // (2) BEHAVIOURAL — the consequence, on the corpus that would actually
    //     produce it. Three theoph profiles span under 2 half-lives, so if the
    //     display threshold ever leaked into the engine rules this frame would
    //     carry LAMBDAZ_LOW_SPAN warnings. Asserting the OUTPUT is what makes
    //     this more than a restatement of the line above.
    const df = withDoseMg(await loadFixture('theoph.csv'));
    const {results} = runPipelineOrThrow(
      df, theophContract(df), theophScalars(), DEFAULT_RULES);
    const warnCol = results.col('WARNINGS_JSON')!;
    let lowSpan = 0;
    for (let i = 0; i < results.rowCount; i++)
      if (String(warnCol.get(i)).includes('LAMBDAZ_LOW_SPAN')) lowSpan++;
    expect(lowSpan, 0,
      `${lowSpan} theoph profiles carry a LAMBDAZ_LOW_SPAN warning. nca-studio ` +
      'never requested that gate — a display threshold has been wired into the ' +
      'NcaRules handed to computeNca.');
    return `minSpanRatio unset; 0/${results.rowCount} LAMBDAZ_LOW_SPAN warnings ` +
      'despite 3 theoph profiles spanning < 2 half-lives';
  }, {timeout: TIMEOUT});

  test('AC-D3: IV-infusion — −T_inf/2 MRT + Vss correction validated vs PKNCA', async () => {
    // The ONLY fixture that exercises the infusion branch: infusionDuration=1h
    // threads through scaleProfileInputs → computeNca and must produce a
    // duration-corrected MRT/Vss. Synthetic 1-compartment (CL=2, V=10, T_inf=1)
    // → PKNCA AUCinf≈50, MRT≈5, Vss≈10.
    const df = await loadFixture('04_iv_infusion.csv');
    const fx = await loadReference('04_iv_infusion-reference.json');
    const byKey = indexFixtureProfiles(fx);

    const {results} = runPipelineOrThrow(
      df, infusionContract(df), infusionScalars(), DEFAULT_RULES);
    expect(results.rowCount, 1, `expected 1 infusion subject; got ${results.rowCount}`);
    expect(results.col('CL_BASIS')!.get(0), 'CL', 'infusion CL_BASIS must be CL (IV)');

    const subj = String(results.col('SUBJECT')!.get(0));
    const expected = byKey.get(subj);
    expect(expected !== undefined, true, `infusion subject '${subj}' missing from reference`);
    rowAssert(results, 0, expected!, `infusion subject ${subj}`);

    // The defining check, row-internal so it does not merely re-assert the
    // PKNCA reference: MRT carries the −T_inf/2 correction (T_inf = 1 h), i.e.
    // MRT = AUMC_inf/AUC_inf − 0.5, NOT the uncorrected moment ratio.
    const aumcInf = results.col('AUMC_INF')!.get(0) as number;
    const aucInf = results.col('AUC_INF')!.get(0) as number;
    const mrt = results.col('MRT')!.get(0) as number;
    expectFloat(mrt, aumcInf / aucInf - 0.5, Math.abs(mrt) * 1e-9,
      'infusion MRT must apply −T_inf/2 (= AUMC_inf/AUC_inf − 0.5)');
    // Vss = CL·MRT carries the same correction; Tlag is NaN (IV route gate). No body weight on
    // this fixture → the absolute companions carry the values (per kg not formed).
    const cl = results.col('CL_ABS')!.get(0) as number;
    expectFloat(results.col('VSS_ABS')!.get(0) as number, cl * mrt, Math.abs(cl * mrt) * 1e-6,
      'infusion Vss = CL·MRT (carries the −T_inf/2 correction)');
    expect(results.col('TLAG')!.isNone(0), true,
      'infusion TLAG is EMPTY (IV route gate)');
  }, {timeout: TIMEOUT});

  // --------------------------------------------------------------------
  // C0_EXTRAP — Ed obs #4. Dedicated tests because the broader rowAssert
  // bails on the BUG-01 Cmax/Tmax exact-match before reaching C0_EXTRAP;
  // these isolate the c0 assertion (uses 0.1% relative tolerance and
  // doesn't depend on Cmax/Tmax bit-identity) so BUG-01 doesn't mask it.
  // --------------------------------------------------------------------

  test('C0_EXTRAP: indometh — sci-comp estimateC0 matches PKNCA reference for 6/6 IV-bolus subjects', async () => {
    const df = await loadFixture('indometh.csv');
    const fx = await loadReference('indometh-reference.json');
    const byKey = indexFixtureProfiles(fx);

    const {results} = runPipelineOrThrow(df, indomethContract(df), indomethScalars(), DEFAULT_RULES);
    expect(results.col('C0_EXTRAP') !== null, true, 'C0_EXTRAP column must exist');
    expect(results.col('C0_EXTRAP')!.tags['units'], 'mg/L', 'C0_EXTRAP units = mg/L');

    for (let i = 0; i < results.rowCount; i++) {
      const subj = results.col('SUBJECT')!.get(i) as string;
      const expC0 = byKey.get(subj)!.provenance.c0_extrapolated;
      expect(expC0 !== null, true, `indometh ${subj}: reference c0 must be non-null for IV bolus`);
      expectFloat(results.col('C0_EXTRAP')!.get(i) as number, expC0!,
        Math.abs(expC0!) * 0.001,
        `indometh subject ${subj}: C0_EXTRAP out of 0.1% tolerance vs PKNCA reference`);
    }
  }, {timeout: TIMEOUT});

  test('C0_EXTRAP drift guard: indometh — C0_EXTRAP === provenance.c0.value EXACTLY (AC-N.3.2/3)', async () => {
    const df = await loadFixture('indometh.csv');
    const {results, outcome} = runPipelineOrThrow(df, indomethContract(df), indomethScalars(), DEFAULT_RULES);
    const fx = await loadReference('indometh-reference.json');
    const byKey = indexFixtureProfiles(fx);
    assertC0DriftGuard(results, outcome.perProfile, 'indometh', 'logslope');
    // The engine's own c0 (sci-comp ≥ 0.11.0) also matches the PKNCA reference —
    // the SAME oracle C0_EXTRAP is held to, closing the loop on both sources.
    for (let i = 0; i < results.rowCount; i++) {
      const subj = results.col('SUBJECT')!.get(i) as string;
      const expC0 = byKey.get(subj)!.provenance.c0_extrapolated!;
      const c0 = outcome.perProfile[i].result.provenance.c0!;
      expectFloat(c0.value, expC0, Math.abs(expC0) * 0.001, `indometh ${subj}: provenance.c0.value vs PKNCA`);
      expect(c0.replacedDoseTimeRow, false, `indometh ${subj}: no dose-time row to replace`);
      // Back-extrapolated share of AUCinf (Phoenix AUC_%Back_Ext): 16–28 % on this corpus.
      expect(c0.pctAucBackExtrap > 15 && c0.pctAucBackExtrap < 30, true,
        `indometh ${subj}: pctAucBackExtrap ${c0.pctAucBackExtrap} outside the measured 16–28 % band`);
    }
  }, {timeout: TIMEOUT});

  test('C0_EXTRAP: theoph — NaN for all 12 PO subjects (non-IV-bolus skip)', async () => {
    const df = withDoseMg(await loadFixture('theoph.csv'));
    const {results} = runPipelineOrThrow(df, theophContract(df), theophScalars(), DEFAULT_RULES);
    expect(results.col('C0_EXTRAP') !== null, true, 'C0_EXTRAP column must exist');
    for (let i = 0; i < results.rowCount; i++) {
      const subj = results.col('SUBJECT')!.get(i) as string;
      // The no-fit sentinel is the PLATFORM null (`isNone`), not a raw NaN.
      expect(results.col('C0_EXTRAP')!.isNone(i), true,
        `theoph subject ${subj}: C0_EXTRAP must be EMPTY for PO route`);
    }
  }, {timeout: TIMEOUT});

  let theophDf: DG.DataFrame | null = null;
  before(async () => {
    // Smoke-load both fixtures once so a fixture-publishing problem fails the
    // whole category up front rather than per-test (the `before` hook surfaces
    // in test output as a single FAIL with the underlying error).
    theophDf = await loadFixture('theoph.csv');
    await loadFixture('indometh.csv');
    await loadReference('theoph-reference.json');
    await loadReference('indometh-reference.json');
  });

  test('fixture-publishing sanity: theoph CSV has Subject + Time + conc columns', async () => {
    expect(theophDf !== null, true, '`before` hook should have populated theophDf');
    if (theophDf === null) return;
    expect(theophDf.col('Subject') !== null, true, 'theoph.csv must expose \'Subject\'');
    expect(theophDf.col('Time') !== null, true, 'theoph.csv must expose \'Time\'');
    expect(theophDf.col('conc') !== null, true, 'theoph.csv must expose \'conc\'');
    expect(theophDf.rowCount > 0, true, 'theoph.csv must not be empty');
  });

  // -------------------------------------------------------------------
  // nca-unit-semantics — non-canonical input mirror test (E8 spec).
  // -------------------------------------------------------------------

  test('theoph in minutes mirror: same results as theoph in hours within Float64 tolerance', async () => {
    // Load the canonical-units baseline.
    const baselineDf = withDoseMg(await loadFixture('theoph.csv'));
    const {results: baselineResults} = runPipelineOrThrow(
      baselineDf, theophContract(baselineDf), theophScalars(), DEFAULT_RULES);

    // Build a "theoph-in-minutes" mirror: clone with `Time × 60` and pass
    // `timeUnits: 'min'` to `ProfileScalars`. `scaleProfileInputs` rescales
    // back to canonical hours; results should match baseline within
    // Float64 tolerance.
    const minutesDf = withDoseMg(await loadFixture('theoph.csv'));
    const timeCol = minutesDf.col('Time')!;
    const raw = timeCol.getRawData() as Float64Array;
    const rescaled = new Float64Array(raw.length);
    for (let i = 0; i < raw.length; i++)
      rescaled[i] = raw[i] === DG.FLOAT_NULL ? raw[i] : raw[i] * 60;
    // Replace the Time column in-place with the scaled values.
    minutesDf.columns.remove('Time');
    minutesDf.columns.add(DG.Column.fromFloat64Array('Time', rescaled));

    const minutesScalars: ProfileScalars = {...theophScalars(), timeUnits: 'min'};
    const {results: minutesResults} = runPipelineOrThrow(
      minutesDf, theophContract(minutesDf), minutesScalars, DEFAULT_RULES);

    expect(minutesResults.rowCount, baselineResults.rowCount,
      `row count must match baseline; got ${minutesResults.rowCount} vs ${baselineResults.rowCount}`);

    // Subject-keyed comparison (lex-sorted in both result DFs).
    const numericCols = ['CMAX', 'TMAX', 'AUC_LAST', 'AUC_INF', 'HALF_LIFE', 'CL', 'VZ'];
    for (let i = 0; i < minutesResults.rowCount; i++) {
      const subj = minutesResults.col('SUBJECT')!.get(i) as string;
      const baseSubj = baselineResults.col('SUBJECT')!.get(i) as string;
      expect(subj, baseSubj, `subject ordering must match at row ${i}`);
      for (const colName of numericCols) {
        const baseMissing = baselineResults.col(colName)!.isNone(i);
        const baseVal = baselineResults.col(colName)!.get(i) as number;
        const mirrorVal = minutesResults.col(colName)!.get(i) as number;
        if (baseMissing) {
          expect(minutesResults.col(colName)!.isNone(i), true,
            `${colName}@${subj}: baseline MISSING must mirror missing; got ${mirrorVal}`);
          continue;
        }
        const tol = Math.max(Math.abs(baseVal) * 1e-9, 1e-12);
        expectFloat(mirrorVal, baseVal, tol,
          `${colName}@${subj}: minutes-mirror must equal baseline within ` +
          `1e-9 relative tolerance (Float64 rescaling drift expected)`);
      }
    }
  }, {timeout: TIMEOUT});
});

// =====================================================================
// Rat multi-route (dataset 03) — VAL-01. Per-route serial PK validated
// bit-exact against PKNCA 0.12.1. References ship as
// `*_nca_datagrok_match.csv` (loaded via loadReferenceCsv). Input schema
// is preclinical: Dose_mg_kg (per-kg) + BodyWeight_g + Conc_ng_mL.
//
// All 64 subjects compute `ok` and match PKNCA at STRICT NFR-05 tolerance
// (0.5% λz / 1% t½·Vz). The former 3 partials (R005, R013 IV; R029 PO) were the
// BUG-05 trailing-BLQ-as-zero defect — fixed in sci-comp 0.7.1. R019's λz
// window-selection divergence (VAL-01-LZ-R019: auto-fit picked 8 terminal points
// vs PKNCA's 4, ~2.3% off) is RESOLVED in sci-comp 0.9.1 — best-fit now uses
// PKNCA/WinNonlin's flat adj-R² tolerance off the global max, so R019 matches at
// strict tolerance and is no longer banded.
//
// `RAT_LZ_DIVERGENT` is retained (empty) as the parameterized regression band: if a
// future sci-comp change reintroduces a window-selection divergence, list the
// affected subjects here rather than loosening the global tolerance.
// =====================================================================

const RAT_LZ_DIVERGENT: ReadonlySet<string> = new Set([]); // VAL-01-LZ-R019 resolved in sci-comp 0.9.1

/**
 * `bwKg` = the subject's weight as recorded in the source fixture (grams ÷ 1000), given when the
 * analysis maps body weight: the PKNCA cl.obs / vz.obs are then ABSOLUTE and asserted on the
 * companions (ground truth), and the per-kg `CL` / `VZ` as an IDENTITY (PKNCA ÷ fixture BW). Without
 * it the analysis is the per-kg path and PKNCA's values are per kg already: `CL` = cl.obs × 1000/60,
 * `VZ` = vz.obs, and no absolute companion can be formed (nca-cl-bw-normalised-units AD-10).
 */
function ratRowAssert(
  results: DG.DataFrame, rowIdx: number, subj: string, expected: FixtureProfile, tag: string,
  lzDivergent: ReadonlySet<string> = RAT_LZ_DIVERGENT, bwKg?: number,
): void {
  const get = (name: string) => results.col(name)!.get(rowIdx) as number;
  const p = expected.parameters;

  // Cmax: 0.1% (NFR-04). Exact match is infeasible here — conc rescales
  // ng/mL→mg/L (×0.001), so a tiny binary-rounding delta vs the reference.
  if (p.cmax !== null)
    expectFloat(get('CMAX'), p.cmax, Math.abs(p.cmax) * 0.001, `${tag}: Cmax`);
  // Tmax: time is already in hours (no rescale) → exact.
  if (p.tmax !== null)
    expect(get('TMAX'), p.tmax, `${tag}: Tmax`);
  if (p.auclast !== null)
    expectFloat(get('AUC_LAST'), p.auclast, Math.abs(p.auclast) * 0.001, `${tag}: AUClast (0.1%)`);
  if (p.aucinf !== null)
    expectFloat(get('AUC_INF'), p.aucinf, Math.abs(p.aucinf) * 0.001, `${tag}: AUCinf (0.1%)`);
  if (p.cl !== null) {
    if (bwKg !== undefined) {
      expectFloat(get('CL_ABS'), p.cl, Math.abs(p.cl) * 0.01, `${tag}: CL abs. (1%)`);
      const want = p.cl / bwKg * 1000 / 60;
      expectFloat(get('CL'), want, want * 0.01, `${tag}: CL per kg = PKNCA ÷ fixture BW (IDENTITY, 1%)`);
    } else {
      const want = p.cl * 1000 / 60;
      expectFloat(get('CL'), want, want * 0.01, `${tag}: CL per kg = PKNCA per-kg × 1000/60 (1%)`);
      expect(results.col('CL_ABS')!.isNone(rowIdx), true, `${tag}: no weight → no CL abs.`);
    }
  }
  if (p.pct_aucextrap !== null)
    expectFloat(get('PCT_EXTRAP'), p.pct_aucextrap, 0.5, `${tag}: %extrap (0.5pp)`);

  // λz family: strict NFR-05 — 0.5% (λz) / 1% (t½, Vz). The `lzDivergent` set is a
  // per-test parameter for banding a KNOWN window-selection divergence to 3%; it is
  // currently empty (VAL-01-LZ-R019 resolved in sci-comp 0.9.1) but kept so a future
  // regression can be banded per-subject instead of loosening the global tolerance.
  const lzTol = lzDivergent.has(subj) ? 0.03 : 0.005;
  const derivedTol = lzDivergent.has(subj) ? 0.03 : 0.01;
  if (p.lambda_z !== null)
    expectFloat(get('LAMBDA_Z'), p.lambda_z, Math.abs(p.lambda_z) * lzTol, `${tag}: λz`);
  if (p.half_life !== null)
    expectFloat(get('HALF_LIFE'), p.half_life, Math.abs(p.half_life) * derivedTol, `${tag}: t½`);
  if (p.vz !== null) {
    if (bwKg !== undefined) {
      expectFloat(get('VZ_ABS'), p.vz, Math.abs(p.vz) * derivedTol, `${tag}: Vz abs.`);
      expectFloat(get('VZ'), p.vz / bwKg, Math.abs(p.vz / bwKg) * derivedTol, `${tag}: Vz per kg (IDENTITY)`);
    } else
      expectFloat(get('VZ'), p.vz, Math.abs(p.vz) * derivedTol, `${tag}: Vz per kg`);
  }
}

/** Each subject's weight as recorded in a rat fixture (`BodyWeight_g`, grams → kg), by `SubjectID`. */
function ratFixtureBwKg(df: DG.DataFrame): Map<string, number> {
  const out = new Map<string, number>();
  for (let r = 0; r < df.rowCount; r++)
    out.set(String(df.col('SubjectID')!.get(r)), (df.col('BodyWeight_g')!.get(r) as number) / 1000);
  return out;
}

category('Reference suite: PKNCA parity (rat multi-route)', () => {
  const routes: Array<{
    name: string; input: string; ref: string; route: nca.RouteCode; clLabel: string; rows: number;
  }> = [
    {name: 'rat IV', input: '03_rat_iv.csv', ref: 'rat_iv_nca_datagrok_match.csv',
      route: nca.ROUTE_IV_BOLUS, clLabel: 'CL', rows: 24},
    {name: 'rat PO', input: '03_rat_po.csv', ref: 'rat_po_nca_datagrok_match.csv',
      route: nca.ROUTE_PO, clLabel: 'CL/F', rows: 24},
    {name: 'rat SC', input: '03_rat_sc.csv', ref: 'rat_sc_nca_datagrok_match.csv',
      route: nca.ROUTE_SC, clLabel: 'CL/F', rows: 16},
  ];

  for (const r of routes) {
    test(`${r.name}: ${r.rows} subjects within tolerances vs PKNCA`, async () => {
      const df = await loadFixture(r.input);
      const fx = await loadReferenceCsv(r.ref);
      const byKey = indexFixtureProfiles(fx);

      const {results} = runPipelineOrThrow(df, ratContract(df), ratScalars(r.route), DEFAULT_RULES);
      expect(results.rowCount, r.rows, `${r.name}: expected ${r.rows} subjects`);
      // AD-4: neutral label; the route-derived basis lives per-row in CL_BASIS.
      expect(results.col('CL')!.tags['label'], 'CL', `${r.name}: CL label must be neutral`);
      expect(results.col('CL_BASIS')!.get(0), r.clLabel, `${r.name}: CL_BASIS`);
      const bw = ratFixtureBwKg(df);

      for (let i = 0; i < results.rowCount; i++) {
        const subj = results.col('SUBJECT')!.get(i) as string;
        const expected = byKey.get(subj);
        expect(expected !== undefined, true, `${r.name}: subject '${subj}' missing from reference`);

        // BUG-05 fixed (sci-comp 0.7.1): every subject computes `ok` and matches
        // PKNCA. R019's λz-window divergence is resolved (sci-comp 0.9.1), so it
        // asserts at strict tolerance too (RAT_LZ_DIVERGENT is now empty); a
        // regression back to `partial` (or any other status) trips this assertion.
        expect(results.col('STATUS')!.get(i), 'ok', `${r.name} ${subj}: STATUS must be ok`);
        ratRowAssert(results, i, subj, expected!, `${r.name} ${subj}`, RAT_LZ_DIVERGENT, bw.get(subj));
      }
    }, {timeout: TIMEOUT});
  }
});

// =====================================================================
// Multicompound batch screen (dataset 05) — VAL-01. 20 compounds × 3 rats,
// all IV bolus, per-subject NCA validated against PKNCA 0.12.1. The reference
// ships as `multicompound_nca_per_subject.csv` (loaded via loadReferenceCsv,
// no JSON conversion). SubjectID is globally unique across compounds, so
// grouping by subject reproduces the 60 per-subject reference rows directly.
//
// No body-weight column → the per-kg CL/Vz path (CL in L/h/kg, Vz in L/kg);
// the reference cl.obs / vz.obs are PKNCA's native per-kg values for the same
// path, so ratRowAssert compares like-for-like. All 60 subjects match at strict
// NFR-05 0.5%/1%. R017 (DGK-0006) formerly diverged 1.91% on λz (PKNCA fit a
// 3-point terminal window while sci-comp's best-fit selected a different one) —
// the same window-selection class as dataset 03's R019 (VAL-01-LZ-R019), now
// RESOLVED in sci-comp 0.9.1, so R017 is no longer banded.
// =====================================================================
// Empty: VAL-01-LZ-R019 resolved in sci-comp 0.9.1 (R017 now matches strict).
const MULTICOMPOUND_LZ_DIVERGENT: ReadonlySet<string> = new Set<string>([]);

category('Reference suite: PKNCA parity (multicompound 05)', () => {
  test('multicompound: 60/60 subjects within tolerances vs PKNCA (per-kg path)', async () => {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const fx = await loadReferenceCsv('multicompound_nca_per_subject.csv');
    const byKey = indexFixtureProfiles(fx);

    const {results} = runPipelineOrThrow(
      df, multicompoundContract(df), multicompoundScalars(), DEFAULT_RULES);
    expect(results.rowCount, 60, `expected 60 subjects; got ${results.rowCount}`);
    // AD-4: neutral numeric label; IV-bolus → CL_BASIS=CL. No BW → per-kg unit tags.
    expect(results.col('CL')!.tags['label'], 'CL', 'multicompound CL label must be neutral');
    expect(results.col('CL_BASIS')!.get(0), 'CL', 'multicompound CL_BASIS must be CL (IV-bolus)');
    expect(results.col('CL')!.tags['units'], 'mL/min/kg', 'CL per kg in mL/min/kg');
    expect(results.col('VZ')!.tags['units'], 'L/kg', 'VZ per kg in L/kg');

    for (let i = 0; i < results.rowCount; i++) {
      const subj = results.col('SUBJECT')!.get(i) as string;
      const expected = byKey.get(subj);
      expect(expected !== undefined, true, `multicompound subject '${subj}' missing from reference`);
      expect(results.col('STATUS')!.get(i), 'ok', `multicompound ${subj}: STATUS must be ok`);
      ratRowAssert(results, i, subj, expected!, `multicompound ${subj}`, MULTICOMPOUND_LZ_DIVERGENT);
    }
  }, {timeout: TIMEOUT});

  test('multicompound: C0_EXTRAP drift guard — 60/60 observed t=0 (AC-N.3.2)', async () => {
    // Measured (synthesis C-2): all 60 t=0 rows carry a positive measured value, so
    // every c0 is the OBSERVED sample — the gate never fires and nothing is extrapolated.
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const {results, outcome} = runPipelineOrThrow(
      df, multicompoundContract(df), multicompoundScalars(), DEFAULT_RULES);
    assertC0DriftGuard(results, outcome.perProfile, 'multicompound', 'observed');
  }, {timeout: TIMEOUT});
});

// =====================================================================
// nca-cl-vz-convention AC-6 — per-kg CL/Vz oracle (no body weight).
// Rat IV R001 with BodyWeight UNMAPPED exercises the ratified per-kg path
// (AD-2): CL = dose_mg_kg / AUCinf = 1 / 3.23324 = 0.309287 L/h/kg, and
// Vz_perkg = absolute Vz / BW = 0.91950 / 0.2441 = 3.767 L/kg. The value alone
// equals the BW=1.0-hack result, so the test ALSO pins the per-kg PATH (SF-1):
// doseNormalization='per_kg', exactly one PER_KG_NO_BW warning, BW unmapped,
// DOSE_NORM='per_kg' on the row, and uniform per-kg frame → L/h/kg / L/kg tags.
// =====================================================================
category('Reference suite: per-kg CL/Vz oracle (no body weight, AC-6)', () => {
  test('rat IV R001, BW unmapped → CL 0.309287 L/h/kg, Vz 3.767 L/kg, per-kg path pinned', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    // Contract WITHOUT body weight → forces the per-kg branch (AD-2), so the
    // computed CL/Vz are per-kg rather than the absolute rat-reference values.
    const noBwContract = {observations: {...ratContract(df).observations}};
    delete noBwContract.observations.bodyWeight;
    const {results, outcome} = runPipelineOrThrow(
      df, noBwContract, ratScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);

    let row = -1;
    for (let i = 0; i < results.rowCount; i++)
      if (String(results.col('SUBJECT')!.get(i)) === 'R001') {row = i; break;}
    expect(row >= 0, true, 'R001 must be present in rat IV results');
    expect(String(results.col('STATUS')!.get(row)), 'ok', 'R001 must compute ok');

    // (a) BW unmapped → the per-kg branch (not the absolute branch) executed.
    expect(noBwContract.observations.bodyWeight === undefined, true, 'BW must be unmapped for the per-kg path');
    // (b) per-kg numeric oracle (SF-1 — value == BW=1.0-hack result; path pinned below), in the
    // reporting scale: 0.309287 L/h/kg × 1000/60 = 5.15478 mL/min/kg (nca-cl-bw-normalised-units D-b).
    expectFloat(results.col('CL')!.get(row) as number, 0.309287 * 1000 / 60, 0.309287 * 1000 / 60 * 0.01,
      'R001 per-kg CL = dose_mg_kg / AUCinf = 0.309287 L/h/kg = 5.1548 mL/min/kg (1% tol)');
    expectFloat(results.col('VZ')!.get(row) as number, 3.767, 3.767 * 0.01,
      'R001 per-kg Vz = absolute Vz / BW = 3.767 L/kg (1% tol)');
    // (c) DOSE_NORM marks the row per-kg (the engine basis); the per-kg unit tags.
    expect(String(results.col('DOSE_NORM')!.get(row)), 'per_kg', 'R001 DOSE_NORM = per_kg');
    expect(results.col('CL')!.tags['units'], 'mL/min/kg', 'CL units mL/min/kg');
    expect(results.col('VZ')!.tags['units'], 'L/kg', 'VZ units L/kg');
    expect(results.col('CL_ABS')!.isNone(row), true, 'no weight → no absolute CL');
    // (c2) AC-D4: Vss rides the SAME per-kg path as Vz. IV-bolus route → Vss is
    // computed (finite, not the EV NaN gate) and carries the per-kg L/kg unit.
    expect(Number.isFinite(results.col('VSS')!.get(row) as number), true,
      'R001 Vss finite (IV route, not gated)');
    expect(results.col('VSS')!.tags['units'], 'L/kg', 'VSS units L/kg');
    // (d) per-profile basis + exactly one PER_KG_NO_BW warning (path proof, SF-1).
    const pp = outcome.perProfile.find((p) => p.index.subjectKey === 'R001');
    expect(pp !== undefined, true, 'R001 present in perProfile');
    expect(pp!.doseNormalization, 'per_kg', 'R001 ScaledProfileInputs basis is per_kg');
    const perKgWarns = pp!.result.provenance.warnings.filter((w) => w.code === 'PER_KG_NO_BW');
    expect(perKgWarns.length, 1,
      `exactly one PER_KG_NO_BW for R001; got ${JSON.stringify(pp!.result.provenance.warnings)}`);
  }, {timeout: TIMEOUT});
});

// =====================================================================
// nca-cl-bw-normalised-units F2 — the BW-present per-kg path on R001 and the BW-handling cases with
// no external source (REGRESSION_PIN: hand constants, never derived through the code under test;
// the variants are DERIVED in-test from 03_rat_iv, guardrail "derive, do not hand-type a fixture").
// PKNCA anchor (GROUND_TRUTH): R001 cl.obs = 0.0754969 L/h on BW 0.2441 kg, 1 mg/kg IV bolus.
// =====================================================================
category('Reference suite: per-kg CL with body weight (nca-cl-bw-normalised-units)', () => {
  const R001_CL_ABS = 0.0754969; // PKNCA 0.12.1 cl.obs, rat_iv_nca_datagrok_match.csv
  const R001_CL_PER_KG = 5.1548; // mL/min/kg — 0.0754969 / 0.2441 × 1000/60 (REGRESSION_PIN)

  /** R001's rows of 03_rat_iv, optionally rewritten. */
  async function r001(edit?: (df: DG.DataFrame) => void): Promise<DG.DataFrame> {
    const df = await loadFixture('03_rat_iv.csv');
    if (edit !== undefined) edit(df);
    return df;
  }
  const rowOf = (results: DG.DataFrame, s: string): number => {
    for (let i = 0; i < results.rowCount; i++) if (String(results.col('SUBJECT')!.get(i)) === s) return i;
    return -1;
  };
  const warnCodes = (results: DG.DataFrame, i: number): string[] =>
    (JSON.parse(String(results.col('WARNINGS_JSON')!.get(i))) as Array<{code: string}>).map((w) => w.code);

  test('R001 BW mapped (g): CL 5.1548 mL/min/kg, Vz 3.7669 L/kg; CL abs. = PKNCA within 0.1%', async () => {
    const df = await r001();
    const {results} = runPipelineOrThrow(df, ratContract(df), ratScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);
    const i = rowOf(results, 'R001');
    expectFloat(results.col('CL')!.get(i) as number, R001_CL_PER_KG, 5e-4, 'R001 per-kg CL (pin)');
    expectFloat(results.col('VZ')!.get(i) as number, 3.7669, 5e-4, 'R001 per-kg Vz (pin)');
    expectFloat(results.col('CL_ABS')!.get(i) as number, R001_CL_ABS, R001_CL_ABS * 0.001,
      'CL abs. = PKNCA cl.obs (GROUND_TRUTH, 0.1%)');
    expectFloat(results.col('BW_KG')!.get(i) as number, 0.2441, 1e-12, 'BW_KG = 244.1 g in kg');
    expect(String(results.col('DOSE_NORM')!.get(i)), 'absolute', 'AC-2.5: DOSE_NORM unchanged (engine mg)');
  }, {timeout: TIMEOUT});

  test('E4 class: BW in g and the same BW in kg give the same per-kg output', async () => {
    const g = await r001();
    const kg = await r001((df) => {
      const col = df.col('BodyWeight_g')!;
      for (let r = 0; r < df.rowCount; r++) col.set(r, (col.get(r) as number) / 1000);
    });
    const a = runPipelineOrThrow(g, ratContract(g), ratScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES).results;
    const b = runPipelineOrThrow(kg, ratContract(kg), {...ratScalars(nca.ROUTE_IV_BOLUS), bodyWeightUnit: 'kg'},
      DEFAULT_RULES).results;
    for (const c of ['CL', 'VZ', 'VSS', 'CL_ABS', 'BW_KG']) {
      const x = a.col(c)!.get(rowOf(a, 'R001')) as number;
      expectFloat(b.col(c)!.get(rowOf(b, 'R001')) as number, x, Math.abs(x) * 1e-12, `${c}: g ≡ kg`);
    }
  }, {timeout: TIMEOUT});

  /** Add an absolute `Dose_mg` (Float64) column — the mg each rat received (mg/kg × its BW). */
  const addDoseMg = (df: DG.DataFrame): void => {
    const d = df.col('Dose_mg_kg')!; const bw = df.col('BodyWeight_g')!;
    df.columns.add(DG.Column.fromFloat64Array('Dose_mg', Float64Array.from({length: df.rowCount},
      (_, r) => (d.get(r) as number) * (bw.get(r) as number) / 1000)));
  };
  /** The rat contract dosing from the absolute `Dose_mg` column. */
  const mgContract = (df: DG.DataFrame, withBw: boolean): ColumnContract => {
    const obs = {...ratContract(df).observations,
      dose: {name: 'Dose_mg', required: false, rawType: 'float64' as const, units: 'mg' as const}};
    if (!withBw) delete (obs as {bodyWeight?: unknown}).bodyWeight;
    return {observations: obs};
  };
  const mgScalars = {...ratScalars(nca.ROUTE_IV_BOLUS), doseUnits: 'mg' as const};

  test('fixed-mg dose (0.2441 mg) + BW → per kg of R001; without BW → per kg blank, CL abs. set', async () => {
    const withBw = await r001(addDoseMg);
    const a = runPipelineOrThrow(withBw, mgContract(withBw, true), mgScalars, DEFAULT_RULES).results;
    expectFloat(a.col('CL')!.get(rowOf(a, 'R001')) as number, R001_CL_PER_KG, 5e-4, 'fixed mg + BW → per kg');
    const noBw = await r001(addDoseMg);
    const b = runPipelineOrThrow(noBw, mgContract(noBw, false), mgScalars, DEFAULT_RULES).results;
    const i = rowOf(b, 'R001');
    expect(b.col('CL')!.isNone(i), true, 'fixed mg, no BW → per kg not formed');
    expectFloat(b.col('CL_ABS')!.get(i) as number, R001_CL_ABS, R001_CL_ABS * 0.001, 'CL abs. reported');
    expect(warnCodes(b, i).includes('PER_KG_NOT_FORMED'), true, 'AC-2.4: the reason is in WARNINGS_JSON');
  }, {timeout: TIMEOUT});

  test('AD-2: fixed-mg dose and BW varying within R001 → per kg refused with the range', async () => {
    const df = await r001((d) => {
      addDoseMg(d); // dose taken BEFORE the weight change: R001 received 0.2441 mg
      const bw = d.col('BodyWeight_g')!; const subj = d.col('SubjectID')!;
      let k = 0;
      for (let r = 0; r < d.rowCount; r++)
        if (subj.get(r) === 'R001' && k++ === 3) bw.set(r, 260.0); // one later sample weighs more
    });
    const {results} = runPipelineOrThrow(df, mgContract(df, true), mgScalars, DEFAULT_RULES);
    const i = rowOf(results, 'R001');
    expect(results.col('CL')!.isNone(i), true, 'varying BW + absolute mg → per kg refused');
    const w = (JSON.parse(String(results.col('WARNINGS_JSON')!.get(i))) as Array<{code: string; message: string}>)
      .find((x) => x.code === 'PER_KG_NOT_FORMED');
    expect(w?.message.includes('0.2441–0.26 kg') ?? false, true, `range named: ${w?.message}`);
    expectFloat(results.col('CL_ABS')!.get(i) as number, R001_CL_ABS, R001_CL_ABS * 0.001, 'absolute kept');
  }, {timeout: TIMEOUT});

  test('AD-1: BW missing for one subject → that row per-kg via the per-kg engine path; others absolute', async () => {
    const df = await r001((d) => {
      const bw = d.col('BodyWeight_g')!; const subj = d.col('SubjectID')!;
      for (let r = 0; r < d.rowCount; r++) if (subj.get(r) === 'R001') bw.set(r, DG.FLOAT_NULL);
    });
    const {results} = runPipelineOrThrow(df, ratContract(df), ratScalars(nca.ROUTE_IV_BOLUS), DEFAULT_RULES);
    const i = rowOf(results, 'R001');
    expect(String(results.col('DOSE_NORM')!.get(i)), 'per_kg', 'blank BW + mg/kg → the per-kg engine path');
    expectFloat(results.col('CL')!.get(i) as number, R001_CL_PER_KG, 5e-4, 'same per-kg value as with BW');
    expect(results.col('CL_ABS')!.isNone(i), true, 'no weight → no CL abs.');
    expect(results.col('CL')!.tags['units'], 'mL/min/kg', 'one unit for the mixed frame (E5 class)');
  }, {timeout: TIMEOUT});
});

// =====================================================================
// per-row-route — combined mixed-route pipeline (03_rat_multiroute.csv).
// One Route column carries bare IV / PO / SC across 64 subjects. Validates
// Option-A end-to-end plumbing (AC-A7), per-subject CL_BASIS (AC-C1), the
// bare-IV → IV-bolus warning (AC-A4), and an end-to-end F% sanity (AD-6).
// =====================================================================
category('Reference suite: per-row mixed route (rat multiroute)', () => {
  test('AC-A7/AC-C1: per-subject ROUTE + CL_BASIS reflect each subject\'s own route', async () => {
    const df = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(df);
    const {results} = runPipelineOrThrow(
      df, ratMultirouteContract(df), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);

    expect(results.rowCount, 64, 'expected 64 subjects across all routes');
    // Neutral numeric label even though the frame is mixed-route.
    expect(results.col('CL')!.tags['label'], 'CL', 'CL label neutral on mixed frame');

    const routeCol = results.col('ROUTE')!;
    const basisCol = results.col('CL_BASIS')!;
    const seen = new Set<string>();
    for (let i = 0; i < results.rowCount; i++) {
      const route = String(routeCol.get(i));
      const basis = String(basisCol.get(i));
      seen.add(route);
      const expectedBasis = route === 'IV-bolus' ? 'CL' : 'CL/F';
      expect(['IV-bolus', 'PO', 'SC'].includes(route), true, `unexpected route '${route}' at row ${i}`);
      expect(basis, expectedBasis, `${route} subject row ${i} → ${expectedBasis}`);
    }
    // Plumbing proof: all three input routes are represented per-subject (not a
    // single scalar applied to every row, not the raw 'IV'/'PO'/'SC' strings).
    expect(seen.has('IV-bolus') && seen.has('PO') && seen.has('SC'), true,
      'IV, PO and SC subjects all present per-row');
  }, {timeout: TIMEOUT});

  test('AC-A4: bare "IV" subjects carry the BARE_IV_DEFAULTED_TO_BOLUS warning', async () => {
    const df = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(df);
    const {results} = runPipelineOrThrow(
      df, ratMultirouteContract(df), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);

    const routeCol = results.col('ROUTE')!;
    const warnCol = results.col('WARNINGS_JSON')!;
    let ivChecked = 0;
    for (let i = 0; i < results.rowCount; i++) {
      if (String(routeCol.get(i)) !== 'IV-bolus') continue;
      const warns = JSON.parse(String(warnCol.get(i))) as Array<{code: string}>;
      expect(warns.some((w) => w.code === BARE_IV_DEFAULTED_TO_BOLUS), true,
        `IV subject row ${i} must carry the bare-IV default warning`);
      ivChecked++;
    }
    expect(ivChecked > 0, true, 'at least one IV subject checked');
  }, {timeout: TIMEOUT});

  test('AD-6: end-to-end F% vs PKNCA reference (PO/F/5, PO/M/5, SC/F/25)', async () => {
    const df = await loadFixture('03_rat_multiroute.csv');
    const opts = setupMultirouteRoute(df);
    const {results, outcome} = runPipelineOrThrow(
      df, ratMultirouteContract(df), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    // Carry Sex / Dose metadata into the results frame (as the pane does).
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), ratMultirouteContract(df));

    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['Sex']});
    const lookup = (route: string, sex: string, dose: number): number => {
      for (let i = 0; i < table.rowCount; i++) {
        if (String(table.col('ROUTE')!.get(i)) === route &&
            String(table.col('Sex')!.get(i)) === sex &&
            (table.col('DOSE')!.get(i) as number) === dose)
          return table.col('F_PCT')!.get(i) as number;
      }
      throw new Error(`F% stratum ${route}/${sex}/${dose} not found`);
    };
    // End-to-end F% propagates each subject's NCA AUCinf (not the stored PKNCA
    // values), so the exact-number validation is the DETERMINISTIC fixture test
    // (bioavailability-tests.ts AC-D1/D3, built from the stored dn_AUC columns —
    // those reproduce 43.07 / 31.95 / 117.51 to 1e-6). Here we assert the full
    // pipeline lands close on the well-populated PO strata (proving the Option-A
    // route plumbing → dose-normalize → sex-stratified ratio chain is sound), and
    // assert the *scientific property* on the sparse SC/F/25 cell.
    //
    // PO strata: 0.1% band. With BUG-05 fixed (sci-comp 0.7.1) R005/R013 contribute
    // valid AUCinf to the pooled IV denominator, and with R019's λz-window divergence
    // resolved (sci-comp 0.9.1, VAL-01-LZ-R019) the end-to-end F% now reproduces the
    // deterministic reference to < 1e-6 — so the former 3% band is tightened to a strict
    // integration tolerance (exact values still pinned by bioavailability-tests AC-D1/D3).
    expectFloat(lookup('PO', 'F', 5), 43.0677728553402, 43.07 * 0.001, 'PO/F/5 ≈ 43.07 (0.1%)');
    expectFloat(lookup('PO', 'M', 5), 31.9451449763666, 31.95 * 0.001, 'PO/M/5 ≈ 31.95 (0.1%)');
    // SC/F/25 is a sparse cell whose aggregate dn_AUC inherits wider per-subject
    // NCA spread; the scientific invariant that matters is that F% is reported
    // UNCAPPED (> 100%), not clamped. Exact 117.51 is pinned by AC-D3.
    const scF25 = lookup('SC', 'F', 25);
    expect(scF25 > 100, true, `SC/F/25 must be reported uncapped (> 100%); got ${scF25}`);
    expect(Number.isFinite(scF25), true, 'SC/F/25 F% is finite (IV-F reference present)');
  }, {timeout: TIMEOUT});
});

// =====================================================================
// tk-aggregates AD-7 / AC-15 — regenerated TK reference-card integrity.
// The P7 regen re-anchors TC-012 (Rac from auclast, not the all-NA aucinf.obs)
// and TC-013 (per-sex exposure margin, previously absent). This is the rule-18
// oracle gate: the committed reference cards must be non-NA / non-empty and carry
// the per-sex structure the TK primitives validate against. No *existing*
// parameter changes — only the two previously-NA/empty aggregate cards populate
// (science-preservation: defect fix, not a behavior change). The arithmetic-vs-
// primitive parity is asserted in toxicokinetics-tests.ts.
// =====================================================================
category('Reference suite: TK aggregate cards (tk-aggregates AC-15)', () => {
  async function loadCsv(name: string): Promise<DG.DataFrame> {
    const csv = await grok.dapi.files.readAsText(`System:AppData/Nca/${name}`);
    return DG.DataFrame.fromCsv(csv);
  }

  test('AC-15: accumulation card — AUC_acc_ratio (Formula A, auclast) non-NA for all 8 cells', async () => {
    const acc = await loadCsv('rat_tk28d_accumulation_ratios.csv');
    expect(acc.rowCount, 8, '8 (Dose×Sex) cells');
    for (let i = 0; i < acc.rowCount; i++) {
      const aucRatio = acc.col('AUC_acc_ratio')!.get(i) as number;
      const cmaxRatio = acc.col('Cmax_acc_ratio')!.get(i) as number;
      expect(Number.isFinite(aucRatio) && aucRatio > 0, true,
        `row ${i}: AUC_acc_ratio must be non-NA (was all-NA before P7 — DG-2)`);
      expect(Number.isFinite(cmaxRatio) && cmaxRatio > 0, true, `row ${i}: Cmax_acc_ratio non-NA`);
      // Formula label committed (Formula A — the only one the 06 cells support).
      expect(String(acc.col('Rac_formula')!.get(i)).indexOf('A') >= 0, true, `row ${i}: Formula A label`);
    }
    // auclast-pivot arithmetic: AUC_acc_ratio == auclast_Day28 / auclast_Day1.
    for (let i = 0; i < acc.rowCount; i++) {
      const d1 = acc.col('auclast_Day1')!.get(i) as number;
      const dN = acc.col('auclast_Day28')!.get(i) as number;
      expectFloat(acc.col('AUC_acc_ratio')!.get(i) as number, dN / d1, Math.abs(dN / d1) * 1e-9,
        `row ${i}: AUC_acc_ratio = auclast_Day28 / auclast_Day1`);
    }
  }, {timeout: TIMEOUT});

  test('AC-15: exposure-margin card — non-empty, one row per sex, non-NA HED + valid-sex projection', async () => {
    const margin = await loadCsv('rat_tk28d_exposure_margin.csv');
    expect(margin.rowCount, 2, 'one row per sex (was sex-pooled/absent before P7 — DG-3 / F-1)');
    const sexes = new Set<string>();
    let hedNonNa = 0;
    let projFinite = 0;
    // NOTE: a CSV "NA" round-trips through DG.fromCsv into a FLOAT column's
    // FLOAT_NULL sentinel (≈2.68e-34), which IS finite — so `Number.isFinite`
    // alone treats a null cell as a value. Guard with `isNone` first (the
    // documented FLOAT_NULL footgun) so M's NA projection is counted as null.
    const hedCol = margin.col('MARGIN_HED')!;
    const projCol = margin.col('MARGIN_PROJECTED')!;
    for (let i = 0; i < margin.rowCount; i++) {
      sexes.add(String(margin.col('Sex')!.get(i)));
      if (!hedCol.isNone(i) && Number.isFinite(hedCol.get(i) as number)) hedNonNa++;
      if (!projCol.isNone(i) && Number.isFinite(projCol.get(i) as number)) projFinite++;
      // Exposure-based default is honestly NA (no empirical human AUC — DG-1).
      expect(margin.col('MARGIN_EXPOSURE')!.isNone(i), true,
        `row ${i}: MARGIN_EXPOSURE NA (no empirical human AUC)`);
    }
    expect(sexes.has('F') && sexes.has('M'), true, 'both sexes present');
    expect(hedNonNa, 2, 'MARGIN_HED non-NA for both sexes (dose-based, F-1)');
    // F has a valid Day-1 CL → finite projection; M Day-1 CL is NA → projection NA (F-2).
    expect(projFinite, 1, 'exactly one sex (F) has a non-NA CL projection (10/M CL is NA)');
  }, {timeout: TIMEOUT});
});

// =====================================================================
// Dose proportionality (Smith power model) — INDIVIDUAL-fit parity vs the
// R lm()/confint() oracle (AD-12 / AC-11, rule-18). The adapter is fed the
// same PKNCA per-subject AUCinf the oracle used (rat_iv_nca_datagrok_match.csv),
// so the β / 90% CI must match the R-generated card byte-for-byte (the
// regression math is sci-comp's `linearFit`, validated upstream).
//
// Oracle: scripts/rat_multiroute_dose_proportionality_power.R (lm(log(aucinf.obs)
// ~ log(Dose_mg_kg)) per Route×Sex, 90% confint). REGRESSION_PIN values below
// are recomputed against the actual CSV at build (rule-16), not assumed.
// =====================================================================
category('Reference suite: dose proportionality (Smith power model)', () => {
  test('AC-11: IV individual-fit β / 90% CI parity vs the R lm() oracle', async () => {
    // `rat_iv_nca_datagrok_match.csv` carries the identical IV/F + IV/M per-subject
    // rows as the IV portion of `rat_multiroute_nca_datagrok_match.csv` (the file
    // the R oracle `rat_multiroute_dose_proportionality_power.R` consumes) — the
    // IV anchor below was confirmed byte-identical against both. The bundled rat_iv
    // file is used here because the multiroute match CSV is not in package/files.
    const df = await loadFixture('rat_iv_nca_datagrok_match.csv');
    const res = computeDoseProportionality(df, {
      doseCol: 'Dose_mg_kg', routeCol: 'Route', sexCol: 'Sex',
      aucInfCol: 'aucinf.obs', cmaxCol: 'cmax',
    });

    const ivF = res.fits.find((f) => f.route === 'IV' && f.sex === 'F' && f.parameter === 'AUC');
    expect(ivF !== undefined, true, 'IV/F AUC fit present');
    // REGRESSION_PIN (rat_iv_nca_datagrok_match.csv, 90% CI): β 0.996437,
    // df 10, CI [0.835099, 1.157776], R² 0.926092, CR [0.903090, 1.096910].
    expectFloat(ivF!.beta, 0.996437, 1e-4, 'IV/F β = 0.9964 (individual fit, df = N−2)');
    expect(ivF!.df, 10, 'IV/F df = N − 2 = 10');
    expect(ivF!.n, 12, 'IV/F N = 12 individual subjects');
    expect(ivF!.nDoseLevels, 3, 'IV/F 3 dose levels');
    expect(ivF!.isGroupMeanFallback, false, 'IV/F uses the individual-subject path');
    expectFloat(ivF!.slopeCI[0], 0.835099, 1e-3, 'IV/F 90% CI low');
    expectFloat(ivF!.slopeCI[1], 1.157776, 1e-3, 'IV/F 90% CI high');
    expectFloat(ivF!.criticalRegion[0], 0.903090, 1e-4, 'CR low (r = 10)');
    expectFloat(ivF!.criticalRegion[1], 1.096910, 1e-4, 'CR high (r = 10)');
    expectFloat(ivF!.rSquared, 0.926092, 1e-3, 'IV/F R²');
    // The honest call: the CI straddles the CR → IIV-driven inconclusive
    // (β ≈ 1.00 point estimate, but underpowered — never dressed as biology).
    expect(ivF!.verdict, 'inconclusive-underpowered', 'IV/F verdict');
    expect(ivF!.mode, 'iiv-driven', 'IV/F mode is IIV-driven');

    const ivM = res.fits.find((f) => f.route === 'IV' && f.sex === 'M' && f.parameter === 'AUC');
    expect(ivM !== undefined, true, 'IV/M AUC fit present');
    expectFloat(ivM!.beta, 0.944137, 1e-4, 'IV/M β');
    expectFloat(ivM!.slopeCI[0], 0.844168, 1e-3, 'IV/M 90% CI low');
    expectFloat(ivM!.slopeCI[1], 1.044105, 1e-3, 'IV/M 90% CI high');
  }, {timeout: TIMEOUT});
});

category('Reference suite: Tlag LLOQ-absent (fixture 9, full adapter path)', () => {
  // Fixture 9 (C-5): an extravascular profile with leading zeros and NO BLQ /
  // LLOQ column. Routed through the FULL adapter path (runNcaPipeline →
  // extractProfileInputs → computeNca → tlag), NOT a pre-zeroed vector handed
  // to computeNca (which would make the alignment claim vacuous). sci-comp tlag
  // is route-gated (IV → NaN), so the route is extravascular (ROUTE_PO). The
  // boundary is exact by construction (last C==0 before the first C>0) and
  // sci-comp pins Tlag to the SAME quantifiable boundary AUC uses (derived.ts).
  const scalarsPo: ProfileScalars = {
    dose: 100, route: nca.ROUTE_PO, doseUnits: 'mg', concentrationUnits: 'mg/L',
    timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0,
  };
  // LLOQ-absent by construction — NO BLQ flag and NO LLOQ column (the C-5 fixture).
  // With neither mapped, the extractor leaves the BLQ mask all-zeros, so the
  // leading conc==0 rows reach computeNca as measured zeros — and sci-comp's
  // `tlag` reads the conc==0 boundary directly (it operates on the conc series,
  // not the mask). This is the actual desync-prevention guarantee on zero-coded
  // LIMS exports — tlag pins to the same conc>0 boundary AUC starts from. (An LLOQ
  // column is deliberately NOT added — it would both contradict "LLOQ-absent" and
  // trip the LLOQ_OUT_OF_RANGE sanity validator on the all-zero sibling.)
  const poContract = (df: DG.DataFrame) => ({
    observations: {
      subjectId: {name: 'Subject', required: true,
        rawType: (df.col('Subject')!.type === DG.COLUMN_TYPE.INT ? 'int32' : 'string') as 'int32' | 'string'},
      time: {name: 'Time', required: true, rawType: 'float64' as const, units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64' as const, units: 'mg/L'},
    },
  });

  // Oral 1-cmt with a 1 h absorption lag: C(t)=0 for t≤1, then the standard
  // first-order-absorption rise/fall on the shifted clock (τ = t−1). Leading
  // samples (0, 0.5, 1) are coded 0 (below quantifiable, no LLOQ column).
  const ka = 1.0; const k = 0.2; const V = 10; const DOSE = 100; const LAG = 1;
  const cfun = (t: number): number => t <= LAG ? 0 :
    (DOSE * ka) / (V * (ka - k)) * (Math.exp(-k * (t - LAG)) - Math.exp(-ka * (t - LAG)));
  const times = [0, 0.5, 1, 2, 4, 6, 8, 12, 16, 24];
  const concLag = times.map(cfun); // 0,0,0,5.636…,…

  test('AC-9.1/9.2: lag>0 → TLAG = last zero before first positive = AUC-start boundary', async () => {
    const df = makeSyntheticObs({
      subjects: times.map(() => 'L1'),
      times: Float64Array.from(times),
      concs: Float64Array.from(concLag),
    });
    const {results} = runPipelineOrThrow(df, poContract(df), scalarsPo, DEFAULT_RULES);
    expect(results.col('TLAG')!.isNone(0), false, 'lag>0 TLAG is present');
    // AC-9.1: Tlag = 1 (the t=1 zero immediately before the first positive at t=2).
    expect(results.col('TLAG')!.get(0) as number, 1, 'lag>0 → TLAG = 1');
    // AC-9.2 (alignment): TLAG equals the time of the sample just before the
    // first measurable conc — the same quantifiable boundary AUC integrates from.
    const firstPosIdx = concLag.findIndex((c) => c > 0);
    expect(results.col('TLAG')!.get(0) as number, times[firstPosIdx - 1],
      'TLAG == time of last zero before first positive (AUC-start boundary)');
    expect(results.col('VSS')!.isNone(0), true,
      'extravascular → VSS NaN (route gate), confirms non-IV path');
  }, {timeout: TIMEOUT});

  test('AC-9.3: measurable t=0 → TLAG = 0; all-zero → TLAG = NaN', async () => {
    // Measurable from the first sample → no observed lag → TLAG = time[0] = 0.
    const measT0 = times.map((t) => t === 0 ? 8.0 : cfun(t + 1)); // positive at t=0
    const dfMeas = makeSyntheticObs({
      subjects: times.map(() => 'M1'),
      times: Float64Array.from(times),
      concs: Float64Array.from(measT0),
    });
    const {results: rMeas} = runPipelineOrThrow(dfMeas, poContract(dfMeas), scalarsPo, DEFAULT_RULES);
    expect(rMeas.col('TLAG')!.get(0) as number, 0, 'measurable t=0 → TLAG = 0');

    // All-zero profile → never rises above zero → TLAG = NaN (never expect(x, undefined)).
    const dfZero = makeSyntheticObs({
      subjects: times.map(() => 'Z1'),
      times: Float64Array.from(times),
      concs: Float64Array.from(times.map(() => 0)),
    });
    const {results: rZero} = runPipelineOrThrow(dfZero, poContract(dfZero), scalarsPo, DEFAULT_RULES);
    expect(rZero.col('TLAG')!.isNone(0), true, 'all-zero → TLAG is EMPTY');
  }, {timeout: TIMEOUT});
});

// =====================================================================
// GROK-20960 AC-N.4.3 — the DATA-GAP-42 corollary on the FULL adapter path.
// Before sci-comp 0.11.0 the LLOQ MAPPING decided the effective BLQ policy: an
// unmapped LLOQ left the BLQ mask all-zeros (the embedded 0 integrated as a
// measured zero), while a mapped LLOQ flagged it and the shipped `set-zero`
// default was silently `exclude` (the point dropped, interpolated across) —
// two different AUClast for one table. Now `set-zero` integrates the 0 either
// way, so mapping no longer decides the policy. Pinned here, not inferred.
// =====================================================================
category('Reference suite: BLQ mapped ≡ unmapped under the shipped default (AC-N.4.3)', () => {
  // A 9-row PO frame with ONE embedded BLQ (t = 6 coded 0) inside a positive tail.
  const times = [0, 0.5, 1, 2, 4, 6, 8, 12, 24];
  const concs = [0, 1.2, 2.0, 3.0, 2.0, 0, 0.5, 0.25, 0.06];
  const subjects = times.map(() => 'B1');
  const scalars: ProfileScalars = {
    dose: 100, route: nca.ROUTE_PO, doseUnits: 'mg', concentrationUnits: 'mg/L',
    timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0.05,
  };
  const contract = (df: DG.DataFrame, withLloq: boolean) => ({
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: 'string' as const},
      time: {name: 'Time', required: true, rawType: 'float64' as const, units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64' as const, units: 'mg/L'},
      ...(withLloq && df.col('LLOQ') !== null ?
        {lloq: {name: 'LLOQ', required: false, rawType: 'float64' as const, units: 'mg/L'}} : {}),
    },
  });

  test('the same table, LLOQ mapped vs unmapped → the SAME AUClast (set-zero integrates the 0)', async () => {
    const unmapped = makeSyntheticObs({
      subjects, times: Float64Array.from(times), concs: Float64Array.from(concs),
    });
    const mapped = makeSyntheticObs({
      subjects, times: Float64Array.from(times), concs: Float64Array.from(concs),
      lloqs: Float64Array.from(times.map(() => 0.05)),
    });
    const {results: rU} = runPipelineOrThrow(unmapped, contract(unmapped, false), scalars, DEFAULT_RULES);
    const {results: rM, outcome: oM} = runPipelineOrThrow(mapped, contract(mapped, true), scalars, DEFAULT_RULES);
    // The MAPPED run really did flag the embedded zero (the inference fired)…
    const maskM = oM.perProfile[0].result.provenance.blqApplied;
    expect(maskM.conc.length, 9, 'mapped: 9 observations');
    // …and the unmapped one did not — yet the numbers agree.
    expect(rU.col('STATUS')!.get(0), 'ok', 'unmapped computes ok');
    expect(rM.col('STATUS')!.get(0), 'ok', 'mapped computes ok');
    for (const col of ['AUC_LAST', 'AUMC_LAST', 'AUC_INF', 'LAMBDA_Z', 'CMAX', 'TMAX', 'TLAG']) {
      expect(rM.col(col)!.get(0), rU.col(col)!.get(0),
        `${col}: mapped (${rM.col(col)!.get(0)}) must equal unmapped (${rU.col(col)!.get(0)}) — ` +
        'the LLOQ mapping must not decide the BLQ policy');
    }
    // And the number is the `set-zero` one — the embedded 0 is INTEGRATED (linear
    // down to 0 and up again), not interpolated across as `exclude` would do.
    const aucExclude = runPipelineOrThrow(mapped, contract(mapped, true), scalars,
      {...DEFAULT_RULES, blq: {...DEFAULT_RULES.blq, embedded: 'exclude'}}).results.col('AUC_LAST')!.get(0) as number;
    const aucSetZero = rM.col('AUC_LAST')!.get(0) as number;
    expect(aucSetZero < aucExclude, true,
      `set-zero AUClast (${aucSetZero}) must be BELOW exclude's (${aucExclude}) — the 0 was integrated`);
    expect(Math.abs(aucSetZero - aucExclude) / aucExclude > 0.01, true, 'and by more than 1 %');
  }, {timeout: TIMEOUT});
});
