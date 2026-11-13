/**
 * Sex-stratifier designation — Slice A science-preservation gates.
 *
 * @see files/03_rat_multiroute_gender.csv (Sex→Gender rename; data byte-identical)
 *
 * The feature lets a scientist DESIGNATE the sex column (any name) instead of
 * binding to one hardcoded `/^sex$/i`. The science-preservation guarantee is
 * that designating a renamed column reproduces the EXACT per-sex F% the
 * `Sex`-named column produces — the generator (`gen_sex_stratifier_fixtures.py`)
 * rewrites ONLY the Sex header, so every numeric cell is identical and the two
 * pipelines must land bit-for-bit equal per stratum (AC-A1).
 *
 * AC-A5 (safe degradation): a sex column that is NOT constant-per-subject is not
 * carried by `appendGroupingColumns`, so the pane resolves the stratifier to
 * null and F% stays absent (never a garbled partial stratification). This test
 * pins the carriage mechanism the pane's `carried.includes(designated)` gate
 * relies on — verified against real generated output (rule 16).
 *
 * **Footgun gate:** `expect(cond === expected, true, msg)` form throughout.
 *
 * AC-A2 (suggestion-not-bind) is observable UX (the static hint) verified at the
 * live `grok test` / visual-review gate, not here.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, ratMultirouteContract, setupMultirouteRoute, ratScalars, runPipelineOrThrow,
} from './test-utils';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {computeBioavailability} from '../adapter/bioavailability';
import {applySexMap, type SexCode} from '../adapter/sex-mapping';
import {
  normalizeEventDrivenObservations, CONC_COL, TAD_COL, SUBJECT_COL,
  type EventContract,
} from '../adapter/event-driven';
import {poolCompositeByCell} from '../adapter/composite-pool';

const TIMEOUT = 60000;

/** Run the full mixed-route pipeline on a fixture and return its sex-stratified
 *  F% table, stratifying on `sexCol` (the DESIGNATED column name). Mirrors the
 *  pane's chain: route-encode → pipeline → appendGroupingColumns → F%. */
async function fpctTable(
  csv: string, sexCol: string, map?: ReadonlyMap<string, SexCode>,
): Promise<DG.DataFrame> {
  const df = await loadFixture(csv);
  // Slice B: apply the value→M/F/U map to the sex column BEFORE the pipeline +
  // carriage, exactly as the pane does on the source frame.
  if (map !== undefined) applySexMap(df, sexCol, map);
  const opts = setupMultirouteRoute(df);
  const {results, outcome} = runPipelineOrThrow(
    df, ratMultirouteContract(df), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
  const carried = appendGroupingColumns(
    results, df, outcome.perProfile.map((p) => p.index), ratMultirouteContract(df));
  expect(carried.includes(sexCol), true, `designated sex column "${sexCol}" must be carried`);
  return computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: [sexCol]}).table;
}

/** Index an F% table by `route|sex|dose` → F_PCT for stratum-wise comparison. */
function indexFpct(table: DG.DataFrame, sexCol: string): Map<string, number> {
  const out = new Map<string, number>();
  for (let i = 0; i < table.rowCount; i++) {
    const route = String(table.col('ROUTE')!.get(i));
    const sex = String(table.col(sexCol)!.get(i));
    const dose = table.col('DOSE')!.get(i) as number;
    out.set(`${route}|${sex}|${dose}`, table.col('F_PCT')!.get(i) as number);
  }
  return out;
}

category('Sex stratifier: designation (Slice A)', () => {
  test('AC-A1: Gender-renamed designation reproduces the Sex per-stratum F% oracle', async () => {
    // Baseline: the conventionally-named `Sex` column (the prior hardcoded path).
    const sexTable = indexFpct(await fpctTable('03_rat_multiroute.csv', 'Sex'), 'Sex');
    // Designated: the SAME study with the header renamed `Sex` → `Gender`.
    const genderTable = indexFpct(await fpctTable('03_rat_multiroute_gender.csv', 'Gender'), 'Gender');

    expect(genderTable.size, sexTable.size, 'same stratum count');
    expect(genderTable.size > 0, true, 'F% strata produced');
    // Byte-identical data ⇒ each stratum's F% must match exactly (science gate).
    for (const [key, sexVal] of sexTable.entries()) {
      const genderVal = genderTable.get(key);
      expect(genderVal !== undefined, true, `stratum ${key} present in Gender run`);
      expectFloat(genderVal!, sexVal, 1e-9, `F% stratum ${key} reproduced via designation`);
    }
    // Pin to the ABSOLUTE PKNCA oracle on the well-populated PO strata (matching
    // reference-suite AD-6's 0.1% band) so the Gender path is verified against the
    // reference, not merely equal to the Sex run. SC/F/25 keeps the uncapped-only
    // assertion — the full pipeline's sparse SC/F/25 lands wide of 117.51 (its
    // exact value is pinned by the deterministic bioavailability-tests AC-D3, not
    // the end-to-end path), so AD-6 asserts only the >100 scientific property.
    // (Band tightened 3%→0.1% with VAL-01-LZ-R019 resolved in sci-comp 0.9.1.)
    expectFloat(genderTable.get('PO|F|5')!, 43.0677728553402, 43.07 * 0.001,
      'PO/F/5 ≈ 43.07 via Gender designation (0.1%)');
    expectFloat(genderTable.get('PO|M|5')!, 31.9451449763666, 31.95 * 0.001,
      'PO/M/5 ≈ 31.95 via Gender designation (0.1%)');
    const scF25 = genderTable.get('SC|F|25')!;
    expect(scF25 > 100, true, `SC/F/25 reported uncapped (> 100%); got ${scF25}`);
  }, {timeout: TIMEOUT});

  test('AC-A5: a non-constant-per-subject designation is not carried → F% degrades safely', async () => {
    const df = await loadFixture('03_rat_multiroute.csv');
    // Make `Sex` vary WITHIN one subject (a mis-designated, non-per-subject
    // column). The constancy gate must then refuse to carry it.
    const sexSrc = df.col('Sex')!;
    const subj = df.col('SubjectID')!;
    const firstSubj = subj.getString(0);
    let flipped = 0;
    for (let i = 0; i < df.rowCount && flipped < 1; i++) {
      if (subj.getString(i) === firstSubj && i > 0) {
        sexSrc.set(i, sexSrc.getString(i) === 'M' ? 'F' : 'M');
        flipped++;
      }
    }
    expect(flipped, 1, 'flipped one row of the first subject to break constancy');

    const opts = setupMultirouteRoute(df);
    const {results, outcome} = runPipelineOrThrow(
      df, ratMultirouteContract(df), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    const carried = appendGroupingColumns(
      results, df, outcome.perProfile.map((p) => p.index), ratMultirouteContract(df));
    // The gate the pane relies on: a non-constant designation is NOT carried, so
    // `carried.includes('Sex')` is false → sexForDecorators === null → F% absent
    // + the needs-Sex banner fires (safe degradation, never a partial stratum).
    expect(carried.includes('Sex'), false,
      'a non-constant-per-subject Sex column must not be carried (safe degradation)');
  }, {timeout: TIMEOUT});
});

category('Sex stratifier: value→M/F/U mapping (Slice B)', () => {
  test('AC-B3: token-collision map (M/Male/MALE→M, F/Female→F) reproduces the clean F% oracle', async () => {
    // Baseline: the clean `Sex` = {M,F} study.
    const baseline = indexFpct(await fpctTable('03_rat_multiroute.csv', 'Sex'), 'Sex');
    // Collision: the SAME study with males in {M,Male,MALE}, females in {F,Female};
    // the explicit map collapses them back to M/F. Byte-identical data ⇒ the F%
    // must equal the baseline exactly (the explicit-merge science guarantee).
    const map = new Map<string, SexCode>([
      ['M', 'M'], ['Male', 'M'], ['MALE', 'M'], ['F', 'F'], ['Female', 'F'],
    ]);
    const mapped = indexFpct(await fpctTable('03_rat_multiroute_sexcollision.csv', 'Sex', map), 'Sex');

    expect(mapped.size, baseline.size, 'same stratum count after collapse');
    expect(mapped.size > 0, true, 'F% strata produced');
    for (const [key, baseVal] of baseline.entries()) {
      const mapVal = mapped.get(key);
      expect(mapVal !== undefined, true, `stratum ${key} present after collapse`);
      expectFloat(mapVal!, baseVal, 1e-9, `F% stratum ${key} unchanged by explicit token merge`);
    }
  }, {timeout: TIMEOUT});

  test('AC-B4: coded-1/2 map (1→M, 2→F) before pooling reproduces the M/F composite cells', async () => {
    // Why this is safe: applying the documented convention (1=M, 2=F) to
    // the source sex column BEFORE `poolCompositeByCell` is a LABEL-ONLY relabel of
    // the cell key — every pooled concentration is unchanged vs the M/F fixture.
    const ec: EventContract = {
      id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE',
      mdv: 'MDV', addl: 'ADDL', ii: 'II', route: 'ROUTE', bodyWeight: 'BW_kg', day: 'DAY',
      covariates: ['Sex', 'DOSE_MGKG'],
    };

    // M/F reference frame → production pooling.
    const mfDf = await loadFixture('06_rat_tk28d_nonmem.csv');
    const mfNorm = normalizeEventDrivenObservations(mfDf, ec);
    expect(mfNorm.errors.length, 0, 'M/F: no blocking normalization errors');
    const mfPooled = poolCompositeByCell(mfNorm.observations, {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex'});

    // Coded-1/2 frame → apply the documented map on the normalized source BEFORE
    // pooling (exactly the pane's application point) → production pooling.
    const codedDf = await loadFixture('06_rat_tk28d_nonmem_sexcoded.csv');
    const codedNorm = normalizeEventDrivenObservations(codedDf, ec);
    expect(codedNorm.errors.length, 0, 'coded: no blocking normalization errors');
    applySexMap(codedNorm.observations, 'Sex', new Map<string, SexCode>([['1', 'M'], ['2', 'F']]));
    const codedPooled = poolCompositeByCell(codedNorm.observations, {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex'});

    expect(codedPooled.cellCount, mfPooled.cellCount, 'same composite cell count');
    // Index each pooled frame by cell|TAD → mean concentration.
    const indexPooled = (p: DG.DataFrame): Map<string, number> => {
      const out = new Map<string, number>();
      const subj = p.col(SUBJECT_COL)!; const tad = p.col(TAD_COL)!; const conc = p.col(CONC_COL)!;
      for (let i = 0; i < p.rowCount; i++)
        out.set(`${subj.getString(i)}|${tad.get(i) as number}`, conc.get(i) as number);
      return out;
    };
    const mfMap = indexPooled(mfPooled.pooled);
    const codedMap = indexPooled(codedPooled.pooled);
    expect(codedMap.size, mfMap.size, 'same (cell × TAD) point count');
    let checked = 0;
    for (const [key, mfConc] of mfMap.entries()) {
      const codedConc = codedMap.get(key);
      expect(codedConc !== undefined, true, `point ${key} present in coded-then-mapped pooling`);
      expectFloat(codedConc!, mfConc, 1e-9, `composite conc at ${key} unchanged by 1/2→M/F relabel`);
      checked++;
    }
    expect(checked > 0, true, 'at least one composite point compared');
    // Confirm the relabel actually happened (keys are M/F, not 1/2).
    expect([...codedMap.keys()].some((k) => k.includes('_M_') || k.includes('_F_')), true,
      'coded cell keys relabelled to M/F');
    expect([...codedMap.keys()].some((k) => k.includes('_1_') || k.includes('_2_')), false,
      'no raw 1/2 survives in the composite cell key');
  }, {timeout: TIMEOUT});

  test('AC-B5: single-sex study maps to one stratum and computes (no error)', async () => {
    // A 1-value sex column is valid, not an error: every subject maps to one
    // stratum. Force the multiroute study to a single sex and confirm F% still
    // computes (one IV reference, EV strata divided by it).
    const df = await loadFixture('03_rat_multiroute.csv');
    applySexMap(df, 'Sex', new Map<string, SexCode>([['M', 'U'], ['F', 'U']])); // collapse to one stratum
    const opts = setupMultirouteRoute(df);
    const {results, outcome} = runPipelineOrThrow(
      df, ratMultirouteContract(df), ratScalars(nca.ROUTE_PO), DEFAULT_RULES, opts);
    appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), ratMultirouteContract(df));
    const {table} = computeBioavailability(results, {doseCol: 'Dose_mg_kg', stratifyCols: ['Sex']});
    expect(table.rowCount > 0, true, 'single-sex F% still produces rows');
    const sexVals = new Set<string>();
    for (let i = 0; i < table.rowCount; i++) sexVals.add(String(table.col('Sex')!.get(i)));
    expect(sexVals.size, 1, 'exactly one sex stratum');
    expect(sexVals.has('U'), true, 'the single stratum is the mapped value U');
  }, {timeout: TIMEOUT});
});
