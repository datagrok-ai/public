/**
 * Adapter — column autodetect category.
 *
 * @see src/adapter/column-autodetect.ts (implementation under test)
 *
 * Coverage:
 * - Reference-fixture happy paths (theoph: 5 fields, indometh: 3 fields).
 * - AD-3 single-letter synonym removal (c / t / id never auto-fill).
 * - AD-4 BLQ INT 0/1 acceptance + INT-with-cardinality-or-value-violation rejection.
 * - AD-5 actualTime vs nominalTime bucket structure (PRD FR-106: nominal time never fills the AUC time slot).
 * - Pure-function contract (reference-equal results across calls; no mutation).
 * - Case-insensitive name match; type-mismatch rejection.
 *
 * Unit tests construct synthetic DataFrames via `DG.DataFrame.fromColumns([...])`;
 * integration tests load the bundled theoph / indometh CSVs via
 * `grok.dapi.files.readAsText` (same path the launcher uses).
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {detectMapping, isBoolLikeNumeric, detectEventDrivenColumns} from '../adapter/column-autodetect';
import {loadFixture} from './test-utils';

function dfFromColumns(...cols: DG.Column[]): DG.DataFrame {
  return DG.DataFrame.fromColumns(cols);
}

category('Adapter: event-driven (NONMEM) detection', () => {
  test('AC-A1: 03 NONMEM → isEventDriven + evid/amt/cmt/rate/mdv resolved, no addl/ii', async () => {
    const df = await loadFixture('03_rat_multiroute_nonmem.csv');
    const d = detectEventDrivenColumns(df);
    expect(d.isEventDriven, true, '03 triggers Event-driven mode');
    expect(d.evid !== undefined, true, 'EVID resolved');
    expect(d.amt !== undefined, true, 'AMT resolved');
    expect(d.cmt !== undefined, true, 'CMT resolved');
    expect(d.rate !== undefined, true, 'RATE resolved');
    expect(d.mdv !== undefined, true, 'MDV resolved');
    expect(d.addl === undefined, true, 'no ADDL in 03');
    expect(d.ii === undefined, true, 'no II in 03');
  }, {timeout: 30000});

  test('AC-A2: 06 NONMEM additionally resolves ADDL + II', async () => {
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const d = detectEventDrivenColumns(df);
    expect(d.isEventDriven, true, '06 triggers Event-driven mode');
    expect(d.addl !== undefined, true, 'ADDL resolved in 06');
    expect(d.ii !== undefined, true, 'II resolved in 06');
  }, {timeout: 30000});

  test('AC-A3: theoph (no EVID/AMT) → isEventDriven=false (dense study byte-unchanged)', async () => {
    const df = await loadFixture('theoph.csv');
    expect(detectEventDrivenColumns(df).isEventDriven, false, 'theoph never triggers NONMEM mode');
  }, {timeout: 30000});

  test('AC-A4: EVID present but no AMT → isEventDriven=false (three-part conjunction)', async () => {
    const df = dfFromColumns(
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ID', ['1', '1']),
      DG.Column.fromInt32Array('EVID', Int32Array.from([1, 0])),
      DG.Column.fromFloat64Array('DV', Float64Array.from([NaN, 9])),
    );
    expect(detectEventDrivenColumns(df).isEventDriven, false, 'EVID without AMT must not trigger');
  });

  test('a coincidental EVID with no valid {0..4} value → isEventDriven=false', async () => {
    const df = dfFromColumns(
      DG.Column.fromInt32Array('EVID', Int32Array.from([7, 9])),
      DG.Column.fromFloat64Array('AMT', Float64Array.from([1, 1])),
    );
    expect(detectEventDrivenColumns(df).isEventDriven, false, 'EVID with no {0..4} row must not trigger');
  });
});

function strCol(name: string, values: ReadonlyArray<string>): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, Array.from(values));
}

function floatCol(name: string, values: ReadonlyArray<number>): DG.Column {
  return DG.Column.fromFloat64Array(name, Float64Array.from(values));
}

function intCol(name: string, values: ReadonlyArray<number>): DG.Column {
  return DG.Column.fromInt32Array(name, Int32Array.from(values));
}

function intColWithNulls(name: string, values: ReadonlyArray<number | null>): DG.Column {
  const buf = new Int32Array(values.length);
  for (let i = 0; i < values.length; i++)
    buf[i] = values[i] === null ? DG.INT_NULL : (values[i] as number);
  return DG.Column.fromInt32Array(name, buf);
}

function boolCol(name: string, values: ReadonlyArray<boolean>): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.BOOL, name, Array.from(values));
}

category('Adapter: column autodetect', () => {
  // -------------------------------------------------------------------
  // Integration tests — bundled reference fixtures
  // -------------------------------------------------------------------

  test('integration: theoph fixture — 5 fields auto-fill (Subject/Time/Conc/Dose/Wt)', async () => {
    const df = await loadFixture('theoph.csv');
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'Subject', 'theoph.Subject → subject');
    expect(detected.time?.name, 'Time', 'theoph.Time → time');
    expect(detected.concentration?.name, 'conc', 'theoph.conc → concentration');
    expect(detected.dose?.name, 'Dose', 'theoph.Dose → dose');
    expect(detected.bodyWeight?.name, 'Wt', 'theoph.Wt → bodyWeight');
    expect(detected.blqFlag === undefined, true, 'theoph has no BLQ column');
    expect(detected.lloq === undefined, true, 'theoph has no LLOQ column');
  });

  test('route autodetect: a unique Route/EXROUTE column fills the route field (AC-A2)', async () => {
    const df = dfFromColumns(
      strCol('SubjectID', ['A', 'B']),
      floatCol('Time_hr', [0, 1]),
      floatCol('Conc_ng_mL', [10, 5]),
      strCol('Route', ['IV', 'PO']),
    );
    const detected = detectMapping(df);
    expect(detected.route?.name, 'Route', 'unique Route column → route field');
    // Ambiguous: two route-synonym columns → none (uniqueMatch semantics).
    const df2 = dfFromColumns(
      strCol('Route', ['IV']), strCol('EXROUTE', ['IV']),
      strCol('SubjectID', ['A']), floatCol('Time_hr', [0]), floatCol('Conc_ng_mL', [1]),
    );
    expect(detectMapping(df2).route === undefined, true, 'two route synonyms → ambiguous → none');
  });

  test('integration: indometh fixture — 3 required fields auto-fill, optionals empty', async () => {
    const df = await loadFixture('indometh.csv');
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'Subject', 'indometh.Subject → subject');
    expect(detected.time?.name, 'time', 'indometh.time → time');
    expect(detected.concentration?.name, 'conc', 'indometh.conc → concentration');
    expect(detected.dose === undefined, true, 'indometh has no Dose column');
    expect(detected.bodyWeight === undefined, true, 'indometh has no Wt column');
    expect(detected.blqFlag === undefined, true, 'indometh has no BLQ column');
    expect(detected.lloq === undefined, true, 'indometh has no LLOQ column');
  });

  // -------------------------------------------------------------------
  // AD-3 — single-letter synonyms `c` / `t` / `id` must NOT auto-fill
  // -------------------------------------------------------------------

  test('AD-3 (E9a): bare `c` does NOT fill concentration', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1']),
      floatCol('Time', [0]),
      floatCol('c', [1.0]),
    );
    const detected = detectMapping(df);
    expect(detected.concentration === undefined, true,
      `bare 'c' must NOT auto-fill concentration (AD-3); got '${detected.concentration?.name}'`);
  });

  test('AD-3 (E9b): bare `t` does NOT fill time (no actualTime peer)', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1']),
      floatCol('t', [0]),
      floatCol('conc', [1.0]),
    );
    const detected = detectMapping(df);
    expect(detected.time === undefined, true,
      `bare 't' must NOT auto-fill time (AD-3); got '${detected.time?.name}'`);
  });

  test('AD-3 (E9c): bare `id` does NOT fill subject', async () => {
    const df = dfFromColumns(
      intCol('id', [1, 2, 3]),
      floatCol('time', [0, 1, 2]),
      floatCol('conc', [10, 5, 2]),
    );
    const detected = detectMapping(df);
    expect(detected.subject === undefined, true,
      `bare 'id' must NOT auto-fill subject (AD-3); got '${detected.subject?.name}'`);
    // `time` and `conc` should still auto-fill (sanity)
    expect(detected.time?.name, 'time', 'time should still auto-fill from `time` column');
    expect(detected.concentration?.name, 'conc', 'conc should still auto-fill');
  });

  // -------------------------------------------------------------------
  // GAP-14 — preclinical species-id subject columns (MouseID / AnimalNo)
  // -------------------------------------------------------------------

  test('GAP-14: `MouseID` auto-fills subject (preclinical species-id form)', async () => {
    const df = dfFromColumns(
      intCol('MouseID', [1, 2, 3]),
      floatCol('time', [0, 1, 2]),
      floatCol('conc', [10, 5, 2]),
    );
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'MouseID',
      `'MouseID' must auto-fill subject (GAP-14); got '${detected.subject?.name}'`);
  });

  test('GAP-14: `AnimalNo` auto-fills subject (non-conventional fixture)', async () => {
    const df = dfFromColumns(
      strCol('AnimalNo', ['A1', 'A2', 'A3']),
      floatCol('time', [0, 1, 2]),
      floatCol('conc', [10, 5, 2]),
    );
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'AnimalNo',
      `'AnimalNo' must auto-fill subject (GAP-14); got '${detected.subject?.name}'`);
  });

  test('GAP-14 guard: `StudyID` / `GroupID` do NOT false-match subject (no bare `id` token)', async () => {
    const df = dfFromColumns(
      strCol('StudyID', ['ST1', 'ST1', 'ST1']),
      strCol('GroupID', ['G1', 'G1', 'G2']),
      floatCol('time', [0, 1, 2]),
      floatCol('conc', [10, 5, 2]),
    );
    const detected = detectMapping(df);
    expect(detected.subject === undefined, true,
      `neither 'StudyID' nor 'GroupID' may auto-fill subject (GAP-14); got '${detected.subject?.name}'`);
  });

  // -------------------------------------------------------------------
  // AD-4 — BLQ INT 0/1 acceptance (NONMEM MDV, SAS, Excel conventions)
  // -------------------------------------------------------------------

  test('AD-4 (E9d): INT 0/1 column named `BLQ` is auto-detected', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1', 'S1', 'S1']),
      floatCol('time', [0, 1, 2, 3]),
      floatCol('conc', [10, 5, 2, 1]),
      intCol('BLQ', [0, 1, 0, 1]),
    );
    const detected = detectMapping(df);
    expect(detected.blqFlag?.name, 'BLQ', `INT 0/1 BLQ must auto-detect (AD-4); got '${detected.blqFlag?.name}'`);
  });

  test('AD-4 guard: INT column named `BLQ` with values [0,1,2] is NOT detected (cardinality > 2)', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1', 'S1']),
      floatCol('time', [0, 1, 2]),
      floatCol('conc', [10, 5, 2]),
      intCol('BLQ', [0, 1, 2]),
    );
    const detected = detectMapping(df);
    expect(detected.blqFlag === undefined, true,
      `INT BLQ with 3+ distinct values must NOT auto-detect; got '${detected.blqFlag?.name}'`);
  });

  test('AD-4 guard: INT column named `BLQ` with values [0,5] is NOT detected (value outside {0,1})', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('time', [0, 1]),
      floatCol('conc', [10, 5]),
      intCol('BLQ', [0, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.blqFlag === undefined, true,
      `INT BLQ with value 5 must NOT auto-detect; got '${detected.blqFlag?.name}'`);
  });

  test('AD-4 nulls (FM-7): INT 0/1 BLQ with null values [0, INT_NULL, 1, 0] IS detected', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1', 'S1', 'S1']),
      floatCol('time', [0, 1, 2, 3]),
      floatCol('conc', [10, 5, 2, 1]),
      intColWithNulls('BLQ', [0, null, 1, 0]),
    );
    const detected = detectMapping(df);
    expect(detected.blqFlag?.name, 'BLQ',
      `INT 0/1 BLQ with nulls must auto-detect (FM-7 nulls don't inflate uniqueCount); ` +
      `got '${detected.blqFlag?.name}'`);
  });

  test('AD-4 legacy: BOOL column named `BLQ` is auto-detected (existing encoding still works)', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1', 'S1', 'S1']),
      floatCol('time', [0, 1, 2, 3]),
      floatCol('conc', [10, 5, 2, 1]),
      boolCol('BLQ', [false, true, false, true]),
    );
    const detected = detectMapping(df);
    expect(detected.blqFlag?.name, 'BLQ', `BOOL BLQ must still auto-detect; got '${detected.blqFlag?.name}'`);
  });

  // -------------------------------------------------------------------
  // AD-5 — actualTime vs nominalTime bucket structure (PRD FR-106: nominal time never fills the AUC slot)
  // -------------------------------------------------------------------

  test('AD-5 (E9f / case 1): `actual_time` + `nominal_time` → time fills from actual_time', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('actual_time', [0, 1.05]),
      floatCol('nominal_time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time?.name, 'actual_time',
      `actualTime wins over nominalTime (PRD FR-106); got '${detected.time?.name}'`);
  });

  test('AD-5 (E9e): `nominal_time` alone → time stays empty (no fallback)', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('nominal_time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time === undefined, true,
      `nominalTime must NOT promote to time calculation slot (AD-5); got '${detected.time?.name}'`);
  });

  test('AD-5 fallback: bare `time` alone (no actual_time, no nominal_time) → time fills', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time?.name, 'time',
      `bare 'time' alone must fill (AD-5 — 'time' is in actualTime bucket); got '${detected.time?.name}'`);
  });

  test('AD-5 (R1 C3): ambiguous actualTime bucket `time_h` + `time_min` → time stays empty', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('time_h', [0, 1]),
      floatCol('time_min', [0, 60]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time === undefined, true,
      `2 actualTime candidates must leave time empty; got '${detected.time?.name}'`);
  });

  test('AD-5 (R1 C3 negative): `time_h` + `time_min` + `nominal_time` → time still empty (no fallback)', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('time_h', [0, 1]),
      floatCol('time_min', [0, 60]),
      floatCol('nominal_time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time === undefined, true,
      `ambiguous actualTime must NOT fall through to nominalTime (AD-5); got '${detected.time?.name}'`);
  });

  test('AD-5 (R2 N2): `actual_time` + bare `time` both → time empty (actualTime bucket ambiguity)', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('actual_time', [0, 1.05]),
      floatCol('time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time === undefined, true,
      `2 actualTime candidates (actual_time + bare time) must leave time empty; got '${detected.time?.name}'`);
  });

  test('AD-5 (R2 N3 — case 3): 1 actualTime + 2 nominalTime → time fills from actualTime', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('time', [0, 1]),
      floatCol('nominal_time', [0, 1]),
      floatCol('nomtime', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time?.name, 'time',
      `nominalTime ambiguity must NOT block time fill from actualTime; got '${detected.time?.name}'`);
  });

  // -------------------------------------------------------------------
  // General contract — case-insensitive, type filter, empty DF, purity
  // -------------------------------------------------------------------

  test('case-insensitive match: `USUBJID` matches `usubjid` synonym', async () => {
    const df = dfFromColumns(
      strCol('USUBJID', ['SUBJ-001']),
      floatCol('time', [0]),
      floatCol('conc', [1.0]),
    );
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'USUBJID',
      `case-insensitive match must accept USUBJID (CDISC PC); got '${detected.subject?.name}'`);
  });

  test('type-mismatch rejection: FLOAT col named `subject` NOT detected (filter is STRING|INT)', async () => {
    const df = dfFromColumns(
      floatCol('subject', [1.5, 2.7]),
      floatCol('time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.subject === undefined, true,
      `FLOAT 'subject' must NOT pass type filter; got '${detected.subject?.name}'`);
  });

  test('empty / no-match DF: all fields undefined', async () => {
    const df = dfFromColumns(
      strCol('Foo', ['x']),
      floatCol('Bar', [1.0]),
      floatCol('Baz', [2.0]),
    );
    const detected = detectMapping(df);
    expect(detected.subject === undefined, true, 'subject empty');
    expect(detected.time === undefined, true, 'time empty');
    expect(detected.concentration === undefined, true, 'concentration empty');
    expect(detected.dose === undefined, true, 'dose empty');
    expect(detected.bodyWeight === undefined, true, 'bodyWeight empty');
    expect(detected.blqFlag === undefined, true, 'blqFlag empty');
    expect(detected.lloq === undefined, true, 'lloq empty');
  });

  test('pure-function contract: two calls return reference-equal columns', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('time', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const a = detectMapping(df);
    const b = detectMapping(df);
    expect(a.subject === b.subject, true, 'subject ref-equal across calls');
    expect(a.time === b.time, true, 'time ref-equal across calls');
    expect(a.concentration === b.concentration, true, 'concentration ref-equal across calls');
  });

  // -------------------------------------------------------------------
  // Unit-suffixed preclinical / CRO column names (e.g. rat dataset 03)
  // -------------------------------------------------------------------

  test('unit-suffixed: Time_hr / Conc_ng_mL auto-fill (rat schema)', async () => {
    const df = dfFromColumns(
      strCol('SubjectID', ['R001', 'R001']),
      floatCol('Time_hr', [0, 1]),
      floatCol('Conc_ng_mL', [100, 50]),
    );
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'SubjectID', 'SubjectID → subject');
    expect(detected.time?.name, 'Time_hr', 'Time_hr → time');
    expect(detected.concentration?.name, 'Conc_ng_mL', 'Conc_ng_mL → concentration');
  });

  test('unit-suffixed: BodyWeight_g / Dose_mg_kg auto-fill (rat schema)', async () => {
    const df = dfFromColumns(
      strCol('SubjectID', ['R001', 'R001']),
      floatCol('Time_hr', [0, 1]),
      floatCol('Conc_ng_mL', [100, 50]),
      floatCol('BodyWeight_g', [244.1, 244.1]),
      floatCol('Dose_mg_kg', [1, 1]),
    );
    const detected = detectMapping(df);
    expect(detected.bodyWeight?.name, 'BodyWeight_g', 'BodyWeight_g → bodyWeight');
    expect(detected.dose?.name, 'Dose_mg_kg', 'Dose_mg_kg → dose');
  });

  test('unit-suffixed: Conc_µg_mL auto-fills concentration', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1']),
      floatCol('time', [0]),
      floatCol('Conc_µg_mL', [1.0]),
    );
    const detected = detectMapping(df);
    expect(detected.concentration?.name, 'Conc_µg_mL', 'Conc_µg_mL → concentration');
  });

  test('SCIENCE-FLAG #1: unit-suffixed NominalTime_hr does NOT auto-fill time', async () => {
    // The unit-suffixed time variants are ACTUAL-only; a nominal-time column —
    // even unit-suffixed — must stay out of the AUC-calculation slot.
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('NominalTime_hr', [0, 1]),
      floatCol('conc', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.time === undefined, true,
      `NominalTime_hr must NOT auto-fill time; got '${detected.time?.name}'`);
  });

  test('integration: 03_rat_iv fixture — Subject/Time/Conc/Dose/BodyWeight auto-fill', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    const detected = detectMapping(df);
    expect(detected.subject?.name, 'SubjectID', 'rat SubjectID → subject');
    expect(detected.time?.name, 'Time_hr', 'rat Time_hr → time');
    expect(detected.concentration?.name, 'Conc_ng_mL', 'rat Conc_ng_mL → concentration');
    expect(detected.dose?.name, 'Dose_mg_kg', 'rat Dose_mg_kg → dose');
    expect(detected.bodyWeight?.name, 'BodyWeight_g', 'rat BodyWeight_g → bodyWeight');
  });

  // -------------------------------------------------------------------
  // Unit autofill (GAP-11) — parse the unit token from the matched column
  // name so the mapping-pane dropdown pre-selects the right unit instead of
  // silently leaving canonical (the 1000× CL/Vz/AUC mis-scaling trap).
  // -------------------------------------------------------------------

  test('unit autofill integration: 03_rat_iv fixture — ng/mL · mg/kg · h · g parsed', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    const detected = detectMapping(df);
    expect(detected.concentrationUnit, 'ng/mL', 'Conc_ng_mL → concentrationUnit ng/mL');
    expect(detected.doseUnit, 'mg/kg', 'Dose_mg_kg → doseUnit mg/kg');
    expect(detected.timeUnit, 'h', 'Time_hr → timeUnit h');
    expect(detected.bodyWeightUnit, 'g', 'BodyWeight_g → bodyWeightUnit g');
  });

  test('unit autofill: unsuffixed columns leave units undefined (canonical default kept)', async () => {
    // theoph columns carry no unit suffix → no unit emitted → pane keeps
    // its canonical mg/L · h · mg defaults (no silent mis-guess).
    const df = await loadFixture('theoph.csv');
    const detected = detectMapping(df);
    expect(detected.concentration?.name, 'conc', 'theoph conc still maps');
    expect(detected.concentrationUnit === undefined, true, 'bare conc → no unit');
    expect(detected.timeUnit === undefined, true, 'bare Time → no unit');
    expect(detected.doseUnit === undefined, true, 'bare Dose → no unit');
  });

  test('unit autofill: µg/mL variant and minutes parse to dropdown values', async () => {
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('Time_min', [0, 30]),
      floatCol('Conc_µg_mL', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.concentrationUnit, 'µg/mL', 'Conc_µg_mL → µg/mL');
    expect(detected.timeUnit, 'min', 'Time_min → min');
  });

  test('unit autofill: unsupported token (mg/mL conc) emits no unit', async () => {
    // mg/mL is a valid synonym match for the column but NOT a Concentration
    // dropdown choice → conservative: leave the unit unset rather than guess.
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('Time', [0, 1]),
      floatCol('Conc_mg_mL', [10, 5]),
    );
    const detected = detectMapping(df);
    expect(detected.concentration?.name, 'Conc_mg_mL', 'Conc_mg_mL still maps as concentration');
    expect(detected.concentrationUnit === undefined, true,
      `mg/mL is not a dropdown choice → no unit; got '${detected.concentrationUnit}'`);
  });

  test('unit autofill: longest-suffix-first ordering — no shorter-suffix shadowing', async () => {
    // Anti-shadow regression: a shorter suffix must NOT pre-empt a longer one.
    // `_mg_kg` must win over `_mg` (else Dose_mg_kg → 'mg', a 0× vs per-kg
    // mis-scale); `_kg` must win over `_g` (else BodyWeight_kg → 'g', a 1000×
    // BW mis-scale). And a column matched by one field's synonym must not pick
    // up another field's suffix (`Dose_mg` is not a BW column at all).
    const df = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('Time', [0, 1]),
      floatCol('Conc', [10, 5]),
      floatCol('Dose_mg_kg', [5, 5]),
      floatCol('BW_kg', [0.3, 0.3]),
    );
    const detected = detectMapping(df);
    expect(detected.doseUnit, 'mg/kg', '_mg_kg must win over _mg (not shadowed to mg)');
    expect(detected.bodyWeightUnit, 'kg', '_kg must win over _g (not shadowed to g)');
    // `_mg` alone still resolves correctly when it is the actual suffix.
    const df2 = dfFromColumns(
      strCol('Subject', ['S1', 'S1']),
      floatCol('Time', [0, 1]),
      floatCol('Conc', [10, 5]),
      floatCol('Dose_mg', [5, 5]),
    );
    const d2 = detectMapping(df2);
    expect(d2.doseUnit, 'mg', 'Dose_mg → mg (bare _mg suffix)');
    expect(d2.bodyWeightUnit === undefined, true, 'Dose_mg is not a BW column → no bodyWeightUnit');
  });

  // -------------------------------------------------------------------
  // isBoolLikeNumeric — shared predicate used by detector + picker filter
  // -------------------------------------------------------------------

  test('isBoolLikeNumeric: BOOL column → true', async () => {
    const col = boolCol('flag', [true, false, true]);
    expect(isBoolLikeNumeric(col), true, 'BOOL must pass');
  });

  test('isBoolLikeNumeric: INT 0/1 column → true', async () => {
    const col = intCol('flag', [0, 1, 0, 1]);
    expect(isBoolLikeNumeric(col), true, 'INT 0/1 must pass');
  });

  test('isBoolLikeNumeric: INT 0/2 column → false', async () => {
    const col = intCol('flag', [0, 2, 0]);
    expect(isBoolLikeNumeric(col), false, 'INT with 2 must fail');
  });

  test('isBoolLikeNumeric: FLOAT column → false', async () => {
    const col = floatCol('flag', [0, 1, 0]);
    expect(isBoolLikeNumeric(col), false, 'FLOAT must fail');
  });

  test('isBoolLikeNumeric: empty INT column → false', async () => {
    const col = DG.Column.fromInt32Array('flag', new Int32Array(0));
    expect(isBoolLikeNumeric(col), false, 'empty INT must fail (no observable 0/1 pattern)');
  });
});
