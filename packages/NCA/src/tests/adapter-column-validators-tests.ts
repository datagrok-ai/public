/**
 * Adapter — column validators category.
 *
 * @see src/adapter/validators.ts (implementation under test)
 *
 * Each validation code is exercised by a synthetic DataFrame that violates
 * exactly one rule, plus a happy-path test that runs the contract end-to-end
 * with clean data.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {
  validateObservations,
  validateScalarsForDose,
  validateScalarsForMolar,
  validateMultiAnalyte,
  validateLloqSanity,
  validateSparseDesign,
} from '../adapter/validators';
import type {ColumnContract, ValidationCode} from '../adapter/column-contract';
import type {ProfileScalars} from '../adapter/profile-inputs';
import {buildProfileIndices} from '../adapter/profile-index';
import {isCompoundIdColumnName} from '../adapter/column-autodetect';
import {computeSparseNca} from '../adapter/sparse-nca';
import {poolCompositeByCell} from '../adapter/composite-pool';
import {normalizeEventDrivenObservations, SUBJECT_COL, type EventContract} from '../adapter/event-driven';
import {loadFixture} from './test-utils';

function baseContract(): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: 'string'},
      time: {name: 'Time', required: true, rawType: 'float64'},
      concentration: {name: 'conc', required: true, rawType: 'float64'},
    },
  };
}

function cleanDataFrame(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S2', 'S2']),
    DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1])),
    DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 8, 4])),
  ]);
}

function codes(errs: ReadonlyArray<{code: ValidationCode}>): ValidationCode[] {
  return errs.map((e) => e.code);
}

category('Adapter: column validators', () => {
  test('happy path: clean DataFrame yields 0 errors', async () => {
    const errs = validateObservations(cleanDataFrame(), baseContract());
    expect(errs.length, 0, `expected zero errors, got ${JSON.stringify(errs)}`);
  });

  test('MISSING_COLUMN: required column absent', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S2']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
    ]);
    const errs = validateObservations(df, baseContract());
    expect(errs.some((e) => e.code === 'MISSING_COLUMN' && e.column === 'conc'), true,
      `expected MISSING_COLUMN on 'conc'; got ${JSON.stringify(codes(errs))}`);
  });

  test('WRONG_TYPE: string column mapped where float64 expected', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S2']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Time', ['0', '1']),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const errs = validateObservations(df, baseContract());
    expect(errs.some((e) => e.code === 'WRONG_TYPE' && e.column === 'Time'), true,
      `expected WRONG_TYPE on 'Time' (STRING mapped to float64); got ${JSON.stringify(codes(errs))}`);
  });

  test('INT column accepted for rawType=float64 (NOTE-04 closure)', async () => {
    // Phase 3 extractor (extractFloat64) widens INT-backed columns to Float64
    // via INT_NULL → NaN translation. Validator must agree — closing the
    // strict-validator / permissive-extractor inconsistency logged as NOTE-04.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S2']),
      DG.Column.fromInt32Array('Time', Int32Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const errs = validateObservations(df, baseContract());
    expect(errs.some((e) => e.code === 'WRONG_TYPE' && e.column === 'Time'), false,
      `INT column on rawType=float64 must NOT raise WRONG_TYPE; got ${JSON.stringify(errs)}`);
  });

  test('NEGATIVE_TIME: row with time < 0', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, -1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const errs = validateObservations(df, baseContract());
    const negs = errs.filter((e) => e.code === 'NEGATIVE_TIME');
    expect(negs.length, 1, `expected 1 NEGATIVE_TIME; got ${JSON.stringify(errs)}`);
    expect(negs[0].row, 1, 'expected NEGATIVE_TIME at row index 1');
  });

  test('NON_NUMERIC_CONC: +Infinity in conc column', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, Infinity])),
    ]);
    const errs = validateObservations(df, baseContract());
    const non = errs.filter((e) => e.code === 'NON_NUMERIC_CONC');
    expect(non.length, 1, `expected 1 NON_NUMERIC_CONC; got ${JSON.stringify(errs)}`);
    expect(non[0].row, 1, 'expected NON_NUMERIC_CONC at row index 1');
  });

  test('DUPLICATE_SUBJECT_TIME: same (subject, time) twice', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S2']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 0, 0])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 11, 8])),
    ]);
    const errs = validateObservations(df, baseContract());
    const dup = errs.filter((e) => e.code === 'DUPLICATE_SUBJECT_TIME');
    expect(dup.length, 1, `expected 1 DUPLICATE_SUBJECT_TIME on S1@0; got ${JSON.stringify(errs)}`);
  });

  test('optional column ignored when absent (no MISSING_COLUMN for bodyWeight)', async () => {
    const contract: ColumnContract = {
      observations: {
        ...baseContract().observations,
        bodyWeight: {name: 'Wt', required: false, rawType: 'float64'},
      },
    };
    const errs = validateObservations(cleanDataFrame(), contract);
    expect(errs.length, 0, `expected 0 errors when optional col absent; got ${JSON.stringify(errs)}`);
  });
});

// =====================================================================
// nca-unit-semantics — new unit-load-bearing validators (E6, AD-4/8/9)
// =====================================================================

function baseScalars(overrides: Partial<ProfileScalars> = {}): ProfileScalars {
  return {
    dose: 100,
    route: nca.ROUTE_IV_BOLUS,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    bodyWeight: null,
    bodyWeightUnit: 'kg',
    infusionDuration: null,
    lloq: 0,
    ...overrides,
  };
}

category('Adapter: validators (unit semantics)', () => {
  // -------------- /kg dose + no body weight (nca-cl-vz-convention AC-1) --------------
  // P-A / AD-1: a /kg dose with no body-weight source is NO LONGER a block. It
  // is the valid per-kg path; validateScalarsForDose returns []. The per-kg
  // advisory (PER_KG_NO_BW) rides provenance.warnings from scaleProfileInputs
  // (asserted in adapter-unit-conversion-tests.ts), not the validator.

  test('AC-1: mg/kg dose + null BW scalar + no BW column → NO LONGER blocks (returns [])', async () => {
    const errs = validateScalarsForDose(
      {observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      }},
      baseScalars({doseUnits: 'mg/kg', bodyWeight: null}),
    );
    expect(errs.length, 0, `mg/kg+no-BW is now the per-kg path, not a block; got ${JSON.stringify(errs)}`);
  });

  test('AC-1: µg/kg dose + null BW scalar + no BW column → NO LONGER blocks (returns [])', async () => {
    const errs = validateScalarsForDose(
      {observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      }},
      baseScalars({doseUnits: 'µg/kg', bodyWeight: null}),
    );
    expect(errs.length, 0, `µg/kg+no-BW is now the per-kg path, not a block; got ${JSON.stringify(errs)}`);
  });

  // The pre-Run per-kg confirmation (`requiresPerKgConfirmation`) was retired with per kg becoming
  // the reported basis everywhere.

  // -------------- UNITS_MISMATCH (molar — Path C) --------------

  test('UNITS_MISMATCH: µmol/L emits with molar discriminator', async () => {
    const errs = validateScalarsForMolar(baseScalars({concentrationUnits: 'µmol/L'}));
    expect(errs.length, 1, `expected 1 UNITS_MISMATCH for µmol/L`);
    expect(errs[0].code, 'UNITS_MISMATCH');
    expect(errs[0].details?.kind, 'molar-not-supported');
    expect(errs[0].details?.unit, 'µmol/L');
  });

  test('UNITS_MISMATCH: nmol/L emits with molar discriminator', async () => {
    const errs = validateScalarsForMolar(baseScalars({concentrationUnits: 'nmol/L'}));
    expect(errs.length, 1);
    expect(errs[0].code, 'UNITS_MISMATCH');
    expect(errs[0].details?.kind, 'molar-not-supported');
    expect(errs[0].details?.unit, 'nmol/L');
  });

  test('UNITS_MISMATCH: mass-based units (mg/L, ng/mL, µg/L, µg/mL) → no emit', async () => {
    for (const u of ['mg/L', 'ng/mL', 'µg/L', 'µg/mL'] as const) {
      const errs = validateScalarsForMolar(baseScalars({concentrationUnits: u}));
      expect(errs.length, 0, `mass unit ${u} must not raise UNITS_MISMATCH`);
    }
  });

  // -------------- MULTI_ANALYTE_NOT_SUPPORTED --------------

  test('MULTI_ANALYTE_NOT_SUPPORTED: Parent_/Metabolite_ columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('Parent_Drug', Float64Array.from([10, 5])),
      DG.Column.fromFloat64Array('Metabolite_M1', Float64Array.from([1, 0.5])),
    ]);
    const errs = validateMultiAnalyte(df, 'Subject');
    expect(errs.length, 1, `expected 1 MULTI_ANALYTE_NOT_SUPPORTED`);
    expect(errs[0].code, 'MULTI_ANALYTE_NOT_SUPPORTED');
  });

  test('MULTI_ANALYTE_NOT_SUPPORTED: per-row analyte categorical col', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 1, 0.5])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'analyte',
        ['Parent', 'Parent', 'M1', 'M1']),
    ]);
    const errs = validateMultiAnalyte(df, 'Subject');
    expect(errs.length, 1, `expected 1 MULTI_ANALYTE_NOT_SUPPORTED for per-row analyte cat col`);
  });

  test('MULTI_ANALYTE_NOT_SUPPORTED: single-analyte CSV → no emit', async () => {
    const errs = validateMultiAnalyte(cleanDataFrame(), 'Subject');
    expect(errs.length, 0, `single-analyte CSV must not false-positive`);
  });

  // -------------- BUG-07 interim: cassette (shared subject across compounds) --------------

  test('BUG-07: the cassette fixture blocks — CompoundID varies within subject R001', async () => {
    const df = await loadFixture('07_cassette_shared_subject.csv');
    const errs = validateMultiAnalyte(df, 'SubjectID');
    expect(errs.length, 1, 'one cassette error');
    expect(errs[0].code, 'MULTI_ANALYTE_NOT_SUPPORTED');
    expect(errs[0].details?.kind, 'cassette-shared-subject');
    expect(errs[0].details?.subject, 'R001');
    expect(errs[0].column, 'CompoundID');
    expect(errs[0].message.includes('DGK-0001') && errs[0].message.includes('DGK-0002'), true,
      'the message names both compounds');
  });

  test('BUG-07: a batch screen with disjoint subjects per compound is NOT blocked', async () => {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    expect(validateMultiAnalyte(df, 'SubjectID').length, 0,
      'CompoundID constant within every subject → no cassette');
  });

  test('BUG-07: only a compound-named column triggers — BLQ / Sex varying within a subject do not', async () => {
    // The structural-only scan reverted 2026-07-14 false-positived on exactly these columns.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1', 'S2', 'S2', 'S2']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 2, 0, 1, 2])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 1, 9, 4, 1])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'BLQ', ['False', 'False', 'True', 'False', 'True', 'True']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'F', 'M', 'F', 'F', 'F']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Test Article', ['A', 'A', 'A', 'B', 'B', 'B']),
    ]);
    expect(validateMultiAnalyte(df, 'Subject').length, 0, 'no compound column varies within a subject');
  });

  test('BUG-07: recognised compound-id names (and only those)', async () => {
    for (const n of ['CompoundID', 'Compound_Name', 'CmpdID', 'Analyte Name', 'AnalyteID', 'Drug', 'drug_name',
      'Test Article', 'Article', 'TestItem', 'Molecule', 'molecule-id', 'Substance', 'Chemical', 'Metabolite',
      'DGK_ID'])
      expect(isCompoundIdColumnName(n), true, `${n} reads as a compound id`);
    for (const n of ['Sex', 'BLQ', 'Dose', 'Route', 'StudyID', 'Compounding', 'DrugLevel', 'SubjectID',
      'Parent_Conc', 'Metabolite_M1', 'Drug_Conc', 'CompoundMW', 'Articles'])
      expect(isCompoundIdColumnName(n), false, `${n} is not a compound id`);
  });

  test('BUG-07: a column already reported by the long-form signal is not reported twice', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 1, 0.5])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'analyte', ['Parent', 'Parent', 'M1', 'M1']),
    ]);
    expect(validateMultiAnalyte(df, 'Subject').length, 1, 'one error, not two');
  });

  // -------------- LLOQ_OUT_OF_RANGE (AD-8 a/c/d/e) --------------

  function lloqDf(opts: {
    conc: ReadonlyArray<number>,
    blqFlag?: ReadonlyArray<boolean>,
    lloqCol?: ReadonlyArray<number>,
  }): DG.DataFrame {
    const n = opts.conc.length;
    const cols: DG.Column[] = [
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', Array(n).fill('S1')),
      DG.Column.fromFloat64Array('Time', Float64Array.from(opts.conc.map((_, i) => i))),
      DG.Column.fromFloat64Array('conc', Float64Array.from(opts.conc)),
    ];
    if (opts.blqFlag !== undefined)
      cols.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'BLQ', Array.from(opts.blqFlag)));
    if (opts.lloqCol !== undefined)
      cols.push(DG.Column.fromFloat64Array('LLOQ', Float64Array.from(opts.lloqCol)));
    return DG.DataFrame.fromColumns(cols);
  }

  function lloqContract(includeBlq = false, includeLloqCol = false): ColumnContract {
    return {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
        ...(includeBlq ?
          {blqFlag: {name: 'BLQ', required: false, rawType: 'uint8' as const}} : {}),
        ...(includeLloqCol ?
          {lloq: {name: 'LLOQ', required: false, rawType: 'float64' as const}} : {}),
      },
    };
  }

  test('LLOQ_OUT_OF_RANGE (a): LLOQ > min(non-BLQ) × 1.5 fires', async () => {
    // lloq=1.0 chosen so condition (e) "all < LLOQ" does NOT pre-empt (a):
    // max(conc)=2.0 ≥ 1.0, so (e) is false. min(non-BLQ)=0.5; 0.5×1.5=0.75;
    // lloq=1.0 > 0.75 → (a) fires with kind 'above-min-non-blq'.
    const df = lloqDf({conc: [0.5, 1.0, 2.0]});
    const indices = buildProfileIndices(df, lloqContract());
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 1.0}), indices);
    expect(errs.length, 1, `expected condition (a) fire; got ${JSON.stringify(errs)}`);
    expect(errs[0].code, 'LLOQ_OUT_OF_RANGE');
    expect(errs[0].details?.kind, 'above-min-non-blq');
  });

  test('LLOQ_OUT_OF_RANGE (a): BLQ-coded-as-zero (no blq flag col) does NOT fire (R1-F3 guard)', async () => {
    // Watson/Thermo Fisher LIMS encode BLQ as conc=0 with no BLQ flag column.
    // "Non-BLQ" must = conc > 0 (excludes zeros from min calculation).
    const df = lloqDf({conc: [0, 0, 0.5, 1.0, 2.0]});
    const indices = buildProfileIndices(df, lloqContract());
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 0.4}), indices);
    expect(errs.length, 0, `BLQ-as-zero must not false-positive-block; got ${JSON.stringify(errs)}`);
  });

  test('LLOQ_OUT_OF_RANGE (a): explicit BLQ flag column at 0 → low-conc min used correctly', async () => {
    // Same data-shape rationale as the previous test: lloq chosen so (e)
    // doesn't pre-empt (a). With blqFlag=[false,false,false] all conc are
    // counted as non-BLQ → min(non-BLQ)=0.5 → (a) fires for lloq>0.75.
    const df = lloqDf({conc: [0.5, 1.0, 2.0], blqFlag: [false, false, false]});
    const indices = buildProfileIndices(df, lloqContract(true));
    const errs = validateLloqSanity(df, lloqContract(true), baseScalars({lloq: 1.0}), indices);
    expect(errs.length, 1, `with explicit blqFlag, min(non-BLQ)=0.5 → lloq=1.0 > 0.75 fires`);
    expect(errs[0].details?.kind, 'above-min-non-blq');
  });

  test('LLOQ_OUT_OF_RANGE (a): predose t=0=0 with blqFlag=false does NOT fire (BUG-04 / Ed #8)', async () => {
    // 04 mouse composite / 06 TK composite encode the predose row as
    // (t=0, conc=0, blqFlag=false) — semantically "no drug administered
    // yet", NOT "below detection". Before the fix the no-flag-branch
    // conc>0 guard didn't apply on the explicit-flag path, so 0 polluted
    // min(non-BLQ) and condition (a) hard-blocked legitimate composites.
    // Unified invariant `NOT(flagged BLQ) AND conc > 0` excludes the
    // predose zero from min(non-BLQ) regardless of flag tagging.
    const df = lloqDf({conc: [0, 1.0, 2.0], blqFlag: [false, false, false]});
    const indices = buildProfileIndices(df, lloqContract(true));
    const errs = validateLloqSanity(df, lloqContract(true), baseScalars({lloq: 1.2}), indices);
    expect(errs.length, 0,
      `predose t=0=0 must be excluded from min(non-BLQ); ` +
      `min(rest)=1.0, lloq=1.2 ≤ 1.5 → no fire; got ${JSON.stringify(errs)}`);
  });

  test('LLOQ_OUT_OF_RANGE (b) DROPPED in v1: high-sensitivity assay NOT hard-blocked', async () => {
    // Modern LC-MS/MS: lloq=0.001 vs min-non-BLQ=5 (ratio 5000) — condition (b)
    // would have fired but is dropped from v1 per R2-finding-4.
    const df = lloqDf({conc: [5, 50, 100]});
    const indices = buildProfileIndices(df, lloqContract());
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 0.001}), indices);
    expect(errs.length, 0, `high-sensitivity assay must NOT be hard-blocked in v1`);
  });

  test('LLOQ_OUT_OF_RANGE (c): LLOQ > max(conc) × 10 fires', async () => {
    const df = lloqDf({conc: [0.5, 2.0]});
    const indices = buildProfileIndices(df, lloqContract());
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 25}), indices);
    expect(errs.length, 1, `condition (c) must fire`);
    expect(errs[0].details?.kind, 'above-max-times-ten');
  });

  test('LLOQ_OUT_OF_RANGE (d): LLOQ ≤ 0 fires', async () => {
    const df = lloqDf({conc: [0.5, 1.0, 2.0]});
    const indices = buildProfileIndices(df, lloqContract(false, true));
    // Provide LLOQ via column so the validator path uses the lloqCol value
    const dfWithLloq = lloqDf({conc: [0.5, 1.0, 2.0], lloqCol: [-1, -1, -1]});
    const indices2 = buildProfileIndices(dfWithLloq, lloqContract(false, true));
    const errs = validateLloqSanity(dfWithLloq, lloqContract(false, true),
      baseScalars({lloq: -1}), indices2);
    expect(errs.length, 1, `LLOQ = -1 must fire condition (d); got ${JSON.stringify(errs)}`);
    expect(errs[0].details?.kind, 'non-positive');
    // Reassure indices1 was used above
    expect(indices.length, 1);
  });

  test('LLOQ_OUT_OF_RANGE (e): all conc < LLOQ fires', async () => {
    const df = lloqDf({conc: [0.1, 0.2, 0.3]});
    const indices = buildProfileIndices(df, lloqContract());
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 10}), indices);
    // Note: condition (c) max*10 = 0.3*10=3 < 10 → fires (c). condition (e) also true.
    // Validator emits the first matching condition (continue after each); (c)
    // fires first, but the test asserts at-least-one LLOQ_OUT_OF_RANGE.
    expect(errs.length >= 1, true, `condition (c) or (e) must fire`);
    expect(errs[0].code, 'LLOQ_OUT_OF_RANGE');
  });

  test('LLOQ (c)/(e) suppressed: sparse all-BLQ, BLQ-coded-as-zero (BUG-03 / Ed #7)', async () => {
    // Ed's real scenario: sparse-destructive design, 1 sample per mouse.
    // A mouse sampled at the terminal tail with conc < LLOQ encoded as 0
    // (Watson/Thermo convention, no blqFlag column). Without the fix:
    // - (c) fires (max=0, lloq > 0×10=0 for any positive lloq).
    // - (e) fires (concVals=[0] all < lloq).
    // - (a) self-skips (conc>0 guard excludes the zero from nonBlqConc).
    // With single-observation exemption: (c)/(e) skipped → 0 errors;
    // compute layer marks the profile `failed` per-subject (no Cmax) and
    // the batch continues for the remaining sparse subjects.
    const df = lloqDf({conc: [0]});
    const indices = buildProfileIndices(df, lloqContract());
    expect(indices[0].observationRows.length, 1, 'fixture must be single-row');
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 10}), indices);
    expect(errs.length, 0,
      `single-observation all-BLQ-as-zero must skip (c)/(e); got ${JSON.stringify(errs)}`);
  });

  test('LLOQ (c)/(e) suppressed: sparse all-BLQ, explicit blqFlag (BUG-03 / Ed #7)', async () => {
    // Alternative encoding of the same scenario: the lone sample carries
    // its assay value but is flagged BLQ via an explicit boolean column.
    // Without the fix: (c) fires (max=0.5, lloq=10 > 5), then (e) would
    // also (0.5 < 10). (a) self-skips because the flag empties nonBlqConc.
    // With single-observation exemption: (c)/(e) skipped → 0 errors.
    const df = lloqDf({conc: [0.5], blqFlag: [true]});
    const indices = buildProfileIndices(df, lloqContract(true));
    expect(indices[0].observationRows.length, 1, 'fixture must be single-row');
    const errs = validateLloqSanity(df, lloqContract(true), baseScalars({lloq: 10}), indices);
    expect(errs.length, 0,
      `single-observation all-BLQ-via-flag must skip (c)/(e); got ${JSON.stringify(errs)}`);
  });

  test('LLOQ sanity: no LLOQ specified (scalar=0, no lloq col) → skip (no emit)', async () => {
    const df = lloqDf({conc: [0.5, 1.0, 2.0]});
    const indices = buildProfileIndices(df, lloqContract());
    const errs = validateLloqSanity(df, lloqContract(), baseScalars({lloq: 0}), indices);
    expect(errs.length, 0, `no LLOQ → skip sanity; got ${JSON.stringify(errs)}`);
  });
});

// --- Sparse-design contract + topology validator (UC-04, Slice A) ----------

function sparseContract(): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'MouseID', required: true, rawType: 'string'},
      time: {name: 'ActualTime', required: true, rawType: 'float64'},
      concentration: {name: 'Conc', required: true, rawType: 'float64'},
      animalId: {name: 'MouseID', required: false, rawType: 'string'},
      nominalTime: {name: 'NomTime', required: false, rawType: 'float64'},
      samplingType: {name: 'SType', required: false, rawType: 'string'},
    },
  };
}

category('Adapter: sparse design validator (UC-04)', () => {
  // AC-A1: a well-formed sparse CSV (pooled animal IDs + nominal time) validates.
  test('AC-A1: a sparse CSV with the new roles yields 0 blocking errors', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'MouseID', ['M1', 'M2', 'M3', 'M4']),
      DG.Column.fromFloat64Array('ActualTime', Float64Array.from([0.25, 0.5, 1, 2])),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from([0.25, 0.5, 1, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([10, 8, 5, 2])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SType', ['Destructive', 'Destructive', 'Destructive', 'Destructive']),
    ]);
    const errs = validateObservations(df, sparseContract());
    expect(errs.length, 0, `sparse CSV should validate; got ${JSON.stringify(errs)}`);
  });

  // AC-A2: declared "destructive" but an animal spans ≥2 nominal times → flag.
  test('AC-A2: declared-destructive with paired animals flags SPARSE_TOPOLOGY_MISMATCH', async () => {
    const df = DG.DataFrame.fromColumns([
      // M1 and M2 each appear at two nominal times → actually batch/serial.
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'MouseID', ['M1', 'M2', 'M1', 'M2']),
      DG.Column.fromFloat64Array('ActualTime', Float64Array.from([1, 1, 2, 2])),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from([1, 1, 2, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([10, 12, 6, 7])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SType', ['Destructive', 'Destructive', 'Destructive', 'Destructive']),
    ]);
    const warns = validateSparseDesign(df, sparseContract());
    expect(warns.some((w) => w.code === 'SPARSE_TOPOLOGY_MISMATCH'), true,
      `expected SPARSE_TOPOLOGY_MISMATCH; got ${JSON.stringify(codes(warns))}`);
  });

  // A genuinely destructive design (1 sample/animal) raises nothing.
  test('a true destructive design (1 sample/animal) does not flag a mismatch', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'MouseID', ['M1', 'M2', 'M3', 'M4']),
      DG.Column.fromFloat64Array('ActualTime', Float64Array.from([1, 1, 2, 2])),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from([1, 1, 2, 2])),
      DG.Column.fromFloat64Array('Conc', Float64Array.from([10, 12, 6, 7])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SType', ['Destructive', 'Destructive', 'Destructive', 'Destructive']),
    ]);
    const warns = validateSparseDesign(df, sparseContract());
    expect(warns.length, 0, `true destructive should not flag; got ${JSON.stringify(codes(warns))}`);
  });

  // No-op when the sparse roles are not all mapped (dense analysis).
  test('no sparse roles mapped → validator is a no-op', async () => {
    const df = cleanDataFrame();
    expect(validateSparseDesign(df, baseContract()).length, 0, 'dense contract → no sparse check');
  });
});

// BUG-07 PATTERN (bug-stress 2026-10-05): every pooling key without a compound dimension. Animals
// carry ONE compound each here, so the cassette check cannot see it — the GROUP is what mixes.
category('Adapter: pooling across compounds (BUG-07 pattern)', () => {
  /** A destructive design: 12 animals, one sample each, at 3 nominal times; compound X on animals
   *  1-6, Y on 7-12; `doseOf` sets each compound's dose. */
  function destructive(doseOf: (compound: string) => number): DG.DataFrame {
    const animal: string[] = []; const t: number[] = []; const conc: number[] = [];
    const dose: number[] = []; const cmp: string[] = [];
    for (let a = 1; a <= 12; a++) {
      const c = a <= 6 ? 'X' : 'Y';
      const time = [0.5, 2, 8][(a - 1) % 3];
      animal.push(`A${a}`); t.push(time); conc.push(10 * Math.exp(-0.2 * time) + a * 0.01);
      dose.push(doseOf(c)); cmp.push(c);
    }
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Animal', animal),
      DG.Column.fromFloat64Array('NomTime', Float64Array.from(t)),
      DG.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
      DG.Column.fromFloat64Array('Dose', Float64Array.from(dose)),
      DG.Column.fromStrings('CompoundID', cmp),
    ]);
  }
  const sparseOpts = {nominalTimeCol: 'NomTime', animalIdCol: 'Animal', concCol: 'Conc', doseCol: 'Dose',
    groupCol: 'Dose', bootstrapIterations: 20};

  test('sparse: a dose group holding two compounds is refused (not averaged)', async () => {
    let msg = '';
    try {computeSparseNca(destructive(() => 10), sparseOpts);} catch (e) {msg = String((e as Error).message);}
    expect(msg.includes('Dose group') && msg.includes('X') && msg.includes('Y'), true,
      `refused with the mixed-compound message (got: "${msg}")`);
  });

  test('sparse: the real two-compound fixture (04c) is refused', async () => {
    // 04_mouse_sparse_destructive with animals alternating DGK-0001 / DGK-0002 inside each dose group.
    const df = await loadFixture('04c_mouse_sparse_two_compounds.csv');
    let msg = '';
    try {
      computeSparseNca(df, {nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', concCol: 'Conc_ng_mL',
        doseCol: 'Dose_mg_kg', groupCol: 'Dose_mg_kg', bootstrapIterations: 20});
    } catch (e) {msg = String((e as Error).message);}
    expect(msg.includes('Dose group \'10\'') && msg.includes('DGK-0001') && msg.includes('DGK-0002'), true,
      `refused naming the dose group and both compounds (got: "${msg}")`);
  });

  test('sparse: no dose-group column — the whole study is the composite and is refused', async () => {
    let msg = '';
    try {
      computeSparseNca(destructive(() => 10),
        {nominalTimeCol: 'NomTime', animalIdCol: 'Animal', concCol: 'Conc', bootstrapIterations: 20});
    } catch (e) {msg = String((e as Error).message);}
    expect(msg.startsWith('The study (no dose-group column) pools'), true,
      `names the study, not "Dose group 'all'" (got: "${msg}")`);
  });

  test('sparse: one compound per dose group computes as before', async () => {
    const r = computeSparseNca(destructive((c) => c === 'X' ? 10 : 50), sparseOpts);
    expect(r.groups.length, 2, 'two single-compound dose groups');
  });

  test('NONMEM composite: a dose×sex×day cell holding two compounds is reported', async () => {
    const EC06: EventContract = {
      id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
      addl: 'ADDL', ii: 'II', route: 'ROUTE', bodyWeight: 'BW_kg', covariates: ['Sex', 'DOSE_MGKG'], day: 'DAY',
    };
    const norm = normalizeEventDrivenObservations(await loadFixture('06_rat_tk28d_nonmem.csv'), EC06).observations;
    const pool = (): ReturnType<typeof poolCompositeByCell> =>
      poolCompositeByCell(norm, {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex'});
    expect(pool().mixedCompound.length, 0, 'the real single-compound TK fixture pools cleanly');
    // Two compounds dosed in the same cells: alternate by animal (the pseudo-subject's animal part).
    const subj = norm.col(SUBJECT_COL)!;
    const animals = [...new Set(Array.from({length: norm.rowCount}, (_, i) => subj.getString(i).split('_')[0]))];
    norm.columns.add(DG.Column.fromStrings('CompoundID', Array.from({length: norm.rowCount},
      (_, i) => animals.indexOf(subj.getString(i).split('_')[0]) % 2 === 0 ? 'DGK-A' : 'DGK-B')));
    const mixed = pool().mixedCompound;
    expect(mixed.length, 1, 'the compound column is reported');
    expect(mixed[0].column, 'CompoundID');
    expect(mixed[0].values.length >= 2, true, 'with the compounds the cell mixes');
  });
});
