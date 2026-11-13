/**
 * Adapter — profile-inputs category.
 *
 * @see src/adapter/profile-inputs.ts (implementation under test)
 *
 * Covers: copy-not-view semantics, Float32 widening warning, BLQ inference
 * from LLOQ column, route per-row INT_NULL fallback, and the route-string
 * → RouteCode mapping.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {extractProfileInputs, routeStringToCode} from '../adapter/profile-inputs';
import type {ProfileScalars} from '../adapter/profile-inputs';
import {buildProfileIndices} from '../adapter/profile-index';
import type {ColumnContract} from '../adapter/column-contract';

function baseScalars(): ProfileScalars {
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
  };
}

category('Adapter: profile inputs', () => {
  test('extracts fresh Float64Arrays (copy, not view)', async () => {
    const timeArr = Float64Array.from([0, 1, 2]);
    const concArr = Float64Array.from([10, 5, 2.5]);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1']),
      DG.Column.fromFloat64Array('Time', timeArr),
      DG.Column.fromFloat64Array('conc', concArr),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const {inputs, warnings} = extractProfileInputs(df, index, contract, baseScalars());

    expect(inputs.time instanceof Float64Array, true);
    expect(inputs.conc instanceof Float64Array, true);
    expectArray(Array.from(inputs.time), [0, 1, 2]);
    expectArray(Array.from(inputs.conc), [10, 5, 2.5]);

    // Mutate the source — the extracted arrays must NOT change (copy invariant).
    const sourceRaw = df.col('Time')!.getRawData() as Float64Array;
    sourceRaw[0] = 999;
    expect(inputs.time[0], 0, 'extracted time[0] must not see source mutation');

    expect(warnings.length, 0, `expected no warnings, got ${JSON.stringify(warnings)}`);
  });

  test('FLOAT32_WIDENED warning when source column is Float32', async () => {
    const f32 = Float32Array.from([0, 1, 2]);
    const concArr = Float64Array.from([10, 5, 2.5]);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1']),
      DG.Column.fromFloat32Array('Time', f32),
      DG.Column.fromFloat64Array('conc', concArr),
    ]);
    // Validator would reject this with WRONG_TYPE (Float32 stored as Datagrok
    // FLOAT column has type FLOAT — same as Float64); the extractor doesn't
    // run the validator, so it sees the Float32 path and widens.
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const {inputs, warnings} = extractProfileInputs(df, index, contract, baseScalars());

    expect(inputs.time instanceof Float64Array, true, 'widened to Float64');
    const flagged = warnings.find((w) => w.code === 'FLOAT32_WIDENED');
    expect(flagged !== undefined, true, `FLOAT32_WIDENED expected; got ${JSON.stringify(warnings)}`);
    expect(flagged!.severity, 'warning');
  });

  test('BLQ inferred from conc=0 when LLOQ column present and no BLQ flag', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 2])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([0, 5, 0])),
      DG.Column.fromFloat64Array('LLOQ', Float64Array.from([0.1, 0.1, 0.1])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
        lloq: {name: 'LLOQ', required: false, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const {inputs} = extractProfileInputs(df, index, contract, baseScalars());
    expectArray(Array.from(inputs.blqMask), [1, 0, 1]);
  });

  test('BLQ all zeros when no flag column AND no LLOQ column', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([0, 5])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const {inputs} = extractProfileInputs(df, index, contract, baseScalars());
    expectArray(Array.from(inputs.blqMask), [0, 0]);
  });

  test('route per-row INT_NULL falls back to scalars.route (no silent IV-bolus mislabel)', async () => {
    const routeArr = Int32Array.from([DG.INT_NULL, DG.INT_NULL]);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
      DG.Column.fromInt32Array('Route', routeArr),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
        route: {name: 'Route', required: false, rawType: 'int32'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const scalars: ProfileScalars = {...baseScalars(), route: nca.ROUTE_PO};
    const {inputs} = extractProfileInputs(df, index, contract, scalars);
    expect(inputs.route, nca.ROUTE_PO, 'null route cell must fall back to scalars.route');
  });

  test('routeStringToCode covers all 6 sci-comp codes; throws on unknown', async () => {
    expect(routeStringToCode('IV-bolus'), nca.ROUTE_IV_BOLUS);
    expect(routeStringToCode('IV-infusion'), nca.ROUTE_IV_INFUSION);
    expect(routeStringToCode('PO'), nca.ROUTE_PO);
    expect(routeStringToCode('SC'), nca.ROUTE_SC);
    expect(routeStringToCode('IM'), nca.ROUTE_IM);
    expect(routeStringToCode('other'), nca.ROUTE_OTHER);
    let threw = false;
    try {
      routeStringToCode('subcu');
    } catch (e) {
      threw = e instanceof Error && e.message.includes('subcu');
    }
    expect(threw, true, 'unknown route string must throw');
  });

  // ---------------------------------------------------------------
  // nca-unit-semantics — extractor delegates rescaling to
  // scaleProfileInputs (AD-1 / AD-12 single-rescaling-point invariant)
  // ---------------------------------------------------------------

  test('canonical output: extractor rescales min→h and ng/mL→mg/L', async () => {
    // Source: time in minutes, conc in ng/mL.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 60, 120])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([1000, 500, 250])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const scalars: ProfileScalars = {
      ...baseScalars(), timeUnits: 'min', concentrationUnits: 'ng/mL',
    };
    const {inputs} = extractProfileInputs(df, index, contract, scalars);
    expect(inputs.timeUnits, 'h', 'output must carry canonical h');
    expect(inputs.concentrationUnits, 'mg/L', 'output must carry canonical mg/L');
    expect(inputs.doseUnits, 'mg', 'output must carry canonical mg');
    expectArray(Array.from(inputs.time), [0, 1, 2]);
    expectArray(Array.from(inputs.conc), [1, 0.5, 0.25]);
  });

  test('mg/kg dose + BW=70 + bwUnit=kg → canonical 70 mg dose', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const scalars: ProfileScalars = {
      ...baseScalars(), dose: 1, doseUnits: 'mg/kg', bodyWeight: 70,
    };
    const {inputs} = extractProfileInputs(df, index, contract, scalars);
    expect(Math.abs(inputs.dose - 70) < 1e-12, true,
      `mg/kg × 70 kg ⇒ 70 mg; got ${inputs.dose}`);
    expect(inputs.bodyWeight, 70);
  });

  // ---------------------------------------------------------------
  // blank-BW-dosenorm-fix — a mapped BW COLUMN with per-subject blanks
  // must classify each profile independently: a blank cell → null BW →
  // per_kg (NOT absolute with a NaN dose). Regression for the bug where
  // readNumericCell's FLOAT_NULL→NaN leaked past scaleProfileInputs'
  // `bodyWeight === null` check, mislabeling DOSE_NORM and NaN-ing the dose.
  // ---------------------------------------------------------------

  test('mg/kg + mapped BW column, blank cell → per_kg with finite dose (no absolute/NaN)', async () => {
    // S1 has BW=250 g; S2's BW cell is blank (null). Same mg/kg dose scalar.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S2', 'S2']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 10, 5])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'BW', [250, 250, null, null]),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
        bodyWeight: {name: 'BW', required: false, rawType: 'float64'},
      },
    };
    const indices = buildProfileIndices(df, contract);
    expect(indices.length, 2, 'two subjects');
    const scalars: ProfileScalars = {
      ...baseScalars(), dose: 5, doseUnits: 'mg/kg', bodyWeightUnit: 'g',
    };

    // S1 — BW present (250 g = 0.25 kg) → absolute, dose = 5 mg/kg × 0.25 kg = 1.25 mg.
    const s1 = extractProfileInputs(df, indices[0], contract, scalars);
    expect(s1.doseNormalization, 'absolute', 'present-BW subject is absolute');
    expect(Math.abs(s1.inputs.dose - 1.25) < 1e-9, true,
      `present BW → 5 mg/kg × 0.25 kg = 1.25 mg; got ${s1.inputs.dose}`);

    // S2 — BW blank → null → per_kg, dose stays per-kg (5 mg/kg, no ×BW), finite.
    const s2 = extractProfileInputs(df, indices[1], contract, scalars);
    expect(s2.doseNormalization, 'per_kg', 'blank-BW subject must be per_kg, not absolute');
    expect(s2.inputs.bodyWeight, null, 'blank BW cell must resolve to null, not NaN');
    expect(Number.isFinite(s2.inputs.dose), true,
      `blank-BW dose must be finite per-kg, not NaN; got ${s2.inputs.dose}`);
    expect(Math.abs(s2.inputs.dose - 5) < 1e-9, true,
      `per-kg dose = 5 mg/kg (no ×BW); got ${s2.inputs.dose}`);
    expect(s2.warnings.some((w) => w.code === 'PER_KG_NO_BW'), true,
      `blank-BW subject must surface PER_KG_NO_BW; got ${JSON.stringify(s2.warnings)}`);
  });

  test('BW grams 250 g with bwUnit=g → canonical 0.25 kg', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const scalars: ProfileScalars = {
      ...baseScalars(), bodyWeight: 250, bodyWeightUnit: 'g',
    };
    const {inputs} = extractProfileInputs(df, index, contract, scalars);
    expect(Math.abs((inputs.bodyWeight ?? NaN) - 0.25) < 1e-12, true,
      `250 g → 0.25 kg; got ${inputs.bodyWeight}`);
  });

  test('mouse 30 with bwUnit=kg surfaces BW_UNIT_HEURISTIC_WARNING in extractor warnings', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const scalars: ProfileScalars = {
      ...baseScalars(), dose: 1, doseUnits: 'mg/kg', bodyWeight: 30, bodyWeightUnit: 'kg',
    };
    const {warnings} = extractProfileInputs(df, index, contract, scalars);
    expect(warnings.some((w) => w.code === 'BW_UNIT_HEURISTIC_WARNING'), true,
      `expected BW_UNIT_HEURISTIC_WARNING in extractor warnings; got ${JSON.stringify(warnings)}`);
  });

  test('canonical units present even when source units are already canonical', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 2])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 2.5])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const {inputs} = extractProfileInputs(df, index, contract, baseScalars());
    expect(inputs.timeUnits, 'h');
    expect(inputs.concentrationUnits, 'mg/L');
    expect(inputs.doseUnits, 'mg');
    expectArray(Array.from(inputs.time), [0, 1, 2]);
    expectArray(Array.from(inputs.conc), [10, 5, 2.5]);
  });

  test('LLOQ scalar rescaled by same factor as concentration (AD-7)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([1000, 500])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'string'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const [index] = buildProfileIndices(df, contract);
    const scalars: ProfileScalars = {
      ...baseScalars(), concentrationUnits: 'ng/mL', lloq: 2,
    };
    const {inputs} = extractProfileInputs(df, index, contract, scalars);
    expect(typeof inputs.lloq, 'number');
    expect(Math.abs((inputs.lloq as number) - 0.002) < 1e-15, true,
      `lloq 2 ng/mL → 0.002 mg/L (same factor as conc); got ${inputs.lloq}`);
  });
});
