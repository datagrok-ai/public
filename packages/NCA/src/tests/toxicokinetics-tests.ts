/**
 * Adapter — toxicokinetics (TK) aggregates (tk-aggregates AD-1..AD-4).
 *
 * @see files/rat_tk28d_nca_per_cell.csv          (real PKNCA per-cell NCA — the input)
 * @see files/rat_tk28d_accumulation_ratios.csv   (regenerated oracle, AC-1/AC-2)
 * @see files/rat_tk28d_exposure_margin.csv       (regenerated oracle, AC-9/AC-10)
 *
 * Validation tests (AC-1/AC-2/AC-9/AC-10) build the primitive's input frame from
 * the **real** PKNCA per-cell NCA output (`rat_tk28d_nca_per_cell.csv`, reparsed
 * to Float64 to avoid the BUG-01 Float32 precision loss) and assert the TS
 * primitive reproduces the **independently** R-computed oracle CSVs — real output
 * vs reference, not a mirror-pattern (rule 16 / rule 18). Unit tests (AC-3..AC-7,
 * AC-11..AC-12) drive synthetic frames to exercise gate logic deterministically.
 *
 * **Footgun gate:** `expect(cond === expected, true, msg)` form throughout — never
 * `expect(x, undefined, msg)` (guardrails).
 */

import {expectNoRawNan} from './test-utils';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';

import {
  computeAccumulationRatios,
  accumulationIndexTheoretical,
  computeExposureMargin,
  computeSexFold,
} from '../adapter/toxicokinetics';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {reparseFloatColumnsToFloat64} from '../adapter/csv-float64';
import {unmappedGridColumns} from '../views/results-view';

const PER_CELL = 'System:AppData/Nca/rat_tk28d_nca_per_cell.csv';
const ACC_REF = 'System:AppData/Nca/rat_tk28d_accumulation_ratios.csv';
const MARGIN_REF = 'System:AppData/Nca/rat_tk28d_exposure_margin.csv';

/** Load a CSV fixture as a Float64-reparsed DataFrame (BUG-01 precision). */
async function loadCsv(path: string): Promise<DG.DataFrame> {
  const csv = await grok.dapi.files.readAsText(path);
  const df = DG.DataFrame.fromCsv(csv);
  reparseFloatColumnsToFloat64(df, csv);
  return df;
}

/** Per-cell accumulation options pinning the per_cell.csv column names. */
const ACC_OPTS = {
  tau: 24,
  doseCol: 'Dose_mg_kg',
  studyDayCol: 'StudyDay',
  stratifyCols: ['Sex'] as const,
  cmaxCol: 'cmax',
  aucTauCol: 'auclast',
  aucInfCol: 'aucinf.obs',
  lambdaZCol: 'lambda.z',
  adjR2Col: 'adj.r.squared',
  pctExtrapCol: 'aucpext.obs',
  tlastCol: 'tlast',
};

const MARGIN_OPTS = {
  noaelMgKg: 100,
  humanDoseMg: 200,
  humanBodyWeightKg: 70,
  clRefDoseMgKg: 10,
  ratBodyWeightKg: 0.25,
  doseCol: 'Dose_mg_kg',
  studyDayCol: 'StudyDay',
  stratifyCols: ['Sex'] as const,
  aucTauCol: 'auclast',
  clCol: 'cl.obs',
};

category('Adapter: toxicokinetics — accumulation (AD-1/AD-2)', () => {
  let perCell: DG.DataFrame;
  let accRef: DG.DataFrame;
  before(async () => {
    perCell = await loadCsv(PER_CELL);
    accRef = await loadCsv(ACC_REF);
  });

  // BUG-12 drift guard: every column of the REAL accumulation table must resolve
  // through `tkColumnDisplay` (fixed keys + the sex-fold prefixes), or
  // configureTkGrid leaves it with a raw header and DG's generic stats tooltip.
  // Only the caller's stratify columns are exempt (user-named, self-describing).
  test('BUG-12: every accumulation-table column has a presentation entry', async () => {
    const {table} = computeAccumulationRatios(perCell, ACC_OPTS);
    // AC-13 (F4 Arm 1): this producer's REAL output must carry no raw IEEE NaN.
    // Asserted here rather than in a duplicate fixture — the frame already exists.
    expectNoRawNan(table, 'toxicokinetics.accumulationTable');
    const unmapped = unmappedGridColumns(table, 'tk', [...ACC_OPTS.stratifyCols]);
    expect(unmapped.length, 0,
      `TK columns unresolved by tkColumnDisplay: [${unmapped.join(', ')}]`);
  });

  test('AC-1/AC-2: RAC_CMAX + RAC_AUC (Formula A) reproduce the regenerated oracle (8/8)', async () => {
    const {table, day1, dayN} = computeAccumulationRatios(perCell, ACC_OPTS);
    expect(day1, 1, 'day1 resolves to min StudyDay = 1');
    expect(dayN, 28, 'dayN resolves to max StudyDay = 28');
    expect(table.rowCount, 8, `expected 8 (Dose×Sex) cells; got ${table.rowCount}`);

    // Index computed output by Sex|Dose.
    const idx = new Map<string, number>();
    for (let i = 0; i < table.rowCount; i++)
      idx.set(`${table.col('Sex')!.get(i)}|${table.col('DOSE')!.get(i)}`, i);

    let checked = 0;
    for (let r = 0; r < accRef.rowCount; r++) {
      const sex = String(accRef.col('Sex')!.get(r));
      const dose = accRef.col('Dose_mg_kg')!.get(r) as number;
      const expCmax = accRef.col('Cmax_acc_ratio')!.get(r) as number;
      const expAuc = accRef.col('AUC_acc_ratio')!.get(r) as number;
      const row = idx.get(`${sex}|${dose}`);
      expect(row !== undefined, true, `computed row for ${sex}/${dose} missing`);
      // Same source NCA numbers → tight tolerance (pure pivot/ratio).
      expectFloat(table.col('RAC_CMAX')!.get(row!) as number, expCmax, Math.abs(expCmax) * 1e-6,
        `RAC_CMAX ${sex}/${dose}`);
      expectFloat(table.col('RAC_AUC')!.get(row!) as number, expAuc, Math.abs(expAuc) * 1e-6,
        `RAC_AUC ${sex}/${dose}`);
      checked++;
    }
    expect(checked, 8, 'all 8 oracle cells checked');
    // Spot-pin the headline cell (synthesis AC-1/AC-2).
    const f10 = idx.get('F|10')!;
    expectFloat(table.col('RAC_AUC')!.get(f10) as number, 5.425, 0.01, '10/F Rac(AUC) = 5.425');
    expectFloat(table.col('RAC_CMAX')!.get(f10) as number, 4.618, 0.01, '10/F Rac(Cmax) = 4.618');
    // Default formula label is A.
    expect(String(table.col('RAC_FORMULA')!.get(f10)), 'A (AUCτ/AUCτ)', '10/F formula = A');
  });

  test('AC-3: Formula B on 06 → all 8 cells fall back to A (label asserts fallback)', async () => {
    const {table} = computeAccumulationRatios(perCell, {...ACC_OPTS, formula: 'B'});
    expect(table.rowCount, 8, '8 cells');
    for (let i = 0; i < table.rowCount; i++) {
      const label = String(table.col('RAC_FORMULA')!.get(i));
      expect(label.indexOf('fallback') >= 0, true,
        `cell ${table.col('Sex')!.get(i)}/${table.col('DOSE')!.get(i)} must fall back to A; got '${label}'`);
      // The fallback value still equals Formula A (AUCτ/AUCτ) — never NaN.
      expect(Number.isFinite(table.col('RAC_AUC')!.get(i) as number), true,
        'fallback RAC_AUC is finite (= Formula A), not NaN');
    }
  });

  test('AC-6: RAC_THEO = NaN for all 8 cells (every cell fails the λz gate)', async () => {
    const {table} = computeAccumulationRatios(perCell, ACC_OPTS);
    for (let i = 0; i < table.rowCount; i++) {
      expect(table.col('RAC_THEO')!.isNone(i), true,
        `cell ${table.col('Sex')!.get(i)}/${table.col('DOSE')!.get(i)} RAC_THEO must be NaN (gate fails)`);
      // SS flag is empty when theo is unavailable (honest uncertainty).
      expect(String(table.col('RAC_SS_FLAG')!.get(i)), '', 'RAC_SS_FLAG empty when theo unavailable');
    }
    // The column carries the "Accumulation Index" label (tooltip surface).
    expect(table.col('RAC_THEO')!.tags['label'], 'Accumulation Index (theoretical)',
      'RAC_THEO labelled Accumulation Index');
  });

  test('AC-7: synthetic pass-cell → RAC_THEO equals the closed form 1/(1−e^(−λz·τ))', async () => {
    // λz present, adj.R² = 0.95, extrap = 10%, τ = 24 → gate passes.
    const lambdaZ = 0.05;
    const tau = 24;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Sex', ['F', 'F']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([10, 10])),
      DG.Column.fromInt32Array('StudyDay', Int32Array.from([1, 28])),
      DG.Column.fromFloat64Array('cmax', Float64Array.from([2, 8])),
      DG.Column.fromFloat64Array('auclast', Float64Array.from([40, 200])),
      DG.Column.fromFloat64Array('aucinf.obs', Float64Array.from([50, 250])),
      DG.Column.fromFloat64Array('lambda.z', Float64Array.from([lambdaZ, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('adj.r.squared', Float64Array.from([0.95, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('aucpext.obs', Float64Array.from([10, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('tlast', Float64Array.from([24, 24])),
    ]);
    const {table} = computeAccumulationRatios(df, {...ACC_OPTS, tau});
    expect(table.rowCount, 1, 'one cell');
    const expectedTheo = 1 / (1 - Math.exp(-lambdaZ * tau));
    expectFloat(table.col('RAC_THEO')!.get(0) as number, expectedTheo, Math.abs(expectedTheo) * 1e-9,
      'RAC_THEO = closed form');
    // Direct helper parity.
    expectFloat(accumulationIndexTheoretical(lambdaZ, tau), expectedTheo, 1e-12, 'helper = closed form');
    // SS flag set (observed Rac 5.0 vs theo ~1.45 ⇒ observed ≫ theoretical).
    expect(String(table.col('RAC_SS_FLAG')!.get(0)).length > 0, true, 'SS flag text set on pass-cell');
  });

  test('AC-4: Day-1 tlast ≠ τ → that cell RAC_AUC = NaN + RAC_TLAST_NE_TAU; others unaffected', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Sex', ['F', 'F', 'M', 'M']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([10, 10, 10, 10])),
      DG.Column.fromInt32Array('StudyDay', Int32Array.from([1, 28, 1, 28])),
      DG.Column.fromFloat64Array('cmax', Float64Array.from([2, 8, 2, 8])),
      DG.Column.fromFloat64Array('auclast', Float64Array.from([40, 200, 40, 200])),
      // F Day-1 stops at 8 h (≠ τ=24); M samples to 24 h.
      DG.Column.fromFloat64Array('tlast', Float64Array.from([8, 24, 24, 24])),
    ]);
    const {table, warnings} = computeAccumulationRatios(df, {
      tau: 24, doseCol: 'Dose_mg_kg', studyDayCol: 'StudyDay', stratifyCols: ['Sex'],
      cmaxCol: 'cmax', aucTauCol: 'auclast', tlastCol: 'tlast',
    });
    const idx = new Map<string, number>();
    for (let i = 0; i < table.rowCount; i++) idx.set(String(table.col('Sex')!.get(i)), i);
    expect(table.col('RAC_AUC')!.isNone(idx.get('F')!), true,
      'F cell RAC_AUC = NaN (tlast ≠ τ)');
    expect(Number.isFinite(table.col('RAC_AUC')!.get(idx.get('M')!) as number), true,
      'M cell RAC_AUC finite (tlast == τ)');
    expect(warnings.some((w) => w.code === 'RAC_TLAST_NE_TAU'), true, 'RAC_TLAST_NE_TAU raised');
  });

  test('AC-5: a cell mixing CL_BASIS → Rac refused (NaN) + MIXED_CL_BASIS_AGGREGATE', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Sex', ['F', 'F']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([10, 10])),
      DG.Column.fromInt32Array('StudyDay', Int32Array.from([1, 28])),
      DG.Column.fromFloat64Array('cmax', Float64Array.from([2, 8])),
      DG.Column.fromFloat64Array('auclast', Float64Array.from([40, 200])),
      DG.Column.fromFloat64Array('tlast', Float64Array.from([24, 24])),
      DG.Column.fromStrings('CL_BASIS', ['CL', 'CL/F']),
    ]);
    const {table, warnings} = computeAccumulationRatios(df, {
      tau: 24, doseCol: 'Dose_mg_kg', studyDayCol: 'StudyDay', stratifyCols: ['Sex'],
      cmaxCol: 'cmax', aucTauCol: 'auclast', tlastCol: 'tlast',
    });
    expect(table.col('RAC_AUC')!.isNone(0), true, 'mixed-basis cell RAC_AUC is EMPTY');
    expect(warnings.some((w) => w.code === 'MIXED_CL_BASIS_AGGREGATE'), true,
      'MIXED_CL_BASIS_AGGREGATE raised (not divided)');
  });

  test('AC-5: a cell mixing DOSE_NORM → Rac refused (NaN) + MIXED_DOSE_NORM_AGGREGATE', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Sex', ['F', 'F']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([10, 10])),
      DG.Column.fromInt32Array('StudyDay', Int32Array.from([1, 28])),
      DG.Column.fromFloat64Array('cmax', Float64Array.from([2, 8])),
      DG.Column.fromFloat64Array('auclast', Float64Array.from([40, 200])),
      DG.Column.fromFloat64Array('tlast', Float64Array.from([24, 24])),
      DG.Column.fromStrings('DOSE_NORM', ['absolute', 'per_kg']),
    ]);
    const {table, warnings} = computeAccumulationRatios(df, {
      tau: 24, doseCol: 'Dose_mg_kg', studyDayCol: 'StudyDay', stratifyCols: ['Sex'],
      cmaxCol: 'cmax', aucTauCol: 'auclast', tlastCol: 'tlast',
    });
    expect(table.col('RAC_AUC')!.isNone(0), true, 'mixed dose-norm cell RAC_AUC is EMPTY');
    expect(warnings.some((w) => w.code === 'MIXED_DOSE_NORM_AGGREGATE'), true,
      'MIXED_DOSE_NORM_AGGREGATE raised (not divided)');
  });
});

category('Adapter: toxicokinetics — exposure margin (AD-3)', () => {
  let perCell: DG.DataFrame;
  let marginRef: DG.DataFrame;
  before(async () => {
    perCell = await loadCsv(PER_CELL);
    // The card's cl.obs is PKNCA CL on a mg/kg dose with no weight — per kg, L/h/kg. The margin
    // refuses a CL that does not say so (nca-cl-bw-normalised-units D-m), so the test declares it.
    perCell.col('cl.obs')!.setTag('units', 'L/h/kg');
    marginRef = await loadCsv(MARGIN_REF);
  });

  test('AC-8/AC-9: per-sex MARGIN_HED reproduces oracle (sex-identical); MARGIN_EXPOSURE NaN', async () => {
    const {table, warnings} = computeExposureMargin(perCell, MARGIN_OPTS);
    expect(table !== null, true, 'margin table produced');
    expect(table!.rowCount, 2, 'one row per sex');
    expect(warnings.some((w) => w.code === 'NO_HUMAN_AUC_INPUT'), true, 'AC-8: NO_HUMAN_AUC_INPUT (DG-1)');

    const idx = new Map<string, number>();
    for (let i = 0; i < table!.rowCount; i++) idx.set(String(table!.col('Sex')!.get(i)), i);
    for (let r = 0; r < marginRef.rowCount; r++) {
      const sex = String(marginRef.col('Sex')!.get(r));
      const row = idx.get(sex)!;
      // AC-8: exposure-based default NaN (no empirical human AUC in 06).
      expect(table!.col('MARGIN_EXPOSURE')!.isNone(row), true,
        `${sex} MARGIN_EXPOSURE = NaN`);
      // AC-9: HED matches oracle (dose-based, sex-identical).
      const expHed = marginRef.col('MARGIN_HED')!.get(r) as number;
      expectFloat(table!.col('MARGIN_HED')!.get(row) as number, expHed, Math.abs(expHed) * 1e-6,
        `${sex} MARGIN_HED`);
    }
    // NI-3: both sexes carry the identical HED value (no sex term).
    expectFloat(table!.col('MARGIN_HED')!.get(idx.get('F')!) as number,
      table!.col('MARGIN_HED')!.get(idx.get('M')!) as number, 1e-9, 'MARGIN_HED sex-identical');
    expectFloat(table!.col('MARGIN_HED')!.get(idx.get('F')!) as number, 5.67568, 0.001,
      'MARGIN_HED ≈ 5.68× (NOAEL 100 × 6/37 ÷ (200/70))');
  });

  test('AC-9: MARGIN_PROJECTED reproduces oracle — F valid, M NaN (10/M CL is NA)', async () => {
    const {table, warnings} = computeExposureMargin(perCell, MARGIN_OPTS);
    const idx = new Map<string, number>();
    for (let i = 0; i < table!.rowCount; i++) idx.set(String(table!.col('Sex')!.get(i)), i);
    const oracle = new Map<string, number>();
    for (let r = 0; r < marginRef.rowCount; r++)
      oracle.set(String(marginRef.col('Sex')!.get(r)), marginRef.col('MARGIN_PROJECTED')!.get(r) as number);

    // F: finite, matches oracle (~10.78).
    const fProj = table!.col('MARGIN_PROJECTED')!.get(idx.get('F')!) as number;
    expectFloat(fProj, oracle.get('F')!, Math.abs(oracle.get('F')!) * 1e-6, 'F MARGIN_PROJECTED');
    expectFloat(fProj, 10.784, 0.01, 'F projected margin ≈ 10.78');
    // M: NaN (10 mg/kg Day-1 cl.obs is NA → no valid reference CL).
    expect(table!.col('MARGIN_PROJECTED')!.isNone(idx.get('M')!), true,
      'M MARGIN_PROJECTED = NaN (no valid Day-1 CL)');
    expect(warnings.some((w) => w.code === 'NO_CL_REFERENCE_FOR_SEX'), true, 'NO_CL_REFERENCE_FOR_SEX for M');
    // Worst-case = F (the only finite, hence smallest, projection).
    expect(table!.col('WORST_CASE')!.get(idx.get('F')!) as boolean, true, 'F flagged worst-case');
  });

  test('AC-10: allometric default; disabling allometry inflates the margin ≈4.1× + warning', async () => {
    const allometric = computeExposureMargin(perCell, MARGIN_OPTS);
    const nonAllo = computeExposureMargin(perCell, {...MARGIN_OPTS, disableAllometry: true});
    const fA = allometric.table!.col('Sex')!;
    const findF = (t: DG.DataFrame): number => {
      for (let i = 0; i < t.rowCount; i++) if (String(t.col('Sex')!.get(i)) === 'F') return i;
      throw new Error('F row missing');
    };
    void fA;
    const aF = allometric.table!.col('MARGIN_PROJECTED')!.get(findF(allometric.table!)) as number;
    const nF = nonAllo.table!.col('MARGIN_PROJECTED')!.get(findF(nonAllo.table!)) as number;
    // (70/0.25)^0.25 = 4.091 → non-allometric inflates the margin by ~4.1×.
    expectFloat(nF / aF, Math.pow(70 / 0.25, 0.25), 0.01, 'non-allometric / allometric ≈ 4.09×');
    expect(allometric.table!.col('PROJECTION_ALLOMETRIC')!.get(findF(allometric.table!)) as boolean, true,
      'default projection is allometric');
    expect(nonAllo.table!.col('PROJECTION_ALLOMETRIC')!.get(findF(nonAllo.table!)) as boolean, false,
      'disabled projection is non-allometric');
    expect(nonAllo.warnings.some((w) => w.code === 'MARGIN_NON_ALLOMETRIC'), true,
      'MARGIN_NON_ALLOMETRIC inflation warning raised');
    expect(allometric.warnings.some((w) => w.code === 'MARGIN_NON_ALLOMETRIC'), false,
      'allometric default carries no inflation warning');
  });

  // --- nca-cl-bw-normalised-units D-m: the margin's CL input contract (AC-6.2 / 6.3 / 6.4). ---

  test('AC-6.4: no BW on the card → labelled 0.25 kg default; MARGIN_PROJECTED = 10.784 unchanged', async () => {
    const {table} = computeExposureMargin(perCell, MARGIN_OPTS);
    const f = [...Array(table!.rowCount).keys()].find((i) => table!.col('Sex')!.get(i) === 'F')!;
    expectFloat(table!.col('MARGIN_PROJECTED')!.get(f) as number, 10.784, 0.01, 'F projected margin');
    expect(table!.col('RAT_REF_BW_KG')!.get(f), 0.25, 'the default anchor, as a number');
    expect(table!.col('RAT_REF_BW_SOURCE')!.get(f), 'default', 'and labelled as the default');
  });

  test('AC-6.2: a CL column that is not per kg (absolute or untagged) is refused, never read as per kg', async () => {
    for (const unit of ['L/h', '']) {
      const df = perCell.clone();
      df.col('cl.obs')!.setTag('units', unit);
      const {table, warnings} = computeExposureMargin(df, MARGIN_OPTS);
      expect(warnings.some((w) => w.code === 'MARGIN_CL_NOT_PER_KG'), true, `refused for units '${unit}'`);
      for (let i = 0; i < table!.rowCount; i++)
        expect(table!.col('MARGIN_PROJECTED')!.isNone(i), true, `no projected margin for units '${unit}'`);
    }
  });

  test('AC-6.2: mL/min/kg CL is converted to L/h/kg before the ratified formula', async () => {
    const df = perCell.clone();
    const cl = df.col('cl.obs')!;
    for (let i = 0; i < df.rowCount; i++) if (!cl.isNone(i)) cl.set(i, (cl.get(i) as number) * 1000 / 60);
    cl.setTag('units', 'mL/min/kg');
    const {table} = computeExposureMargin(df, MARGIN_OPTS);
    const f = [...Array(table!.rowCount).keys()].find((i) => table!.col('Sex')!.get(i) === 'F')!;
    expectFloat(table!.col('MARGIN_PROJECTED')!.get(f) as number, 10.784, 0.01, 'same margin in mL/min/kg');
  });

  test('AC-6.3: a cohort body weight is the anchor (number + source); rows without CL are counted', async () => {
    // Synthetic F cell: three Day-1 10 mg/kg profiles, one without a CL; NOAEL AUC on Day 28.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Sex', ['F', 'F', 'F', 'F']),
      DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([10, 10, 10, 100])),
      DG.Column.fromInt32Array('StudyDay', Int32Array.from([1, 1, 1, 28])),
      DG.Column.fromFloat64Array('auclast', Float64Array.from([40, 41, 42, 2000])),
      DG.Column.fromFloat64Array('cl.obs', Float64Array.from([0.06, 0.07, DG.FLOAT_NULL, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('BW_KG', Float64Array.from([0.2, 0.24, 0.3, 0.25])),
    ]);
    df.col('cl.obs')!.setTag('units', 'L/h/kg');
    const {table, warnings} = computeExposureMargin(df, MARGIN_OPTS);
    // Anchor = mean BW over the rows that gave a CL (0.2, 0.24), not the CL-less one.
    expectFloat(table!.col('RAT_REF_BW_KG')!.get(0) as number, 0.22, 1e-12, 'cohort mean over the CL rows');
    expect(table!.col('RAT_REF_BW_SOURCE')!.get(0), 'cohort mean of 2 profiles', 'source states n');
    const dropped = warnings.find((w) => w.code === 'MARGIN_CL_ROWS_DROPPED');
    expect(dropped?.message.includes('2 of 3') ?? false, true, `n used of N stated: ${dropped?.message}`);
    // The ratified formula, recomputed in the test from its inputs (not through the function).
    const expected = 2000 / ((200 / 70) / (0.065 * Math.pow(70 / 0.22, -0.25)));
    expectFloat(table!.col('MARGIN_PROJECTED')!.get(0) as number, expected, expected * 1e-9, 'margin');
  });

  test('AC-11: no NOAEL input → NO_NOAEL_INPUT, no margin computed (no auto-derivation)', async () => {
    const {table, warnings} = computeExposureMargin(perCell, {...MARGIN_OPTS, noaelMgKg: undefined});
    expect(table === null, true, 'no margin table without NOAEL');
    expect(warnings.some((w) => w.code === 'NO_NOAEL_INPUT'), true, 'NO_NOAEL_INPUT raised');
  });
});

category('Adapter: toxicokinetics — sex fold (AD-4)', () => {
  test('AC-12: FOLD_CMAX = max/min of per-sex MEAN_CMAX; direction + N/%CV present; no significance', async () => {
    // Build a tiny per-subject frame: 2 F + 2 M, then aggregate by Sex.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('SEX', ['F', 'F', 'M', 'M']),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([10, 14, 20, 24])),
    ]);
    const agg = aggregateNcaResults(results, 'SEX', {params: ['CMAX']});
    const {table} = computeSexFold(agg, {sexCol: 'SEX', params: ['CMAX']});
    expect(table.rowCount, 1, 'one fold cell (no other group key)');

    const meanF = 12; // (10+14)/2
    const meanM = 22; // (20+24)/2
    expectFloat(table.col('FOLD_CMAX')!.get(0) as number, meanM / meanF, 1e-9, 'FOLD_CMAX = max/min');
    expect(String(table.col('FOLD_DIR_CMAX')!.get(0)), 'M>F', 'direction M>F (males higher)');
    // N and %CV the fold rests on are carried.
    expect(table.col('F_N_CMAX') !== null && table.col('M_N_CMAX') !== null, true, 'per-sex N carried');
    expect(table.col('F_CV_PCT_CMAX') !== null && table.col('M_CV_PCT_CMAX') !== null, true, 'per-sex %CV carried');
    // No significance / p-value column exists.
    let hasSig = false;
    for (const c of table.columns as Iterable<DG.Column>) {
      const n = c.name.toLowerCase();
      if (n.indexOf('pval') >= 0 || n.indexOf('p_value') >= 0 || n.indexOf('signif') >= 0) hasSig = true;
    }
    expect(hasSig, false, 'no significance/p-value column (composite N — thesis #4)');
  });
});
