// Tests for ANOVA

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {oneWayAnova, FactorizedData} from '../anova/anova-tools';
import {FIXTURES, Tol} from './anova-fixtures';

const ROWS_M = 1;
const M = 1000000;
const TIMEOUT = 4000;
const NIST_TIMEOUT = 30000;
const ALPHA = 0.05;
const CATEGORIES = 'race';
const FEATURES = 'height';
const TO_VALIDATE = false;
const ERR = 0.01;

/** Validation features*/
const FEATURES_COL = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'features', [
  9, 12, 4, 8, 7, 4, 6, 8, 2, 10, 1, 3, 4, 5, 2,
]);

/** Validation categories */
const CATEGORIES_COL = DG.Column.fromStrings('features', [
  'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'C',
]);

/** Expected ANOVA results for the validation data */
enum EXPECTED {
  DF_BN = 2,
  DF_TOT = 14,
  DF_WN = 12,
  SS_BN = 63.333,
  SS_TOT = 147.333,
  SS_WN = 84,
  MS_BN = 31.666,
  MS_WN = 7,
  F_STAT = 4.523,
  F_CRIT = 3.885,
  P_VAL = 0.034,
};

/** Approximate-equality check with separate relative + absolute tolerances.
 *  Threshold = atol + rtol * |expected|. NaN==NaN is treated as a match. */
function expectClose(actual: number, expected: number, tol: Tol, label?: string): void {
  if (Number.isNaN(actual) && Number.isNaN(expected)) return;
  if (!Number.isFinite(actual) || !Number.isFinite(expected)) {
    expect(actual, expected, label ?? 'non-finite mismatch');
    return;
  }
  const diff = Math.abs(actual - expected);
  const threshold = tol.atol + tol.rtol * Math.abs(expected);
  expect(diff <= threshold, true,
    `${label ?? 'value'}: |${actual} - ${expected}| = ${diff} > ${threshold} ` +
    `(rtol=${tol.rtol}, atol=${tol.atol})`);
}

/** Build a string category column from integer-encoded factor labels.
 *  The 'g' prefix is required: DG.Column.fromStrings auto-promotes purely
 *  numeric strings to an INT column, which makes the underlying setStats
 *  reinterpret factor values as indices and run off the end of its buffers. */
function factorCol(name: string, data: readonly number[]): DG.Column {
  return DG.Column.fromStrings(name, data.map((n) => `g${n}`));
}

/** Build a double-precision (Float64) value column.
 *  Critical for NIST stiff datasets (SmLs07/SmLs09) where values are ~1e12 + 0.4:
 *  Float32 truncates those to a single constant and stdev collapses to 0. */
function valueCol(name: string, data: readonly number[]): DG.Column<DG.COLUMN_TYPE.FLOAT> {
  return DG.Column.fromFloat64Array(name, Float64Array.from(data)) as unknown as DG.Column<DG.COLUMN_TYPE.FLOAT>;
}

category('ANOVA', () => {
  test(`Performance: ${ROWS_M}M rows demog`, async () => {
    const df = grok.data.demo.demog(ROWS_M * M);
    const categories = df.col(CATEGORIES);
    const features = df.col(FEATURES);

    const factorized = new FactorizedData(categories!, features!, categories!.stats.uniqueCount);
    factorized.areVarsEqual(ALPHA);

    oneWayAnova(categories!, features!, ALPHA, {method: 'Fisher', toValidate: TO_VALIDATE});
  }, {timeout: TIMEOUT, benchmark: true});

  test(`Correctness`, async () => {
    const analysis = oneWayAnova(CATEGORIES_COL, FEATURES_COL, ALPHA, {
      method: 'Fisher',
      toValidate: TO_VALIDATE,
    });
    if (analysis.method !== 'Fisher')
      throw new Error('Expected Fisher report');
    const anova = analysis.anovaTable;

    // check degrees of freedom (df-s)
    expect(anova.dfBn, EXPECTED.DF_BN, 'Incorrect degrees of freedom: dfBn');
    expect(anova.dfTot, EXPECTED.DF_TOT, 'Incorrect degrees of freedom: dfTot');
    expect(anova.dfWn, EXPECTED.DF_WN, 'Incorrect degrees of freedom: dfWn');

    const eq = (x: number, y: number) => Math.abs(x - y) < ERR;

    // check sum of squares (ss-s)
    expect(eq(anova.ssBn, EXPECTED.SS_BN), true, 'Incorrect sum of squares: ssBn');
    expect(eq(anova.ssTot, EXPECTED.SS_TOT), true, 'Incorrect sum of squares: ssTot');
    expect(eq(anova.ssWn, EXPECTED.SS_WN), true, 'Incorrect sum of squares: ssWn');

    // check mean squares (ms-s)
    expect(eq(anova.msBn, EXPECTED.MS_BN), true, 'Incorrect mean squares: msBn');
    expect(eq(anova.msWn, EXPECTED.MS_WN), true, 'Incorrect mean squares: msWn');

    // check F-statistics
    expect(eq(anova.fStat, EXPECTED.F_STAT), true, 'Incorrect F-statistics value');

    // check p-value
    expect(eq(anova.pValue, EXPECTED.P_VAL), true, 'Incorrect p-value');

    // check F-critical
    expect(eq(analysis.fCritical, EXPECTED.F_CRIT), true, 'Incorrect F-critical');
  }, {timeout: TIMEOUT});

  test(`Correctness: Welch on validation data`, async () => {
    const analysis = oneWayAnova(CATEGORIES_COL, FEATURES_COL, ALPHA, {
      method: 'Welch',
      toValidate: false,
    });
    if (analysis.method !== 'Welch')
      throw new Error('Expected Welch report');
    const w = analysis.anovaTable;
    const f = FIXTURES.validation_simple.welch;

    expectClose(w.fStat, f.fStat, f.tol, 'Welch fStat');
    expect(w.dfBn, f.dfBn, 'Welch dfBn');
    expectClose(w.dfWn, f.dfWn, f.tol, 'Welch dfWn');
    expectClose(w.pValue, f.pValue, f.tol, 'Welch pValue');
    expectClose(analysis.fCritical, f.fCritical, f.tol, 'Welch fCritical');
  }, {timeout: TIMEOUT});

  // Canary tests — these protect the two numerical fixes. If somebody
  // reverts the Welford accumulator or brings back `1 - jStat.centralF.cdf(...)`,
  // these will fail.

  test('Numerical: variance stable on shifted data (Welford guard)', async () => {
    // Offset values by 1e9. Naive Σx² − (Σx)²/n catastrophically loses precision;
    // Welford preserves the variance. Float64 storage is required — Float32 (the
    // default for DG.COLUMN_TYPE.FLOAT) would itself truncate the offset.
    const offset = 1e9;
    const cats = DG.Column.fromStrings('cat', [
      'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',
    ]);
    const features = valueCol('feature', [
      offset + 0.0, offset + 1.0, offset + 2.0, offset + 1.0, offset + 0.5,
      offset + 10.0, offset + 11.0, offset + 12.0, offset + 11.0, offset + 10.5,
    ]);

    const r = oneWayAnova(cats, features, 0.05, {method: 'Welch', toValidate: false});
    if (r.method !== 'Welch')
      throw new Error('Expected Welch report');
    expect(Number.isFinite(r.anovaTable.fStat), true, 'F-statistic is NaN/Infinity');
    expect(r.anovaTable.fStat > 100, true, 'F-statistic too small — variance lost precision');
  }, {timeout: TIMEOUT});

  test('Numerical: p-value preserves precision in extreme tail (ibeta guard)', async () => {
    // Two extremely well-separated groups. With `1 - cdf` the cdf rounds to 1
    // and p collapses to 0 exactly; with ibeta-via-survival the p stays positive.
    const cats = DG.Column.fromStrings('cat', [
      'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
      'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B',
    ]);
    const features = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'feature', [
      0, 0.01, 0.02, 0.01, 0, 0.03, 0.02, 0.01, 0.02, 0.01,
      1000, 1000.1, 999.9, 1000.2, 999.8, 1000.1, 999.9, 1000.0, 1000.1, 999.9,
    ]);

    const r = oneWayAnova(cats, features, 0.05, {method: 'Fisher', toValidate: false});
    if (r.method !== 'Fisher')
      throw new Error('Expected Fisher report');
    expect(r.anovaTable.pValue > 0, true, 'p-value collapsed to 0 in tail');
    expect(r.anovaTable.pValue < 1e-30, true, 'p-value too large for this separation');
  }, {timeout: TIMEOUT});
});

category('ANOVA: NIST StRD', () => {
  for (const datasetName of ['SiRstv', 'AtmWtAg', 'SmLs04', 'SmLs07', 'SmLs09'] as const) {
    const fixture = FIXTURES.nist[datasetName];

    test(`Fisher: ${datasetName}`, async () => {
      const cats = factorCol('group', fixture.factor);
      const vals = valueCol('value', fixture.value);
      const r = oneWayAnova(cats, vals, 0.05, {method: 'Fisher', toValidate: false});
      if (r.method !== 'Fisher')
        throw new Error('Expected Fisher report');

      const f = fixture.fisher;
      expect(r.anovaTable.dfBn, f.dfBn, `${datasetName} Fisher dfBn`);
      expect(r.anovaTable.dfWn, f.dfWn, `${datasetName} Fisher dfWn`);
      expect(r.anovaTable.dfTot, f.dfTot, `${datasetName} Fisher dfTot`);
      expectClose(r.anovaTable.fStat, f.fStat, f.tol, `${datasetName} Fisher fStat`);
      expectClose(r.anovaTable.ssBn, f.ssBn, f.tol, `${datasetName} Fisher ssBn`);
      expectClose(r.anovaTable.ssWn, f.ssWn, f.tol, `${datasetName} Fisher ssWn`);
      expectClose(r.anovaTable.ssTot, f.ssTot, f.tol, `${datasetName} Fisher ssTot`);
      expectClose(r.anovaTable.msBn, f.msBn, f.tol, `${datasetName} Fisher msBn`);
      expectClose(r.anovaTable.msWn, f.msWn, f.tol, `${datasetName} Fisher msWn`);
      expectClose(r.fCritical, f.fCritical, f.tol, `${datasetName} Fisher fCritical`);
    }, {timeout: NIST_TIMEOUT});

    test(`Welch: ${datasetName}`, async () => {
      const cats = factorCol('group', fixture.factor);
      const vals = valueCol('value', fixture.value);
      const r = oneWayAnova(cats, vals, 0.05, {method: 'Welch', toValidate: false});
      if (r.method !== 'Welch')
        throw new Error('Expected Welch report');

      const w = fixture.welch;
      expect(r.anovaTable.dfBn, w.dfBn, `${datasetName} Welch dfBn`);
      expectClose(r.anovaTable.fStat, w.fStat, w.tol, `${datasetName} Welch fStat`);
      expectClose(r.anovaTable.dfWn, w.dfWn, w.tol, `${datasetName} Welch dfWn`);
      expectClose(r.fCritical, w.fCritical, w.tol, `${datasetName} Welch fCritical`);
    }, {timeout: NIST_TIMEOUT});
  }
});
