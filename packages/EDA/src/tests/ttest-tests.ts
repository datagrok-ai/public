// Tests for the two-sample t-test.

import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {twoSampleTTest} from '../ttest/ttest-tools';
import {oneWayAnova} from '../anova/anova-tools';

const ALPHA = 0.05;

/** Build a string category column from integer-encoded labels.
 *  The 'g' prefix is required: DG.Column.fromStrings auto-promotes purely numeric
 *  strings to an INT column, which breaks the factorization indexing. */
function factorCol(name: string, data: readonly number[]): DG.Column {
  return DG.Column.fromStrings(name, data.map((n) => `g${n}`));
}

/** Build a double-precision value column. */
function valueCol(name: string, data: readonly number[]): DG.Column<DG.COLUMN_TYPE.FLOAT> {
  return DG.Column.fromFloat64Array(name, Float64Array.from(data)) as unknown as DG.Column<DG.COLUMN_TYPE.FLOAT>;
}

/** Approximate-equality check with relative + absolute tolerances. */
function close(actual: number, expected: number, label: string, rtol = 1e-4, atol = 1e-6): void {
  const diff = Math.abs(actual - expected);
  const threshold = atol + rtol * Math.abs(expected);
  expect(diff <= threshold, true,
    `${label}: |${actual} - ${expected}| = ${diff} > ${threshold} (rtol=${rtol}, atol=${atol})`);
}

category('t-test: two-sample', () => {
  // Group g0 = [1..5] (mean 3, var 2.5), g1 = [6..10] (mean 8, var 2.5); n0 = n1 = 5.
  const equalCats = factorCol('group', [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]);
  const equalVals = valueCol('feature', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);

  test('Student: known values', async () => {
    const r = twoSampleTTest(equalCats, equalVals, equalCats.categories.length,
      {method: 'Student', alpha: ALPHA});

    expect(r.method, 'Student');
    close(r.meanDiff, 5, 'meanDiff');
    close(r.t, 5, 't');
    close(r.df, 8, 'df');
    close(r.pValue, 0.0010520, 'pValue', 1e-2);
    close(r.ciLow, 2.693996, 'ciLow', 1e-3);
    close(r.ciHigh, 7.306004, 'ciHigh', 1e-3);
    close(r.cohenD, 3.1622777, 'cohenD', 1e-4);
    close(r.hedgesG, 2.8565, 'hedgesG', 1e-3);
  }, {timeout: 10000});

  test('Welch on equal variances == Student df', async () => {
    // With equal variances and equal sizes, Welch's df reduces to n0 + n1 − 2.
    const r = twoSampleTTest(equalCats, equalVals, equalCats.categories.length,
      {method: 'Welch', alpha: ALPHA});

    expect(r.method, 'Welch');
    close(r.meanDiff, 5, 'meanDiff');
    close(r.t, 5, 't');
    close(r.df, 8, 'df (Welch reduces to 8)');
    close(r.cohenD, 3.1622777, 'cohenD');
  }, {timeout: 10000});

  test('t² = F identity (Student vs Fisher ANOVA)', async () => {
    const r = twoSampleTTest(equalCats, equalVals, equalCats.categories.length,
      {method: 'Student', alpha: ALPHA});
    const anova = oneWayAnova(equalCats, equalVals, ALPHA, {method: 'Fisher', toValidate: false});
    if (anova.method !== 'Fisher')
      throw new Error('Expected Fisher report');

    close(anova.anovaTable.fStat, r.t * r.t, 'F == t²');
    expect(anova.anovaTable.dfBn, 1, 'Fisher dfBn must be 1 for k=2');
    expect(anova.anovaTable.dfWn, 8, 'Fisher dfWn');
    close(anova.anovaTable.pValue, r.pValue, 'p-value (Fisher == Student t)');
  }, {timeout: 10000});

  test('Welch: unequal variances + t² = W identity', async () => {
    // g0 = [1..5] (var 2.5), g1 = [10,20,30,40,50] (var 250); fractional Welch df.
    const cats = factorCol('group', [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]);
    const vals = valueCol('feature', [1, 2, 3, 4, 5, 10, 20, 30, 40, 50]);

    const r = twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});

    close(r.meanDiff, 27, 'meanDiff');
    close(r.t, 3.799408, 't', 1e-4);
    close(r.df, 4.079995, 'df (Satterthwaite, fractional)', 1e-3);
    close(r.cohenD, 2.402992, 'cohenD (d_s, non-pooled SD)', 1e-4);
    close(r.hedgesG, 2.170459, 'hedgesG (g_s*)', 1e-3);

    // Welch one-way ANOVA with k=2 must reproduce the Welch t-test: W = t², matching df and p.
    const anova = oneWayAnova(cats, vals, ALPHA, {method: 'Welch', toValidate: false});
    if (anova.method !== 'Welch')
      throw new Error('Expected Welch report');
    close(anova.anovaTable.fStat, r.t * r.t, 'W == t²');
    close(anova.anovaTable.dfWn, r.df, 'Welch dfWn == t-test df');
    close(anova.anovaTable.pValue, r.pValue, 'p-value (Welch ANOVA == Welch t)');
  }, {timeout: 10000});

  test('Sign of mean difference (group 1 − group 0)', async () => {
    const cats = factorCol('group', [0, 0, 0, 1, 1, 1]);
    const vals = valueCol('feature', [10, 11, 12, 1, 2, 3]);
    const r = twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});

    close(r.meanDiff, -9, 'meanDiff must be negative');
    expect(r.t < 0, true, 't must be negative when group 1 < group 0');
  }, {timeout: 10000});

  test('Guard: more than two groups throws', async () => {
    const cats = factorCol('group', [0, 0, 1, 1, 2, 2]);
    const vals = valueCol('feature', [1, 2, 3, 4, 5, 6]);

    let threw = false;
    try {
      twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});
    } catch (e) {
      threw = true;
      expect((e as Error).message.includes('exactly 2'), true, 'message should mention exactly 2 groups');
    }
    expect(threw, true, 'three groups must throw');
  }, {timeout: 10000});

  test('Guard: degenerate group (n < 2) throws', async () => {
    const cats = factorCol('group', [0, 1, 1, 1]);
    const vals = valueCol('feature', [1, 2, 3, 4]);

    let threw = false;
    try {
      twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});
    } catch (e) {
      threw = true;
      expect((e as Error).message.includes('2 observations'), true, 'message should mention 2 observations');
    }
    expect(threw, true, 'a singleton group must throw');
  }, {timeout: 10000});
});
