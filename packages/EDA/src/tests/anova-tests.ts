// Tests for ANOVA

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';

import {oneWayAnova, FactorizedData} from '../anova/anova-tools';

const ROWS_M = 1;
const M = 1000000;
const TIMEOUT = 4000;
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

category('ANOVA', () => {
  test(`Performance: ${ROWS_M}M rows demog`, async () => {
    const df = grok.data.demo.demog(ROWS_M * M);
    const categories = df.col(CATEGORIES);
    const features = df.col(FEATURES);

    const factorized = new FactorizedData(categories!, features!, categories!.stats.uniqueCount);
    factorized.areVarsEqual(ALPHA);

    oneWayAnova(categories!, features!, ALPHA, TO_VALIDATE);
  }, {timeout: TIMEOUT, benchmark: true});

  test(`Correctness`, async () => {
    const analysis = oneWayAnova(CATEGORIES_COL, FEATURES_COL, ALPHA, TO_VALIDATE);
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
});
