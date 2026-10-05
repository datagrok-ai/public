// Tests for missing values imputation

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {MetricInfo, DISTANCE_TYPE, impute} from '../missing-values-imputation/knn-imputer';
import {getFeatureInputSettings, imputeColumns} from '../missing-values-imputation/impute-columns';
import {dataWithMissingVals} from './utils';

const ROWS_K = 100;
const K = 1000;
const INT_COLS = 5;
const FLOAT_COLS = 5;
const STRING_COLS = 5;
const MIS_VALS_COUNT = 5;
const NEIGHBORS = 5;
const TIMEOUT = 10000;
const TOTAL_COLS = INT_COLS + FLOAT_COLS + STRING_COLS;

const testKNN = (dist: DISTANCE_TYPE) => {
  test(`${dist} dist, ${ROWS_K}K rows, ${TOTAL_COLS} cols, ${MIS_VALS_COUNT * TOTAL_COLS} missing vals`, async () => {
    // Data
    const data = dataWithMissingVals(ROWS_K * K, INT_COLS, FLOAT_COLS, STRING_COLS, MIS_VALS_COUNT);
    const df = data.df;
    const cols = df.columns;

    // Inputs for kNN imputer
    const targetColNames = cols.names();
    const featuresMetrics = new Map<string, MetricInfo>();
    const missingValsIndices = data.misValsIds;

    // Imputation settings
    for (const col of df.columns) {
      const settings = getFeatureInputSettings(col.type as DG.COLUMN_TYPE);
      featuresMetrics.set(col.name, {
        weight: settings.defaultWeight,
        type: settings.defaultMetric,
      });
    }

    // Impute missing values & get fails
    const failedToImput = impute(df, targetColNames, featuresMetrics, missingValsIndices, dist, NEIGHBORS, true);

    // Check fails
    let fails = 0;
    failedToImput.forEach((inds, _) => fails += inds.length);
    expect(fails, 0, `Failed to impute ${fails} missing values`);
  }, {timeout: TIMEOUT, benchmark: true});
};

category(`Missing values imputation`, () => {
  testKNN(DISTANCE_TYPE.EUCLIDEAN);
  testKNN(DISTANCE_TYPE.MANHATTAN);

  test('Impute columns in place', async () => {
    const rows = 30;
    const allNullRow = 12;
    const values = (f: (i: number) => number, nullRows: number[]) =>
      Array.from({length: rows}, (_, i) => nullRows.includes(i) ? null : f(i));
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', values((i) => i, [3, allNullRow])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', values((i) => 2 * i, [7, allNullRow])),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'c', values((i) => rows - i, [allNullRow])),
    ]);
    const nullRows = (name: string) => Array.from({length: rows}, (_, i) => i).filter((i) => df.getCol(name).isNone(i));

    const failed = imputeColumns(df, ['a', 'b'], ['a', 'b', 'c'], 4, DISTANCE_TYPE.EUCLIDEAN);

    expect(df.columns.names().join(), 'a,b,c');
    expect(df.get('a', 3), 3);
    expect(df.get('b', 7), 14);
    expect(nullRows('a').join(), `${allNullRow}`);
    expect(nullRows('b').join(), `${allNullRow}`);
    expect(nullRows('c').join(), `${allNullRow}`);
    expect([...failed.keys()].join(), 'a,b');
    expect(failed.get('a')?.join(), `${allNullRow}`);
    expect(failed.get('b')?.join(), `${allNullRow}`);
  });
});
