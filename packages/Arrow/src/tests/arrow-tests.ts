// Feather round-trip tests have moved to LibTests
// (packages/LibTests/src/tests/compute-utils/fitting/arrow-{roundtrip,titanic}.ts)
// since the conversion code now lives in @datagrok-libraries/arrow.
// Parquet tests stay here because parquet-wasm's WASM init is bound to this
// package's webRoot.

import {before, category, expect, expectTable, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import {default as init} from "parquet-wasm";
import {fromParquet} from "../api/api";

const expectedColumns = ['pclass', 'survived', 'name', 'sex', 'age',
  'sibsp', 'parch', 'ticket', 'fare', 'cabin',
  'embarked', 'boat', 'body', 'home.dest'
];

category('Parquet', () => {
  let dfFromParquet: DG.DataFrame | null;

  before(async () => {
    await init(_package.webRoot + 'dist/parquet_wasm_bg.wasm');
    const bytesParquet = await _package.files.readAsBytes('titanic.parquet');
    dfFromParquet = fromParquet(bytesParquet);
  });

  test('fromParquet: column names', async () => {
    for (const colName of expectedColumns)
      expect(dfFromParquet?.columns.contains(colName), true);
  });

  test('fromParquet: column types', async () => {
    expect(dfFromParquet?.getCol('pclass').type, DG.COLUMN_TYPE.INT);
    expect(dfFromParquet?.getCol('name').type, DG.COLUMN_TYPE.STRING);
    expect(dfFromParquet?.getCol('age').type, DG.COLUMN_TYPE.FLOAT);
  });

  test('fromParquet: number of rows and columns', async () => {
    expect(dfFromParquet?.columns.length, 14);
    expect(dfFromParquet?.rowCount, 1311);
  });
  
  test('fromParquet: serialization', async () => {
    expectTable(dfFromParquet!, DG.DataFrame.fromByteArray((dfFromParquet ?? DG.DataFrame.create()).toByteArray()));
  });
});

// Regression tests for null fidelity and Date32 handling on the Parquet read
// path. Fixture null-fidelity.parquet (3 rows, null in row 1 of every column):
//   i32 int32 [1, null, 3] · i64 int64 [10, null, 30] · i64big int64 [2^40, null, -2^40]
//   f64 double [1.5, null, 2.5] · ts timestamp[ms] · d32 date32[day]
//   cat dictionary<utf8> · b bool · s utf8
category('Parquet: null fidelity', () => {
  let df: DG.DataFrame | null;

  before(async () => {
    await init(_package.webRoot + 'dist/parquet_wasm_bg.wasm');
    const bytes = await _package.files.readAsBytes('null-fidelity.parquet');
    df = fromParquet(bytes);
  });

  test('nulls are nulls, not zeros', async () => {
    for (const colName of ['i32', 'i64', 'i64big', 'f64', 'ts', 'd32', 'cat', 's'])
      expect(df?.getCol(colName).isNone(1), true, `${colName}[1] should be null`);
  });

  test('non-null values survive around nulls', async () => {
    expect(df?.getCol('i32').get(0), 1);
    expect(df?.getCol('i32').get(2), 3);
    expect(df?.getCol('i64').get(0), 10);
    expect(df?.getCol('f64').get(0), 1.5);
    expect(df?.getCol('cat').get(2), 'b');
  });

  test('int64 beyond int32 range keeps values and nulls', async () => {
    const c = df?.getCol('i64big');
    expect(c?.type, DG.COLUMN_TYPE.BIG_INT);
    expect(c?.get(0)?.toString(), '1099511627776');
    expect(c?.get(2)?.toString(), '-1099511627776');
  });

  test('date32 converts day units, not raw day counts', async () => {
    const d = df?.getCol('d32');
    expect(d?.type, DG.COLUMN_TYPE.DATE_TIME);
    // 2020-01-15: raw value is 18276 days; a broken path yields a 1970 epoch-adjacent date
    expect(new Date(d?.get(0)?.valueOf()).getUTCFullYear(), 2020);
  });

  test('timestamp nulls do not become epoch', async () => {
    const t = df?.getCol('ts');
    expect(t?.isNone(1), true);
    expect(new Date(t?.get(0)?.valueOf()).getUTCFullYear(), 2020);
  });
});

