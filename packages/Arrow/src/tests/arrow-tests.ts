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

