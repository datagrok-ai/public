import {before, category, expect, expectTable, test} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';
import {default as init} from "parquet-wasm";
import {fromFeather, fromParquet} from "../api/api";

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

category('Feather', () => {
  let dfFromArrow: DG.DataFrame | null;

  before(async () => {
    const bytesArrow = await _package.files.readAsBytes('titanic.feather');
    dfFromArrow = fromFeather(bytesArrow);
  });

  test('fromFeather: column names', async () => {
    for (const colName of expectedColumns)
      expect(dfFromArrow?.columns.contains(colName), true);
  });

  test('fromFeather: column types', async () => {
    expect(dfFromArrow?.getCol('pclass').type, DG.COLUMN_TYPE.INT);
    expect(dfFromArrow?.getCol('name').type, DG.COLUMN_TYPE.STRING);
    expect(dfFromArrow?.getCol('age').type, DG.COLUMN_TYPE.FLOAT);
  });

  test('fromFeather: number of rows and columns', async () => {
    expect(dfFromArrow?.columns.length, 14);
    expect(dfFromArrow?.rowCount, 1311);
  });

  test('fromFeather: serialization', async () => {
    expectTable(dfFromArrow!, DG.DataFrame.fromByteArray((dfFromArrow ?? DG.DataFrame.create()).toByteArray()));
  });
});
