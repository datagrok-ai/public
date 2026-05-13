// Titanic-based Feather round-trip tests for @datagrok-libraries/arrow.
// Parquet round-trip tests stay in packages/Arrow because parquet-wasm's
// WASM init is bound to that package's webRoot.

import * as DG from 'datagrok-api/dg';
import {before, category, expect, expectTable, test} from '@datagrok-libraries/test/src/test';
import {fromFeather, toFeather} from '@datagrok-libraries/arrow';
import {_package} from '../../package-test';

const expectedColumns = ['pclass', 'survived', 'name', 'sex', 'age',
  'sibsp', 'parch', 'ticket', 'fare', 'cabin',
  'embarked', 'boat', 'body', 'home.dest'];

category('Arrow / Feather titanic', () => {
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

  test('toFeather: titanic roundtrip', async () => {
    const bytes = toFeather(dfFromArrow!, true);
    expect(bytes !== null, true, 'toFeather returned null');
    const back = fromFeather(bytes!)!;
    expect(back.rowCount, dfFromArrow!.rowCount);
    expect(back.columns.length, dfFromArrow!.columns.length);
  });
});
