import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {tableFromIPC} from 'apache-arrow'
import { Buffer } from 'buffer';
import { _package } from '../package-test';
import { default as init, readParquet, writeParquet, WriterPropertiesBuilder, Compression } from '../arrow1';


category('Arrow', () => {
  let bytes;
  let arrow;
  let table;

  before(async () => {
    bytes = await _package.files.readAsBytes('titanic.feather');
    arrow = Buffer.from(bytes);
    table = tableFromIPC(arrow);
  });

  test('feather file column names', async () => {
    const expectedTables = ['pclass', 'survived', 'name', 'sex', 'age', 
                            'sibsp', 'parch', 'ticket', 'fare', 'cabin',
                            'embarked', 'boat', 'body', 'home.dest'
                          ];
    expect(table.schema.fields.map((table) => table.name).toString(), expectedTables.toString());
  });

  test('feather file column types', async () => {
    expect(table.getChild('pclass').type.toString(), 'Int32');
    expect(table.getChild('name').type.toString(), 'Dictionary<Int32, Utf8>');
    expect(table.getChild('age').type.toString(), 'Float64');
  });

  test('feather file number of rows and columns', async () => {
    expect(table.numCols, 14);
    expect(table.numRows, 1310);
  });

});
