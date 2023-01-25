import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {Table, tableFromIPC} from 'apache-arrow'
import { Buffer } from 'buffer';
import { _package } from '../package-test';
//@ts-ignore
import { default as init, readParquet, writeParquet, WriterPropertiesBuilder, Compression } from '../arrow1';

category('Arrow', () => {
  let bytes_arrow;
  let arrow;
  let table_arrow: any;
  let bytes_parquet;
  let table_parquet: any;

  before(async () => {
    bytes_arrow = await _package.files.readAsBytes('titanic.feather');
    arrow = Buffer.from(bytes_arrow);
    table_arrow = tableFromIPC(arrow);
    await init(_package.webRoot + 'dist/arrow1_bg.wasm');
    bytes_parquet = await _package.files.readAsBytes('titanic.parquet');
    table_parquet = tableFromIPC(readParquet(bytes_parquet));
  });

  test('feather file column names', async () => {
    const expectedTables = ['pclass', 'survived', 'name', 'sex', 'age', 
                            'sibsp', 'parch', 'ticket', 'fare', 'cabin',
                            'embarked', 'boat', 'body', 'home.dest'
                          ];
    expect(table_arrow.schema.fields.map((table_arrow: any) => table_arrow.name).toString(), expectedTables.toString());
  });

  test('feather file column types', async () => {
    expect(table_arrow.getChild('pclass')?.type.toString(), 'Int32');
    expect(table_arrow.getChild('name')?.type.toString(), 'Dictionary<Int32, Utf8>');
    expect(table_arrow.getChild('age')?.type.toString(), 'Float64');
  });

  test('feather file number of rows and columns', async () => {
    expect(table_arrow.numCols, 14);
    expect(table_arrow.numRows, 1310);
  });
  
  test('parquet file column names', async () =>{
    const expectedTables = ['pclass', 'survived', 'name', 'sex', 'age', 
                            'sibsp', 'parch', 'ticket', 'fare', 'cabin',
                            'embarked', 'boat', 'body', 'home.dest'
                          ];
    expect(table_parquet.schema.fields.map((table_parquet: any) => table_parquet.name).toString(), expectedTables.toString());
  });

  test('parquet file column types', async () => {
    expect(table_parquet.getChild('pclass')?.type.toString(), 'Int32');
    expect(table_parquet.getChild('name')?.type.toString(), 'Dictionary<Int32, Utf8>');
    expect(table_parquet.getChild('age')?.type.toString(), 'Float64');
  });

  test('parquet file number of rows and columns', async () => {
    expect(table_parquet.numCols, 14);
    expect(table_parquet.numRows, 1310);
  });
});
