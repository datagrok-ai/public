import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {searchSubstructure} from '../package';

export async function loadFileAsText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}

export async function loadFileAsBytes(name: string): Promise<Uint8Array> {
  return await _package.files.readAsBytes(name);
}

export async function dfFromColWithOneCategory(colName: string, value: string, length: number): Promise<DG.DataFrame> {
  const col = DG.Column.fromType(DG.COLUMN_TYPE.STRING, colName, length);
  col.init((i) => value);
  return DG.DataFrame.fromColumns([col]);
}
export async function createTableView(tableName: string): Promise<DG.TableView> {
  const df = await readDataframe(tableName);
  df.name = tableName.replace('.csv', '');
  await grok.data.detectSemanticTypes(df);
  const view = grok.shell.addTableView(df);
  return view;
}

export async function readDataframe(tableName: string): Promise<DG.DataFrame> {
  const file = await loadFileAsText(tableName);
  const df = DG.DataFrame.fromCsv(file);
  df.name = tableName.replace('.csv', '');
  return df;
}

export async function _testSearchSubstructure(df: DG.DataFrame, colName: string,
  pattern: string, trueIndices: number[]): Promise<void> {
  const col = df.columns.byName(colName);
  const bitset: DG.BitSet = (await searchSubstructure(col, pattern, '')).get(0);
  const bitsetString = bitset.toBinaryString();
  const bitsetArray = [...bitsetString];
  for (let k = 0; k < trueIndices.length; ++k) {
    expect(bitsetArray[trueIndices[k]] === '1', true);
    bitsetArray[trueIndices[k]] = '0';
  }
  for (let i = 0; i < bitsetArray.length; ++i)
    expect(bitsetArray[i] === '0', true);
}

export async function _testSearchSubstructureSARSmall(): Promise<void> {
  const file = await loadFileAsText('sar-small.csv');

  const df = DG.DataFrame.fromCsv(file);
  const col = df.columns.byIndex(0);
  const bitset: DG.BitSet = (await searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', '')).get(0);
  const countDataframe = col.length;
  const countResult = bitset.trueCount;
  expect(countDataframe, 200);
  expect(countResult, 90);
}

export async function _testSearchSubstructureAllParameters(foo: any): Promise<void> {
  await foo({});
  await foo();
}
