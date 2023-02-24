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
  const bitset: DG.BitSet = (await searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)c2ccccc2N1', '')).get(0);
  const countDataframe = col.length;
  const countResult = bitset.trueCount;
  expect(countDataframe, 200);
  expect(countResult, 90);
}

export async function _testSearchSubstructureAllParameters(foo: any): Promise<void> {
  await foo({});
  await foo();
}

export const molV2000 = `NVS-LGIVUK
RDKit          2D 0   0.00000     0.00000     0

19 20  0  0  0  0              1 V2000 
-0.6933    2.3553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
0.4388    1.3714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-0.6581    0.3483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-2.0926    0.7867    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
-3.1896   -0.2364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-3.0066   -1.7252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-4.3660   -2.3593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-5.3891   -1.2623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
-4.6620    0.0497    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
-5.2961    1.4091    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
1.2991    2.6002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
2.7336    2.1617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
2.7599    0.6620    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
3.9887   -0.1983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
5.3481    0.4358    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
3.8581   -1.6926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
5.0869   -2.5528    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
2.4987   -2.3266    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
1.3417    0.1735    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
1  2  1  0  0  0  0
2  3  1  0  0  0  0
3  4  1  0  0  0  0
4  5  1  0  0  0  0
5  6  1  0  0  0  0
6  7  1  0  0  0  0
7  8  1  0  0  0  0
8  9  1  0  0  0  0
9 10  2  0  0  0  0
2 11  1  0  0  0  0
11 12  1  0  0  0  0
12 13  1  0  0  0  0
13 14  1  0  0  0  0
14 15  2  0  0  0  0
14 16  1  0  0  0  0
16 17  1  0  0  0  0
16 18  1  0  0  0  0
13 19  1  0  0  0  0
19  2  1  0  0  0  0
9  5  1  0  0  0  0
M  END`;

export const molV3000 = `Azaguanine-8
Actelion Java MolfileCreator 2.0

  0  0  0  0  0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 12 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0 0.0 0.0 0
M  V30 2 C -0.0119 -1.5014 0.0 0
M  V30 3 N -1.3227 -2.2522 0.0 0
M  V30 4 C 1.275 -2.2522 0.0 0
M  V30 5 C -1.3346 -3.7418 0.0 0
M  V30 6 N 2.705 -1.7755 0.0 0
M  V30 7 C 1.2631 -3.7418 0.0 0
M  V30 8 N -0.0476 -4.4806 0.0 0
M  V30 9 N -2.6812 -4.4806 0.0 0
M  V30 10 N 3.5749 -2.9791 0.0 0
M  V30 11 N 2.6931 -4.1946 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 3 5
M  V30 5 1 4 6
M  V30 6 2 4 7
M  V30 7 2 5 8
M  V30 8 1 5 9
M  V30 9 2 6 10
M  V30 10 1 7 11
M  V30 11 1 7 8
M  V30 12 1 10 11
M  V30 END BOND
M  V30 END CTAB
M  END`;
