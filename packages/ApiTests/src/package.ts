/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {expect, expectTable as _expectTable} from '@datagrok-libraries/utils/src/test';
import dayjs from "dayjs";

export const _package = new DG.Package();

//name: getTable
//input: string name
//input: string path {optional: true}
//output: dataframe result
export async function getTable(name: string, path: string): Promise<DG.DataFrame> {
  const file = (await grok.dapi.files.list(path ? `system:appdata/${path}/` : 'System:DemoFiles/', true, name))[0];
  const str = await file.readAsString();
  const result = DG.DataFrame.fromCsv(str);
  return result;
}

//name: getColumn
//input: dataframe table
//input: string columnName
//output: column col
export function getColumn(table: DG.DataFrame, columnName: string): DG.Column {
  const col = table.getCol(columnName);
  return col;
}

//name: getDT [get demo table]
//input: int rows {optional: true}
//input: string name {optional: true}
//output: dataframe df
export function getDT(rows: number = 20, name: any = 'demog'): DG.DataFrame {
  const df = grok.data.demo.getDemoTable(name, rows);
  return df;
}

//name: getCell [get table cell]
//input: dataframe table
//input: int rowIndex
//input: string columnName
//output: cell x
export function getCell(table: DG.DataFrame, rowIndex: number, columnName: string): DG.Cell {
  return table.cell(rowIndex, columnName);
}

//name: expectTable
//input: dataframe actual
//input: dataframe expected
//output: bool result
export function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): boolean {
  _expectTable(actual, expected);
  return true;
}

//name: dummyPackageFunctionWithDefaultValue
//input: string a = "test"
//output: int c
export function dummyPackageFunctionWithDefaultValue(a: string) {
  const c = a.length;
  return c;
}

//name: dummyPackageFunction
//input: int a
//input: int b
//output: int c
export function dummyPackageFunction(a: number, b: number) {
  const c = a + b;
  return c;
}

//name: dummyDataFrameFunction
//input: dataframe table [Data table]
//output: int count 
//output: dataframe tableOut
export function dummyDataFrameFunction(table: DG.DataFrame) {
  return {'tableOut': table, 'count': table.rowCount};
}

//name: testStringAsync
//input: int a
//output: int b
export async function testIntAsync(a: number): Promise<number> {
  return new Promise((r) => r(a + 10));
}

//name: CustomStringInput
//input: object params
//output: object input
export function CustomStringInput(params: any) {
  const defaultInput = ui.input.string('Custom input', {value: ''});
  defaultInput.root.style.backgroundColor = 'aqua';
  defaultInput.input.style.backgroundColor = 'aqua';
  return defaultInput;
}

//name: testOutputAnnotationJoinDf
//input: dataframe data
//input: column col
//output: dataframe res {action:join(data)}
export function testOutputAnnotationJoinDf(data: DG.DataFrame, col: DG.Column<string>): DG.DataFrame {
  const colRes = DG.Column.string('joined', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return DG.DataFrame.fromColumns([colRes]);
}

//name: testOutputAnnotationJoinCol
//input: dataframe data
//input: column col
//output: column res {action:join(data)}
export function testOutputAnnotationJoinCol(data: DG.DataFrame, col: DG.Column<string>): DG.Column {
  const colRes = DG.Column.string('joined', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return colRes;
}

//name: testOutputAnnotationJoinColList
//input: dataframe data
//input: column col
//output: column res {action:join(data)}
export function testOutputAnnotationJoinColList(data: DG.DataFrame, col: DG.Column<string>): DG.ColumnList {
  const colRes = DG.Column.string('joined', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return DG.DataFrame.fromColumns([colRes]).columns;
}

//name: testOutputAnnotationReplaceDf
//input: dataframe data
//input: column col
//output: dataframe res {action:replace(data)}
export function testOutputAnnotationReplaceDf(data: DG.DataFrame, col: DG.Column<string>): DG.DataFrame {
  const colRes = DG.Column.string('val', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return DG.DataFrame.fromColumns([colRes]);
}

//name: testOutputAnnotationReplaceCol
//input: dataframe data
//input: column col
//output: column res {action:replace(data)}
export function testOutputAnnotationReplaceCol(data: DG.DataFrame, col: DG.Column<string>): DG.Column {
  const colRes = DG.Column.string('val', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return colRes;
}

//name: testOutputAnnotationReplaceColList
//input: dataframe data
//input: column col
//output: column res {action:replace(data)}
export function testOutputAnnotationReplaceColList(data: DG.DataFrame, col: DG.Column<string>): DG.ColumnList {
  const colRes = DG.Column.string('val', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return DG.DataFrame.fromColumns([colRes]).columns;
}

//name: testOutputWithoutAction
//input: dataframe data
//input: column col
//output: column res
export function testOutputWithoutAction(data: DG.DataFrame, col: DG.Column<string>): DG.Column {
  const res = DG.Column.string('val', data.rowCount).init((i) => `${col.get(i)}_abc`);
  return res;
}

//name: expectDate
//input: datetime actual
//input: datetime expected
export function expectDate(actual: dayjs.Dayjs, expected: dayjs.Dayjs): void {
  expect(actual.valueOf(), expected.valueOf());
}


//name: testVectorFunc
//meta.vectorFunc: true
//input: column col
//input: string prefix
//output: column res
export function testVectorFunc(col: DG.Column, prefix: string) {
  const res = DG.Column.string('result', col.length).init((i) => `${prefix}_${col.get(i)}`);
  return res;
}


//name: testVectorFuncNonVectorizableParam
//meta.vectorFunc: true
//input: column col
//input: string prefix
//input: column postfix
//output: column res
export function testVectorFuncNonVectorizableParam(col: DG.Column, prefix: string, postfix: DG.Column) {
  const res = DG.Column.string('result', col.length).init((i) => `${prefix}_${col.get(i)}_${postfix.get(i)}`);
  return res;
}