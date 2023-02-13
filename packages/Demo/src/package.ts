/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoView} from './demo-app/demo-app';
// eslint-disable-next-line no-unused-vars
import {expect, expectArray, expectFloat, expectObject} from '@datagrok-libraries/utils/src/test';


export const _package = new DG.Package();


//name: Demo
//tags: app
//description: Interactive demo of major Datagrok capabilities
export function demoApp() {
  grok.shell.addView(new DemoView());
}

//name: getTable
//input: string name
//input: string folder {optional: true}
//output: dataframe result
export async function getTable(name: string, folder: string = ''): Promise<DG.DataFrame> {
  const file = (await grok.dapi.files.list(folder ? `Demo:Files/${folder}/` : 'Demo:Files/', true, name))[0];
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

//name: getColumns
//input: dataframe table
//input: string_list columnNames
//output: column_list cols
export function getColumns(table: DG.DataFrame, columnNames: string[]): DG.Column[] {
  const cols = columnNames.map((c) => table.getCol(c));
  return cols;
}

//name: compareObjects
//input: map actual
//input: map expected
//output: bool equals
export function compareObjects(actual: {[key: string]: any}, expected: {[key: string]: any}): boolean {
  expectObject(actual, expected);
  const equals = true;
  return equals;
}
