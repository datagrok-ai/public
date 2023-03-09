/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TestViewerForProperties} from './viewers/test-viewer-for-properties';
import {expectTable as _expectTable} from '@datagrok-libraries/utils/src/test';
// import {hashDataFrame} from '@datagrok-libraries/utils/src/dataframe-utils';


//name: TestViewerForProperties
//description: Viewer to test properties and others
//tags: viewer, panel
//output: viewer result
export function testViewerForProperties() {
  return new TestViewerForProperties();
}

//name: getTable
//input: string name
//input: string path {optional: true}
//output: dataframe result
export async function getTable(name: string, path: string): Promise<DG.DataFrame> {
  const file = (await grok.dapi.files.list(path ? `system:appdata/${path}/` : 'Demo:Files/', true, name))[0];
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

//name: expectTable
//shortName: expectTable
//input: dataframe actual
//input: dataframe expected
//output: bool result
export function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): boolean {
  _expectTable(actual, expected);
  return true;
}
