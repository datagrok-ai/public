/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TestViewerForProperties} from './viewers/test-viewer-for-properties';

//name: TestViewerForProperties
//description: Viewer to test properties and others
//tags: viewer, panel
//output: viewer result
export function testViewerForProperties() {
  return new TestViewerForProperties();
}

//name: compareTables
//input: dataframe actual
//input: string name
//input: string path
//output: bool result
export async function compareTables(actual: DG.DataFrame, name:string, path: string): Promise<boolean> {
  const file = (await grok.dapi.files.list(path ? `system:appdata/${path}/` : 'Demo:Files/', true, name))[0];
  const expectedText = await file.readAsString();
  const actualText = actual.toCsv();
  return expectedText === actualText;
}

//name: getTable
//input: string name
//input: string path {optional: true}
//output: dataframe result
export async function getTable(name: string, path: string): Promise<DG.DataFrame> {
  const file = (await grok.dapi.files.list(path ? `system:appdata/${path}/` : 'Demo:Files/', true, name))[0];
  const str = await file.readAsString();
  return DG.DataFrame.fromCsv(str);
}
