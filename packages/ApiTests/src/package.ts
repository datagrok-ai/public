/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TestViewerForProperties} from './viewers/test-viewer-for-properties';
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

//name: compareTables
//input: dataframe actual
//input: dataframe expected
//output: bool result
export function compareTables(actual: DG.DataFrame, expected: DG.DataFrame): boolean {
  const result = actual.toCsv() === expected.toCsv();
  return result;
  /*
  let result: boolean = true;
  const uint8ArrayActual = hashDataFrame(actual);
  const uint8ArrayExpected = hashDataFrame(expected);
  if (uint8ArrayExpected.length !== uint8ArrayActual.length) {
    result = false;
    return result;
  }
  for (let i = 0; i < uint8ArrayExpected.length; i++) {
    if (uint8ArrayExpected[i] !== uint8ArrayActual[i]) {
      result = false;
      break;
    }
  }
  return result;
  */
}
