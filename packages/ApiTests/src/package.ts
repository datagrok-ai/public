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

//name: expectTable
//shortName: expectTable
//input: dataframe actual
//input: dataframe expected
//output: bool result
export function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): boolean {
  _expectTable(actual, expected);
  return true;
}

//name: getOutput
//input: string func
//input: string param
//output: map res
export async function getOutput(func: string, param: string): Promise<object> {
  const f = await grok.functions.eval(func);
  const params: {[key: string]: any} = {};
  f.inputs.forEach((i: any) => params[i.name] = typeof i.defaultValue === 'string' ?
    i.defaultValue.slice(1, -1) : i.defaultValue);
  const call = await f.prepare(params);
  await call.call();
  const res = {val: call.getParamValue(param)};
  return res;
}
