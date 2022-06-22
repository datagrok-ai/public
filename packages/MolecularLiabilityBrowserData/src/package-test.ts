import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {runTests} from '@datagrok-libraries/utils/src/test';

export const _package = new DG.Package();

//name: test
//output: dataframe result
export async function test() {
  const data = await runTests();
  return DG.DataFrame.fromObjects(data);
}
