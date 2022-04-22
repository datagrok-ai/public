import * as DG from 'datagrok-api/dg';
import {runTests, tests} from '@datagrok-libraries/utils/src/test';

import './tests/sequence-validation';

export const _package = new DG.Package();
export {tests};

//name: testOligoBatchCalculator
//output: dataframe result
export async function testOligoBatchCalculator(): Promise<DG.DataFrame> {
  const data = await runTests();
  return DG.DataFrame.fromObjects(data)!;
}
