import * as DG from 'datagrok-api/dg';
import {runTests, tests} from '@datagrok-libraries/utils/src/test';

import './tests/drugbank-tests';

export const _package = new DG.Package();
export {tests}

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
export async function test(category: string, test: string): Promise<DG.DataFrame> {
  const data = await runTests({category, test});
  return DG.DataFrame.fromObjects(data)!;
}
