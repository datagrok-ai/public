import * as DG from 'datagrok-api/dg';

import {runTests, TestContext, tests} from '@datagrok-libraries/test/src/test';

import './tests/numbering-tests';
import './tests/performance-tests';
import './tests/search-tests';

export const _package = new DG.Package();
export {tests};

/** For the 'test' function argument names are fixed as 'category' and 'test' because of way it is called. */
//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool stressTest {optional: true}
//output: dataframe result
export async function test(
  category: string, test: string, testContext: TestContext, stressTest?: boolean,
): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, stressTest});
  return DG.DataFrame.fromObjects(data)!;
}
