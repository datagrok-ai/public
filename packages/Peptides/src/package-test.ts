import * as DG from 'datagrok-api/dg';
import {runTests, TestContext, tests} from '@datagrok-libraries/utils/src/test';

import './tests/core';
import './tests/benchmarks';
import './tests/viewers';
import './tests/widgets';
import './tests/table-view';
import './tests/model';
import './tests/misc';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
/**
 *
 * @param category
 * @param test
 * @param testContext
 */
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}
