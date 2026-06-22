import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';

import './tests/rich-function-view-tests';
import './tests/rich-function-view-direct-tests';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
