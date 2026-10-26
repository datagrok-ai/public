import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {runTests, tests, TestContext, initAutoTests as initTests } from '@datagrok-libraries/test/src/test';
import {_package} from './package';
import './tests/campaigns';

export {_package, tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  // verbose: true - for tests returning dataframe
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
