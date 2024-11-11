import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import '@datagrok-libraries/bio/src/types/ngl'; // To enable import from the NGL module declared in bio lib
import {runTests, tests, TestContext, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';

import './tests/autodock-tests';

export let _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool stressTest {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext, stressTest?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, stressTest});
  return DG.DataFrame.fromObjects(data)!;
}


//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
