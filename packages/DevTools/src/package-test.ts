import * as DG from 'datagrok-api/dg';
// import * as grok from 'datagrok-api/grok';
import {runTests, TestContext, tests, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';

import './tests/function-signature-editor-test';
// import './tests/dev-panel-test';
import './tests/test-manager-tests';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool skipCore {optional: true}
//input: bool verbose {optional: true}
//output: dataframe result
export async function test(category: string, test: string,
  testContext: TestContext, skipCore: boolean = false, verbose = true): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, exclude: skipCore ? ['Core'] : undefined, verbose});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
