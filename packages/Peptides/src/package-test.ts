import * as DG from 'datagrok-api/dg';
import {runTests, TestContext, tests, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';

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
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  // const helmInit = DG.Func.find({name: 'initHelm'})[0];
  // if (helmInit)
  //   await helmInit.apply();
  testContext?.catchUnhandled && (testContext.catchUnhandled = false);
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
