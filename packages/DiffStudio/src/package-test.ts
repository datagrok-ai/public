import * as DG from 'datagrok-api/dg';
import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/utils/src/test';
import './tests/numerical-methods-tests';
import './tests/features-tests';
import './tests/platform-funcs-tests';
import './tests/pipeline-tests';

export const _package = new DG.Package();
export {tests};


//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  // verbose: true - for tests returning dataframe
  const data = await runTests({category, test, testContext, verbose: true});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
