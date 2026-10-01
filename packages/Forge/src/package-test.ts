import * as DG from 'datagrok-api/dg';
import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';
import './tests/engines-tests';
import './tests/storage-tests';
import './tests/ui-tests';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data) ?? DG.DataFrame.create();
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
