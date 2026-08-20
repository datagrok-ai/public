import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import './tests/versioning-tests';
import './tests/clone-tests';
import './tests/curation-tests';
import './tests/security-tests';
import './tests/integration-tests';
import './tests/platform-mechanics-tests';
import './tests/function-run-tests';
import './tests/ui-forms-tests';
import './tests/restricted-publish-tests';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool excludeNodeTests {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext,
  excludeNodeTests?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, excludeNodeTests});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
