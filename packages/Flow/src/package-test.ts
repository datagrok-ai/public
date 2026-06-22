import { runTests, tests, TestContext , initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';

import './tests/type-map-tests';
import './tests/node-factory-tests';
import './tests/compiler-tests';
import './tests/serializer-tests';
import './tests/creation-script-import-tests';
import './tests/panel-tests';
import './tests/layout-tests';

export let _package = new DG.Package();
export { tests };

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({ category, test, testContext });
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
