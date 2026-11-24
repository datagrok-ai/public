import * as DG from 'datagrok-api/dg';
import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';
import './tests/engines-tests';
import './tests/metrics-tests';
import './tests/training-tests';
import './tests/search-tests';
import './tests/preparation-tests';
import './tests/storage-tests';
import './tests/apply-tests';
import './tests/check-method-tests';
import './tests/catalog-tests';
import './tests/ui-tests';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: string skipToCategory {optional: true}
//input: string skipToTest {optional: true}
//input: bool returnOnFail {optional: true}
//input: list<string> exclude {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext, skipToCategory?: string,
  skipToTest?: string, returnOnFail?: boolean, exclude?: string[]): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, skipToCategory, skipToTest, returnOnFail, exclude});
  return DG.DataFrame.fromObjects(data) ?? DG.DataFrame.create();
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
