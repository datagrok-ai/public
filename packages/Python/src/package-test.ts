import { runTests, tests, TestContext , initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';
import { DataFrame } from 'datagrok-api/dg';

export let _package = new DG.Package();
export { tests };

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DataFrame> {
  const data = await runTests({ category, test, testContext });
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
