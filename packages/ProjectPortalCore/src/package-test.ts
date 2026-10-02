import { runTests, tests, TestContext , initAutoTests as initTests } from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import './tests/model-tests';
import './tests/portal-tests';
import './tests/schema-tests';
import './tests/readiness-tests';
import './tests/work-unit-tests';
import './tests/workspace-growth-tests';
import './tests/compound-tests';
import './tests/portfolio-table-tests';
import './tests/model-run-tests';
import './tests/catalogue-tests';
import './tests/save-run-tests';
import './tests/run-selection-tests';
import './tests/program-data-tests';
import './tests/attention-tests';
import './tests/sync-run-tests';
import './tests/auto-advance-tests';
import './tests/profile-tests';
import './tests/recommend-tests';

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
