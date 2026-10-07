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
import './tests/new-data-tests';
import './tests/persona-tests';
import './tests/program-role-tests';
import './tests/evidence-tests';
import './tests/snapshot-tests';
import './tests/document-tests';
import './tests/provenance-tests';
import './tests/feed-tests';
import './tests/prediction-tests';
import './tests/approval-tests';
import './tests/parking-tests';
import './tests/evidence-table-tests';
import './tests/pinning-tests';
import './tests/flow-folder-tests';
import './tests/data-package-tests';
import './tests/compare-runs-tests';
import './tests/observation-tests';
import './tests/design-set-tests';
import './tests/dmta-tests';
import './tests/panel-tests';
import './tests/write-grant-tests';
import './tests/program-row-tests';
import './tests/empty-filter-tests';
import './tests/find-tests';
import './tests/decision-record-tests';
import './tests/data-access-tests';
import './tests/read-scope-tests';

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
