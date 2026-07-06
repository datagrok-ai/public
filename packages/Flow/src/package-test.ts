import { runTests, tests, TestContext , initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';

import './tests/type-map-tests';
import './tests/node-factory-tests';
import './tests/compiler-tests';
import './tests/serializer-tests';
import './tests/creation-script-import-tests';
import './tests/panel-tests';
import './tests/layout-tests';
import './tests/order-edge-tests';
import './tests/minimap-tests';
import './tests/creation-script-emit-tests';
import './tests/function-browser-tests';
import './tests/inspect-tests';
import './tests/test-ids-tests';
import './tests/guide-tests';
import './tests/summary-tests';
import './tests/files-tree-tests';
import './tests/execution-preview-tests';
import './tests/viewer-tests';
import './tests/column-picker-tests';
import './tests/connect-interaction-tests';
import './tests/string-list-tests';
import './tests/rerun-node-tests';
import './tests/func-editor-tests';
import './tests/entity-tests';
import './tests/editor-bridge-tests';

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
