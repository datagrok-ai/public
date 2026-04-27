import * as DG from 'datagrok-api/dg';
import {TestContext, runTests, tests, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';

import './tests/compute-api/rich-function-view-tests';
import './tests/utils/expect-tests';
import './tests/utils/json-serialization-tests';
import './tests/compute-utils/rich-function-view-tests';
import './tests/compute-utils/reactive-tree-driver/config-processing';
import './tests/compute-utils/reactive-tree-driver/instance-init';
import './tests/compute-utils/reactive-tree-driver/instance-persistence';
import './tests/compute-utils/reactive-tree-driver/instance-mutations';
import './tests/compute-utils/reactive-tree-driver/instance-readonly';
import './tests/compute-utils/reactive-tree-driver/links-matching';
import './tests/compute-utils/reactive-tree-driver/links-reactivity-data';
import './tests/compute-utils/reactive-tree-driver/links-reactivity-validators';
import './tests/compute-utils/reactive-tree-driver/links-reactivity-meta';
import './tests/compute-utils/reactive-tree-driver/links-batching';
import './tests/compute-utils/reactive-tree-driver/links-reactivity-actions';
import './tests/compute-utils/reactive-tree-driver/instance-additional';
import './tests/compute-utils/reactive-tree-driver/funcall-wrappers';
import './tests/compute-utils/reactive-tree-driver/instance-bridge';
import './tests/compute-utils/reactive-tree-driver/links-additional-states-propagation';
import './tests/compute-utils/reactive-tree-driver/step-deps-tracking';
import './tests/compute-utils/reactive-tree-driver/obsolete-meta';
import './tests/compute-utils/reactive-tree-driver/oninit-hook';
import './tests/compute-utils/reactive-tree-driver/links-workflow';
import './tests/compute-utils/reactive-tree-driver/sequence-run';
import './tests/compute-utils/reactive-tree-driver/links-retention';
import './tests/compute-utils/reactive-tree-driver/structurecheck-hook';
import './tests/compute-utils/reactive-tree-driver/pipeline-granular-mutations';
import './tests/compute-utils/reactive-tree-driver/error-handling';
import './tests/compute-utils/reactive-tree-driver/advanced-coverage';
import './tests/compute-utils/reactive-tree-driver/buffer-keys-during-lock';
import './tests/compute-utils/fitting/nm-pure-math';
import './tests/compute-utils/fitting/worker-dg-shim';
import './tests/compute-utils/fitting/cost-functions';
import './tests/compute-utils/fitting/end-to-end';
import './tests/compute-utils/fitting/cross-executor-parity';
import './tests/arrow/roundtrip';
import './tests/arrow/titanic';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool stressTest {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext, stressTest?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, stressTest});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
