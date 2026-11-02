import * as DG from 'datagrok-api/dg';
import {TestContext, runTests, tests, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';

import './test/viewers-hook';
import './test/app-launch';
import './test/rfv';
import './test/navigation';
import './test/float-display';
import './test/custom-export';
import './test/url-inputs';
import './test/project-export-data';
import './test/run-comparison-extraction';
import './test/run-comparison-views';
import './test/run-comparison-matching';
import './test/run-comparison-selection';
import './test/run-comparison-builders';
import './test/run-comparison-target-selection';
import './test/run-comparison-chart-spec';
import './test/run-comparison-component';
import './test/viewer-frame-swap-contract';
import './test/graphics-output';
import './test/dock-layout';
import './test/validation-actions';
import './test/input-toggles';
import './test/tree-node';
import './test/export-summary';
import './test/headless-run';
import './test/run-storage';
import './test/space-index';
import './test/help-service';
import './test/space-walk';
import './test/custom-function-view';
import './test/feedback';
import './test/tree-utils';
import './test/inspector-annotation-links';

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
  await DG.Func.byName('WebComponents:init').prepare().call();
  await initTests(_package, _package.getModule('package-test.js'));
}
