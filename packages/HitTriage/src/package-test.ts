import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {runTests, tests, TestContext, initAutoTests as initTests } from '@datagrok-libraries/test/src/test';
import {_package} from './package';
import './tests/campaigns';
import './tests/compat';
import './tests/templates';
import './tests/template-editor';
import './tests/save-lock';
import './tests/campaign-session';
import './tests/compute';
import './tests/merge';
import './tests/campaign-layout';
import './tests/triage-start';
import './tests/vid';
import './tests/storage';
import './tests/catalog';
import './tests/ai';
import './tests/u2-guards';

export {_package, tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  // verbose: true - for tests returning dataframe
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
