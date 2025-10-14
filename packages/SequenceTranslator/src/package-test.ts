import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/utils/src/test';

import './tests/formats-to-helm';
import './tests/helm-to-nucleotides';
import './tests/formats-support';
import './tests/files-tests';
import './tests/polytool-detectors-custom-notation-test';
import './tests/polytool-convert-tests';
import './tests/polytool-unrule-tests';
import './tests/polytool-enumerate-tests';
import './tests/polytool-enumerate-breadth-tests';
import './tests/polytool-chain-parse-notation-tests';
import './tests/polytool-chain-from-notation-tests';
import './tests/toAtomicLevel-tests';

import {OligoToolkitTestPackage} from './tests/utils';

export const _package = new OligoToolkitTestPackage();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, verbose: true});
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package as any, _package.getModule('package-test.js'));
}
