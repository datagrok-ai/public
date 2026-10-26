import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {_package} from './package';
import './tests/priors-tests';
import './tests/plan-tests';
import './tests/analyze-tests';
import './tests/band-chart-tests';

// the platform initialises the _package this bundle exports: the same object priors.ts and crux.ts read
export {_package};
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool excludeNodeTests {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext,
  excludeNodeTests?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, excludeNodeTests});
  return DG.DataFrame.fromObjects(data)!;
}

/** Headless entry for the `grok test` Node pass — runs only tests marked {node: true}. */
export async function testNode(pkg: DG.Package,
  options: {category?: string, test?: string, stressTest?: boolean, verbose?: boolean}): Promise<any[]> {
  return await runTests({category: options.category, test: options.test, stressTest: options.stressTest,
    verbose: options.verbose, nodeOnly: true, nodeOptions: {package: pkg}});
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
