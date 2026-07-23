import * as DG from 'datagrok-api/dg';

// import './scripts/scripts';
import './shell/ml';
import './scripts/scripts_tests';
import './celery/celery_tests';
import './celery/node_celery_tests';
import './celery/node_celery_benchmark';
import './docker/docker';
import './files/files';
// import './gui/dialogs'; To fix!

import {runTests, tests, TestContext, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool stressTest {optional: true}
//input: bool excludeNodeTests {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext, stressTest?: boolean,
    excludeNodeTests?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext, stressTest, excludeNodeTests});
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
