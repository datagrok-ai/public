import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {runTests, TestContext, tests} from '@datagrok-libraries/utils/src/test';

import './tests/grid-with-tree-tests';
import './tests/hierarchical-clustering-tests';
import './tests/tree-cut-tests';
import './tests/tree-for-grid-tests';
import './tests/tree-helpers-tests';

export const _package = new DG.Package();
export {tests};

/*
Entry point 'test' is required in webpack.config.js

entry: {
  test: {
    filename: 'package-test.js',
    library: {type: 'var', name: `${packageName}_test`},
    import: './src/package-test.ts',
  },
  package: './src/package.ts',
}
*/


//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}
