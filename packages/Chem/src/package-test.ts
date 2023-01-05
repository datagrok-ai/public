import * as DG from 'datagrok-api/dg';
import {runTests, TestContext, tests} from '@datagrok-libraries/utils/src/test';

import './tests/menu-tests-chem-space';
import './tests/menu-tests-cliffs';
import './tests/menu-tests-similarity-diversity';
import './tests/menu-tests-script-based';
import './tests/menu-tests-rgroups';

import './tests/col-panel-tests';
import './tests/cell-panel-tests';

import './tests/substructure-search-tests';
import './tests/rendering-tests';
import './tests/sketcher-tests';

// import './tests/ui-tests';
import './tests/detector-tests';
import './tests/api-based-tests';
import './tests/notation-converter-tests';
import './tests/screening-tools';

import './tests/save-as-sdf-tests';
import './tests/substructure-filter-tests';

import './tests/viewers';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

