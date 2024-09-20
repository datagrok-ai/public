import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {runTests, TestContext, tests, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';

import './tests/_first-tests';
import './tests/Palettes-test';
import './tests/detectors-tests';
import './tests/detectors-weak-and-likely-tests';
import './tests/detectors-benchmark-tests';
import './tests/msa-tests';
import './tests/splitters-test';
import './tests/monomer-libraries-tests';
import './tests/renderers-test';
import './tests/renderers-monomer-placer-tests';
import './tests/converters-test';
import './tests/fasta-handler-test';
import './tests/fasta-export-tests';
import './tests/bio-tests';
import './tests/WebLogo-positions-test';
import './tests/WebLogo-project-tests';
import './tests/WebLogo-layout-tests';
import './tests/checkInputColumn-tests';
import './tests/similarity-diversity-tests';
import './tests/substructure-filters-tests';
import './tests/pepsea-tests';
import './tests/viewers';
import './tests/seq-handler-tests';
import './tests/seq-handler-splitted-tests';
import './tests/seq-handler-get-region-tests';
import './tests/seq-handler-get-helm-tests';
import './tests/to-atomic-level-tests';
import './tests/to-atomic-level-ui-tests';
import './tests/mm-distance-tests';
import './tests/activity-cliffs-tests';
import './tests/sequence-space-test';
import './tests/scoring';


export const _package = new DG.Package();
export {tests};

/** For the 'test' function argument names are fixed as 'category' and 'test' because of way it is called. */
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
