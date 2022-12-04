import * as DG from 'datagrok-api/dg';

import {runTests, TestContext, tests} from '@datagrok-libraries/utils/src/test';

import './tests/Palettes-test';
import './tests/detectors-tests';
import './tests/detectors-benchmark-tests';
import './tests/msa-tests';
import './tests/sequence-space-test';
import './tests/activity-cliffs-tests';
import './tests/splitters-test';
import './tests/renderers-test';
import './tests/converters-test';
import './tests/fasta-handler-test';
import './tests/fasta-export-tests';
import './tests/bio-tests';
import './tests/WebLogo-positions-test';
import './tests/checkInputColumn-tests';
import './tests/similarity-diversity-tests';
import './tests/substructure-filters-tests';

export const _package = new DG.Package();
export {tests};

/** For the 'test' function argument names are fixed as 'category' and 'test' because of way it is called. */
//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}
