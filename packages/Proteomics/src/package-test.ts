import * as DG from 'datagrok-api/dg';
import {runTests, TestContext, tests, initAutoTests as initTests} from '@datagrok-libraries/test/src/test';

import './tests/parsers';
import './tests/analysis';
import './tests/generic-parser';
import './tests/qc-dashboard';
import './tests/enrichment';
import './tests/enrichment-visualization';
import './tests/spectronaut-parser';
import './tests/spectronaut-candidates-parser';
import './tests/spectronaut-candidates-e2e';
import './tests/fragpipe-parser';
import './tests/fragpipe-e2e';
import './tests/subcellular-location';
import './tests/volcano';
import './tests/gene-label-resolver';
import './tests/smart-pathway-filter';
import './tests/uniprot-panel';
import './tests/group-mean-correlation';
import './tests/publish-spike';
import './tests/publish-roundtrip';
import './tests/spc-formula-lines-spike';
import './tests/spc';

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

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}
