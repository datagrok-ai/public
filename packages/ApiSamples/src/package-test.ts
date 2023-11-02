import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {runTests, tests, TestContext, category, test as _test, delay} from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();
export {tests};

import './tests/test';

const skip = [
  'function-events', 'demo', 'ui-events', 'last-error', 'chem-benchmark', 'menu-customization',
  'README', '10k-columns-updates',
  // To fix
  'scatter-plot-3d', 'network-diagram', // Break
  'files', 'open-table-by-id', 'all-chembl-structures', 'charts-in-cells',
  'property-grid', 'tree-view-adv', 'attached-properties', 'all-input-types',
  'add-single-filter', 'custom-filters', 'filter-group', 'dynamic-loading'
];

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//tags: init
export async function initTests() {
  const scripts = await grok.dapi.scripts.filter('package.shortName = "ApiSamples"').list();
  for (const script of scripts) {
    category(('Scripts:' + script.options.path as string).replaceAll('/', ':'), () => {
      _test(script.friendlyName, async () => {
        await script.apply();
        await delay(300);
        // if (grok.shell.lastError) {
        //   const err = grok.shell.lastError;
        //   grok.shell.lastError = '';
        //   throw new Error(err);
        // }
      }, skip.includes(script.friendlyName) ? {skipReason: 'skip'} : undefined);
    });
  }
}
