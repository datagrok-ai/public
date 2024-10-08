import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import './dataframe/dataframe';
import './dataframe/detector';
import './dataframe/calculated-columns';
import './dataframe/events';
import './dataframe/datetime-columns-join';
import './dataframe/dataframe-join';
import './dataframe/dataframe-link';
// import './dataframe/add-new-column';
import './functions/functions';
import './functions/conversion-functions';
import './functions/date-functions';
import './functions/logical-functions';
import './functions/math-functions';
import './functions/stats-functions';
import './functions/text-functions';
import './functions/cache';
import './shell/shell';
import './shell/ml';
import './shell/settings';
import './dapi/files';
import './dapi/files-list';
import './dapi/functions';
import './dapi/fetch';
import './dapi/groups';
import './dapi/dapi';
import './dapi/connection';
import './dapi/entities';
import './dapi/layouts';
import './dapi/packages';
import './dapi/projects';
import './dapi/tables';
import './dapi/sticky_meta';
import './dapi/users';
import './dapi/benchmarks';
import './dapi/functions-annotations';
import './widgets/files-widget';
import './widgets/legend';
import './widgets/tree-view';
import './utils/color';
// import './package/upload';
import './packages/properties';
import './packages/docker';
import './packages/user-settings-storage';
import './grid/grid';
import './grid/color-coding';
import './grid/multi-value-column';
import './stats/stats';
// import './bitset/bitset';
import './valuematcher/valuematcher';
import './property/property';
import './widgets/input-form';

import { runTests, tests, TestContext, initAutoTests as initTests } from '@datagrok-libraries/utils/src/test';

export const _package = new DG.Package();
export { tests };

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//input: bool stressTest {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext, stressTest?: boolean): Promise<DG.DataFrame> {
  const data = await runTests({ category, test, testContext, stressTest });
  return DG.DataFrame.fromObjects(data)!;
}

//name: testPlatform
//output: dataframe result
//top-menu: Tools | Dev | Test | Platform
export async function testPlatform(): Promise<DG.DataFrame> {
  const skip = ['Benchmarks', 'Packages: Docker', 'Functions: Client-side cache'];
  const data = await runTests({ exclude: skip });
  return DG.DataFrame.fromObjects(data)!;
}

//name: testPackages
//output: dataframe result
//top-menu: Tools | Dev | Test | Packages
export async function testPackages(): Promise<DG.DataFrame> {
  const funcs = DG.Func.find({ name: 'test' });
  const dfs: DG.DataFrame[] = [];
  for (const f of funcs) {
    if (f.package?.name != null) {
      grok.shell.closeAll();
      grok.shell.info(`Testing ${f.package.name}`);
      const df = await f.apply();
      if (df == null) {
        grok.shell.error(`Failed to fetch test results from ${f.package.name}`);
        continue;
      }
      const packageColumn = DG.Column.string('package', df.rowCount);
      packageColumn.init((n) => f.package.name);
      df.columns.insert(packageColumn, 0);
      dfs.push(df);
      grok.shell.closeAll();
    }
  }

  let result: DG.DataFrame | null = null;
  for (const df of dfs) {
    if (result == null)
      result = df;
    else result.append(df, true);
  }

  return result!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
} 
