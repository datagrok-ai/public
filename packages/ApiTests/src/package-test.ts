import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import './dataframe/dataframe';
import './dataframe/detector';
import './dataframe/calculated-columns';
import './dataframe/events';
import './dataframe/datetime-columns-join';
import './functions/functions';
import './shell/shell';
import './shell/windows';
import './views/docking';
import './views/events';
import './views/layouts';
import './views/files-view';
import './dapi/files';
import './dapi/functions';
import './dapi/fetch';
import './dapi/groups';
import './ui/inputs';
import './ui/forms';
import './dapi/dapi';
import './dapi/connection';
import './dapi/entities';
import './dapi/layouts';
import './dapi/packages';
import './dapi/projects';
import './dapi/tables';
import './dapi/user-data-storage';
import './dapi/users';
import './dapi/benchmarks';
import './shell/ml';
import './shell/settings';
import './ui/divs';
import './ui/buttons';
import './widgets/files-widget';
import './widgets/legend';
import './widgets/tree-view';
import './ui/icons';
import './ui/tables';
import './ui/range-slider';
import './ui/accordion';
import './ui/tab-control';
import './ui/list';
import './ui/image';
import './ui/users';
import './ui/groups';
import './ui/tags';
import './ui/sharing';
import './package/upload';
import './viewers/viewers';
import './grid/grid';
import './grid/color-coding';
import './connections/queries-test';
import './connections/get-all-top100';
import './scripts/scripts-params';
import './gui/dialogs';
import './gui/files';
import './gui/grid';
import './gui/project-upload';

import {runTests, tests, TestContext} from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();
export {tests};


//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
//top-menu: Tools | Dev | JS API Tests
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}

//name: testPackages
//output: dataframe result
//top-menu: Tools | Dev | Test Packages
export async function testPackages(): Promise<DG.DataFrame> {
  const funcs = DG.Func.find({name: 'test'});
  const dfs:DG.DataFrame[] = [];
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
