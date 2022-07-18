import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import './dataframe/dataframe';
import './dataframe/calculated-columns';
import './dataframe/events';
import './functions/functions';
import './shell/shell';
import './shell/windows';
import './viewer/viewer';
import './views/docking';
import './views/events';
import './views/layouts';
import './dapi/files';
import './dapi/fetch';
import './dapi/admin';
import './dapi/groups';
import './ui/inputs';
import './ui/forms';
import './dapi/dapi';
import './dapi/connection';
import './dapi/entities';
import './dapi/layouts';
import './dapi/projects';
import './dapi/tables';
import './dapi/user-data-storage';
import './dapi/users';
import './shell/ml';
import './ui/divs';
import './ui/buttons';
import './widgets/legend';
import './ui/icons';
import './ui/tables';
import './ui/range-slider';
import './ui/accordion';
import './ui/tab-control';
import './ui/list';
import './ui/image';
import './gui/anonimize';
import './gui/cluster';
import './ui/viewers-adding';
import './grid/grid';
import './gui/random-data';
import './gui/missing-values-imputation';
import './connections/demo_queries';
import './gui/aggregate-rows';
import './gui/multivariate-analysis';
import './gui/pca';
import './gui/spe';
import './gui/split-columns';
import './gui/add-new-column';
import './gui/data-search';
import './gui/delete-rows-columns';
import './gui/filters';
import './gui/unpivot';
import './gui/scatter-plot'
import './gui/bar-chart'
import './gui/pie-chart'
import './gui/histogram'
import './gui/matrix-plot'
import './gui/pc-plot'
import './gui/density-plot'
import './gui/scatter-plot-3d.ts'
import './gui/word-cloud'
import './gui/network-diagram'
import './gui/box-plot'

import {runTests, tests} from '@datagrok-libraries/utils/src/test';
export const _package = new DG.Package();
export {tests};


//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
//top-menu: Tools | Dev | JS API Tests
export async function test(category: string, test: string): Promise<DG.DataFrame> {
  const data = await runTests({category, test});
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
