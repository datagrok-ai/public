import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import './dataframe/dataframe';
import './dataframe/detector';
import './dataframe/calculated-columns';
import './dataframe/events';
import './dataframe/datetime-columns-join';
import './dataframe/dataframe-join';
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
import './utils/color';
import './package/upload';
import './viewers/viewers';
import './viewers/filters';
import './grid/grid';
import './grid/color-coding';
// import './connections/queries-test';
// import './connections/perf-tests';
import './connections/get-all-top100';
import './scripts/scripts-params';
import './gui/dialogs';
import './gui/files';
import './gui/grid';
import './gui/data-science';
import './gui/project-upload';
// import './gui/viewers/bar-chart';
// import './gui/viewers/box-plot';
// import './gui/viewers/density-plot';
// import './gui/viewers/form';
// import './gui/viewers/histogram';
// import './gui/viewers/leaflet';
// import './gui/viewers/line-chart';
// import './gui/viewers/matrix-plot';
// import './gui/viewers/network-diagram';
// import './gui/viewers/pc-plot';
// import './gui/viewers/pie-chart';
// import './gui/viewers/scatter-plot';
// import './gui/viewers/scatter-plot-3d';
// import './gui/viewers/word-cloud';

import {runTests, tests, TestContext} from '@datagrok-libraries/utils/src/test';
import {Column, DataFrame, DataQuery, FuncCall} from 'datagrok-api/dg';
import {filter} from 'rxjs/operators';

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


//name: testConnections
//output: dataframe result
//top-menu: Tools | Dev | Test Connections
export async function testConnections(): Promise<DG.DataFrame> {
  const connections: string[] = ['PostgreSQLApiTests', 'SnowflakeApiTests', 'MSSQLApiTests', 'OracleApiTests'];
  const tables: string[] = ['Long', 'Normal', 'Wide', 'Tiny'];
  const fetchSizes: string[] = ['big', 'dynamic', 'low'];

  // const queriesFriendlyNames: string[] = ['PostgresNormal', 'PostgresLong', 'PostgresWide'];
  const l = connections.length * tables.length * fetchSizes.length;

  const df = DataFrame.fromColumns([Column.string('type', l), Column.string('fetch', l), 
    Column.string('db', l), Column.int('TTFR', l), Column.int('TTC', l)]);

  let startTime: number;
  let ttfr: number;

  let callCheck: (value: FuncCall) => boolean;
  let ttfrSet = false;
  grok.functions.onParamsUpdated.pipe(filter((c) => callCheck(c) && !ttfrSet)).subscribe(() => {
    ttfr = Date.now() - startTime;
    df.columns.byName('TTFR').set(row, ttfr);
    ttfrSet = true; 
  });

  let row = 0;
  for (const con of connections) {
    for (const table of tables) {
      for (const fetchSize of fetchSizes) {
        if (table == 'Long' && fetchSize == 'low')
          continue;  

        const connection = await grok.dapi.connections.filter(`name = "${con}"`).first();
        ttfrSet = false;

        df.columns.byName('type').set(row, table); 
        df.columns.byName('fetch').set(row, fetchSize); 
        df.columns.byName('db').set(row, con); 

        callCheck = (c: FuncCall) => c.aux.get('fetchSize') == fetchSize && 
          (c.func as DataQuery).connection.name == con;

        const preTable = con.startsWith('Snowflake') ? 'TEST.' : '';

        let sql;

        if (table == 'Tiny') {
          if (con.startsWith('Oracle'))
            sql = `select * from Test_Long WHERE ROWNUM = 1`;
          else if (con.startsWith('MS'))
            sql = `select TOP 1 * from Test_Long`;
          else 
            sql = 'select 1';
        } else 
          sql = `select * from ${preTable}Test_${table}`;
        
        const query = `--fetchSize: ${fetchSize}\n${sql}\n--end`;

        startTime = Date.now();

        console.log('executing' + query);
        const q = connection.query('adhoc', sql);
        const call = q.prepare();
        call.setAuxValue('fetchSize', fetchSize);
        await call.call();
        // await grok.data.db.query(connection.id, query);
        console.log('executed');

        df.columns.byName('TTC').set(row, Date.now() - startTime);

        row++;
      };
    };
  };
  df;
  grok.shell.addTableView(df);
  return df;
}
