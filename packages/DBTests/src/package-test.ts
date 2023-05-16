import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {runTests, tests, TestContext} from '@datagrok-libraries/utils/src/test';
import {Column, DataFrame, DataQuery, FuncCall} from 'datagrok-api/dg';
import './connections/queries-test';
import './connections/get-all-top100';

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

//name: testConnections
//output: dataframe result
//top-menu: Tools | Dev | Test Connections
export async function testConnections(): Promise<DG.DataFrame> {
  const connections: string[] = ['PostgreSQLDBTests', 'SnowflakeDBTests', 'MSSQLDBTests', 'OracleDBTests'];
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
  // @ts-ignore
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
            // @ts-ignore
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
        // @ts-ignore
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
