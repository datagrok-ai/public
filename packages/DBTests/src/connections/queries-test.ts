import {before, category, expect, test, timeout} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {TableInfo} from "datagrok-api/dg";

category('Connections', () => {
  test('External Provider Chembl Perf', async () => {
    if (!DG.Test.isInBenchmark) return;
    const query = await grok.dapi.queries.filter(`friendlyName = "ChemblPerfGenerated"`).include('params').first();
    const lim = 500000;
    const call = query.prepare({'num': lim});
    await call.call();
    const t = call.getOutputParamValue() as DG.DataFrame;
    if (t.rowCount != lim)
      // eslint-disable-next-line no-throw-literal
      throw 'Rows number in' + query.name + 'table is not as expected';
  }, {benchmark: true});

  test('External Provider: Columns in empty result', async () => {
    const query = await grok.dapi.queries.filter(`friendlyName = "TestForColumnsOnEmptyResult"`).include('params').first();
    const call = query.prepare();
    await call.call();
    const t = call.getOutputParamValue() as DG.DataFrame;
    expect(t.columns.length, 10);
    expect(t.columns.contains('first_name'), true);
    expect(t.columns.byName('first_name').length, 0);
  }, {stressTest: true});

  test('External Provider: grok.data.query no params', async () => {
    const result: DG.DataFrame = await grok.data.query('DbTests:PostgresqlPatternsAll', null, true);
    expect(result?.rowCount ?? 0, 30);
  });

  test('External Provider: grok.data.query with params', async () => {
    const result: DG.DataFrame = await grok.data.query('DbTests:PostgresqlPatternsAllParams', {'first_name': 'starts with p', 'id': '>1', 'bool': false,
      'email': 'contains com', 'some_number': '>20', 'country': 'in (Indonesia)', 'date': 'before 1/1/2022'});
    expect(result?.rowCount ?? 0, 1);
  });
});

category('Docker connection', () => {
  let testConnection: DG.DataConnection | null;

  before(async () => {
    testConnection = await grok.functions.eval('DbTests:PostgresDocker');
  });

  test('Connection test', async () => {
    await testConnection!.test();
  }, {timeout: 120000 /* on demand start */});

  test('Connection getSchemas', async () => {
    const schemas: string[] = await grok.dapi.connections.getSchemas(testConnection!);
    expect(schemas.includes('public'));
    expect(schemas.includes('information_schema'));
  });

  test('Connection getSchema', async () => {
    const schema: TableInfo[] = await grok.dapi.connections.getSchema(testConnection!, 'public');
    expect(schema.some((ti) => ti.name === 'City'));
    expect(schema.some((ti) => ti.name === 'Country'));
    expect(schema.some((ti) => ti.name === 'Countrylanguage'));
    expect(schema.find((ti) => ti.name === 'City')!.columns.some((ci) => ci.name === 'id'));
  });

  test('Query', async () => {
    const query: DG.DataQuery = await testConnection!.query('test', 'select * from city');
    const call: DG.FuncCall = await query.prepare().call();
    const result: DG.DataFrame = call.getOutputParamValue();
    expect(result.rowCount, 4079);
  });
});
