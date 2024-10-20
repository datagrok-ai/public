import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

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
});
