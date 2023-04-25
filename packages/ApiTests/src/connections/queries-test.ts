import {category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

category('Connections', () => {
  
  test('queriesTest', async () => {
    const queries = await grok.dapi.queries
      .filter(`options.testExpectedRows != null and package.shortName = "${_package.name}"`).include('params').list();
    if (queries.length == 0)
      throw new Error('test queries were not found!');
    for (const query of queries) {
      const call = query.prepare();
      for (const property of query.inputs) {
        property.set(call, property.defaultValue == null ? property.defaultValue :
          await grok.functions.eval(`${property.defaultValue}`));
      }
      await call.call();
      const t = call.getOutputParamValue() as DG.DataFrame;
      if (t == null)
        throw new Error('Result of ' + query.name + 'is not a DataFrame');
      if (t.rowCount.toString() != query.options.testExpectedRows)
        // eslint-disable-next-line no-throw-literal
        throw 'Rows number in' + query.name + 'table is not as expected';
    }
  });

  test('perfGen', async () => {
    const query = await grok.dapi.queries.include('params,connection').filter(`friendlyName="Perf"`).first();
    const call = query.prepare();
    await call.call();
    const t = call.getOutputParamValue() as DG.DataFrame;
    console.log(t);
    console.log(query);
  });

  test('ScalarQueryTest', async () => {
    const query = await grok.dapi.queries.filter(`friendlyName = "Postgre Scalar Output"`).include('params').first();
    const call = query.prepare();
    await call.call();  
    const t = call.getOutputParamValue() as number;
    console.log(t);
    if (t != 830)
      // eslint-disable-next-line no-throw-literal
      throw 'Rows number in' + query.name + 'table is not as expected';
  });
  
  test('External Provider Chembl Perf', async () => {
    const query = await grok.dapi.queries.filter(`friendlyName = "ChemblPerfGenerated"`).include('params').first();
    const lim = 5000000;
    const call = query.prepare({'num': lim});
    await call.call();  
    const t = call.getOutputParamValue() as DG.DataFrame;
    if (t.rowCount != lim)
      // eslint-disable-next-line no-throw-literal
      throw 'Rows number in' + query.name + 'table is not as expected';
  });

  test('External Provider First part', async () => {
    const query = await grok.dapi.queries.filter(`friendlyName = "Compounds"`).include('params').first();
    const call = query.prepare();
    call.call();
    setTimeout(() => {
      const t = call.getOutputParamValue() as DG.DataFrame;
      if (t == null)
      // eslint-disable-next-line no-throw-literal
        throw 'First rows await time exceeded';
    }, 2000);
  });
});
