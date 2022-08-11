import { after, before, category, delay, expect, test } from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import {DataFrame} from 'datagrok-api/dg';

category('Connections: Demo Queries', () => {

  test('connections.demoQueries', async () => {
    let queries = await grok.dapi.queries.filter('#unit-test').include('params').list();

    for (let query of queries) {
      let call = query.prepare();

      for (let property of query.inputs) {
        property.set(call, await grok.functions.eval(property.defaultValue));
      }

      await call.call()

      let t = call.getOutputParamValue() as DataFrame;

      if (t.rowCount != query.options.testExpected) {
        throw ' Rows number in' + query.name + 'table is not as expected';
      }
    }
  });
});