// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

// import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';


// category('Demo: Coffee_company', () => {
//     test('project.coffee_company', async () => {
//       const p = await grok.dapi.projects.open('coffee_company')
//       expect(p instanceof DG.Project,true)
    
//     });
// });

import {after, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package-test';

category('Connections', () => {
  test('queriesTest', async () => {
    const queries = await grok.dapi.queries.filter(`options.testExpectedRows != null and package.shortName = "${_package.name}"`).include('params').list();
   
    for (const query of queries) {
      const call = query.prepare();

      console.log("NAME: " + query.name);      

      for (const property of query.inputs)
        property.set(call, property.defaultValue == null ? property.defaultValue : await grok.functions.eval(`${property.defaultValue}`));

      await call.call();

      const t = call.getOutputParamValue() as DG.DataFrame;
      if (t == null)
        throw 'Result of ' + query.name + 'is not a DataFrame';

      if (t.rowCount != query.options.testExpectedRows)
        // eslint-disable-next-line no-throw-literal
        throw 'Rows number in' + query.name + 'table is not as expected';
    }
  });
});
