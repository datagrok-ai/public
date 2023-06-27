import {category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

category('Dapi: layouts', () => {
  test('get applicable', async () => {
    await grok.dapi.layouts.getApplicable(grok.data.demo.demog());
  });
});
