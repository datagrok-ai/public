import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: layouts', () => {
  test('Dapi: layouts - get applicable', async () => {
    const layouts = await grok.dapi.layouts.getApplicable(grok.data.demo.demog());
  });
});
