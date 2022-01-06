import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: admin', () => {

  test('Dapi: admin - getServiceInfos', async () => {
    if((await grok.dapi.admin.getServiceInfos()).length == 0)
      throw 'No services';
  });

});
