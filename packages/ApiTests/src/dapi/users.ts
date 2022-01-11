import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: users', () => {

  test('Dapi: users - current', async () => {
    await grok.dapi.users.current();
  });

  test('Dapi: users - current session', async () => {
    await grok.dapi.users.currentSession();
  });

});
