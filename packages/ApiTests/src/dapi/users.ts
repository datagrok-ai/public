import type * as _grok from 'datagrok-api/grok';
declare let grok: typeof _grok;

import {category, test} from '@datagrok-libraries/test/src/test';

category('Dapi: users', () => {
  test('current', async () => {
    await grok.dapi.users.current();
  }, {stressTest: true});

  test('current session', async () => {
    await grok.dapi.users.currentSession();
  }, {stressTest: true});
}, {owner: 'aparamonov@datagrok.ai'});
