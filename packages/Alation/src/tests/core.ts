import * as grok from 'datagrok-api/grok';

import {category, test} from '@datagrok-libraries/utils/src/test';

category('Alation', () => {
  test('App test', async () => {
    await grok.functions.call('Alation:Alation');
  });
});
