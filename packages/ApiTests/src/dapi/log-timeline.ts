import type * as _grok from 'datagrok-api/grok';
declare let grok: typeof _grok;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: log timeline', () => {
  test('timeline of the current session', async () => {
    const session = await grok.dapi.users.currentSession();
    const rows = await grok.dapi.log.getTimeline({session: session.id, limit: 20});
    expect(Array.isArray(rows), true);
    expect(rows.length <= 20, true);
    for (const r of rows) {
      expect(typeof r.time, 'string');
      expect(['client', 'server'].includes(r.source), true);
    }
  }, {owner: 'aparamonov@datagrok.ai'});
});
