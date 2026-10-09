import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: log timeline', () => {
  test('timeline of the current session', async () => {
    const session = (await grok.dapi.users.currentSession()).id;
    const from = new Date(Date.now() - 60000);
    let rows = await grok.dapi.log.getTimeline({session, from});
    // Requests are written every few seconds; this call and the one above are requests of the session.
    for (let i = 0; i < 30 && !rows.some((r) => r.kind === 'request'); i++) {
      await DG.delay(1000);
      rows = await grok.dapi.log.getTimeline({session, from});
    }
    expect(rows.some((r) => r.kind === 'request'), true, 'no request of the current session in the last minute');
    expect(rows.every((r, i) => i === 0 || rows[i - 1].time <= r.time), true, 'rows are not in time order');
  }, {owner: 'aparamonov@datagrok.ai'});
});
