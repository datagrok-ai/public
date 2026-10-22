import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, waitUntilAsync} from '@datagrok-libraries/test/src/test';

category('Dapi: log timeline', () => {
  test('timeline of the current session', async () => {
    const session = (await grok.dapi.users.currentSession()).id;
    const from = new Date(Date.now() - 60000);
    let rows: _DG.TimelineRow[] = [];
    // Requests are written every few seconds; this call and the one above are requests of the session.
    await waitUntilAsync(async () => {
      rows = await grok.dapi.log.getTimeline({session, from});
      return rows.some((r) => r.kind === 'request');
    }, 'no request of the current session in the last minute', 30000, 1000);
    expect(rows.every((r, i) => i === 0 || Date.parse(rows[i - 1].time) <= Date.parse(r.time)), true, 'rows are not in time order');
  }, {owner: 'aparamonov@datagrok.ai'});
});
