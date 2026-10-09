import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: log timeline', () => {
  test('timeline of a request', async () => {
    const script = await grok.dapi.scripts.save(DG.Script.create(
      `//name: ApiTestsLogTimeline${Date.now()}\n//language: javascript\n//output: int x\nx = 1;`));
    try {
      let event: _DG.LogEvent | undefined;
      for (let i = 0; i < 30 && !event; i++) {
        event = (await grok.dapi.log.where({entityId: script.id}).list({pageSize: 10})).find((e) => e.requestId != null);
        if (!event)
          await DG.delay(1000);
      }
      expect(event != null, true, 'no event of the saved script carries a request id');
      const rows = await grok.dapi.log.getTimeline({request: event!.requestId!});
      expect(rows.some((r) => r.requestId === event!.requestId), true);
    }
    finally {
      await grok.dapi.scripts.delete(script);
    }
  }, {owner: 'aparamonov@datagrok.ai'});
});
