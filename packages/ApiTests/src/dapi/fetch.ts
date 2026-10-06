import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../test-package';


category('Dapi: fetch', () => {
  test('post', async () => {
    const url = 'https://jsonplaceholder.typicode.com/posts';
    const data = {name: 'username', password: 'password'};

    const res = await grok.dapi.fetchProxy(url, {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(data),
    });

    if (!res.ok)
      throw new Error('Post failed');
  });  // not stressTest: POSTs to an external echo service (jsonplaceholder);
       // under load it measures that host + the internet, not the stand. The
       // 'get' test (stand-local webRoot) covers the fetchProxy path under stress.

  test('get', async () => {
    // Stand-local target so the stress signal reflects the fetchProxy round-trip, not an
    // external host (was dev.datagrok.ai). Use the stand's own API root (the same
    // /info/server the deploy health-checks) — a debug-published package's webRoot is
    // served by the publish client, which has exited by test time, so webRoot fetches
    // come back empty in the stress harness.
    const url = `${grok.dapi.root}/info/server`;

    const res = await grok.dapi.fetchProxy(url);
    const resText = await res.text();

    if (resText.length == 0)
      throw new Error('Response text is empty');
  }, {stressTest: true});
}, {owner: 'aparamonov@datagrok.ai'});
