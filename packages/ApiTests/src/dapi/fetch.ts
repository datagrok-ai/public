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
  }, {stressTest: true});

  test('get', async () => {
    // Stand-local asset (the package's own webRoot) so the stress signal reflects the
    // proxy round-trip, not the load/latency of an external host (was dev.datagrok.ai).
    const url = `${_package.webRoot}files/cars.csv`;

    const res = await grok.dapi.fetchProxy(url);
    const resText = await res.text();

    if (resText.length == 0)
      throw new Error('Response text is empty');

    DG.DataFrame.fromCsv(resText);
  }, {stressTest: true});
}, {owner: 'aparamonov@datagrok.ai'});
