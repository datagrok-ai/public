import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi', () => {
  // One test per data source: the old combined test summed 14 sequential
  // round trips, so at high concurrency it drew a single misleading 3s+ band
  // in the stress report instead of showing which endpoint is actually slow.
  const sources: {[name: string]: () => Promise<any>} = {
    'queries': () => grok.dapi.queries.first(),
    'connections': () => grok.dapi.connections.first(),
    'credentials': () => grok.dapi.credentials.first(),
    'jobs': () => grok.dapi.jobs.first(),
    'notebooks': () => grok.dapi.notebooks.first(),
    'models': () => grok.dapi.models.first(),
    'packages': () => grok.dapi.packages.first(),
    'layouts': () => grok.dapi.layouts.first(),
    'tables': () => grok.dapi.tables.first(),
    'users': () => grok.dapi.users.first(),
    'groups': () => grok.dapi.groups.first(),
    'scripts': () => grok.dapi.scripts.first(),
    'projects': () => grok.dapi.projects.first(),
    'environments': () => grok.dapi.environments.first(),
  };
  for (const name of Object.keys(sources)) {
    test(`all data sources: ${name}`, async () => {
      await sources[name]();
    }, {owner: 'aparamonov@datagrok.ai', stressTest: true});
  }

  test('logging', async () => {
    const logger = new DG.Logger((m) => (m.params as {[key: string]: any})['jsApiTest2'] = 'jsApiTest3');
    const jsApiTestType = 'jsApiTestType';
    logger.log('jsApiTest0', {jsApiTest1: 'jsApiTest2'}, jsApiTestType);
    for (let i = 0; i < 100; i++) {
        await DG.delay(100);
        if ((await grok.dapi.logTypes.list({filter: jsApiTestType}))[0]?.name == jsApiTestType)
            return;
    }
    throw new Error('Log type not found');
    //TODO: find log
    // console.log(await grok.dapi.log.list({filter: 'jsApiTest1 = "jsApiTest2"'}));
  }, {owner: 'aparamonov@datagrok.ai', skipReason: typeof process !== 'undefined' ? 'NodeJS environment' : undefined});
}, {node: true});
