import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi', () => {
  test('all data sources', async () => {
    await grok.dapi.queries.first();
    await grok.dapi.connections.first();
    await grok.dapi.credentials.first();
    await grok.dapi.jobs.first();
    await grok.dapi.notebooks.first();
    await grok.dapi.models.first();
    await grok.dapi.packages.first();
    await grok.dapi.layouts.first();
    await grok.dapi.tables.first();
    await grok.dapi.users.first();
    await grok.dapi.groups.first();
    await grok.dapi.scripts.first();
    await grok.dapi.projects.first();
    await grok.dapi.environments.first();
  }, {owner: 'aparamonov@datagrok.ai', stressTest: true});

  test('logging', async () => {
    const logger = new DG.Logger((m) => (m.params as {[key: string]: any})['jsApiTest2'] = 'jsApiTest3');
    const jsApiTestType = 'jsApiTestType';
    logger.log('jsApiTest0', {jsApiTest1: 'jsApiTest2'}, jsApiTestType);
    await DG.delay(1000);
    expect((await grok.dapi.logTypes.list({filter: jsApiTestType}))[0]?.name, jsApiTestType);
    //TODO: find log
    // console.log(await grok.dapi.log.list({filter: 'jsApiTest1 = "jsApiTest2"'}));
  }, {owner: 'aparamonov@datagrok.ai', skipReason: typeof process !== 'undefined' ? 'NodeJS environment' : undefined});
});
