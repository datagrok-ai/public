import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

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
  }, {skipReason: 'GROK-11670',  owner: 'aparamonov@datagrok.ai'});

  test('logging', async () => {
    const logger = new DG.Logger((m) => (m.params as {[key: string]: any})['jsApiTest2'] = 'jsApiTest3');
    const jsApiTestType = 'jsApiTestType';
    logger.log('jsApiTest0', {jsApiTest1: 'jsApiTest2'}, jsApiTestType);
    expect((await grok.dapi.logTypes.list({filter: jsApiTestType}))[0]?.name, jsApiTestType);
    //TODO: find log
    // console.log(await grok.dapi.log.list({filter: 'jsApiTest1 = "jsApiTest2"'}));
  }, {skipReason: 'GROK-11670',  owner: 'aparamonov@datagrok.ai'});
});
