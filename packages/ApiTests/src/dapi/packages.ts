import {category, expect, expectObject, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: packages', () => {
  test('list', async () => {
    const list = (await grok.dapi.packages.list());
    expect(list.some((pack) => pack.name === 'Demo'), true);
  });
  
  test('find', async () => {
    const apiTestsPack = 
      DG.toJs((await grok.dapi.packages.list()).filter((pack) => pack.name === 'ApiTests')[0]) as DG.Package;
    expectObject((await grok.dapi.packages.find(apiTestsPack.id)), apiTestsPack);
  });
});
