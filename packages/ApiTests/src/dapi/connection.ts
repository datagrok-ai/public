import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expectObject, expect, test} from '@datagrok-libraries/utils/src/test';

const GDC = grok.dapi.connections;

category('Dapi: connection', () => {
  const dcParams = {
    dataSource: 'PostgresDart', server: 'localhost:5432', db: 'datagrok_dev', login: 'datagrok_dev', password: '123'};

  test('Create, save, delete, share', async () => {
    let dc = DG.DataConnection.create('Local DG Test', dcParams);
    dc = await GDC.save(dc);
    expectObject(dc.parameters, {server: 'localhost:5432', db: 'datagrok_dev'});
    expect(dc.friendlyName, 'Local D G Test');
    expect((await GDC.find(dc.id)).id, dc.id);
    await GDC.delete(dc);
    expect(await GDC.find(dc.id), undefined);
  });
});
