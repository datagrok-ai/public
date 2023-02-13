import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';

const GDC = grok.dapi.connections;

category('Dapi: connection', () => {
  const dcParams = {
    dataSource: 'PostgresDart', server: 'localhost:5432', db: 'datagrok_dev', login: 'datagrok_dev', password: '123'};

  test('Create, save, delete, share', async () => {
    let dc = DG.DataConnection.create('Local DG Test', dcParams);
    dc = await GDC.save(dc);
    expect((dc.parameters as any)['schema'], null);
    expect((dc.parameters as any)['db'], dcParams.db);
    expect(dc.friendlyName, 'Local DG Test');
    expect((await GDC.find(dc.id)).id, dc.id);
    await GDC.delete(dc);
    expect(await GDC.find(dc.id), undefined);
  });
});
