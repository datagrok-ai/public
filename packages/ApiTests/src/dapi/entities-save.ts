import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, category, expect, test} from '@datagrok-libraries/test/src/test';

const created: _DG.Entity[] = [];

category('Dapi: entities.save (polymorphic)', () => {

  test('Project: save via entities.save round-trips via projects.find', async () => {
    const project = DG.Project.create();
    project.name = `entities_save_project_${DG.Utils.randomString(6)}`;
    const saved = await grok.dapi.entities.save(project) as _DG.Project;
    created.push(saved);
    expect(!!saved.id, true);

    const reloaded = await grok.dapi.projects.find(saved.id);
    expect(reloaded.id, saved.id);
    expect(reloaded.name, project.name);
  }, {owner: 'aparamonov@datagrok.ai'});

  test('DataConnection: save via entities.save round-trips via connections.find', async () => {
    const dcParams = {
      dataSource: 'PostgresDart', server: 'localhost:5432', db: 'datagrok_dev',
      login: 'datagrok_dev', password: '123',
    };
    const dc = DG.DataConnection.create(`entities_save_conn_${DG.Utils.randomString(6)}`, dcParams);
    const saved = await grok.dapi.entities.save(dc) as _DG.DataConnection;
    created.push(saved);
    expect(!!saved.id, true);

    const reloaded = await grok.dapi.connections.find(saved.id);
    expect(reloaded.id, saved.id);
    expect(reloaded.friendlyName, dc.friendlyName);
  }, {owner: 'aparamonov@datagrok.ai'});

  // The bare EntityRecord path is exercised end-to-end by the existing
  // 'Dapi: entities' tests in entities.ts (which call entities.save with a
  // bare row through the `setProperties` flow). We don't repeat it here.

  after(async () => {
    for (const e of created) {
      try {
        if (e instanceof DG.Project) await grok.dapi.projects.delete(e);
        else if (e instanceof DG.DataConnection) await grok.dapi.connections.delete(e);
        else await grok.dapi.entities.delete(e);
      } catch (_) {}
    }
    created.length = 0;
  });

}, {owner: 'aparamonov@datagrok.ai', node: true});
