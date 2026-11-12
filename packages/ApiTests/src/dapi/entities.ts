import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, waitUntilAsync} from '@datagrok-libraries/test/src/test';

//@ts-ignore
import { _package } from '../test-package';


category('Dapi: entities', () => {
  test('getProperties', async () => {
    let group: _DG.Group | undefined;
    try {
      group = await createTestGroup();
      const props = await group.getProperties();
      expect(typeof props === 'object', true);
      expect(Object.keys(props).length, 1);
    } finally {
      await safeDeleteGroup(group);
    }
  });

  test('setProperties', async () => {
    let group: _DG.Group | undefined;
    try {
      group = await createTestGroup();
      await group.setProperties({testProp1: 'prop1', testProp2: 'prop2'});
      expect(Object.keys(await group.getProperties()).length, 3);
    } finally {
      await safeDeleteGroup(group);
    }
  });

}, { owner: 'oserhiienko@datagrok.ai'});

category('Dapi: entities: smart search', () => {
  test('users', async () => {
    expect((await grok.dapi.users.filter('firstName = "admin"').list()).length, 1);
    expect((await grok.dapi.users.filter('status = "active" && firstName = "admin"').list({pageSize: 5})).length, 1);
    expect((await grok.dapi.users.filter('status != "active" && firstName = "admin"').list({pageSize: 5})).length, 0);
    expect((await grok.dapi.users.filter('id = "878c42b0-9a50-11e6-c537-6bf8e9ab02ee"').list()).length, 1);
  }, {stressTest: true});

  test('groups', async () => {
    expect((await grok.dapi.groups.filter('nistratora').list()).length, 1);
    expect((await grok.dapi.groups.filter('friendlyName = "all users"').list()).length, 1);
    expect((await grok.dapi.groups.filter('id = "1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5"').list()).length, 1);
  }, {stressTest: true});

  test('packages', async () => {
    expect((await grok.dapi.packages.filter('name="Api Tests"').list()).length > 0, true);
    expect((await grok.dapi.packages.filter(`name="Api Tests" and version = "${_package.version}"`).list()).length > 0, true);
  }, {stressTest: true});
}, {owner: 'aparamonov@datagrok.ai'});

// Recent is read from the account's audit log (a project's creation and edits; opening it is a usage
// event), which the server writes a moment after the save.
category('Dapi: entities: recent', () => {
  test('getRecentEntities holds a project after project.open()', async () => {
    const df = grok.data.demo.demog(20);
    const project = DG.Project.create();
    project.name = `apitests-recent-${DG.Utils.randomString(8)}`;
    const tableInfo = df.getTableInfo();
    project.addChild(tableInfo);
    try {
      await grok.dapi.tables.uploadDataFrame(df);
      await grok.dapi.tables.save(tableInfo);
      await grok.dapi.projects.save(project);
      await (await grok.dapi.projects.find(project.id)).open();
      await waitUntilAsync(async () => (await grok.dapi.entities.getRecentEntities()).some((e) => e.id === project.id),
        `the project ${project.name} is not among the recent entities after it was opened`, 30000, 200);
    }
    finally {
      grok.shell.closeAll();
      const saved = await grok.dapi.projects.find(project.id).catch(() => null);
      if (saved)
        await grok.dapi.projects.delete(saved);
      const table = await grok.dapi.tables.find(tableInfo.id).catch(() => null);
      if (table)
        await grok.dapi.tables.delete(table);
    }
  }, {timeout: 60000});
}, {owner: 'aparamonov@datagrok.ai'});

async function createTestGroup(): Promise<_DG.Group> {
  let group = DG.Group.create(DG.Utils.randomString(6));
  group = await grok.dapi.groups.save(group);
  const propName = DG.Utils.randomString(5);
  const properties: any = {};
  properties[propName] = 'value';
  await group.setProperties(properties);
  return group;
}

async function safeDeleteGroup(group?: _DG.Group): Promise<void> {
  try {
    if (group)
      await grok.dapi.groups.delete(group);
  } catch (_) {}
}
