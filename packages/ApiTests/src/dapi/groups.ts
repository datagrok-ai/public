import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, test, expect} from '@datagrok-libraries/test/src/test';

category('Dapi: groups', () => {
  test('find group', async () => {
    expect((await grok.dapi.groups.filter('unexisting group').first()) == undefined);
  }, {stressTest: true});

  test('create group', async () => {
    let localTestGroup = null as any;
    try {
      const localTestGroupName = DG.Utils.randomString(7);
      localTestGroup = await grok.dapi.groups.createNew(localTestGroupName);
    } finally {
      try {
        if (localTestGroup)
          await grok.dapi.groups.delete(localTestGroup);
      } catch (_) {}
    }
  }, {stressTest: true});

  test('create subgroup', async () => {
    let group: _DG.Group | undefined;
    let subgroup: _DG.Group | undefined;
    try {
      const groupName = DG.Utils.randomString(7);
      const subgroupName = DG.Utils.randomString(7);
      group = await grok.dapi.groups.createNew(groupName);
      subgroup = await grok.dapi.groups.createNew(subgroupName);
      await grok.dapi.groups.addMember(group, subgroup);
      group = await grok.dapi.groups.find(group.id);
      expect(group.members.some((g) => g.id === subgroup!.id));
    } finally {
      if (group) {
        try {
          await grok.dapi.groups.delete(group);
        } catch (_) {}
      }
      if (subgroup) {
        try {
          await grok.dapi.groups.delete(subgroup);
        } catch (_) {}
      }
    }
  }, {stressTest: true});

  test('include member', async () => {
    let subgroup: _DG.Group | null = null;
    let demoGroup: _DG.Group | null = null;
    try {
      const localTestGroupName = `js-api-test-group1_${DG.Utils.randomString(6)}`;
      const localTestGroup2Name = `js-api-test-group2_${DG.Utils.randomString(6)}`;
      demoGroup = await grok.dapi.groups.createNew(localTestGroup2Name);
      subgroup = DG.Group.create(localTestGroupName);
      subgroup.includeTo(demoGroup);
      const adminUser = await grok.dapi.users.filter('login = "admin"').first();
      subgroup.addAdminMember(adminUser.group);
      await grok.dapi.groups.saveRelations(subgroup);
      subgroup = await grok.dapi.groups.include('children.child').filter(`shortName="${localTestGroupName}"`).first();

      let hasAdmin = false;
      for (const m of subgroup.adminMembers) {
        if (m.friendlyName.toLowerCase() == 'admin')
          hasAdmin = true;
      }

      if (!hasAdmin)
        throw new Error('Member not added');
    } finally {
      try {
        if (demoGroup)
          await grok.dapi.groups.delete(demoGroup);
      } catch (_) {}
      try {
        if (subgroup)
          await grok.dapi.groups.delete(subgroup);
      } catch (_) {}
    }
  }, {stressTest: true});

  test('delete group', async () => {
    const localTestGroupName = 'js-api-test-group1';
    const countBefore = await grok.dapi.groups.filter(localTestGroupName).count;
    const localTestGroup = await grok.dapi.groups.createNew(localTestGroupName);
    expect((await grok.dapi.groups.filter(localTestGroupName).first())?.name, localTestGroupName);
    await grok.dapi.groups.delete(localTestGroup);
    expect((await grok.dapi.groups.filter(localTestGroupName)).count == countBefore);
  }, {stressTest: true});

}, {owner: 'aparamonov@datagrok.ai'});
