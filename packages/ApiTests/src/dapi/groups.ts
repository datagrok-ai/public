import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, test, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';

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

  test('request membership', async () => {
    let target: _DG.Group | undefined;
    try {
      target = await grok.dapi.groups.createNew(`js-api-test-target_${DG.Utils.randomString(6)}`);
      await grok.dapi.groups.requestMembership(target, grok.shell.user.group);
    } finally {
      try {
        if (target)
          await grok.dapi.groups.delete(target);
      } catch (_) {}
    }
  });

  test('current user groups', async () => {
    const groups = await grok.dapi.groups.currentUserGroups();
    expect(Array.isArray(groups) && groups.length > 0, true);
    expect(groups.some((g) => g.friendlyName === 'All users'), true);
  });

  test('delete group', async () => {
    // Random name + exact-name filter so concurrent stress runs never collide
    // on the group name (a hard-coded name caused "already exists" / wrong-group
    // assertions under load); always delete it afterwards.
    const localTestGroupName = `js-api-test-delete_${DG.Utils.randomString(8)}`;
    let localTestGroup: _DG.Group | null = await grok.dapi.groups.createNew(localTestGroupName);
    try {
      expect((await grok.dapi.groups.filter(`shortName="${localTestGroupName}"`).first())?.name, localTestGroupName);
      await grok.dapi.groups.delete(localTestGroup);
      localTestGroup = null;
      expect((await grok.dapi.groups.filter(`shortName="${localTestGroupName}"`).first()) == undefined);
    } finally {
      try {
        if (localTestGroup)
          await grok.dapi.groups.delete(localTestGroup);
      } catch (_) {}
    }
  }, {stressTest: true});

  test('a name a group or role already has is refused', async () => {
    const name = `js-api-test-taken_${DG.Utils.randomString(8)}`;
    await grok.dapi.groups.createNew(name);
    try {
      await expectExceptionAsync(() => grok.dapi.groups.createNew(name),
        (e) => String(e?.message ?? e).includes('A group or role with the same name already exists'));
      expect((await grok.dapi.groups.filter(`shortName="${name}"`).list()).length, 1);
    } finally {
      for (const group of await grok.dapi.groups.filter(`shortName="${name}"`).list())
        await grok.dapi.groups.delete(group);
    }
  });

  // the account that creates a group is its only admin
  test('the only admin cannot leave a group or lose Admin', async () => {
    const me = grok.shell.user.group;
    const group = await grok.dapi.groups.createNew(`js-api-test-admin_${DG.Utils.randomString(8)}`);
    const refused = (e: any) => String(e?.message ?? e).includes('should have at least one admin user');
    try {
      await expectExceptionAsync(() => grok.dapi.groups.removeMember(group, me), refused);
      await expectExceptionAsync(() => grok.dapi.groups.addMember(group, me), refused);
      const saved = await grok.dapi.groups.include('children.child').find(group.id);
      expect(saved.adminMembers.some((m) => m.id === me.id), true);
    } finally {
      await grok.dapi.groups.delete(group);
    }
  });
}, {owner: 'aparamonov@datagrok.ai'});
