import {after, before, category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Dapi: groups', () => {
  const testGroupName = 'js-api-test-group';
  let testGroup : DG.Group;

  before(async () => {
    testGroup = await grok.dapi.groups.createNew(testGroupName);
  });

  test('Dapi: groups - find group', async () => {
    if ((await grok.dapi.groups.filter(testGroupName).first()) === undefined)
      throw 'Group doesn\'t exist';
  });

  test('Dapi: groups - create group', async () => {
    let localTestGroup = null as any;
    try {
      const localTestGroupName = 'js-api-test-group1';
      localTestGroup = await grok.dapi.groups.createNew(localTestGroupName);
    } finally {
      await grok.dapi.groups.delete(localTestGroup);
    }
  });

  test('Dapi: groups - create subgroup', async () => {
    let subgroup = null as any;
    try {
      const subgroupName = 'js-api-test-group1';
      subgroup = DG.Group.create(subgroupName);
      testGroup.includeTo(subgroup);
    } finally {
      await grok.dapi.groups.delete(subgroup);
    }
  });

  test('Dapi: groups - include member', async () => {
    let subgroup = null as any;
    let demoGroup = null as any;
    try {
      const localTestGroupName = 'js-api-test-group1';
      const localTestGroup2Name = 'js-api-test-group2';
      demoGroup = await grok.dapi.groups.createNew(localTestGroup2Name);
      subgroup = DG.Group.create(localTestGroupName);
      subgroup.includeTo(demoGroup);
      const adminUser = await grok.dapi.users.filter('login = "admin"').first();
      subgroup.addAdminMember(adminUser.group);
      await grok.dapi.groups.saveRelations(subgroup);

      subgroup = await grok.dapi.groups.include('children.child').filter(localTestGroupName).first();

      let hasAdmin = false;
      for (const m of subgroup.adminMembers) {
        if (m.friendlyName.toLowerCase() == 'admin')
          hasAdmin = true;
      }

      if (!hasAdmin)
        throw 'Member not added';
    } finally {
      await grok.dapi.groups.delete(demoGroup);
      await grok.dapi.groups.delete(subgroup);
    }
  });

  test('Dapi: groups - delete group', async () => {
    let localTestGroup = null as any;
    try {
      const localTestGroupName = 'js-api-test-group1';

      localTestGroup = await grok.dapi.groups.createNew(localTestGroupName);
      await grok.dapi.groups.delete(localTestGroup);

      if (await grok.dapi.groups.filter(localTestGroupName).first() !== undefined)
        throw 'Group not deleted';
    } finally {
      await grok.dapi.groups.delete(localTestGroup);
    }
  });

  after(async () => {
    await grok.dapi.groups.delete(testGroup);
  });
});
