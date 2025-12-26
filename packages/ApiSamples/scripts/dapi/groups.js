//help-url: https://datagrok.ai/help/govern/group
async function groupsTest() {

  // Get list of groups named "demo"
  let demoGroup = await grok.dapi.groups.filter('demo').first();

  // Create a new Demo group, if doesn't exist and save
  if (demoGroup == null)
    demoGroup = await grok.dapi.groups.createNew('Demo Group');

  // Create a subgroup
  let demoSubGroup = DG.Group.create('Demo Sub Group');

  // Add to "Demo" group
  demoSubGroup.includeTo(demoGroup);

  // Find "admin" user
  let adminUser = await grok.dapi.users.filter('login = "admin"').first();

  // Add to subgroup and save
  demoSubGroup.addAdminMember(adminUser.group);
  await grok.dapi.groups.saveRelations(demoSubGroup);

  // Delete subgroup
  await grok.dapi.groups.delete(demoSubGroup);

  // Delete Demo group
  await grok.dapi.groups.delete(demoGroup);
}

groupsTest();