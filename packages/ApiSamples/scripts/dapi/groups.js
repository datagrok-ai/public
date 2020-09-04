async function groupsTest() {
//get list of groups named "demo"
    let demoGroup = await grok.dapi.groups.filter('demo').first();

//create a new Demo group, if doesn't exist and save
    if (demoGroup == null)
        demoGroup = await grok.dapi.groups.createNew("Demo");

//create a subgroup
    let demoSubGroup = DG.Group.create(`Demo ${Math.floor(Math.random() * 100)}`);

//add to "Demo" group
    demoSubGroup.includeTo(demoGroup);

//find "admin" user
    let adminUser = await grok.dapi.users.filter('login = "admin"').first();

//add to subgroup and save
    demoSubGroup.addAdminMember(adminUser.group);
    await grok.dapi.groups.saveRelations(demoSubGroup);

}

groupsTest();