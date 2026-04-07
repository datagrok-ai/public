//help-url: https://datagrok.ai/help/develop/how-to/layouts
async function layoutsTest() {
  let allLayouts = await grok.dapi.layouts.list();
  grok.shell.info(`all layouts: ${allLayouts.length}`);

  await grok.dapi.entities.saveProperties([
    {'entityId': allLayouts[0].id, 'property': 'prop', 'value': 'mark'},
    {'entityId': allLayouts[1].id, 'property': 'prop', 'value': 'mark2'},
  ]);

  await grok.dapi.layouts.save(allLayouts[1]);

  let markedLayouts = await grok.dapi.layouts.filter('prop = "mark2"').list();
  grok.shell.info(`marked layouts: ${markedLayouts.length}`);

  let adminUser = await grok.dapi.users.filter('login = "admin"').first();
  await grok.dapi.permissions.grant(allLayouts[1], adminUser.group, true);
  console.log(await grok.dapi.permissions.get(allLayouts[1]));
}

layoutsTest();