async function layoutsTest() {
    let allLayouts = await grok.dapi.layouts.list();
    grok.shell.info(`all layouts: ${allLayouts.length}`);

    await grok.dapi.entities.saveProperties([
        {'entityId': allLayouts[0].id, 'property': 'prop', 'value': 'mark'},
        {'entityId': allLayouts[1].id, 'property': 'prop', 'value': 'mark2'},
    ]);

    allLayouts[1].id = 'cddcabeb-70ad-42da-a526-dfb187784ffd';

    await grok.dapi.layouts.save(allLayouts[1]);

    let markedLayouts = await grok.dapi.layouts.filter('prop = "mark2"').list();
    grok.shell.info(`marked layouts: ${markedLayouts.length}`);

}

layoutsTest();