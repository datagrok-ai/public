// TabControl

let v = grok.shell.newView('Tabs');

// Fluent API
v.append(ui.h1('Horizontal'));
v.append(
    ui.tabControl({
        'FIRST': () => ui.divText('First panel'),
        'SECOND': () => ui.divText('Second panel')
    })
);

// Detailed API
v.append(ui.h1('Vertical'));
let tabsV = ui.tabControl(null, true);
tabsV.addPane('FIRST', () => ui.divText('First panel'));
tabsV.addPane('SECOND', () => ui.divText('Second panel'));
v.append(tabsV);
