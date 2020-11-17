// TabControl

let v = grok.shell.newView('Tabs');

// Simple API
v.append(ui.h1('Horizontal'));
v.append(
    ui.tabControl({
        'FIRST': ui.divText('First panel'),           // an element
        'SECOND': () => ui.divText('Second panel')    // a function that returns an element
    })
);

// Detailed API
v.append(ui.h1('Vertical'));
let tabsV = ui.tabControl(null, true);
tabsV.addPane('FIRST', () => ui.divText('First panel'));
tabsV.addPane('SECOND', () => ui.divText('Second panel'));
v.append(tabsV);
