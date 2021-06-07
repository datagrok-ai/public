// TabControl

let v = grok.shell.newView('Tabs');

// Simple API
v.append(ui.h1('Horizontal'));
v.append(
    ui.tabControl({
        'FIRST': ui.panel('First panel', {style :{backgroundColor: '#EFE'}}),           // an element
        'SECOND': () => ui.panel('Second panel', {style :{backgroundColor: '#EEF'}})    // a function that returns an element
    })
);

// Detailed API
v.append(ui.h1('Vertical'));
let tabsV = ui.tabControl(null, true);
tabsV.addPane('FIRST', () => ui.panel('First panel', {style :{backgroundColor: '#EFE'}}));
tabsV.addPane('SECOND', () => ui.panel('Second panel', {style :{backgroundColor: '#EEF'}}));
v.append(tabsV);
