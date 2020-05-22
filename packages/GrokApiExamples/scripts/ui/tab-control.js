// TabControl

let v = grok.shell.newView('Tabs');

v.root.appendChild(ui.h1('Horizontal'));
let tabs = ui.tabControl();
tabs.addPane('FIRST', () => ui.divText('First panel'));
tabs.addPane('SECOND', () => ui.divText('Second panel'));
v.root.appendChild(tabs.root);

v.root.appendChild(ui.h1('Vertical'));
let tabsV = ui.tabControl(true);
tabsV.addPane('FIRST', () => ui.divText('First panel'));
tabsV.addPane('SECOND', () => ui.divText('Second panel'));
v.root.appendChild(tabsV.root);