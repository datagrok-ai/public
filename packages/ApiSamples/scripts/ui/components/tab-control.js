// TabControl

let v = grok.shell.newView('Tabs');

let tabControl = ui.tabControl({
  'FIRST': ui.panel('First panel', {style: {backgroundColor: '#EFE'}}), // an element
  'SECOND': () => ui.panel('Second panel', {style: {backgroundColor: '#EEF'}}) // a function that returns an element
});

// Simple API
v.append(ui.h1('Horizontal'));
v.append(tabControl);

// Detailed API
v.append(ui.h1('Vertical'));
let tabsV = ui.tabControl(null, true);
tabsV.addPane('FIRST', () => ui.panel('First panel', {style: {backgroundColor: '#EFE'}}));
tabsV.addPane('SECOND', () => ui.panel('Second panel', {style: {backgroundColor: '#EEF'}}));
v.append(tabsV);

// events
tabControl.onTabChanged.subscribe((_) => grok.shell.info(tabControl.currentPane.name));
