// Custom view toolbox

let v = grok.shell.newView('toolbox demo', [ui.divText('See custom toolbox on the left.')]);

let acc = ui.accordion();
acc.addPane('header 1', () => ui.divText('Dynamic content'));
acc.addPane('header 2', () => ui.divText('More dynamic content'));

v.toolbox = acc.root;
