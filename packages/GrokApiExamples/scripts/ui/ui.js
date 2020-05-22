// Accordion allows you to add panels which content is created
// dynamically during the expansion.

let acc = ui.accordion();
acc.addPane('pane 1', () => ui.divText(new Date().toTimeString()));
acc.addPane('pane 2', () => ui.divText(new Date().toTimeString()));

grok.shell.newView('tables', [acc.root]);