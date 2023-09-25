grok.shell.addTableView(grok.data.demo.demog());
grok.shell.addTableView(grok.data.demo.randomWalk());

let tables = grok.shell.tables;
let input = ui.tableInput('Table', tables[0], tables, (t) => grok.shell.info(t.name));

grok.shell.newView('tableInput', [input.root]);