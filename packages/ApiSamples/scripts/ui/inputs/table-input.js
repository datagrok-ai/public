grok.shell.addTableView(grok.data.demo.demog());
grok.shell.addTableView(grok.data.demo.randomWalk());

let tables = grok.shell.tables;
let input = ui.input.table('Table', {value: tables[0], items: tables, onValueChanged: (t) => grok.shell.info(t.name)});

grok.shell.newView('tableInput', [input.root]);