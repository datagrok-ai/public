// Hiding columns

let table = grok.data.demo.demog();
let view = grok.shell.addTableView(table);

view.grid.columns.setVisible(['age', 'site', 'race']);

// Hiding columns by adding '~' prefix to column names affects all views

table.columns.byName('sex').name = '~sex';
grok.shell.addTableView(table);
