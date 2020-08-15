// Resize a column.

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.columns.byName('age').width = 200;
