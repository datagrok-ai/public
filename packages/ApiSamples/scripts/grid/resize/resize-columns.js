// Resize a column.

let view = grok.shell.addTableView(grok.data.demo.demog());

view.grid.columns.byName('age').width = 200;
view.grid.columns.byIndex(4).width = 300;
view.grid.columns.rowHeader.width = 100;
