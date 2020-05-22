// Resize a column.

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));
view.grid.columns.byName('age').width = 200;