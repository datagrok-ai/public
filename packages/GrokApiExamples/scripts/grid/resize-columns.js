// Resize a column.

let view = grok.addTableView(grok.testData('demog', 5000));
view.grid.columns.byName('age').width = 200;