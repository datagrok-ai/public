// Hiding columns

let table = grok.testData('demog', 5000);
let view = grok.addTableView(table);

view.grid.columns.setVisible(['age', 'sex', 'race']);

// or hiding by adding '~' prefix to column name

table.columns.byName('sex').name = '~sex';
