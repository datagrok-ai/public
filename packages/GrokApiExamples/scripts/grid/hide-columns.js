// Hiding columns

let table = gr.testData('demog', 5000);
let view = gr.addTableView(table);

view.grid.columns.setVisible(['age', 'sex', 'race']);

// or hiding by adding '~' prefix to column name

table.columns.byName('sex').name = '~sex';
