// Resize a column.

let view = gr.addTableView(gr.testData('demog', 5000));
view.grid.columns.byName('age').width = 200;