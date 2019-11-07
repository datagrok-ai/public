// Dynamically generated HTML-based cells

let view = gr.addTableView(gr.testData('demog', 5000));
let col = view.grid.columns.byName('disease');
view.grid.options({'rowHeight': 100});
col.width = 200;
col.cellType = 'html';

view.grid.onCellPrepare(function (gc) {
    if (gc.isTableCell && gc.gridColumn.name == 'disease') {
        //debugger;
        gc.style.element = ui.divV([
            ui.h1(gc.tableRow.subj),
            ui.button('CONTACT')
        ]);
    }
});