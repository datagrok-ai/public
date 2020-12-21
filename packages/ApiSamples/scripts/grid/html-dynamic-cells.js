// Dynamically generated HTML-based cells

let view = grok.shell.addTableView(grok.data.demo.demog());
let col = view.grid.columns.byName('disease');
view.grid.setOptions({'rowHeight': 100});
col.width = 200;
col.cellType = 'html';

view.grid.onCellPrepare(function (gc) {
  if (gc.isTableCell && gc.gridColumn.name === 'disease') {
    //debugger;
    gc.style.element = ui.divV([
      ui.h1(gc.tableRow.subj),
      ui.button('CONTACT')
    ]);
  }
});
