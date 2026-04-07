// Dynamically generated HTML-based cells

let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);
for (let c of t.columns.toList())
  view.grid.columns.byName(c.name).cellType = 'html';

view.grid.onCellPrepare(function(gc) {
  if (gc.isTableCell) 
    gc.style.element = ui.divText(gc.cell.value.toString());
});