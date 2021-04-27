// Handles onCellTooltip event to customize the tooltip for "age" column and all column headers

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.onCellTooltip(function (cell, x, y) {
  if (cell.isTableCell && cell.tableColumn.name === 'age') {
    ui.tooltip.show(ui.divV([
      ui.h1('Custom tooltip'),
      ui.divText(cell.tableRow.race)
    ]), x, y);
    return true;
  }
  if (cell.isColHeader) {
    ui.tooltip.show(ui.divV([ui.h1(`Custom header for ${cell.tableColumn.name}`)]), x, y);
    return true;
  }
});
