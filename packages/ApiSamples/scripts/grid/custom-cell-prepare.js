// Handles onCellPrepare event to set GridCell's parameters

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.onCellPrepare(function (gc) {
  // let's not process row headers, column headers, and other special cells
  if (!gc.isTableCell)
    return;

  if (gc.gridColumn.idx > 4 && (gc.gridRow + gc.gridColumn.idx) % 2 === 0)
    gc.style.backColor = 0xff1f77b4;

  if (gc.gridColumn.name === 'race')
    gc.style.font = '14px MarkerFelt-Thin, Comic Sans MS';

  if (gc.cell.value === 'M')
    gc.customText = 'male';

  if (gc.gridColumn.name === 'age' && gc.cell.value > 30)
    gc.style.textColor = 0xFFFF0000;
});
