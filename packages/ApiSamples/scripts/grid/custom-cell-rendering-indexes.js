// Grid row indexes vs table row indexes

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.onCellRender.subscribe(function (args) {
  args.g.fillStyle = (args.cell.isColHeader ? 'red' : (args.cell.isRowHeader ? 'green' : 'blue'));
  args.g.fillText(args.cell.isColHeader ? args.cell.gridColumn.name : args.cell.tableRowIndex, args.bounds.x + args.bounds.width / 2, args.bounds.y + args.bounds.height / 2);
  args.preventDefault();
});