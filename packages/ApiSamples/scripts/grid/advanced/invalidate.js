let counts = {};
let grid = grok.shell.addTableView(grok.data.demo.demog()).grid;
grid.onCellPrepare(function(gc) {
  let xy = `${gc.gridColumn.name}_${gc.gridRow}`;
  counts[xy] = counts[xy] == null ? 1 : counts[xy] + 1;
  gc.customText = `${counts[xy]}`;
});

setInterval(() => {grid.invalidate();}, 1000);