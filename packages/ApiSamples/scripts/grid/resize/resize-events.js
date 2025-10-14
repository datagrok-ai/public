// Handling of column and row resizing

let grid = grok.shell.addTableView(grok.data.demo.demog()).grid;

grid.onRowsResized.subscribe((ev) => {
  grok.shell.info("Resizing row height: " + (ev.args.dragging ? "in progress" : "done"));
});

grid.onColumnResized.subscribe((ev) => {
  grok.shell.info(`Resizing ${ev.args.column.name}: ` + (ev.args.dragging ? "in progress" : "done"));
});