//name: event-on-grid-resize
//tags: demo
//language: javascript
//
let t = grok.data.demo.demog()
let view = grok.shell.addTableView(t);

view.grid.onRowResized.subscribe(function (ev) {
  if (ev.args.dragging)
    console.log("Row resize in progress");
  else
    grok.shell.info("Row Resized!");
});


view.grid.onColumnResized.subscribe(function (ev) {
  if (ev.args.dragging)
    console.log("Col resize in progress");
  else
    grok.shell.info("Column resized!");
});


