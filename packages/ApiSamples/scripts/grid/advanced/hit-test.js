// Use Grid.hitTest to find a cell at the specified coordinates

let table = grok.data.demo.demog();
let view = grok.shell.addTableView(table);
rxjs.fromEvent(view.grid.overlay, 'mousemove').subscribe((mm) => {
  let cell = view.grid.hitTest(mm.offsetX, mm.offsetY);
  grok.shell.o = cell.gridColumn;
});