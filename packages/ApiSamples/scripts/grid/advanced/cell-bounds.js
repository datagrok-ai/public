let table = grok.data.demo.demog();
let view = grok.shell.addTableView(table);
rxjs.interval(1000).pipe(rxjs.operators.startWith(0)).subscribe((i) => {
  let r = view.grid.cell('site', i).bounds;
  view.grid.canvas.getContext('2d').fillRect(r.x, r.y, r.width, r.height);
});