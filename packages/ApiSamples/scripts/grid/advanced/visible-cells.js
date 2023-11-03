// grid.getVisibleBounds

let grid = grok.shell.addTableView(grok.data.demo.demog()).grid;
grid.onAfterDrawContent.subscribe(() => setTimeout(() => {
  let g = grid.canvas.getContext('2d');

  g.strokeStyle = 'green';
  for (let cell of grid.getVisibleCells())
    g.strokeRect(cell.bounds.x, cell.bounds.y, cell.bounds.width, cell.bounds.height);

  g.strokeStyle = 'blue';
  for (let cell of grid.columns.byName('age').getVisibleCells())
    g.strokeRect(cell.bounds.x, cell.bounds.y, cell.bounds.width, cell.bounds.height);
  
}, 300));