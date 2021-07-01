// grid.getVisibleBounds

let grid = grok.shell.addTableView(grok.data.demo.demog()).grid;

setInterval(() => {
  let g = grid.canvas.getContext('2d');

  g.strokeStyle = 'green';
  for (let cell of grid.getVisibleCells())
    g.strokeRect(cell.bounds.x, cell.bounds.y, cell.bounds.width, cell.bounds.height);

  g.strokeStyle = 'blue';
  for (let cell of grid.getVisibleCells(grid.columns.byName('age')))
    g.strokeRect(cell.bounds.x, cell.bounds.y, cell.bounds.width, cell.bounds.height);

}, 1000);