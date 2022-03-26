// Reusing core grid renderers

let view = grok.shell.addTableView(grok.data.demo.molecules());
view.dataFrame.col('logD').colors.setLinear();
let canvas = ui.canvas(300, 200);

view.grid.onCurrentCellChanged.subscribe((gridCell) =>
  gridCell.renderer.render(canvas.getContext('2d'), 20, 20, 200, 100, gridCell, gridCell.style));

ui.dialog()
  .add('Click on a cell to render it in the canvas below')
  .add(canvas)
  .show();