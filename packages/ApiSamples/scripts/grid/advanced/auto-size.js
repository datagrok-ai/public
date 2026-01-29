// Grid.autoSize resizes the grid to fit the content

const grid = grok.data.demo.demog(5).plot.grid();
grid.autoSize(300, 300);
ui.dialog()
  .add(grid.root)
  .show();
