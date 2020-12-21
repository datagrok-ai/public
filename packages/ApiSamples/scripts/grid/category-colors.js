// Assigning custom category colors

let view = grok.shell.addTableView(grok.data.demo.demog());

view.grid.col('sex').categoryColors = {
  'M': 0xFF0000FF,
  'F': 0xFF800080
};
