// Adds a virtual column with the renderer defined in the PowerGrid package (make sure the package is installed).

var grid = grok.shell.addTableView(grok.data.demo.demog()).grid;
grid.columns.add({cellType: 'sparkline'});