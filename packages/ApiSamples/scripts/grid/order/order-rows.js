// Note that you can also filter visible rows using this method.

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.onRowsSorted.subscribe((ev) => grok.shell.info('Sorted'));
view.grid.setRowOrder([1, 56, 3, 6, 4]);

