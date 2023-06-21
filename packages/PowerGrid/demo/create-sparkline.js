// Adds a "radar" sparkline
// Requires PowerGrid package

const tv = grok.shell.addTableView(grok.data.demo.biosensor());
const gc = tv.grid.columns.add({gridColumnName: 'radar', cellType: 'radar'});
gc.settings = {columnNames: ['x', 'y', 'z']};
