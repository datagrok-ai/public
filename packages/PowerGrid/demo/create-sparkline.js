// Adds a "radar" sparkline
// Requires PowerGrid package

let tv = grok.shell.addTableView(grok.data.demo.biosensor());
let gc = tv.grid.columns.add({gridColumnName: 'radar', cellType: 'radar'});
gc.settings = {columnNames: ['x', 'y', 'z']};