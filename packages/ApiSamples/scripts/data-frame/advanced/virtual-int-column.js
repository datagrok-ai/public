// This virtual column behaves like a normal int column.
// In this particular case, it is auto-calculated on the fly

let table = grok.data.demo.demog();
table.columns.addNewVirtual('car', (i) => i * 2, DG.TYPE.INT);
grok.shell.addTableView(table);
