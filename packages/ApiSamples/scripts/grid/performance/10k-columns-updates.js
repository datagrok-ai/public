let table = grok.data.testData('random walk', 100, 10000);
let grid = grok.shell.addTableView(table).grid;

for (let i = 0; i < table.columns.length; i++)
  grid.columns.byIndex(i).width = 10;

