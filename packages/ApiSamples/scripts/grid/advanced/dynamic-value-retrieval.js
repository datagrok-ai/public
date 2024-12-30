// Adds a virtual column with the renderer defined in the PowerGrid package (make sure the package is installed).

const grid = grok.shell.addTableView(grok.data.demo.demog()).grid;
grid.columns.byName('age').onPrepareValueScript = 'return 42;';
grid.columns.byName('sex').onPrepareValueScript = 'return "42";';

//grid.columns.add({gridColumnName: 'table', cellType: 'dataframe'})
//  .onPrepareValueScript = 'return await grok.dapi.files.readCsv("System:DemoFiles/demog.csv")';

grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
  .onPrepareValueScript = 'return await grok.dapi.files.readCsv("System:DemoFiles/hts/plate-96-1.csv")';