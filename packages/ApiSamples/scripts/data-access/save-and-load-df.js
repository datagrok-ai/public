// Saves a data frame remotely, loads it using id, and shows the data frame info

let df = DG.DataFrame.fromColumns([DG.Column.fromStrings('myColumn', ['first row', 'second row', 'third row'])]);
df.name = 'my data frame';
let id = await grok.dapi.tables.uploadDataFrame(df);

grok.shell.info('Data frame id: ' + id);
let loadedDf = await grok.dapi.tables.getTable(id);
grok.shell.addTableView(loadedDf);
let meta = await grok.dapi.tables.find(id);
grok.shell.info('Data frame name: ' + meta.name);
