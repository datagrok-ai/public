// Saves a data frame remotely, loads it using id, and shows the data frame info

let df = DG.DataFrame.fromColumns([DG.Column.fromStrings('myColumn', ['first row', 'second row', 'third row'])]);
df.name = 'my data frame';
grok.dapi.tables.uploadDataFrame(df).then((id) => {
    grok.shell.info('Data frame id: ' + id);
    grok.dapi.tables.getTable(id).then((loadedDf) => grok.shell.addTableView(loadedDf));

    grok.dapi.tables.find(id).then((meta) => grok.shell.info('Data frame name: ' + meta.name));
});
