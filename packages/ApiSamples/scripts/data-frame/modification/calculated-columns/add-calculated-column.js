let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);

// Creates a calculated column programmatically
let col = await df.columns.addNewCalculated('new', '${x}+${y}-${z}');
grok.shell.addTableView(df);

// Opens a dialog to create a calculated column
df.dialogs.addNewColumn();
