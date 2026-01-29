let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);
let cols = df.columns;
// Removes the second column, named 'y'
cols.remove('y');
let rows = df.rows;
// Removes the second row
rows.removeAt(1);
let view = grok.shell.addTableView(df);