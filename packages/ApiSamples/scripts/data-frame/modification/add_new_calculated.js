let df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);
df.columns.addNewCalculated('new', '${x}+${y}-${z}').then((_) => {
    let view = grok.shell.addTableView(df);
}) ;
