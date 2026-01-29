let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);

await df.columns.addNewCalculated('new', '${x}+${y}-${z}');
grok.shell.addTableView(df);

df.onColumnsRemoved
  .subscribe((data) => {
    data.args.columns
      .filter((col) => col.meta.formula != null)
      .forEach((col) => grok.shell.info(`Removed calculated column ${col.name}`));
  });

df.columns.remove('x');
df.columns.remove('new');
