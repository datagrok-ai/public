let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);

df.onColumnsAdded
  .subscribe((data) => {
    data.args.columns
      .filter((col) => col.meta.formula != null)
      .forEach((col) => grok.shell.info(`Added ${col.name}`));
  });

df.columns.addNewInt('regular column').init(5);
await df.columns.addNewCalculated('calculated column', '${x}+${y}-${z}');
grok.shell.addTableView(df);
