let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);
let col = await df.columns.addNewCalculated('new', '${x}+${y}-${z}');
let view = grok.shell.addTableView(df);

// Listen to formula changes of a calculated column
df.onMetadataChanged
  .pipe(rxjs.operators.filter(data => data.args.key === DG.TAGS.FORMULA))
  .subscribe(async (data) => {
    grok.shell.info(`${data.args.change} - ${data.args.key} - ${data.args.value}`);
    if (data.args.change === 'set' && data.args.source === col) {
      // Apply a formula to the column (returns a new column instance preserving the metadata)
      col = await data.args.source.applyFormula(data.args.value);
    }
  });

// Change the formula either programmatically or from the interface (tab 'Formula' in the context panel)
col.meta.formula = 'Avg($[x])';
