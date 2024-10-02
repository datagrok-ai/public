//tags: DataFrame, Column, modification
let df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.FLOAT, 'x', [1, 2, 3]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'y', [4, 5, 6]),
  DG.Column.fromList(DG.TYPE.FLOAT, 'z', [7, 8, 9])
]);

let col = await df.columns.addNewCalculated('new', '0');
grok.shell.addTableView(df);

// Changes a column's formula programmatically. Note that setting a new value
// to column.meta.formula will not trigger re-calculation of values in the column.
// See also:
// https://dev.datagrok.ai/js/samples/data-frame/modification/calculated-columns/apply-formula
// https://dev.datagrok.ai/js/samples/data-frame/events/calculated-columns/update
col = await col.applyFormula('${x}+${y}-${z}');

// Opens a dialog to edit a calculated column
col.meta.dialogs.editFormula();
