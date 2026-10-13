//api: DG.ColumnList.addNewCalculated
// What a calculated column does with the rows its formula fails on.
// By default they stay empty. The behavior is kept with the column: layouts and recalculations reuse it.
let df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('sample date', ['2023-01-05', 'n/a', '2023-03-07', 'pending']),
]);

// Fills the failed rows with a value, and adds 'parsed errors' with each failed row's message
await df.columns.addNewCalculated('parsed', 'DateParse(${sample date})', {
  type: 'datetime',
  onError: {mode: 'value', value: dayjs.utc('1900-01-01'), errorColumn: true},
});
grok.shell.addTableView(df);

// Rejects on the first failed row, and adds no column
try {
  await df.columns.addNewCalculated('strict', 'DateParse(${sample date})', {onError: {mode: 'stop'}});
}
catch (e) {
  grok.shell.warning(`Not added: ${e}`);
}
