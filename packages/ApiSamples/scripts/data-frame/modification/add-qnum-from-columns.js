// Construct a qnum (qualified number) column from a qualifier column and a value column.

const t = DG.DataFrame.fromCsv(
`qualifier, value
<,         10
,          20
>,         30
=,         40`);

// Either column references or column names work.
t.columns.addNewQnum('result', {
  qualifierColumn: 'qualifier',
  valueColumn: t.col('value'),
});

grok.shell.addTableView(t);
