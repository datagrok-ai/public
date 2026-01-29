// API References:
//
// https://datagrok.ai/js-api/dg/classes/ColumnList
// https://datagrok.ai/js-api/dg/classes/DataFrame
// https://datagrok.ai/js-api/dg/classes/Column#fromList
// https://datagrok.ai/js-api/dg/classes/Column#qnum
//
// Note that the "population" data type becomes int

let t = DG.DataFrame.create(3);
t.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
t.columns.add(DG.Column.fromStrings('population', ['321', '35', '121']));
grok.shell.addTableView(t);

let t2 = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']),
  DG.Column.fromStrings('population', ['321', '35', '121']),
  DG.Column.fromType(DG.COLUMN_TYPE.FLOAT, 'float', 3),
  DG.Column.qnum('qnum', 3),
]);
grok.shell.addTableView(t2);
