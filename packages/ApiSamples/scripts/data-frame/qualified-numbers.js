let qCol = DG.Column.qnum(DG.COLUMN_TYPE.QNUM, 3);
let t = DG.DataFrame.fromColumns([ qCol ]);
qCol.set(0, DG.Qnum.greater(5));

grok.shell.addTableView(t);
