let manual = DG.Column.qnum('manual', 3);
manual.set(0, DG.Qnum.greater(5));
manual.set(1, DG.Qnum.exact(5));
manual.set(1, DG.Qnum.less(5));

// qualifier is stripped
let fromExact = DG.Column.qnum('exact', 3, [1.123, 2.2345, 3.4567], true);

// the platform expects qualified numbers
let fromQnum = DG.Column.qnum('qnum', 3, [1.123, 2.2345, 3.4567], false);

let t = DG.DataFrame.fromColumns([ manual, fromExact, fromQnum ]);
grok.shell.addTableView(t);
