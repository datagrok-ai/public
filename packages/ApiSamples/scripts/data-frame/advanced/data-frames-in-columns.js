let dataFrameColumn = DG.Column.dataFrame('tables', 5);
let t = DG.DataFrame.fromColumns([dataFrameColumn]);
dataFrameColumn.set(1, grok.data.demo.demog());
grok.shell.addTableView(t);