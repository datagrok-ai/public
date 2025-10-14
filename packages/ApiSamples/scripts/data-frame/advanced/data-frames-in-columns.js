let dataFrameColumn = DG.Column.dataFrame('tables', 5);
let t = DG.DataFrame.fromColumns([dataFrameColumn]);
dataFrameColumn.set(1, grok.data.demo.demog());
let grid = grok.shell.addTableView(t).grid;

grid.col('tables').width = 500;
grid.props.rowHeight = 300;