// A new instance of the column is created when column type is changed.

// Metod 1: column.convertTo(newType)
let strCol = DG.Column.fromList('string', 'x', ['5', '4', '3']);
let intCol = strCol.convertTo(DG.COLUMN_TYPE.INT);
intCol.name = 'foo';
let t = DG.DataFrame.fromColumns([strCol, intCol]);

// Metod 2: table.changeColumnType(columnIt, newType)
t.changeColumnType('foo', 'qnum');

grok.shell.addTableView(t);