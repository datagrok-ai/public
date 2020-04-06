let t = DataFrame.fromColumns([
    Column.fromList('int', 'int', [1, 2, 3]),
    Column.fromList('double', 'double', [1.1, 2.1, 3.1]),
    Column.fromList('string', 'string', ["a", "b", "c"])
]);

grok.addTableView(t);
