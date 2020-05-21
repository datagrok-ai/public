let t = grok.DataFrame.fromColumns([
    grok.Column.fromList('int', 'int', [1, 2, 3]),
    grok.Column.fromList('double', 'double', [1.1, 2.1, 3.1]),
    grok.Column.fromList('string', 'string', ["a", "b", "c"])
]);

grok.addTableView(t);
