// Note that the "population" data type becomes int

let t = DG.DataFrame.create(3);
t.columns.add(DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
t.columns.add(DG.Column.fromStrings('population', ['321', '35', '121']));
grok.shell.addTableView(t);

let t2 = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']),
    DG.Column.fromStrings('population', ['321', '35', '121'])
]);
grok.shell.addTableView(t2);
