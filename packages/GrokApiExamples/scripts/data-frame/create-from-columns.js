// Note that the "population" data type becomes int

let t = grok.DataFrame.create(3);
t.columns.add(grok.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
t.columns.add(grok.Column.fromStrings('population', ['321', '35', '121']));
grok.addTableView(t);

let t2 = grok.DataFrame.fromColumns([
    grok.Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']),
    grok.Column.fromStrings('population', ['321', '35', '121'])
]);
grok.addTableView(t2);
