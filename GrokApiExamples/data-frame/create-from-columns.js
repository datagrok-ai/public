// Note that the "population" data type becomes int

let t = DataFrame.create(3);
t.cols.add(Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']));
t.cols.add(Column.fromStrings('population', ['321', '35', '121']));
gr.addTableView(t);

let t2 = DataFrame.fromColumns([
    Column.fromStrings('countries', ['USA', 'Canada', 'Mexico']),
    Column.fromStrings('population', ['321', '35', '121'])
]);
gr.addTableView(t2);
