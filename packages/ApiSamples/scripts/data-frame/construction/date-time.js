// Datetime columns return [Dayjs object](https://github.com/iamkun/dayjs/) but accept multiple types:

let t = DG.DataFrame.create(10)

// Batch initialization using functions
t.columns.addNewDateTime('fun dayJs').init((i) => dayjs());
t.columns.addNewDateTime('fun Date').init((i) => Date.now());
t.columns.addNewDateTime('fun ms since epoch').init((i) => 42);
t.columns.addNewDateTime('fun string').init((i) => '2024/10/30');

// Batch initialization with the same value
t.columns.addNewDateTime('dayJs').init(dayjs());
t.columns.addNewDateTime('Date').init(Date.now());
t.columns.addNewDateTime('ms since epoch').init(42);
t.columns.addNewDateTime('string').init('2024/10/30');

// Index-based initialization
let dc = t.columns.addNewDateTime('new');
dc.set(0, dayjs());
dc.set(1, Date.now());
dc.set(2, 42);
dc.set(3, '2024/10/30');

grok.shell.addTableView(t);