// Different ways to add columns

let t = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);

// specify type identifier as a parameter
t.columns.addNew('foo', DG.TYPE.INT);

// Use type-specific methods if you know the data type in advance.
// Also, check out different ways of column initialization.
t.columns.addNewInt('int').init(5); // scalar initializer
t.columns.addNewFloat('float').init((i) => i / 10); // function initializer
t.columns.addNewBool('bool');
t.columns.addNewString('string');
t.columns.addNewDateTime('datetime');
t.columns.addNewQnum('qnum');

grok.shell.addTableView(t);