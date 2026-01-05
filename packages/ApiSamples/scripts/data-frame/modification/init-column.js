// Different ways to initialize the content of the column

let t = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);

// initialize to scalar
t.col('make').init('VAZ');
t.col('cylinders').init(4);

// index-based initializer
t.columns.addNewInt('index').init((i) => i * 2);

grok.shell.addTableView(t);
