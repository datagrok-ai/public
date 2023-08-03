let t = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);
t.rows.removeWhere((row) => row.get('make') === 'Honda');
grok.shell.addTableView(t);