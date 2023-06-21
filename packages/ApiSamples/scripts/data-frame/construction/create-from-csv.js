//tags: DataFrame, construction
// Creating a DataFrame from a CSV string

let table = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
BMW,   328i,     4,         1.7,    60000        
BMW,   535i,     6,         1.5,    35000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`, {columnImportOptions: [{name: 'model', semType: 'carModel'}]});

grok.shell.addTableView(table);
