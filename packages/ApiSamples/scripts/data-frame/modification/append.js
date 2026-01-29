// Appends two tables. Note that the missing columns are created automatically.

let t1 = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);

let t2 = DG.DataFrame.fromCsv(
  `make,   volume, model,  price
BMW,     1.7,    328i,   60000        
BMW,     1.5,    535i,   35000`);

grok.shell.addTableView(t1.append(t2));
