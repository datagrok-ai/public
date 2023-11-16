// Hides completely empty columns

let table = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    ,         1.4,    15000
Honda, Accord,   ,         1.8,    20000
BMW,   328i,     ,         1.7,    60000        
BMW,   535i,     ,         1.5,    35000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);

let tv = grok.shell.addTableView(table);

setTimeout(() => {
  for (let c of tv.dataFrame.columns) {
    if (c.stats.missingValueCount === tv.dataFrame.rowCount)
      tv.grid.col(c.name).visible = false;
  }
}, 1000);