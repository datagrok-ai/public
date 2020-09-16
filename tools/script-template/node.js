//name: Template
//description: Following template calculates number of cells in table
//language: nodejs
//tags: template, demo
//sample: cars.csv
//input: dataframe table [Data table]
//output: int count [Number of cells in table]
let count = table.dim()[0] * table.dim()[1];
