//name: #{NAME}
//description: Following template calculates number of cells in table
//language: nodejs
//sample: cars.csv
//input: dataframe table [Data table]
//output: int count [Number of cells in table]
const count = table.rowCount * table.columns.length;
