#name: #{NAME}
#description: Following template calculates number of cells in table
#language: r
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
count <- nrow(table) * ncol(table)
