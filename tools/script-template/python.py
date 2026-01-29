#name: #{NAME}
#description: Following template calculates number of cells in table
#language: python
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
count = table.shape[0] * table.shape[1]
