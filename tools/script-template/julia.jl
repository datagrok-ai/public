#name: #{NAME}
#description: Following template calculates number of cells in table
#language: julia
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
count = size(table)[1] * size(table)[2]
