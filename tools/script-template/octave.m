#name: #{NAME}
#description: Following template calculates number of cells in table
#language: octave
#tags: template, demo
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
count = rows(table) * columns(table)
