#name: Sample pandas input
#description: Calculates number of cells in the table
#language: python
#input: dataframe table [Data table, in pandas]
#output: int count [Number of cells in table]
#meta.queueName: python_docker

count = table.shape[0] * table.shape[1]
