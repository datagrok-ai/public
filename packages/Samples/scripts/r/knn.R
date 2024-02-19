#name: KNN R
#description: Imputes (numerical) missing values using the kNN algorithm
#reference: https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm
#language: r
#tags: demo
#sample: demog.csv
#input: dataframe data [Input data table with NA elements]
#input: column_list columns [Input data table columns]
#input: int neighbours = 5 [Number of Nearest Neighbours used]
#output: dataframe data_out {action:replace(data)} [Output data table without NA elements]

require(VIM)

data_out = kNN(data[columns], k=neighbours)[,columns,drop=FALSE]
