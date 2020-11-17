#name: knnImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: column_list columns
#input: int k = 10 [number of nearest neighbours]
#output: dataframe imputedDF [imputed dataset]

require(VIM)
require(gdata)
require(dplyr)

data<-data[rowSums(is.na(data)) != ncol(data), ]
Xcolnames <- colnames(data)
imputedDF <- VIM::kNN(data = data, variable = Xcolnames, k = k, trace = F, imp_var = F, numFun = median)

imputedDF <- imputedDF[,columns]

