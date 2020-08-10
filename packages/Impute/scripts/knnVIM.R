#name: knnImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int k = 10 [number of nearest neighbours]
#output: dataframe outputDF [imputed dataset]

require(VIM)

Xcolnames <- colnames(X)
outputDF <- VIM::kNN(data = X, variable = Xcolnames, k = 10, trace = F, imp_var = F)

