#name: knnImpl
#description: uses mice alrorithm to perform multiple imputation
#language: r
#tags: template, demo
#input: dataframe data [Input data table]
#input: int k = 10 [number of nearest neighbours]
#output: dataframe imputedDF [imputed dataset]

require(VIM)

Xcolnames <- colnames(data)
imputedDF <- VIM::kNN(data = data, variable = Xcolnames, k = k, trace = F, imp_var = F, numFun = 'medianMixed')

